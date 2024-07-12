""" 
Calcutes scattered intensities (in electron units) using Debye equation

This provides:
  - a very good approaximation for molecules in gas phases
  - an OK approaximation for proteins in solutions (more specialized programs
    exist (e.g. crysol, freesaxs, etc.)
      - note that excluded volume effects can be taken into account (if solvent_electron_density>0)
        but is an approximate way (CH, CH2, CH3 are all considered C in terms of excluded volumes)
  - for condensed phases (liquid, glasses, etc.) the so called structure factor
      must be included
    - this can be done by calculating a large "droplets" or by more advanced
      techniques (using partials g(r) from MD) that are outside the scope of this
      very simple module

It uses xraylib to calculates atomic form factors

"""
#! /usr/bin/python
from functools import lru_cache
import numpy as np
import scipy.spatial
import xraylib
import pathlib

import copy

from datastorage import DataStorage as ds
import datastorage

PATH_SCRIPT = pathlib.Path(__file__).parent

ATOM_VOLUME = dict()
ATOM_VOLUME["H"] = 5.15
ATOM_VOLUME["O"] = 9.13
ATOM_VOLUME["C"] = 16.44
ATOM_VOLUME["N"] = 2.49
ATOM_VOLUME["S"] = 19.86


def list_database():
    folder = PATH_SCRIPT / "database/"
    files = list(folder.glob("*.h5"))
    files.sort()
    names = [f.stem for f in files]
    return names


def read_from_database(name):
    folder = PATH_SCRIPT / "database"
    name = pathlib.Path(name).stem
    fname = folder / f"{name}.h5"
    return from_file(fname)


def from_file(fname):
    keys = "atoms", "q", "intensity", "intensity_compton"
    d = datastorage.read(fname)
    return Debye(**{k: d[k] for k in keys})


def atomic_volume(atom):
    if isinstance(atom,int): atom = xraylib.AtomicNumberToSymbol(atom)
    if atom in ATOM_VOLUME:
        return ATOM_VOLUME[atom]
    else:
        return 20


@lru_cache
def _fq(Z, q, solvent_electron_density=0):
    """note that for X-ray lib drops the 4π factor in q"""
    q_xraylib = q / 4 / np.pi
    fq = xraylib.FF_Rayl(Z, q_xraylib)

    if solvent_electron_density != 0:
        v = atomic_volume(Z)
        fq -= solvent_electron_density * v * np.exp(-(q**2) * v ** (2 / 3) / 4 / np.pi)
    return fq


@lru_cache
def _compton(Z, q):
    """note that for X-ray lib drops the 4π factor in q"""
    q_xraylib = q / 4 / np.pi
    try:
        return xraylib.SF_Compt(Z, q_xraylib)
    except ValueError:  # ( for q<0.1)
        return 0


def _to_str(a):
    if isinstance(a, str):
        return a
    else:
        return a.decode()


class Debye:
    def __init__(self, atoms, q, intensity, intensity_compton):
        atoms = [_to_str(a) for a in atoms]
        self.atoms = np.asarray(atoms).astype(str)
        self.q = q
        self.intensity_compton = intensity_compton
        self.intensity = intensity
        self.intensity_total = intensity + intensity_compton

    def __add__(self, other):
        if np.all(self.q == other.q):
            i = self.intensity + other.intensity
            c = self.intensity_compton + other.intensity_compton
        else:
            i = self.intensity + np.interp(self.q, other.q, other.intensity)
            c = self.intensity_compton + np.interp(
                self.q, other.q, other.intensity_compton
            )
        atoms = list(self.atoms) + list(other.atoms)
        return Debye(
            atoms=atoms,
            q=self.q,
            intensity=i,
            intensity_compton=c,
        )

    def __mul__(self, f):
        return Debye(
            atoms=[f"{f}*{a}" for a in self.atoms],
            q=self.q,
            intensity=f * self.intensity,
            intensity_compton=f * self.intensity_compton,
        )

    __rmul__ = __mul__

    def save(self, fname):
        d = ds(
            atoms=self.atoms,
            q=self.q,
            intensity=self.intensity,
            intensity_compton=self.intensity_compton,
            intensity_total=self.intensity_total,
        )
        d.save(fname)

    def save_in_database(self, key):
        fname = PATH_SCRIPT / "database" / f"{key}.h5"
        return self.save(fname)

    def __repr__(self):
        a_types = np.unique(self.atoms)
        N = {}
        for a in a_types:
            N[a] = np.count_nonzero(self.atoms == a)
        return f"Debye {str(N)}"


def ff_coherent(symbol, q, solvent_electron_density=0):
    z = xraylib.SymbolToAtomicNumber(symbol)
    if isinstance(q, (float, int)):
        return _fq(z, q, solvent_electron_density=solvent_electron_density)
    else:
        return np.asarray(
            [_fq(z, qi, solvent_electron_density=solvent_electron_density) for qi in q]
        )


def iq_compton(symbol, q):
    z = xraylib.SymbolToAtomicNumber(symbol)
    if isinstance(q, (float, int)):
        return _compton(z, q)
    else:
        return np.asarray([_compton(z, qi) for qi in q])


def mysinc(x):
    x = np.asarray(x)
    sinc = np.empty_like(x)
    idx = x == 0.0
    sinc[~idx] = np.sin(x[~idx]) / x[~idx]
    sinc[idx] = 1.0
    return sinc


def apply_Sq(debye, sq_file):
    debye = copy.deepcopy(debye)
    if isinstance(sq_file, str):
        sq_file = pathlib.Path(sq_file)
    if not sq_file.is_file():
        sq_file = PATH_SCRIPT / "database" / f"{sq_file}_Sq.txt"
    print("Reading Sq file", str(sq_file))
    q_sq, sq = np.loadtxt(sq_file, usecols=(0, 1), unpack=True)
    sq = np.interp(debye.q, q_sq, sq)
    debye.intensity *= sq
    debye.intensity_total *= sq
    return debye


def coords_to_distance_matrix(posList, posList2=None):
    """0.24 sec with 4400 atoms, 2 sec for 15000"""

    posList = np.asarray(posList)
    if posList2 is None:
        dist = scipy.spatial.distance.cdist(posList, posList)
    else:
        dist = scipy.spatial.distance.cdist(posList, posList2)
    return dist


def calc_debye(atoms, coords, q=np.arange(0, 0.5, 0.01), solvent_electron_density=0):
    atoms = np.asarray(atoms)
    atoms_unique = np.unique(atoms)
    coords = np.asarray(coords)
    if coords.ndim == 1:
        coords = coords[np.newaxis, :]
    fq_coh = dict(
        (a, ff_coherent(a, q, solvent_electron_density=solvent_electron_density))
        for a in atoms_unique
    )
    iq_com = dict((a, iq_compton(a, q)) for a in atoms_unique)
    nAtoms = len(atoms)
    intensity = np.zeros_like(q)
    intensity_compton = np.zeros_like(q)
    dist_matrix = coords_to_distance_matrix(coords, coords)
    for i in range(nAtoms):
        for j in range(nAtoms):
            r_ij = dist_matrix[i, j]
            if r_ij == 0:
                sinc = 1
            else:
                sinc = mysinc(r_ij * q)
            intensity += fq_coh[atoms[i]] * fq_coh[atoms[j]] * sinc
            if i == j:
                intensity_compton += iq_com[atoms[i]]  # * fq_com[atoms[j]]
    return Debye(
        atoms=atoms,
        q=q,
        intensity=intensity,
        intensity_compton=intensity_compton,
    )


def read_qi_file(fname, atoms):
    """read file with two columns, q,i"""
    q, intensity = np.loadtxt(fname, usecols=(0, 1), unpack=True)
    return ds(
        q=q,
        intensity=intensity,
        intensity_compton=np.zeros_like(intensity),
        atoms=atoms,
        intensity_total=intensity,
    )


def calc_debye_distribution(
    atoms,
    coords,
    q=np.arange(0, 0.5, 0.01),
    order=2,
    N=1000,
    dmax=100,
    solvent_electron_density=0,
):
    atoms = np.asarray(atoms)
    atoms_unique = np.unique(atoms)
    coords = np.asarray(coords)
    if coords.ndim == 1:
        coords = coords[np.newaxis, :]
    fq_coh = dict(
        (a, ff_coherent(a, q, solvent_electron_density=solvent_electron_density))
        for a in atoms_unique
    )
    iq_com = dict((a, iq_compton(a, q)) for a in atoms_unique)
    # write d(n) = a*n**order
    a = dmax / (N**order)
    n = np.arange(N)
    dbin = a * n**order
    d = (dbin[:-1] + dbin[1:]) / 2.0
    # prepare matrix of q * d
    qd = [qi * d for qi in q]
    qd = np.asarray(qd)
    sync = mysinc(qd)
    intensity = []
    intensity_compton = []
    for i1 in range(len(atoms_unique)):
        for i2 in range(i1, len(atoms_unique)):
            a1 = atoms_unique[i1]
            a2 = atoms_unique[i2]
            idx1 = atoms == a1
            idx2 = atoms == a2
            dMat = coords_to_distance_matrix(coords[idx1], coords[idx2])
            idx = np.digitize(dMat.ravel(), dbin)
            Nbins = np.bincount(idx, minlength=len(d))[: len(d)]
            I12 = np.sum(Nbins * sync, axis=1)
            if i1 != i2:
                I12 *= 2  # because i1>=i2 (we are not calculating simmetric)
            else:
                intensity_compton.append(iq_com[a1] * len(a1))
            # multiply by form factor
            I12 *= fq_coh[a1] * fq_coh[a2]
            intensity.append(I12)
            # print time.time()-t0
    intensity = np.array(intensity)
    intensity_compton = np.array(intensity_compton)
    intensity = intensity.sum(axis=0)
    intensity_compton = intensity_compton.sum(axis=0)
    return Debye(
        atoms=atoms,
        q=q,
        intensity=intensity,
        intensity_compton=intensity_compton,
    )


def make_simple_molecules():
    pass
    # air


if __name__ == "__main__":
    pass
