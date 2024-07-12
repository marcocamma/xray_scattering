import pathlib
import numpy as np
import urllib3
import os

from pdb_reader import download_pdb

try:
    import pdb_reader
except ImportError:
    print("Could not import pdb_reader")


PATH_SCRIPT = pathlib.Path(__file__).parent


def list_xyz():
    folder = PATH_SCRIPT / "xyz"
    files = list(folder.glob("*.xyz"))
    files.sort()
    names = [f.stem for f in files]
    return names


def read_xyz(fname, skip_lines=1):
    fname = pathlib.Path(fname)
    if not fname.is_file():
        fname = PATH_SCRIPT / "xyz" / f"{fname.stem}.xyz"
    with open(fname, "r") as f:
        lines = f.readlines()
    atom_list = []
    coords = []
    for line in lines[skip_lines:]:
        temp = line.split()
        atom_list.append(temp[0])
        coords.append([float(x) for x in temp[1:4]])
    return atom_list, np.array(coords)


def read_pdb(fname):
    fname = pathlib.Path(fname)
    if not fname.is_file():
        fname = PATH_SCRIPT / "pdb" / f"{fname.stem}.pdb"
    pdb = pdb_reader.PDB(fname)
    atom_list = np.asarray(pdb.get_atom_types()).astype(str)
    coords = pdb.get_coords()
    return atom_list, coords


if __name__ == "__main__":
    pass
