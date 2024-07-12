#! /usr/bin/python
import copy
from datastorage import DataStorage as ds
import numpy as np
from scipy import ndimage
import pathlib
try:
    import sr
    HAS_SR = True
except ImportError:
    HAS_SR = False

# constants
c_NA = 6.022e23
# Avogadro number
c_re = 2.817e-15
# classical electron radius, in meter
g_kb = 1.3806485e-23
# elementary charge
g_e = 1.6e-19

PATH_SCRIPT = pathlib.Path(__file__).parent


class Detector(object):
    def __init__(
        self,
        pixel_size=100e-6,
        npixel=1000,
        material="Si",
        thickness=300e-6,
        center="auto",
    ):
        """
        pixel_size in m
        """
        self.pixel_size = pixel_size
        self.npixel = npixel
        if isinstance(center, str) and center == "auto":
            xc = yc = npixel / 2
        else:
            xc, yc = center
        tempx = (np.arange(npixel) - xc) * pixel_size
        tempy = (np.arange(npixel) - yc) * pixel_size
        self.x, self.y = np.meshgrid(tempx, tempy)
        self.r = np.sqrt(self.x**2 + self.y**2)
        self.azimuth = np.arctan2(self.y, self.x)
        self.thickness = thickness
        self.material = material

        self.x1D = np.arange(npixel / 2) * pixel_size

        self._precomputed = dict()

    def prepare(self, wavelength=1, dist=-0.1, do_2d=True):
        """
        dist = position from detector (negative means detector is downbeam with respect to the sample) [ in m ]
        """
        temp = ds()
        self.wavelength = wavelength
        self.energy = 12.398 / wavelength
        if do_2d:
            angle = np.arctan2(self.r, -dist)
        else:
            angle = np.arctan2(self.x1D, -dist)
        temp.scattering_angle = angle
        temp.q = 4 * np.pi / self.wavelength * np.sin(angle / 2)
        # cos**3
        temp.geom_factor = np.abs(np.cos(angle)) ** 3
        if do_2d:
            temp.pol_factor = (
                0.5
                - 0.5 * np.sin(angle) ** 2 * np.cos(2 * self.azimuth)
                + 0.5 * np.cos(angle) ** 2
            )
        else:
            temp.pol_factor = 1
        temp.total_i_factor = (
            (self.pixel_size / np.abs(dist)) ** 2 * temp.geom_factor * temp.pol_factor
        )
        temp.thickness = self.thickness / np.abs(np.cos(angle))
        if HAS_SR:
            temp.det_absorption = sr.materials.absorption(
                material=self.material, thickness=temp.thickness, energy=self.energy
            )
        else:
            temp.det_absorption = 1
        self.data = temp
        return temp


def photons_to_data(
    debye,
    det,
    wavelength=1,
    dist=-0.1,
    do_2d=True,
    consider_det_abs=True,
    beam_size=None,
):
    """
    I, intensity (in e2)
    dist = position from detector [in m]
         negative means that detector is downbeam with respect to the sample)
    if consider_det_abs, it takes into account finite (and angular dependent) absorption
    beamsize is the fwhm beamsize. if not None, convolution is used to mimic finite beamsize
    """
    det_data = det.prepare(dist=dist, wavelength=wavelength, do_2d=do_2d)
    data = np.interp(det_data.q, debye.q, debye.intensity_photons)
    data *= det_data.total_i_factor
    if consider_det_abs:
        data *= det_data.det_absorption

    if beam_size is not None:
        beam_size = beam_size / det.pixel_size
        if beam_size > 0.5:
            print(
                "Convoluting to mimic finite beamsize; kernel size in pixels", beam_size
            )
            data = ndimage.filters.uniform_filter(data, beam_size)
    return data


def gas_scattering(
    debye,
    detector,
    I0=1e12,
    wavelength=1,
    T=300,
    P=10e5,
    xrange=(-30, -0.5),
    N=100,
    do_2d=True,
    consider_det_abs=True,
    beam_size=None,
):
    """
    P is pressure in Pascal (1atm = 10^5 Pa)
    T is the temperature (in K)
    xrange = beginning and end of region with gas [in m]
        negative values mean that detector is downbeam with respect to the sample)
    NOTE: it does *not* include transmission losses
    """
    # in molecule per m3
    density = P / (g_kb * T)
    x = np.linspace(xrange[0], xrange[1], N)
    dx = abs(x[1] - x[0])
    if do_2d:
        Isum = np.zeros_like(detector.x)
    else:
        Isum = np.zeros_like(detector.x1D)
    debye.intensity_photons = I0 * c_re**2 * density * dx * debye.intensity_total
    for i, _x in enumerate(x):
        temp = photons_to_data(
            debye,
            detector,
            wavelength=wavelength,
            dist=_x,
            do_2d=do_2d,
            consider_det_abs=consider_det_abs,
            beam_size=beam_size,
        )
        Isum += temp
    return Isum


def liquid_scattering(
    debye,
    detector,
    I0=1e12,
    wavelength=1,
    density_g_cm3=None,
    molar_density=None,
    sample_det_dist=-0.1,
    sample_thickness=1e-3,
    do_2d=True,
    consider_det_abs=True,
    beam_size=None,
):
    """
    if molar_density is given, density_gm_cm3 is discarded
    sample_det_dist = position from detector [in m]
        negative means that detector is downbeam with respect to the sample)
    sample_thickness = thickness [in m]
    """
    try:
        atoms_string = "".join(debye.atoms)
    except TypeError:
        atoms_string = "".join(debye.atoms.astype(str))

    if molar_density is None and HAS_SR:
        molar_density = sr.materials.get_molar_density(
            atoms_string, density=density_g_cm3
        )
        t = sr.materials.transmission(
            material=atoms_string,
            thickness=sample_thickness,
            energy=12.398/wavelength,
            density=density_g_cm3,
        )
#        print("Transmission", t)

    else:
        t=1
#        print("Neglecting transmission")

    rho_mol_m3 = molar_density * 1e3 * c_NA  # 1e3 to go from L to m3

    debye.intensity_photons = (
        t * I0 * c_re**2 * rho_mol_m3 * sample_thickness * debye.intensity_total
    )
    return photons_to_data(
        debye,
        detector,
        wavelength=wavelength,
        dist=sample_det_dist,
        do_2d=do_2d,
        consider_det_abs=consider_det_abs,
        beam_size=beam_size,
    )


def photons_to_current(I, energy=10, bandgap=0.0036):
    """I is in (absorbed)photons/sec"""
    return I * energy / bandgap * g_e
