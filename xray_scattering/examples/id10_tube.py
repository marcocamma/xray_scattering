#! /usr/bin/python

from matplotlib import pyplot as plt

from xray_scattering import debye
from xray_scattering import photonize


center = (1700, 1160)
det = photonize.Detector(pixel_size=75e-6, npixel=2000, center=center)


def id10_tube(N=100, P_mbar=1, length=5.5, dist_min=0.3, wavelength=1):
    air = debye.read_from_database("air")
    P = P_mbar * 1e2
    return photonize.gas_scattering(
        air,
        det,
        I0=1,
        T=300,
        P=P,
        xrange=(-length, -dist_min),
        N=N,
        wavelength=wavelength,
    )

image = id10_tube(N=30,P_mbar=10)
plt.imshow(1e12*image)
plt.title("gas scattering (1mbar, 5 m long tube, 10^{12} ph on incident beam)")
