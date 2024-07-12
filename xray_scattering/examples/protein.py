import numpy as np
from matplotlib import pyplot as plt

from xray_scattering import io
from xray_scattering import debye
from xray_scattering import photonize



water = debye.read_from_database("water_narten")
air = debye.read_from_database("air")

q = np.linspace(0,3,1001)

atoms,coords=io.read_pdb("1bbb")
protein=debye.calc_debye_distribution(atoms,coords,q, solvent_electron_density=0)#.334)



det = photonize.Detector(npixel=2000, pixel_size=75e-6)

sample_det_dist = -0.15

water_c = photonize.liquid_scattering(
    water,
    det,
    sample_det_dist=sample_det_dist,
    sample_thickness=1e-3,
    wavelength=1,
    density_g_cm3=1,
    do_2d=False,
)

prot_c = photonize.liquid_scattering(
    protein,
    det,
    sample_det_dist=sample_det_dist,
    sample_thickness=1e-3,
    wavelength=1,
    molar_density=1e-3,
    do_2d=False,
)


air_c = photonize.gas_scattering(
    air, det, xrange=(-0.05+sample_det_dist, 0.05+sample_det_dist), N=20, wavelength=1, do_2d=False
)

plt.title("1e12 ph on incident beam, {sample_det_dist} detector distance")
plt.semilogy(prot_c, label="protein (1mM, 1mm)")
plt.semilogy(air_c, label="air (10cm)")
plt.semilogy(water_c, label="water (1mm)")
plt.grid(True)
plt.legend()
plt.xlabel("pixel")
plt.ylabel("counts / pixel")
