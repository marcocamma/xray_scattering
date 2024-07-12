
# xray_scattering

This simple project aims at calculating the scattering from liquids, gases (including distributed sources) and molecules in solutions. In particular by using so called "electron units", the number of photons expected on a detector can be calculated.
This can be used to estimate signal to noise, level of background signals, etc.

*Be aware*: this project is un-tested

The overall idea is as follow:

1. The calculated (or measured) intensity is read or calculated in so called (electron units). Two modules are used for that:

  - io (to read xyz and pdb files)
  - debye (to calculate the intensity from the atom list and the coordinates)

2. A "detector" is defined (pixel size, wavelength, etc.)

3. The expected signal on the detector is calculated for a given scatterer (i.e. debye output) and detector configuration

## Create virtual environment (optional)
```
# better to create a virtual environment ?
cd your_favourite_folder
python -m venv xray_scattering_venv
source xray_scattering_venv/bin/activate
```

## Install
```
pip install git+https://github.com/marcocamma/xray_scattering.git
```

## Usage

### run examples
```py
from xray_scattering import examples
```

### xray_scattering.io

```py
from xray_scattering import io

# read atoms and coordinates
atoms,coordinates = io.read_xyz(filename)

# list pre-defined molecules
io.list_xyz()

# read pdb file
io.download("1bbb")
atoms,coordinates = io.read_pdb("1bbb") # pdbid or filename
```

### xray_scattering.debye

```py
# for small molecule
data = debye.calc_debye(atoms, coords)

# for larger cluster, calculates pair distribution funtions (much faster)
          
data = debye.calc_debye_distribution(atoms,coords)

# interact with 'database'
debye.list_database()
debye.save_in_database(data,name)
data = debye.read_from_database(name)

# include effect of Sq
data_sq = debye.apply_Sq(data,fname_with_sq)

# include (approximate) excluded volume effect (important for proteins)
data = debye.calc_debye(atoms, coords, solvent_electron_density=0.334)

# "sum" scattering
n2 = debye.read_from_database("nitrogen")
o2 = debye.read_from_database("oxygen")
argon = debye.read_from_database("argon")

air = 0.79*n2+0.2*o2+0.01*argon
air.save_in_database("air")

```

### xray_scattering.photonize
```py
detector = photonize.Detector(npixel=2000,pixel_size=100e-6)

image = photonize.liquid_scattering(
def liquid_scattering(
    data,
    detector,
    I0=1e12,
    wavelength = 1,
    density_g_cm3=None,
    molar_density=None,
    sample_det_dist =-0.1, # neg means sample before det
    sample_thickness=1e-3,
    do_2d=True,
    consider_det_abs=True,
    beam_size=None,
)


# note that gas scattering is considered as 'extended source'.
# the scattering volme is sub-divided into N equally spaces scatterers
image = photonize.gas_scattering(
    data,
    detector,
    I0=1e12,
    wavelength = 1,
    T=300,
    P=10e5,              # pressure in Pa
    xrange =(-30,-0.5),
    N=100,
    do_2d=True,
    consider_det_abs=True,
    beam_size=None,
)
```
