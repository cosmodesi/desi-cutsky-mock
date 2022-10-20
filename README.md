# Generate Survey Mocks

## Table of Contents

-   [Introduction](#introduction)
-   [Dependencies](#Dependencies)
-   [Usage](#Usage)
-   [References](#references)
-   [Contributors][#Contributors]
-   [Acknowledgements](#Acknowledgements)

## Introduction

**G**enerate **S**urvey **M**ocks (GSM) is a tool for generating Light-Cone like mocks starting from boxes. The principle behind the code is to repeat the box, convert from (x,y,z) coordinates to (RA, DEC, Z) and then cut the necessary survey volume. Of course, for cubic mocks that do not have a large enough volume to encapsulate the entire survey volume, some of the regions of the box are used multiple times. However, if the box is sufficiently large and it is rotated accordingly, the survey volume can be composed from unique regions of the box. 

The code accepts a rotation matrix which is defined in ([rotation_matrix.py](rotation_matrix.py)).


If you use this tool in research work that results in publications, please cite the following paper:

> Variu et al., in preparation.

## Dependencies

The dependencies of this code are as follows:

-   [Python3](https://www.python.org/)  (>= 3.6)
-   [NumPy](https://numpy.org/)
-   [argparse]
-   [configparser]
-   [sys]
-   [os]
-   [itertools]
-   [glob]
-   [multiprocessing]
-   [astropy](https://www.astropy.org/)
-   [healpy](https://healpy.readthedocs.io/en/latest/index.html)
-   [numexpr](https://github.com/pydata/numexpr)
-   [h5py](https://docs.h5py.org/en/stable/)
-   [fitsio](https://github.com/esheldon/fitsio)
-   [desimodel](https://github.com/desihub/desimodel)
-   [camb](https://github.com/cmbant/CAMB)

## Usage

After downloading the code, one has to fill the configuration file and use ore create a function in ([main.py](main.py)) that is specific to data structure of the boxes. There are already some written functions in ([main.py](main.py)) that are dedicated to the First Generation of Mocks for DESI based on EZmocks and ABACUS simulations. Finally, one runs:

```bash
python main.py
```


## Configuration parameters

Part of the configuration parameters can be given as arguments when running the code. They can be displayed via the `-h` option as follows:

```bash
python main.py -h
```

Most of them, however, should be given in the configuration file. Examples of configuration files can be found in the ([EZmock/config](EZmock/config)) and ([ABACUS/config](ABACUS/config)) folders.

## Contributors

I thank [Dr. Shadab Alam](https://github.com/shadaba) for his support and suggestions. 

## Acknowledgements
This code is inspired by [prcaetano's](https://github.com/prcaetano/gallightcone) and [Dr. Chia-Hsun Chuang's](https://github.com/chia-hsun-chuang/apply_desifootprint_nz) packages.