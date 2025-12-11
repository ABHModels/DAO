# DAO: A New and Public Non-Relativistic Reflection Model

We present a new non-relativistic reflection model, DAO, designed to calculate reflection spectra in
the rest frame of accretion disks in X-ray binaries and active galactic nuclei. The model couples the
XSTAR code, which treats atomic processes, with the Feautrier method for solving the radiative transfer
equation. A key feature of DAO is the incorporation of a high-temperature corrected cross section and
an exact redistribution function to accurately treat Compton scattering. Furthermore, the model
accommodates arbitrary illuminating spectra, enabling applications across diverse physical conditions.
We investigate the spectral dependence on key physical parameters and benchmark the results against
the widely used reflionx and xillver codes.

[![Documentation Status](https://readthedocs.org/projects/dao-model/badge/?version=latest)](https://dao.readthedocs.io/en/latest/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

## 📚 Documentation

For detailed installation instructions, usage guides, and theoretical background, please visit our official documentation:

👉 **[DAO Documentation](https://dao.readthedocs.io/en/latest/)**

---

## 🚀 Quick Start

### 1. Clone the Repository

```bash
git clone [https://github.com/ABHModels/DAO.git](https://github.com/ABHModels/DAO.git)
cd DAO
```

### 2.Prepare the Build Directory
```bash
mkdir mods
```

### 3.Configure HEASoft Paths

Open the `Makefile` and locate the HEADS variable. You must update this path to point to your local HEASoft installation directory.

*Example*:
```bash
HEADS = /path/to/heasoft/heasoft-6.31.1/x86_64-pclinux-gnu-libc2.17
```
### 4.Complie the Model
```bash
make
```