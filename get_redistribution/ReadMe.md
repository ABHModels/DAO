# Physics Simulation Code

## Overview
This code implements numerical routines for the **Compton scattering kernel** originally developed by García (2020).  

We extend the implementation to include a **Gaussian-approximated kernel** following Ross (1978).  

Different kernel types can be selected at runtime.

## Usage
Compile with `make`:
```bash
make
```

**Please update the HEADS variable in the Makefile to match your local HEASoft installation path.**

```
./drive_SRF.x [KERNEL_TYPE]
```

**KERNEL_TYPE**

11: Exact kernel function (García 2020)

22: Gaussian-approximated kernel (Ross 1978)