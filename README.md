# dft-mixing-plugin
PSI4 Plugin to Support Functional Mixing

This repository consists of two independent plugins made for PSI4, the functional database management system (for better replication), and linear combination support for PSI4.

## Installation Instructions

Source `scf_lcom/pathconfig.sh` to apply the changes from the plugin to PSI4. To make this permanent, you would need to add this to your `~/.bashrc` file. 
```
source scf_lcom/pathconfig.sh
```

## Plugins Overview

### 1. Linear Combination Functional and Dispersion Mixing (LCOM)

This plugin was adapted from the previous work, but made more adaptive to PSI4 by using existing constructs to mix density functionals.

In short, the plugin performs linear combination mixing of functionals, and dispersion (as post-processing after the functional mixing). 

Mathematically, given functional coefficients $c_1, c_2, ..., c_n$, and density functionals $\rho_1, \rho_2, ..., \rho_n$, we could linearly mix these density functionals, and thus the energy contribution by these functionals can be written:
$$E[\rho] = \sum_{k = 1}^n c_k E[\rho_k]$$

Likewise, in the plugin, the dispersion is computed similarly.

The plugin has the following interface:
```python
# input.dat

import scf_lcom  # Imports this plugin to enable LCOM

# Initialization code
...

# Example functional
YOUR_FUNCTIONAL = {
    "name": "your_functional_name" (str)
    "lcom_functionals": {
        # You can specify/use pre-existing PSI4 functionals as such:
        "BLYP": 1.0,

        # or define your own functional, indexed by the key.
        "custom_functional": {
            "name": "custom_functional",
            "coef": 1.0,

            # The rest is identical to a normal PSI4 functional
            # description.
            # Please refer to the PSI4 documentation for the
            # formatting.
            #
            # https://github.com/psi4/psi4/blob/master/psi4/driver/procrouting/dft/dft_builder.py
        } 
    }

    "lcom_dispersion":[
        {    
            "type": "d2",

            # Dispersion parameters, need to specify manually
            # unfortunately.

            "params": {}, 
            "lcom_coef": 0.1 # (1.0 by default)
        }
    ]
}

# Actual Call
energy('scf_lcom', dft_functional=YOUR_FUNCTIONAL)
```

Additionally, you may specify to break down LibXC XC functionals to enable individual X and C mixing for more efficient calculations (if XC functions overlap a lot). This is illustrated below:

```
energy('scf_lcom', dft_functional=..., decompose_xc=True)
```
(note the use of the `decompose_xc` flag).

### 2. Functional Database Interface

This addon works independently from PSI4, since it relies on a precompiled database of functionals. The addon aims to streamline the interface to specifying functionals, and allows for flexible user defined functionals to be used, and adds another layer of infrastructure to the system. 






