# dft-mixing-plugin
PSI4 Plugin to Support Functional Mixing

This repository consists of two independent plugins made for PSI4, the functional database management system (for better replication), and linear combination support for PSI4.

## Prerequisites
### For Anaconda
------
Install the following prerequisites
```
conda install -c conda-forge \
    sqlalchemy \
    sqlite3 \
    pylibxc \
    psi4 \
```

## Installation Instructions
1. Clone this repository `git clone https://github.com/UBCC3/dft-mixing-plugin`.
2. Source `install.sh` in the project root directory. 
```
source install.sh
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

#### Initializing the Database
---

To use initialize database, you need a config file `db_config.init.yaml` (see [sample](./example_usage/compile_database/db_init_config.init.yaml))

If there is no database at that file's location, the module will attempt to load in functionals from PSI4 and DFTD4 as a preset, and generate a new `.db` file for the database. Otherwise, the module
will just load the preexisting database at the specified path.

Afterwards, you can load in multifunctionals, base dispersion parameters, and also dispersion mixing configurations in the `preload` field of the configuration.

You may also call the adapter's method to dynamically
add data to the database.

Source Resolution
---
For any dispersion configuration/multifunctional, all sources
for composing dispersions/functionals are assumed to be the same 
as the mixed one.

By default, the source resolution order are illustrated by the diagram below:
- placeholder for user defined sources
- dftd3
- psi4 

Note that two different user defined sources cannot reference each other's parameter resolution.

How to use the database
---

### Populating the Database (Before Running PSI4)
The database supports three input files/formats:
1. Multifunctional Data (contains composing functionals and coefficients)
2. Base Dispersion Data (internal dispersion parameters)
3. Multi Dispersion Data (coefficients for multidispersion)

See example usages/format JSON files [here](../example_usage). The JSON files contain the structure
that the database expects from the data. Note that you may need to change some paths in the 
yaml configuration files to match the files. 

Populating the database is in `compile_database`,
database is stored in `database`, and a query
example is in `query_database`

To run the example, run `python3 init_db.py` from within `compile_database`  for populating the database (this should put the database
in the `database` folder. Then to run a query example, run `psi4 input.dat` in the `query_database` folder.

### Querying from the Database

In Python/PSIthon:
```python
...
import psi4
from database import Psi4DbAdapter

# Construct configuration
config_path = "path/to/config.yaml"

db = Psi4DbAdapter()
db.bind_database(config_path)

# Get functional (with dispersion) from database
functional_dict = db.get_functional_dict('BLYP', 'dispersion_config', 'func_source', 'disp_source')

# To use within PSI4
scf_energy = psi4.energy('scf_lcom', dft_functionals=functional_dict)
```
This queries for both functional (with source func_source) and dispersion (disp_source), and if the database
cannot find the respective sources, it will resort to the backup sources in order.
(By default, the order is `dftd3`, then `psi4`)







