Database Interface for Functionals and Multifunctional
---

To use the database, you need a config file `db_config.yaml`, with the following information:
```
db_path: path/to/sqlitedb/file

# Sources with lower index means higher priority.
external_resol: ['src1', 'src2']

# Optional if not going to use dftd3.
dftd3_config: "/home/kenrickmh/CHEM/dft-mixing-plugin/database/sdftd3_ref_params.toml"
load_dftd3: True
```

If there is no database at that file's location, the module will attempt to load in functionals from PSI4 and DFTD4 as a preset, and generate a new `.db` file for the database. Otherwise, the module
will just load the preexisting database at the specified path.

Afterwards, you can load in multifunctionals, base dispersion parameters, and also dispersion mixing configurations. 

You may also query a functional 


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


How to use this module
---

### Populating the Database (Before Running PSI4)
```python
...
import psi4
from database import Psi4DbAdapter

# Construct configuration
config = "path/to/config.yaml"

db = Psi4DbAdapter(config)
multifunc_dict = {
    "multifunc1": {
        "functionals": {
            "blyp": 0.25,
            "HCTH": 0.1,
            "PBE": 0.65
        }  
    },

    "multifunc2": {
        "functionals": {
            "blyp": 0.25,
            "B3LYP": 0.75
        },
        "citation": "Citation",
        "description": "description"
    },

    "multifunc3": {
        "functionals": {
            "BP86": 0.25,
            "PW91": 0.75
        }
    }    
}

db.load_multi_functional_data(multifunc_dict, "src1")
```
### Querying from the Database

In Python/PSIthon:
```python
...
import psi4
from database import Psi4DbAdapter

# Construct configuration
config = "path/to/config.yaml"

db = Psi4DbAdapter(config)

# Get functional (with dispersion) from database
functional_dict = db.get_functional_dict('BLYP', 'dispersion_config', 'func_source', 'disp_source')

# To use within PSI4
scf_energy = psi4.energy('scf_lcom', dft_functionals=functional_dict)
```
This queries for both functional (with source func_source) and dispersion (disp_source), and if the database
cannot find the respective sources, it will resort to the backup sources in order.
(By default, the order is `dftd3`, then `psi4`)

