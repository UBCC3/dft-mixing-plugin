Database Interface for Functionals and Multifunctional
---

To use the database, you need a config file `db_config.yaml`, with the following information:
```
db_config: path/to/sqlitedb/file
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
- placeholder for user def sources
- dftd4
- psi4 

Note that two different user defined sources cannot reference each other's parameter resolution.


How to use this module
---

### Populating the Database (Before Running PSI4)
```python
...
import psi4
from scf_lcom.database import FunctionalDatabase

# Construct configuration
config = "path/to/config.yaml"

db = FunctionalDatabase(config)

# Get functional (with dispersion) from database
functional_dict = db.query_functional_disp('BLYP', 'dispersion_config', 'source', format='psi4')

# To use within PSI4
scf_energy = psi4.energy('scf_lcom', dft_functionals=functional_dict)
```



### Querying from the Database

In Python/PSIthon:
```python
...
import psi4
from scf_lcom.database import FunctionalDatabase

# Construct configuration
config = "path/to/config.yaml"

db = FunctionalDatabase(config)

# Get functional (with dispersion) from database
functional_dict = db.query_functional_disp('BLYP', 'dispersion_config', 'source', format='psi4')

# To use within PSI4
scf_energy = psi4.energy('scf_lcom', dft_functionals=functional_dict)
```

