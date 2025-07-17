
from database import Psi4DbAdapter
import json
import pprint
import os

file_dir = os.path.dirname(__file__) 
path = os.path.join(file_dir, "database_config.yaml")
psi_adapter = Psi4DbAdapter(path)

multifunc_file = os.path.join(file_dir, "multifunc.json")
# Add funcitonal to database
with open(multifunc_file, "r") as f: 
    multifunc_dict = json.load(f)

base_disp_file = os.path.join(file_dir, "base_disp.json")
# Add funcitonal to database
with open(base_disp_file, "r") as f: 
    base_disp_dict = json.load(f)
    
disp_config_file = os.path.join(file_dir, "disp_config.json")
# Add funcitonal to database
with open(disp_config_file, "r") as f: 
    disp_config_dict = json.load(f)

psi_adapter.load_multi_functional_data(multifunc_dict, "src1")
psi_adapter.load_base_dispersion_data(base_disp_dict, 'src1')
psi_adapter.load_dispersion_config_data(disp_config_dict, 'src1')

# Multifunctional query
print("=" * 40 + " MULTIFUNCTIONAL QUERY " + "=" * 40)
print("MULTIFUNC1")
pprint.pprint(psi_adapter.get_functional_dict("multifunc1"))

print("MULTIFUNC2")
pprint.pprint(psi_adapter.get_functional_dict("multifunc2"))

print("MULTIFUNC3")
pprint.pprint(psi_adapter.get_functional_dict("multifunc3"))


print("MULTIFUNC1 w DISP1")
pprint.pprint(psi_adapter.get_functional_dict("multifunc1", "disp_mix1"))
