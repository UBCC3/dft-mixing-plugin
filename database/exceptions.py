'''
    List of exceptional functionals that 
    somehow **DO NOT** have a parent DFT
    but do have dispersion models. 
    
    This has to be patched in manually...
'''

import json, os

# For now, just blacklist them, still need to 
# patch them manually...

# And for aliases, just insert "duplicate" names into 
# the mix

json_path = os.path.join(os.path.dirname(__file__), 'problematic.json')
with open(json_path, 'r') as f:
    blacklist_funcs = json.load(f) 