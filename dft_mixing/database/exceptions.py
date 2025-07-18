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
    
# These are functionals whose canonical name cannot be separated into 
# dashcoeff nomeclature
weird_dashcoeff_functionals : dict = {
    "b973c": {
        "name": "B973c",
        "alias": [
            "B97-3c"
        ],
        "xc_functionals": {
            "GGA_XC_B97_3C": {}
        },
        "description": "    B97-3c GGA-based 3C composite method with a TZ basis set, D3 and short-range basis set correction.\n",
        "citation": "     J. G. Brandenburg, C.Bannwarth, A. Hansen, S. Grimme J. Chem. Phys. 148, 064104, 2018\n",
        "doi": "10.1063/1.5012601",
        "dispersion": {
            "type": "d3bjatm",
            "params": {
                "s6": 1.0,
                "s8": 1.5,
                "a1": 0.37,
                "a2": 4.1,
                "s9": 1.0
            },
            "citation": False
        }
    },
    
    "b97-3c": {
        "name": "B973c",
        "alias": [
            "B97-3c"
        ],
        "xc_functionals": {
            "GGA_XC_B97_3C": {}
        },
        "description": "    B97-3c GGA-based 3C composite method with a TZ basis set, D3 and short-range basis set correction.\n",
        "citation": "     J. G. Brandenburg, C.Bannwarth, A. Hansen, S. Grimme J. Chem. Phys. 148, 064104, 2018\n",
        "doi": "10.1063/1.5012601",
        "dispersion": {
            "type": "d3bjatm",
            "params": {
                "s6": 1.0,
                "s8": 1.5,
                "a1": 0.37,
                "a2": 4.1,
                "s9": 1.0
            },
            "citation": False
        }
    },
    
    "wb97x3c": {
        "name": "wB97X3c",
        "alias": [
            "wB97X-3C"
        ],
        "xc_functionals": {
            "HYB_GGA_XC_WB97X_V": {}
        },
        "dispersion": {
            "type": "d4bjeeqatm",
            "nlc": False,
            "params": {
                "a1": 0.2464,
                "a2": 4.737,
                "s6": 1.0,
                "s8": 0.0,
                "s9": 1.0,
                "alp": 16.0,
                "ga": 3.0,
                "gc": 2.0,
                "wf": 6.0
            },
            "citation": False
        },
        "description": "    wB97X basied 3C composite method with a small basis set, gCP and D4\n",
        "citation": "    M. Muller, A. Hansen, S. Grimme, J. Chem. Phys. 158, 014103 (2023)\n",
        "doi": "10.1063/5.0133026"
    },
    
    
}

# Dummy functionals 
# 
dummy_functionals : dict = {
    
    
    
    
    
}