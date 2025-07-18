'''
    Sample data for testing the database functionality.
    The samples consist of two sets, the raw data that is going 
    to be inserted into the databse, and also a set of queries and 
    correct answer to those queries
'''
from psi4.driver.procrouting.dft.dft_builder import functionals

lcom_func_dataset = {   
    "singlefunc1": {
        "functionals": {
            "BLYP": 1.0
        }   
    },  
    
    "singlefunc2": {
        "functionals": {
            "pbe0": 1.0
        }   
    },  
    
    "multifunc2": {
        "functionals": {
            "blyp": 0.25,
            "B3LYP": 0.75
        },
        "citation": "ABCD",
        "description": "ABCD",
    },

    "multifunc3": {
        "functionals": {
            "BP86": 0.25,
            "PW91": 0.75,
        }
    },    
}

fail_dataset = {
    "fail1": {
        "functionals": {
            "NOT_A_FUNC": 1.0,
        }
    },  
    
    "fail2": {
        "functionals": {
            "BLYP": 0.25,
            "NOT_A_FUNC": 1.0,
        }
    },  
}

func_dataset_ref = {
    "singlefunc1": {
        "name": "singlefunc1",
        "lcom_functionals": {
            "BLYP": 1.0
        }
    },  
    
    "singlefunc2": {
        "name": "singlefunc2",
        "lcom_functionals": {
            "pbe0": 1.0
        }   
    },  
    
    "multifunc2": {
        "name": "multifunc2",
        "lcom_functionals": {
            "BLYP": 0.25,
            "B3LYP": 0.75,
        },
    },

    "multifunc3": {
        "name": "multifunc3",
        "lcom_functionals": {
            "BP86": 0.25,
            
            # Aliased name
            "PW91": 0.75
        }
    },    
}

lcom_psi_dispconfig_dataset : dict = {
    "TPSS": {
        "disp_mix1": {
            "coeffs": {
                "d2": 0.25,
                "d3zero": 0.25,
            }
        },
        
        "disp_mix2": {
            "coeffs": {
                "d2": 0.55,
                "d3zero": 0.25,
                "d3bj": 0.15,
            }   
        },
    },  
    
    "singlefunc1": {
        "disp_mix1": {
            "coeffs": {
                "d2": 0.25,
                "d3": 0.25,
                "d4": 0.5,
            }   
        },
        
        "disp_mix2": {
            "coeffs": {
                "d2": 0.1,
                "d3": 0.2,
                "d4": 0.8,
            },
        }, 
    },  
    
    "multifunc2": {
        "disp_mix1": {
            "coeffs": {
                "d2": 0.25,
                "d3": 0.25,
                "d4": 0.5,
            },
        },
        
        "disp_mix2": {
            "coeffs": {
                "d2": 0.1,
                "d3": 0.2,
                "d4": 0.7,
            },
        }, 
    },  
}


lcom_base_disp_dataset : dict = {
    "singlefunc1": {
        "d2": {
            "type": "d2",
            "params": {
                "test_param": "test_param"
            }
        },
        
        "d3": {
            "type": "d3",
            "params": {
                "test_param": "test_param"
            }
        },
        
        "d4": {
            "type": "d4",
            "params": {
                "test_param": "test_param"
            }
        }
    },  
    
    "multifunc2": {
        "d2": {
            "type": "d2",
            "params": { "test_param": "test_param2" }
        },
        
        "d3": {
            "type": "d3",
            "params": {
                "test_param": "test_param2"
            }
        },
        
        "d4": {
            "type": "d4",
            "params": {
                "test_param": "test_param2"
            }
        }   
    },  
}

multi_disp_ans = {    
    ("multifunc2", "disp_mix1"): {
        **func_dataset_ref['multifunc2'],
        "lcom_dispersion":[
            {    
                "type": "d2",
                "params": {"test_param": "test_param2"},
                "lcom_coef": 0.25 # (1.0 by default)
            },
            
            {    
                "type": "d3",
                "params": {"test_param": "test_param2"},
                "lcom_coef": 0.25 # (1.0 by default)
            },
            
            {    
                "type": "d4",
                "params": {"test_param": "test_param2"},
                "lcom_coef": 0.5 # (1.0 by default)
            },
        
        ]
    },
    
    ("multifunc2", "disp_mix2"): {
        **func_dataset_ref['multifunc2'],
        "lcom_dispersion":[
            {    
                "type": "d2",
                "params": {"test_param": "test_param2"},
                "lcom_coef": 0.1 # (1.0 by default)
            },
            
            {    
                "type": "d3",
                "params": {"test_param": "test_param2"},
                "lcom_coef": 0.2 # (1.0 by default)
            },
            
            {    
                "type": "d4",
                "params": {"test_param": "test_param2"},
                "lcom_coef": 0.7 # (1.0 by default)
            },
        
        ]
    },
    
    
    ("singlefunc1", "disp_mix1"): {
        **func_dataset_ref['singlefunc1'],
        "lcom_dispersion":[
            {    
                "type": "d2",
                "params": {"test_param": "test_param"},
                "lcom_coef": 0.25 # (1.0 by default)
            },
            
            {    
                "type": "d3",
                "params": {"test_param": "test_param"},
                "lcom_coef": 0.25 # (1.0 by default)
            },
            
            {    
                "type": "d4",
                "params": {"test_param": "test_param"},
                "lcom_coef": 0.5 # (1.0 by default)
            },
        
        ]
    },
    
    ("singlefunc1", "disp_mix2"): {
        **func_dataset_ref['singlefunc1'],
        "lcom_dispersion":[
            {    
                "type": "d2",
                "params": {"test_param": "test_param"},
                "lcom_coef": 0.1 # (1.0 by default)
            },
            
            {    
                "type": "d3",
                "params": {"test_param": "test_param"},
                "lcom_coef": 0.2 # (1.0 by default)
            },
            
            {    
                "type": "d4",
                "params": {"test_param": "test_param"},
                "lcom_coef": 0.8 # (1.0 by default)
            },
        
        ]
    },   
}