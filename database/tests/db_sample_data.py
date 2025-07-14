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
            "PWPW": 0.75,
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
            "BLYP": {
                'coef': 1.0,
                **functionals['blyp']
            }
        }
    },  
    
    "singlefunc2": {
        "name": "singlefunc2",
        "lcom_functionals": {
            "pbe0": {
                'coef': 1.0,
                **functionals['pbe0']
            }
        }   
    },  
    
    "multifunc2": {
        "name": "multifunc2",
        "lcom_functionals": {
            "BLYP": {
                'coef': 0.25,
                **functionals['blyp']
            },
            "B3LYP": {
                'coef': 0.75,
                **functionals['b3lyp']
            }
        },
    },

    "multifunc3": {
        "name": "multifunc3",
        "lcom_functionals": {
            "BP86": {
                'coef': 0.25,
                **functionals['bp86']
            },
            
            # Aliased name
            "PWPW": {
                'coef': 0.75,
                **functionals['pwpw']
            },
        }
    },    
}

lcom_psi_dispconfig_dataset : dict = {
    "TPSS": {
        "disp_mix1": {
            "disp1": 0.25,
            "disp2": 0.25,
            "disp3": 0.5,
        },
        
        "disp_mix2": {
            "d2": 0.55,
            "d3": 0.25,
            "d3bj": 0.15,
        },
    },  
    
    "singlefunc1": {
        "disp_mix1": {
            "disp2": 0.25,
            "disp3": 0.25,
            "disp4": 0.5,
        },
        
        "disp_mix2": {
            "d2": 0.1,
            "d3": 0.2,
            "d3bj": 0.8,
        }, 
    },  
    
    "multifunc2": {
        "disp_mix1": {
            "d2": 0.25,
            "d3": 0.25,
            "d4": 0.5,
        },
        
        "disp_mix2": {
            "d2": 0.1,
            "d3": 0.2,
            "d3bj": 0.8,
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
            "params": {
                "test_param": "test_param2"
            }
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