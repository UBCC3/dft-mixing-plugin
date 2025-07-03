'''
    Sample data for testing the database functionality.
    The samples consist of two sets, the raw data that is going 
    to be inserted into the databse, and also a set of queries and 
    correct answer to those queries
'''

lcom_func_dataset = {   
    "multifunc1": {
        "functionals": {
            "BLYP": 1.0
        }   
    },  
    
    "multifunc2": {
        "functionals": {
            "BLYP": 0.25,
            "B3LYP": 0.75
        },
        "citation": "ABCD",
        "description": "ABCD",
    },

    "multifunc3": {
        "functionals": {
            "BP86": 0.25,
            "BPBE": 0.75,
        }
    },    
    
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

dispersion_dataset = {   
    "multifunc1": {
        "functionals": {
            "BLYP": 1.0
        }   
    },  
    
    "multifunc2": {
        "functionals": {
            "BLYP": 0.25,
            "B3LYP": 0.75
        },
        "citation": "ABCD",
        "description": "ABCD",
    },

    "multifunc3": {
        "functionals": {
            "BP86": 0.25,
            "BPBE": 0.75,
        }
    },    
    
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


































