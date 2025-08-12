
# "Bogus" parent functionals will have _D suffixes instead of -D to 
# differentiate from actual dispersions
missing_parent_functionals = [
    {
        "name": "b97_d",
        "xc_functionals": {
            "GGA_XC_B97_D": {}
        }
    }, 
    
    {
        "name": "opwlyp_d",
        "xc_functionals": {
            "GGA_XC_OPWLYP_D": {}
        }
    }, 
    
    {
        "name": "oblyp_d",
        "xc_functionals": {
            "GGA_XC_OBLYP_D": {}
        }
    }, 
    
    {
        "name": "opbe_d",
        "xc_functionals": {
            "GGA_XC_OPBE_D": {}
        }
    }, 
    
    {
        "name": "otpss_d",
        "xc_functionals": {
            "GGA_XC_OPBE_D": {}
        }
    }, 
    
    {
        "name": "DSD-SVWN",
        "x_functionals": {
            "LDA_X": {
                "alpha": 0.29
            }
        },
        "x_hf": {
            "alpha": 0.71
        },
        "c_functionals": {
            "LDA_C_VWN": {
                "alpha": 0.34
            }
        },
        "c_mp2": {
            "os": 0.58,
            "ss": 0.11
        },
        "citation": '    S. Kozuch, J.M.L. Martin, J. Comp. Chem., 34, 2327-2344, 2013\n',
        "description": '    DSD-SVWN5 SCS Double Hybrid XC Functional\n',
    },
    
    {
        "name": "DSD-BP86",
        "x_functionals": {
            "GGA_X_B88": {
                "alpha": 0.33
            }
        },
        "x_hf": {
            "alpha": 0.67
        },
        "c_functionals": {
            "GGA_C_P86": {
                "alpha": 0.49
            }
        },
        "c_mp2": {
            "os": 0.49,
            "ss": 0.24
        },
        "citation": '    S. Kozuch, J.M.L. Martin, J. Comp. Chem., 34, 2327-2344, 2013\n',
        "description": '    DSD-BP86 SCS Double Hybrid XC Functional\n'
    }
]

additional_aliases = {
    # Alias -> canonical name
    "b97_1": "b97-1",
    "b97_2": "b97-2",
    "b98": "sb98-2c", 
    "b97d": "b97_d",
    
    
    
    
}


def insert_missing_functionals():

















