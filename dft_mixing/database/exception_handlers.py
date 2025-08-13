
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
        "name": "mp2",
        "x_hf": {
            "alpha": 1.0
        },
        "c_functionals": {},
        "c_mp2": {"alpha": 1.0}, 
        'citation': '    Rezac, J.; Greenwell, C.; Beran, G. (2018), J. Chem. Theory Comput., 14: 4711-4721\n'
    },
    
    {
        "name": "b97_3c",
        "xc_functionals": {
            "GGA_XC_B97_3C": {}
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

dashcoeff_mapping = {
    "b97-3c": ("b97_3c", "d3bjatm"),
    "b973c": ("b97_3c", "d3bjatm"),
    
    # PSI4, Names whose canonical names cannot be broken to dash_coeffs
    # For this, only need the name of the canonical parent functional name.
    # But also need to add the alternate dispersion name as well.
    "r2scan3c": ("r2scan", "d4bjeeqatm"), 
    "r2scan-3c":  ("r2scan", "d4bjeeqatm"),
    "hf+d": ("hf", "das2010"),
    "hf-d": ("hf", "das2010"),
    "hf3c": ("hf", "d3bj2b"),
    "hf-3c": ("hf", "d3bj2b"),
    "pbeh3c": ("pbe0", "d3bj2b"),
    "pbeh-3c": ("pbe0", "d3bj2b"),
    "dldf+d09": ("dldf", "das2009"),
    "dldf-d09": ("dldf", "das2009"),
    "dldf+d10": ("dldf", "das2010"),
    "dldf-d10": ("dldf", "das2010"),
    "mp2d": ("mp2", "dmp2")
}
