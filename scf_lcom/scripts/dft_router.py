import psi4
from psi4.driver.procrouting.dft.dft_builder import (
    functionals, _dispersion_aliases, build_superfunctional_from_dictionary
)

from lcom_dispersion import EmpiricalDispersion, LcomDispersion
from psi4.driver.p4util.exceptions import ValidationError
from qcengine.programs.empirical_dispersion_resources import dashcoeff, get_dispersion_aliases, new_d4_api


def check_consistency(func_dictionary):
    """
    This checks the consistency of the definitions of exchange and correlation components
    of the functional, including detecting duplicate requests for LibXC params, inconsistent
    requests for HF exchange and missing correlation. It also makes sure that names of methods
    passed in using dft_functional={} syntax have a non-implemented name.
    """
    # 0a) make sure method name is set:
    if "name" not in func_dictionary:
        raise ValidationError("SCF: No method name was specified in functional dictionary.")
    else:
        name = func_dictionary["name"]
        # 0b) make sure provided name is unique:
        if (name.lower() in functionals.keys()) and (func_dictionary not in functionals.values()):
            raise ValidationError(
                "SCF: Provided name for a custom dft_functional matches an already defined one: %s." % (name))

    # 1a) sanity checks definition of lcom_functionals
    if "lcom_functionals" in func_dictionary:
        if "x_functionals" in func_dictionary or "x_hf" in func_dictionary:
            raise ValidationError("SCF: Duplicate specification of exchange (LCOM + X) in functional %s." % (name))
        elif "c_functionals" in func_dictionary or "c_mp2" in func_dictionary:
            raise ValidationError("SCF: Duplicate specification of correlation (LCOM + C) in functional %s." % (name))
        elif "xc_functionals" in func_dictionary:
            raise ValidationError("SCF: Duplicate specification of correlation (LCOM + XC) in functional %s." % (name))
        elif not isinstance(func_dictionary["lcom_functionals"], dict):
            raise ValidationError("SCF: Invalid lcom_functionals definition in functional %s." % (name))
        for val in func_dictionary["lcom_functionals"].values():
            if (not isinstance(val, dict)) and (not isinstance(val, float)) and (not isinstance(val, int)):
                raise ValidationError("SCF: lcom_functionals definition in functional %s has invalid entry." % (name))

    # 1a) sanity checks definition of xc_functionals
    elif "xc_functionals" in func_dictionary:
        if "x_functionals" in func_dictionary or "x_hf" in func_dictionary:
            raise ValidationError("SCF: Duplicate specification of exchange (XC + X) in functional %s." % (name))
        elif "c_functionals" in func_dictionary or "c_mp2" in func_dictionary:
            raise ValidationError("SCF: Duplicate specification of correlation (XC + C) in functional %s." % (name))

    # 1b) require at least an empty exchange functional entry or X_HF
    elif "x_functionals" not in func_dictionary and "x_hf" not in func_dictionary:
        raise ValidationError("SCF: No exchange specified in functional %s." % (name))

    # 1c) require at least an empty correlation functional entry or C_MP2
    elif "c_functionals" not in func_dictionary and "c_mp2" not in func_dictionary:
        raise ValidationError("SCF: No correlation specified in functional %s." % (name))

    # 2) use_libxc handling:
    use_libxc = 0
    if "x_functionals" in func_dictionary:
        for item in func_dictionary["x_functionals"]:
            if "use_libxc" in func_dictionary["x_functionals"][item] and \
            func_dictionary["x_functionals"][item]["use_libxc"]:
                use_libxc += 1

    # 2a) only 1 component in x_functionals can have "use_libxc": True to prevent libxc conflicts
    if use_libxc > 1:
        raise ValidationError("SCF: Duplicate request for libxc exchange parameters in functional %s." % (name))

    # 2b) if "use_libxc" is defined in x_functionals, there shouldn't be an "x_hf" key
    elif use_libxc == 1 and "x_hf" in func_dictionary:
        raise ValidationError("SCF: Inconsistent definition of exchange in functional %s." % (name))

    # 2c) ensure libxc params requested in "x_hf" are for a functional that is included in "x_functionals"
    elif "x_hf" in func_dictionary and "use_libxc" in func_dictionary["x_hf"] \
    and func_dictionary["x_hf"]["use_libxc"] not in func_dictionary["x_functionals"]:
        raise ValidationError(
            "SCF: Libxc parameters requested for an exchange functional not defined as a component of %s." % (name))

    # 3) checks would be caught at runtime or involve only formatting.
    #    included here to preempt driver definition problems, if specific fctl not in tests.
    # 3a) check formatting for citation
    if "citation" in func_dictionary:
        cit = func_dictionary["citation"]
        if cit and not (cit.startswith('    ') and cit.endswith('\n')):
            raise ValidationError(
                f"SCF: All citations should have the form '    A. Student, B. Prof, J. Goodstuff Vol, Page, Year\n', not : {cit}"
            )
    if "dispersion" in func_dictionary:
        disp = func_dictionary["dispersion"]
        # 3b) check dispersion type present and known
        if "type" not in disp or disp["type"] not in _dispersion_aliases:
            raise ValidationError(
                f"SCF: Dispersion type ({disp['type']}) should be among ({_dispersion_aliases.keys()})")
    # 3c) check dispersion params complete
        allowed_params = sorted(dashcoeff[_dispersion_aliases[disp["type"]]]["default"].keys())
        if "params" not in disp or sorted(disp["params"].keys()) != allowed_params:
            raise ValidationError(
                f"SCF: Dispersion params for {name} ({list(disp['params'].keys())}) must include all ({allowed_params})")
    # 3d) check formatting for dispersion citation
        if "citation" in disp:
            cit = disp["citation"]
            if cit and not ((cit.startswith('    ') and cit.endswith('\n')) or re.match(r"^10.\d{4,9}/[-._;()/:A-Z0-9]+$", cit)):
                raise ValidationError(
                    f"SCF: All citations should have the form '    A. Student, B. Prof, J. Goodstuff Vol, Page, Year\n', not : {cit}"
                )
                
def merge_superfunctionals(parent: psi4.core.SuperFunctional, child: psi4.core.SuperFunctional, 
                            coef : float = 1.0) -> None:
    '''
        Helper function to mix a child superfunctional into a parent superfunctional.
        Modifies the parent functional.
    '''
    
    # Assign functional groups to functional adders.
    func_handlers = [(child.x_functionals, parent.add_x_functional),
                     (child.c_functionals, parent.add_c_functional)]
    
    # Merge functionals
    for child_fnctls, parent_add_func in func_handlers:
        for child_fnctl in child_fnctl:
            child_fnctl.set_alpha(coef * child_fnctl.alpha())    
            parent_add_func(child_fnctl)
            
    # After merging parent, need to merge global exchange (HF) params
    # For now, only support alpha merging and cannot merge diffeerent omegas

    
    # Merge correlation mp2 params
    
    
def build_lcom_functional(func_dict, restricted) -> tuple[psi4.core.SuperFunctional, list[dict]]:
    '''
        Builds a SuperFunctional object (compatible with current version of PSI4)  
        that mixes multiple DFT functionals as a linear combination.
        
        Format of a lcom_functional (RECCOMENDED) is:
        {
            "lcom_functional":
                {
                    <fnctl_name>: {
                        "coef": <leading coefficient for functional>
                        ... 
                        <other superfunctional fields> 
                    }
                    
                    ...
                    
                }
        }     
        
        You can also map the functional name to a number (as the lcom_coefficient)
        for backwards compatibility.
    '''
    check_consistency(func_dict)
    dispersion = []

    # Special case if dict has lcom_functionals (overrides other functions)
    if "lcom_functionals" in func_dict:
        sup = psi4.core.SuperFunctional.blank()
        lcom_funcs = func_dict["lcom_functionals"]
        
        for (name, func_params) in lcom_funcs.items():
            if (isinstance(func_params, (int, float))):
                # Must also be able to find the functional in the functional list 
                # (should have been checked) before
                lcom_func = functionals[name.lower()]
                lcom_coef = func_params
            
            else:
                lcom_coef = func_params["coef"]
                lcom_func = func_params
                lcom_func["name"] = name
            
            # "Compile functional" as a superfunctional. 
            child_sup, disp = build_lcom_functional(func_dict, restricted)
            
            # Merge child superfunctional into parent's        
            merge_superfunctionals(sup, child_sup, lcom_coef)
            
    
    # Otherwise, direct this to the normal PSI4 call
    return build_superfunctional_from_dictionary(func_dict, )


    
    









