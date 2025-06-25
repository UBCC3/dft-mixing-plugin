import psi4
from psi4.driver.procrouting.dft.dft_builder import (
    functionals, _dispersion_aliases, build_superfunctional_from_dictionary
)

from psi4.driver.procrouting.dft import build_superfunctional
from .lcom_dispersion import EmpiricalDispersion, LcomDispersion
from psi4.driver.p4util.exceptions import ValidationError
from qcengine.programs.empirical_dispersion_resources import dashcoeff, get_dispersion_aliases, new_d4_api

# libxc for XC decomposition
import pylibxc as xc

# Legacy build_from_dict
legacy_build_superfunc_from_dict = build_superfunctional_from_dictionary

def _check_consistency(func_dictionary):
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
        # 3b.1) Case where dispersion is a list
        if isinstance(disp, list):
            for lcom_disp in disp:
                if ("type" not in disp or disp["type"] not in _dispersion_aliases):
                    raise ValidationError(
                        f"SCF: Dispersion type ({disp['type']}) should be among ({_dispersion_aliases.keys()})")      
        
        # 3b.2) check dispersion type present and known (if it is a dict)
        elif ("type" not in disp or disp["type"] not in _dispersion_aliases):
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
                
def _merge_superfunctionals(parent: psi4.core.SuperFunctional, child: psi4.core.SuperFunctional, 
                            xc_func_dict: dict, coef : float = 1.0, ) -> None:
    '''
        Helper function to mix a child superfunctional into a parent superfunctional.
        Modifies the parent functional (must be the master functional)
        
        Notes:
        - coef is the absolute weight with relation to the parent.
        
        Also updates the global parameters for a more consistent handling. 
    '''
    
    # EARLY STOPPING
    
    if child.needs_vv10(): raise RuntimeError("SCF: Sorry, I don't know how to combine VV10 yet!")
    # GRAC fnctls are not added by the build func; They shouldn't be in the SF
    if child.needs_grac(): raise RuntimeError("SCF: Somehow, the provided functional has GRAC. This shouldn't happen!")
    # Only one MP2 parameter set; if the LR is different, we're a bit screwed
    if parent.is_c_hybrid() and child.is_c_hybrid() and parent.c_omega() != child.c_omega():
        raise RuntimeError("SCF: Sorry, I don't know how to combine MP2 with different omega yet!")
    # For this version, just stop early if we encounter a range separated hybrid
    if child.x_omega() != 0:
        raise RuntimeError("SCF: Don't really know how to handle multi range dispersion at this point!")
    
    # @TODO: Update this to use RangeStep in a later version for range separated hybrid
    
    # Assign functional groups to functional adders.
    func_handlers = [(child.x_functionals(), parent.add_x_functional),
                     (child.c_functionals(), parent.add_c_functional)]
   
    # Combine LibXC functionals 
    for child_fnctls, parent_add_func in func_handlers:
        for child_fnctl in child_fnctls:
            child_fnctl.set_alpha(coef * child_fnctl.alpha())        
            
            # TODO: Optimize this
            # This could have been optimized :P, right now it is O(n^2) w.r.t number of functionals
            
            # parent_masterlist = parent.x_functionals() + parent.c_functionals()
            
            # Next, check if we can find a compatible functional to add to
            # (Needs to be the same fnctl (name) and have the same LR cutoff (omega))        
            # needs_addition = True
            
            # Now just simply check if the functional has an entry in
            # the parent dict. (No need to loop over everything)
            
            # !TODO: Need to implement support for tweaks, though, it is a private field
            # I don't think making this a plugin is no longer sustainable
            # probably going back to fork : (
            
            xc_func_key = (child_fnctl.name(), child_fnctl.omega(), frozenset({}))
            if xc_func_key in xc_func_dict:
                # Add the coefficient
                target_fnctl = xc_func_dict[xc_func_key]
                target_fnctl.set_alpha(target_fnctl.alpha() + child_fnctl.alpha())
                
            else: 
                xc_func_dict[xc_func_key] = child_fnctl
                parent_add_func(child_fnctl) 
                 
            # for potential_tgt in parent_masterlist:
            #     if potential_tgt.name() != child_fnctl.name(): continue
            #     if potential_tgt.omega() != child_fnctl.omega(): continue
            #     if potential_tgt.tweaks() != child_fnctl.tweaks(): continue
            #     needs_addition = False
            #     print(f"PAR COEF {potential_tgt.alpha()} CHILD COEF{child_fnctl.alpha()}")
            #     potential_tgt.set_alpha(potential_tgt.alpha() + child_fnctl.alpha())
            #     break
            
            # if needs_addition: 
            #     parent_add_func(child_fnctl)    
             
    # Merge x_hf coeffcients, total x_hf_alpha = linear combination sum of children
    parent.set_x_alpha(parent.x_alpha() + coef * child.x_alpha())
    parent.set_x_beta(parent.x_beta() + coef * child.x_beta())
    
    # Note: to combine css and cos, we need to sum the multiplication of alpha and c_os_alpha
    # For modifying alpha, it will come in the end instead (set it to 1.0)
    # Only modify if child is not scs.
    if (coef * child.c_alpha() != 0):
        parent.set_c_os_alpha(parent.c_os_alpha() + coef * child.c_alpha()* child.c_os_alpha())
        parent.set_c_ss_alpha(parent.c_os_alpha() + coef * child.c_alpha()* child.c_os_alpha()) 
        parent.set_c_alpha(1.0)
        
    # Merge x_hf coeffcients, total x_hf_alpha = linear combination sum of children
    parent.set_x_alpha(parent.x_alpha() + coef * child.x_alpha())
    
# def _print_out_lcom(xc_func_map):
#     '''
#         Prints out the final LCOM parameters (already combined), as well as the 
#         final HF parameters.
#     '''




# def _decompose_xc_helper():
#     '''
#         Helper function to decompose XC functionals into X and C parts.
#     '''
    
    


def build_lcom_helper(  func_dict: dict,
                        npoints, deriv, 
                        restricted,
                        xc_func_map: dict,
                        master_sup: psi4.core.SuperFunctional, 
                        sup_coef : float = 1.0, ) \
                    -> None:
    '''
        Recursive helper to build a lcom functional.
        Modifies a shared map between all build calls that maps 
        functional LIBXC name to the PSI4 Functional object.
        
        This function only **retrieves** individual LibXC PSI4 functionals,
        and does not create any new superfunctional, unless they are the 
        base PSI4 versions.
    '''
    _check_consistency(func_dict)
    print(func_dict)
    
    # Base Case: normal version, call legacy version (BOTH X, C and XC)
    if "lcom_functionals" not in func_dict:
        
        # Now do the superfunctional, ignore dispersions since it is already done in the top level          
        leaf_sup, _ = legacy_build_superfunc_from_dict(func_dict, npoints, deriv, restricted)

        # All we do now is just adding coefficients to the functional map and also add it 
        # to the superfunctional.        
        _merge_superfunctionals(master_sup, leaf_sup, xc_func_map, sup_coef)    

        return
            
    # ============= RECURSIVE PART ======================

    # Special case if dict has lcom_functionals (overrides other functions)
    # sup = psi4.core.SuperFunctional.blank()
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
        
        # "Compile functional" as a superfunctional. (is a list)
        child_coef = lcom_coef * sup_coef
        print("CHILD FUNCS",name, lcom_func)
        build_lcom_helper(lcom_func, npoints, deriv, restricted, xc_func_map, master_sup, sup_coef=child_coef)
    return


def build_lcom_functional(func_dict: dict, npoints, deriv, restricted) -> tuple[psi4.core.SuperFunctional, list[dict]]:
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
        
        Returns a tuple of superfunctional and a list of dispersion represented as a dict or string
        with dispersion name.
        If there are no dispersion, return an empty list.
    '''
    
    # Not really sure why PSI4 doesn't store the functionals as a map instead...
    # (why use a list)
    
    # Do a legacy build if the function does not have lcom functional.
    if "lcom_functionals" not in func_dict:
        disp_list = []
        sup, disp = legacy_build_superfunc_from_dict(func_dict, npoints, deriv, restricted) 
        
        if disp:
            disp_list.append(disp)
            
        return sup, disp_list
    
    
    # Now get PSI4 environment variable to determine whether to decompose XC functionals or not

    # Build the actual superfunctional
    master_sup = psi4.core.SuperFunctional.blank()
    
    # Dynamic programming to map the tuple (fnctl_name, omega, tweaks) -> lcom_coef 
    # in the final superfunctional.
    xc_func_map : dict[tuple[str, float, frozenset], float] = {}
    build_lcom_helper(func_dict, npoints, deriv, restricted, xc_func_map, master_sup, 1.0)

    # Master should be built, now set c_alpha to 1.0
    master_sup.set_c_alpha(1.0)

    # After collecting all functionals, print them out

    master_sup.set_max_points(npoints)
    master_sup.set_deriv(deriv)
    master_sup.allocate()
    master_sup.set_name(func_dict["name"].upper())
    master_sup.set_xclib_description(psi4.core.LibXCFunctional.xclib_description())
    
    dispersions = []
    
    # Now handle dispersions (only read top level), 
    # Assume the format of a list of dicts with a coef field inside them.
    if "lcom_dispersion" in func_dict:
        dispersions : list[dict ]= func_dict["lcom_dispersion"]
        
        # Separate coefficients from dictionary
        for dispersion in dispersions:
            if "citation" not in dispersion:
                dispersion["citation"] = False
                
            if "nlc" in dispersion:
                sup.set_vv10_b(-1.0)
                sup.set_do_vv10(dispersion["nlc"])
                
            if dispersion["type"] == 'nl':
                sup.set_vv10_b(dispersion["params"]["b"])
                sup.set_vv10_c(dispersion["params"]["c"])
            
        # Now validate vv10 and kill if there's more from vv10
        n_disp_vv10 = sum(map(lambda d: "nlc" in d or d["type"] == 'nl', dispersions))
        if n_disp_vv10 > 1:
            raise RuntimeError('Multiple dispersion coefficients with VV10. Cannot combine VV10')
        
    return master_sup, dispersions        
        

# PATCHED CODE TO SUPPORT COMBINING DISPERSION
def lcom_build_functional_and_disp(name, restricted, save_pairwise_disp=False, **kwargs):

    if psi4.core.has_option_changed("SCF", "DFT_DISPERSION_PARAMETERS"):
        modified_disp_params = psi4.core.get_option("SCF", "DFT_DISPERSION_PARAMETERS")
    else:
        modified_disp_params = None

    # Check if we should decompose_xc, this is a bit hacky
    if "decompose_xc" in kwargs and kwargs["decompose_xc"]:
        name["decompose_xc"] = True

    # Figure out functional
    # Returns with a list of dispersion_type, or None at all.
    superfunc, disp_type = build_superfunctional(name, restricted)

    if len(disp_type) == 0:
        return superfunc, None
    
    _disp_functor = LcomDispersion()
    
    # Handle dispersion
    for disp in disp_type:
        
        coef = disp.pop('lcom_coef', 1.0)
        if isinstance(name, dict):
            emp_disp = EmpiricalDispersion(name_hint='',
                                            level_hint=disp_type["type"],
                                            param_tweaks=disp_type["params"],
                                            save_pairwise_disp=save_pairwise_disp,
                                            engine=kwargs.get('engine', None))
        else:
            emp_disp = EmpiricalDispersion(name_hint=superfunc.name(),
                                            level_hint=disp_type["type"],
                                            param_tweaks=modified_disp_params,
                                            save_pairwise_disp=save_pairwise_disp,
                                            engine=kwargs.get('engine', None))
            
        _disp_functor.add_empirical_disp(emp_disp, coef)

    _disp_functor.print_out()
    return superfunc, _disp_functor    
    

