import psi4
import pylibxc

# Constants from pylibxc (as defined in the source)
# (https://gitlab.com/libxc/libxc/-/blob/devel/src/xc.h#L37-40)

XC_EXCHANGE                 = 0
XC_CORRELATION              = 1
XC_EXCHANGE_CORRELATION     = 2  
XC_KINETIC                  = 3  

from psi4.driver.procrouting.dft.dft_builder import (
    functionals, _dispersion_aliases, build_superfunctional_from_dictionary
)
legacy_build_superfunc_from_dict = build_superfunctional_from_dictionary

class LCOMSuperFunctionalBuilder:
    def __init__(self, func_dict: dict, npoints, deriv, restricted):
        self._npoints = npoints
        self._deriv = deriv
        self._restricted = restricted

        self._master_sup = psi4.core.SuperFunctional.blank()
        self._disps = []
        
        # Maps (parent_func, xc_func) -> tweak
        self._tweak_mapper : dict[tuple[str, str], dict] = {}
        self._xc_func_map: dict[tuple[str, float, frozenset], psi4.core.LibXCFunctional] = {}
        
        self.build(func_dict)

    def get_psi_func_dispersion(self) -> tuple[psi4.core.SuperFunctional, list[dict]]:
        return self._master_sup, self._disps
        
    def _build(self, func_dict: dict) -> None:
        
        # Not really sure why PSI4 doesn't store the functionals as a map instead...
        # (why use a list)
    
        # Do a legacy build if the function does not have lcom functional.
        if "lcom_functionals" not in func_dict:
            disp_list = []
            sup, disp = legacy_build_superfunc_from_dict(
                func_dict, self._npoints, self._deriv, self._restricted
            )
            if disp:
                disp_list.append(disp)
                
            self._master_sup = sup
            self._disps = disp_list
            return

        self._build_lcom_helper(func_dict, 1.0)

        self._master_sup.set_c_alpha(1.0)
        self._master_sup.set_max_points(self._npoints)
        self._master_sup.set_deriv(self._deriv)
        self._master_sup.allocate()
        self._master_sup.set_name(func_dict["name"].upper())
        self._master_sup.set_xclib_description(psi4.core.LibXCFunctional.xclib_description())

        dispersions = []
        # Now handle dispersions (only read top level), 
        # Assume the format of a list of dicts with a coef field inside them.
        if "lcom_dispersion" in func_dict:
            dispersions = func_dict["lcom_dispersion"]
            
            # Separate coefficients from dictionary
            for dispersion in dispersions:
                if "citation" not in dispersion:
                    dispersion["citation"] = False

                if "nlc" in dispersion:
                    self._master_sup.set_vv10_b(-1.0)
                    self._master_sup.set_do_vv10(dispersion["nlc"])

                if dispersion["type"] == 'nl':
                    self._master_sup.set_vv10_b(dispersion["params"]["b"])
                    self._master_sup.set_vv10_c(dispersion["params"]["c"])

            # Now validate vv10 and kill if there's more from vv10
            n_disp_vv10 = sum("nlc" in d or d["type"] == 'nl' for d in dispersions)
            if n_disp_vv10 > 1:
                raise RuntimeError('Multiple dispersion coefficients with VV10. Cannot combine VV10')

        # set dispersions
        self._disps = dispersions

    def _build_lcom_helper(self, func_dict: dict, sup_coef: float):
        if "lcom_functionals" not in func_dict:
            leaf_sup, _ = legacy_build_superfunc_from_dict(func_dict, self._npoints, self._deriv, self._restricted)
            self._merge_superfunctionals(leaf_sup, sup_coef)
            
            master_xc_funcs = {}
            master_xc_funcs.update(func_dict.get("x_functionals"), {})
            master_xc_funcs.update(func_dict.get("c_functionals"), {})
            master_xc_funcs.update(func_dict.get("xc_functionals"), {})
            
            # Get tweakers
            for func_name, args in master_xc_funcs.items():
                if "tweak" in args:
                    tweak = args['tweak']
                    self._tweak_mapper[(leaf_sup.name(), func_name)] = tweak
                    
            return

        # Special case if dict has lcom_functionals (overrides other functions)
        for name, func_params in func_dict["lcom_functionals"].items():
            if isinstance(func_params, (int, float)):
            # Must also be able to find the functional in the functional list 
            # (should have been checked) before
                lcom_func = functionals[name.lower()]
                lcom_coef = func_params
            else:
                lcom_coef = func_params["coef"]
                lcom_func = func_params
                lcom_func["name"] = name

            child_coef = lcom_coef * sup_coef
            self._build_lcom_helper(lcom_func, child_coef)

    def _merge_superfunctionals(self, child: psi4.core.SuperFunctional, coef: float):
        parent = self._master_sup

        if child.needs_vv10():
            raise RuntimeError("SCF: Sorry, I don't know how to combine VV10 yet!")
        
        # GRAC fnctls are not added by the build func; They shouldn't be in the SF
        if child.needs_grac():
            raise RuntimeError("SCF: Somehow, the provided functional has GRAC. This shouldn't happen!")
        
        # Only one MP2 parameter set; if the LR is different, we're a bit screwed
        if parent.is_c_hybrid() and child.is_c_hybrid() and parent.c_omega() != child.c_omega():
            raise RuntimeError("SCF: Sorry, I don't know how to combine MP2 with different omega yet!")
        
        # For this version, just stop early if we encounter a range separated hybrid  
        if child.x_omega() != 0:
            raise RuntimeError("SCF: Don't really know how to handle multi range dispersion at this point!")

        # Merge x_hf coeffcients, total x_hf_alpha = linear combination sum of children
        parent.set_x_alpha(parent.x_alpha() + coef * child.x_alpha())
        parent.set_x_beta(parent.x_beta() + coef * child.x_beta())

        # Note: to combine css and cos, we need to sum the multiplication of alpha and c_os_alpha
        # For modifying alpha, it will come in the end instead (set it to 1.0)
        # Only modify if child is not scs.
        if coef * child.c_alpha() != 0:
            parent.set_c_os_alpha(parent.c_os_alpha() + coef * child.c_alpha() * child.c_os_alpha())
            parent.set_c_ss_alpha(parent.c_ss_alpha() + coef * child.c_alpha() * child.c_ss_alpha())
            parent.set_c_alpha(1.0)

        x_funcs = []
        c_funcs = []
        
        # Special case where child is a xc_func, decompose XC functionals.
        if child.is_libxc_func():
            
            xc_psi_func = child_fnctl.c_functionals()[0]
            xc_name = xc_psi_func.name()
            libxc_func = pylibxc.LibXCFunctional(xc_name, self._restricted)
        
            xc_coef = xc_psi_func.alpha()
            xc_decomp : list[tuple[float, str]] = libxc_func.aux_funcs()
            
            if xc_decomp is None:
                c_funcs.append(xc_psi_func)
                
            else:    
                for internal_coef, func_name in xc_decomp:
                    psi_fname = ("XC_" + func_name).upper()
                    internal_libxc_func = pylibxc.LibXCFunctional(func_name, self._restricted)
                    psi_internal_func = psi4.core.LibXCFunctional(psi_fname, self._restricted)
                    
                    combined_coef = xc_coef * internal_coef
                    psi_internal_func.set_alpha(psi_internal_func.alpha() * combined_coef)

                    internal_func_kind = internal_libxc_func.get_kind()
                    
                    # Check the type of fnctl              
                    if internal_func_kind == XC_EXCHANGE:
                        x_funcs.append(psi_internal_func)
                    elif internal_func_kind == XC_CORRELATION:
                        c_funcs.append(psi_internal_func)
                    
                    elif (internal_func_kind == XC_EXCHANGE_CORRELATION
                            or internal_libxc_func == XC_KINETIC):
                        raise RuntimeError(f"ERROR: Cannot decompose functional {func_name}")    

        else:
            x_funcs.extend(child.x_functionals())
            c_funcs.extend(child.c_functionals())


        func_handlers = [
            (x_funcs, parent.add_x_functional),
            (c_funcs, parent.add_c_functional),
        ]

        for child_fnctls, parent_add_func in func_handlers:
            for child_fnctl in child_fnctls:
                child_fnctl.set_alpha(coef * child_fnctl.alpha())

                # Check for tweaks: 
                tweak = self._tweak_mapper.get((child.name(), child_fnctl.name()), {})

                key = (child_fnctl.name(), child_fnctl.omega(), frozenset(tweak))
                if key in self._xc_func_map:
                    target_fnctl = self._xc_func_map[key]
                    target_fnctl.set_alpha(target_fnctl.alpha() + child_fnctl.alpha())
                else:
                    self._xc_func_map[key] = child_fnctl
                    parent_add_func(child_fnctl)

    