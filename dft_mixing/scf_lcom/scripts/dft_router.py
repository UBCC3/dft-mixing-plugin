import psi4
from psi4.driver.procrouting.dft.dft_builder import (
    functionals, _dispersion_aliases, build_superfunctional_from_dictionary
)

from psi4.driver.procrouting.dft import build_superfunctional
from .lcom_dispersion import EmpiricalDispersion, LcomDispersion
from psi4.driver.p4util.exceptions import ValidationError
from qcengine.programs.empirical_dispersion_resources import dashcoeff, get_dispersion_aliases, new_d4_api
from psi4.driver.procrouting.proc import _set_external_potentials_to_wavefunction

from .lcom_superfunctional import LCOMSuperFunctionalBuilder
# libxc for XC decomposition
import pylibxc as xc
import re



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
                
def build_lcom_from_dict(func_dict: dict, npoints, deriv, restricted):
    decomp_xc : bool = func_dict.pop('decompose_xc', False)
    return LCOMSuperFunctionalBuilder(func_dict,
                                      npoints,
                                      deriv, restricted, decomp_xc).get_psi_func_dispersion()        

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
    
    _disp_functor = LcomDispersion([])
    
    # Handle dispersion
    for disp in disp_type:
        
        coef = disp.pop('lcom_coef', 1.0)
        if isinstance(name, dict):
            emp_disp = EmpiricalDispersion(name_hint='',
                                            level_hint=disp["type"],
                                            param_tweaks=disp["params"],
                                            save_pairwise_disp=save_pairwise_disp,
                                            engine=kwargs.get('engine', None))
        else:
            emp_disp = EmpiricalDispersion(name_hint=superfunc.name(),
                                            level_hint=disp["type"],
                                            param_tweaks=modified_disp_params,
                                            save_pairwise_disp=save_pairwise_disp,
                                            engine=kwargs.get('engine', None))
            
        _disp_functor.add_empirical_disp(emp_disp, coef)

    _disp_functor.print_out()
    return superfunc, _disp_functor    

# PATCHED CODE TO CREATE A MODIFIED SCF WAVEFUNCTION
# Replaces psi4.driver.procrouting.proc.scf_wavefunction_factory
def lcom_scf_wavefunction_factory(name, ref_wfn, reference, **kwargs):
    """Builds the correct (R/U/RO/CU HF/KS) wavefunction from the
    provided information, sets relevant auxiliary basis sets on it,
    and prepares any empirical dispersion.

    """
    # Figure out functional and dispersion
    superfunc, _disp_functor = lcom_build_functional_and_disp(name, restricted=(reference in ["RKS", "RHF"]), **kwargs)

    # Build the wavefunction
    psi4.core.prepare_options_for_module("SCF")
    if reference in ["RHF", "RKS"]:
        wfn = psi4.core.RHF(ref_wfn, superfunc)
    elif reference == "ROHF":
        wfn = psi4.core.ROHF(ref_wfn, superfunc)
    elif reference in ["UHF", "UKS"]:
        wfn = psi4.core.UHF(ref_wfn, superfunc)
    elif reference == "CUHF":
        wfn = psi4.core.CUHF(ref_wfn, superfunc)
    else:
        raise ValidationError("SCF: Unknown reference (%s) when building the Wavefunction." % reference)

    if _disp_functor:
        wfn._disp_functor = _disp_functor.removeNlDisp()

    # Set the DF basis sets
    df_needed = psi4.core.get_global_option("SCF_TYPE") in ["DF", "MEM_DF", "DISK_DF" ]
    df_needed |= "DFDIRJ" in psi4.core.get_global_option("SCF_TYPE")
    df_needed |= (psi4.core.get_global_option("SCF_TYPE") == "DIRECT" and core.get_option("SCF", "DF_SCF_GUESS"))
    if df_needed:
        aux_basis = psi4.core.BasisSet.build(wfn.molecule(), "DF_BASIS_SCF",
                                        psi4.core.get_option("SCF", "DF_BASIS_SCF"),
                                        "JKFIT", psi4.core.get_global_option('BASIS'),
                                        puream=wfn.basisset().has_puream())
        wfn.set_basisset("DF_BASIS_SCF", aux_basis)
    else:
        wfn.set_basisset("DF_BASIS_SCF", psi4.core.BasisSet.zero_ao_basis_set())

    # Set the relativistic basis sets
    if psi4.core.get_global_option("RELATIVISTIC") in ["X2C", "DKH"]:
        decon_basis = psi4.core.BasisSet.build(wfn.molecule(), "BASIS_RELATIVISTIC",
                                        psi4.core.get_option("SCF", "BASIS_RELATIVISTIC"),
                                        "DECON", psi4.core.get_global_option('BASIS'),
                                        puream=wfn.basisset().has_puream())
        wfn.set_basisset("BASIS_RELATIVISTIC", decon_basis)

    # Set the multitude of SAD basis sets
    if (psi4.core.get_option("SCF", "GUESS") in ["SAD", "SADNO", "HUCKEL", "MODHUCKEL"]):
        sad_basis_list = psi4.core.BasisSet.build(wfn.molecule(), "ORBITAL",
                                             psi4.core.get_global_option("BASIS"),
                                             puream=wfn.basisset().has_puream(),
                                             return_atomlist=True)
        wfn.set_sad_basissets(sad_basis_list)

        if ("DF" in psi4.core.get_option("SCF", "SAD_SCF_TYPE")):
            # We need to force this to spherical regardless of any user or other demands.
            optstash = psi4.p4util.OptionsState(['PUREAM'])
            psi4.core.set_global_option('PUREAM', True)
            sad_fitting_list = psi4.core.BasisSet.build(wfn.molecule(), "DF_BASIS_SAD",
                                                   psi4.core.get_option("SCF", "DF_BASIS_SAD"),
                                                   puream=True,
                                                   return_atomlist=True)
            wfn.set_sad_fitting_basissets(sad_fitting_list)
            optstash.restore()

    if psi4.core.get_option("SCF", "GUESS") == "SAPGAU":
        # Populate sapgau basis
        sapgau = psi4.core.BasisSet.build(wfn.molecule(), "SAPGAU_BASIS", core.get_global_option("SAPGAU_BASIS"))
        wfn.set_basisset("SAPGAU", sapgau)

    if hasattr(psi4.core, "EXTERN") and 'external_potentials' in kwargs:
        psi4.core.print_out("\n  Warning! Both an external potential EXTERN object and the external_potential" +
                       " keyword argument are specified. The external_potentials keyword argument will be ignored.\n")
        raise ValidationError("double extern")

    ep = kwargs.get("external_potentials", None)
    if ep is not None:
        _set_external_potentials_to_wavefunction(ep, wfn)

    return wfn

    

