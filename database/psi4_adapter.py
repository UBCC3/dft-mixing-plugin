
from .database_v2 import FunctionalDatabase
from .db_models import (
    Sources,
    Functional,
    FunctionalAlias,
    FunctionalCoeffs,
    DispersionBase,
    DispersionConfig,
    DispersionAlias
)

import logging
from typing import Any

logger = logging.getLogger(__name__)

class Psi4DbAdapter:
    '''
        PSI4 database w
    '''

    def __init__(self, config_path):
        self.db = FunctionalDatabase(config_path)
    
    
    def _import_psi4_data(self):
        pass
    
    def load_base_functional_data(self):
        pass
    
    def load_multi_functional_data(self):
        pass
    
    def load_base_dispersion_data(self):
        pass
    
    def load_dispersion_config_data(self):
        pass
        
    def get_functional_dict(self, 
                              functional_name : str,
                              dispersion_name : str|None = None,
                              functional_source : str|None = None,
                              dispersion_source : str|None = None):
        '''
            Queries for a functional and an optional dispersion
            configuration for that functional, and formats it
            to a PSI4 compatible dictionary. 
            
            Important Note:
                For dispersions, the adapater queries for the 
                dashcoeff pair "{functional_name}-{dispersion_name}".
        ''' 
        
        dashcoeff_name = f"{functional_name}-{dispersion_name}"

        # Query for functional and dispersion components
        with self.db.get_session() as session:
            func, func_coeffs = self.db.get_functional(
                session, 
                functional_name,
                functional_source,
            )
            
            disp_coeffs = self.db.get_multi_dispersion(
                session, dashcoeff_name, functional_source,
                dispersion_source
            )
            session.commit()
    
        return self._format_to_psi4(func, func_coeffs, disp_coeffs)
    
    def _format_to_psi4(self, 
                        par_functional : Functional,
                        functional_coeffs: list[tuple[Functional, float]],
                        dispersion_coeffs: list[tuple[DispersionBase, float]] = []):
        '''
            Formats query results to PSI4 lcom plugin format.
        '''    
        
        func_dict : dict[str, Any] = {}
        func_dict["name"] = par_functional.fnctl_name
        func_dict["citation"] = par_functional.citation
        func_dict["description"] = par_functional.description
        
        if len(functional_coeffs) > 1:
            lcom_functionals : dict[str, Any] = {}
            for (child_fnctl, coef) in functional_coeffs:
                child_name = str(child_fnctl.fnctl_name)
                child_data = child_fnctl.fnctl_data.copy()
                lcom_functionals[child_name] = child_data
                lcom_functionals[child_name]["coef"] = coef    
            
            func_dict["lcom_functionals"] = lcom_functionals
        
        else:
            func_dict = par_functional.fnctl_data        
            
            
        logger.warning(f"Dipsersion Coeffs: {dispersion_coeffs}") 
            
        if len(dispersion_coeffs) > 0:
            lcom_disps = []
            for (base_disp, coef) in dispersion_coeffs:
                disp_dict = {}
                disp_dict["params"] = base_disp.disp_params
                disp_dict["type"] = base_disp.subdisp_name
                disp_dict["citation"] = base_disp.citation
                disp_dict["description"] = base_disp.description
                disp_dict["lcom_coef"] = coef
                lcom_disps.append(disp_dict)
            
            if len(dispersion_coeffs) == 1:
                func_dict["dispersion"] = lcom_disps[0]
            else:
                func_dict["lcom_dispersion"] = lcom_disps
        
        if func_dict is None:
            logger.error(f"FUNC DICT FORMATTING IS NONE {func_dict}")
        
        return func_dict        
        

    def get_functional_dict(self, 
                              functional_name : str,
                              dispersion_name : str|None = None,
                              functional_source : str|None = None,
                              dispersion_source : str|None = None):
        '''
            Queries for a functional and an optional dispersion
            configuration for that functional, and formats it
            to a PSI4 compatible dictionary. 
        ''' 
        pass   



