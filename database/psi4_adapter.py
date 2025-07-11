
from database import FunctionalDatabase
from db_models import (
    Sources,
    Functional,
    FunctionalAlias,
    FunctionalCoeffs,
    DispersionBase,
    DispersionConfig,
    DispersionAlias
)

from db_errors import (
    DBDuplicateError,
    DBNotFoundError
)

import logging
from typing import Any

import re

# PSI4 internal functionals
from psi4.driver.procrouting.dft.dft_builder import functionals

# DFTD4 internal functionals
logger = logging.getLogger(__name__)

class Psi4DbAdapter:
    '''
        PSI4 database wrapper for a Database handler
    '''

    def __init__(self, config_path):
        self.db = FunctionalDatabase(config_path)
        self._import_psi4_data()
    
    def _import_psi4_data(self):
        
        parent_functionals = {
            func_name : functionals[func_name] \
                for func_name in functionals \
                if "dispersion" not in functionals[func_name]
        }
        
        dispersion_functionals = {
            func_name : functionals[func_name] \
                for func_name in functionals \
                if "dispersion" in functionals[func_name]
        }
        
        
        for func_alias, func_dict in parent_functionals.items():
            canon_fname = func_dict["name"]
            func_citation = func_dict.get("citation", "No citations available")
            func_description = func_dict.get("description", "")
            func_data = func_dict
            
            #   For any insertion phase, we need to check whether 
            #   the functional already exists or not, 
            with self.db.get_session() as session:
                try:
                    self.db.get_single_functional(session, canon_fname, "psi4")
                
                # Have not exist, actually add functional + canon alias
                except DBNotFoundError as e:
                    session.rollback()
                    
                    self.db.add_functional_alias(session, canon_fname, canon_fname)
                    self.db.insert_base_functional(
                        session, 
                        canon_fname,
                        func_citation,
                        func_description,
                        func_data,
                        "psi4",
                    )
                    
                # Now add the alias (if its not the same as canon fname)
                finally:
                    if (canon_fname.lower() != func_alias.lower()):
                        self.db.add_functional_alias(session, canon_fname, func_alias)
                    session.commit()
                
            # Just to make sure functional is in database.
            with self.db.get_session() as session:
                # Insert base functional into PSI4
                # logger.warning(f"TEST QUERY {func_alias} -> {canon_fname} into db")
                self.db.get_functional(session, func_alias, "psi4") 
                session.commit()
                
            
        # Insert dispersion
        for func_dashcoeff, func_dict in dispersion_functionals.items():
            canon_fdashcoeff = func_dict["name"]
            pattern = r"([-\d\w()]+)-([\w\d()]+)\)?$"
            match = re.match(pattern, canon_fdashcoeff)
            
            if not match: 
                logger.error(f"Dispersion causing error, Skipping functional {func_dashcoeff}:{func_dict} ")
                continue
                
            parent_func, _ = match.groups()
            disp_dict = func_dict["dispersion"]
            canon_dname = disp_dict["type"]
            
            # Check if parent functional exists or not
            with self.db.get_session() as session:
                try:                   
                    func_id = self.db.get_single_functional(session,
                        functional_name= parent_func, source="psi4")
                except DBNotFoundError as e:
                    logger.error(f"Dispersion {func_dashcoeff} does not have parent functional")
                    continue                
        
            with self.db.get_session() as session:
                try:
                    # Check if CANONICAL dispersion exists first
                    self.db.get_base_dispersion(session, f'{parent_func}-{canon_dname}', "psi4", "psi4")

                except DBNotFoundError as e:
                    session.rollback()
                    
                    # Insert canonical alias
                    self.db.add_dash_coeff_mapping(session, parent_func, canon_dname, f'{parent_func}-{canon_dname}')
                    self.db.insert_single_disp(
                        session, 
                        parent_func,
                        canon_dname,
                        disp_dict.get("citation", "No citation."),
                        disp_dict.get("description", "No description."),
                        disp_dict["params"],
                        "psi4",
                        "psi4",
                    )
                    
                finally:
                    if (f'{parent_func}-{canon_dname}'.lower() != func_dashcoeff.lower()):
                        self.db.add_dash_coeff_mapping(session, parent_func, canon_dname, func_dashcoeff)
                    session.commit()
                    
            
            with self.db.get_session() as session:
                # Just check if its actually there...
                self.db.get_base_dispersion(session, func_dashcoeff, "psi4", "psi4") 
                session.commit()
    
    
    def load_base_functional_data(self, func_dict):
        raise NotImplementedError("Not implemented yet!")
    
    def load_multi_functional_data(self, multi_func_dict, source: str):
        '''
            Loads in multi-functional data.
        '''
        error_list = []
        
        for multifunc_alias, multi_args in multi_func_dict.items():
            try:
                components = multi_args["functionals"]
                canon_name = multi_args.get("name", multifunc_alias)
                
                with self.db.get_session() as session:                
                    self.db.insert_multi_functional(
                        session,
                        canon_name,
                        multi_args.get("citation", "No citation."),
                        multi_args.get("description", "No description."),
                        components,
                        source
                    )
                    
                    session.commit()
            except Exception as e:
                error_list.append(e)
                
        if len(error_list) > 0:
            raise ExceptionGroup("Some errors were encountered:", error_list)            
        return
        
    
    def load_base_dispersion_data(self, base_disp_dict, src: str):
        error_list = []
        
        for parent_func, avail_disps in base_disp_dict.items():
            for disp_alias, disp_info in avail_disps.items():
                try:
                    
                    dash_coeff_name = f'{parent_func}-{disp_alias}'
                    with self.db.get_session() as session:
                    
                        self.db.insert_single_disp(
                            session,
                            dash_coeff_name,
                            parent_func,
                            disp_alias,
                            disp_info.get("citation", "No Citation"),
                            disp_info.get("description", "No Description"),
                            disp_info["disp_params"],
                            src,                            
                        )
                    
                        session.commit()
                except Exception as e:
                    error_list.append(e)
                    
        if len(error_list) > 0:
            raise ExceptionGroup("Some errors were encountered:", error_list)            
        return
    
    def load_dispersion_config_data(self, disp_config_dict, src: str):
        error_list = []
        
        for parent_func, avail_configs in disp_config_dict.items():
            for disp_config_name, config_args in avail_configs.items():
                try:
                    dash_coeff_name = f'{parent_func}-{disp_config_name}'
                    with self.db.get_session() as session:
                        self.db.insert_disp_config(
                            session,
                            dash_coeff_name,
                            parent_func,
                            disp_config_name,
                            config_args.get("citation", "No Citation"),
                            config_args.get("description", "No Description"),
                            config_args["coeffs"],
                            src,                            
                        )
                    
                        session.commit()
                except Exception as e:
                    error_list.append(e)
                    
        if len(error_list) > 0:
            raise ExceptionGroup("Some errors were encountered:", error_list)            
        return
        
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
        
        # Query for functional first
        with self.db.get_session() as session:
            par_functional, func_coeffs = self.db.get_functional(
                session, functional_name, functional_source
            )
            
            if dispersion_name is not None:
                dashcoeff_name = f'{par_functional.fnctl_name}-{dispersion_name}'
                dispersion_coeffs = self.db.get_multi_dispersion(
                    session, dashcoeff_name, functional_source, dispersion_source
                )
            else:
                dispersion_coeffs = []
            
            session.commit()
            
        return self._format_to_psi4(par_functional, func_coeffs, dispersion_coeffs)
        