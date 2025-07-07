import logging.config
import os
from pathlib import Path

import sqlalchemy.orm
from sqlalchemy import Column, String, Double, ForeignKey, UUID, JSON, UniqueConstraint
from sqlalchemy.orm import declarative_base, relationship
from sqlalchemy.ext.hybrid import hybrid_property

# For functionals
from psi4.driver.procrouting.dft.dft_builder import functionals

# Default dispersions
from qcengine.programs.empirical_dispersion_resources import dashcoeff, get_dispersion_aliases, new_d4_api

# blacklisted functionals
from .exceptions import blacklist_funcs

# Error types
from db_errors import (
    SourceResolutionError,
    DBNotFoundError,
    DBDuplicateError
)

from .db_models import (
    get_base,
    Functional,
    FunctionalAlias,
    FunctionalCoeffs,
    DispersionAlias,
    DispersionBase, 
    DispersionConfig, 
    Sources
)


import uuid
import re
import warnings
import json

import logging

# Config parser
import yaml

logger = logging.getLogger(__name__) 


class FunctionalDatabase:
    
    def __init__(self, database_config : Path):
        '''
            Creates a functional database from scratch, or reads from
            an existing database configuration (in a yaml file.)
            
            Yaml contains:
            1. Path to database
            
        '''     
        self.source_fallback_stack : list[str] = ["dftd4", "psi4"]
        
        # Import from existing database if path is specified 
        if database_config is None:
            raise ValueError("Error: database configurations cannot be None")
      
        # Import configs
        with open(database_config) as f:
            config = yaml.safe_load(f)

        DB_PATH = config["db_path"]
        
        # Check if db exists
        db_exists = os.path.exists(DB_PATH)
        
        db_uri = f'sqlite:///{DB_PATH}'
        # Create a SQL database for this
        engine = sqlalchemy.create_engine(db_uri)
        
        Base = get_base()
        Base.metadata.create_all(engine)
        
        self.Session = sqlalchemy.orm.sessionmaker(bind=engine, expire_on_commit=False)
        
        # If database does not exist, create a new database
        # with PSI4 columns filled in
        if not db_exists:
            self._import_psi4_fnctls_disp()

    def _get_session(self) -> sqlalchemy.orm.Session:
        return self.Session()
    
    def insert_base_functional(self, 
                                f_alias: str,
                                f_name: str,
                                f_citation : str,
                                f_desc: str,
                                f_data: dict,
                                src: str) -> None:
        
        # Check if master name already there
        
        
        
        
        
        pass
        
    def insert_multi_functional(self, 
                                f_alias: str,
                                f_name: str,
                                f_citation : str,
                                f_desc: str,
                                f_data: dict,
                                src: str) -> None:
        pass
    
    
    
    def _resolve_source(self, err_source: str) -> str:
        '''
            This functional is called when the query
            cannot find something in ERR_SOURCE, and 
            returns the next best source to search in.
        '''
        
        lower_resol_arr = list(map(lambda s: s.lower(), self.source_fallback_stack))
        err_source = err_source.lower()
        
        resol_index = -1
        if err_source in lower_resol_arr:
            resol_index = lower_resol_arr.index(err_source)
        
        next_source = resol_index + 1
        if (next_source == len(lower_resol_arr)):
            raise SourceResolutionError("ERROR: Cannot resolve source!")

        return lower_resol_arr[next_source]            
    
    def _add_dispersion_alias(self, dispersion_name: str) -> str:
        '''
            Adds an alias of a dispersion to the database.
        '''
    
    def _resolve_dispersion_alias(self, dispersion_name: str) -> str:
        '''
            Returns the canonical name of a dispersion.
        '''
        logger.warning("Not implemented yet")
        
    def _add_functional_alias(self, functional_name: str) -> str:
        '''
            Adds an alias to the database.
        '''
    
    def _resolve_functional_alias(self, functional_name: str) -> str:
        '''
            Returns the canonical name of a functional.
        '''
        logger.warning("Not implemented yet")
        return functional_name
    
    def get_dispersion(self,
                       functional_name: str,
                       dispersion_name: str,
                       source: str):
        '''
            Returns a list of base dispersions associated with a
            single dispersion model
        '''
        
        # Resolve both functional and dispersion alias first
        proper_fname = self._resolve_functional_alias(functional_name)
        proper_dname = self._resolve_dispersion_alias(dispersion_name)
        
        with self._get_session() as session:
            disp_list : list[tuple[DispersionBase, float]] = (
                session.query(DispersionBase, DispersionConfig.subdisp_coef)
                    .join(Functional, DispersionConfig.fnctl_id == Functional.fnctl_id)
                    .join(Sources, DispersionConfig.source == Sources.id) 
                    .join(DispersionBase, DispersionBase.subdisp_id == DispersionConfig.subdisp_id)
                    # .all()
                    .filter(
                        sqlalchemy.func.lower(Functional.fnctl_name) == functional_name.lower(),
                        sqlalchemy.func.lower(Sources.name) == source.lower(),
                        sqlalchemy.func.lower(DispersionConfig.disp_name) == dispersion_name.lower()                        
                    ).all()       
            )
            
            logger.warning(disp_list)
            
            # Check for results first
            
            # Try to resolve to another source if not found
            if len(disp_list) == 0:
                try:
                    next_source = self._resolve_source(source)
                except SourceResolutionError as e:
                    raise DBNotFoundError("Cannot find functional!") from e
                
                
            
            
            
            
            session.commit()
        
        

    def get_single_functional(self, 
                            functional_name: str, 
                            source: str) -> Functional:
        '''
            Queries only a single functional (without getting 
            all the coeffs).
            
            If the source is not 
            found, it just resolves to the fallback sources, 
            in priority.
            
            Raises:
                DBNotFoundError: If the functional is not present in any source.
        '''
        
        # Resolve alias first (canonical name is always in alias list)
        res_functional_name = self._resolve_functional_alias(functional_name)
        
        # Find funcitonal
        with self._get_session() as session:
            query_funcs : Functional = (
                session.query(Functional)
                .join(Sources)
                .filter(sqlalchemy.func.lower(Functional.fnctl_name) == functional_name.lower(),
                        Sources.name == source)  
                .all()
            )
            
            if len(query_funcs) > 1:
                raise DBDuplicateError("Error: More than one functional found on a single source")
            
            if len(query_funcs) == 0:
                try:
                    next_source = self._resolve_source(source)
                except SourceResolutionError as e:
                    raise DBNotFoundError("Cannot find functional!") from e
                
                return self.get_single_functional(functional_name, next_source)
                
            return query_funcs                    
    
    def get_functional(self, 
                             functional_name: str, 
                             source: str) -> (
                                tuple[Functional | None, list[tuple[Functional, float]]]
                             ):
        '''
            Queries a functional and returns the functional data
            and (if any) returns a list of tuples containing 
            the subcomponent functional and coefficients.
            
            The subcomponent functional is guranteed to be 
            base functionals.
            
            If queried functional is a base functional, 
            the function returns an empty list.
        '''
        
        # Resolve alias first
        
        # Query base functional 
        query_func = self.get_single_functional(functional_name, source)
    
        # If not lcom, we are done!
        if not query_func.is_lcom:
            return query_func, []
        
        parent_func = query_func
        parent_func_id = query_func.fnctl_id
        # If it is lcom, we need more work!
        with self._get_session() as session:        
            
            func_coefs : list[tuple[Functional, float]] = (
                            session.query(Functional, FunctionalCoeffs.coef)
                            .join(FunctionalCoeffs, FunctionalCoeffs.parent_fnctl_id == Functional.fnctl_id)
                            .filter(
                                FunctionalCoeffs.parent_fnctl_id == parent_func_id,
                            
                            )  
                            .all()
                        )
            
            # Validate, cannot have a lcom functional in the result
            if any(func.is_lcom for (func, _) in func_coefs):
                raise RuntimeError("Error: Somehow we have a multifunctional here")

            session.commit()

        return parent_func, func_coefs
    
    
    
    

