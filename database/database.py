import logging.config
import os
from pathlib import Path

import sqlalchemy.orm
from sqlalchemy import Column, String, Double, ForeignKey, UUID, JSON, UniqueConstraint
from sqlalchemy.orm import declarative_base, relationship, Session
from sqlalchemy.ext.hybrid import hybrid_property

# Default dispersions
from qcengine.programs.empirical_dispersion_resources import dashcoeff, get_dispersion_aliases, new_d4_api

# Error types
from .db_errors import (
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
        
        # Havent imported dftd4 dispersions yet    
        self.source_fallback_stack : list[str] = ["psi4"]
        self._is_empty = True
        
        # Import from existing database if path is specified 
        if database_config is None:
            raise ValueError("Error: database configurations cannot be None")
      
        # Import configs
        with open(database_config) as f:
            config = yaml.safe_load(f)

        if (config.get("load_dftd3", False)):
            self.source_fallback_stack.insert(0, "dftd3")

        # Configure source resolution
        DB_EXTERNAL_RESOL : list = config.get('external_resol', [])
        if (len(DB_EXTERNAL_RESOL) > 0):
            while len(DB_EXTERNAL_RESOL) != 0:
                last_priority = DB_EXTERNAL_RESOL.pop(-1)
                self.source_fallback_stack.insert(0, last_priority)

        DB_PATH = config["db_path"]
        
        # Check if db exists
        db_exists = os.path.exists(DB_PATH)
        if (db_exists): 
            self._is_empty = False
        
        # Configure base        
        db_uri = f'sqlite:///{DB_PATH}'
        # Create a SQL database for this
        engine = sqlalchemy.create_engine(db_uri)
        
        Base = get_base()
        Base.metadata.create_all(engine)
        
        self.Session = sqlalchemy.orm.sessionmaker(bind=engine, expire_on_commit=False)

    def is_empty(self):
        return self._is_empty

    def get_source_name(self, session: sqlalchemy.orm.Session,
                        source_id):
        src_name = (session.query(Sources.name)
                    .filter(Sources.id == source_id)
                    .first())
        
        if src_name is None:
            raise DBNotFoundError("Cannot find source")

        return src_name[0]
    
    def get_source_atomic(self, source_id):
        with self.get_session() as session:
            return self.get_source_name(session, source_id)


    def get_session(self) -> sqlalchemy.orm.Session:
        return self.Session()
    
    def _query_source(self, session: sqlalchemy.orm.Session, 
                        source_name: str) -> uuid.UUID:
        '''
            Queries for a source within the databse
            
            Returns:
                `src_id` id of source
            
            Raises:
                DBNotFoundError: if the source data is not found.
        '''
            
        src_id = (session.query(Sources.id)
                    .filter(sqlalchemy.func.lower(Sources.name) == source_name.lower())
                    .first())  
    
        if src_id is None:
            raise DBNotFoundError(f"Error, source is missing for {source_name}.")
        
        return src_id[0]
        
    def _create_source(self,
                        session: sqlalchemy.orm.Session, 
                        source_name: str) -> uuid.UUID:
        '''
            Create a new source/namespace for functionals/
            dispersions. 
            
            Returns:
                `src_id` id of source
                
            Raises:
        '''

        src_id = str(uuid.uuid4())
        new_src = Sources(
            id = src_id,
            name = source_name
        )
        
        session.add(new_src)
        return src_id
    
    def insert_base_functional(self, 
                                session: sqlalchemy.orm.Session,
                                canon_fname: str,
                                f_citation : str,
                                f_desc: str,
                                f_data: dict,
                                src: str) -> uuid.UUID:
        '''
            Creates a new base functional and returns the
            functional id for the base functional.
        '''
        
        try:
            src_id = self._query_source(session, src)
            
        # Cannot find source, create a new source
        # just for this functional.
        except DBNotFoundError as e:
            src_id = self._create_source(session, src)
        
        new_functional_id = str(uuid.uuid4())
        # logger.warning(f'ACTUALLY INSERT {f_name} into DB')
        new_functional = Functional(
            fnctl_id = new_functional_id,
            source= src_id,
            fnctl_data = f_data, 
            fnctl_name = canon_fname,
            citation = f_citation,
            description = f_desc    
        )
        
        session.add(new_functional)
        # logger.warning(f'{f_name} into DB OK')
        
        return new_functional_id      
        
    def insert_multi_functional(self, 
                                session: sqlalchemy.orm.Session,
                                multi_name: str,
                                multi_citation : str,
                                multi_desc: str,
                                multi_coeffs: dict[str, float],
                                src: str) -> uuid.UUID:
        
        # This is a functional consisting of its dispersion

        # Create new functional entry for the multifunctional
        multi_id = self.insert_base_functional(
            session, multi_name, multi_citation, 
            multi_desc, None, src
        )
    
        # Now, loop through the coefficients and lookup
        # functionals
        for child_name, coef in multi_coeffs.items():
            child_func = self.get_single_functional(session, child_name, src)
            
            # Add new coef entry
            coef_entry = FunctionalCoeffs(
                parent_fnctl_id = multi_id,
                child_fnctl_id = child_func.fnctl_id, 
                coef=coef
            )
            
            session.add(coef_entry)

        return multi_id
                
    def insert_single_disp(self, 
                           session: sqlalchemy.orm.Session,
                            func_name: str,
                            disp_type: str,
                            disp_citation : str,
                            disp_desc: str,
                            disp_params: dict,
                            disp_src: str,
                            func_src: str | None = None) -> uuid.UUID:
    
        try:
            src_id = self._query_source(session, disp_src)
            
        # Cannot find source, create a new source
        # just for this functional.
        except DBNotFoundError as e:
            src_id = self._create_source(session, disp_src)
    
        # Check if functional exists
        func_id = self.get_single_functional(session, 
                                             func_name, 
                                             func_src).fnctl_id
    
        base_disp_id = uuid.uuid4()
        
        # If functional does exist (does not throw error), insert the dispersion info
        base_disp = DispersionBase(
            subdisp_id = str(base_disp_id),
            fnctl_id = str(func_id),
            disp_params = disp_params,
            subdisp_name = disp_type,
            disp_base_source = str(src_id),
            citation = disp_citation,
            description = disp_desc,
        )
        
        disp_config_id = uuid.uuid4()
        
        # Also add a new dispersion config
        disp_config = DispersionConfig(
            disp_id = str(disp_config_id),
            fnctl_id = str(func_id),
            disp_config_source = str(src_id), 
            citation = disp_citation,
            description = disp_desc,
            disp_name = disp_type,
            subdisp_id = str(base_disp_id)
        )
        
        session.add(base_disp)
        session.add(disp_config)
        
        return base_disp_id

    def insert_disp_config(self, 
                           session: sqlalchemy.orm.Session,
                            func_name: str,
                            disp_name: str,
                            disp_citation : str,
                            disp_desc: str,
                            disp_coeffs: dict[str, float],
                            disp_src: str,
                            func_src: str | None = None) -> None:
        
        try:
            src_id = self._query_source(session, disp_src)
            
        # Cannot find source, create a new sourceg
        # just for this functional.
        except DBNotFoundError as e:
            src_id = self._create_source(session, disp_src)
        
        # Get functional id
        functional = self.get_single_functional(session, 
                                            func_name, None)

        func_id = functional.fnctl_id

        # Add new dispersion entries
        for subdisp_name, coef in disp_coeffs.items():
            
            subdisp_alias = f'{func_name}-{subdisp_name}'
            # logger.warning(subdisp_alias)
            
            # Assume dispersion part of the same configuration,
            # otherwise, fallback to existing dispersion.
            subdisp_id : str = self.get_base_dispersion(session, subdisp_alias,
                                          disp_source=disp_src,
                                          func_source=func_src).subdisp_id
            
            disp_coef_id = str(uuid.uuid4())
            disp_entry = DispersionConfig(
                disp_id = str(disp_coef_id),
                fnctl_id = str(func_id), 
                disp_config_source = src_id,
                citation = disp_citation,
                description = disp_desc,
                disp_name = disp_name, 
                subdisp_id = str(subdisp_id),
                subdisp_coef = coef 
            ) 
            
            session.add(disp_entry)
        
        
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
    
    def add_dash_coeff_mapping(self,
                              session : sqlalchemy.orm.Session, 
                              functional_name: str, 
                              dispersion_name: str,
                              dash_coeff_name: str):
        '''
            Adds an alias of a dispersion to the database.
        '''
        new_alias = DispersionAlias(
            func_name = functional_name,
            disp_name = dispersion_name, 
            func_dash_name = dash_coeff_name
        )
        
        session.add(new_alias)        
    
    
    def resolve_dash_coeff(self,
                            session: sqlalchemy.orm.Session, 
                            dash_coeff_name: str) -> tuple[str, str]:
        '''
            Resolve dash coefficient name to resolve
            into their canonical functional name, and
            canonical dispersion name.
            
            Returns:
                `(canon_func_name, canon_disp_name)`, canonical
                functional name and dispersion respectively.
        '''
        
        res = (session.query(DispersionAlias.func_name, DispersionAlias.disp_name)
                .filter(sqlalchemy.func.lower(DispersionAlias.func_dash_name)
                         == dash_coeff_name.lower())
                .first())
        
        if res is None:
            raise DBNotFoundError("Error: Cannot resolve dashcoeff")
        
        return (res[0], res[1])        
    
    def add_functional_alias(self, 
                              session: sqlalchemy.orm.Session,
                              functional_name: str, 
                              functional_alias : str):
        '''
            Adds a functional alias to the database,
            maps it to the canonical name.
        '''
        # logger.warning(f"INSERTING ALIAS {functional_alias} -> {functional_name}")
        new_alias = FunctionalAlias(
            func_name = functional_name,
            alias_name = functional_alias
        )
        
        session.add(new_alias)      
        
    
    def resolve_functional_alias(self, 
                                  session,
                                  functional_name: str) -> str:
        '''
            Returns the canonical name of a functional.
        '''
        canon_name = ( session.query(FunctionalAlias.func_name)
                        .filter(sqlalchemy.func.lower(FunctionalAlias.alias_name) == functional_name.lower())
                        .first())

        if canon_name is None:
            raise DBNotFoundError(f"Error: Cannot resolve alias {functional_name}!")
        
        return canon_name[0]   
    
    # ==================================
    #          GETTERS (BY NAME)
    #
    # ==================================
    def get_base_dispersion(self,
                       session: sqlalchemy.orm.Session,
                       dashcoeff_name: str,
                       func_source: str | None = None,
                        disp_source: str | None = None) -> DispersionBase:
        '''
            Returns a list of base dispersions associated with a
            single dispersion model
            
            Source Resolution Notes:
                Here, the functional source resolution takes 
                takes priority over the dispersion source resolution.
        '''
        
        # Resolve dashcoeff name
        canon_fname, canon_dname = self.resolve_dash_coeff(session, dashcoeff_name)
        
        # Resolve functional first
        functional = self.get_single_functional(session, canon_fname,
                                      func_source)
        
        
        canon_func_source = str(functional.source)
        func_id = functional.fnctl_id
        
        if disp_source is None:
            disp_source = self.source_fallback_stack[0]  
            
        source_id = self._query_source(session, disp_source)
        
        base_disp = (
            session.query(DispersionBase)
            .filter(
                DispersionBase.fnctl_id == func_id, 
                sqlalchemy.func.lower(DispersionBase.subdisp_name) == canon_dname.lower(),
                DispersionBase.disp_base_source == source_id
            )
            .first() 
        )
        
        if base_disp is None:
            try:
                next_source = self._resolve_source(disp_source)
            except SourceResolutionError as e:
                raise DBNotFoundError("Cannot find functional!") from e
            
            return self.get_base_dispersion(session,
                                            dashcoeff_name,
                                            canon_func_source,
                                            next_source)
             
        return base_disp
    
    def get_multi_dispersion(self,
                       session: sqlalchemy.orm.Session,
                       dashcoeff_name: str,
                       func_source: str | None = None,
                        disp_source: str | None = None
                        ) -> list[tuple[DispersionBase, float]]:
        '''
            Returns a list of base dispersions associated with a
            single dispersion model
        '''
        
        # Resolve dashcoeff name
        canon_fname, canon_dname = self.resolve_dash_coeff(session, dashcoeff_name)
        
        # Resolve functional first
        functional = self.get_single_functional(session, canon_fname,
                                      func_source)
        
        logger.warning(f"QUERY DISP? f:{canon_fname} d:{canon_dname}")
        
        canon_func_source = str(functional.source)
        func_id = functional.fnctl_id
        
        if disp_source is None:
            logger.warning(f"FIRST OPTION {self.source_fallback_stack[0]}")
            disp_source = self.source_fallback_stack[0]  
            
        source_id = self._query_source(session, disp_source)
        
        disp_coeffs = (
            session.query(DispersionBase, DispersionConfig.subdisp_coef)
            .join(
                DispersionConfig, 
                DispersionConfig.subdisp_id == DispersionBase.subdisp_id
            )
            .filter(
                DispersionConfig.fnctl_id == func_id, 
                sqlalchemy.func.lower(DispersionConfig.disp_name) == canon_dname.lower(),
                DispersionConfig.disp_config_source == source_id
            )
            .all() 
        )
        
        if len(disp_coeffs) == 0:
            try:
                next_source = self._resolve_source(disp_source)
            except SourceResolutionError as e:
                raise DBNotFoundError("Cannot find dispersion coefficients!") from e
            
            return self.get_multi_dispersion(session,
                                            dashcoeff_name,
                                            canon_func_source,
                                            next_source)
            
        # logger.warning(disp_coeffs)
        return disp_coeffs
    
    def get_single_functional(self, 
                            session: sqlalchemy.orm.Session,
                            functional_name: str, 
                            source: str | None = None) -> Functional:
        '''
            Queries only a single functional (without getting 
            all the coeffs).
            
            If the source is not 
            found, it just resolves to the fallback sources, 
            in priority.
            
            Raises:
                DBNotFoundError: If the functional is not present in any source.
        '''
        
        if source is None: 
            source = self.source_fallback_stack[0]
        
        # Resolve alias first (canonical name is always in alias list)
        res_functional_name = self.resolve_functional_alias(session, 
                                                            functional_name)
        
        # Find funcitonal
    
        query_funcs : list[Functional] = (
            session.query(Functional)
            .join(Sources)
            .filter(sqlalchemy.func.lower(Functional.fnctl_name) == res_functional_name.lower(),
                    sqlalchemy.func.lower(Sources.name) == source.lower())  
            .all()
        )
        
        if len(query_funcs) > 1:
            raise DBDuplicateError("Error: More than one functional found on a single source")
        
        if len(query_funcs) == 0:
            try:
                next_source = self._resolve_source(source)
            except SourceResolutionError as e:
                raise DBNotFoundError("Cannot find functional!") from e
            
            return self.get_single_functional(session, functional_name, next_source)
            
        return query_funcs[0]

    def get_functional(self, 
                        session: sqlalchemy.orm.Session, 
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
        
        if source is None:
            source = self.source_fallback_stack[0]
    
        # Query base functional 
        query_func = self.get_single_functional(session, functional_name, source)
        logger.info(f"QUERY_FUNC LCOM? {query_func.is_lcom}")
    
        # If not lcom, we are done!
        if not query_func.is_lcom:
            return query_func, []
        
        parent_func = query_func
        parent_func_id = query_func.fnctl_id
        # If it is lcom, we need more work!
                   
        func_coefs : list[tuple[Functional, float]] = (
                        session.query(Functional, FunctionalCoeffs.coef)
                        .join(FunctionalCoeffs, FunctionalCoeffs.child_fnctl_id == Functional.fnctl_id)
                        .filter(
                            FunctionalCoeffs.parent_fnctl_id == parent_func_id,        
                        )  
                        .all()
                    )
            
        # Validate, cannot have a lcom functional in the result
        if any(func.is_lcom for (func, _) in func_coefs):
            raise RuntimeError("Error: Somehow we have a multifunctional here")

        return parent_func, func_coefs
    
    
    
    # ==================================
    #          GETTERS (BY ID)
    #
    # ==================================
    def get_single_functional_by_id(self, 
                            session: sqlalchemy.orm.Session,
                            functional_id: uuid.UUID) -> Functional:
        '''
            Queries only a single functional (without getting 
            all the coeffs) by id.
            
            Raises:
                DBNotFoundError: If the functional is not present in any source.
        '''
        # Find funcitonal
    
        query_func: Functional | None = (
            session.query(Functional)
            .filter(Functional.fnctl_id == functional_id)  
            .first()
        )
    
        if query_func is None:
            raise DBNotFoundError("Cannot find functional!")
            
        return query_func                    
    
    def get_functional_by_id(self, 
                            session: sqlalchemy.orm.Session, 
                            functional_id: uuid.UUID) -> (
                                tuple[Functional | None, list[tuple[Functional, float]]]
                            ):
        '''
            Queries a functional by ID and returns the functional data
            and (if any) returns a list of tuples containing 
            the subcomponent functional and coefficients.
            
            The subcomponent functional is guranteed to be 
            base functionals.
            
            If queried functional is a base functional, 
            the function returns an empty list.
        '''
        
        # Query base functional 
        query_func = self.get_single_functional_by_id(session, functional_id)
    
        # If not lcom, we are done!
        if not query_func.is_lcom:
            return query_func, []
        
        parent_func = query_func
        parent_func_id = query_func.fnctl_id
        # If it is lcom, we need more work!
                   
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

        return parent_func, func_coefs