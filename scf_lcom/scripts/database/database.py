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

import uuid
import re
import warnings
import json

import logging

# Config parser
import yaml

logger = logging.getLogger(__name__)

Base = declarative_base() 
class Sources(Base):
    __tablename__ = "sources"

    id = Column(String(36), primary_key=True, default=lambda: str(uuid.uuid4()))
    name = Column(String(40), unique=True, nullable=False)

    # Reverse access to related rows
    functionals = relationship("Functional", back_populates="source_ref")
    dispersion_configs = relationship("DispersionConfig", back_populates="source_ref")
    dispersion_bases = relationship("DispersionBase", back_populates="source_ref")

class Functional(Base):
    __tablename__ = "functional"

    fnctl_id = Column(String(36), primary_key=True, default=lambda: str(uuid.uuid4()))
    source = Column(String(36), ForeignKey('sources.id'), nullable=False)
    fnctl_data = Column(JSON(), nullable=True, default=None)
    fnctl_name = Column(String(40), nullable=False, unique=True)
    citation = Column(String(512), nullable=False, default="")
    description = Column(String(512), nullable=False, default="")

    # Navigation
    source_ref = relationship("Sources", back_populates="functionals")
    child_coeffs = relationship("FunctionalCoeffs", foreign_keys="[FunctionalCoeffs.parent_fnctl_id]", back_populates="parent")
    parent_coeffs = relationship("FunctionalCoeffs", foreign_keys="[FunctionalCoeffs.child_fnctl_id]", back_populates="child")
    dispersion_configs = relationship("DispersionConfig", back_populates="functional")
    dispersion_bases = relationship("DispersionBase", back_populates="functional")

    @hybrid_property
    def is_lcom(self):
        return self.fnctl_data is None

class FunctionalAlias(Base):
    __tablename__ = "functionalalias"
    alias_id = Column(String(36), primary_key=True, default=lambda: str(uuid.uuid4()))
    fnctl_id = Column(String(36), ForeignKey('functional.fnctl_id'), nullable=False)
    alias_name = Column(String(40), nullable=False)

class FunctionalCoeffs(Base):
    __tablename__ = "functionalcoeffs"

    parent_fnctl_id = Column(String(36), ForeignKey('functional.fnctl_id'), primary_key=True)
    child_fnctl_id = Column(String(36), ForeignKey('functional.fnctl_id'), primary_key=True)
    coef = Column(Double(), nullable=False, default=1.0)

    # Navigation
    parent = relationship("Functional", foreign_keys=[parent_fnctl_id], back_populates="child_coeffs")
    child = relationship("Functional", foreign_keys=[child_fnctl_id], back_populates="parent_coeffs")

class DispersionBase(Base):
    __tablename__ = "dispersionbase"

    subdisp_id = Column(String(36), primary_key=True, default=lambda: str(uuid.uuid4()))
    fnctl_id = Column(String(36), ForeignKey('functional.fnctl_id'), nullable=False)
    disp_params = Column(JSON(), nullable=True)
    
    # Note that this cannot be the aliased name, must be the master name.
    subdisp_name = Column(String(40), nullable=False)
    source = Column(String(36), ForeignKey('sources.id'), nullable=False)
    citation = Column(String(512), nullable=False, default="")
    description = Column(String(512), nullable=False, default="")

    # Navigation
    source_ref = relationship("Sources", back_populates="dispersion_bases")
    functional = relationship("Functional", back_populates="dispersion_bases")
    used_in_configs = relationship("DispersionConfig", back_populates="subdisp")

class DispersionConfig(Base):
    __tablename__ = "dispersionconfig"

    disp_id = Column(String(36), primary_key=True, default=lambda: str(uuid.uuid4()))
    fnctl_id = Column(String(36), ForeignKey('functional.fnctl_id'), primary_key=True)
    source = Column(String(36), ForeignKey('sources.id'), nullable=False)
    citation = Column(String(512), nullable=False, default="")
    description = Column(String(512), nullable=False, default="")
    disp_name = Column(String(40), nullable=False)
    subdisp_id = Column(String(36), ForeignKey('dispersionbase.subdisp_id'), nullable=False)
    subdisp_coef = Column(Double(), nullable=False, default=1.0)

    __table_args__ = (
        UniqueConstraint('disp_name', 'subdisp_id', name='uq_disp_subdisp'),
    )

    # Navigation
    source_ref = relationship("Sources", back_populates="dispersion_configs")
    functional = relationship("Functional", back_populates="dispersion_configs")
    subdisp = relationship("DispersionBase", back_populates="used_in_configs")

class DispersionAlias(Base):
    __tablename__ = "dispersionalias"
    alias_id = Column(String(36), primary_key=True, default=lambda: str(uuid.uuid4()))
    disp_config_id = Column(String(36), ForeignKey('dispersionconfig.fnctl_id'), nullable=False)
    alias_name = Column(String(40), nullable=False)

class FunctionalDatabase:
    
    def __init__(self, database_config : Path):
        '''
            Creates a functional database from scratch, or reads from
            an existing database configuration (in a yaml file.)
            
            Yaml contains:
            1. Path to database
            
        '''     
        self.source_resol_stack : list[str] = []
        
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
        
        Base.metadata.create_all(engine)
        
        self.Session = sqlalchemy.orm.sessionmaker(bind=engine, expire_on_commit=False)
        
        # If database does not exist, create a new database
        # with PSI4 columns filled in
        if not db_exists:
            self._import_psi4_fnctls_disp()
    
    def get_source_resolution_stack(self) -> list[str]:
        return self.source_resol_stack.copy()
    
    def add_json_source(self, source: str, 
                        fnctl_base_json : dict,
                        fnctl_coef_json : dict,  
                        disp_param_json : dict = None, 
                        disp_mixing_json: dict = None):        
        
        # Need to add feedback if process went successfully, or
        # some functionals 
        pass
    
    def query_functional_disp(self, func_name : str, 
                            disp_config_name : str = None,
                            source: str = None, 
                            fmt: str = "psi4") -> dict:
        '''
            Queries for a functional, and returns it in a specified format
            (PSI4 JSON format by default). 
            
            If `disp_config_name` is given as an argument, it will
            add that dispersion configuration to the queried
            functional.
            
            By default, the database queries for the top of the 
            source stack by default, and keep progressing
            
            For example, if the source stack was specified as:
                - sourceA
                - dftd4     (hidden)
                - psi4      (hidden)
                
            it will query for sourceA first, then dftd4 if not found, then psi4.
            
            If `source` was given (reccomended), it will immediately search for
            that source. 
        '''
        
        par_func, functional_coefs = self._get_only_functional(func_name, source)
        
        if disp_config_name is not None:
            disp_coefs = self._get_dispersion(func_name, disp_config_name, source)
        else:
            disp_coefs = []
            
        # Assemble functional and dispersion
        if fmt.lower() == "psi4":
            return self._format_to_psi4(par_func, functional_coefs, disp_coefs)

        raise RuntimeError("Error: Invalid format!")


    def _get_session(self) -> sqlalchemy.orm.Session:
        return self.Session()
    
    def _resolve_fnctl_coefficients(self, source: str, multi_fnctl_components: dict):
        raise NotImplementedError("Have not implemented this yet!")        
    
    # Stores
    def _store_fnctl_base_json(self, 
                        source: str,
                        fnctl_base_json : dict,) -> None:
        raise NotImplementedError("Error: Have not implemented this yet...")

    def _store_disp_json(self, 
                        source: str,
                        disp_param_json : dict = None) -> None:
        
        error_list = []
        
        # Add all base dispersions first
        for parent_func, avail_disps in disp_param_json.items():
            for disp_name, disp_info in avail_disps.items():
                try:
                    with self._get_session() as session:
                        # Query for functional id
                        fnctl_id, source_id = (
                            session.query(Functional.fnctl_id, Sources.id)
                            .join(Sources)
                            .filter(Functional.fnctl_name == parent_func,
                                    Sources.name == source)
                            .first()
                        )

                        # Create new Base Disperison
                        base_disp = DispersionBase(
                            fnctl_id = fnctl_id, 
                            subdisp_name = disp_info["type"],
                            disp_params = disp_info["disp_params"],
                            citation = disp_info.get("citation", ""),
                            description = disp_info.get("citation", "")
                        )
                        
                        # Create new DispersionConfig
                        base_disp_cfg = DispersionBase(
                            fnctl_id = fnctl_id, 
                            subdisp_name = disp_info["type"],
                            citation = disp_info.get("citation", ""),
                            description = disp_info.get("citation", ""),
                            disp_name = disp_info["type"],
                            subdisp_id = base_disp.subdisp_id
                        )
                        
                        session.add(base_disp)
                        session.add(base_disp_cfg)
                        session.commit()
                        
                except Exception as e:
                    error_list.append(e)    
                    
        if len(error_list > 1):
            raise ExceptionGroup("[ERROR] Several errors were found: ", error_list)
        
    def _store_disp_mixing(self, 
                            source: str,
                            disp_coeffs_json : dict = None) -> None:
            
            error_list = []
            
            # Add all base dispersions first
            for parent_func, avail_disps in disp_coeffs_json.items():
                for disp_name, disp_info in avail_disps.items():
                    try:
                        with self._get_session() as session:
                            # Query for functional id
                            fnctl_id, source_id = (
                                session.query(Functional.fnctl_id, Sources.id)
                                .join(Sources)
                                .filter(Functional.fnctl_name == parent_func,
                                        Sources.name == source)
                                .first()
                            )

                            # Create new Base Disperison
                            base_disp = DispersionBase(
                                fnctl_id = fnctl_id, 
                                subdisp_name = disp_info["type"],
                                disp_params = disp_info["disp_params"],
                                citation = disp_info.get("citation", ""),
                                description = disp_info.get("citation", "")
                            )
                            
                            # Create new DispersionConfig
                            base_disp_cfg = DispersionBase(
                                fnctl_id = fnctl_id, 
                                subdisp_name = disp_info["type"],
                                citation = disp_info.get("citation", ""),
                                description = disp_info.get("citation", ""),
                                disp_name = disp_info["type"],
                                subdisp_id = base_disp.subdisp_id
                            )
                            
                            session.add(base_disp)
                            session.add(base_disp_cfg)
                            session.commit()
                            
                    except Exception as e:
                        error_list.append(e) 
    
    def _store_fnctl_coef_json(self, 
                               source: str,
                               fnctl_coef_json: dict) -> None:
        '''
            Imports functionals from a JSON file and also validates
            the JSON file format.
            

            ```
        '''
        # Since we are doing this as a single transaction,  
        # need to rollback the entire process
        error_lists : list[Exception] = []

        for multi_fnctl_name, multi_args in fnctl_coef_json.items():
            # Typecheck the types
            assert type(multi_fnctl_name) == str, "Error: incompatible type for fnctl_name"
            
            components = multi_args["fnctls"]
            
            # Create the new functional first, then queue up.
            new_fnctl : Functional = Functional(
                    source=source,
                    fnctl_name=multi_fnctl_name, 
                    citation=multi_args.get("citation", ""),
                    description=multi_args.get("description", "")
                )
                
            multi_fnctl_id = new_fnctl.fnctl_id
                        
            for child_fnctl_name, coef in components.items():
                child_fnctl_name = child_fnctl_name.lower()
                try:
                    with self._get_session() as session:
                        
                        # Lookup existing functionals.
                        child_fnctl : Functional = (session.query(Functional)
                                .filter(Functional.fnctl_name.lower() == child_fnctl_name)
                                .first()
                        )
                        
                        if child_fnctl is None:
                            raise RuntimeError(f"Error: Cannot find functional {child_fnctl_name} in database!")
                        
                        # !TODO: for future support, resolve coefficients.                 
                        if child_fnctl.is_lcom:
                            raise RuntimeError(f"Error: Cannot create multi-functional by specifying multi-functional, need base functionals.")
                        
                        child_id = child_fnctl.fnctl_id
                        new_coef_entry = FunctionalCoeffs(
                            parent_fnctl_id = multi_fnctl_id,
                            child_fnctl_id = child_id, 
                            coef=coef
                        )
                        
                        session.add(new_coef_entry)
                
                except Exception as e:
                    error_lists.append(e)
        
        if len(error_lists) > 0:
            raise ExceptionGroup("Some errors were encountered:", error_lists)            
        return

    def _import_psi4_fnctls_disp(self) -> None:
        '''
            Helper function to import PSI4 dispersions and functions
            into the database.
        '''
        self.source_resol_stack.append("psi4")
        with self._get_session() as session:
            src = Sources(name="psi4")
            session.add(src)
            session.commit()
            src_id = src.id
          
        logger.info(f"src_id: {src_id}")
    
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
        

        func_master_set = {}
        dispersion_master_set = {}
        
        # Add parent functionals        
        for func_name, func_dict in parent_functionals.items():
            alias_name = func_name
            master_name = func_dict["name"]
            func_citation = func_dict.get("citation", "No citations available")
            func_description = func_dict.get("description", "")
            func_data = func_dict
            
            if master_name.lower() not in func_master_set:
                with self._get_session() as session:

                    fnctl = Functional(
                        source = src_id,
                        fnctl_data=master_name,
                        fnctl_name=func_name,
                        citation=func_citation,
                        description=func_description,
                    )
                    
                    session.add(fnctl)
                    session.commit()

                func_master_set[master_name.lower()] = fnctl.fnctl_id
        
            # Add an alias
            else: 
                with self._get_session() as session:
                    alias = FunctionalAlias(
                        fnctl_id = func_master_set[master_name.lower()],
                        alias_name = alias_name
                    )
                    session.add(alias)
                    session.commit()
            
        error_fnctls = {}
            
        # Now add dispersions, only add those whose parent functionals we added before
        for func_name, func_dict in dispersion_functionals.items():
            
            # Skip blacklisted functionals for now
            if func_name.lower() in blacklist_funcs:
                continue
            
            # For dispersion, ALWAYS link towards master name
            # func_name is with dashcoeff
            func_name = func_dict["name"]
            
            # Extract the parent part of fnctl
            pattern = r"([-\d\w()]+)-([\w\d()]+)\)?$"
            match = re.match(pattern, func_name)
            
            # END immediately if in blacklisted functionals
            # if func_name.lower() in blacklist_funcs:
            #     logger.error(f'Skipping functional {func_name} since the dispersion is weird')
            #     continue
            
            if not match: 
                logger.error(f"Dispersion causing error, Skipping functional {func_name}:{func_dict} ")
                error_fnctls[func_name] = func_dict
                continue
                # raise RuntimeError(f"Error: somehow encountered a non dash-coeff functional {func_name}")
            
            parent_func, disp_alias = match.groups()
            # Search for parent_func id
            with self._get_session() as session:
                fnctl : Functional = (session.query(Functional)
                            .filter(sqlalchemy.func.lower(Functional.fnctl_name) == parent_func.lower())
                            .first())
                
                if fnctl is None:
                    logger.error(f"Skipping functional {func_name}:{func_dict}")
                    session.commit()
                    error_fnctls[func_name] = func_dict
                    continue
                    raise RuntimeError(f"Error: no parent functional corresponds to the functional {func_name} : {func_dict}")
                
                func_id = fnctl.fnctl_id
                session.commit()
            
            func_disp = func_dict["dispersion"]
            disp_name = func_disp["type"]
            disp_citation = func_disp.get("citation", "No citations available")
            disp_description = func_disp.get("description", "")
            disp_params = func_disp.get("params", {})
            
            with self._get_session() as session:
                
                subdisp_id = str(uuid.uuid4())
                
                # Add base dispersion
                subdisp = DispersionBase(
                    subdisp_id=subdisp_id,
                    fnctl_id=func_id,
                    disp_params=disp_params,
                    subdisp_name=disp_name,
                    source = src_id,
                    citation=disp_citation,
                    description=disp_description,
                )
                
                # Add a configuration consisting of only ONE dispersion.
                # under the name of the dispersion itself so PSI4 can run 
                # singular dispersion model
                
                disp_config = DispersionConfig(
                    fnctl_id=func_id,
                    source=src_id,
                    citation=disp_citation,
                    description=disp_description,
                    disp_name=disp_name,
                    subdisp_id=subdisp_id,    
                )
                
                session.add(subdisp)
                session.add(disp_config)
                session.commit()
        
    def _get_only_functional(self, 
                             functional_name: str, 
                             source: str) -> (
                                tuple[Functional, list[tuple[Functional, float]]]
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
        
        with self._get_session() as session:        
            query_func : Functional = (
                            session.query(Functional)
                            .join(Sources)
                            .filter(sqlalchemy.func.lower(Functional.fnctl_name) == functional_name.lower(),
                                    Sources.name == source)  
                            .first()
                        )
            
            session.commit()
    
        if query_func is None:
            raise RuntimeError('Error: Cannot find functional.')
    
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
    
    def _get_dispersion(self, functional_name: str, dispersion_name: str, source: str):
            
        logger.warning(f"Querying: fname: {functional_name}, disp: {dispersion_name}, source: {source} ")
        
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
            
            session.commit()
            
        if len(disp_list) == 0:
            raise RuntimeError("Error: Dispersion Not Found!")
             
        return disp_list
            
    def _format_to_psi4(self, 
                        par_functional : Functional,
                        functional_coeffs: list[tuple[Functional, float]],
                        dispersion_coeffs: list[tuple[DispersionBase, float]] = []):
        '''
            Formats query results to PSI4 lcom plugin format.
        '''    
        
        func_dict = {}
        func_dict["name"] = par_functional.fnctl_name
        func_dict["citation"] = par_functional.citation
        func_dict["description"] = par_functional.description
        
        if len(functional_coeffs) > 1:
            lcom_fucntionals = {}
            for (child_fnctl, coef) in functional_coeffs:
                child_name = child_fnctl.fnctl_name
                child_data = child_fnctl.fnctl_data.copy()
                lcom_fucntionals[child_name] = child_data
                lcom_fucntionals[child_name]["coef"] = coef    
            
            func_dict["lcom_functionals"] = lcom_fucntionals
        
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
    
   