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

import uuid
import re
import warnings

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

class FunctionalDatabase:
    
    def __init__(self, database_config : Path):
        '''
            Creates a functional database from scratch, or reads from
            an existing database configuration (in a yaml file.)
        '''     
        self.source_resol_stack : list[str] = []
        
        # Import from existing database if path is specified 
        if database_config is None:
            raise ValueError("Error: database configurations cannot be None")
        
        # Otherwise, load from PSI4, then followed by sdftd3
        self._import_psi4_fnctls_disp()
        
        pass
    
    def get_source_resolution_stack(self) -> list[str]:
        return self.source_resol_stack.copy()
    
    def add_json_source(source: str, 
                        fnctl_base_json : dict,
                        fnctl_coef_json : dict,  
                        disp_param_json : dict = None, 
                        disp_mixing_json: dict = None):        
        
        # Need to add feedback if process went successfully, or
        # some functionals 
        pass
    
    
    def _get_session(self) -> sqlalchemy.orm.Session:
        return self.Session()
    
    def _resolve_fnctl_coefficients(self, source: str, multi_fnctl_components: dict):
        raise NotImplementedError("Have not implemented this yet!")
        
    
    # Stores
    def _store_fnctl_base_json(self, 
                        source: str,
                        fnctl_base_json : dict,) -> None:
        pass

    def _store_disp_json(self, 
                        source: str,
                        disp_param_json : dict = None, 
                        disp_mixing_json: dict = None) -> None:
        pass
    
    def _store_fnctl_coef_json(self, 
                               source: str,
                               fnctl_coef_json: dict) -> None:
        '''
            Imports functionals from a JSON file and also validates
            the JSON file format.
            
            We accept the following JSON format
            ```json
                data: 
                    {
                        <multi_fnctl_name>:
                            // Mandatory
                            // This is a dictionary of functional components and their coeffcients
                            "fnctls": {
                                <fnctl_name>: coef,
                                ...
                            },    
                            
                            // Optional, citation of functional
                            // Defaults to an empty string
                            "citation": str,
                            
                            // Optional, description of functional
                            // Defaults to an empty string
                            "description": str
                    }
                    
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
            src_id = src.id
            session.add(src)
            session.commit()
          
        parent_functionals = {
            functionals[func_name] \
                for func_name in functionals \
                if "dispersion" not in functionals[func_name]
        }
        
        dispersion_functionals = {
            functionals[func_name] \
                for func_name in functionals \
                if "dispersion" in functionals[func_name]
        }
        
        # Add parent functionals        
        for func_name, func_dict in parent_functionals.items():
            func_name = func_dict["name"]
            func_citation = func_dict.get("citation", "No citations available")
            func_description = func_dict.get("description", "")
            func_data = func_dict
            
            with self._get_session() as session:
                fnctl = Functional(
                    source = src_id,
                    fnctl_data=func_data,
                    fnctl_name=func_name,
                    citation=func_citation,
                    description=func_description,
                )
                
                session.add(fnctl)
                session.commit()
            
        # Now add dispersions, only add those whose parent functionals we added before
        for func_name, func_dict in dispersion_functionals.items():
            func_name = func_dict["name"]
            
            # Extract the parent part of fnctl
            pattern = r"([-\d\w()]+)-([\w\d()]+)\)?$"
            match = re.match(pattern, func_name)
            if not match: 
                raise RuntimeError(f"Error: somehow encountered a non dash-coeff functional {func_name}")
            
            parent_func, _ = match.groups()
            # Search for parent_func id
            with self._get_session() as session:
                fnctl : Functional = (session.query(Functional)
                            .filter(Functional.fnctl_name.lower() == parent_func.lower())
                            .first())
                
                if fnctl is None:
                    raise RuntimeError(f"Error: no parent functional corresponds to the functional {func_name}")
                
                func_id = fnctl.fnctl_id
                session.commit()
            
            func_disp = func_dict["dispersion"]
            disp_name = func_disp["type"]
            disp_citation = func_disp.get("citation", "No citations available")
            disp_description = func_disp.get("description", "")
            disp_params = func_disp.get("params", {})
            
            with self._get_session() as session:
                
                # Add base dispersion
                subdisp = DispersionBase(
                    fnctl_id=func_id,
                    disp_params=disp_params,
                    subdisp_name=disp_name,
                    source = src_id,
                    citation=disp_citation,
                    description=disp_description,
                )
                
                subdisp_id = subdisp.subdisp_id
                
                # Add a configuration consisting of only ONE dispersion.
                # under the name of the dispersion itself so PSI4 can run 
                # singular dispersion model
                
                disp_config = DispersionConfig(
                    fnctl_id=func_id,
                    source=src_id,
                    citation=disp_citation,
                    func_description=disp_description,
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
                            .filter(Functional.fnctl_name == functional_name,
                                    Sources.name == source)  
                            .first()
                        )
            
            session.commit()
        pass
    
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
            
        with self._get_session() as session:
            disp_list : list[tuple[DispersionBase, float]] = (
                session.query(DispersionBase, DispersionConfig.subdisp_coef)
                    .join(Functional, DispersionConfig.fnctl_id == Functional.fnctl_id)
                    .join(Sources, DispersionConfig.source == Sources.name) 
                    .join(DispersionBase, DispersionBase.subdisp_id == DispersionConfig.subdisp_id)
                    .filter(
                        Functional.fnctl_name == functional_name,
                        Sources.name == source,
                        DispersionConfig.disp_name == dispersion_name                        
                    ).all()       
            )
            
            session.commit()
            
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
            
            
        if len(dispersion_coeffs) > 1:
            lcom_disps = []
            for (base_disp, coef) in dispersion_coeffs:
                disp_dict = {}
                disp_dict["params"] = base_disp.disp_params
                disp_dict["type"] = base_disp.subdisp_name
                disp_dict["citation"] = base_disp.citation
                disp_dict["description"] = base_disp.description
                disp_dict["lcom_coef"] = coef
                lcom_disps.append(disp_dict)
            
            func_dict["lcom_dispersion"] = lcom_disps
        
        return func_dict
    
   