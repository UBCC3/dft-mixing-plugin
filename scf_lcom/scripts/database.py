import os
from pathlib import Path

from sqlalchemy import Column, String, Double, ForeignKey, UUID, JSON, UniqueConstraint
from sqlalchemy.orm import declarative_base, relationship
from sqlalchemy.ext.hybrid import hybrid_property

# For functionals
from psi4.driver.procrouting.dft.dft_builder import functionals

# Default dispersions
from qcengine.programs.empirical_dispersion_resources import dashcoeff, get_dispersion_aliases, new_d4_api

import uuid
import re

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
    
    def __init__(self, database_config : str | Path = None):
        '''
            Creates a functional database from scratch, or reads from
            an existing database configuration (in a yaml file.)
        '''     
        self.source_resol_stack : list[str] = []
        
        # Import from existing database if path is specified 
        if database_config is not None:
            return
        
        # Otherwise, load from PSI4, then followed by sdftd3
        self._import_psi4_fnctls_disp()
        
        pass
    
    def get_source_resolution_stack(self) -> list[str]:
        return self.source_resol_stack.copy()
    
    def add_json_source(fnctl_coef_json, disp_param_json=None, disp_mixing_json=None):
        pass
    
    def _import_psi4_fnctls_disp(self) -> None:
        '''
            Helper function to import PSI4 dispersions and functions
            into the database.
        '''
        self.source_resol_stack.append("psi4")
        with self.Session() as session:
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
            
            with self.Session() as session:
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
            with self.Session() as session:
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
            
            with self.Session() as session:
                
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
                

    
    def query_functional(self, functional_name: str, disp_config: str | None = None) -> dict:
        '''
            Query a functional within a database. 
            Returns a dictionary that is compatible with the 
            PSI4 dft_functionals interface or the scf_lcom
            plugin interface.
        '''
        self.source_resol_stack.append("psi4")

        func_dict = {}
        func_dict["name"]

        
        
        
        
        if disp_config is None:
            return func_dict
        
        
        




        pass