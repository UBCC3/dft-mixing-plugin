import logging.config
import os
from pathlib import Path

import sqlalchemy.orm
from sqlalchemy import Column, String, Double, ForeignKey, UUID, JSON, UniqueConstraint
from sqlalchemy.orm import declarative_base, relationship
from sqlalchemy.ext.hybrid import hybrid_property

import uuid

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
    func_name = Column(String(36), nullable=False)
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
    disp_base_source = Column(String(36), ForeignKey('sources.id'), nullable=False)
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
    disp_config_source = Column(String(36), ForeignKey('sources.id'), nullable=False)
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
    func_name =  Column(String(36), nullable=False)
    disp_name = Column(String(36), nullable=False)
    alias_name = Column(String(40), nullable=False)
    
def get_base():
    return Base