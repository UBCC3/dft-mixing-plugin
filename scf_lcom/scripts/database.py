from sqlalchemy import Column, String, Double, ForeignKey, UUID, JSON
from sqlalchemy.orm import declarative_base, relationship
from sqlalchemy.ext.hybrid import hybrid_property
import uuid

Base = declarative_base()

class Source(Base):
    source_id = Column(String(36), primary_key=True, default=lambda: str(uuid.uuid4()))
    name = Column(String, unique=True, nullable=False) 
    
class Functional(Base):
    '''
        Functional table, contains the PSI4 compatible functional data
    ''' 
    source = Column(String(36), ForeignKey('sources.id'), nullable=False)
    fnctl_id = Column(String(36), primary_key=True, default=lambda: str(uuid.uuid4()))
    fnctl_data = Column(JSON(), nullable=True, default= None)
    fnctl_name = Column(String(40), nullable=False, unique=True)
    citation = Column(String(512), nullable=True, default="")

    @hybrid_property
    def is_lcom(self):
        return self.fnctl_data is None

class FunctionalCoeffs(Base):
    '''
        Table containing functional linear combination coefficients.
    ''' 
    parent_fnctl_id = Column(String(36), ForeignKey('functional.fnctl_id'), primary_key=True)
    child_fnctl_id = Column(String(36), ForeignKey('functional.fnctl_id'), primary_key=True)
    coef = Column(Double(), nullable=False)

class Dispersion(Base):
    source = Column(String(36), ForeignKey('sources.id'), nullable=False)
    citation = Column(String(512), nullable=True, default="")
    fnctl_id = Column(String(36), ForeignKey('functional.fnctl_id'), primary_key=True)
    disp_id = Column(String(36), primary_key=True, default=lambda: str(uuid.uuid4()))
    disp_params = Column(JSON(), nullable=True)
    disp_name = Column(String(40), nullable=False)

class Source(Base):
    id = Column(String(36), primary_key=True, default=lambda: str(uuid.uuid4()))
    name = Column(String(40), unique=True, nullable=False)


