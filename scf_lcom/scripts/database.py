from sqlalchemy import Column, String, Double, ForeignKey, UUID, JSON
from sqlalchemy.orm import declarative_base, relationship
import uuid

Base = declarative_base()

class Source(Base):
    source_id = Column(String(36), primary_key=True, default=lambda: str(uuid.uuid4()))
    name = Column(String, unique=True, nullable=False) 
    
class Functional(Base):
    '''
        Functional table, contains the PSI4 compatible functional data
    ''' 
    source = Column(String(40), nullable=False)
    func_id = Column(String(36), primary_key=True, default=lambda: str(uuid.uuid4()))

class FunctionalCoeffs(Base):
    '''
        Table containing functional linear combination coefficients.
    '''

class Dispersion(Base):
    source = Column(String(40), primary_key=True)
    fnctl_id = Column(primary_key=True, ForeignKey())


class FunctionalDatabase(Base):
    '''
        Database wrapper for querying functionals and dispersion.
    '''    
    


