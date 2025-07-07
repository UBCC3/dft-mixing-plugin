import pytest
import os
import yaml

import logging

from psi4.driver.procrouting.dft.dft_builder import functionals

from scripts.database import FunctionalDatabase
from .db_sample_data import lcom_func_dataset

logger = logging.getLogger(__name__)

class TestBase:
    @pytest.fixture
    def initialize_db(self, tmp_path):
        # Initializes the database for the first time
        config_path = tmp_path / 'test_data_config.yaml'
        db_path = tmp_path / 'test.db'        

        config = {"db_path": str(db_path)}
        
        # Load in database path
        with open(config_path, 'w') as f:
            config = yaml.safe_dump(config, f)
        
        func_db = FunctionalDatabase(config_path)

        yield func_db
        
    def test_base_queries(self, initialize_db):        
        func_db : FunctionalDatabase = initialize_db

        # Test queries
        func_dict = func_db.query_functional_disp('b3lyp', source='psi4')
        
        logging.info(func_dict)
        
        ref_dict = functionals['b3lyp']
        for attr in ref_dict:
            if attr != 'dispersion':
                assert (func_dict[attr] == func_dict[attr]), f"Error, mismatch at {attr}"    
        
        print(func_dict) 
        
    def test_base_w_dispersion(self, initialize_db):
        
        func_db : FunctionalDatabase = initialize_db

        # Test queries
        disp = 'd2'
        
        ref_dict = functionals['b3lyp-d2']
        logger.warning(ref_dict)
        
        func_dict = func_db.query_functional_disp('b3lyp', disp, source="psi4")
        logger.warning(func_dict)

        assert ('dispersion' in func_dict), "No dispersion found!"

        for attr in ref_dict:
            if attr != 'dispersion':
                assert (func_dict[attr] == func_dict[attr]), f"Error, mismatch at {attr}"    
            
            else: 
                disp = func_dict['dispersion']
                for disp_attr in ref_dict['dispersion']:
                    assert(disp[disp_attr] == ref_dict[disp_attr])
                    
class TestFunctional:
    
    @pytest.fixture
    def test_load_fully_success(tmp_path):
        # Initializes the database for the first time
        config_path = tmp_path / 'test_data_config.yaml'
        db_path = tmp_path / 'test.db'        

        config = {"db_path": str(db_path)}
        
        # Load in database path
        with open(config_path, 'w') as f:
            config = yaml.safe_dump(config, f)
        
        func_db = FunctionalDatabase(config_path)
        func_db._store_fnctl_coef_json("test", lcom_func_dataset)

        yield func_db
    
    def test_queries(test_load_fully_success):
        
        func_db = test_load_fully_success
        
    
        
        
        
        
        pass

class TestDispersion:
    
    @pytest.fixture
    def load_functional_dispersion():
        pass
    
    def test_queries():
        pass



    
    

