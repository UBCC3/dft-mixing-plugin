import pytest
import os
import yaml

from scripts.database import FunctionalDatabase

class TestBase:
    @pytest.fixture
    def initialize_db(self, tmp_path):
        # Initializes the database for the first time
        config_path = tmp_path / 'test_data_config.yaml'
        db_path = tmp_path / 'test.db'        

        config = {"db_path": db_path}
        
        # Load in database path
        with open(config_path, 'w') as f:
            config = yaml.safe_dump(config)
        
        func_db = FunctionalDatabase(config_path)

        yield func_db
        
    def test_base_queries(self, initialize_db):
        
        func_db : FunctionalDatabase = initialize_db

        # Test queries
        func_db.
        
        return     
    

# class TestFunctional:
    
#     @pytest.fixture
#     def test_load_functional_exceptions():
#         pass
    
#     def test_load_fully_success():
#         pass
    
#     def test_queries():
#         pass

# class TestDispersion:
    
#     @pytest.fixture
#     def load_functionalf_dispersion():
#         pass
    
#     def test_load_fully_success():
#         pass
    
#     def test_queries():
#         pass



    
    

