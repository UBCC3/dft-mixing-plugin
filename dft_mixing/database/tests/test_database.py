import pytest
import os
import yaml
import re
import pprint

import logging

from psi4.driver.procrouting.dft.dft_builder import functionals

from dft_mixing.database import Psi4DbAdapter
from db_sample_data import (
    lcom_func_dataset,
    func_dataset_ref,
    lcom_psi_dispconfig_dataset,
    lcom_base_disp_dataset,
    multi_disp_ans
)

logger = logging.getLogger(__name__)

def compare_psi_disp(dict, ref_disp):
    for attr in ref_disp:
        if attr in {'citation', 'description', 'type', 'alias'}: continue
        if attr == 'name':
            assert dict[attr].lower() == ref_disp[attr].lower(), f"Mismatch at attribute {attr}"
        
        
        assert dict[attr] == ref_disp[attr], f"Mismatch at attribute {attr}"
    

def compare_lcom_dispersion(test_dict, ref_disp):
    normalized_test = normalize_name_dict_lowercase(test_dict)
    
    for subdisp_name in ref_disp:
        assert (subdisp_dict.lower() in normalized_test), f"Missing dispersion name {subdisp_name}"
        subdisp_dict = normalized_test[subdisp_name.lower()]
        compare_psi_disp(subdisp_dict, ref_disp[subdisp_name])
    
def compare_psi_functionals(test_dict, ref_dict):
    for attr in ref_dict:
        if attr in {'name', 'citation', 'description', 'alias'}: continue
        
        if attr == 'dispersion':
            ref_disp = ref_dict[attr]
            disp_dict = test_dict[attr]
            compare_psi_disp(disp_dict, ref_disp)        
            continue
        
        if attr == 'lcom_dispersion':
            ref_lcom_disp = ref_dict[attr]
            lcom_disp_dict = test_dict[attr]
            compare_lcom_dispersion(lcom_disp_dict, ref_lcom_disp)
            continue
        
        if (type(ref_dict[attr]) == str):
            assert test_dict[attr].lower() == ref_dict[attr].lower(), f"Mismatch at attribute {attr}"
            continue
        
        assert test_dict[attr] == ref_dict[attr], f"Mismatch at attribute {attr}"

def normalize_name_dict_lowercase(name_dict: dict):
    '''
        Returns a copy of name_dict with normalized 
        keys (lowercase).
    '''
    norm_dict = {}
    for name, val in name_dict.items():
        assert type(name) == str
        norm_dict[name.lower()] = val
        
    return norm_dict

def compare_lcom_functionals(func_dict, ref_dict):
    for attr in ref_dict:
        if attr == 'dispersion':
            ref_disp = ref_dict[attr]
            disp_dict = func_dict[attr]
            compare_psi_disp(disp_dict, ref_disp)        
            continue
    
        if attr == 'name':
            assert ref_dict['name'].lower() == func_dict['name'].lower()
        
        if attr != 'lcom_functionals': continue
        
        all_lcom_dict = normalize_name_dict_lowercase(func_dict['lcom_functionals'])
        for ref_fname, ref_fdict in ref_dict['lcom_functionals'].items():
            target_fdict = all_lcom_dict[ref_fname.lower()]
            if type(ref_fdict) == dict:
                assert (target_fdict['coef'] == pytest.approx(ref_fdict['coef']))
                compare_psi_functionals(target_fdict, ref_fdict)
            else:
                # Raw coeff
                assert(target_fdict == pytest.approx(ref_fdict))
    
class TestBase:
    @pytest.fixture(scope="class")
    def initialize_db(self, tmp_path_factory):
        
        tmp_path = tmp_path_factory.mktemp('db_test_data')
        # Initializes the database for the first time
        config_path = tmp_path / 'test_data_config.yaml'
        db_path = tmp_path / 'test.db'        

        config = {
            "db_path": str(db_path),
            "load_dftd3": False,
            "external_resol": [],
            "preload": []
        }
        
        # Load in database path
        with open(config_path, 'w') as f:
            config = yaml.safe_dump(config, f)
        
        func_db = Psi4DbAdapter()
        func_db.initialize_database(config_path)

        yield func_db
    
    @pytest.mark.parametrize(
        "func_name", [
            'b3lyp',        
            'op-pbe',
            'oppbe',
            'pbeop',
            'DSDPBEPBE',
            'DSD-PBEPBE'
        ]
    )       
    def test_base_queries(self, initialize_db, func_name):        
        func_db : Psi4DbAdapter = initialize_db
        func_name = func_name.lower()
        func_dict = func_db.get_functional_dict(func_name, functional_source='psi4')    
        logger.warning(func_dict)
        compare_psi_functionals(func_dict, functionals[func_name])

    @pytest.mark.parametrize(
        "dashcoeff_name", [
            'b3lyp-d',
            'b3lyp-d2',
            'b3lyp-d3bj',
            'b3lyp-d3zero',
            'b3lyp-d3', 
            'bp86-nl',
            'bp86-d',
            'bp86-d3bj',        
        ]
    )       
    def test_base_w_dispersion(self, initialize_db, dashcoeff_name):
        func_db : Psi4DbAdapter = initialize_db

        pattern = r"([-\d\w()]+)-([\w\d()]+)\)?$"
        match = re.match(pattern, dashcoeff_name)
        
        fname, dname = match.groups()
        
        ref_dict = functionals[dashcoeff_name]
        logger.warning(f"REF  {ref_dict}")
        
        func_dict = func_db.get_functional_dict(fname, dname, "psi4", "psi4")
        logger.warning(f"ACTUAL  {func_dict}")

        assert ('dispersion' in func_dict), "No dispersion found!"
        compare_psi_functionals(func_dict, ref_dict)
        
class TestFunctional:
    
    @pytest.fixture(scope="class")
    def initialize_db(self, tmp_path_factory):
        
        tmp_path = tmp_path_factory.mktemp('db_test_data')
        # Initializes the database for the first time
        config_path = tmp_path / 'test_data_config.yaml'
        db_path = tmp_path / 'test.db'        

        config = {
            "db_path": str(db_path),
            "load_dftd3": False,
            "external_resol": [],
            "preload": []
        }
        
        # Load in database path
        with open(config_path, 'w') as f:
            config = yaml.safe_dump(config, f)
        
        func_db = Psi4DbAdapter()
        func_db.initialize_database(config_path)
    
        func_db.load_multi_functional_data(lcom_func_dataset, 'test')
        yield func_db
    
    @pytest.mark.parametrize(
        'fname', list(func_dataset_ref.keys())
    )
    def test_queries(self, initialize_db, fname):
        ref_ans = func_dataset_ref[fname]
        logger.warning(f'REF: {pprint.pformat(ref_ans)}')
        func_db : Psi4DbAdapter = initialize_db
        
        ans =  func_db.get_functional_dict(fname, None, 'test')
        logger.warning(f'ANS: {ans}')
        compare_lcom_functionals(ans, ref_ans)
    
    # def test_not_found(initialize_db, fnames):

        
class TestDispersion:
    
    @pytest.fixture(scope="class")
    def initialize_db(self, tmp_path_factory):
        
        tmp_path = tmp_path_factory.mktemp('db_test_data')
        # Initializes the database for the first time
        config_path = tmp_path / 'test_data_config.yaml'
        db_path = tmp_path / 'test.db'        

        
        config = {
            "db_path": str(db_path),
            "load_dftd3": False,
            "external_resol": ["test"],
            "preload": []
        }
        
        # Load in database path
        with open(config_path, 'w') as f:
            config = yaml.safe_dump(config, f)
        
        func_db = Psi4DbAdapter()
        func_db.initialize_database(config_path)

        func_db.load_multi_functional_data(lcom_func_dataset, 'test')
        func_db.load_base_dispersion_data(lcom_base_disp_dataset, 'test')
        func_db.load_dispersion_config_data(lcom_psi_dispconfig_dataset, 'test')
        
        yield func_db
    
    @pytest.mark.parametrize(
        'fname, disp_name', 
        [
            # ("TPSS", "disp_mix1"),
            # ("TPSS", "disp_mix2"),
            ("singlefunc1", "disp_mix1"),
            ("singlefunc1", "disp_mix2"),
            ("multifunc2", "disp_mix1"),
            ("multifunc2", "disp_mix2"),
        ]
    )
    def test_multidisp_query(self, initialize_db, fname, disp_name):
        ref_ans = multi_disp_ans[(fname, disp_name)]
        func_db : Psi4DbAdapter = initialize_db
        ans = func_db.get_functional_dict(fname, disp_name, 'test', 'test')
        logger.warning(f'DISP DICT {pprint.pformat(ans)}')
        compare_lcom_functionals(ans, ref_ans)
    

