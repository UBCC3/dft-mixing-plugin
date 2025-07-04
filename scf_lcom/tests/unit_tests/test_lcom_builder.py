
from scripts.lcom_superfunctional import LCOMSuperFunctionalBuilder
import pytest


# Smoketests

def test_create_single_xc_existing():
    pass

def test_create_single_manual():
    pass

def test_create_single_xc():
    pass


# basic GGA, non hybrid
def test_lcom_manual_def():
    pass

def test_lcom_gga_basic():
    pass

def test_lcom_tweaks():
    pass

# hybrid, with HF component
def test_lcom_hybrid_gga():
    pass


# XC Related tests

# Important: test whether it computes functionals twice
def test_lcom_duplicate_xc_functional():
    pass

# Minor feature: decompose xc components
def test_lcom_duplicate_xc_decomp():
    pass

# Minor feature, not decompose xc with tweaks, just a performance test
# since this is experimental
def test_lcom_duplicate_xc_decomp():
    pass




