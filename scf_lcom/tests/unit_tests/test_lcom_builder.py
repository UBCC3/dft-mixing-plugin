# Isolated builder into a class, so unit testing is way more easier
from scripts.lcom_superfunctional import LCOMSuperFunctionalBuilder
import psi4
import pytest


from psi4.driver.procrouting.dft.dft_builder import functionals

def compare_libxc_func(funcA, funcB):
    '''
        Compares two LibXC functionals if they are equal.
        (Right now, hard to compare if tweaks are equal)
    '''
    
    assert (funcA.name() == funcB.name())
    assert (funcA.alpha() == pytest.approx(funcB.alpha(), abs=1e-6))
    assert (funcA.omega() == pytest.approx(funcB.omega(), abs=1e-6))

# Tweaks can only be tested externally for now.


# Smoketests, on single, make sure single input still works
# in PSI4

# Single LCOM functional
def test_create_single_xc_existing():    
    
    functional = functionals['blyp']

    # Compare energies
    lcom_energy = psi4.energy('scf_lcom', dft_functional=functional)
    ref_energy = psi4.energy('scf', dft_functional=functional)
    
    assert (lcom_energy == pytest.approx(ref_energy, abs=1e-6))

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

# Minor feature, nested LCOMs,
# probably not will be supported since DB is 
# already taking care of it
def test_nested_lcom():
    pass

