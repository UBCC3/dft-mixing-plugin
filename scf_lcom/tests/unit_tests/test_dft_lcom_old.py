import numpy as np

import pytest

import psi4

from psi4.driver.procrouting.dft.dft_builder import build_superfunctional_base
from psi4.driver.procrouting.dft.hyb_functionals import get_pw6b95_tweaks

pytestmark = [pytest.mark.psi, pytest.mark.api]

def compare_sf_hf(a, b):
    """ Returns true if two Superfunctionals have the same HF exchange """
    if a.x_alpha() != b.x_alpha():                             return False
    if a.x_omega() == 0 and b.x_omega() == 0:                  return True # If long range is turned off, we're done
    if a.x_omega() != b.x_omega() or a.x_beta() != b.x_beta(): return False
    return True
def compare_sf_mp2(a, b):
    """ Returns true if two Superfunctionals have the same MP2 correlation """
    if a.c_alpha() != b.c_alpha():                             return False
    if a.c_ss_alpha() != b.c_ss_alpha():                       return False
    if a.c_os_alpha() != b.c_os_alpha():                       return False
    if a.c_omega() != b.c_omega():                             return False # Psi4 actually never sets this; Should never be different
    return True
def compare_fnctl(a, b):
    """ Returns true if two functionals are identical """
    if a.name() != b.name():                                   return False
    if a.alpha() != b.alpha() or a.omega() != b.omega():       return False
    # TODO: Tweaks!!!
    return True
    
def compare_fnctls_uni(a, b):
    """ Returns true if all functionals in a appear exactly once in b """
    for af in a.functionals():
        matching = 0
        for bf in b.functionals():
            if compare_fnctl(af, bf): matching += 1
        if matching != 1: return False
    return True
def compare_fnctls(a, b):
    """ Returns true if all functionals in a and b Superfunctionals are identical """
    return compare_fnctls_uni(a, b) and compare_fnctls_uni(b, a)

def real_print_superfunctional(f):
    """ Prints a superfunctional. Actually. """
    print(f.name())
    print("hf_alpha =", f.x_alpha())
    print("hf_omega =", f.x_omega())
    print("hf_beta =", f.x_beta())
    print("mp2_alpha =", f.c_alpha())
    print("mp2_omega =", f.c_omega())
    print("mp2_ss_alpha =", f.c_ss_alpha())
    print("mp2_os_alpha =", f.c_os_alpha())
    for fn in f.functionals():
        print(fn.name() + ":")
        print("   ", "alpha =", fn.alpha())
        print("   ", "omega =", fn.omega())
        for (k, v) in fn.tweaks().items():
            print("   ", "tweak" + k, "=", v)
    print()

def make_sf(sf_def): return build_superfunctional_base(sf_def, False)
def add_fnctls_to_first(first, *nextsf):
    (first, d) = make_sf(first)
    for sf in nextsf:
        for fnctl in make_sf(sf)[0].functionals():
            first.add_functional(fnctl)
    return (first, d)
def add_lxc_fnctls(fbase = { "name": "fnctl", "x_hf": { "alpha": 0 }, "c_mp2": { "alpha": 0 } }, **lxcs):
    (sf, d) = make_sf(fbase)
    for (name, extra) in lxcs.items():
        (sfn, d) = make_sf({
            "name": name,
            "xc_functionals": { name: extra }
        })
        fn = sfn.functionals()[0]
        fn.set_alpha(extra["alpha"] if "alpha" in extra else 1.0)
        if "omega" in extra: fn.set_omega(extra["omega"])
        sf.add_functional(fn)
    return (sf, d)

@pytest.mark.scf
@pytest.mark.dft
@pytest.mark.parametrize("name,fa,fb", [
    # Copy tests... Just make sure we didn't really mess anything up
    pytest.param(
        "COPY_HF",
        {
            "name": "test-a",
            "x_hf": { "alpha": 1.0, },
            "c_functionals": {},
        },
        {
            "name": "test-b",
            "lcom_functionals": { "HF": 1.0, },
        }
    ),
    pytest.param(
        "COPY_HF_MP2",
        {
            "name": "test-a",
            "x_hf": { "alpha": 1.0, },
            "c_mp2": { "alpha": 1.0, },
        },
        {
            "name": "test-b",
            "lcom_functionals": { "MP2MP2": 1.0, }
        }
    ),
    pytest.param(
        "COPY_LXC",
        {
            "name": "test-a",
            "xc_functionals": { "HYB_GGA_XC_B97": {} }
        },
        {
            "name": "test-b",
            "lcom_functionals": { "B97-0": 1.0, }
        }
    ),
    pytest.param(
        "COPY_X_C",
        {
            "name": "test-a",
            "x_functionals": { "GGA_X_B88": {} },
            "c_functionals": { "GGA_C_LYP": {} },
        },
        {
            "name": "test-b",
            "lcom_functionals": { "BLYP": 1.0, }
        }
    ),

    # Make sure coefs greater than one work. I think there's a print bug only making these all =1 in the out file
    pytest.param(
        "GT1_HF",
        {
            "name": "test-a",
            "x_hf": { "alpha": -10.0, },
            "c_functionals": {},
        },
        {
            "name": "test-b",
            "lcom_functionals": { "HF": -10.0, },
        }
    ),
    pytest.param(
        "GT1_HF_MP2",
        {
            "name": "test-a",
            "x_hf": { "alpha": 5.0, "beta": 0.0 },
            "c_mp2": { "alpha": -10.0, },
        },
        {
            "name": "test-b",
            "lcom_functionals": { "HF": 15.0, "MP2MP2": -10.0, }
        }
    ),
    pytest.param(
        "GT1_X_C",
        {
            "name": "test-a",
            "x_functionals": { "GGA_X_B88": { "alpha": -10.0 } },
            "c_functionals": { "GGA_C_LYP": { "alpha": -10.0 } },
        },
        {
            "name": "test-b",
            "lcom_functionals": { "BLYP": -10.0, }
        }
    ),

    # Basic tests just merge two known functionals
    pytest.param(
        "BASIC_HF_MP2",
        {
            "name": "test-a",
            "x_hf": { "alpha": 1.0, },
            "c_mp2": { "alpha": 0.75, },
        },
        {
            "name": "test-b",
            "lcom_functionals": { "HF": 0.25, "MP2MP2": 0.75 }
        }
    ),
    pytest.param(
        "BASIC_X_C",
        {
            "name": "test-a",
            "x_functionals": { "GGA_X_B88": {} },
            "c_functionals": {
                "GGA_C_LYP": { "alpha": 0.25 },
                "GGA_C_OP_B88": { "alpha": 0.75 }
            },
        },
        {
            "name": "test-b",
            "lcom_functionals": { "BLYP": 0.25, "BOP": 0.75 }
        }
    ),
    pytest.param(
        "BASIC_LXC",
        add_lxc_fnctls(
            fbase = { 
                "name": "fnctl",
                "x_hf": { "alpha": 0, "beta": 0.75, "omega": 0.4 },
                "c_mp2": { "alpha": 0 }
            },
            MGGA_XC_ZLP = { "alpha": 0.25 },
            HYB_GGA_XC_WB97 = { "alpha": 0.75 }
        ),
        {
            "name": "test-b",
            "lcom_functionals": { "ZLP": 0.25, "wB97": 0.75 }
        }
    ),

    pytest.param(
        "TWEAK_NO_OVERWRITE",
        add_fnctls_to_first(
            {
                "name": "mPWPW_with_PW6B95_hf",
                "x_functionals": { "GGA_X_mPW91": { "alpha": 0.5 } },
                "c_functionals": { "GGA_C_PW91": { "alpha": 0.5 } },
                "x_hf": { "alpha": 0.14 },
            },
            {
                "name": "PW6B95_no_hf",
                "x_functionals": {
                    "GGA_X_MPW91": {
                        "tweak": get_pw6b95_tweaks(),
                        "alpha": 0.36
                    }
                },
                "c_functionals": {
                    "MGGA_C_BC95": {
                        "alpha": 0.5,
                        "tweak": {
                            "_css": 0.03668,
                            "_copp": 0.00262,
                        },
                    }
                },
            },
        ),
        {
            "name": "test-b",
            "lcom_functionals": { "PW6B95": 0.5, "mPWPW": 0.5 }
        }
    ),
    pytest.param(
        "TWEAK_MERGE",
        {
            "name": "test-a",
            "x_functionals": {
                "GGA_X_MPW91": {
                    "tweak": get_pw6b95_tweaks(),
                    "alpha": 0.36
                }
            },
            "x_hf": {
                "alpha": 0.14
            },
            "c_functionals": {
                "MGGA_C_BC95": {
                    "tweak": {
                        "_css": 0.03668,
                        "_copp": 0.00262,
                    },
                }
            },
        },
        {
            "name": "test-b",
            "lcom_functionals": {
                "PW6B95": 0.5,
                "PW6B95_c": {
                    "coef": 0.5,
                    "x_functionals": {},
                    "c_functionals": {
                        "MGGA_C_BC95": {
                            "tweak": {
                                "_css": 0.03668,
                                "_copp": 0.00262,
                            },
                        }
                    },
                }
            }
        }
    ),
])
def test_lcom_hf_basic(name, fa, fb, request):
    (fa, da) = fa if isinstance(fa, tuple) else make_sf(fa)
    (fb, db) = fb if isinstance(fb, tuple) else make_sf(fb)

    if da or db: raise "Uhhh... None of these should be generating dispersion"

    # Debugging doesn't actually work here, since Psi4 redirects this print god knows where
    #fa.print_detail(10)
    #fb.print_detail(10)

    real_print_superfunctional(fa)
    real_print_superfunctional(fb)

    assert(compare_sf_hf(fa, fb))
    assert(compare_sf_mp2(fa, fb))
    assert(compare_fnctls(fa, fb))

@pytest.fixture
def dft_bench_systems():
    ang = np.array([
      [ -1.551007,  -0.114520,   0.000000],
      [ -1.934259,   0.762503,   0.000000],
      [ -0.599677,   0.040712,   0.000000]])
    oldang = ang * 0.52917721067 / 0.52917720859
    oldmass = [15.99491461956, 1.00782503207, 1.00782503207]

    h2o = psi4.core.Molecule.from_arrays(geom=oldang, elez=[8, 1, 1], units='Angstrom', mass=oldmass)

    psi4.set_options({
       'dft_radial_points': 200,
       'dft_spherical_points': 590,
       'guess': 'sad',
       'e_convergence': 9,
       'd_convergence': 9,
       'basis': '6-31G',
    }, verbose=10)

    return h2o

# Here I just double check that beta in fact doesn't matter when omega is zero
@pytest.mark.scf
@pytest.mark.dft
def test_lcom_prereq_beta_when_omega_zero(dft_bench_systems, request):
    fa = {
        "name": "test-a",
        "x_hf": { "alpha": 1.0, "beta": 0.0, "omega": 0.0 },
        "c_mp2": { "alpha": 0, },
    }
    fb = {
        "name": "test-b",
        "x_hf": { "alpha": 1.0, "beta": 1.0, "omega": 0.0 },
        "c_mp2": { "alpha": 0, },
    }
    a = psi4.energy('scf', dft_functional=fa, molecule=dft_bench_systems)
    b = psi4.energy('scf', dft_functional=fb, molecule=dft_bench_systems)
    assert compare_values(b, a, 4, request.node.name), (a - b)

# After talking in the Psi4 dev channel, it appears that LibXC doesn't actually spec the meaning
# of the alpha parameter. If a functional is used in a LCOM, it **absolutely must** behave consistently
# with respect to alpha and omega. Omega is kindof guaranteeded by how Psi4 copies the omega param;
# But we need to make sure the alpha parameter is used as a master coefficient. This test verifies that
# by setting alpha to zero and verifying that this is the same as not having the fnctl altogether.
# It's faster and more reliable than carefully checking all of LibXC.
@pytest.mark.scf
@pytest.mark.dft
@pytest.mark.parametrize("name,dest", [
    pytest.param("HYB_MGGA_X_M06_HF", None),
    pytest.param("MGGA_C_M06_HF", None),
    pytest.param("HYB_MGGA_X_M06_2X", None),
    pytest.param("MGGA_C_M06_2X", None),
    pytest.param("HYB_GGA_X_N12_SX", None),
    pytest.param("GGA_C_N12_SX", None),
    pytest.param("MGGA_X_MN15_L", None),
    pytest.param("MGGA_C_MN15_L", None),
    pytest.param("MGGA_X_MVS", None),
    pytest.param("GGA_C_REGTPSS", None),
    pytest.param("MGGA_X_MN12_L", None),
    pytest.param("MGGA_C_MN12_L", None),
    pytest.param("GGA_X_GAM", None),
    pytest.param("GGA_C_GAM", None),
    pytest.param("HYB_MGGA_X_MN12_SX", None),
    pytest.param("MGGA_C_MN12_SX", None),
    pytest.param("MGGA_X_M11_L", None),
    pytest.param("MGGA_C_M11_L", None),
    pytest.param("HYB_MGGA_X_MN15", None),
    pytest.param("MGGA_C_MN15", None),
    pytest.param("MGGA_X_TM", None),
    pytest.param("MGGA_C_TM", None),
    pytest.param("HYB_MGGA_X_M08_HX", None),
    pytest.param("MGGA_C_M08_HX", None),
    pytest.param("HYB_MGGA_X_M06", None),
    pytest.param("MGGA_C_M06", None),
    pytest.param("HYB_GGA_XC_CAM_B3LYP", None),
    pytest.param("MGGA_X_MS2", None),
    pytest.param("HYB_MGGA_X_M11", None),
    pytest.param("MGGA_C_M11", None),
    pytest.param("MGGA_X_PKZB", None),
    pytest.param("MGGA_C_PKZB", None),
    pytest.param("HYB_MGGA_X_M08_SO", None),
    pytest.param("MGGA_C_M08_SO", None),
    pytest.param("MGGA_X_BLOC", None),
    pytest.param("MGGA_C_TPSSLOC", None),
    pytest.param("HYB_MGGA_X_M05_2X", None),
    pytest.param("MGGA_C_M05_2X", None),
    pytest.param("HYB_MGGA_X_BMK", None),
    pytest.param("GGA_C_BMK", None),
    pytest.param("GGA_X_N12", None),
    pytest.param("GGA_C_N12", None),
    pytest.param("GGA_X_PBE", None),
    pytest.param("GGA_C_PBE", None),
    pytest.param("MGGA_X_TPSS", None),
    pytest.param("MGGA_C_TPSS", None),
    pytest.param("HYB_MGGA_X_M05", None),
    pytest.param("MGGA_C_M05", None),
    pytest.param("GGA_X_OPTX", None),
    pytest.param("GGA_C_LYP", None),
    pytest.param("HYB_GGA_XC_WB97X", None),
    pytest.param("MGGA_X_M06_L", None),
    pytest.param("MGGA_C_M06_L", None),
    pytest.param("HYB_MGGA_X_MVSH", None),
    pytest.param("HYB_GGA_XC_LC_WPBE08_WHS", None),
    pytest.param("MGGA_X_MBEEF", None),
    pytest.param("GGA_C_PBE_SOL", None),
    pytest.param("HYB_GGA_XC_B97_2", None),
    pytest.param("HYB_GGA_XC_B97_3", None),
    pytest.param("MGGA_X_MS0", None),
    pytest.param("GGA_X_PW91", None),
    pytest.param("GGA_C_PW91", None),
    pytest.param("HYB_GGA_XC_WB97", None),
    pytest.param("HYB_GGA_X_SOGGA11_X", None),
    pytest.param("GGA_C_SOGGA11_X", None),
    pytest.param("GGA_XC_HCTH_93", None),
    pytest.param("GGA_XC_HCTH_407", None),
    pytest.param("GGA_X_MPW91", None),
    pytest.param("MGGA_C_BC95", None),
    pytest.param("HYB_GGA_XC_B97_K", None),
    pytest.param("LDA_C_PW", None),
    pytest.param("HYB_MGGA_XC_PWB6K", None),
    pytest.param("HYB_GGA_XC_LRC_WPBE", None),
    pytest.param("MGGA_X_REVTPSS", None),
    pytest.param("MGGA_C_REVTPSS", None),
    pytest.param("MGGA_X_TAU_HCTH", None),
    pytest.param("GGA_C_TAU_HCTH", None),
    pytest.param("GGA_X_MPW91", None),
    pytest.param("HYB_GGA_XC_HFLYP", None),
    pytest.param("HYB_MGGA_X_SCAN0", None),
    pytest.param("MGGA_C_SCAN", None),
    pytest.param("HYB_MGGA_X_TAU_HCTH", None),
    pytest.param("GGA_C_HYB_TAU_HCTH", None),
    pytest.param("HYB_MGGA_XC_TPSSH", None),
    pytest.param("HYB_MGGA_XC_REVTPSSH", None),
    pytest.param("GGA_X_B88", None),
    pytest.param("GGA_C_OP_B88", None),
    pytest.param("GGA_C_P86", None),
    pytest.param("HYB_GGA_XC_LRC_WPBEH", None),
    pytest.param("MGGA_X_SCAN", None),
    pytest.param("MGGA_X_MS1", None),
    pytest.param("HYB_MGGA_X_MS2H", None),
    pytest.param("GGA_X_RPW86", None),
    pytest.param("GGA_XC_HCTH_120", None),
    pytest.param("GGA_XC_HCTH_147", None),
    pytest.param("LDA_X", "x_functionals"),
    pytest.param("LDA_C_VWN", None),
    pytest.param("GGA_X_PBE_SOL", None),
    pytest.param("HYB_GGA_XC_B3LYP", None),
    pytest.param("HYB_GGA_XC_B3PW91", None),
    pytest.param("GGA_C_OP_PBE", None),
    pytest.param("HYB_GGA_XC_B3P86", None),
    pytest.param("GGA_X_SOGGA", None),
    pytest.param("GGA_X_PBE_R", None),
    pytest.param("HYB_GGA_XC_B97_1", None),
    pytest.param("HYB_GGA_XC_B97", None),
    pytest.param("HYB_GGA_XC_PBEH", None),
    pytest.param("GGA_X_RPBE", None),
])
def test_lcom_prereq_fnctl_alpha_turns_off(name, dest, dft_bench_systems, request):
    fnctl_blank = { "name": "blank", "x_hf": {}, "c_mp2": {} }
    fnctl_def = { "name": "base" }

    if dest is None: # Infer destination
        if "_XC_" in name: dest = "xc_functionals"
        elif "_X_" in name: dest = "x_functionals"
        elif "_C_" in name: dest = "c_functionals"
        else: raise "Unable to infer destination for functional "+name

    if dest != "xc_functionals":
        fnctl_def["x_functionals"] = {}
        fnctl_def["c_functionals"] = {}
        fnctl_def[dest] = {
            name: { "alpha": 0.0 }
        }
    else:
        fnctl_def[dest] = {
            name: { "_override_alpha": 0.0 }
        }

    (param_ref, d) = make_sf(fnctl_def)
    print(param_ref)

    # We're not looking at dispersion fnctls yet
    if d: raise "I didn't think about dispersion when I wrote this (TODO)"

    fnctl_blank["x_hf"]["alpha"] = param_ref.x_alpha()
    fnctl_blank["x_hf"]["beta"] = param_ref.x_beta()
    fnctl_blank["x_hf"]["omega"] = param_ref.x_omega()

    if param_ref.is_c_scs_hybrid():
        fnctl_blank["c_mp2"]["ss"] = param_ref.c_ss_alpha()
        fnctl_blank["c_mp2"]["os"] = param_ref.c_os_alpha()
    else:
        fnctl_blank["c_mp2"]["alpha"] = param_ref.c_alpha()
    
    real_print_superfunctional(make_sf(fnctl_blank)[0])
    real_print_superfunctional(make_sf(fnctl_def)[0])
    
    a = psi4.energy('scf', dft_functional=fnctl_blank, molecule=dft_bench_systems)
    b = psi4.energy('scf', dft_functional=fnctl_def, molecule=dft_bench_systems)
    assert compare_values(b, a, 4, request.node.name + ": " + name), (a - b)
