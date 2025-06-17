#
# @BEGIN LICENSE
#
# scf_lcom by Psi4 Developer, a plugin to:
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2024 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

import psi4
import psi4.driver.p4util as p4util
import psi4.driver.procrouting.dft.dft_builder as dft_builder
from psi4.driver.procrouting import proc_util
from scripts.dft_router import build_lcom_functional

run_scf = psi4.driver.procedures['energy']['scf']

def run_scf_lcom(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    scf_lcom can be called via :py:func:`~driver.energy`. For scf plugins.

    >>> energy('scf_lcom')
    
    This plugin only patches the dft builder to support linear combination
    of DFT superfunctionals.
    """
    
    # Monkey patch the scf function (for lcom dft suppport, then call run_scf) to avoid reimplementing
    # everything for scf energy calculation
    original_func = dft_builder.build_superfunctional_from_dictionary
    
    try:
        # replace the psi4 builder with our implementation
        dft_builder.build_superfunctional_from_dictionary = build_lcom_functional
        
        # Now just run scf normally
        return run_scf(name, kwargs)

    except Exception as e:
        raise e
        
    finally:
        # Regardless, restore the original function
        dft_builder.build_superfunctional_from_dictionary = original_func

# Integration with driver routines
psi4.driver.procedures['energy']['scf_lcom'] = run_scf_lcom
