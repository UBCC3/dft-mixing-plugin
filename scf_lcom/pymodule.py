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
from psi4.driver.procrouting import proc_util

def run_scf_lcom(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    scf_lcom can be called via :py:func:`~driver.energy`. For scf plugins.

    >>> energy('scf_lcom')
    
    This plugin creates a LCOMSuperFunctional and computes its energy

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    # Your plugin's psi4 run sequence goes here
    psi4.core.set_local_option('SCF_LCOM', 'PRINT', 1)

    # Build a new blank wavefunction to pass into scf
    scf_lcom_molecule = kwargs.get('molecule', psi4.core.get_active_molecule())

    if scf_lcom_molecule.schoenflies_symbol() != 'c1':
        psi4.core.print_out("This SCF code must be run in c1 symmetry, switching\n")
        scf_lcom_molecule = scf_lcom_molecule.clone()
        scf_lcom_molecule.reset_point_group('c1')
        scf_lcom_molecule.update_geometry()

    new_wfn = psi4.core.Wavefunction.build(scf_lcom_molecule, psi4.core.get_global_option('BASIS'))

    scf_lcom_wfn = psi4.core.plugin('scf_lcom.so', new_wfn)
    psi4.set_variable('CURRENT ENERGY', scf_lcom_wfn.energy())

    if kwargs.get('ref_wfn', False):
        return (scf_lcom_wfn, scf_lcom_wfn.energy())
    else:
        return scf_lcom_wfn.energy()

# Integration with driver routines
psi4.driver.procedures['energy']['scf_lcom'] = run_scf_lcom


def exampleFN():
    # Your Python code goes here
    pass
