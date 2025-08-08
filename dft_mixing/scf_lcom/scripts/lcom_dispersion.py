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
from psi4.driver.procrouting.empirical_dispersion import EmpiricalDispersion
import logging
logger = logging.getLogger(__name__)
logging.disable(logging.WARNING)
class LcomDispersion(EmpiricalDispersion):
    """
    A type of dispersion that implements a linear combination of parent empirical dispersion algorithms.
    
    Parameters
    ----------
    `parents_coefs`:
    
        A list of tuples of the form of `(EmpiricalDispersion, dispersion_coefficient)`.
        Each tuple represents a parent dispersion scheme, which is
        linearly combined using the coefficient contained in the last float.
    """
    def __init__(self, parents_coefs: list[tuple[EmpiricalDispersion, float]]):
        self.parents_coefs = parents_coefs
        for (parent, coef) in parents_coefs:
            if not (isinstance(coef, float) or isinstance(coef, int)):
                raise RuntimeError("Invalid LCOM coefficient in dispersion: " + str(type(coef)))
            if not issubclass(type(parent), EmpiricalDispersion):
                raise RuntimeError("LCOM parent dispersion must be an instance of EmpiricalDispersion")

    def add_empirical_disp(self, disp: EmpiricalDispersion, coef: float) -> None:
        self.parents_coefs.append((disp, coef))

    def print_out(self):
        should_add_coef = len(self.parents_coefs) != 1 or self.parents_coefs[0][1] != 1.0
        for (p, c) in self.parents_coefs:
            # TODO: Make this less hack
            if should_add_coef: psi4.core.print_out("    With a coefficient of %.5f...\n" % c)
            p.print_out()
                
    def removeNlDisp(self):
        return LcomDispersion(
            list(filter(lambda p: p[0].engine != 'nl', self.parents_coefs))
        )

    def compute_energy(self, molecule: psi4.core.Molecule, wfn: psi4.core.Wavefunction = None) -> float:
        return sum((coef * dispersion.compute_energy(molecule, wfn)
                        for (dispersion, coef) in self.parents_coefs))

    def compute_gradient(self,
                         molecule: psi4.core.Molecule,
                         wfn: psi4.core.Wavefunction = None) -> psi4.core.Matrix:
        
        # Not sure if this is really correct under all dispersions, but
        # if the gradient is referring to first derrivatives, then
        # this should be a sum of all gradient computations

    
        raise RuntimeError("Not yet implemented")

    def compute_hessian(self,
                        molecule: psi4.core.Molecule,
                        wfn: psi4.core.Wavefunction = None) -> psi4.core.Matrix:
        raise RuntimeError("Not yet implemented")



