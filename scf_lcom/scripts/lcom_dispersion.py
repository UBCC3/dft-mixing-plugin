import psi4
from psi4.driver.procrouting.empirical_dispersion import EmpiricalDispersion

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

    def print_out(self):
        should_add_coef = len(self.parents_coefs) != 1 or self.parents_coefs[0][1] != 1.0
        for (p, c) in self.parents_coefs:
            # TODO: Make this less hack
            if should_add_coef: psi4.core.print_out("    With a coefficient of %.5f...\n" % c)
            p.print_out()
    
    
    def removeNlDisp(self):
        return LcomDispersion(
            *filter(lambda p: p[0].engine != 'nl', self.parents_coefs)
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



