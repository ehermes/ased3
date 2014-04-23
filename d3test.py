from ase.lattice.surface import fcc111
from ase.data.s22 import create_s22_system
from ase.calculators.d3 import D3
from ase.calculators.test import numeric_forces

fortran = True

try:
    from ase.calculators.d3ef import d3ef
except ImportError:
    fortran = False

def testforce(system, fortran=False, bj=False):
    system.set_calculator(D3(fortran=fortran, bj=bj))
    f1 = system.get_forces()
    f2 = numeric_forces(system, d=0.0001)
    
    assert abs(f1 - f2).max() < 1e-6

waterdim = create_s22_system('Water_dimer')
slab = fcc111('Au', (1, 1, 4), vacuum=7.5)

testforce(waterdim, fortran=False, bj=False)
#testforce(slab, fortran=False, bj=False) # This is too slow!

testforce(waterdim, fortran=False, bj=True)
#testforce(slab, fortran=False, bj=True) # This is too slow!

if fortran:
    testforce(waterdim, fortran=True, bj=False)
    testforce(slab, fortran=True, bj=False)
    
    testforce(waterdim, fortran=True, bj=True)
    testforce(slab, fortran=True, bj=True)
