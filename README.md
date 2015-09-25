### ASED3

This project implements a Grimme D3 Calculator for ASE.

Usage:

```
from ased3 import D3

calc = D3(
    xc='PBE',       # Which exchange-correlation functional to use
    bj=True,        # Whether to use Becke-Johnson damping, instead of zero-damping
    threebody=True, # Whether to include the Axilrod-Teller-Muto three-body terms
    calculator=dft, # A DFT calculator instance, e.g. VASP or GPAW. Optional.
    )

my_molecule.set_calculator(calc)

energy = my_molecule.get_potential_energy()
forces = my_molecule.get_forces()
stress = my_molecule.get_stress()
```
