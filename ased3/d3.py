from __future__ import division, print_function

import numpy as np

from ase.units import Bohr
from ase.calculators.calculator import Calculator, all_changes
from ased3.d3params import damp, dampbj, alp
from ased3.d3_fort import d3_fort

d3_calc = d3_fort.d3_calc

class D3(Calculator):
    implemented_properties = ['energy', 'forces', 'stress']
    default_parameters = {'xc': 'pbe',
                          'bj': True,
                          'threebody': True,
                          'rcut': 95 * Bohr,
                          'rcutcn': 40 * Bohr,
                          'rs6': None,
                          's18': None,
                          'rs18': None,
                          's6': None,
                          'calculator': None}

    nolabel = True
    def __init__(self, **kwargs):
        Calculator.__init__(self, **kwargs)

    def calculate(self, atoms=None, properties=['energy'],
                   system_changes=all_changes):
        Calculator.calculate(self, atoms, properties, system_changes)

        xc = self.parameters.xc
        bj = self.parameters.bj
        threebody = self.parameters.threebody
        rcut = self.parameters.rcut
        rcutcn = self.parameters.rcutcn
        calculator = self.parameters.calculator

        if bj:
            rs6, s18, rs18, s6 = dampbj[xc]
        else:
            rs6, s18, rs18, s6 = damp[xc]

        if self.parameters.s6 is not None:
            s6 = self.parameters.s6
        if self.parameters.s18 is not None:
            s18 = self.parameters.s18
        if self.parameters.rs6 is not None:
            rs6 = self.parameters.rs6
        if self.parameters.rs18 is not None:
            rs18 = self.parameters.rs18

        energy, forces, stress = d3_calc(
                self.atoms.get_atomic_numbers(),
                self.atoms.get_cell(),
                self.atoms.get_positions().T,
                rcut=rcut,
                rcutcn=rcutcn,
                s6=s6,
                s18=s18,
                rs6=rs6,
                rs18=rs18,
                alp6=alp,
                alp8=alp + 2,
                pbc=self.atoms.get_pbc(),
                bj=bj,
                threebody=threebody)
        
        self.results['energy'] = energy
        self.results['forces'] = forces.T
        self.results['stress'] = stress.flat[[0, 4, 8, 5, 2, 1]]

        if calculator is not None:
            calculate(self.atoms)
            self.results['energy'] += calculator.results['energy']
            self.results['forces'] += calculator.results['forces']
            self.results['stress'] += calculator.results['stress']
