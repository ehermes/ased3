"""Grimme D3(BJ) correction scheme for DFT"""

import itertools
import numpy as np
from ase.units import *
from ase.calculators.d3params import *
from ase.calculators.general import Calculator
#from ase.calculators.d3ef import d3ef
from d3ef import d3ef

class D3(Calculator):
    """D3(BJ) correction of Grimme et al"""
    def __init__(self, bj=True, xc='pbe', rmax=30., calculator=None):

        Calculator.__init__(self)

        self.bj = bj
        if self.bj:
            self.damp = dampbj
        else:
            self.damp = damp
        # Atomic parameters
        self.rcov = None
        self.r2r4 = None
        self.r0 = None
        # Scaling parameters
        self.k1 = k1
        self.k2 = k2
        self.k3 = k3
        # Damping parameters
        self.rs6, self.s18, self.rs18, self.s6 = self.damp[xc.lower()]
        if self.bj:
            self.a1 = self.rs6
            self.a2 = self.rs18 * Bohr
        self.rs8 = self.rs18
        self.rs10 = self.rs18
        self.alp = alp
        self.alp6 = alp
        self.alp8 = alp+2
        # Cutoff distance
        self.rmax = rmax
        self.atoms = None
        self.calculator = calculator
        # Coordination number
        self.cn = None
        # Number of images
        self.nimages = np.array([0.,0.,0.],dtype=np.int16)
        # C6, C8, C9 parameters
        self.c6 = None
        self.c8 = None
        self.c9 = None
        # List of all atoms in periodic cells for easier iteration
        self.allatoms = []

    def set_atoms(self, atoms):
        self.atoms = atoms
        self.energy = 0.
        self.forces = None

    def updateparams(self, atoms=None):
        if atoms == None:
            atoms = self.atoms
        if self.rcov == None:
            self.rcov = np.zeros(len(self.atoms))
            for i, atom in enumerate(self.atoms):
                self.rcov[i] = rcov[atom.number - 1]
        if self.r2r4 == None:
            self.r2r4 = np.zeros(len(self.atoms))
            for i, atom in enumerate(self.atoms):
                self.r2r4[i] = r2r4[atom.number - 1]
        if self.r0 == None:
            self.r0 = np.zeros((len(self.atoms),len(self.atoms)))
            for i, atomi in enumerate(self.atoms):
                for j, atomj in enumerate(self.atoms):
                    self.r0[i,j] = r0ab[atomi.number - 1,atomj.number - 1]
        # BJ damping stuff
        if self.bj:
            dmp = self.a1 * np.sqrt(3.0 * np.outer(self.r2r4,self.r2r4)) + self.a2
            self.dmp6 = dmp**6
            self.dmp8 = dmp**8
        self.cn = np.zeros(len(atoms))
        cell = self.atoms.get_cell()
        # Calculate unit cell surface normal vectors
        surf = np.cross(np.roll(cell,2,axis=0),np.roll(cell,1,axis=0))
        # Dot unit vector into corresponding surface normal to find
        # "length" of each vector in x, y, and z
        vert = (np.diag(np.dot(cell,surf.T)) 
                / np.array([np.linalg.norm(subvec) for subvec in surf]))
        # Find the minimum number of images in each direction that contain
        # all atoms within a cutoff radius, so long as the system is periodic
        # in that direction
        self.nimages = np.int16(np.ceil(np.abs(self.rmax/vert))) * self.atoms.pbc
        # Iterate through all images
        for (i,j,k) in itertools.product(
                *[xrange(-self.nimages[ii],self.nimages[ii]+1) for ii in xrange(3)]):
            # Get translation vector for atom b
            tvec = np.dot(np.array([i,j,k]),cell)
            # Iterate through all pairs of atoms in unit cell
            added = np.zeros(len(atoms))
            for a, atoma in enumerate(atoms):
                for b, atomb in enumerate(atoms):
                    # Don't calculate interaction of an atom with itself in
                    # the central image
                    if (a == b) and (i == j == k == 0):
                        continue
                    # Get distance
                    rab = np.linalg.norm(atomb.position - atoma.position + tvec)
                    # If it's within the cutoff, calculate the contribution
                    # to CN^A *only*
                    if rab < self.rmax:
                        if not added[b]:
                            #print "Adding:", b, (i,j,k), atomb.position + tvec
                            self.allatoms.append([b,(i,j,k),atomb.position + tvec])
                            added[b] += 1
                        rcovab = self.rcov[a] + self.rcov[b]
                        # CN^A = sum(1/(1 + e^(-k1(k2(RcovA + RcovB)/rAB)-1)))
                        self.cn[a] += 1./(1. + np.exp(-self.k1 * (rcovab / rab - 1.)))
        # C6 terms are set individually for each pair, as D3 does not use simple 
        # combining rules for crossterms
        self.c6 = np.zeros((len(atoms),len(atoms)))
        for a, atoma in enumerate(atoms):
            for b, atomb in enumerate(atoms):
                na = atoma.number - 1
                nb = atomb.number - 1
                # Lij weights the c6 parameters listed in c6ab
                lij = np.zeros((numcn[na],numcn[nb]))
                c6abij = np.zeros((numcn[na],numcn[nb]))
                for i in xrange(numcn[na]):
                    for j in xrange(numcn[nb]):
                        c6abij[i,j] = c6ab[na,i,nb,j]
                        # Lij = e^(-k3*[(CNA - CNAi)^2 + (CNB - CNBj)^2])
                        lij[i,j] = np.exp(-k3*((self.cn[a] - cn[na,i])**2
                            + (self.cn[b] - cn[nb,j])**2))
                self.c6[a,b] = np.average(c6abij,weights=lij)
        self.c8 = 3.0 * self.c6 * np.outer(self.r2r4,self.r2r4)
        self.c9 = np.zeros((len(atoms),len(atoms),len(atoms)))
        for (a,b,c) in itertools.product(*[xrange(len(atoms)) for i in xrange(3)]):
            self.c9[a,b,c] = -np.sqrt(self.c6[a,b] * self.c6[b,c] * self.c6[c,a] / Hartree)

    def E2E3(self, atoms=None):
        if atoms == None:
            atoms = self.atoms
        e2 = 0.
        e3 = 0.
        for a, atoma in enumerate(atoms):
            for (bnum, (b, imageb, xyzb)) in enumerate(self.allatoms):
                ecalc = True
                if (imageb == (0,0,0)):
                    if (b == a):
                        continue
                    elif (b < a):
                        ecalc = False
                xyzab = xyzb - atoma.position
                # We need rab^2, so this is more computationally efficient
                # than using np.linalg.norm
                rab2 = np.dot(xyzab,xyzab)
                rab = np.sqrt(rab2)
                if self.bj:
                    dedc6 = -1./(rab**6 + self.dmp6[a,b])
                    dedc8 = -1./(rab**8 + self.dmp8[a,b])
                else:
                    dedc6 = -1./((rab**6)*(1. + 6.*(self.rs6*self.r0[a,b]/rab)**self.alp6))
                    dedc8 = -1./((rab**8)*(1. + 6.*(self.rs8*self.r0[a,b]/rab)**self.alp8))
                if ecalc:
                    e2 += self.s6 * self.c6[a,b] * dedc6 + self.s18 * self.c8[a,b] * dedc8
                for c, imagec, xyzc in self.allatoms[bnum+1:]:
                    if (imagec == (0,0,0)):
                        if (c == a):
                            continue
                        elif (c <= a):
                            ecalc = False
                    xyzac = xyzc - atoma.position
                    rac2 = np.dot(xyzac,xyzac)
                    rac = np.sqrt(rac2)
                    xyzbc = xyzc - xyzb
                    rbc2 = np.dot(xyzbc,xyzbc)
                    rbc = np.sqrt(rbc2)
                    cosalpha = (rab2 + rac2 - rbc2) / (2 * rab * rac)
                    cosbeta  = (rab2 + rbc2 - rac2) / (2 * rab * rbc)
                    cosgamma = (rac2 + rbc2 - rab2) / (2 * rac * rbc)
                    rav = 4./3. * (self.r0[a,b] * self.r0[b,c] * self.r0[a,c] 
                            / (rab * rbc * rac))**(1./3.)
                    dedc9 = -(3. * cosalpha * cosbeta * cosgamma + 1) / (rab * rac * rbc *
                            rab2 * rac2 * rbc2 * (1. + 6. * rav ** self.alp8))
                    if ecalc:
                        e3 += self.c9[a,b,c] * dedc9
        return e2 + e3

    def update(self, atoms=None):
        if atoms is None:
            atoms = self.calculator.get_atoms()
#        if (self.atoms and self.atoms.get_positions() == atoms.get_positions()).all():
#            return
        if self.calculator:
            self.energy = self.calculator.get_potential_energy(atoms)
            self.forces = self.calculator.get_forces(atoms)
        else:
            self.energy = 0.
            self.forces = np.zeros((len(atoms),3))
        self.atoms = atoms.copy()
        if self.cn == None:
            self.updateparams(atoms)
#        self.energy += self.E2E3(atoms)
        self.energy += d3ef.e2e3(
                image=np.array([elem[1] for elem in self.allatoms],dtype=np.int16),
                atomindex=np.array([elem[0] for elem in self.allatoms],dtype=np.int16) + 1,
                xyz=self.atoms.get_positions(),
                xyzall=np.array([elem[2] for elem in self.allatoms]),
                dmp6=self.dmp6,
                dmp8=self.dmp8,
                r0=self.r0,
                rmax=self.rmax,
                c6=self.c6,
                c8=self.c8,
                c9=self.c9,
                s6=self.s6,
                s18=self.s18,
                rs6=self.rs6,
                rs8=self.rs8,
                alp6=self.alp6,
                alp8=self.alp8,
                bj=self.bj,
                )

    def get_potential_energy(self, atoms):
        self.update(atoms)
        return self.energy
