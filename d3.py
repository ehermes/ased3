"""Grimme D3(BJ) correction scheme for DFT"""

import itertools
import numpy as np
from ase.units import *
from ase.calculators.d3params import k1, k2, k3, alp, damp, dampbj, numcn, cn, rcov, r2r4, r0ab, c6ab
from ase.calculators.general import Calculator
from ase.calculators.d3ef import d3ef

class D3(Calculator):
    """D3(BJ) correction of Grimme et al"""
    def __init__(self, bj=True, xc='pbe', rcut=95.*Bohr, rcutcn=40.*Bohr, calculator=None):

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
        self.dmp6 = None
        self.dmp8 = None

        # Cutoff distances
        self.rcut = rcut
        if self.rcut < rcutcn:
            self.rcutcn = self.rcut
        else:
            self.rcutcn = rcutcn

        self.atoms = None
        self.calculator = calculator

        # Coordination number
        self.cn = None
        self.dcn = None

        # Number of images
        self.nimages = np.array([0,0,0],dtype=np.int16)

        # C6, C8, C9 parameters, and their gradient
        self.c6 = None
        self.dc6 = None

        self.c8 = None
        self.dc8 = None

        self.c9 = None
        self.dc9 = None

        # List of all atoms in periodic cells for easier iteration
        self.allatomindex = None
        self.allatomimage = None
        self.allatomxyz = None

    def set_atoms(self, atoms):
        self.atoms = atoms
        self.energy = 0.
        self.forces = None

    def calccn(self, atoms):
        cell = atoms.get_cell()
        allatomindex = []
        allatomimage = []
        allatomxyz = []

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
                    # to CN^A *only* (not CN^B, we do that later)
                    if rab < self.rcut:
                        if not added[b]:
                            allatomindex.append(b)
                            allatomimage.append([i,j,k])
                            allatomxyz.append(atomb.position + tvec)
                            added[b] += 1

                        if rab < self.rcutcn:
                            rcovab = self.rcov[a] + self.rcov[b]
                            # CN^A = sum(1/(1 + e^(-k1(k2(RcovA + RcovB)/rAB)-1)))
                            self.cn[a] += 1./(1. + np.exp(-self.k1 * (rcovab / rab - 1.)))
        self.allatomindex = np.array(allatomindex,dtype=np.int16)
        self.allatomimage = np.array(allatomimage,dtype=np.int16)
        self.allatomxyz = np.array(allatomxyz,dtype=np.float64)


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
        self.dmp6 = np.zeros((len(atoms),len(atoms)))
        self.dmp8 = np.zeros((len(atoms),len(atoms)))
        if self.bj:
            dmp = self.a1 * np.sqrt(3.0 * np.outer(self.r2r4,self.r2r4)) + self.a2
            self.dmp6 = dmp**6
            self.dmp8 = dmp**8

        self.cn = np.zeros(len(atoms))
        self.dcn = np.zeros((len(atoms),len(atoms),3))
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
        self.nimages = np.int16(np.ceil(np.abs(self.rcut/vert))) * self.atoms.pbc

        # Calculate coordination number and all image atom positions with fortran code
        nallatoms = np.prod(2*self.nimages+2) * len(self.atoms)
        imageall, indexall, xyzall, self.cn, self.dcn = d3ef.cncalc(
                nallatoms=nallatoms,
                imagelist=self.nimages,
                k1=self.k1,
                cell=cell,
                xyz=self.atoms.get_positions(),
                rcut=self.rcut,
                rcutcn=self.rcutcn,
                rcov=self.rcov,
                )
        indices = []
        for i, index in enumerate(indexall):
            if (index >= 0):
                indices.append(i)
        self.allatomindex = np.zeros(len(indices),dtype=np.int16)
        self.allatomimage = np.zeros((len(indices),3),dtype=np.int16)
        self.allatomxyz = np.zeros((len(indices),3),dtype=np.float64)
        for i, j in enumerate(indices):
            self.allatomindex[i] = indexall[j]
            self.allatomimage[i] = imageall[j]
            self.allatomxyz[i] = xyzall[j]
#        for i, index in enumerate(indexall):
#            if (index >= 0):
#                self.allatoms.append([index, tuple(imageall[i]), xyzall[i]])

        ## Calculate coordination number and all image atom positions with python code
        #self.calccn(self.atoms)
        # C6 terms are set individually for each pair, as D3 does not use simple 
        # combining rules for crossterms
        self.c6 = np.zeros((len(atoms),len(atoms)))
        self.dc6 = np.zeros((len(atoms),len(atoms),len(atoms),3))

        for a, atoma in enumerate(atoms):
            na = atoma.number - 1
            ncna = numcn[na]
            for b, atomb in enumerate(atoms):
                nb = atomb.number - 1
                ncnb = numcn[nb]

                # Lij weights the c6 parameters listed in c6ab
                c6abij = c6ab[na,:ncna,nb,:ncnb]

                lij = np.exp(-k3*(((self.cn[a] - cn[na,:ncna])**2)[:,np.newaxis]
                    + (self.cn[b] - cn[nb,:ncnb])**2))

                dlij = -2 * lij[:,:,np.newaxis] * self.k3 \
                        * ((self.cn[a] - cn[na,:ncna])[:,np.newaxis,np.newaxis] \
                        * self.dcn[:,a][:,np.newaxis,np.newaxis] \
                        + (self.cn[b] - cn[nb,:ncnb])[:,np.newaxis] \
                        * self.dcn[:,b][:,np.newaxis,np.newaxis])
                        
                self.c6[a,b] = np.average(c6abij,weights=lij)
                self.dc6[:,a,b] = (-self.c6[a,b] * np.sum(dlij,axis=(1,2)) \
                        + np.sum(c6abij[:,:,np.newaxis] * dlij,axis=(1,2))) \
                        / np.sum(lij)
        self.c8 = 3.0 * self.c6 * np.outer(self.r2r4,self.r2r4)
        self.dc8 = 3.0 * self.dc6 * np.outer(self.r2r4,self.r2r4)[:,:,np.newaxis]

        self.c9 = np.zeros((len(atoms),len(atoms),len(atoms)))
        self.dc9 = np.zeros((len(atoms),len(atoms),len(atoms),len(atoms),3))
        self.c9 = -np.sqrt(self.c6[:,:,np.newaxis] * self.c6 * self.c6[:,np.newaxis,:] 
                / Hartree)
        self.dc9 = (self.dc6[:,:,:,np.newaxis] 
                * self.c6[:,:,np.newaxis] 
                * self.c6[:,np.newaxis,:,np.newaxis]
                + self.c6[:,:,np.newaxis,np.newaxis] 
                * self.dc6[:,np.newaxis] 
                * self.c6[:,np.newaxis,:,np.newaxis]
                + self.c6[:,:,np.newaxis,np.newaxis] 
                * self.c6[:,:,np.newaxis] 
                * self.dc6[:,:,np.newaxis]) \
                        / (2 * Hartree * self.c9[:,:,:,np.newaxis])

    def E2E3(self, atoms=None):
        if atoms == None:
            atoms = self.atoms

        e6 = 0.
        f6 = np.zeros((natoms,3))

        e8 = 0.
        e8 = np.zeros((natoms,3))

        e9 = 0.
        e9 = np.zeros((natoms,3))

        for a, atoma in enumerate(atoms):
            for bnum, b in enumerate(self.allatomindex):
                imageb = self.allatomimage[bnum]
                xyzb = self.allatomxyz[bnum]
                ecalc = True
                # Don't calculate interaction if atom b == atom a (only central unit cell)
                if (imageb == (0,0,0)):
                    if (b == a):
                        continue

                xyzab = xyzb - atoma.position
                rab2 = np.dot(xyzab,xyzab)
                rab = np.sqrt(rab2)
                uxyzab = xyzab / rab

                if rab > self.rcut:
                    continue

                # Choose which type of damping to use
                if self.bj:
                    dedc6 = -1./(rab**6 + self.dmp6[a,b])
                    dfdc6 = -6 * dedc6**2 * rab**5 * uxyzab

                    dedc8 = -1./(rab**8 + self.dmp8[a,b])
                    dfdc8 = -8 * dedc8**2 * rab**7 * uxyzab
                else:
                    rav = (self.rs6 * self.r0[a,b] / rab)**self.alp6
                    damp6 = 1. / (1. + 6. * rav)
                    dedc6 = -damp6 / rab**6
                    dfdc6 = -6. * uxyzab * (1. - self.alp6 * rav + 6. * rav) \
                            / (rab**7 * damp6**2)

                    rav = (self.rs8 * self.r0[a,b]/rab)**self.alp8
                    damp8 = 1. / (1. + 6. * rav)
                    dedc8 = -damp8 / rab**8
                    dfdc8 = -uxyzab * (8. - 6. * self.alp8 * rav + 48. * rav) \
                            / (rab**9 * damp8**2)

                e6 += self.s6 * self.c6[a,b] * dedc6 
                f6[a] += self.s6 * self.c6[a,b] * dfdc6

                e8 += self.s18 * self.c8[a,b] * dedc8
                f8[a] += self.s18 * self.c8[a,b] * dedc8

                if (rab > self.rcutcn):
                    continue

                for cnum, c in enumerate(self.allatomindex):
                    imagec = self.allatomimage[cnum]
                    xyzc = self.allatomxyz[cnum]
                    # Don't calculate interaction if atom c == atom a, or atom c == atom b
                    if (imagec == (0,0,0)):
                        if (c == a):
                            continue

                    if (cnum == bnum):
                        continue

                    xyzac = xyzc - atoma.position
                    rac2 = np.dot(xyzac,xyzac)
                    rac = np.sqrt(rac2)

                    if (rac > self.rcutcn):
                        continue

                    xyzbc = xyzc - xyzb
                    rbc2 = np.dot(xyzbc,xyzbc)
                    rbc = np.sqrt(rbc2)

                    if (rbc > self.rcutcn):
                        continue

                    rab3 = rab * rab2
                    rac3 = rac * rac2
                    rbc3 = rbc * rbc2
                    xyzab2 = xyzab**2
                    xyzac2 = xyzac**2
                    xyzbc2 = xyzbc**2

                    cosalpha = (rab2 + rac2 - rbc2) / (2 * rab * rac)
                    cosbeta  = (rab2 + rbc2 - rac2) / (2 * rab * rbc)
                    cosgamma = (rac2 + rbc2 - rab2) / (2 * rac * rbc)

                    dcosalpha = cosalpha * (xyzac/rac2 + xyzab/rab2) \
                            - (xyzab + xyzac) / (rab*rac)
                    dcosbeta = (-self.xyz[a] * (rab*rbc*cosbeta + xyzab*xyzbc) \
                            - self.xyz[b] * (rab*rac*cosalpha - xyzab*xyzac) \
                            + self.xyz[c] * (rab2 - xyzab2))/(rab3*rbc)
                    dcosgamma = (-self.xyz[a] * (rac*rbc*cosgamma - xyzac*xyzbc) \
                            - self.xyz[c] * (rab*rac*cosalpha - xyzab*xyzac) \
                            + self.xyz[b] * (rac2 - xyzac2))/(rac3*rbc3)

                    angles = 3. * cosalpha * cosbeta * cosgamma + 1.
                    dangles = 3. * (cosalpha * cosbeta * dcosgamma \
                            + cosalpha * cosgamma * dcosbeta \
                            + cosbeta * cosgamma * dcosalpha)

                    rav = 4./3. * (self.r0[a,b] * self.r0[b,c] * self.r0[a,c] 
                            / (rab * rbc * rac))**(1./3.)
                    drav = (4./9.) * rav * (xyzab/rab2 + xyzac/rac2)

                    damp9 = 1./(1. + 6. * rav**self.alp8)
                    ddamp9 = -6. * rav**(self.alp8-1) * drav * damp9**2

                    r9 = 1. / (rab3 * rac3 * rbc3)
                    dr9 = 3. * r9 * (xyzab/rab2 + xyzac/rac2)

                    dedc9 = -angles * damp9 * r9
                    dfdc9 = -dangles * damp9 * r9 \
                            - angles * ddamp9 * r9 \
                            - angles * damp8 * dr9

                    e9 += self.c9[a,b,c] * dedc9
                    f9[a] += self.c9[a,b,c] * dfdc9
        # We double counted many interactions in the name of easier force calculations, so let's fix that here
        return (e6 + e8)/2. + e3/6.

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
        e, f = d3ef.e2e3(
                image=self.allatomimage,
                atomindex=self.allatomindex + 1,
                xyz=self.atoms.get_positions(),
                xyzall=self.allatomxyz,
                dmp6=self.dmp6,
                dmp8=self.dmp8,
                r0=self.r0,
                rcut=self.rcut,
                rcutcn=self.rcutcn,
                c6=self.c6,
                dc6=self.dc6,
                c8=self.c8,
                dc8=self.dc8,
                c9=self.c9,
                dc9=self.dc9,
                s6=self.s6,
                s18=self.s18,
                rs6=self.rs6,
                rs8=self.rs8,
                alp6=self.alp6,
                alp8=self.alp8,
                bj=self.bj,
                )
        self.energy += e
        self.forces += f

    def get_potential_energy(self, atoms):
        self.update(atoms)
        return self.energy
