"""Grimme D3(BJ) correction scheme for DFT"""

import itertools
import numpy as np
from ase.units import Hartree, Bohr
from ase.calculators.d3params import k1, k2, k3, alp, damp, dampbj, \
        numcn, cn, rcov, r2r4, r0ab, c6ab

usefortran = True
try:
    from ase.calculators.d3ef import d3ef
except ImportError:
    usefortran = False

class D3(object):
    """D3(BJ) correction of Grimme et al"""
    def __init__(self,
            bj=True,
            xc='pbe',
            rcut=95. * Bohr,
            rcutcn=40. * Bohr,
            calculator=None,
            fortran=usefortran):

        # Whether we use faster Fortran code, or slower native python
        self.fortran = fortran

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
        self.alp8 = alp + 2
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
        self.nimages = np.zeros(3, dtype=np.int16)

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

        # Energy, forces, and stress
        self.energy = None
        self.forces = None
        self.stress = None

    def calccn(self, atoms):
        """Calculates the coordination number of each atom in the system
        and its gradient."""

        cell = atoms.get_cell()
        allatomindex = []
        allatomimage = []
        allatomxyz = []

        self.cn = np.zeros(len(atoms))
        self.dcn = np.zeros((len(atoms), len(atoms), 3))

        # Iterate through all images
        for (i, j, k) in itertools.product(
                *[xrange(-self.nimages[ii], self.nimages[ii] + 1) \
                        for ii in xrange(3)]):
            # Get translation vector for atom b
            tvec = np.dot(np.array([i, j, k]), cell)

            # Iterate through all pairs of atoms in unit cell
            added = np.zeros(len(atoms))
            for a, atoma in enumerate(atoms):
                for b, atomb in enumerate(atoms):
                    # Don't calculate interaction of an atom with itself in
                    # the central image
                    if (a == b) and (i == j == k == 0):
                        continue
                    # Get distance
                    xyzab = atomb.position - atoma.position + tvec
                    rab2 = np.dot(xyzab, xyzab)
                    rab = np.sqrt(rab2)

                    # If it's within the cutoff, calculate the contribution
                    # to CN^A *only* (not CN^B, we do that later)
                    if rab < self.rcut:
                        if not added[b]:
                            allatomindex.append(b)
                            allatomimage.append([i, j, k])
                            allatomxyz.append(atomb.position + tvec)
                            added[b] += 1

                        if rab < self.rcutcn:
                            uxyzab = xyzab / rab
                            rcovab = self.rcov[a] + self.rcov[b]
                            cnexp = np.exp(-k1 * (rcovab / rab - 1.))
                            cnab = 1. / (1. + cnexp)
                            dcnab = cnexp * k1 * rcovab * cnab**2 \
                                    * uxyzab / rab2
                            self.cn[a] += cnab
                            self.dcn[a, b, :] = dcnab
                            self.dcn[a, a, :] += dcnab
        self.allatomindex = np.array(allatomindex, dtype=np.int16)
        self.allatomimage = np.array(allatomimage, dtype=np.int16)
        self.allatomxyz = np.array(allatomxyz, dtype=np.float64)

    def updateparams(self, atoms=None):
        """Calculate dispersion parameters and their gradients"""

        if atoms is None:
            atoms = self.atoms

        self.rcov = np.zeros(len(self.atoms))
        self.r2r4 = np.zeros(len(self.atoms))
        self.r0 = np.zeros((len(self.atoms), len(self.atoms)))

        for i, atomi in enumerate(self.atoms):
            self.rcov[i] = rcov[atomi.number - 1]
            self.r2r4[i] = r2r4[atomi.number - 1]
            for j, atomj in enumerate(self.atoms):
                self.r0[i, j] = r0ab[atomi.number - 1, atomj.number - 1]

        # BJ damping stuff
        self.dmp6 = np.zeros((len(atoms), len(atoms)))
        self.dmp8 = np.zeros((len(atoms), len(atoms)))
        if self.bj:
            dmp = self.a1 * np.sqrt(3.0 \
                    * np.outer(self.r2r4, self.r2r4)) + self.a2
            self.dmp6 = dmp**6
            self.dmp8 = dmp**8

        self.cn = np.zeros(len(atoms))
        self.dcn = np.zeros((len(atoms), len(atoms), 3))
        cell = self.atoms.get_cell()

        # Calculate unit cell surface normal vectors
        surf = np.cross(np.roll(cell, 2, axis=0), \
                np.roll(cell, 1, axis=0))

        # Dot unit vector into corresponding surface normal to find
        # 'length' of each vector in x, y, and z
        vert = (np.diag(np.dot(cell, surf.T))
                / np.array([np.linalg.norm(subvec) for subvec in surf]))

        # Find the minimum number of images in each direction that contain
        # all atoms within a cutoff radius, so long as the system is periodic
        # in that direction
        self.nimages = np.int16(np.ceil(np.abs(self.rcut / vert))) \
                * self.atoms.pbc

        # Calculate coordination number and all image atom positions
        # with fortran code
        if self.fortran:
            nallatoms = np.prod(2 * self.nimages + 1) * len(self.atoms)
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
                if index >= 0:
                    indices.append(i)
            self.allatomindex = np.zeros(len(indices), dtype=np.int16)
            self.allatomimage = np.zeros((len(indices), 3), dtype=np.int16)
            self.allatomxyz = np.zeros((len(indices), 3), dtype=np.float64)
            for i, j in enumerate(indices):
                self.allatomindex[i] = indexall[j]
                self.allatomimage[i] = imageall[j]
                self.allatomxyz[i] = xyzall[j]
        else:
            # Calculate coordination number and all image atom
            # positions with python code
            self.calccn(self.atoms)
        # C6 terms are set individually for each pair, as D3 does not
        # use simple combining rules for crossterms
        self.c6 = np.zeros((len(atoms), len(atoms)))
        self.dc6 = np.zeros((len(atoms), len(atoms), len(atoms), 3))

        for a, atoma in enumerate(atoms):
            na = atoma.number - 1
            ncna = numcn[na]
            for b, atomb in enumerate(atoms):
                nb = atomb.number - 1
                ncnb = numcn[nb]

                # Lij weights the c6 parameters listed in c6ab
                c6abij = c6ab[na, :ncna, nb, :ncnb]

                lij = np.exp(-k3 * (((self.cn[a] - cn[na, :ncna])**2)[:, \
                        np.newaxis] + (self.cn[b] - cn[nb, :ncnb])**2))

                dlij = -2 * lij[:, :, np.newaxis] * self.k3 * ((self.cn[a] \
                        - cn[na, :ncna])[:, np.newaxis, np.newaxis] \
                        * self.dcn[:, a][:, np.newaxis, np.newaxis] \
                        + (self.cn[b] - cn[nb, :ncnb])[:, np.newaxis] \
                        * self.dcn[:, b][:, np.newaxis, np.newaxis])

                self.c6[a, b] = np.average(c6abij, weights=lij)
                self.dc6[:, a, b] = (-self.c6[a, b] * np.sum(dlij, \
                        axis=(1, 2)) + np.sum(c6abij[:, :, np.newaxis] \
                        * dlij, axis=(1, 2))) / np.sum(lij)
        self.c8 = 3.0 * self.c6 * np.outer(self.r2r4, self.r2r4)
        self.dc8 = 3.0 * self.dc6 * np.outer(self.r2r4, \
                self.r2r4)[:, :, np.newaxis]

        self.c9 = np.zeros((len(atoms), len(atoms), len(atoms)))
        self.dc9 = np.zeros((len(atoms), len(atoms), len(atoms), \
                len(atoms), 3))
        self.c9 = -np.sqrt(self.c6[:, :, np.newaxis] * self.c6 * \
                self.c6[:, np.newaxis, :] / Hartree)
        self.dc9 = (self.dc6[:, :, :, np.newaxis] \
                * self.c6[:, :, np.newaxis] \
                * self.c6[:, np.newaxis, :, np.newaxis] \
                + self.c6[:, :, np.newaxis, np.newaxis] \
                * self.dc6[:, np.newaxis] \
                * self.c6[:, np.newaxis, :, np.newaxis] \
                + self.c6[:, :, np.newaxis, np.newaxis] \
                * self.c6[:, :, np.newaxis] \
                * self.dc6[:, :, np.newaxis]) \
                / (2 * Hartree * self.c9[:, :, :, np.newaxis])

    def efcalc(self, atoms=None):
        """Calculate the system energy and atomic forces"""

        if atoms is None:
            atoms = self.atoms
        natoms = len(atoms)

        e6 = 0.
        f6 = np.zeros((natoms, 3))

        e8 = 0.
        f8 = np.zeros((natoms, 3))

        eabc = 0.
        fabc = np.zeros((natoms, 3))

        # Atom a loops over atoms in the central unit cell
        for a, atoma in enumerate(atoms):
            # Atom b loops over all image atoms
            for bnum, b in enumerate(self.allatomindex):
                imageb = self.allatomimage[bnum]
                xyzb = self.allatomxyz[bnum]

                # Self-interaction scaling parameter
                sself = 1.

                if b < a:
                    continue

                # Don't calculate the interaction if atom b <= atom a
                # if atom b is in the central unit cell
                if b == a:
                    if np.all(imageb == np.zeros(3, dtype=np.int16)):
                        continue
                    else:
                        sself = 1. / 2.

                # Set some variables for interaction between a and b
                xyzab = xyzb - atoma.position
                rab2 = np.dot(xyzab, xyzab)
                rab = np.sqrt(rab2)

                # Don't calculate the interaction if rab > rcut
                if rab > self.rcut:
                    continue

                # Choose which type of damping to use, and calculate
                # the damping factor
                # Becke-Johnson damping
                if self.bj:
                    dedc6 = -1. / (rab**6 + self.dmp6[a, b])
                    dfdc6 = -6. * dedc6**2 * rab**4 * xyzab

                    dedc8 = -1. / (rab**8 + self.dmp8[a, b])
                    dfdc8 = -8. * dedc8**2 * rab**6 * xyzab
                # 'Zero-damping'
                else:
                    rav = (self.rs6 * self.r0[a, b] / rab)**self.alp6
                    drav = -xyzab * self.alp6 * rav / rab2
                    damp6 = 1. / (1. + 6. * rav)
                    ddamp6 = -6. * damp6**2 * drav
                    dedc6 = -damp6 / rab**6
                    dfdc6 = 6. * xyzab * dedc6 / rab2 + ddamp6 / rab**6

                    rav = (self.rs8 * self.r0[a, b] / rab)**self.alp8
                    drav = -xyzab * self.alp8 * rav / rab2
                    damp8 = 1. / (1. + 6. * rav)
                    ddamp8 = -6. * damp8**2 * drav
                    dedc8 = -damp8 / rab**8
                    dfdc8 = 8. * xyzab * dedc8 / rab2 + ddamp8 / rab**8

                # C6 energy and force contributions
                e6 += self.s6 * self.c6[a, b] * dedc6 * sself
                if b != a:
                    f6[a] -= self.s6 * self.c6[a, b] * dfdc6
                    f6[b] += self.s6 * self.c6[a, b] * dfdc6
                f6 -= self.s6 * self.dc6[:, a, b] * dedc6 * sself

                e8 += self.s18 * self.c8[a, b] * dedc8 * sself
                if b != a:
                    f8[a] -= self.s18 * self.c8[a, b] * dfdc8
                    f8[b] += self.s18 * self.c8[a, b] * dfdc8
                f8 -= self.s18 * self.dc8[:, a, b] * dedc8 * sself

                # Don't calculate 3-body term if rab > rcutcn
                if rab > self.rcutcn:
                    continue

                for cnum, c in enumerate(self.allatomindex):
                    imagec = self.allatomimage[cnum]
                    xyzc = self.allatomxyz[cnum]

                    # 3-body self interaction scaling term. Can be
                    # either 1, 1/2, or 1/6.
                    sself = 1.

                    # Don't calculate interaction if c < a
                    if c < a:
                        continue

                    # Don't calculate interaction if c is in the
                    # central unit cell and
                    # a == c
                    if c == a:
                        if np.all(imagec == np.zeros(3, dtype=np.int16)):
                            continue

                    # Don't calculate interaction if c < b
                    if c < b:
                        continue

                    # Don't calculate interaction if c and b are in the
                    # same unit cell and
                    # c == b
                    if c == b:
                        if np.all(imagec == imageb):
                            continue

                    # Figure out if we're calculating a self energy.
                    # If exactly two of a, b, and c are the same
                    # index, then the system is doubly degenerate
                    # and the energy must be divided by 2. If all three
                    # indices are the same, then the system is 6-fold
                    # degenerate, and the energy must be divided by 6.
                    if (a == c) or (a == b) or (b == c):
                        if (a == b) and (b == c):
                            sself = 1. / 6.
                        else:
                            sself = 1. / 2.

                    # Set soem variables for a interacting with c
                    xyzac = xyzc - atoma.position
                    rac2 = np.dot(xyzac, xyzac)
                    rac = np.sqrt(rac2)

                    # Don't calculate the interaction if rac > rcutcn
                    if rac > self.rcutcn:
                        continue

                    # Set some variables for b interacting with c
                    xyzbc = xyzc - xyzb
                    rbc2 = np.dot(xyzbc, xyzbc)
                    rbc = np.sqrt(rbc2)

                    # Don't calculate the interaction if rbc > rcutcn
                    if rbc > self.rcutcn:
                        continue

                    # Pre-calculate some powers of rab that we're going
                    # to need later
                    rab3 = rab * rab2
                    rac3 = rac * rac2
                    rbc3 = rbc * rbc2

                    # Use dot product instead of law of cosines to
                    # calculate angle terms to be consistent with
                    # derivatives below (it doesn't actually matter)
                    cosalpha = np.dot(xyzab, xyzac) / (rab * rac)
                    cosbeta = -np.dot(xyzab, xyzbc) / (rab * rbc)
                    cosgamma = np.dot(xyzac, xyzbc) / (rac * rbc)

                    # Gradient of cosalpha, cosbeta, cosgamma.  Very
                    # complicated... Figured this all out using
                    # Mathematica and defining cosalpha, cosbeta,
                    # cosgamma as above.
                    dacosalpha = cosalpha * (xyzac / rac2 + xyzab / rab2) \
                            - (xyzac + xyzab) / (rab * rac)
                    dacosbeta = xyzbc / (rab * rbc) + xyzab * cosbeta / rab2
                    dacosgamma = -xyzbc / (rac * rbc) + xyzac * cosgamma \
                            / rac2

                    dbcosalpha = xyzac / (rab * rac) - xyzab * cosalpha / rab2
                    dbcosbeta = cosbeta * (xyzbc / rbc2 - xyzab / rab2) \
                            + (xyzab - xyzbc) / (rab * rbc)
                    dbcosgamma = -xyzac / (rac * rbc) + xyzbc * cosgamma \
                            / rbc2

                    dccosalpha = xyzab / (rab * rac) - xyzac * cosalpha / rac2
                    dccosbeta = -xyzab / (rab * rbc) - xyzbc * cosbeta / rbc2
                    dccosgamma = -cosgamma * (xyzac / rac2 + xyzbc / rbc2) \
                            + (xyzac + xyzbc) / (rac * rbc)

                    # I have no idea what 'rav' stands for, but that's
                    # what Grimme called this variable.  Cube root of
                    # the product of the ratios of r0ab/rab, times 4/3
                    # for some reason. I don't know.
                    rav = 4. / 3. * (self.r0[a, b] * self.r0[b, c] \
                            * self.r0[a, c] / (rab * rbc * rac))**(1. / 3.)
                    darav = (rav / 3.) * (xyzab / rab2 + xyzac / rac2)
                    dbrav = (rav / 3.) * (-xyzab / rab2 + xyzbc / rbc2)
                    dcrav = -(rav / 3.) * (xyzac / rac2 + xyzbc / rbc2)

                    # Three-body term *always* uses 'zero' damping,
                    # even if we are using the BJ version of DFT-D3
                    damp9 = 1. / (1. + 6. * rav**self.alp8)
                    ddamp9 = -6. * self.alp8 * rav**(self.alp8 - 1) * damp9**2

                    # Three-body depends on 'average' r^9
                    r9 = 1. / (rab3 * rac3 * rbc3)
                    dar9 = 3. * r9 * (xyzab / rab2 + xyzac / rac2)
                    dbr9 = 3. * r9 * (-xyzab / rab2 + xyzbc / rbc2)
                    dcr9 = -3. * r9 * (xyzac / rac2 + xyzbc / rbc2)

                    # The derivatives change if two or more of the
                    # atoms are the same in different unit cells.
                    # These are not obvious, but mathematica/sympy/etc
                    # will confirm that it is correct.
                    if (a == b) and (b != c) and (a != c):
                        dacosalpha = xyzac * cosalpha / rac2 \
                                - xyzab / (rab * rac)
                        dacosbeta = xyzbc * cosbeta / rbc2 \
                                + xyzab / (rab * rac)
                        dacosgamma = -dccosgamma

                        dbcosalpha = dacosalpha
                        dbcosbeta = dacosbeta
                        dbcosgamma = dacosgamma

                        darav = -dcrav
                        dbrav = -dcrav

                        dar9 = -dcr9
                        dbr9 = -dcr9
                    elif (a != b) and (b == c) and (a != c):
                        dbcosalpha = -dacosalpha
                        dbcosbeta = -xyzab * cosbeta / rab2 \
                                - xyzbc / (rab * rbc)
                        dbcosgamma = -xyzac * cosgamma / rac2 \
                                + xyzbc / (rac * rbc)

                        dccosalpha = dbcosalpha
                        dccosbeta = dbcosbeta
                        dccosgamma = dbcosgamma

                        dbrav = -darav
                        dcrav = -darav

                        dbr9 = -dar9
                        dcr9 = -dcr9
                    elif (a != b) and (b != c) and (a == c):
                        dacosalpha = xyzab * cosalpha / rab2 \
                                - xyzac / (rac * rab)
                        dacosbeta = -dbcosbeta
                        dacosgamma = -xyzbc * cosgamma / rbc2 \
                                + xyzac / (rbc * rac)

                        dccosalpha = dacosalpha
                        dccosbeta = dacosbeta
                        dccosgamma = dacosgamma

                        darav = -dbrav
                        dcrav = -dbrav

                        dar9 = -dbr9
                        dcr9 = -dbr9
                    elif (a == b) and (b == c) and (a == c):
                        dacosalpha = 0.
                        dacosbeta = 0.
                        dacosgamma = 0.

                        dbcosalpha = 0.
                        dbcosbeta = 0.
                        dbcosgamma = 0.

                        dccosalpha = 0.
                        dccosbeta = 0.
                        dccosgamma = 0.

                        darav = 0.
                        dbrav = 0.
                        dcrav = 0.

                        dar9 = 0.
                        dbr9 = 0.
                        dcr9 = 0.

                    # Angle terms of the three body energy, and its gradient
                    angles = 3. * cosalpha * cosbeta * cosgamma + 1.
                    daangles = 3. * (dacosalpha * cosbeta * cosgamma
                            + cosalpha * dacosbeta * cosgamma
                            + cosalpha * cosbeta * dacosgamma)
                    dbangles = 3. * (dbcosalpha * cosbeta * cosgamma
                            + cosalpha * dbcosbeta * cosgamma
                            + cosalpha * cosbeta * dbcosgamma)
                    dcangles = 3. * (dccosalpha * cosbeta * cosgamma
                            + cosalpha * dccosbeta * cosgamma
                            + cosalpha * cosbeta * dccosgamma)

                    # Damping derivatives
                    dadamp9 = ddamp9 * darav
                    dbdamp9 = ddamp9 * dbrav
                    dcdamp9 = ddamp9 * dcrav

                    # Three-body energy
                    dedc9 = -angles * damp9 * r9

                    # Product rule
                    dafdc9 = -daangles * damp9 * r9 \
                            - angles * dadamp9 * r9 \
                            - angles * damp9 * dar9
                    dbfdc9 = -dbangles * damp9 * r9 \
                            - angles * dbdamp9 * r9 \
                            - angles * damp9 * dbr9
                    dcfdc9 = -dcangles * damp9 * r9 \
                            - angles * dcdamp9 * r9 \
                            - angles * damp9 * dcr9

                    # C9 energies and force contributions
                    eabc += self.c9[a, b, c] * dedc9 * sself
                    fabc[a] -= self.c9[a, b, c] * dafdc9
                    fabc[b] -= self.c9[a, b, c] * dbfdc9
                    fabc[c] -= self.c9[a, b, c] * dcfdc9
                    if a == b:
                        fabc[a] += self.c9[a, b, c] * dbfdc9
                        fabc[b] += self.c9[a, b, c] * dafdc9
                    if a == c:
                        fabc[a] += self.c9[a, b, c] * dcfdc9
                        fabc[c] += self.c9[a, b, c] * dafdc9
                    if b == c:
                        fabc[b] += self.c9[a, b, c] * dcfdc9
                        fabc[c] += self.c9[a, b, c] * dbfdc9
                    fabc -= self.dc9[:, a, b, c] * dedc9 * sself
        self.energy += e6 + e8 + eabc
        self.forces += f6 + f8 + fabc

    def update(self, atoms=None):
        """Update system parameters and calculate energy, forces,
        and stress"""

        if atoms is None:
            atoms = self.calculator.get_atoms()
#        if (self.atoms and
#            (self.atoms.get_positions() == atoms.get_positions()).all()):
#            return
        if self.calculator is not None:
            self.energy = self.calculator.get_potential_energy(atoms)
            self.forces = self.calculator.get_forces(atoms)
            self.stress = self.calculator.get_stress(atoms)
        else:
            self.energy = 0.
            self.forces = np.zeros((len(atoms), 3))
#            self.stress = np.zeros(6)
        self.atoms = atoms.copy()
        self.updateparams(atoms)
        if self.fortran:
            e, f = d3ef.efcalc(
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
        else:
            self.efcalc(atoms)
        cell = atoms.get_cell()
        # Not sure this is correct, leaving it unimplemented for the time being
        #stress = np.sum(self.forces[:, :, np.newaxis] \
        #        * atoms.get_positions()[:, np.newaxis, :], axis=0) \
        #        / atoms.get_volume()
        #self.stress = stress.flat[[0, 4, 8, 5, 2, 1]]

    def get_potential_energy(self, atoms):
        self.update(atoms)
        return self.energy

    def get_forces(self, atoms):
        self.update(atoms)
        return self.forces

    def get_stress(self, atoms):
        raise NotImplementedError
