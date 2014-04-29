"""Grimme D3(BJ) correction scheme for DFT"""

import itertools
import numpy as np
from ase.units import Hartree, Bohr
from ase.calculators.d3params import k1, k2, k3, alp, damp, dampbj, \
        numcn, cn, rcov, r2r4, r0ab, c6ab
from ase.calculators.general import Calculator


class D3:
    """D3(BJ) correction of Grimme et al"""
    def __init__(self,
            bj=True,
            xc='pbe',
            rcut=95. * Bohr,
            rcutcn=40. * Bohr,
            calculator=None,
            ):

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
        self.images = None

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
        pos = atoms.get_positions()

        natoms = len(atoms)

        # Create a list of all images in which there are atoms that we will
        # calculate the interaction with atoms in the central unit cell
        images = []
        for i in xrange(-self.nimages[0], self.nimages[0] + 1):
            for j in xrange(-self.nimages[1], self.nimages[1] + 1):
                for k in xrange(-self.nimages[2], self.nimages[2] + 1):
                    images.append(np.array([i, j, k], dtype=np.int16))
        self.images = np.array(images, dtype=np.int16)

        # Create natoms x natoms x nimages array of displacement vectors
        xyzab = np.array([[[bpos - apos + np.dot(i, cell) for i in images] \
                for b, bpos in enumerate(pos)] \
                for a, apos in enumerate(pos)])

        # Calculate distance squared between all pairs of atoms
        rab2 = np.array([[[np.dot(xyzi,xyzi) for xyzi in xyz] \
                for xyz in xyza] for xyza in xyzab])

        # Mask to exclude atoms interacting with themselves, and atoms
        # interacting with atoms beyond the cutoff radius rcutcn
        cnmask = np.logical_or(rab2 < 1e-6, rab2 > self.rcutcn**2)
        rab2 = np.ma.array(rab2, mask=cnmask)

        # Absolute displacement matrix
        rab = np.sqrt(rab2)

        # Unit vector pointing in the direction of the displacement between
        # atoms.  We will get some divide-by-zeros here, so let's ignore those
        # particular warnings...
        with np.errstate(invalid='ignore', divide='ignore'):
            uxyzab = xyzab / rab[:,:,:,np.newaxis]

        # natoms x natoms x nimages array of rcova + rcovb
        rcovab = np.array([[[self.rcov[a] + self.rcov[b] for i in self.images] \
                for b in xrange(natoms)] for a in xrange(natoms)])

        # Calculate the coordination number and its derivative.  Again,
        # we will get some divide-by-zeros, so we want to ignore it.
        with np.errstate(divide='ignore'):
            cnexp = np.exp(-self.k1 * (rcovab / rab - 1.))

        cnab = 1. / (1. + cnexp)

        with np.errstate(invalid='ignore'):
            dcnab = cnexp[:,:,:,np.newaxis] * self.k1 \
                    * rcovab[:,:,:,np.newaxis] \
                    * cnab[:,:,:,np.newaxis]**2 \
                    * uxyzab / rab2[:,:,:,np.newaxis]

        # Sum over images to get the total coordination number of each atom
        self.cn = cnab.sum((1, 2))
        self.dcn = dcnab.sum(2)
        self.dcn[np.diag_indices(natoms)] = self.dcn.sum(1)

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
            self.r0 = np.zeros((len(self.atoms), len(self.atoms)))
            for i, atomi in enumerate(self.atoms):
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
        if atoms == None:
            atoms = self.atoms
        natoms = len(atoms)
        nimages = len(self.images)
        cell = atoms.get_cell()
#        nallatoms = len(self.allatomindex)

        e6 = 0.
        f6 = np.zeros((natoms, 3))

        e8 = 0.
        f8 = np.zeros((natoms, 3))

        eabc = 0.
        fabc = np.zeros((natoms, 3))

        # Figure out which image is the central image.
        for i, image in enumerate(self.images):
            if np.all(image == np.zeros(3, dtype=np.int16)):
                selfimage = i

        # Set self-interaction correcting term and mask to avoid double
        # counting interactions.
        xyzmask = np.zeros((natoms, natoms, nimages, 3), dtype=np.bool)
        countmask = np.zeros((natoms, natoms, nimages), dtype=np.bool)
        xyzmaskbc = np.zeros((natoms, natoms, nimages, nimages, 3), \
                dtype=np.bool)
        countmaskbc = np.zeros((natoms, natoms, nimages, nimages), \
                dtype=np.bool)

        imagediag = np.arange(nimages)

        # Mask out entries where rab == 0
        for i in xrange(natoms):
            for j in xrange(natoms):
                if i == j:
                    countmask[i, j, selfimage] = True
                    xyzmask[i, j, selfimage, :] = True
                    countmaskbc[i, j, imagediag, imagediag] = True
                    xyzmaskbc[i, j, imagediag, imagediag, :] = True

        # Array of lattice vectors, one for each image we're interacting with
        tvec = np.array([np.dot(image, cell) for image in self.images])

        # Displacement vector matrix (n x n x nimages x 3)
        pos = atoms.get_positions()
        xyzab = np.array(pos[np.newaxis, :, np.newaxis, :] \
                + tvec[np.newaxis, np.newaxis, :, :] \
                - pos[:, np.newaxis, np.newaxis, :])
        xyzab = np.ma.array(xyzab, mask=xyzmask)

        # Distance and distance^2 matrices (n x n x nimages)
        rab2 = np.ma.array([[[np.dot(xyzi, xyzi) \
                for xyzi in xyz] \
                for xyz in xyza] \
                for xyza in xyzab], \
                mask=countmask)
        rab = np.ma.array(np.sqrt(rab2), mask=countmask)

        # Displacement vector matrix (n x n x nimages x nimages x 3)
        xyzbc = np.array(pos[np.newaxis, :, np.newaxis, np.newaxis, :] \
                + tvec[np.newaxis, np.newaxis, :, np.newaxis, :] \
                - pos[:, np.newaxis, np.newaxis, np.newaxis, :] \
                - tvec[np.newaxis, np.newaxis, np.newaxis, :, :])
        xyzbc = np.ma.array(xyzbc, mask=xyzmaskbc)

        # Distance and distance^2 matrices (n x n x nimages x nimages)
        rbc2 = np.ma.array([[[[np.dot(xyzij, xyzij) \
                for xyzij in xyzi] \
                for xyzi in xyz] \
                for xyz in xyzb] \
                for xyzb in xyzbc], \
                mask=countmaskbc)
        rbc = np.ma.array(np.sqrt(rbc2), mask=countmaskbc)


        # Update mask to exclude interactions with distances greater
        # than the cutoff distance
        countmask = np.logical_or(countmask, rab > self.rcut)

        # 3-body interaction has a lower cutoff.
        countmask3b = np.logical_or(countmask, rab > self.rcutcn)
        countmaskbc = np.logical_or(countmaskbc, rbc > self.rcutcn)

        # Choose which type of damping to use, and calculate the damping factor
        # Becke-Johnson damping
        if self.bj:
            with np.errstate(divide='ignore'):
                dedc6 = -1. / (rab**6 + self.dmp6[:, :, np.newaxis])
                dfdc6 = -6. * dedc6[:, :, :, np.newaxis]**2 \
                        * rab[:, :, :, np.newaxis]**4 * xyzab

                dedc8 = -1. / (rab**8 + self.dmp8[:, :, np.newaxis])
                dfdc8 = -8. * dedc8[:, :, :, np.newaxis]**2 \
                        * rab[:, :, :, np.newaxis]**6 * xyzab
        # 'Zero-damping'
        else:
            rav = (self.rs6 * self.r0[:, :, np.newaxis] / rab)**self.alp6
            drav = -xyzab * self.alp6 \
                    * rav[:, :, :, np.newaxis] \
                    / rab2[:, :, :, np.newaxis]
            damp6 = 1. / (1. + 6. * rav)
            ddamp6 = -6. * damp6[:, :, :, np.newaxis]**2 * drav
            dedc6 = -damp6 / rab**6
            dfdc6 = 6. * xyzab \
                    * (dedc6 / rab2)[:, :, :, np.newaxis] \
                    + ddamp6 / rab[:, :, :, np.newaxis]**6

            rav = (self.rs8 * self.r0[:, :, :, np.newaxis] / rab)**self.alp8
            drav = -xyzab * self.alp8 \
                    * rav[:, :, :, np.newaxis] \
                    / rab2[:, :, :, np.newaxis]
            damp8 = 1. / (1. + 6. * rav)
            ddamp8 = -6. * damp8[:, :, :, np.newaxis]**2 * drav
            dedc8 = -damp8 / rab**8
            dfdc8 = 8. * xyzab \
                    * (dedc8 / rab2)[:, :, :, np.newaxis] \
                    + ddamp8 / rab[:, :, :, np.newaxis]**8

        # C6 energy and force contributions
        e6 = self.s6 * np.sum(self.c6[:, :, np.newaxis] * dedc6)

        f6ab = self.c6[:, :, np.newaxis, np.newaxis] * dfdc6
        f6 = self.s6 * (-np.sum(f6ab, axis=(1, 2)) \
               - 0.5 * np.sum(self.dc6[:, :, :, np.newaxis, :] \
                * dedc6[np.newaxis, :, :, :, np.newaxis], \
                axis=(1,2,3)))

        # C8 energy and force contributions
        e8 = self.s18 * np.sum(self.c8[:, :, np.newaxis] * dedc8)

        f8ab = self.c8[:, :, np.newaxis, np.newaxis] * dfdc8
        f8 = self.s18 * (-np.sum(f8ab, axis=(1, 2)) \
                - 0.5 * np.sum(self.dc8[:, :, :, np.newaxis, :] \
                * dedc8[np.newaxis, :, :, :, np.newaxis], \
                axis=(1,2,3)))
        
        # Move on to calculating 3-body interactions, update masks with
        # the new cutoff radius
        rab.mask = countmask3b
        rab2.mask = countmask3b
        
        # Align xyzab, rab, and rab2 for xyzac, rac, and rac2
        # Essentially, b -> c and insert new axes for atom b
        # and its image.
        xyzac = xyzab[:, np.newaxis, :, np.newaxis, :, :]
        rac = rab[:, np.newaxis, :, np.newaxis, :]
        rac2 = rab2[:, np.newaxis, :, np.newaxis, :]
        
        # Insert new axes for atom c and its image
        xyzab = xyzab[:, :, np.newaxis, :, np.newaxis, :]
        rab = rab[:, :, np.newaxis, :, np.newaxis]
        rab2 = rab2[:, :, np.newaxis, :, np.newaxis]
        
        # Insert new axis for atom a
        xyzbc = xyzbc[np.newaxis, :, :, :, :, :]
        rbc = rbc[np.newaxis, :, :, :, :]
        rbc2 = rbc2[np.newaxis, :, :, :, :]

        # Use law of cosines for the sake of simplicity.
        cosalpha = (rab2 + rac2 - rbc2) / (2. * rab * rac)
        cosbeta = (rab2 + rbc2 - rac2) / (2. * rab * rbc)
        cosgamma = (rac2 + rbc2 - rab2) / (2. * rac * rbc)

        # na will add a new axis (x, y, z) to an array
        na = (Ellipsis, np.newaxis)

        # Gradient of cosalpha, cosbeta, and cosgamma.
        dcosalpha = cosalpha[na] * (xyzac / rac2[na] + xyzab / rab2[na]) \
                - (xyzac + xyzab) / (rab * rac)[na]
        dcosbeta = xyzbc / (rab * rbc)[na] + xyzab * (cosbeta / rab2)[na]
        dcosgamma = -xyzbc / (rac * rbc)[na] + xyzac * (cosgamma / rac2)[na]

        # Angle term and its gradient
        angles = 3. * cosalpha * cosbeta * cosgamma + 1.
        dangles = 3. * (dcosalpha * cosbeta[na] * cosgamma[na] \
                + cosalpha[na] * dcosbeta * cosgamma[na] \
                + cosalpha[na] * cosbeta[na] * dcosgamma)

        # rav for the damping of the three-body term and its gradient
        with np.errstate(divide='ignore'):
            rav = 4. / 3.* (self.r0[:, :, np.newaxis, np.newaxis, np.newaxis] \
                    * self.r0[np.newaxis, :, :, np.newaxis, np.newaxis] \
                    * self.r0[:, np.newaxis, :, np.newaxis, np.newaxis] \
                    / (rab * rbc * rac))**(1. / 3.)
        drav = (rav / 3.)[na] * (xyzab / rab2[na] + xyzac / rac2[na])

        # Damping for the three-body term and its gradient
        damp9 = 1. / (1. + 6. * rav**self.alp8)
        ddamp9 = -6. * self.alp8 * (rav**(self.alp8 - 1) * damp9**2)[na] * drav

        # 'Average' r^9 and its gradient
        with np.errstate(divide='ignore'):
            r9 = 1. / (rab * rab2 * rac * rac2 * rbc * rbc2)
        dr9 = 3. * r9[na] * (xyzab / rab2[na] + xyzac / rac2[na])

        dedc9 = -angles * damp9 * r9

        dfdc9 = -dangles * damp9[na] * r9[na] \
                - angles[na] * ddamp9 * r9[na] \
                - angles[na] * damp9[na] * dr9
        
        # Three-body contribution to the energy and forces
        eabc = np.sum(self.c9[:, :, :, np.newaxis, np.newaxis] \
                * dedc9)
        fabc = -np.sum(self.c9[:, :, :, np.newaxis, np.newaxis, np.newaxis] \
                * dfdc9, axis=(1, 2, 3, 4)) \
                - np.sum(self.dc9[:, :, :, :, np.newaxis, np.newaxis, :] \
                * dedc9[np.newaxis, ..., np.newaxis], axis=(1, 2, 3, 4, 5)) / 3.
        
        # We overcounted some interactions, so these factors
        # are to correct for that.
        etot = (e6 + e8) / 2. + eabc / 6.
        ftot = f6 + f8 + fabc / 2.

        # We used masked arrays, and this is just to make sure there
        # aren't any components of the force that are erroneously masked
        assert (ftot.mask == False).all()

        # ftot.data is to strip the mask (which we determined above is
        # all False) from the data.
        self.energy = etot
        self.forces = ftot.data

    def update(self, atoms=None):
        if atoms is None:
            atoms = self.calculator.get_atoms()
#        if (self.atoms and
#            (self.atoms.get_positions() == atoms.get_positions()).all()):
#            return
        if self.calculator:
            self.energy = self.calculator.get_potential_energy(atoms)
            self.forces = self.calculator.get_forces(atoms)
            self.stress = self.calculator.get_stress(atoms)
        else:
            self.energy = 0.
            self.forces = np.zeros((len(atoms), 3))
            self.stress = np.zeros(6)
        self.atoms = atoms.copy()
        if self.cn == None:
            self.updateparams(atoms)

        self.efcalc(atoms)

        cell = atoms.get_cell()
        stress = np.sum(self.forces[:, :, np.newaxis] \
                * atoms.get_positions()[:, np.newaxis, :], axis=0) \
                / (len(atoms) * np.dot(cell[0], np.cross(cell[1], cell[2])))
        self.stress += np.array([stress[0][0], stress[1][1], stress[2][2], \
                stress[1][2], stress[0][2], stress[0][1]])

    def get_potential_energy(self, atoms=None):
        self.update(atoms)
        return self.energy

    def get_forces(self, atoms):
        self.update(atoms)
        return self.forces

    def get_stress(self, atoms):
        self.update(atoms)
        return self.stress
