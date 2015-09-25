module d3

   use d3params, only: sp, dp, qp, Bohr, Hartree, k1, k2, k3, alp, max_elem, &
      max_cn, numcn, cntab, r2r4a, rcova, r0

   implicit none

   private
   public d3_calc

   contains

      subroutine d3_calc(natoms, atomnumber, cell, xyz, rcutin, rcutcn, s6, s18, &
            rs6, rs18, alp6, alp8, pbc, bj, threebody, energy, forces, stress)

         use d3params, only: initialize_c6, cross_prod, outer_prod

         implicit none
         !f2py threadsafe

         integer, intent(in) :: natoms, atomnumber(natoms)

         real(dp), intent(in) :: cell(3, 3), xyz(3, natoms), rcutin, rcutcn
         real(dp), intent(in) :: s6, s18, rs6, rs18, alp6, alp8

         logical, intent(in) :: pbc(3), bj, threebody
         
         integer :: a, b, c, i, j, k
         integer :: na, nb, nc, ncna, ncnb, bnum, cnum
         integer :: elemab, elemac, elembc
         integer :: images_2b(3), images_3b(3), natoms_pbc
         integer :: nimages_2b, nimages_3b
         integer, external :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM

         real(dp) :: rcut
         real(dp) :: rcut2, rcutcn2, rcovab
         real(dp) :: cn(natoms)
         real(dp) :: cnexp, cnab, dcnab(3)
         real(dp) :: lij, lijsum, dlij(3, natoms), slij(3, 3)
         real(dp) :: xyza(3), xyzb(3), xyzc(3), xyzab(3), xyzac(3), xyzbc(3)
         real(dp) :: uxyzab(3) 
         real(dp) :: self, lattlen(3), V
         real(dp) :: c8, c9, e6, e8, e9, f6ab(3), f8ab(3)
         real(dp) :: a1, a2, dmp
         real(dp) :: rab, rab2, rab3, rab6, rab8
         real(dp) :: rac, rac2, rac3
         real(dp) :: rbc, rbc2, rbc3
         real(dp) :: dmp6, dmp8, rav6, rav8
         real(dp) :: dmp9, ddmp9, dadmp9(3), dbdmp9(3)
         real(dp) :: r9, dar9(3), dbr9(3)
         real(dp) :: fa9(3), fb9(3), fc9(3)
         real(dp) :: angles, daangles(3), dbangles(3)
         real(dp) :: rav9, darav9(3), dbrav9(3)
         real(dp) :: c6ab, c6ac, c6bc
         real(dp) :: alph, beta, gamm
         real(dp) :: daalph(3), dabeta(3), dagamm(3)
         real(dp) :: dbalph(3), dbbeta(3), dbgamm(3)
         real(dp) :: sc6ab(3, 3), sc6ac(3, 3), sc6bc(3, 3), sc8(3, 3)
         real(dp) :: sc9(3, 3)
         real(dp) :: dc6ab(3, natoms), dc6ac(3, natoms), dc6bc(3, natoms)
         real(dp) :: dc8(3, natoms), dc9(3, natoms)
         real(dp) :: sc_tmp(3, 3), c6abij
         real(dp) :: c6(natoms, natoms), dc6(3, natoms, natoms, natoms)
         real(dp) :: sc6(3, 3, natoms, natoms)
         real(dp) :: r0ab, r0ac, r0bc

         logical, allocatable :: atom_3b(:), central(:)

         integer, allocatable :: atom_pbc(:)
         real(dp), allocatable :: images(:, :)
         real(dp), allocatable :: xyz_pbc(:, :)
         real(dp), allocatable :: dcn(:, :, :), scn(:, :, :)
         real(dp), allocatable :: c6_ref(:,  :, :)

         real(dp), intent(out) :: energy, forces(3, natoms), stress(3, 3)

         if (rcutcn > rcutin) then
            rcut = rcutcn
         else
            rcut = rcutin
         endif

         if (bj) then
            a1 = rs6
            a2 = rs18 * Bohr
         endif

         rcut2 = rcut**2
         rcutcn2 = rcutcn**2

         do i = 1, 3
         lattlen(i) = sqrt(dot_product(cell(i, :), cell(i, :)))
         enddo

         V = dot_product(cell(1, :), cross_prod(cell(2, :), cell(3, :)))

         ! Number of images within 2body cutoff
         images_2b = ceiling(rcut * product(lattlen) / (V * lattlen))
         do i = 1, 3
         if (.not. pbc(i)) images_2b(i) = 0
         enddo
         nimages_2b = product(2 * images_2b + 1)
         natoms_pbc = natoms * nimages_2b

         ! Number of images within 3body cutoff
         images_3b = ceiling(rcutcn * product(lattlen) / (V * lattlen))
         do i = 1, 3
         if (.not. pbc(i)) images_3b(i) = 0
         enddo
         nimages_3b = product(2 * images_3b + 1)

         allocate(images(3, product(2 * images_2b + 1)))

         b = 1
         images(:, 1) = 0

         do i = -images_3b(1), images_3b(1)
         do j = -images_3b(2), images_3b(2)
         do k = -images_3b(3), images_3b(3)
         if ((i == 0) .and. (j == 0) .and. (k == 0)) cycle
         b = b + 1
         images(:, b) = matmul((/i, j, k/), cell)
         enddo
         enddo
         enddo

         do i = -images_2b(1), images_2b(1)
         do j = -images_2b(2), images_2b(2)
         do k = -images_2b(3), images_2b(3)
         if ((abs(i) <= images_3b(1)) .and. &
            (abs(j) <= images_3b(2)) .and. &
            (abs(k) <= images_3b(3))) cycle
         b = b + 1
         images(:, b) = matmul((/i, j, k/), cell)
         enddo
         enddo
         enddo

         allocate(atom_3b(natoms_pbc))
         allocate(central(natoms_pbc))
         allocate(atom_pbc(natoms_pbc))
         allocate(xyz_pbc(3, natoms_pbc))
         
         b = 0
         do a = 1, natoms
         do i = 1, nimages_2b
         b = b + 1
         atom_pbc(b) = a
         xyz_pbc(:, b) = xyz(:, a) + images(:, i)
         if (i == 0) then
            central(b) = .true.
         else
            central(b) = .false.
         endif

         if (i <= nimages_3b) then
            atom_3b(b) = .true.
         else
            atom_3b(b) = .false.
         endif
         enddo
         enddo

         deallocate(images)

         ! Calculate number of 2body and 3body interactions, and calculate
         ! cn & dcn
         allocate(dcn(3, natoms, natoms))
         allocate(scn(3, 3, natoms))

!$OMP PARALLEL default(private) shared(natoms, images_3b, xyz, atomnumber) &
!$OMP shared(cell, rcutcn2, cn, dcn, scn, natoms_pbc, atom_pbc, xyz_pbc, atom_3b) &
!$OMP shared(nimages_2b)
!$OMP WORKSHARE
         cn = 0
         dcn = 0
         scn = 0
!$OMP END WORKSHARE
!$OMP DO reduction(+: cn, dcn, scn) schedule(dynamic, natoms) collapse(2)
         do a = 1, natoms
         do bnum = 1, natoms_pbc

         if (bnum <= nimages_2b * (a - 1) + 1) cycle
         if (.not. atom_3b(bnum)) cycle

         b = atom_pbc(bnum)

         xyza = xyz(:, a)
         na = atomnumber(a)

         xyzb = xyz_pbc(:, bnum)
         nb = atomnumber(b)
         
         xyzab = xyzb - xyza

         rab2 = dot_product(xyzab, xyzab)
         
         if (rab2 < rcutcn2) then
            rab = sqrt(rab2)
            uxyzab = xyzab / rab
            rcovab = rcova(na) + rcova(nb)
            cnexp = exp(-k1 * (rcovab / rab - 1.0_dp))
            cnab = 1.0_dp / (1.0_dp + cnexp)
            dcnab = cnexp * k1 * rcovab * cnab**2 * uxyzab / rab2
            cn(a) = cn(a) + cnab
            sc_tmp = spread(xyzab, 1, 3) * spread(dcnab, 2, 3)
            scn(:, :, a) = scn(:, :, a) + sc_tmp
            if (b > a) then
               cn(b) = cn(b) + cnab

               dcn(:, b, a) = dcn(:, b, a) - dcnab
               dcn(:, a, a) = dcn(:, a, a) + dcnab
               dcn(:, a, b) = dcn(:, a, b) + dcnab
               dcn(:, b, b) = dcn(:, b, b) - dcnab

               scn(:, :, b) = scn(:, :, b) + sc_tmp
            endif
         endif
         enddo
         enddo
!$OMP END DO
!$OMP END PARALLEL

         allocate(c6_ref(max_cn, max_cn, (max_elem + 1) * max_elem / 2))

         call initialize_c6(c6_ref)

!!$OMP PARALLEL default(private) shared(natoms, atomnumber, numcn, cn, cntab) &
!!$OMP shared(dcn, scn, c6_ref, c6, dc6, sc6)
!!$OMP WORKSHARE
         c6 = 0
         dc6 = 0
         sc6 = 0
!!$OMP END WORKSHARE
!!$OMP BARRIER
!!$OMP DO schedule(dynamic, natoms) collapse(2)
         do a = 1, natoms
         do b = a, natoms
         na = atomnumber(a)
         nb = atomnumber(b)
         ncna = numcn(na)
         ncnb = numcn(nb)

         if (na < nb) then
            elemab = max_elem * (na - 1) - (na * (na - 1))/2 + nb
         else
            elemab = max_elem * (nb - 1) - (nb * (nb - 1))/2 + na
         endif

         if ((ncna == 1) .and. (ncnb == 1)) then
            c6(b, a) = c6_ref(1, 1, elemab)
            dc6(:, :, b, a) = 0
            sc6(:, :, b, a) = 0
         else
            lijsum = 0
            dlij = 0
            slij = 0
            c6ab = 0
            dc6ab = 0
            sc6ab = 0
            do i = 1, ncna
            do j = 1, ncnb
            if (na < nb) then
               c6abij = c6_ref(j, i, elemab)
            else
               c6abij = c6_ref(i, j, elemab)
            endif
            lij = exp(-k3 * ((cn(a) - cntab(i, na))**2 + (cn(b) - cntab(j, nb))**2))
            dlij = dlij - 2.0_dp * k3 * lij &
               * ((cn(a) - cntab(i, na)) * dcn(:, :, a) &
               + (cn(b) - cntab(j, nb)) * dcn(:, :, b))
            slij = slij - 2.0_dp * k3 * lij &
               * ((cn(a) - cntab(i, na)) * scn(:, :, a) &
               + (cn(b) - cntab(j, nb)) * scn(:, :, b))
            lijsum = lijsum + lij
            c6ab = c6ab + c6abij * lij
            dc6ab = dc6ab - c6abij &
               * 2.0_dp * k3 * lij &
               * ((cn(a) - cntab(i, na)) * dcn(:, :, a) &
               + (cn(b) - cntab(j, nb)) * dcn(:, :, b))
            sc6ab = sc6ab - c6abij &
               * 2.0_dp * k3 * lij &
               * ((cn(a) - cntab(i, na)) * scn(:, :, a) &
               + (cn(b) - cntab(j, nb)) * scn(:, :, b))
            enddo
            enddo

            c6ab = c6ab / lijsum
            c6(b, a) = c6ab
            dc6(:, :, b, a) = (dc6ab - c6ab * dlij) / lijsum
            sc6(:, :, b, a) = (sc6ab - c6ab * slij) / lijsum
            if (a /= b) then
               c6(a, b) = c6(b, a)
               dc6(:, :, a, b) = dc6(:, :, b, a)
               sc6(:, :, a, b) = sc6(:, :, a, b)
            endif
         endif
         enddo
         enddo
!!$OMP END DO
!!$OMP END PARALLEL

         deallocate(c6_ref)
         deallocate(dcn)
         deallocate(scn)

!$OMP PARALLEL default(private) shared(natoms, xyz, atomnumber) &
!$OMP shared(rcut2, rcutcn2, c6, dc6, sc6, atom_pbc, xyz_pbc) &
!$OMP shared(a1, a2, rs6, alp6, rs18, alp8, natoms_pbc) &
!$OMP shared(bj, threebody, s6, s18, energy, forces, stress, nimages_2b) &
!$OMP shared(atom_3b)
!$OMP WORKSHARE
         energy = 0
         forces = 0
         stress = 0
!$OMP END WORKSHARE
!$OMP DO reduction(+: energy, forces, stress) schedule(dynamic, natoms) collapse(2)
         do a = 1, natoms
         do bnum = 1, natoms_pbc

         if (bnum <= nimages_2b * (a - 1) + 1) cycle

         xyza = xyz(:, a)
         na = atomnumber(a)

         xyzb = xyz_pbc(:, bnum)
         b = atom_pbc(bnum)
         nb = atomnumber(b)

         xyzab = xyzb - xyza
         rab2 = dot_product(xyzab, xyzab)
         if (rab2 > rcut2) cycle

         rab = sqrt(rab2)
         rab2 = rab**2
         rab6 = rab2**3
         rab8 = rab6 * rab2
         uxyzab = xyzab / rab

         if (a == b) then
            self = 0.5_dp
         else
            self = 1.0_dp
         endif

         if (na < nb) then
            elemab = max_elem * (na - 1) - (na * (na - 1))/2 + nb
         else
            elemab = max_elem * (nb - 1) - (nb * (nb - 1))/2 + na
         endif
         c6ab = c6(b, a)
         dc6ab = dc6(:, :, b, a)
         sc6ab = sc6(:, :, b, a)
         r0ab = r0(elemab)

         e6 = 0
         e8 = 0

         c8 = 3.0_dp * c6ab * r2r4a(na) * r2r4a(nb)
         dc8 = 3.0_dp * dc6ab * r2r4a(na) * r2r4a(nb)
         sc8 = 3.0_dp * sc6ab * r2r4a(na) * r2r4a(nb)

         if (bj) then
            dmp = a1 * sqrt(3.0_dp * r2r4a(na) * r2r4a(nb)) + a2
            e6 = -1.0_dp / (rab6 + dmp**6)
            f6ab = -6.0_dp * e6**2 * rab2**2 * xyzab

            e8 = -1.0_dp / (rab8 + dmp**8)
            f8ab = -8.0_dp * e8**2 * rab6 * xyzab
         else
            rav6 = (rs6 * r0ab / rab)**alp6
            dmp6 = 1.0_dp / (1.0_dp + 6.0_dp * rav6)
            e6 = -dmp6 / rab6
            f6ab = 6.0_dp * xyzab * (e6 &
               + dmp6**2 * alp6 * rav6 / rab6) / rab2

            rav8 = (rs18 * r0ab / rab)**alp8
            dmp8 = 1.0_dp / (1.0_dp + 6.0_dp * rav8)
            e8 = -dmp8 / rab8
            f8ab = xyzab * (8.0_dp * e8 &
               + 6.0_dp * dmp8**2 * alp8 * rav8 / rab8) / rab2
         endif

         e6 = e6 * s6 * self
         e8 = e8 * s18 * self

         f6ab = c6ab * f6ab * s6 * self
         f8ab = c8 * f8ab * s18 * self

         energy = energy + c6ab * e6 + c8 * e8

         if (a /= b) then
            forces(:, a) = forces(:, a) - f6ab - f8ab
            forces(:, b) = forces(:, b) + f6ab + f8ab
         endif
         forces = forces - dc6ab * e6 - dc8 * e8

         stress = stress + sc6ab * e6 + sc8 * e8 &
            + outer_prod(3, xyzab, f6ab + f8ab)
         
         if (threebody .and. (rab2 < rcutcn2)) then
            do cnum = bnum + 1, natoms_pbc
            if (.not. atom_3b(cnum)) cycle

            c = atom_pbc(cnum)

            xyzc = xyz_pbc(:, cnum)

            xyzac = xyzc - xyza
            rac2 = dot_product(xyzac, xyzac)

            if (rac2 > rcutcn2) cycle

            xyzbc = xyzc - xyzb
            rbc2 = dot_product(xyzbc, xyzbc)

            if (rbc2 > rcutcn2) cycle

            rab3 = rab2 * rab
            rac = sqrt(rac2)
            rac3 = rac2 * rac
            rbc = sqrt(rbc2)
            rbc3 = rbc2 * rbc

            nc = atomnumber(c)
   
            if (na < nc) then
               elemac = max_elem * (na - 1) - na * (na - 1)/2 + nc
            else
               elemac = max_elem * (nc - 1) - nc * (nc - 1)/2 + na
            endif
   
            if (nb < nc) then
               elembc = max_elem * (nb - 1) - nb * (nb - 1)/2 + nc
            else
               elembc = max_elem * (nc - 1) - nc * (nc - 1)/2 + nb
            endif

            r0ac = r0(elemac)
            r0bc = r0(elembc)

            c6ac = c6(c, a)
            dc6ac = dc6(:, :, c, a)
            sc6ac = sc6(:, :, c, a)

            c6bc = c6(c, b)
            dc6bc = dc6(:, :, c, b)
            sc6bc = sc6(:, :, c, b)

            if (a == b) then
               if (a == c) then
                  self = 1.0_dp / 3.0_dp
               else
                  self = 1.0_dp / 2.0_dp
               endif
            else
               self = 1.0_dp
            endif

            alph = dot_product(xyzab, xyzac) / (rab * rac)
            beta = -dot_product(xyzab, xyzbc) / (rab * rbc)
            gamm = dot_product(xyzac, xyzbc) / (rac * rbc)
   
            ! Gradient of alph, beta, gamm. Very complicated...
            ! Figured this all out using Mathematica and defining
            ! alph = dot_product(xyzab,xyzac)/(rab * rac), etc.
            daalph = alph * (xyzac / rac2 + xyzab / rab2) &
               - (xyzac + xyzab) / (rab * rac)
            dabeta = xyzbc / (rab * rbc) + xyzab * beta / rab2
            dagamm = -xyzbc / (rac * rbc) + xyzac * gamm / rac2
   
            dbalph = xyzac / (rab * rac) - xyzab * alph / rab2
            dbbeta = beta * (xyzbc / rbc2 - xyzab / rab2) &
               + (xyzab - xyzbc) / (rab * rbc)
            dbgamm = -xyzac / (rac * rbc) + xyzbc * gamm / rbc2
   
            ! I have no idea what 'rav' stands for, but that's what Grimme
            ! called this variable.  Cube root of the product of the
            ! ratios of r0ab/rab, times 4/3 for some reason. I don't know.
            rav9 = (4.0_dp/3.0_dp) * (r0ab * r0ac * r0bc &
               / (rab * rbc * rac))**(1.0_dp/3.0_dp)
            darav9 = (rav9/3.0_dp) * (xyzab / rab2 + xyzac / rac2)
            dbrav9 = (rav9/3.0_dp) * (-xyzab / rab2 + xyzbc / rbc2)
   
            ! Three-body term *always* uses "zero" damping, even if
            ! we are using the BJ version of DFT-D3
            dmp9 = 1.0_dp/(1.0_dp + 6.0_dp * rav9**alp8)
            ddmp9 = -6.0_dp * alp8 * rav9**(alp8-1) * dmp9**2
   
            ! Three-body depends on "average" r^9
            r9 = 1.0_dp / (rab3 * rac3 * rbc3)
            dar9 = 3.0_dp * r9 * (xyzab / rab2 + xyzac / rac2)
            dbr9 = 3.0_dp * r9 * (-xyzab / rab2 + xyzbc / rbc2)
   
            ! Angle term of the three body energy, and its gradient
            angles = 3.0_dp * alph * beta * gamm + 1.0_dp
            daangles = 3.0_dp * (daalph * beta * gamm &
               + alph * dabeta * gamm &
               + alph * beta * dagamm)
            dbangles = 3.0_dp * (dbalph * beta * gamm &
               + alph * dbbeta * gamm &
               + alph * beta * dbgamm)
   
            ! Damping derivatives
            dadmp9 = ddmp9 * darav9
            dbdmp9 = ddmp9 * dbrav9
   
            ! Three-body energy
            e9 = -angles * dmp9 * r9 * self

            c9 = -sqrt(c6ab * c6ac * c6bc / Hartree)

            ! Forces
            fa9 = c9 * (-daangles * dmp9 * r9 &
               - angles * dadmp9 * r9 &
               - angles * dmp9 * dar9)

            fb9 = c9 * (-dbangles * dmp9 * r9 &
               - angles * dbdmp9 * r9 &
               - angles * dmp9 * dbr9)

            fc9 = -(fa9 + fb9)

            dc9 = (dc6ab * c6ac * c6bc &
               + c6ab * dc6ac * c6bc &
               + c6ab * c6ac * dc6bc) &
               / (2 * Hartree * c9)

            sc9 = (sc6ab * c6ac * c6bc &
               + c6ab * sc6ac * c6bc &
               + c6ab * c6ac * sc6bc) &
               / (2 * Hartree * c9)

            energy = energy + c9 * e9

            forces(:, a) = forces(:, a) - fa9 * self
            forces(:, b) = forces(:, b) - fb9 * self
            forces(:, c) = forces(:, c) - fc9 * self

            forces = forces - dc9 * e9

            stress = stress + sc9 * e9 &
               - spread(xyza, 1, 3) * spread(fa9, 2, 3) * self &
               - spread(xyzb, 1, 3) * spread(fb9, 2, 3) * self &
               - spread(xyzc, 1, 3) * spread(fc9, 2, 3) * self
            enddo
         endif
         enddo
         enddo
!$OMP END DO
!$OMP END PARALLEL

         deallocate(atom_pbc)
         deallocate(xyz_pbc)

         stress = -(stress + transpose(stress)) / (2.0_dp * V)
      end subroutine d3_calc
end module d3
