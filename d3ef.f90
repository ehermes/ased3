module d3ef

   implicit none

   contains

      subroutine cncalc(natoms,nallatoms,imagelist,k1,cell,xyz,rcut,rcutcn,rcov,image, &
            atomindex,xyzall,cn,dcn)

         implicit none

         integer, intent(in) :: natoms, nallatoms, imagelist(3), k1

         real*8, intent(in) :: cell(3,3), xyz(natoms,3)
         real*8, intent(in) :: rcut, rcutcn, rcov(natoms)

         integer :: i, j, k, nadded, a, b

         real*8 :: tvec(3), rcut2, rcutcn2, rcovab
         real*8 :: rab, rab2, xyzab(3), uxyzab(3)
         real*8 :: cnexp, cnab, dcnab(3)

         logical :: added(natoms)

         integer, intent(out) :: image(nallatoms,3), atomindex(nallatoms)
         real*8, intent(out) :: xyzall(nallatoms,3)
         real*8, intent(out) :: cn(natoms)
         real*8, intent(out) :: dcn(natoms,natoms,3)

         rcut2 = rcut**2
         rcutcn2 = rcutcn**2
         nadded = 1

         image = 0
         atomindex = -1
         xyzall = 0.d0
         cn = 0.d0
         dcn = 0.d0

         ! Iterate over unit cells
         do i = -imagelist(1), imagelist(1)
         do j = -imagelist(2), imagelist(2)
         do k = -imagelist(3), imagelist(3)
            ! Calculate translation vector
            tvec = matmul((/i,j,k/),cell)
            
            ! Iterate over pairs of atoms
            added = .FALSE.
            do a = 1, natoms
               do b = 1, natoms
                  ! Don't calculate anything if atom a == atom b
                  if ((a .eq. b) .and. (i .eq. 0) .and. (j .eq. 0) .and. (k .eq. 0)) cycle

                  xyzab = xyz(b,:) + tvec - xyz(a,:)
                  rab2 = dot_product(xyzab,xyzab)

                  ! If rab < rcut, add it to the list of image atoms that we are
                  ! going to calculate interactions between for later
                  if (rab2 < rcut2) then
                     if (.not. added(b)) then
                        atomindex(nadded) = b - 1
                        image(nadded,:) = (/i,j,k/)
                        xyzall(nadded,:) = xyz(b,:) + tvec
                        added(b) = .TRUE.
                        nadded = nadded + 1
                     endif
                  endif

                  ! If rab < rcutcn, add to the coordination number of atom a
                  ! (not b, we will do in a different part of the loop)
                  if (rab2 < rcutcn2) then
                     rab = sqrt(rab2)
                     uxyzab = xyzab / rab
                     rcovab = rcov(a) + rcov(b)
                     cnexp = exp(-k1 * (rcovab / rab - 1.d0))
                     cnab = 1.d0 / (1.d0 + cnexp)
                     dcnab = cnexp * k1 * rcovab * cnab**2 * uxyzab / rab2
                     cn(a) = cn(a) + cnab
                     dcn(a,b,:) = dcnab
                     dcn(a,a,:) = dcn(a,a,:) + dcnab
                  endif
               enddo
            enddo
         enddo
         enddo
         enddo

!         lij = 0.d0
!
!         do a = 1, natoms
!         do b = 1, natoms
!            do i = 1, numcn(a)
!            do j = 1, numcn(b)
!               lij = 0.d0
!               lij(i,j) = exp(-k3*((cn(a) - cni(a,i))**2 + (cn(b) - cni(b,j))**2))
!               dlij(:,i,j,:) = -2 * lij(i,j) * k3 * ((cn(a) - cni(a,i)) * dcn(:,a,:) &
!                  + (cn(b) - cni(b,j)) * dcn(:,b,:)
!            enddo
!            enddo
!            c6(a,b) = sum(c6abij * lij) / sum(lij)
!            dc6(:,a,b,:) = (-c6(a,b) * 


      end subroutine cncalc

      subroutine e2e3(natoms,nallatoms,image,atomindex,xyz,xyzall,dmp6,dmp8,r0,rcut, &
            rcutcn,c6,dc6,c8,dc8,c9,dc9,s6,s18,rs6,rs8,alp6,alp8,bj,etot,ftot)
         
         implicit none

         integer, intent(in) :: natoms, nallatoms

         integer, intent(in) :: image(nallatoms,3), atomindex(nallatoms)

         real*8, intent(in) :: xyz(natoms,3), xyzall(nallatoms,3)
         real*8, intent(in) :: dmp6(natoms,natoms), dmp8(natoms,natoms)
         real*8, intent(in) :: r0(natoms,natoms), rcut, rcutcn
         real*8, intent(in) :: s6, s18, rs6, rs8, alp6, alp8
         real*8, intent(in) :: c6(natoms,natoms), c8(natoms,natoms)
         real*8, intent(in) :: c9(natoms,natoms,natoms)
         real*8, intent(in) :: dc6(natoms,natoms,natoms,3)
         real*8, intent(in) :: dc8(natoms,natoms,natoms,3)
         real*8, intent(in) :: dc9(natoms,natoms,natoms,natoms,3)

         logical, intent(in) :: bj

         integer :: a, b, c, bnum, cnum

         real*8 :: xyzab(3), xyzac(3), xyzbc(3)
         real*8 :: xyzab2(3), xyzac2(3), xyzbc2(3)
         real*8 :: uxyzab(3), uxyzac(3), uxyzbc(3)

         real*8 :: rab, rab2, rab3, rac, rac2, rac3
         real*8 :: rbc, rbc2, rbc3

         real*8 :: dedc6, dedc8, dedc9, damp6, damp8, damp9
         real*8 :: ddamp6(3), ddamp8(3)
         real*8 :: ddamp9(3), dadamp9(3), dbdamp9(3), dcdamp9(3)

         real*8 :: dfdc6(3), dfdc8(3)
         real*8 :: dfdc9(3), dafdc9(3), dbfdc9(3), dcfdc9(3)

         real*8 :: cosalpha, cosbeta, cosgamma
         real*8 :: dacosalpha(3), dacosbeta(3), dacosgamma(3)
         real*8 :: dbcosalpha(3), dbcosbeta(3), dbcosgamma(3)
         real*8 :: dccosalpha(3), dccosbeta(3), dccosgamma(3)

         real*8 :: angles, rav, r9, drav(3)
         real*8 :: daangles(3), dbangles(3), dcangles(3)
         real*8 :: darav(3), dbrav(3), dcrav(3)
         real*8 :: dar9(3), dbr9(3), dcr9(3)

         real*8 :: e6, e8, eabc 
         real*8 :: f6(natoms,3), f8(natoms,3), fabc(natoms,3)

         real*8, intent(out) :: etot
         real*8, intent(out) :: ftot(natoms,3)

         e6 = 0.d0
         f6 = 0.d0
         e8 = 0.d0
         f8 = 0.d0
         eabc = 0.d0
         fabc = 0.d0

         !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(e6,f6,e8,f8,eabc,fabc) &
         !$OMP SHARED(natoms,nallatoms,image,atomindex,xyz,xyzall,dmp6,dmp8) &
         !$OMP SHARED(r0,rcut,rcutcn,c6,dc6,c8,dc8,c9,dc9,s6,s18,rs6,rs8) &
         !$OMP SHARED(alp6,alp8,bj,etot,ftot)

         ! Atom a loops over atoms in the central unit cell
         !$OMP DO REDUCTION(+:e6,e8,eabc,f6,f8,fabc)
         do a = 1, natoms
            ! Atom b loops over all image atoms
            do bnum = 1, nallatoms
               b = atomindex(bnum)

               ! Don't calculate the interaction if atom b <= atom a
               ! if atom b is in the central unit cell
               if (all(image(bnum,:) .eq. (/0,0,0/))) then
                  if (b .le. a)  cycle
               endif

               xyzab = xyzall(bnum,:) - xyz(a,:)
               rab2 = dot_product(xyzab,xyzab)
               rab = sqrt(rab2)
               uxyzab = xyzab / rab
               
               ! Don't calculate the interaction if rab > rcut
               if (rab .gt. rcut)  cycle

               ! Choose which type of damping to use, and calculate the damping
               ! factor
               if (bj) then
                  dedc6 = -1.d0 / (rab**6 + dmp6(a,b))
                  dfdc6 = -6.d0 * dedc6**2 * rab**5 * uxyzab

                  dedc8 = -1.d0 / (rab**8 + dmp8(a,b))
                  dfdc8 = -8.d0 * dedc8**2 * rab**7 * uxyzab
               else
                  rav = (rs6 * r0(a,b) / rab)**alp6
                  drav = -uxyzab * alp6 * rav / rab
                  damp6 = 1.d0 / (1.d0 + 6.d0 * rav)
                  ddamp6 = -6.d0 * damp6**2 * drav
                  dedc6 = -damp6 / rab**6
                  dfdc6 = 6.d0 * uxyzab * dedc6 / rab &
                     + ddamp6 / rab**6

                  rav = (rs8 * r0(a,b)/rab)**alp8
                  drav = -uxyzab * alp8 * rav / rab
                  damp8 = 1.d0 / (1.d0 + 6.d0 * rav)
                  ddamp8 = -6.d0 * damp8**2 * drav
                  dedc8 = -damp8 / rab**8
                  dfdc8 = 8.d0 * uxyzab * dedc8 / rab &
                     + ddamp8 / rab**8
               endif

               ! C6 energy and force contributions
               e6 = e6 + s6 * c6(a,b) * dedc6
               f6(a,:) = f6(a,:) + s6 * c6(a,b) * dfdc6
               f6(b,:) = f6(b,:) - s6 * c6(a,b) * dfdc6
               f6 = f6 + s6 * dc6(:,a,b,:) * dedc6

               ! C8 energy and force contributions
               e8 = e8 + s18 * c8(a,b) * dedc8
               f8(a,:) = f8(a,:) + s18 * c8(a,b) * dfdc8
               f8(b,:) = f8(b,:) - s18 * c8(a,b) * dfdc8
               f8 = f8 + s18 * dc8(:,a,b,:) * dedc8

               ! Don't calculate 3-body term if rab > rcutcn
               if (rab .gt. rcutcn)  cycle
               
               ! Atom c loops over all image atoms starting with bnum+1
               do cnum = bnum + 1, nallatoms
                  c = atomindex(cnum)

                  ! Don't calculate the interaction if atom c <= atom a
                  ! within the central unit cell
                  if (all(image(cnum,:) .eq. (/0,0,0/))) then
                     if (c .le. a)  cycle
                  endif

                  xyzac = xyzall(cnum,:) - xyz(a,:)
                  rac2 = dot_product(xyzac,xyzac)
                  rac = sqrt(rac2)

                  ! Don't calculate the interaction if rac > rcutcn
                  if (rac .gt. rcutcn)  cycle

                  xyzbc = xyzac - xyzab
                  rbc2 = dot_product(xyzbc,xyzbc)
                  rbc = sqrt(rbc2)

                  ! Don't calculate the interaction if rbc > rcutcn
                  if (rbc .gt. rcutcn)  cycle

                  ! Pre-calculate some powers of rab and its vector
                  ! that we're going to need later
                  rab3 = rab * rab2
                  rac3 = rac * rac2
                  rbc3 = rbc * rbc2
                  xyzab2 = xyzab**2
                  xyzac2 = xyzac**2
                  xyzbc2 = xyzbc**2
   
                  ! Use dot product instead of law of cosines to calculate
                  ! angle terms to be consistent with derivatives below.
                  cosalpha = dot_product(xyzab,xyzac)/(rab*rac)
                  cosbeta = -dot_product(xyzab,xyzbc)/(rab*rbc)
                  cosgamma = dot_product(xyzac,xyzbc)/(rac*rbc)

                  ! Gradient of cosalpha, cosbeta, cosgamma. Very complicated...
                  ! Figured this all out using Mathematica and defining
                  ! cosalpha = dot_product(xyzab,xyzac)/(rab * rac), etc.
                  dacosalpha = cosalpha * (xyzac / rac2 + xyzab / rab2) &
                     - (xyzab + xyzac) / (rab * rac)
                  dacosbeta = (-xyz(a,:) * (rab * rbc * cosbeta + xyzab * xyzbc) &
                     - xyz(b,:) * (rab * rac * cosalpha - xyzab * xyzac) &
                     + xyz(c,:) * (rab2 - xyzab2))/(rab3 * rbc)
                  dacosgamma = (-xyz(a,:) * (rac * rbc * cosgamma - xyzac * xyzbc) &
                     - xyz(c,:) * (rab * rac * cosalpha - xyzab * xyzac) &
                     + xyz(b,:) * (rac2 - xyzac2))/(rac3 * rbc)

                  dbcosalpha = (-xyz(b,:) * (rab * rac * cosalpha + xyzab * xyzac) &
                     - xyz(a,:) * (rab * rbc * cosbeta - xyzab * xyzbc) &
                     + xyz(c,:) * (rab2 + xyzab2))/(rab3 * rac)
                  dbcosbeta = cosbeta * (xyzbc / rbc2 - xyzab / rab2) &
                     - (xyzbc - xyzab) / (rbc * rab)
                  dbcosgamma = (-xyz(c,:) * (rab * rbc * cosbeta + xyzab * xyzbc) &
                     - xyz(b,:) * (rbc * rac * cosgamma - xyzbc * xyzac) &
                     + xyz(a,:) * (rbc2 - xyzbc2))/(rbc3 * rac)

                  dccosalpha = (-xyz(a,:) * (rac * rbc * cosgamma - xyzac * xyzbc) &
                     - xyz(c,:) * (rab * rac * cosalpha - xyzab * xyzac) &
                     + xyz(b,:) * (rac2 - xyzac2))/(rac3 * rab)
                  dccosbeta = (-xyz(c,:) * (rab * rbc * cosbeta + xyzab * xyzbc) &
                     - xyz(b,:) * (rac * rbc * cosgamma - xyzac * xyzbc) &
                     + xyz(a,:) * (rbc2 - xyzbc2))/(rbc3 * rab)
                  dccosgamma = -cosgamma * (xyzbc / rbc2 + xyzac / rac2) &
                     + (xyzac + xyzbc) / (rac * rbc)

                  angles = 3.d0 * cosalpha * cosbeta * cosgamma + 1.d0
                  daangles = 3.d0 * (dacosalpha * cosbeta * cosgamma &
                     + cosalpha * dacosbeta * cosgamma &
                     + cosalpha * cosbeta * dacosgamma)
                  dbangles = 3.d0 * (dbcosalpha * cosbeta * cosgamma &
                     + cosalpha * dbcosbeta * cosgamma &
                     + cosalpha * cosbeta * dbcosgamma)
                  dcangles = 3.d0 * (dccosalpha * cosbeta * cosgamma &
                     + cosalpha * dccosbeta * cosgamma &
                     + cosalpha * cosbeta * dccosgamma)

                  ! I have no idea what 'rav' stands for, but that's what Grimme
                  ! called this variable.  Cube root of the product of the
                  ! ratios of r0ab/rab, times 4/3 for some reason. I don't know.
                  rav = (4.d0/3.d0) * (r0(a,b) * r0(b,c) * r0(a,c) &
                     / (rab * rbc * rac))**(1.d0/3.d0)
                  darav = (rav/3.d0) * (xyzab / rab2 + xyzac / rac2)
                  dbrav = (rav/3.d0) * (-xyzab / rab2 + xyzbc / rbc2)
                  dcrav = -(rav/3.d0) * (xyzac / rac2 + xyzbc / rbc2)

                  ! Three-body term *always* uses "zero" damping, even if
                  ! we are using the BJ version of DFT-D3
                  damp9 = 1.d0/(1.d0 + 6.d0 * rav**alp8)
                  ddamp9 = -6.d0 * alp8 * rav**(alp8-1) * damp9**2
                  dadamp9 = ddamp9 * darav
                  dbdamp9 = ddamp9 * dbrav
                  dcdamp9 = ddamp9 * dcrav

                  ! Three-body depends on "average" r^9
                  r9 = 1.d0 / (rab3 * rac3 * rbc3)
                  dar9 = 3.d0 * r9 * (xyzab / rab2 + xyzac / rac2)
                  dbr9 = 3.d0 * r9 * (-xyzab / rab2 + xyzbc / rbc2)
                  dcr9 = -3.d0 * r9 * (xyzac / rac2 + xyzbc / rbc2)

                  ! Product rule
                  dedc9 = -angles * damp9 * r9
                  
                  dafdc9 = -daangles * damp9 * r9 &
                     - angles * dadamp9 * r9 &
                     - angles * damp9 * dar9
                  dbfdc9 = -dbangles * damp9 * r9 &
                     - angles * dbdamp9 * r9 &
                     - angles * damp9 * dbr9
                  dcfdc9 = -dcangles * damp9 * r9 &
                     - angles * dcdamp9 * r9 &
                     - angles * damp9 * dcr9

                  ! C9 energy and force contributions
                  eabc = eabc + c9(a,b,c) * dedc9
                  fabc(a,:) = fabc(a,:) + c9(a,b,c) * dafdc9
                  fabc(b,:) = fabc(b,:) + c9(a,b,c) * dbfdc9
                  fabc(c,:) = fabc(c,:) + c9(a,b,c) * dcfdc9
                  fabc = fabc + dc9(:,a,b,c,:) * dedc9
               enddo
            enddo
         enddo
         !$OMP END DO
         !$OMP END PARALLEL
         ! No more double counting!
         etot = e6 + e8 + eabc
         ftot = f6 + f8 + fabc
         return
      end subroutine e2e3

end module d3ef
