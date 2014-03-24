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
         added = .FALSE.

         ! Iterate over unit cells
         do i = -imagelist(1), imagelist(1)
         do j = -imagelist(2), imagelist(2)
         do k = -imagelist(3), imagelist(3)
            ! Calculate translation vector
            tvec = matmul((/i,j,k/),cell)

            added = .FALSE.
            
            ! Iterate over pairs of atoms
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

      end subroutine cncalc

      subroutine efcalc(natoms,nallatoms,image,atomindex,xyz,xyzall,dmp6,dmp8,r0,rcut, &
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

         real*8 :: self

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

               ! Self-interaction scaling parameter
               self = 1.d0

               if (b .lt. a)  cycle

               ! Don't calculate the interaction if atom b <= atom a
               ! if atom b is in the central unit cell
               if (b .eq. a) then
                  if (all(image(bnum,:) .eq. (/0,0,0/))) then
                     cycle
                  else
                     self = 0.5d0
                  endif
               endif

               ! Set some variables for interaction between a and b
               xyzab = xyzall(bnum,:) - xyz(a,:)
               rab2 = dot_product(xyzab,xyzab)
               rab = sqrt(rab2)
               uxyzab = xyzab / rab
               
               ! Don't calculate the interaction if rab > rcut
               if (rab .gt. rcut)  cycle

               ! Choose which type of damping to use, and calculate the damping
               ! factor
               ! Becke-Johnson damping
               if (bj) then
                  dedc6 = -1.d0 / (rab**6 + dmp6(a,b))
                  dfdc6 = -6.d0 * dedc6**2 * rab**4 * xyzab

                  dedc8 = -1.d0 / (rab**8 + dmp8(a,b))
                  dfdc8 = -8.d0 * dedc8**2 * rab**6 * xyzab
               ! "Zero-damping"
               else
                  rav = (rs6 * r0(a,b) / rab)**alp6
                  drav = -xyzab * alp6 * rav / rab2
                  damp6 = 1.d0 / (1.d0 + 6.d0 * rav)
                  ddamp6 = -6.d0 * damp6**2 * drav
                  dedc6 = -damp6 / rab**6
                  dfdc6 = 6.d0 * xyzab * dedc6 / rab2 &
                     + ddamp6 / rab**6

                  rav = (rs8 * r0(a,b)/rab)**alp8
                  drav = -xyzab * alp8 * rav / rab2
                  damp8 = 1.d0 / (1.d0 + 6.d0 * rav)
                  ddamp8 = -6.d0 * damp8**2 * drav
                  dedc8 = -damp8 / rab**8
                  dfdc8 = 8.d0 * xyzab * dedc8 / rab2 &
                     + ddamp8 / rab**8
               endif

               ! C6 energy and force contributions
               e6 = e6 + s6 * c6(a,b) * dedc6 * self
               if (b .ne. a) then
                  f6(a,:) = f6(a,:) - s6 * c6(a,b) * dfdc6
                  f6(b,:) = f6(b,:) + s6 * c6(a,b) * dfdc6
               endif
               f6 = f6 - s6 * dc6(:,a,b,:) * dedc6 * self

               ! C8 energy and force contributions
               e8 = e8 + s18 * c8(a,b) * dedc8 * self
               if (b .ne. a) then
                  f8(a,:) = f8(a,:) - s18 * c8(a,b) * dfdc8
                  f8(b,:) = f8(b,:) + s18 * c8(a,b) * dfdc8
               endif
               f8 = f8 - s18 * dc8(:,a,b,:) * dedc8 * self

               ! Don't calculate 3-body term if rab > rcutcn
               if (rab .gt. rcutcn)  cycle
               
               ! Atom c loops over all image atoms starting with bnum+1
               do cnum = 1, nallatoms
                  c = atomindex(cnum)

                  ! 3-body self interaction scaling term.  Can be either 1, 1/2,
                  ! or 1/6
                  self = 1.d0

                  ! Don't calculate interaction if c < a
                  if (c .lt. a)  cycle

                  ! Don't calculate interaction if c is in the central unit cell
                  ! and a == c
                  if (c .eq. a) then
                     if (all(image(cnum,:) .eq. (/0,0,0/)))  cycle
                  endif

                  ! Don't calculate interaction if c < b
                  if (c .lt. b)  cycle

                  ! Don't calculate interaction if c and b are in the same unit
                  ! cell and c == b
                  if (c .eq. b) then
                     if (all(image(cnum,:) .eq. image(bnum,:)))  cycle
                  endif

                  ! Figure out if we're calculating a self energy. If exactly
                  ! two of a, b, and c are the same index, then the system is
                  ! doubly degenerate and the energy must be divided by 2.  If
                  ! all three indices are the same, then the system is 6-fold
                  ! degenerate, and the energy must be divided by 6.
                  if ((a .eq. c) .or. (a .eq. b) .or.  (b .eq. c)) then
                     if ((a .eq. b) .and. (b .eq. c)) then
                        self = 1.d0 / 6.d0
                     else
                        self = 1.d0 / 2.d0
                     endif
                  endif

                  ! Set some variables for a interacting with c
                  xyzac = xyzall(cnum,:) - xyz(a,:)
                  rac2 = dot_product(xyzac,xyzac)
                  rac = sqrt(rac2)

                  ! Don't calculate the interaction if rac > rcutcn
                  if (rac .gt. rcutcn)  cycle

                  ! Set some variables for b interacting with c
                  xyzbc = xyzac - xyzab
                  rbc2 = dot_product(xyzbc,xyzbc)
                  rbc = sqrt(rbc2)

                  ! Don't calculate the interaction if rbc > rcutcn
                  if (rbc .gt. rcutcn)  cycle

                  ! Pre-calculate some powers of rab that we're going to need
                  ! later
                  rab3 = rab * rab2
                  rac3 = rac * rac2
                  rbc3 = rbc * rbc2
   
                  ! Use dot product instead of law of cosines to calculate
                  ! angle terms to be consistent with derivatives below.
                  ! (It doesn't actually matter)
                  cosalpha = dot_product(xyzab,xyzac)/(rab*rac)
                  cosbeta = -dot_product(xyzab,xyzbc)/(rab*rbc)
                  cosgamma = dot_product(xyzac,xyzbc)/(rac*rbc)

                  ! Gradient of cosalpha, cosbeta, cosgamma. Very complicated...
                  ! Figured this all out using Mathematica and defining
                  ! cosalpha = dot_product(xyzab,xyzac)/(rab * rac), etc.
                  dacosalpha = cosalpha * (xyzac / rac2 + xyzab / rab2) &
                     - (xyzac + xyzab) / (rab * rac)
                  dacosbeta = xyzbc / (rab * rbc) + xyzab * cosbeta / rab2
                  dacosgamma = -xyzbc / (rac * rbc) + xyzac * cosgamma / rac2

                  dbcosalpha = xyzac / (rab * rac) - xyzab * cosalpha / rab2
                  dbcosbeta = cosbeta * (xyzbc / rbc2 - xyzab / rab2) &
                     + (xyzab - xyzbc) / (rab * rbc)
                  dbcosgamma = -xyzac / (rac * rbc) + xyzbc * cosgamma / rbc2

                  dccosalpha = xyzab / (rab * rac) - xyzac * cosalpha / rac2
                  dccosbeta = -xyzab / (rab * rbc) - xyzbc * cosbeta / rbc2
                  dccosgamma = -cosgamma * (xyzac / rac2 + xyzbc / rbc2) &
                     + (xyzac + xyzbc) / (rac * rbc)

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

                  ! Three-body depends on "average" r^9
                  r9 = 1.d0 / (rab3 * rac3 * rbc3)
                  dar9 = 3.d0 * r9 * (xyzab / rab2 + xyzac / rac2)
                  dbr9 = 3.d0 * r9 * (-xyzab / rab2 + xyzbc / rbc2)
                  dcr9 = -3.d0 * r9 * (xyzac / rac2 + xyzbc / rbc2)

                  ! The derivatives change if two or more of the atoms are the
                  ! same in different unit cells.  These are not obvious, but
                  ! mathematica/sympy/etc will confirm that this is correct.
                  if ((a .eq. b) .and. (b .ne. c) .and. (a .ne. c)) then
                     dacosalpha = xyzac * cosalpha / rac2 - xyzab / (rab * rac)
                     dacosbeta = xyzbc * cosbeta / rbc2 + xyzab / (rab * rbc)
                     dacosgamma = -dccosgamma

                     dbcosalpha = dacosalpha
                     dbcosbeta = dacosbeta
                     dbcosgamma = dacosgamma

                     darav = -dcrav
                     dbrav = -dcrav

                     dar9 = -dcr9
                     dbr9 = -dcr9
                  elseif ((a .ne. b) .and. (b .eq. c) .and. (a .ne. c)) then
                     dbcosalpha = -dacosalpha
                     dbcosbeta = -xyzab * cosbeta / rab2 - xyzbc / (rab * rbc)
                     dbcosgamma = -xyzac * cosgamma / rac2 + xyzbc / (rac * rbc)

                     dccosalpha = dbcosalpha
                     dccosbeta = dbcosbeta
                     dccosgamma = dbcosgamma

                     dbrav = -darav
                     dcrav = -darav

                     dbr9 = -dar9
                     dcr9 = -dcr9
                  elseif ((a .ne. b) .and. (b .ne. c) .and. (a .eq. c)) then
                     dacosalpha = xyzab * cosalpha / rab2 - xyzac / (rac * rab)
                     dacosbeta = - dbcosbeta
                     dacosgamma = -xyzbc * cosgamma / rbc2 + xyzac / (rbc * rac)

                     dccosalpha = dacosalpha
                     dccosbeta = dacosbeta
                     dccosgamma = dacosgamma

                     darav = -dbrav
                     dcrav = -dbrav

                     dar9 = -dbr9
                     dcr9 = -dbr9
                  elseif ((a .eq. b) .and. (b .eq. c) .and. (a .eq. c)) then
                     dacosalpha = 0.d0
                     dacosbeta = 0.d0
                     dacosgamma = 0.d0
                     
                     dbcosalpha = 0.d0
                     dbcosbeta = 0.d0
                     dbcosgamma = 0.d0

                     dccosalpha = 0.d0
                     dccosbeta = 0.d0
                     dccosgamma = 0.d0

                     darav = 0.d0
                     dbrav = 0.d0
                     dcrav = 0.d0

                     dar9 = 0.d0
                     dbr9 = 0.d0
                     dcr9 = 0.d0
                  endif

                  ! Angle term of the three body energy, and its gradient
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


                  ! Damping derivatives
                  dadamp9 = ddamp9 * darav
                  dbdamp9 = ddamp9 * dbrav
                  dcdamp9 = ddamp9 * dcrav

                  ! Three-body energy
                  dedc9 = -angles * damp9 * r9
                  
                  ! Product rule
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
                  eabc = eabc + c9(a,b,c) * dedc9 * self
                  fabc(a,:) = fabc(a,:) - c9(a,b,c) * dafdc9
                  fabc(b,:) = fabc(b,:) - c9(a,b,c) * dbfdc9
                  fabc(c,:) = fabc(c,:) - c9(a,b,c) * dcfdc9
                  if (a .eq. b) then
                     fabc(a,:) = fabc(a,:) + c9(a,b,c) * dbfdc9
                     fabc(b,:) = fabc(b,:) + c9(a,b,c) * dafdc9
                  endif
                  if (a .eq. c) then
                     fabc(a,:) = fabc(a,:) + c9(a,b,c) * dcfdc9
                     fabc(c,:) = fabc(c,:) + c9(a,b,c) * dafdc9
                  endif
                  if (b .eq. c) then
                     fabc(b,:) = fabc(b,:) + c9(a,b,c) * dcfdc9
                     fabc(c,:) = fabc(c,:) + c9(a,b,c) * dbfdc9
                  endif
                  fabc = fabc - dc9(:,a,b,c,:) * dedc9 * self
               enddo
            enddo
         enddo
         !$OMP END DO
         !$OMP END PARALLEL
         ! No more double counting!
         etot = e6 + e8 + eabc
         ftot = f6 + f8 + fabc
         return
      end subroutine efcalc

end module d3ef
