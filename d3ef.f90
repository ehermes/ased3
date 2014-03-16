module d3ef

   implicit none

   contains

      subroutine e2e3(natoms,nallatoms,image,atomindex,xyz,xyzall,dmp6,dmp8,r0,rmax, &
            c6,c8,c9,s6,s18,rs6,rs8,alp6,alp8,bj,etot)
         
         implicit none

         integer, intent(in) :: natoms, nallatoms

         integer, intent(in) :: image(nallatoms,3), atomindex(nallatoms)

         real*8, intent(in) :: xyz(natoms,3), xyzall(nallatoms,3)
         real*8, intent(in) :: dmp6(natoms,natoms), dmp8(natoms,natoms)
         real*8, intent(in) :: r0(natoms,natoms), rmax
         real*8, intent(in) :: s6, s18, rs6, rs8, alp6, alp8
         real*8, intent(in) :: c6(natoms,natoms), c8(natoms,natoms)
         real*8, intent(in) :: c9(natoms,natoms,natoms)

         logical, intent(in) :: bj

         integer :: a, b, c, bnum, cnum, nint
         logical :: ecalc

         real*8 :: xyzab(3), xyzac(3), xyzbc(3)
         real*8 :: rab, rab2, rac, rac2, rbc, rbc2
         real*8 :: dedc6, dedc8, dedc9
         real*8 :: cosalpha, cosbeta, cosgamma, rav
         real*8 :: e6, e8, eabc

         real*8, intent(out) :: etot

         !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(e6, e8, eabc, nint) &
         !$OMP FIRSTPRIVATE(natoms,nallatoms,image,atomindex,xyz,xyzall) &
         !$OMP FIRSTPRIVATE(dmp6,dmp8,r0,rmax,s6,s18,rs6,rs8,alp6,alp8) &
         !$OMP FIRSTPRIVATE(c6,c8,c9,bj)
         e6 = 0.d0
         e8 = 0.d0
         eabc = 0.d0

         nint = 0

         print *, "hi!"

         !$OMP DO REDUCTION(+:e6,e8,eabc,nint)
         do a = 1, natoms
            do bnum = 1, nallatoms
               b = atomindex(bnum)
               ecalc = .TRUE.
               if (all(image(bnum,:) .eq. (/0,0,0/))) then
                  if (b .eq. a) then
                     cycle
                  endif
               endif
               xyzab = xyzall(bnum,:) - xyz(a,:)
               rab2 = dot_product(xyzab,xyzab)
               rab = sqrt(rab2)
               if (rab .gt. rmax) then
                  cycle
               endif
               if (bj) then
                  dedc6 = -1.d0 / (rab**6 + dmp6(a,b))
                  dedc8 = -1.d0 / (rab**8 + dmp8(a,b))
               else
                  dedc6 = -1.d0 / ((rab**6) * (1.d0 + 6.d0 * (rs6 * r0(a,b)/rab)**alp6))
                  dedc8 = -1.d0 / ((rab**8) * (1.d0 + 6.d0 * (rs8 * r0(a,b)/rab)**alp8))
               endif
               if (ecalc) then
                  nint = nint + 1
                  e6 = e6 + s6 * c6(a,b) * dedc6
                  e8 = e8 + s18 * c8(a,b) * dedc8
               endif
               do cnum = 1, nallatoms
                  c = atomindex(cnum)
                  if (all(image(cnum,:) .eq. (/0,0,0/))) then
                     if (c .eq. a) then
                        cycle
                     endif
                  endif
                  if (cnum .eq. bnum)  cycle
                  xyzac = xyzall(cnum,:) - xyz(a,:)
                  rac2 = dot_product(xyzac,xyzac)
                  rac = sqrt(rac2)
                  if (rac .gt. rmax) then
                     cycle
                  endif
                  xyzbc = xyzac - xyzab
                  rbc2 = dot_product(xyzbc,xyzbc)
                  rbc = sqrt(rbc2)
                  if (rbc .gt. rmax) then
                     cycle
                  endif
                  cosalpha = (rab2 + rac2 - rbc2) / (2.d0 * rab * rac)
                  cosbeta  = (rab2 + rbc2 - rac2) / (2.d0 * rab * rbc)
                  cosgamma = (rac2 + rbc2 - rab2) / (2.d0 * rac * rbc)
                  rav = (4.d0/3.d0) * (r0(a,b) * r0(b,c) * r0(a,c) &
                     / (rab * rbc * rac))**(1.d0/3.d0)
                  
                  dedc9 = -(3.d0 * cosalpha * cosbeta * cosgamma + 1.d0) &
                     / (rab * rac * rbc * rab2 * rac2 * rbc2 &
                     * (1.d0 + 6.d0 * rav ** alp8))
                  if (ecalc) then
                     eabc = eabc + c9(a,b,c) * dedc9
                  endif
               enddo
            enddo
         enddo
         !$OMP END DO
         !$OMP END PARALLEL
         etot = (e6 + e8)/2.d0 + eabc/6.d0
         print *, "Number of interactions:", nint
         return
      end subroutine e2e3

end module d3ef
