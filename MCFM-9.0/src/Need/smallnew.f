!  Copyright (C) 2019, respective authors of MCFM.
!
!  This program is free software: you can redistribute it and/or modify it under
!  the terms of the GNU General Public License as published by the Free Software
!  Foundation, either version 3 of the License, or (at your option) any later
!  version.
!
!  This program is distributed in the hope that it will be useful, but WITHOUT ANY
!  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!  PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License along with
!  this program. If not, see <http://www.gnu.org/licenses/>
 
      subroutine smallnew(p,npart,*)
! JC, Sep. 2018: applies a cut on dimensionless variables
!       --- cut on Ei/100 GeV [typical hard scale]
!       --- cut on sij/Ei/Ej, i.e. dot product of directions
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'cutoff.f'
      include 'energy.f'
      include 'first.f'
      include 'plabel.f'
      include 'runstring.f'
      include 'kpart.f'
      include 'taucut.f'
      integer:: npart,j,k
      integer, save :: ipp
!$omp threadprivate(ipp)
      real(dp):: p(mxpart,4),dot,pttwo,ptfour,ptsix,checkmin,Eref

!--- on the first call determine beginning of parton entries in plabel
      if (first) then
        first=.false.
        ipp=3
        do while ((ipp < mxpart ) .and. (plabel(ipp) .ne. 'pp') .and. (plabel(ipp) .ne. ''))
          ipp=ipp+1
        enddo
        if (ipp == mxpart) then
          write(6,*) 'Could not identify partons in smallnew.f'
          stop
        endif
!        write(6,*) 'found ipp=',ipp 
      endif

      checkmin=1.e9_dp    ! a big number
      Eref=1.e2_dp

! loops over particles in positions ipp and above since
! they are (usually) the partons that cause problems
!      do j=ipp,npart+2
! Changed to start at 3 for additional safety
      do j=3,npart+2
!        write(6,*) 'Ej',j,abs(p(j,4)/1.e2_dp)
!        write(6,*) 'j1',j,abs(dot(p,j,1)/p(j,4)/p(1,4))
!        write(6,*) 'j2',j,abs(dot(p,j,2)/p(j,4)/p(2,4))
        checkmin=min(checkmin,abs(p(j,4)/Eref))
        checkmin=min(checkmin,abs(dot(p,j,1)/p(j,4)/p(1,4)))
        checkmin=min(checkmin,abs(dot(p,j,2)/p(j,4)/p(2,4)))
        if (j == npart+2) cycle
        do k=j+1,npart+2
!           write(6,*) 'jk',j,k,abs(dot(p,j,k)/p(j,4)/p(k,4))
          checkmin=min(checkmin,abs(dot(p,j,k)/p(j,4)/p(k,4)))
        enddo
      enddo
!      stop

! for processes such as X+jet, ensure that pt(X) is not abnormally small
! (leading to Gram determinant issues in virtual)
!      if ((ipp == 5) .and. (npart > 2)) then
!        checkmin=min(checkmin,pttwo(3,4,p)/Eref)
!      endif

! for the virtual contribution, always cut below 10^-7 for stability,
      if (kpart == kvirt) then
        if (checkmin < 1.e-7_dp) return 1
      endif

! make cut
      if (checkmin < cutoff) return 1
! for non-SCET calculations do not need such a stringent cutoff
      if (.not. usescet) then
        if (checkmin < 1.e-6_dp) return 1
      endif
      

      return
      end

