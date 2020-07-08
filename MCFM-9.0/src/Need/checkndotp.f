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
 
      subroutine checkndotp(p,n,ig)
      implicit none
      include 'types.f'
c--- routine to check that the vector "n" contracted with the LO
c--- matrix elements, is properly defined such that n.p(ig)=0
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: p(mxpart,4),n(4),test,tolerance
      integer:: ig
      parameter (tolerance=1.e-2_dp)
      logical failedndp
      common/failedndp/failedndp 
!$omp threadprivate(/failedndp/)
      
      failedndp=.false.
      test=n(4)*p(ig,4)-n(1)*p(ig,1)-n(2)*p(ig,2)-n(3)*p(ig,3)
      test=test/max(abs(n(4)),abs(n(1)),abs(n(2)),abs(n(3)))
     &         /max(abs(p(ig,4)),abs(p(ig,1)),abs(p(ig,2)),abs(p(ig,3)))
      
      if (test > tolerance) then
        write(6,*) 'warning: tolerance for n.p check in gvec routine'
     & //' exceeded: ',tolerance,test
c        write(*,*) 'tolerance: ',tolerance
c        write(6,*) 'test: ',test
c        write(*,*) 'index of gluon: ',ig
c        write(6,*) 'n: ',n(1),n(2),n(3),n(4)
c        write(6,*) 'p: ',p(ig,1),p(ig,2),p(ig,3),p(ig,4)
c        write(6,*)
c        call writeout(p)
        call flush(6)
        failedndp=.true.
!        stop
      endif
      
      return
      end
      
