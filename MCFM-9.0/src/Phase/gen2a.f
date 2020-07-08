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
 
      subroutine gen2a(r,p,wt2,*)
      implicit none
      include 'types.f'
c----1+2 --> 3+4
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'mxdim.f'
      include 'phasemin.f'
      include 'x1x2.f'
      integer:: nu
      real(dp):: r(mxdim),p(mxpart,4),wt2
      real(dp):: p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),p7(4)
      include 'energy.f'

      p(:,:)=0._dp
      xx(1)=xmin+(1._dp-xmin)*r(1)
      xx(2)=xmin+(1._dp-xmin)*r(2)

      p1(4)=-0.5_dp*xx(1)*sqrts
      p1(1)=0._dp
      p1(2)=0._dp
      p1(3)=-0.5_dp*xx(1)*sqrts
      
      p2(4)=-0.5_dp*xx(2)*sqrts
      p2(1)=0._dp
      p2(2)=0._dp
      p2(3)=+0.5_dp*xx(2)*sqrts

      call phase2(r,p1,p2,p3,p4,wt2,*999)

c---5,6,7=dummies -- just so nothing crashes
      do nu=1,4
      p(1,nu)=p1(nu)
      p(2,nu)=p2(nu)
      p(3,nu)=p3(nu)
      p(4,nu)=p4(nu)
      enddo
      wt2=wt2*(1._dp-xmin)**2
      return

 999  continue
      wt2=0._dp
      p(:,:)=0._dp
      return
      end
