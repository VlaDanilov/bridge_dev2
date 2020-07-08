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
 
      subroutine itransform(p,tp,x,ip,jp,kp)
      implicit none
      include 'types.f'
************************************************************************
*     Author: R.K. Ellis                                               *
*     July, 1998.                                                      *
*     Given p ((n+1)-phase space) produce tp (n-phase space)           *
*     by creating jp                                                   *
************************************************************************
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'npart.f'
      real(dp):: p(mxpart,4),tp(mxpart,4),k(4),kt(4),ks(4),
     & kDp(3:mxpart),ksDp(3:mxpart),kDk,ksDks,x
      integer:: ip,kp,j,nu,jp

      do nu=1,4
      tp(ip,nu)=x*p(ip,nu)
      tp(kp,nu)=p(kp,nu)
c---just so it is non-zero
      tp(jp,nu)=p(jp,nu)

      k(nu) =-p(ip,nu)-p(kp,nu)-p(jp,nu)
      kt(nu) =-tp(ip,nu)-tp(kp,nu)
      ks(nu)=k(nu)+kt(nu)
      enddo

      kDk=k(4)**2-k(1)**2-k(2)**2-k(3)**2
      ksDks=ks(4)**2-ks(1)**2-ks(2)**2-ks(3)**2


      do j=3,npart+2
      if (j == jp) goto 20
      kDp(j)=k(4)*p(j,4)-k(1)*p(j,1)-k(2)*p(j,2)-k(3)*p(j,3)
      ksDp(j)=ks(4)*p(j,4)-ks(1)*p(j,1)-ks(2)*p(j,2)-ks(3)*p(j,3)

      do nu=1,4
      tp(j,nu)=p(j,nu)-two*ksDp(j)*ks(nu)/ksDks+two*kDp(j)*kt(nu)/kDk
      enddo
 20   continue
      enddo

      return
      end


