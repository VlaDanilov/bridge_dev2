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
 
      subroutine catani(mhsq,k1,k2,k3,k4,k5,k6,coeff2,coeff1)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'momwbbm.f'
      include 'scale.f'
      integer:: k1,k2,k3,k4,k5,k6,nu
      real(dp):: mhsq,p2(4),p3(4),p1Dp2,p3Dp4
      complex(dp):: xl12,xl34,coeff2,coeff1
      
      do nu=1,4
      p2(nu)=bp*mom(k2,nu)+bm*mom(k3,nu)
      p3(nu)=bp*mom(k3,nu)+bm*mom(k2,nu)
      enddo
      p1Dp2=mom(k1,4)*p2(4)-mom(k1,1)*p2(1)
     &     -mom(k1,2)*p2(2)-mom(k1,3)*p2(3)
      p3Dp4=mom(k4,4)*p3(4)-mom(k4,1)*p3(1)
     &     -mom(k4,2)*p3(2)-mom(k4,3)*p3(3)
      xl12=cplx1(log(abs(two*p1Dp2)/musq))
      xl34=cplx1(log(abs(two*p3Dp4)/musq))
      
      coeff2=-cone
      coeff1=xl12+xl34+cplx1(8._dp/3._dp-log(mhsq/musq))
      
      return
      end


c--- Note: this is not Catani's formula, this is my reinterpretation 
c---       of the result that we have found      
      subroutine catanisl(mhsq,k1,k2,k3,k4,k5,k6,coeff2,coeff1)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'momwbbm.f'
      include 'scale.f'
      integer:: k1,k2,k3,k4,k5,k6,nu
      real(dp):: mhsq,p2(4),p3(4),be,p2Dp3,p1Dp4
      complex(dp):: xl14,coeff2,coeff1,Lnrat
      
      do nu=1,4
      p2(nu)=bp*mom(k2,nu)+bm*mom(k3,nu)
      p3(nu)=bp*mom(k3,nu)+bm*mom(k2,nu)
      enddo
      p2Dp3=p2(4)*p3(4)-p2(1)*p3(1)
     &     -p2(2)*p3(2)-p2(3)*p3(3)
      p1Dp4=mom(k1,4)*mom(k4,4)-mom(k1,1)*mom(k4,1)
     &     -mom(k1,2)*mom(k4,2)-mom(k1,3)*mom(k4,3)
      xl14=Lnrat(musq,-two*p1Dp4)

      coeff2=-cone
      be=bp-bm
      coeff1=-cone-xl14-cplx1(half*(one+be**2)/be)
     & *Lnrat(mhsq,-two*(p2Dp3+mhsq)*bp**2)
      
      return
      end
     
 
