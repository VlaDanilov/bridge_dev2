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
 
      subroutine dksrwzc(p1,p2,p3,p4,p5,p6,p7,cc11,cc12)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer:: p1,p2,p3,p4,p5,p6,p7
      real(dp):: s12,s34,s126,s345
      complex(dp):: iza,izb,cc11,cc12
      iza(p1,p2)=cone/za(p1,p2)
      izb(p1,p2)=cone/zb(p1,p2)
      s126=s(p1,p2)+s(p1,p6)+s(p2,p6)
      s345=s(p3,p4)+s(p3,p5)+s(p4,p5)
      s34=s(p3,p4)
      s12=s(p1,p2)
      cc11= + s12**(-1)*s34**(-1) * (  - za(p1,p2)*za(p3,p5)*za(p3,p7)*
     &    zb(p1,p6)**2*zb(p3,p4)*izb(p6,p7)*s126**(-1)*s345**(-1) + za(
     &    p1,p2)*za(p3,p5)*za(p5,p7)*zb(p1,p6)**2*zb(p4,p5)*izb(p6,p7)*
     &    s126**(-1)*s345**(-1) + za(p1,p2)*za(p3,p5)*zb(p1,p4)*zb(p1,
     &    p6)*zb(p5,p6)*izb(p5,p7)*izb(p6,p7)*s126**(-1) - za(p1,p2)*
     &    za(p3,p7)*zb(p1,p4)*zb(p1,p6)*izb(p5,p7)*s126**(-1) + za(p2,
     &    p6)*za(p3,p5)*zb(p1,p6)*zb(p4,p6)*zb(p5,p6)*izb(p5,p7)*izb(p6
     &    ,p7)*s126**(-1) - za(p2,p6)*za(p3,p7)*zb(p1,p6)*zb(p4,p6)*
     &    izb(p5,p7)*s126**(-1) )

      cc12= + s12**(-1)*s34**(-1) * ( za(p1,p2)*za(p3,p5)**2*zb(p1,p6)*
     &    zb(p1,p7)*zb(p3,p4)*iza(p5,p7)*s126**(-1)*s345**(-1) - za(p2,
     &    p3)*za(p3,p5)*za(p5,p6)*zb(p1,p6)*zb(p3,p4)*iza(p5,p7)*iza(p6
     &    ,p7)*s345**(-1) - za(p2,p3)*za(p3,p5)*zb(p1,p7)*zb(p3,p4)*
     &    iza(p6,p7)*s345**(-1) + za(p2,p5)*za(p3,p5)*za(p5,p6)*zb(p1,
     &    p6)*zb(p4,p5)*iza(p5,p7)*iza(p6,p7)*s345**(-1) + za(p2,p5)*
     &    za(p3,p5)*zb(p1,p7)*zb(p4,p5)*iza(p6,p7)*s345**(-1) - za(p2,
     &    p6)*za(p3,p5)**2*zb(p1,p6)*zb(p3,p4)*zb(p6,p7)*iza(p5,p7)*
     &    s126**(-1)*s345**(-1) )

      return
      end
