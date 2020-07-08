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
 
!===== T. Dennen, May 2014
!===== Bubble coefficients for 
!===== q(i1,-)qb(i2,+)ga(i3,-)ga(i4,+)ga(i5,+)ga(i6,+)
      subroutine aaaa_MHV_b_init(i1,i2,i3,i4,i5,i6,za,zb,
     & aaaa_MHV_b)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: i1,i2,i3,i4,i5,i6
      complex(dp):: aaaa_MHV_b(25), zab2

      zab2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)

      aaaa_MHV_b(:) = czip

      aaaa_MHV_b(8) = (
     &  (-2*za(i1,i2)*za(i1,i4)*za(i3,i4)*zab2(i3,i1,i5,i4))/
     &  (za(i1,i5)*za(i2,i4)*za(i2,i6)*za(i4,i5)*za(i4,i6)*
     &    zab2(i4,i1,i5,i4)) + 
     & (2*za(i1,i2)*za(i1,i5)*za(i3,i5)*zab2(i3,i1,i4,i5))/
     &  (za(i1,i4)*za(i2,i5)*za(i2,i6)*za(i4,i5)*za(i5,i6)*
     &    zab2(i5,i1,i4,i5)) + 
     & (za(i1,i6)*za(i3,i6)**2*zab2(i1,i4,i5,i6)**2)/
     &  (za(i1,i4)*za(i1,i5)*za(i2,i6)*za(i4,i6)*za(i5,i6)*
     &    zab2(i6,i2,i3,i6)**2) + 
     & (za(i1,i5)*za(i3,i4)**2*zb(i5,i4)**2)/
     &  (za(i2,i6)*za(i4,i5)*za(i4,i6)*zab2(i4,i1,i5,i4)**2) - 
     & (za(i1,i4)*za(i3,i5)**2*zb(i5,i4)**2)/
     &  (za(i2,i6)*za(i4,i5)*za(i5,i6)*zab2(i5,i1,i4,i5)**2) - 
     & (2*za(i1,i2)**3*za(i2,i3)*za(i3,i6)*zb(i6,i2))/
     &  (za(i1,i4)*za(i1,i5)*za(i2,i4)*za(i2,i5)*za(i2,i6)**2*
     &    zab2(i2,i3,i6,i2)) - 
     & (2*za(i1,i2)*za(i1,i6)**2*za(i2,i3)*za(i3,i6)*zb(i6,i2))/
     &  (za(i1,i4)*za(i1,i5)*za(i2,i6)**2*za(i4,i6)*za(i5,i6)*
     &    zab2(i6,i2,i3,i6))
     & )


      aaaa_MHV_b(9) = (
     &  (2*za(i1,i2)*za(i1,i4)*za(i3,i4)*zab2(i3,i1,i6,i4))/
     &  (za(i1,i6)*za(i2,i4)*za(i2,i5)*za(i4,i5)*za(i6,i4)*
     &    zab2(i4,i1,i6,i4)) + 
     & (za(i1,i5)*za(i3,i5)**2*zab2(i1,i6,i4,i5)**2)/
     &  (za(i1,i4)*za(i1,i6)*za(i2,i5)*za(i4,i5)*za(i6,i5)*
     &    zab2(i5,i2,i3,i5)**2) - 
     & (2*za(i1,i2)*za(i1,i6)*za(i3,i6)*zab2(i3,i1,i4,i6))/
     &  (za(i1,i4)*za(i2,i5)*za(i2,i6)*za(i6,i4)*za(i6,i5)*
     &    zab2(i6,i1,i4,i6)) - 
     & (za(i1,i6)*za(i3,i4)**2*zb(i4,i6)**2)/
     &  (za(i2,i5)*za(i4,i5)*za(i6,i4)*zab2(i4,i1,i6,i4)**2) + 
     & (za(i1,i4)*za(i3,i6)**2*zb(i4,i6)**2)/
     &  (za(i2,i5)*za(i6,i4)*za(i6,i5)*zab2(i6,i1,i4,i6)**2) - 
     & (2*za(i1,i2)**3*za(i2,i3)*za(i3,i5)*zb(i5,i2))/
     &  (za(i1,i4)*za(i1,i6)*za(i2,i4)*za(i2,i5)**2*za(i2,i6)*
     &    zab2(i2,i3,i5,i2)) - 
     & (2*za(i1,i2)*za(i1,i5)**2*za(i2,i3)*za(i3,i5)*zb(i5,i2))/
     &  (za(i1,i4)*za(i1,i6)*za(i2,i5)**2*za(i4,i5)*za(i6,i5)*
     &    zab2(i5,i2,i3,i5))
     & )


      aaaa_MHV_b(10) = (
     &  (za(i1,i4)*za(i3,i4)**2*zab2(i1,i5,i6,i4)**2)/
     &  (za(i1,i5)*za(i1,i6)*za(i2,i4)*za(i5,i4)*za(i6,i4)*
     &    zab2(i4,i2,i3,i4)**2) - 
     & (2*za(i1,i2)*za(i1,i5)*za(i3,i5)*zab2(i3,i1,i6,i5))/
     &  (za(i1,i6)*za(i2,i4)*za(i2,i5)*za(i5,i4)*za(i5,i6)*
     &    zab2(i5,i1,i6,i5)) + 
     & (2*za(i1,i2)*za(i1,i6)*za(i3,i6)*zab2(i3,i1,i5,i6))/
     &  (za(i1,i5)*za(i2,i4)*za(i2,i6)*za(i5,i6)*za(i6,i4)*
     &    zab2(i6,i1,i5,i6)) - 
     & (2*za(i1,i2)**3*za(i2,i3)*za(i3,i4)*zb(i4,i2))/
     &  (za(i1,i5)*za(i1,i6)*za(i2,i4)**2*za(i2,i5)*za(i2,i6)*
     &    zab2(i2,i3,i4,i2)) - 
     & (2*za(i1,i2)*za(i1,i4)**2*za(i2,i3)*za(i3,i4)*zb(i4,i2))/
     &  (za(i1,i5)*za(i1,i6)*za(i2,i4)**2*za(i5,i4)*za(i6,i4)*
     &    zab2(i4,i2,i3,i4)) + 
     & (za(i1,i6)*za(i3,i5)**2*zb(i6,i5)**2)/
     &  (za(i2,i4)*za(i5,i4)*za(i5,i6)*zab2(i5,i1,i6,i5)**2) - 
     & (za(i1,i5)*za(i3,i6)**2*zb(i6,i5)**2)/
     &  (za(i2,i4)*za(i5,i6)*za(i6,i4)*zab2(i6,i1,i5,i6)**2)
     & )


      aaaa_MHV_b(15) = (
     &  (2*za(i1,i2)**3*za(i2,i3)*za(i3,i6)*zb(i6,i2))/
     & (za(i1,i4)*za(i1,i5)*za(i2,i4)*za(i2,i5)*za(i2,i6)**2*
     &   zab2(i2,i3,i6,i2))
     & )


      aaaa_MHV_b(16) = (
     &  (2*za(i1,i2)**3*za(i2,i3)*za(i3,i5)*zb(i5,i2))/
     & (za(i1,i4)*za(i1,i6)*za(i2,i4)*za(i2,i5)**2*za(i2,i6)*
     &   zab2(i2,i3,i5,i2))
     & )


      aaaa_MHV_b(17) = (
     &  (2*za(i1,i2)**3*za(i2,i3)*za(i3,i4)*zb(i4,i2))/
     & (za(i1,i5)*za(i1,i6)*za(i2,i4)**2*za(i2,i5)*za(i2,i6)*
     &   zab2(i2,i3,i4,i2))
     & )


      aaaa_MHV_b(22) = (
     &  (2*za(i1,i2)*za(i1,i3)*za(i3,i4))/
     &  (za(i1,i4)*za(i2,i5)*za(i2,i6)*za(i4,i5)*za(i4,i6)) - 
     & (2*za(i1,i2)*za(i1,i5)*za(i3,i5)*zab2(i3,i1,i4,i5))/
     &  (za(i1,i4)*za(i2,i5)*za(i2,i6)*za(i4,i5)*za(i5,i6)*
     &    zab2(i5,i1,i4,i5)) + 
     & (2*za(i1,i2)*za(i1,i6)*za(i3,i6)*zab2(i3,i1,i4,i6))/
     &  (za(i1,i4)*za(i2,i5)*za(i2,i6)*za(i4,i6)*za(i5,i6)*
     &    zab2(i6,i1,i4,i6)) + 
     & (za(i1,i4)*za(i3,i5)**2*zb(i5,i4)**2)/
     &  (za(i2,i6)*za(i4,i5)*za(i5,i6)*zab2(i5,i1,i4,i5)**2) - 
     & (za(i1,i4)*za(i3,i6)**2*zb(i6,i4)**2)/
     &  (za(i2,i5)*za(i4,i6)*za(i5,i6)*zab2(i6,i1,i4,i6)**2)
     & )


      aaaa_MHV_b(23) = (
     &  (za(i1,i2)**2*za(i1,i3)**2)/
     &  (za(i1,i4)*za(i1,i5)*za(i1,i6)*za(i2,i4)*za(i2,i5)*za(i2,i6)) - 
     & (za(i1,i4)*za(i3,i4)**2*zab2(i1,i2,i3,i4)**2)/
     &  (za(i1,i5)*za(i1,i6)*za(i2,i4)*za(i4,i5)*za(i4,i6)*
     &    zab2(i4,i2,i3,i4)**2) + 
     & (za(i1,i5)*za(i3,i5)**2*zab2(i1,i2,i3,i5)**2)/
     &  (za(i1,i4)*za(i1,i6)*za(i2,i5)*za(i4,i5)*za(i5,i6)*
     &    zab2(i5,i2,i3,i5)**2) - 
     & (za(i1,i6)*za(i3,i6)**2*zab2(i1,i2,i3,i6)**2)/
     &  (za(i1,i4)*za(i1,i5)*za(i2,i6)*za(i4,i6)*za(i5,i6)*
     &    zab2(i6,i2,i3,i6)**2) + 
     & (2*za(i1,i2)*za(i1,i4)**2*za(i2,i3)*za(i3,i4)*zb(i4,i2))/
     &  (za(i1,i5)*za(i1,i6)*za(i2,i4)**2*za(i4,i5)*za(i4,i6)*
     &    zab2(i4,i2,i3,i4)) - 
     & (2*za(i1,i2)*za(i1,i5)**2*za(i2,i3)*za(i3,i5)*zb(i5,i2))/
     &  (za(i1,i4)*za(i1,i6)*za(i2,i5)**2*za(i4,i5)*za(i5,i6)*
     &    zab2(i5,i2,i3,i5)) + 
     & (2*za(i1,i2)*za(i1,i6)**2*za(i2,i3)*za(i3,i6)*zb(i6,i2))/
     &  (za(i1,i4)*za(i1,i5)*za(i2,i6)**2*za(i4,i6)*za(i5,i6)*
     &    zab2(i6,i2,i3,i6))
     & )


      aaaa_MHV_b(24) = (
     &  (2*za(i1,i2)*za(i1,i3)*za(i3,i5))/
     &  (za(i1,i5)*za(i2,i4)*za(i2,i6)*za(i5,i4)*za(i5,i6)) - 
     & (2*za(i1,i2)*za(i1,i4)*za(i3,i4)*zab2(i3,i1,i5,i4))/
     &  (za(i1,i5)*za(i2,i4)*za(i2,i6)*za(i4,i6)*za(i5,i4)*
     &    zab2(i4,i1,i5,i4)) + 
     & (2*za(i1,i2)*za(i1,i6)*za(i3,i6)*zab2(i3,i1,i5,i6))/
     &  (za(i1,i5)*za(i2,i4)*za(i2,i6)*za(i4,i6)*za(i5,i6)*
     &    zab2(i6,i1,i5,i6)) + 
     & (za(i1,i5)*za(i3,i4)**2*zb(i4,i5)**2)/
     &  (za(i2,i6)*za(i4,i6)*za(i5,i4)*zab2(i4,i1,i5,i4)**2) - 
     & (za(i1,i5)*za(i3,i6)**2*zb(i6,i5)**2)/
     &  (za(i2,i4)*za(i4,i6)*za(i5,i6)*zab2(i6,i1,i5,i6)**2)
     & )


      aaaa_MHV_b(25) = (
     &  (2*za(i1,i2)*za(i1,i3)*za(i3,i6))/
     &  (za(i1,i6)*za(i2,i4)*za(i2,i5)*za(i6,i4)*za(i6,i5)) - 
     & (2*za(i1,i2)*za(i1,i4)*za(i3,i4)*zab2(i3,i1,i6,i4))/
     &  (za(i1,i6)*za(i2,i4)*za(i2,i5)*za(i4,i5)*za(i6,i4)*
     &    zab2(i4,i1,i6,i4)) + 
     & (2*za(i1,i2)*za(i1,i5)*za(i3,i5)*zab2(i3,i1,i6,i5))/
     &  (za(i1,i6)*za(i2,i4)*za(i2,i5)*za(i4,i5)*za(i6,i5)*
     &    zab2(i5,i1,i6,i5)) + 
     & (za(i1,i6)*za(i3,i4)**2*zb(i4,i6)**2)/
     &  (za(i2,i5)*za(i4,i5)*za(i6,i4)*zab2(i4,i1,i6,i4)**2) - 
     & (za(i1,i6)*za(i3,i5)**2*zb(i5,i6)**2)/
     &  (za(i2,i4)*za(i4,i5)*za(i6,i5)*zab2(i5,i1,i6,i5)**2)
     & )


      return
      end
