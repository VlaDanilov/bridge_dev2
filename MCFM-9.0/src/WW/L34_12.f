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
 
      function L34_12(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: L34_12
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6
      integer:: i1,i2,i3,i4
      complex(dp):: z2,bit1,bit2

      real(dp):: t,t134,t234,del3,del12,del56
  
      z2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)

      t134=t(j1,j3,j4)
      t234=t(j2,j3,j4)
      del12=s(j1,j2)-s(j3,j4)-s(j5,j6)
      del56=s(j5,j6)-s(j1,j2)-s(j3,j4)
      Del3=s(j1,j2)**2+s(j3,j4)**2+s(j5,j6)**2
     & -2._dp*(s(j1,j2)*s(j3,j4)+s(j1,j2)*s(j5,j6)+s(j3,j4)*s(j5,j6))

      bit1=
     & +zb(j4,j1)*za(j1,j3)*zb(j3,j5)+zb(j4,j2)*za(j2,j3)*zb(j3,j5)
     & +zb(j4,j1)*za(j1,j4)*zb(j4,j5)+zb(j4,j2)*za(j2,j4)*zb(j4,j5)
      bit2=
     & +zb(j4,j5)*za(j5,j1)*zb(j1,j5)+zb(j4,j6)*za(j6,j1)*zb(j1,j5)
     & +zb(j4,j5)*za(j5,j2)*zb(j2,j5)+zb(j4,j6)*za(j6,j2)*zb(j2,j5)

c-----first line
      L34_12=1.5_dp*del56*(t134-t234)*z2(j3,j1,j2,j4)*z2(j6,j1,j2,j5)
     & /(z2(j2,j5,j6,j1)*Del3**2)
      L34_12=L34_12+1.5_dp*za(j3,j6)*bit1/(z2(j2,j5,j6,j1)*Del3)

c----second line
      L34_12=L34_12+0.5_dp*za(j3,j4)*zb(j4,j5)*bit2
     & /(zb(j5,j6)*z2(j2,j5,j6,j1)*del3)
      L34_12=L34_12+zb(j1,j4)*za(j2,j6)*t134
     & *(za(j3,j6)*del12-2._dp*za(j3,j4)*zb(j4,j5)*za(j5,j6))
     & /(za(j5,j6)*z2(j2,j5,j6,j1)**2*del3)

c-----third line
      L34_12=L34_12+0.5_dp*t134/(z2(j2,j5,j6,j1)*del3)
     & *( +za(j3,j4)*zb(j4,j5)**2/zb(j5,j6)
     &    +zb(j3,j4)*za(j3,j6)**2/za(j5,j6)-2._dp*za(j3,j6)*zb(j4,j5)) 

c--line before penultimate line
      L34_12=L34_12
     & +(+z2(j3,j1,j4,j5)/zb(j5,j6)
     &   -za(j3,j4)*zb(j1,j4)*za(j2,j6)/z2(j2,j5,j6,j1))
     &  *(zb(j4,j5)*del12
     &  -2._dp*zb(j4,j3)*za(j3,j6)*zb(j6,j5))/(z2(j2,j5,j6,j1)*Del3)

c--penultimate line
      L34_12=L34_12+4._dp
     & *(za(j3,j4)*zb(j4,j5)*z2(j6,j1,j3,j4)
     &  +za(j6,j3)*zb(j3,j4)*z2(j3,j2,j4,j5))/(z2(j2,j5,j6,j1)*Del3)

c--last line
      L34_12=L34_12+2._dp*del12/(z2(j2,j5,j6,j1)*Del3)
     & *(+zb(j4,j5)*z2(j3,j2,j4,j5)/zb(j5,j6)
     &   -za(j3,j6)*z2(j6,j1,j3,j4)/za(j5,j6))      

      return
      end


