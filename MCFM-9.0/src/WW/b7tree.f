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
 
      function B7treea(j1,j2,j3,j4,j5,j6,j7,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: B7treea
C---This function in the notation of DKS
c---Multiplied by a factor of (-i)
C---in the natural way this function refers to the process
c---u(2)+dbar(1)-->nu(3)+e^+(3)+b(5)+bbar(6)+g(7)
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6,j7
      integer:: i1,i2,i3,i4
      complex(dp):: z2
      real(dp):: s127,s356
      z2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
      s127=s(j1,j2)+s(j1,j7)+s(j2,j7)
      s356=s(j3,j5)+s(j3,j6)+s(j5,j6)
      B7treea=z2(j1,j2,j7,j4)*z2(j1,j3,j6,j5)*za(j3,j6)
     & /(za(j1,j7)*za(j2,j7)*s(j5,j6)*s127*s356)
      return 
      end


      function B7treeb(j1,j2,j3,j4,j5,j6,j7,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: B7treeb
C---This function in the notation of DKS
c---Multiplied by a factor of (-i)
C---in the natural way this function refers to the process
c---u(2)+dbar(1)-->nu(3)+e^+(3)+b(5)+bbar(6)+g(7)
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6,j7
      integer:: i1,i2,i3,i4
      complex(dp):: z2
      real(dp):: s127,s345
      z2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
      s127=s(j1,j2)+s(j1,j7)+s(j2,j7)
      s345=s(j3,j4)+s(j4,j5)+s(j3,j5)

      B7treeb=-za(j2,j6)*zb(j4,j5)
     & *(za(j1,j2)*z2(j3,j4,j5,j1)-za(j2,j7)*z2(j3,j4,j5,j7))
     & /(za(j1,j7)*za(j2,j7)*s(j3,j4)*s127*s345)
      return 
      end
