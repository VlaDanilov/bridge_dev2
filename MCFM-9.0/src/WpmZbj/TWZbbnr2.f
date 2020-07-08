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
 
      function TWZbbnr2(p1,p2,p3,p4,p5,p6,p7,p8)
      implicit none
      include 'types.f'
      complex(dp):: TWZbbnr2
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_com.f'
      include 'sprods_com.f'
C     Author: R.K. Ellis Feb, 2013
C     written by program WZbbdiags.frm
C     These are the non-resonant diagrams
C     Calculation is performed for LH light-line (perforce because of W)
C     Calculation is performed for LH bbbar-line
C     Calculation is performed for LH Z-dcay line
      integer:: p1,p2,p3,p4,p5,p6,p7,p8
      real(dp):: s3,s4,s56,s78,s3456,s356,s178,s278
C     statement functions
      s3(p1,p2,p3)=s(p1,p2)+s(p1,p3)+s(p2,p3)
      s4(p1,p2,p3,p4)=s(p1,p2)+s(p1,p3)+s(p1,p4)
     &               +s(p2,p3)+s(p2,p4)+s(p3,p4)
C     end statement functions
      s56=s(p5,p6)
      s78=s(p7,p8)
      s178=s3(p1,p7,p8)
      s278=s3(p2,p7,p8)
      s356=s3(p3,p5,p6)
      s3456=s4(p3,p4,p5,p6)
      TWZbbnr2= + s3456**(-1)*s56**(-1)*s78**(-1)*s356**(-1)*s178**(-1)
     &  * ( za(p1,p7)*za(p2,p3)*za(p3,p5)*zb(p1,p4)*zb(p1,p8)*zb(p3,p6)
     &     + za(p1,p7)*za(p2,p5)*za(p3,p5)*zb(p1,p4)*zb(p1,p8)*zb(p5,p6
     &    ) + za(p2,p3)*za(p3,p5)*za(p7,p8)*zb(p1,p8)*zb(p3,p6)*zb(p4,
     &    p8) + za(p2,p5)*za(p3,p5)*za(p7,p8)*zb(p1,p8)*zb(p4,p8)*zb(p5
     &    ,p6) )
      TWZbbnr2 = TWZbbnr2 + s3456**(-1)*s56**(-1)*s78**(-1)*s356**(-1)*
     & s278**(-1) * (  - za(p2,p3)*za(p2,p7)*za(p3,p5)*zb(p1,p4)*zb(p2,
     &    p8)*zb(p3,p6) - za(p2,p5)*za(p2,p7)*za(p3,p5)*zb(p1,p4)*zb(p2
     &    ,p8)*zb(p5,p6) + za(p2,p7)*za(p3,p5)*za(p3,p7)*zb(p1,p4)*zb(
     &    p3,p6)*zb(p7,p8) + za(p2,p7)*za(p3,p5)*za(p5,p7)*zb(p1,p4)*
     &    zb(p5,p6)*zb(p7,p8) )

      return
      end
