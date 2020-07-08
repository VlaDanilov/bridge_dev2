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
 
      function TWbbmZab(p1,p2,p3,p4,p5,p6,p7,p8)
      implicit none
      include 'types.f'
      complex(dp):: TWbbmZab
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_com.f'
      include 'sprods_com.f'
C     Author: R.K. Ellis Feb, 2013
C     written by program WbbZdiags.frm
C     These are the Abelian diagrams with both the W and Z coming
C     off light line, only three diagrams with Z emitted before W
C     Calculation is performed for LH light-line (perforce because of W)
C     Calculation is performed for LH bbbar-line
C     Calculation is performed for LH Z-dcay line
      integer:: p1,p2,p3,p4,p5,p6,p7,p8
      real(dp):: s3,s4,s34,s56,s134,s234,s567,s568,s5678
      complex(dp):: d9m,d10m
C     statement functions
      s3(p1,p2,p3)=s(p1,p2)+s(p1,p3)+s(p2,p3)
      s4(p1,p2,p3,p4)=s(p1,p2)+s(p1,p3)+s(p1,p4)
     &               +s(p2,p3)+s(p2,p4)+s(p3,p4)
C     end statement functions
      s34=s(p3,p4)
      s56=s(p5,p6)
      s134=s3(p1,p3,p4)
      s234=s3(p2,p3,p4)
      s567=s3(p5,p6,p7)
      s568=s3(p5,p6,p8)
      s5678=s4(p5,p6,p7,p8)
      d9m= + s34**(-1)*s56**(-1)*s134**(-1)*s567**(-1) * (  - za(p1,p3)
     &    *za(p2,p5)*za(p5,p7)*zb(p1,p4)*zb(p1,p8)*zb(p5,p6)*
     &    s5678**(-1) + za(p1,p3)*za(p2,p7)*za(p5,p7)*zb(p1,p4)*zb(p1,
     &    p8)*zb(p6,p7)*s5678**(-1) + za(p2,p5)*za(p3,p4)*za(p5,p7)*zb(
     &    p1,p4)*zb(p4,p8)*zb(p5,p6)*s5678**(-1) - za(p2,p7)*za(p3,p4)*
     &    za(p5,p7)*zb(p1,p4)*zb(p4,p8)*zb(p6,p7)*s5678**(-1) )
      d9m = d9m + s34**(-1)*s56**(-1)*s134**(-1)*s568**(-1) * (  - za(
     &    p1,p3)*za(p2,p7)*za(p5,p6)*zb(p1,p4)*zb(p1,p6)*zb(p6,p8)*
     &    s5678**(-1) - za(p1,p3)*za(p2,p7)*za(p5,p8)*zb(p1,p4)*zb(p1,
     &    p8)*zb(p6,p8)*s5678**(-1) + za(p2,p7)*za(p3,p4)*za(p5,p6)*zb(
     &    p1,p4)*zb(p4,p6)*zb(p6,p8)*s5678**(-1) + za(p2,p7)*za(p3,p4)*
     &    za(p5,p8)*zb(p1,p4)*zb(p4,p8)*zb(p6,p8)*s5678**(-1) )

      d10m= + s34**(-1)*s56**(-1)*s234**(-1)*s567**(-1) * ( za(p2,p3)*
     &    za(p2,p5)*za(p5,p7)*zb(p1,p8)*zb(p2,p4)*zb(p5,p6)*s5678**(-1)
     &     - za(p2,p3)*za(p2,p7)*za(p5,p7)*zb(p1,p8)*zb(p2,p4)*zb(p6,p7
     &    )*s5678**(-1) + za(p2,p3)*za(p3,p5)*za(p5,p7)*zb(p1,p8)*zb(p3
     &    ,p4)*zb(p5,p6)*s5678**(-1) - za(p2,p3)*za(p3,p7)*za(p5,p7)*
     &    zb(p1,p8)*zb(p3,p4)*zb(p6,p7)*s5678**(-1) )
      d10m = d10m + s34**(-1)*s56**(-1)*s234**(-1)*s568**(-1) * ( za(p2
     &    ,p3)*za(p2,p7)*za(p5,p6)*zb(p1,p6)*zb(p2,p4)*zb(p6,p8)*
     &    s5678**(-1) + za(p2,p3)*za(p2,p7)*za(p5,p8)*zb(p1,p8)*zb(p2,
     &    p4)*zb(p6,p8)*s5678**(-1) + za(p2,p3)*za(p3,p7)*za(p5,p6)*zb(
     &    p1,p6)*zb(p3,p4)*zb(p6,p8)*s5678**(-1) + za(p2,p3)*za(p3,p7)*
     &    za(p5,p8)*zb(p1,p8)*zb(p3,p4)*zb(p6,p8)*s5678**(-1) )

      TWbbmZab=d9m+d10m
      return
      end
