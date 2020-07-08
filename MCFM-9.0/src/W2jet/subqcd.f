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
 
      subroutine subqcd(i1,i2,i3,i4,i5,i6,za,zb,amp)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
c*******************************************************************
c     the matrix elements of the  
C     helicity amplitudes for the QCD process
c     q(-p1)+qbar(-p2) --> l(p3)+abar(p4)+g(p5)+g(p6)
c     multiplied by ((a+l)^2-M**2)/(a+l)^2/g^4
c     one colour ordering only
c     left-on quark line only
c*******************************************************************
      integer:: i1,i2,i3,i4,i5,i6
      real(dp):: s156,s56,s256,s34
      complex(dp):: amp(-1:1,-1:1),bmp(-1:1,-1:1)
      complex(dp):: p1,p2,p3,b1,b2,zba2,zba3
!begin statement functions
      zba2(i1,i2,i3,i4)=zb(i1,i2)*za(i2,i4)+zb(i1,i3)*za(i3,i4)
      zba3(i1,i2,i3,i4,i5)=
     & zb(i1,i2)*za(i2,i5)+zb(i1,i3)*za(i3,i5)+zb(i1,i4)*za(i4,i5)
!end statement functions

      s156=s(i1,i5)+s(i1,i6)+s(i5,i6)
      s256=s(i2,i5)+s(i2,i6)+s(i5,i6)
      s56=s(i5,i6)
      s34=s(i3,i4)

!      amp(1,1)=four*za(i2,i3)**2
!     & /(za(i5,i6)*za(i2,i6)*za(i1,i5)*za(i4,i3))

!      b1=za(i2,i3)*zb(i2,i5)+za(i6,i3)*zb(i6,i5)
!      b2=za(i1,i6)*zb(i1,i4)+za(i5,i6)*zb(i5,i4)
!      p1=four*za(i2,i6)*zb(i2,i5)*zb(i1,i4)*b1
!     & /(zb(i2,i6)*s256*s56*s34)
!      p2=four*za(i1,i6)*zb(i1,i5)*za(i2,i3)*b2
!     & /(za(i1,i5)*s156*s56*s34)
!      p3=four*b1*b2/(zb(i2,i6)*za(i1,i5)*s56*s34)
!      amp(1,-1)=p1+p2+p3

!      b1=za(i1,i5)*zb(i1,i4)+za(i6,i5)*zb(i6,i4)
!      b2=za(i2,i3)*zb(i2,i6)+za(i3,i5)*zb(i6,i5)
!      p1=-four*zb(i1,i6)**2*za(i3,i2)*b1
!     & /(zb(i1,i5)*s156*s56*s34)
!      p2=+four*za(i2,i5)**2*zb(i1,i4)*b2
!     & /(za(i2,i6)*s256*s56*s34)
!      p3=four*zb(i1,i6)*za(i2,i5)*zb(i1,i4)*za(i2,i3)
!     & /(za(i2,i6)*zb(i1,i5)*s56*s34)
!      amp(-1,1)=p1+p2+p3

!      amp(-1,-1)=four*zb(i4,i1)**2
!     & /(zb(i5,i6)*zb(i2,i6)*zb(i1,i5)*zb(i4,i3))

      amp(-1,-1)=-four/s(i3,i4)
     & *zb(i1,i4)*zba3(i1,i2,i6,i5,i3)/zb(i1,i5)/zb(i2,i6)/zb(i6,i5)
      amp(+1,+1)=+four/s(i3,i4)
     & *za(i2,i3)*zba3(i4,i1,i6,i5,i2)/za(i1,i5)/za(i2,i6)/za(i6,i5)
      amp(+1,-1)=four/(s(i3,i4)*s(i5,i6))
     & *(zba2(i4,i1,i5,i6)*zba2(i5,i2,i6,i3)/za(i1,i5)/zb(i2,i6)
     & -za(i1,i6)*za(i2,i3)*zb(i1,i5)*zba2(i4,i1,i5,i6)/za(i1,i5)/s156
     & -za(i2,i6)*zb(i1,i4)*zb(i2,i5)*zba2(i5,i2,i6,i3)/zb(i2,i6)/s256)
      amp(-1,+1)=four/(s(i3,i4)*s(i5,i6))
     & *(za(i2,i3)*za(i2,i5)*zb(i1,i4)*zb(i1,i6)/za(i2,i6)/zb(i1,i5)
     & -za(i2,i3)*zb(i1,i6)**2*zba2(i4,i1,i6,i5)/zb(i1,i5)/s156
     & -za(i2,i5)**2*zb(i1,i4)*zba2(i6,i2,i5,i3)/za(i2,i6)/s256)


      return
      end

