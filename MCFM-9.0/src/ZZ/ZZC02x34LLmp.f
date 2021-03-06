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
 
      subroutine ZZC02x34LLmp(j1,j2,j3,j4,j5,j6,za,zb,mt,Xmp,Xpm)
      implicit none
      include 'types.f'
      
c--- Author: J. M. Campbell, October 2013
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: k1,k2,k3,k4,k5,k6
      integer:: h3,h5,j1,j2,j3,j4,j5,j6
      real(dp):: mt,mtsq,t,s134,s234
      complex(dp):: zab2,amp2,Xmp(2,2),Xpm(2,2)

C---statement functions
      zab2(k1,k2,k3,k4)=za(k1,k2)*zb(k2,k4)+za(k1,k3)*zb(k3,k4)
C---end statement functions

      s134=t(j1,j3,j4)
      s234=t(j2,j3,j4)
      k1=j1
      k2=j2
      mtsq=mt**2

      do h3=1,2
      do h5=1,2
      if (h3 == 1) then
        k3=j3
        k4=j4
      elseif (h3 == 2) then
        k3=j4
        k4=j3
      endif
      if (h5 == 1) then
        k5=j5
        k6=j6
      elseif (h5 == 2) then
        k5=j6
        k6=j5
      endif

c---- NB: amp0=0      
c      amp0=1._dp/zb(k5,k6)/zab2(k2,k3,k4,k1)**2*(
c     & -half*(s234+s(k3,k4))*za(k2,k3)**2*zb(k2,k6)**2/za(k3,k4)
c     & -s234*za(k2,k3)*zb(k2,k4)*zb(k2,k6)*zab2(k2,k3,k4,k6)
c     &  /zab2(k2,k3,k4,k2)
c     & +za(k1,k3)*zb(k2,k4)/za(k1,k2)*zab2(k2,k3,k4,k6)**2
c     & -half*s(k3,k4)*za(k2,k3)*zb(k2,k6)/za(k3,k4)*zab2(k3,k2,k4,k6)
c     & +half*za(k1,k3)*zb(k4,k6)/za(k1,k2)
c     &  *zab2(k2,k3,k4,k2)*zab2(k2,k3,k4,k6)
c     & -half*za(k2,k3)*zb(k2,k6)/za(k3,k4)
c     &  *zab2(k2,k3,k4,k2)*zab2(k3,k2,k4,k6)  
c     & -za(k2,k3)*zb(k2,k6)*za(k2,k3)*zb(k3,k4)/za(k1,k2)
c     &  *zab2(k2,k3,k4,k6)*zab2(k1,k3,k4,k2)/zab2(k2,k3,k4,k2)   
c     & +half*za(k2,k3)*zb(k4,k6)/za(k1,k2)
c     &  *zab2(k1,k3,k4,k6)*zab2(k2,k3,k4,k2)
c     & -half*zb(k2,k4)*zab2(k2,k3,k4,k6)*zab2(k3,k2,k4,k6))

      amp2=1._dp/s(k1,k2)/s(k3,k4)/s(k5,k6)*(
      
     & -two*zab2(k1,k3,k4,k2)/zab2(k2,k3,k4,k1)**3*s234
     &  *za(k2,k3)**2*za(k2,k5)*zb(k1,k6)*zb(k3,k4)*zb(k1,k2)

     & +za(k2,k3)/zab2(k2,k3,k4,k1)**2*(
     &  +zab2(k1,k5,k6,k1)*zb(k2,k4)*zab2(k2,k3,k4,k6)*zab2(k5,k3,k4,k2)
     &  +zab2(k1,k5,k6,k1)*zb(k2,k4)*s(k3,k4)*za(k2,k5)*zb(k2,k6)
     &  -za(k2,k3)*za(k1,k5)*zb(k1,k2)**2*zb(k3,k4)*zab2(k2,k3,k4,k6)
     &   *zab2(k1,k3,k4,k2)/zab2(k2,k3,k4,k2)
     &  +za(k2,k3)*za(k2,k5)*zb(k1,k6)/za(k1,k2)*zb(k3,k4)
     &   *zab2(k1,k3,k4,k2)**2
     &  +2*zb(k3,k4)*za(k2,k3)*za(k2,k5)*zb(k2,k6)*zb(k1,k2)
     &   *zab2(k1,k3,k4,k2)*s234/zab2(k2,k3,k4,k2))
   
     & +1._dp/zab2(k2,k3,k4,k1)*za(k2,k3)*(
     &  +(s134+s(k3,k4))*za(k1,k5)*zb(k2,k4)*zb(k2,k6)
     &  +za(k2,k5)*zb(k2,k4)*zb(k2,k6)*zab2(k1,k3,k4,k2)
     &  -za(k2,k3)*za(k2,k5)*zb(k2,k6)*zb(k3,k4)*zab2(k1,k3,k4,k2)**2
     &   /za(k1,k2)/zab2(k2,k3,k4,k2)
     &  +2*zb(k3,k4)*za(k3,k4)*za(k1,k2)*zb(k2,k4)*zb(k2,k6)**2
     &   *za(k5,k6)/zab2(k2,k3,k4,k2))
   
     & -za(k1,k2)*za(k3,k4)*zb(k2,k4)**2*zb(k2,k6)**2
     &   /zb(k1,k2)*za(k5,k6)/zab2(k2,k3,k4,k2) )
  
      Xmp(h3,h5)=amp2*mtsq

      enddo
      enddo

c--- obtain remaining coefficients by c.c.
      do h3=1,2
      do h5=1,2
      Xpm(h3,h5)=conjg(Xmp(3-h3,3-h5))
      enddo
      enddo

      return
      end
