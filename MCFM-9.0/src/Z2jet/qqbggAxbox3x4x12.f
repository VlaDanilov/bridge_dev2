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
 
      subroutine qqbggAxbox3x4x12(p1,p2,p3,p4,p5,p6,za,zb,coeff0,coeff2)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer::p1,p2,p3,p4,p5,p6
      complex(dp)::coeff0(2,2),coeff2(2,2),bit0(2,2),bit2(2,2)
      
      call qqbggAxbox3x4x12sub(p1,p2,p3,p4,p5,p6,za,zb,coeff0,coeff2)
c---- flip5 operation is this exchange and a minus sign
      call qqbggAxbox3x4x12sub(p2,p1,p3,p4,p6,p5,zb,za,bit0,bit2)
      coeff0(1,1)=-bit0(2,2)
      coeff0(2,1)=-bit0(1,2)
      coeff2(1,1)=-bit2(2,2)
      coeff2(2,1)=-bit2(1,2)
      
      return
      end
      
      subroutine qqbggAxbox3x4x12sub(p1,p2,p3,p4,p5,p6,za,zb,coeff0,coeff2)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer::p1,p2,p3,p4,p5,p6
      complex(dp)::coeff0(2,2),coeff2(2,2),zab2
      real(dp)::s3,s12,s34,s56,s124

C--- begin statement function
      s3(p1,p2,p3)=s(p1,p2)+s(p2,p3)+s(p3,p1)
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
C--- end statement functions

      s12=s(p1,p2)
      s34=s(p3,p4)
      s56=s(p5,p6)
      s124=s3(p1,p2,p4)

      coeff0=czip
      coeff2=czip

C     Indices are the polarizations of the gluons, 1^+  6^+ are fixed 
      coeff0(2,2) =  0  
      coeff2(2,2) =  
     & -za(p2,p5)*zb(p1,p2)*zb(p3,p4)/(2*za(p3,p4)*s12*s56)*(
     &  s124*za(p2,p3)*zb(p6,p4)/(zab2(p3,p1,p2,p4)) 
     & +s124*za(p2,p4)*zb(p6,p3)/(zab2(p4,p1,p2,p3))
     & + zab2(p2,p1,p4,p6)  - za(p2,p3)*zb(p3,p6) )

      coeff0(1,2) =   -1._dp/4._dp*s34*s124*( 
     &(zab2(p2,p1,p4,p3)**2*zab2(p4,p1,p2,p6)**2-za(p2,p4)**2*zb(p3,p6)**2*s124**2)
     &  /(za(p1,p2)*zb(p5,p6)*zab2(p4,p1,p2,p3)**4) ) 
      coeff2(1,2) =   
     &-zb(p1,p4)*s124/(2*zab2(p4,p1,p2,p3)*s12*s56)*( 
     &  za(p1,p2)*za(p3,p5)*zb(p1,p6) + za(p2,p3)*za(p5,p4)*zb(p4,p6) )  
     &  - 3*(s124*zb(p1,p3)*zab2(p3,p1,p2,p4)*za(p4,p5)*zab2(p2,p1,p4,p6) )/(2*zab2(p4,p1,p2,p3)**2*s12*s56)  
     &  - (zab2(p3,p1,p2,p4))/(2*zab2(p4,p1,p2,p3)*s12*s56)* (
     &      zab2(p2,p1,p4,p6)*(za(p3,p5)*zb(p1,p3)-za(p4,p5)*zb(p1,p4))
     &    - za(p2,p5)*zb(p1,p6)*s124    
     &    + za(p2,p4)*za(p3,p5)*zb(p1,p4)*zb(p3,p6) )
     & + za(p2,p3)*za(p3,p5)*zb(p1,p4)*zb(p4,p6)/(2*s12*s56)

      return
      end
