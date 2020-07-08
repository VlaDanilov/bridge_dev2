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
 
      subroutine TWZggSRL(p1,p2,p3,p4,p5,p6,p7,p8,srl)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_com.f'
      include 'sprods_com.f'
C     Author: R.K. Ellis Feb, 2013
C     written by program WZggdiags.frm
C     These are the singly resonant diagrams with the
C     the Z coming off the lepton
C     Calculation is performed for LH light-line (perforce because of W)
C     Calculation is performed for LH Z-dcay line
C     The two indices of srl are the gluon helicities
      integer:: p1,p2,p3,p4,p5,p6,p7,p8
      real(dp):: s3,s4,s56,s78,s278,s178,s356,s3456
      complex(dp):: zba2,srl(2,2),iza,izb
C     statement functions
      s3(p1,p2,p3)=s(p1,p2)+s(p1,p3)+s(p2,p3)
      s4(p1,p2,p3,p4)=s(p1,p2)+s(p1,p3)+s(p1,p4)
     &               +s(p2,p3)+s(p2,p4)+s(p3,p4)
      zba2(p1,p2,p3,p4)=zb(p1,p2)*za(p2,p4)+zb(p1,p3)*za(p3,p4)
      iza(p1,p2)=cone/za(p1,p2)
      izb(p1,p2)=cone/zb(p1,p2)
C     end statement functions
      s56=s(p5,p6)
      s78=s(p7,p8)
      s356=s3(p3,p5,p6)
      s278=s3(p2,p7,p8)
      s178=s3(p1,p7,p8)
      s3456=s4(p3,p4,p5,p6)
      srl(1,1)= + s78**(-1)*s56**(-1)*s3456**(-1)*s356**(-1)*s278**(-1)
     &  * ( za(p2,p7)*za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p1,p4)*izb(p1,p8
     &    )*zba2(p6,p3,p5,p2) + za(p2,p7)*za(p3,p5)*za(p7,p8)*zb(p1,p4)
     &    *zb(p1,p7)*izb(p1,p8)*zba2(p6,p3,p5,p7) + za(p2,p7)*za(p3,p5)
     &    *za(p7,p8)*zb(p1,p4)*zb(p1,p8)*izb(p1,p8)*zba2(p6,p3,p5,p8)
     &     + za(p2,p8)*za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p1,p4)*izb(p1,
     &    p7)*zba2(p6,p3,p5,p2) + za(p2,p8)*za(p3,p5)*za(p7,p8)*zb(p1,
     &    p4)*zb(p1,p7)*izb(p1,p7)*zba2(p6,p3,p5,p7) + za(p2,p8)*za(p3,
     &    p5)*za(p7,p8)*zb(p1,p4)*zb(p1,p8)*izb(p1,p7)*zba2(p6,p3,p5,p8
     &    ) )
      srl(1,1) = srl(1,1) + s56**(-1)*s3456**(-1)*s356**(-1)*s278**(-1)
     &  * (  - za(p2,p8)*za(p3,p5)*zb(p1,p2)**2*zb(p1,p4)*izb(p1,p7)*
     &    izb(p1,p8)*izb(p2,p7)*zba2(p6,p3,p5,p2) - za(p2,p8)*za(p3,p5)
     &    *zb(p1,p2)*zb(p1,p4)*zb(p1,p7)*izb(p1,p7)*izb(p1,p8)*izb(p2,
     &    p7)*zba2(p6,p3,p5,p7) - za(p2,p8)*za(p3,p5)*zb(p1,p2)*zb(p1,
     &    p4)*zb(p1,p8)*izb(p1,p7)*izb(p1,p8)*izb(p2,p7)*zba2(p6,p3,p5,
     &    p8) - za(p3,p5)*za(p7,p8)*zb(p1,p2)*zb(p1,p4)*izb(p1,p8)*izb(
     &    p2,p7)*zba2(p6,p3,p5,p2) - za(p3,p5)*za(p7,p8)*zb(p1,p4)*zb(
     &    p1,p7)*izb(p1,p8)*izb(p2,p7)*zba2(p6,p3,p5,p7) - za(p3,p5)*
     &    za(p7,p8)*zb(p1,p4)*zb(p1,p8)*izb(p1,p8)*izb(p2,p7)*zba2(p6,
     &    p3,p5,p8) )

      srl(2,2)= + s78**(-1)*s56**(-1)*s3456**(-1)*s356**(-1)*s178**(-1)
     &  * ( za(p1,p2)*za(p3,p5)*zb(p1,p4)*zb(p1,p7)*zb(p7,p8)*iza(p2,p8
     &    )*zba2(p6,p3,p5,p2) + za(p1,p2)*za(p3,p5)*zb(p1,p4)*zb(p1,p8)
     &    *zb(p7,p8)*iza(p2,p7)*zba2(p6,p3,p5,p2) + za(p3,p5)*za(p7,p2)
     &    *zb(p1,p7)*zb(p7,p4)*zb(p7,p8)*iza(p2,p8)*zba2(p6,p3,p5,p2)
     &     + za(p3,p5)*za(p7,p2)*zb(p1,p8)*zb(p7,p4)*zb(p7,p8)*iza(p2,
     &    p7)*zba2(p6,p3,p5,p2) + za(p3,p5)*za(p8,p2)*zb(p1,p7)*zb(p7,
     &    p8)*zb(p8,p4)*iza(p2,p8)*zba2(p6,p3,p5,p2) + za(p3,p5)*za(p8,
     &    p2)*zb(p1,p8)*zb(p7,p8)*zb(p8,p4)*iza(p2,p7)*zba2(p6,p3,p5,p2
     &    ) )
      srl(2,2) = srl(2,2) + s56**(-1)*s3456**(-1)*s356**(-1)*s178**(-1)
     &  * (  - za(p1,p2)**2*za(p3,p5)*zb(p1,p4)*zb(p1,p7)*iza(p1,p8)*
     &    iza(p2,p7)*iza(p2,p8)*zba2(p6,p3,p5,p2) - za(p1,p2)*za(p3,p5)
     &    *za(p7,p2)*zb(p1,p7)*zb(p7,p4)*iza(p1,p8)*iza(p2,p7)*iza(p2,
     &    p8)*zba2(p6,p3,p5,p2) - za(p1,p2)*za(p3,p5)*za(p8,p2)*zb(p1,
     &    p7)*zb(p8,p4)*iza(p1,p8)*iza(p2,p7)*iza(p2,p8)*zba2(p6,p3,p5,
     &    p2) - za(p1,p2)*za(p3,p5)*zb(p1,p4)*zb(p7,p8)*iza(p1,p8)*iza(
     &    p2,p7)*zba2(p6,p3,p5,p2) - za(p3,p5)*za(p7,p2)*zb(p7,p4)*zb(
     &    p7,p8)*iza(p1,p8)*iza(p2,p7)*zba2(p6,p3,p5,p2) - za(p3,p5)*
     &    za(p8,p2)*zb(p7,p8)*zb(p8,p4)*iza(p1,p8)*iza(p2,p7)*zba2(p6,
     &    p3,p5,p2) )

      srl(1,2)= + s56**(-1)*s3456**(-1)*s356**(-1)*s278**(-1) * ( za(p2
     &    ,p7)*za(p3,p5)*zb(p1,p4)*zb(p2,p8)*zb(p8,p2)*iza(p7,p8)*izb(
     &    p2,p7)*izb(p7,p8)*zba2(p6,p3,p5,p2) + za(p2,p7)*za(p3,p5)*zb(
     &    p1,p4)*zb(p2,p8)*zb(p8,p7)*iza(p7,p8)*izb(p2,p7)*izb(p7,p8)*
     &    zba2(p6,p3,p5,p7) )
      srl(1,2) = srl(1,2) + s56**(-1)*s3456**(-1)*s356**(-1)*s178**(-1)
     &  * (  - za(p1,p7)**2*za(p3,p5)*zb(p1,p4)*zb(p1,p8)*iza(p1,p8)*
     &    iza(p7,p8)*izb(p7,p8)*zba2(p6,p3,p5,p2) - za(p1,p7)*za(p3,p5)
     &    *za(p8,p7)*zb(p1,p8)*zb(p8,p4)*iza(p1,p8)*iza(p7,p8)*izb(p7,
     &    p8)*zba2(p6,p3,p5,p2) )
      srl(1,2) = srl(1,2) + s56**(-1)*s3456**(-1)*s356**(-1) * (  - za(
     &    p1,p7)*za(p3,p5)*zb(p1,p4)*zb(p2,p8)*iza(p1,p8)*iza(p7,p8)*
     &    izb(p2,p7)*izb(p7,p8)*zba2(p6,p3,p5,p2) - za(p1,p7)*za(p3,p5)
     &    *zb(p1,p4)*iza(p1,p8)*iza(p7,p8)*izb(p2,p7)*zba2(p6,p3,p5,p7)
     &     + za(p3,p5)*zb(p2,p8)*zb(p8,p4)*iza(p1,p8)*izb(p2,p7)*izb(p7
     &    ,p8)*zba2(p6,p3,p5,p2) + za(p3,p5)*zb(p8,p4)*iza(p1,p8)*izb(
     &    p2,p7)*zba2(p6,p3,p5,p7) )

      srl(2,1)= + s56**(-1)*s3456**(-1)*s356**(-1)*s278**(-1) * ( za(p2
     &    ,p8)**2*za(p3,p5)*zb(p1,p4)*zb(p7,p2)*iza(p2,p7)*iza(p7,p8)*
     &    izb(p7,p8)*zba2(p6,p3,p5,p2) + za(p2,p8)**2*za(p3,p5)*zb(p1,
     &    p4)*zb(p7,p8)*iza(p2,p7)*iza(p7,p8)*izb(p7,p8)*zba2(p6,p3,p5,
     &    p8) )
      srl(2,1) = srl(2,1) + s56**(-1)*s3456**(-1)*s356**(-1)*s178**(-1)
     &  * (  - za(p1,p8)*za(p3,p5)*zb(p1,p4)*zb(p1,p7)**2*iza(p7,p8)*
     &    izb(p1,p8)*izb(p7,p8)*zba2(p6,p3,p5,p2) - za(p3,p5)*za(p7,p8)
     &    *zb(p1,p7)**2*zb(p7,p4)*iza(p7,p8)*izb(p1,p8)*izb(p7,p8)*
     &    zba2(p6,p3,p5,p2) )
      srl(2,1) = srl(2,1) + s56**(-1)*s3456**(-1)*s356**(-1) * (  - za(
     &    p2,p8)*za(p3,p5)*zb(p1,p4)*zb(p1,p7)*iza(p2,p7)*iza(p7,p8)*
     &    izb(p1,p8)*izb(p7,p8)*zba2(p6,p3,p5,p2) )

      return
      end
