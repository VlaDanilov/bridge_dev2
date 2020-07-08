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
 
      function A51mmppp(j1,j2,j3,j4,j5,za,zb)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'epinv.f'
      include 'scale.f'
      complex(dp):: A51mmppp
      complex(dp):: Vf,Vs,Ff,Fs,L0,L2,miA5tree,lnrat
      integer:: j1,j2,j3,j4,j5
C----Eq.(8,10) of hep-ph/9302280v1 of BDK multiplied by 16*pi^2*(-i)
C--- to give (16*pi^2)*(-i)*A^{[1/2]}_{5;1}
      miA5tree=za(j1,j2)**4
     & /(za(j1,j2)*za(j2,j3)*za(j3,j4)*za(j4,j5)*za(j5,j1))
      Vf=cplx1(-five/two*epinv-two)
     & -half*(lnrat(musq,-s(j2,j3))+lnrat(musq,-s(j5,j1)))
      Vs=-Vf/three+cplx1(two/nine)
      Ff=-za(j1,j2)**2*L0(-s(j2,j3),-s(j5,j1))
     & *(za(j2,j3)*zb(j3,j4)*za(j4,j1)+za(j2,j4)*zb(j4,j5)*za(j5,j1))
     & /(two*za(j2,j3)*za(j3,j4)*za(j4,j5)*za(j5,j1)*s(j5,j1))
      Fs=-Ff/three
     & -zb(j3,j4)*za(j4,j1)*za(j2,j4)*zb(j4,j5)*L2(-s(j2,j3),-s(j5,j1))
     & *(za(j2,j3)*zb(j3,j4)*za(j4,j1)+za(j2,j4)*zb(j4,j5)*za(j5,j1))
     & /(three*za(j3,j4)*za(j4,j5)*s(j5,j1)**3)
     & -za(j3,j5)*zb(j3,j5)**3
     & /(three*zb(j1,j2)*zb(j2,j3)*za(j3,j4)*za(j4,j5)*zb(j5,j1))
     & +za(j1,j2)*zb(j3,j5)**2
     & /(three*zb(j2,j3)*za(j3,j4)*za(j4,j5)*zb(j5,j1))
     & +za(j1,j2)*zb(j3,j4)*za(j4,j1)*za(j2,j4)*zb(j4,j5)
     & /(six*s(j2,j3)*za(j3,j4)*za(j4,j5)*s(j5,j1))
      A51mmppp=-(Vf+Vs)*miA5tree-(Ff+Fs)
      
c      write(6,*) 'j1,j2,j3,j4,j5,miA5tree',j1,j2,j3,j4,j5,miA5tree
      return
      end 
