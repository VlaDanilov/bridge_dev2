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
 
      subroutine qqbZggtree(j1,j2,j3,j4,j5,j6,za,zb,ampTree)
      implicit none
      include 'types.f'
c--- Tree-level amplitudes for the process
c----    0 -> q(j1) g(j2) g(j3) qb(j4) e+(j5) e-(j6) from
c--- hep-ph/9708239, Eqs(8.4), (8.8), (8.14) multiplied by (-i) 
      integer:: j1,j2,j3,j4,j5,j6
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      real(dp):: t
      complex(dp):: ampTree(2,2)

!      if(st=='q+g-g+qb-') then
      ampTree(1,2)= 
     .((za(j2,j4)*za(j4,j5)*zb(j1,j3)*zb(j1,j6))/
     .(s(j2,j3)*s(j5,j6)*za(j3,j4)*zb(j1,j2))-
     .(za(j4,j5)*zb(j1,j3)**2*(-(za(j1,j2)*zb(j1,j6))+za(j2,j3)*zb(j3,j6
     .)))/
     .(s(j2,j3)*s(j5,j6)*zb(j1,j2)*t(j1,j2,j3))+
     .(za(j2,j4)**2*zb(j1,j6)*(-(za(j2,j5)*zb(j2,j3))+za(j4,j5)*zb(j3,j4
     .)))/
     .(s(j2,j3)*s(j5,j6)*za(j3,j4)*t(j2,j3,j4)))

!      elseif(st=='q+g+g-qb-') then
      ampTree(2,1)= 
     .(-(((za(j3,j5)*zb(j2,j3)+za(j4,j5)*zb(j2,j4))*
     .(-(za(j1,j3)*zb(j1,j6))-za(j2,j3)*zb(j2,j6)))/
     .(s(j2,j3)*s(j5,j6)*za(j1,j2)*zb(j3,j4)))-
     .(za(j1,j3)*za(j4,j5)*zb(j1,j2)*
     .(-(za(j1,j3)*zb(j1,j6))-za(j2,j3)*zb(j2,j6)))/
     .(s(j2,j3)*s(j5,j6)*za(j1,j2)*t(j1,j2,j3))+
     .(za(j3,j4)*zb(j1,j6)*zb(j2,j4)*
     .(za(j3,j5)*zb(j2,j3)+za(j4,j5)*zb(j2,j4)))/
     .(s(j2,j3)*s(j5,j6)*zb(j3,j4)*t(j2,j3,j4)))

!      elseif(st=='q+g+g+qb-') then
c---This amplitude corresponds to
c    q(4)+l(5) --> q_R(1)+l_R(6)+g_R(2)+g_R(3)
      ampTree(2,2)= 
     .(-za(j4,j5)**2)/(za(j1,j2)*za(j2,j3)*za(j3,j4)*za(j5,j6))

 
!      if(st=='q+g-g-qb-') then
!        a6treeg1=a6treeg('q+g+g+qb-',j4,j3,j2,j1,j6,j5,zb,za)
      ampTree(1,1)= 
     .(-zb(j1,j6)**2)/(zb(j4,j3)*zb(j3,j2)*zb(j2,j1)*zb(j6,j5))

      return
      end
