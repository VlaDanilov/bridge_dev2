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
 
      function Faxsl(st,j1,j2,j3,j4,j5,j6,za,zb) 
      implicit none
      include 'types.f'
      complex(dp):: Faxsl
      
      integer:: j1,j2,j3,j4,j5,j6
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      character*9 st
      complex(dp):: L1
      real(dp):: t,mtsq         
      mtsq=mt**2
      if(st=='q+qb-g+g-') then
      Faxsl=
     .(2._dp*(1/(24._dp*mtsq)-L1(-t(j1,j2,j4),-s(j5,j6))/(2._dp*s(j5,j6)))*zb(
     .j1,j3)*
     .(za(j2,j5)*zb(j1,j2)+za(j4,j5)*zb(j1,j4))*zb(j3,j6))/
     .(s(j5,j6)*zb(j1,j4)*zb(j2,j4))+
     .(2._dp*(1/(24._dp*mtsq)-L1(-t(j1,j2,j3),-s(j5,j6))/(2._dp*s(j5,j6)))*za(
     .j2,j4)*za(j4,j5)*
     .(-(za(j1,j2)*zb(j1,j6))+za(j2,j3)*zb(j3,j6)))/
     .(s(j5,j6)*za(j1,j3)*za(j2,j3))
      elseif(st=='q+qb-g+g+') then
      Faxsl= 
     .(2._dp*(1/(24._dp*mtsq)-L1(-t(j1,j2,j4),-s(j5,j6))/(2._dp*s(j5,j6)))*za(
     .j2,j5)*
     .(-(za(j1,j2)*zb(j1,j3))-za(j2,j4)*zb(j3,j4))*zb(j3,j6))/
     .(s(j5,j6)*za(j1,j4)*za(j2,j4))+
     .(2._dp*(1/(24._dp*mtsq)-L1(-t(j1,j2,j3),-s(j5,j6))/(2._dp*s(j5,j6)))*za(
     .j2,j5)*
     .(-(za(j1,j2)*zb(j1,j4))+za(j2,j3)*zb(j3,j4))*zb(j4,j6))/
     .(s(j5,j6)*za(j1,j3)*za(j2,j3))
      endif
      return
      end
