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
 
      function Fax(st,j1,j2,j3,j4,j5,j6,za,zb) 
      implicit none
      include 'types.f'
      complex(dp):: Fax
      
      integer:: j1,j2,j3,j4,j5,j6
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      character*9 st
      complex(dp):: L0,L1,Lsm1_2mh,Lsm1_2me,I3m,Lnrat
      real(dp):: t,mtsq  
      mtsq=mt**2

      if(st=='q+qb-g-g+') then
      Fax=
     .(za(j2,j3)**2*za(j3,j5)*zb(j1,j6))/(12._dp*mtsq*s(j5,j6)*za(j2,j4)*z
     .a(j3,j4))-
     .(Lnrat(-s(j5,j6),-s(j3,j4))*(za(j2,j5)*zb(j1,j2)+za(j3,j5)*zb(j1,j
     .3))**2)/
     .(za(j5,j6)*zb(j1,j2)*(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))*
     .*2)-
     .(za(j2,j3)**2*za(j4,j5)*zb(j1,j6))/
     .(s(j5,j6)*za(j2,j4)*za(j3,j4)*
     .(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3)))-
     .(L0(-t(j1,j2,j3),-s(j1,j2))*za(j2,j3)*
     .(za(j2,j5)*zb(j1,j2)+za(j3,j5)*zb(j1,j3))*
     .(-(za(j1,j5)*zb(j1,j3))-za(j2,j5)*zb(j2,j3)))/
     .(s(j1,j2)*za(j5,j6)*(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))**
     .2)+
     .(L0(-t(j1,j2,j3),-s(j5,j6))*(za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb(j1,j
     .3))*
     .(za(j2,j5)*zb(j1,j2)+za(j3,j5)*zb(j1,j3))*zb(j1,j4)*
     .(-(za(j1,j5)*zb(j1,j3))-za(j2,j5)*zb(j2,j3)))/
     .(s(j5,j6)*za(j5,j6)*zb(j1,j2)*zb(j1,j3)*
     .(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))**2)+
     .(za(j2,j5)*zb(j1,j4)**2*zb(j3,j6))/
     .(s(j5,j6)*zb(j1,j3)*(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))*
     .zb(j3,j4))-(za(j3,j5)*zb(j1,j4)*
     .(((s(j1,j2)-s(j3,j4)-s(j5,j6))*za(j2,j3)*za(j4,j5))/
     .((s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-
     .2._dp*s(j1,j2)*s(j5,j6)-2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*za(j3,j4)
     .*
     .za(j5,j6))+((-s(j1,j2)+s(j3,j4)-s(j5,j6))*za(j4,j5)*
     .zb(j1,j4))/
     .((s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-
     .2._dp*s(j1,j2)*s(j5,j6)-2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*za(j5,j6)
     .*
     .zb(j1,j2))-(2._dp*za(j2,j3)*zb(j3,j6))/
     .(s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s(j5,j
     .6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)-
     .(zb(j1,j4)*zb(j3,j6))/(s(j5,j6)*zb(j1,j2)*zb(j3,j4))+
     .((-s(j1,j2)-s(j3,j4)+s(j5,j6))*zb(j1,j4)*zb(j3,j6))/
     .((s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s(j5,
     .j6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*zb(j1,j2)*zb(j3,j4))))/
     .(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))-
     .(za(j2,j5)*zb(j1,j4)**2*zb(j4,j6))/(12._dp*mtsq*s(j5,j6)*zb(j1,j3)*z
     .b(j3,j4))-
     .(za(j2,j3)*zb(j4,j6)*(-((za(j2,j3)*za(j4,j5))/
     .(s(j5,j6)*za(j1,j2)*za(j3,j4)))+
     .((-s(j1,j2)-s(j3,j4)+s(j5,j6))*za(j2,j3)*za(j4,j5))/
     .((s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-
     .2._dp*s(j1,j2)*s(j5,j6)-2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*za(j1,j2)
     .*
     .za(j3,j4))-(2._dp*za(j4,j5)*zb(j1,j4))/
     .(s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s(j5,j
     .6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)+
     .((-s(j1,j2)+s(j3,j4)-s(j5,j6))*za(j2,j3)*zb(j3,j6))/
     .((s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-
     .2._dp*s(j1,j2)*s(j5,j6)-2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*za(j1,j2)
     .*
     .zb(j5,j6))+((s(j1,j2)-s(j3,j4)-s(j5,j6))*zb(j1,j4)*zb(j3,j6))/
     .((s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s(j5,
     .j6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*zb(j3,j4)*zb(j5,j6))))/
     .(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))+
     .(L0(-t(j1,j2,j4),-s(j1,j2))*zb(j1,j4)*
     .(-(za(j1,j4)*zb(j1,j6))-za(j2,j4)*zb(j2,j6))*
     .(-(za(j1,j2)*zb(j1,j6))+za(j2,j4)*zb(j4,j6)))/
     .(s(j1,j2)*(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))**2*zb(j5,j6
     .))+
     .(L0(-t(j1,j2,j4),-s(j5,j6))*za(j2,j3)*
     .(-(za(j1,j4)*zb(j1,j6))-za(j2,j4)*zb(j2,j6))*
     .(-(za(j1,j2)*zb(j1,j3))-za(j2,j4)*zb(j3,j4))*
     .(-(za(j1,j2)*zb(j1,j6))+za(j2,j4)*zb(j4,j6)))/
     .(s(j5,j6)*za(j1,j2)*za(j2,j4)*
     .(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))**2*zb(j5,j6))-
     .(Lnrat(-s(j5,j6),-s(j3,j4))*(-(za(j1,j2)*zb(j1,j6))+za(j2,j4)*zb(j
     .4,j6))**2)/
     .(za(j1,j2)*(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))**2*zb(j5,j
     .6))+
     .(L1(-s(j5,j6),-t(j1,j2,j3))*za(j4,j5)*
     .(za(j2,j5)*zb(j1,j2)+za(j3,j5)*zb(j1,j3))*zb(j1,j4)**2)/
     .(za(j5,j6)*zb(j1,j2)*zb(j1,j3)*
     .(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))*t(j1,j2,j3))+
     .(Lsm1_2mh(s(j3,j4),t(j1,j2,j3),s(j1,j2),s(j5,j6))*
     .((za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb(j1,j3))**2*
     .(-(za(j1,j5)*zb(j1,j3))-za(j2,j5)*zb(j2,j3))**2-
     .za(j4,j5)**2*zb(j1,j3)**2*t(j1,j2,j3)**2))/
     .(2._dp*za(j5,j6)*zb(j1,j2)*(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j
     .3))**4)


      Fax=Fax
     .-Lnrat(-s(j1,j2),-s(j3,j4))*((za(j2,j3)*za(j3,j5)**2*zb(j1,j4))/
     .((s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s(j5,
     .j6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*za(j3,j4)*za(j5,j6))+
     .(za(j2,j3)*za(j4,j5)*(za(j2,j5)*zb(j1,j2)+za(j3,j5)*zb(j1,j3)))/
     .(za(j3,j4)*za(j5,j6)*(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))*
     .*
     .2)+(6._dp*za(j1,j2)*(za(j2,j5)*zb(j1,j2)+za(j3,j5)*zb(j1,j3))*
     .(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4))*
     .((-s(j1,j2)+s(j3,j4)-s(j5,j6))*zb(j1,j6)-
     .2._dp*za(j2,j5)*zb(j1,j2)*zb(j5,j6)))/
     .((s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s(j5,
     .j6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)**2*
     .(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3)))-
     .(za(j2,j4)*(za(j2,j5)*zb(j1,j2)+za(j3,j5)*zb(j1,j3))*
     .(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4))*
     .(3._dp*(-(za(j1,j4)*za(j3,j5)*zb(j1,j3))-
     .za(j2,j4)*za(j3,j5)*zb(j2,j3))-
     .za(j4,j5)*(t(j1,j2,j3)-t(j1,j2,j4))))/
     .((s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s(j5,
     .j6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*za(j3,j4)*za(j5,j6)*
     .(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))**2))


      Fax=Fax
     .-I3m(s(j1,j2),s(j3,j4),s(j5,j6))*
     .(-((za(j2,j3)*za(j3,j5)*zb(j1,j4)*zb(j4,j6))/
     .(s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s(j5,j
     .6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2))-
     .(3._dp*(-s(j1,j2)+s(j3,j4)-s(j5,j6))*
     .(za(j2,j5)*zb(j1,j2)+za(j3,j5)*zb(j1,j3))*
     .(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4))*
     .(-((s(j1,j2)-s(j3,j4)-s(j5,j6))*za(j1,j2)*zb(j1,j6))-
     .(-s(j1,j2)-s(j3,j4)+s(j5,j6))*za(j2,j5)*zb(j5,j6)))/
     .((s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s(j5,
     .j6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)**2*
     .(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3)))-
     .(3._dp*(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4))*
     .(-(za(j1,j2)*za(j2,j5)*zb(j1,j2)*zb(j1,j6))-
     .za(j2,j3)*za(j4,j5)*zb(j1,j4)*zb(j3,j6)-
     .za(j2,j4)*za(j3,j5)*zb(j1,j3)*zb(j4,j6)-
     .za(j2,j5)*za(j5,j6)*zb(j1,j6)*zb(j5,j6)))/
     .(2._dp*(s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s
     .(j5,j6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*
     .(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3)))-
     .(za(j2,j3)*(za(j2,j5)*zb(j1,j2)+za(j3,j5)*zb(j1,j3))*zb(j4,j6))/
     .(2._dp*(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))*t(j1,j2,j3))+
     .(za(j2,j4)*(za(j2,j5)*zb(j1,j2)+za(j3,j5)*zb(j1,j3))*
     .(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4))*zb(j3,j6)*
     .(t(j1,j2,j3)-t(j1,j2,j4)))/
     .((s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s(j5,
     .j6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*
     .(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))**2))


      Fax=Fax
     .+(L1(-s(j5,j6),-t(j1,j2,j4))*za(j2,j3)**2*zb(j3,j6)*
     .(-(za(j1,j2)*zb(j1,j6))+za(j2,j4)*zb(j4,j6)))/
     .(za(j1,j2)*za(j2,j4)*(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))*
     .zb(j5,j6)*t(j1,j2,j4))


      Fax=Fax
     .+(Lsm1_2mh(s(j3,j4),t(j1,j2,j4),s(j1,j2),
     .s(j5,j6))*((-(za(j1,j4)*zb(j1,j6))-za(j2,j4)*zb(j2,j6))**2*
     .(-(za(j1,j2)*zb(j1,j3))-za(j2,j4)*zb(j3,j4))**2-
     .za(j2,j4)**2*zb(j3,j6)**2*t(j1,j2,j4)**2))/
     .(2._dp*za(j1,j2)*(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))**4*zb(
     .j5,j6))


      Fax=Fax
     .-I3m(s(j1,j2),s(j3,j4),s(j5,j6))*
     .(-((za(j2,j3)*za(j3,j5)*zb(j1,j4)*zb(j4,j6))/
     .(s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s(j5,j
     .6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2))-
     .(3._dp*(-s(j1,j2)+s(j3,j4)-s(j5,j6))*
     .((s(j1,j2)-s(j3,j4)-s(j5,j6))*za(j2,j5)*zb(j1,j2)+
     .(-s(j1,j2)-s(j3,j4)+s(j5,j6))*za(j5,j6)*zb(j1,j6))*
     .(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4))*
     .(-(za(j1,j2)*zb(j1,j6))+za(j2,j4)*zb(j4,j6)))/
     .((s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s(j5,
     .j6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)**2*
     .(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3)))-
     .(3._dp*(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4))*
     .(-(za(j1,j2)*za(j2,j5)*zb(j1,j2)*zb(j1,j6))-
     .za(j2,j3)*za(j4,j5)*zb(j1,j4)*zb(j3,j6)-
     .za(j2,j4)*za(j3,j5)*zb(j1,j3)*zb(j4,j6)-
     .za(j2,j5)*za(j5,j6)*zb(j1,j6)*zb(j5,j6)))/
     .(2._dp*(s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s
     .(j5,j6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*
     .(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3)))-
     .(za(j3,j5)*zb(j1,j4)*(-(za(j1,j2)*zb(j1,j6))+za(j2,j4)*zb(j4,j6)))
     ./
     .(2._dp*(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))*t(j1,j2,j4))+
     .(za(j4,j5)*zb(j1,j3)*(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4))*
     .(-(za(j1,j2)*zb(j1,j6))+za(j2,j4)*zb(j4,j6))*
     .(-t(j1,j2,j3)+t(j1,j2,j4)))/
     .((s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s(j5,
     .j6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*
     .(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))**2))

      
      Fax=Fax
     .-Lnrat(-s(j1,j2),-s(j3,j4))*((-6._dp*zb(j1,j2)*
     .((-s(j1,j2)+s(j3,j4)-s(j5,j6))*za(j2,j5)-
     .2._dp*za(j1,j2)*za(j5,j6)*zb(j1,j6))*
     .(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4))*
     .(-(za(j1,j2)*zb(j1,j6))+za(j2,j4)*zb(j4,j6)))/
     .((s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s(j5,
     .j6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)**2*
     .(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3)))+
     .(za(j2,j3)*zb(j1,j4)*zb(j4,j6)**2)/
     .((s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s(j5,
     .j6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*zb(j3,j4)*zb(j5,j6))+
     .(zb(j1,j4)*zb(j3,j6)*(-(za(j1,j2)*zb(j1,j6))+za(j2,j4)*zb(j4,j6)))
     ./
     .((-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))**2*zb(j3,j4)*
     .zb(j5,j6))-(zb(j1,j3)*
     .(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4))*
     .(-(za(j1,j2)*zb(j1,j6))+za(j2,j4)*zb(j4,j6))*
     .(3._dp*(-(za(j1,j4)*zb(j1,j3)*zb(j4,j6))-
     .za(j2,j4)*zb(j2,j3)*zb(j4,j6))-
     .zb(j3,j6)*(-t(j1,j2,j3)+t(j1,j2,j4))))/
     .((s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s(j5,
     .j6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*
     .(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))**2*zb(j3,j4)*zb(j5,j6
     .)
     .))


      Fax=Fax
     .-Lnrat(-s(j5,j6),-s(j3,j4))*(-((za(j2,j3)**2*za(j3,j5)*zb(j4,j6)
     .)/
     .((s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s(j5,
     .j6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*za(j1,j2)*za(j3,j4)))-
     .(za(j2,j4)*za(j3,j5)*(za(j2,j3)*zb(j3,j6)+za(j2,j5)*zb(j5,j6)))/
     .(za(j1,j2)*za(j3,j4)*(-(za(j4,j5)*zb(j3,j5))-za(j4,j6)*zb(j3,j6))*
     .*
     .2)-(6._dp*za(j5,j6)*(-(za(j3,j5)*zb(j4,j5))-za(j3,j6)*zb(j4,j6))*
     .(za(j2,j3)*zb(j3,j6)+za(j2,j5)*zb(j5,j6))*
     .(-((-s(j1,j2)+s(j3,j4)-s(j5,j6))*zb(j1,j6))+
     .2._dp*za(j2,j5)*zb(j1,j2)*zb(j5,j6)))/
     .((s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s(j5,
     .j6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)**2*
     .(-(za(j4,j5)*zb(j3,j5))-za(j4,j6)*zb(j3,j6)))-
     .(za(j4,j5)*(-(za(j3,j5)*zb(j4,j5))-za(j3,j6)*zb(j4,j6))*
     .(za(j2,j3)*zb(j3,j6)+za(j2,j5)*zb(j5,j6))*
     .(3._dp*(za(j2,j3)*za(j4,j5)*zb(j3,j5)+
     .za(j2,j3)*za(j4,j6)*zb(j3,j6))+
     .za(j2,j4)*(t(j3,j5,j6)-t(j4,j5,j6))))/
     .((s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s(j5,
     .j6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*za(j1,j2)*za(j3,j4)*
     .(-(za(j4,j5)*zb(j3,j5))-za(j4,j6)*zb(j3,j6))**2))


      Fax=Fax
     .-Lnrat(-s(j5,j6),-s(j3,j4))*(-((za(j3,j5)*zb(j1,j4)**2*zb(j4,j6))/
     .((s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s(j5,
     .j6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*zb(j1,j2)*zb(j3,j4)))-
     .(zb(j1,j3)*(za(j4,j5)*zb(j1,j4)-za(j5,j6)*zb(j1,j6))*zb(j4,j6))/
     .(zb(j1,j2)*zb(j3,j4)*(-(za(j4,j5)*zb(j3,j5))-za(j4,j6)*zb(j3,j6))*
     .*2)
     .+(6._dp*(za(j4,j5)*zb(j1,j4)-za(j5,j6)*zb(j1,j6))*
     .(-((-s(j1,j2)+s(j3,j4)-s(j5,j6))*za(j2,j5))+
     .2._dp*za(j1,j2)*za(j5,j6)*zb(j1,j6))*
     .(-(za(j3,j5)*zb(j4,j5))-za(j3,j6)*zb(j4,j6))*zb(j5,j6))/
     .((s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s(j5,
     .j6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)**2*
     .(-(za(j4,j5)*zb(j3,j5))-za(j4,j6)*zb(j3,j6)))-
     .((za(j4,j5)*zb(j1,j4)-za(j5,j6)*zb(j1,j6))*zb(j3,j6)*
     .(-(za(j3,j5)*zb(j4,j5))-za(j3,j6)*zb(j4,j6))*
     .(3._dp*(za(j4,j5)*zb(j1,j4)*zb(j3,j5)+
     .za(j4,j6)*zb(j1,j4)*zb(j3,j6))+
     .zb(j1,j3)*(-t(j3,j5,j6)+t(j4,j5,j6))))/
     .((s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s(j5,
     .j6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*zb(j1,j2)*zb(j3,j4)*
     .(-(za(j4,j5)*zb(j3,j5))-za(j4,j6)*zb(j3,j6))**2))


      elseif(st=='q+qb-g+g-') then
      Fax= 
     .(Lnrat(-s(j5,j6),-s(j3,j4))*(za(j2,j5)*zb(j1,j2)+za(j4,j5)*zb(j1,j
     .4))**2)/
     .(za(j5,j6)*zb(j1,j2)*(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4))*
     .*2)+
     .(L0(-t(j1,j2,j4),-s(j1,j2))*za(j2,j4)*
     .(za(j2,j5)*zb(j1,j2)+za(j4,j5)*zb(j1,j4))*
     .(-(za(j1,j5)*zb(j1,j4))-za(j2,j5)*zb(j2,j4)))/
     .(s(j1,j2)*za(j5,j6)*(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4))**
     .2)-
     .(L0(-t(j1,j2,j4),-s(j5,j6))*zb(j1,j3)*
     .(za(j2,j5)*zb(j1,j2)+za(j4,j5)*zb(j1,j4))*
     .(-(za(j1,j5)*zb(j1,j4))-za(j2,j5)*zb(j2,j4))*
     .(-(za(j1,j3)*zb(j1,j2))-za(j3,j4)*zb(j2,j4)))/
     .(s(j5,j6)*za(j5,j6)*zb(j1,j2)*zb(j2,j4)*
     .(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4))**2)+
     .(zb(j1,j3)*(-(za(j2,j5)*zb(j2,j3))+za(j4,j5)*zb(j3,j4))*zb(j3,j6))
     ./
     .(12._dp*mtsq*s(j5,j6)*zb(j2,j4)*zb(j3,j4))-
     .(za(j2,j4)*za(j4,j5)*(-(za(j1,j4)*zb(j1,j6))-za(j3,j4)*zb(j3,j6)))
     ./
     .(12._dp*mtsq*s(j5,j6)*za(j1,j3)*za(j3,j4))+
     .(za(j2,j4)*za(j3,j5)*(-(za(j1,j4)*zb(j1,j6))-za(j3,j4)*zb(j3,j6)))
     ./
     .(s(j5,j6)*za(j1,j3)*za(j3,j4)*
     .(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4)))-
     .(zb(j1,j3)*(-(za(j2,j5)*zb(j2,j3))+za(j4,j5)*zb(j3,j4))*zb(j4,j6))
     ./
     .(s(j5,j6)*zb(j2,j4)*(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4))*
     .zb(j3,j4))+(za(j4,j5)*zb(j1,j3)*
     .(-(((s(j1,j2)-s(j3,j4)-s(j5,j6))*za(j2,j4)*za(j3,j5))/
     .((s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-
     .2._dp*s(j1,j2)*s(j5,j6)-2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*
     .za(j3,j4)*za(j5,j6)))+
     .((-s(j1,j2)+s(j3,j4)-s(j5,j6))*za(j3,j5)*zb(j1,j3))/
     .((s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-
     .2._dp*s(j1,j2)*s(j5,j6)-2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*za(j5,j6)
     .*
     .zb(j1,j2))-(2._dp*za(j2,j4)*zb(j4,j6))/
     .(s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s(j5,j
     .6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)+
     .(zb(j1,j3)*zb(j4,j6))/(s(j5,j6)*zb(j1,j2)*zb(j3,j4))-
     .((-s(j1,j2)-s(j3,j4)+s(j5,j6))*zb(j1,j3)*zb(j4,j6))/
     .((s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s(j5,
     .j6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*zb(j1,j2)*zb(j3,j4))))/
     .(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4))+
     .(za(j2,j4)*zb(j3,j6)*((za(j2,j4)*za(j3,j5))/
     .(s(j5,j6)*za(j1,j2)*za(j3,j4))-
     .((-s(j1,j2)-s(j3,j4)+s(j5,j6))*za(j2,j4)*za(j3,j5))/
     .((s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-
     .2._dp*s(j1,j2)*s(j5,j6)-2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*za(j1,j2)
     .*
     .za(j3,j4))-(2._dp*za(j3,j5)*zb(j1,j3))/
     .(s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s(j5,j
     .6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)+
     .((-s(j1,j2)+s(j3,j4)-s(j5,j6))*za(j2,j4)*zb(j4,j6))/
     .((s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-
     .2._dp*s(j1,j2)*s(j5,j6)-2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*za(j1,j2)
     .*
     .zb(j5,j6))-((s(j1,j2)-s(j3,j4)-s(j5,j6))*zb(j1,j3)*zb(j4,j6))/
     .((s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s(j5,
     .j6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*zb(j3,j4)*zb(j5,j6))))/
     .(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4))-
     .(L0(-t(j1,j2,j3),-s(j1,j2))*zb(j1,j3)*
     .(-(za(j1,j3)*zb(j1,j6))-za(j2,j3)*zb(j2,j6))*
     .(-(za(j1,j2)*zb(j1,j6))+za(j2,j3)*zb(j3,j6)))/
     .(s(j1,j2)*(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4))**2*zb(j5,j6
     .))-
     .(L0(-t(j1,j2,j3),-s(j5,j6))*za(j2,j4)*
     .(-(za(j1,j3)*zb(j1,j6))-za(j2,j3)*zb(j2,j6))*
     .(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))*
     .(-(za(j1,j2)*zb(j1,j6))+za(j2,j3)*zb(j3,j6)))/
     .(s(j5,j6)*za(j1,j2)*za(j1,j3)*
     .(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4))**2*zb(j5,j6))+
     .(Lnrat(-s(j5,j6),-s(j3,j4))*(-(za(j1,j2)*zb(j1,j6))+za(j2,j3)*zb(j
     .3,j6))**2)/
     .(za(j1,j2)*(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4))**2*zb(j5,j
     .6))-
     .(L1(-s(j5,j6),-t(j1,j2,j3))*za(j1,j4)*za(j2,j4)*
     .(-(za(j1,j2)*zb(j1,j6))+za(j2,j3)*zb(j3,j6))*zb(j4,j6))/
     .(za(j1,j2)*za(j1,j3)*(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4))*
     .zb(j5,j6)*t(j1,j2,j3))-(Lsm1_2mh(s(j3,j4),t(j1,j2,j3),s(j1,j2),
     .s(j5,j6))*((-(za(j1,j3)*zb(j1,j6))-za(j2,j3)*zb(j2,j6))**2*
     .(-(za(j1,j2)*zb(j1,j4))+za(j2,j3)*zb(j3,j4))**2-
     .za(j2,j3)**2*zb(j4,j6)**2*t(j1,j2,j3)**2))/
     .(2._dp*za(j1,j2)*(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4))**4*zb(
     .j5,j6))


      Fax=Fax
     .+Lnrat(-s(j1,j2),-s(j3,j4))*((-6._dp*zb(j1,j2)*
     .((-s(j1,j2)+s(j3,j4)-s(j5,j6))*za(j2,j5)-
     .2._dp*za(j1,j2)*za(j5,j6)*zb(j1,j6))*
     .(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))*
     .(-(za(j1,j2)*zb(j1,j6))+za(j2,j3)*zb(j3,j6)))/
     .((s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s(j5,
     .j6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)**2*
     .(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4)))-
     .(za(j2,j4)*zb(j1,j3)*zb(j3,j6)**2)/
     .((s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s(j5,
     .j6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*zb(j3,j4)*zb(j5,j6))-
     .(zb(j1,j3)*(-(za(j1,j2)*zb(j1,j6))+za(j2,j3)*zb(j3,j6))*zb(j4,j6))
     ./
     .((-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4))**2*zb(j3,j4)*
     .zb(j5,j6))+(zb(j1,j4)*
     .(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))*
     .(-(za(j1,j2)*zb(j1,j6))+za(j2,j3)*zb(j3,j6))*
     .(3._dp*(-(za(j1,j3)*zb(j1,j4)*zb(j3,j6))-
     .za(j2,j3)*zb(j2,j4)*zb(j3,j6))-
     .zb(j4,j6)*(t(j1,j2,j3)-t(j1,j2,j4))))/
     .((s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s(j5,
     .j6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*
     .(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4))**2*zb(j3,j4)*zb(j5,j6
     .)
     .))


      Fax=Fax
     .+I3m(s(j1,j2),s(j3,j4),s(j5,j6))*
     .(-((za(j2,j4)*za(j4,j5)*zb(j1,j3)*zb(j3,j6))/
     .(s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s(j5,j
     .6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2))-
     .(3._dp*(-s(j1,j2)+s(j3,j4)-s(j5,j6))*
     .((s(j1,j2)-s(j3,j4)-s(j5,j6))*za(j2,j5)*zb(j1,j2)+
     .(-s(j1,j2)-s(j3,j4)+s(j5,j6))*za(j5,j6)*zb(j1,j6))*
     .(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))*
     .(-(za(j1,j2)*zb(j1,j6))+za(j2,j3)*zb(j3,j6)))/
     .((s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s(j5,
     .j6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)**2*
     .(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4)))-
     .(3._dp*(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))*
     .(-(za(j1,j2)*za(j2,j5)*zb(j1,j2)*zb(j1,j6))-
     .za(j2,j3)*za(j4,j5)*zb(j1,j4)*zb(j3,j6)-
     .za(j2,j4)*za(j3,j5)*zb(j1,j3)*zb(j4,j6)-
     .za(j2,j5)*za(j5,j6)*zb(j1,j6)*zb(j5,j6)))/
     .(2._dp*(s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s
     .(j5,j6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*
     .(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4)))-
     .(za(j4,j5)*zb(j1,j3)*(-(za(j1,j2)*zb(j1,j6))+za(j2,j3)*zb(j3,j6)))
     ./
     .(2._dp*(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4))*t(j1,j2,j3))+
     .(za(j3,j5)*zb(j1,j4)*(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))*
     .(-(za(j1,j2)*zb(j1,j6))+za(j2,j3)*zb(j3,j6))*
     .(t(j1,j2,j3)-t(j1,j2,j4)))/
     .((s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s(j5,
     .j6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*
     .(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4))**2))


      Fax=Fax
     .-(L1(-s(j5,j6),-t(j1,j2,j4))*za(j3,j5)*zb(j1,j3)*
     .(za(j2,j5)*zb(j1,j2)+za(j4,j5)*zb(j1,j4))*zb(j2,j3))/
     .(za(j5,j6)*zb(j1,j2)*zb(j2,j4)*
     .(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4))*t(j1,j2,j4))


      Fax=Fax
     .-(Lsm1_2mh(s(j3,j4),t(j1,j2,j4),s(j1,j2),s(j5,j6))*
     .((za(j2,j3)*zb(j1,j2)-za(j3,j4)*zb(j1,j4))**2*
     .(-(za(j1,j5)*zb(j1,j4))-za(j2,j5)*zb(j2,j4))**2-
     .za(j3,j5)**2*zb(j1,j4)**2*t(j1,j2,j4)**2))/
     .(2._dp*za(j5,j6)*zb(j1,j2)*(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j
     .4))**4)


      Fax=Fax
     .+I3m(s(j1,j2),s(j3,j4),s(j5,j6))*
     .(-((za(j2,j4)*za(j4,j5)*zb(j1,j3)*zb(j3,j6))/
     .(s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s(j5,j
     .6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2))-
     .(3._dp*(-s(j1,j2)+s(j3,j4)-s(j5,j6))*
     .(za(j2,j5)*zb(j1,j2)+za(j4,j5)*zb(j1,j4))*
     .(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))*
     .(-((s(j1,j2)-s(j3,j4)-s(j5,j6))*za(j1,j2)*zb(j1,j6))-
     .(-s(j1,j2)-s(j3,j4)+s(j5,j6))*za(j2,j5)*zb(j5,j6)))/
     .((s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s(j5,
     .j6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)**2*
     .(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4)))-
     .(3._dp*(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))*
     .(-(za(j1,j2)*za(j2,j5)*zb(j1,j2)*zb(j1,j6))-
     .za(j2,j3)*za(j4,j5)*zb(j1,j4)*zb(j3,j6)-
     .za(j2,j4)*za(j3,j5)*zb(j1,j3)*zb(j4,j6)-
     .za(j2,j5)*za(j5,j6)*zb(j1,j6)*zb(j5,j6)))/
     .(2._dp*(s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s
     .(j5,j6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*
     .(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4)))-
     .(za(j2,j4)*(za(j2,j5)*zb(j1,j2)+za(j4,j5)*zb(j1,j4))*zb(j3,j6))/
     .(2._dp*(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4))*t(j1,j2,j4))+
     .(za(j2,j3)*(za(j2,j5)*zb(j1,j2)+za(j4,j5)*zb(j1,j4))*
     .(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))*zb(j4,j6)*
     .(-t(j1,j2,j3)+t(j1,j2,j4)))/
     .((s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s(j5,
     .j6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*
     .(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4))**2))


      Fax=Fax
     .+Lnrat(-s(j1,j2),-s(j3,j4))*(-((za(j2,j4)*za(j4,j5)**2*zb(j1,j3))/
     .((s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s(j5,
     .j6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*za(j3,j4)*za(j5,j6)))-
     .(za(j2,j4)*za(j3,j5)*(za(j2,j5)*zb(j1,j2)+za(j4,j5)*zb(j1,j4)))/
     .(za(j3,j4)*za(j5,j6)*(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4))*
     .*
     .2)+(6._dp*za(j1,j2)*(za(j2,j5)*zb(j1,j2)+za(j4,j5)*zb(j1,j4))*
     .(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))*
     .((-s(j1,j2)+s(j3,j4)-s(j5,j6))*zb(j1,j6)-
     .2._dp*za(j2,j5)*zb(j1,j2)*zb(j5,j6)))/
     .((s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s(j5,
     .j6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)**2*
     .(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4)))+
     .(za(j2,j3)*(za(j2,j5)*zb(j1,j2)+za(j4,j5)*zb(j1,j4))*
     .(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))*
     .(3._dp*(-(za(j1,j3)*za(j4,j5)*zb(j1,j4))-
     .za(j2,j3)*za(j4,j5)*zb(j2,j4))-
     .za(j3,j5)*(-t(j1,j2,j3)+t(j1,j2,j4))))/
     .((s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s(j5,
     .j6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*za(j3,j4)*za(j5,j6)*
     .(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4))**2))


      Fax=Fax
     .+Lnrat(-s(j5,j6),-s(j3,j4))*((za(j4,j5)*zb(j1,j3)**2*zb(j3,j6))/
     .((s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s(j5,
     .j6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*zb(j1,j2)*zb(j3,j4))+
     .(zb(j1,j4)*(za(j3,j5)*zb(j1,j3)-za(j5,j6)*zb(j1,j6))*zb(j3,j6))/
     .(zb(j1,j2)*zb(j3,j4)*(-(za(j3,j5)*zb(j4,j5))-za(j3,j6)*zb(j4,j6))*
     .*
     .2)+(6._dp*(za(j3,j5)*zb(j1,j3)-za(j5,j6)*zb(j1,j6))*
     .(-((-s(j1,j2)+s(j3,j4)-s(j5,j6))*za(j2,j5))+
     .2._dp*za(j1,j2)*za(j5,j6)*zb(j1,j6))*
     .(-(za(j4,j5)*zb(j3,j5))-za(j4,j6)*zb(j3,j6))*zb(j5,j6))/
     .((s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s(j5,
     .j6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)**2*
     .(-(za(j3,j5)*zb(j4,j5))-za(j3,j6)*zb(j4,j6)))+
     .((za(j3,j5)*zb(j1,j3)-za(j5,j6)*zb(j1,j6))*
     .(-(za(j4,j5)*zb(j3,j5))-za(j4,j6)*zb(j3,j6))*zb(j4,j6)*
     .(3._dp*(za(j3,j5)*zb(j1,j3)*zb(j4,j5)+
     .za(j3,j6)*zb(j1,j3)*zb(j4,j6))+
     .zb(j1,j4)*(t(j3,j5,j6)-t(j4,j5,j6))))/
     .((s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s(j5,
     .j6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*zb(j1,j2)*zb(j3,j4)*
     .(-(za(j3,j5)*zb(j4,j5))-za(j3,j6)*zb(j4,j6))**2))


      Fax=Fax
     .+Lnrat(-s(j5,j6),-s(j3,j4))*((za(j2,j4)**2*za(j4,j5)*zb(j3,j6))/
     .((s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s(j5,
     .j6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*za(j1,j2)*za(j3,j4))+
     .(za(j2,j3)*za(j4,j5)*(za(j2,j4)*zb(j4,j6)+za(j2,j5)*zb(j5,j6)))/
     .(za(j1,j2)*za(j3,j4)*(-(za(j3,j5)*zb(j4,j5))-za(j3,j6)*zb(j4,j6))*
     .*2)
     .-(6._dp*za(j5,j6)*(-(za(j4,j5)*zb(j3,j5))-za(j4,j6)*zb(j3,j6))*
     .(za(j2,j4)*zb(j4,j6)+za(j2,j5)*zb(j5,j6))*
     .(-((-s(j1,j2)+s(j3,j4)-s(j5,j6))*zb(j1,j6))+
     .2._dp*za(j2,j5)*zb(j1,j2)*zb(j5,j6)))/
     .((s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s(j5,
     .j6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)**2*
     .(-(za(j3,j5)*zb(j4,j5))-za(j3,j6)*zb(j4,j6)))+
     .(za(j3,j5)*(-(za(j4,j5)*zb(j3,j5))-za(j4,j6)*zb(j3,j6))*
     .(za(j2,j4)*zb(j4,j6)+za(j2,j5)*zb(j5,j6))*
     .(3._dp*(za(j2,j4)*za(j3,j5)*zb(j4,j5)+
     .za(j2,j4)*za(j3,j6)*zb(j4,j6))+
     .za(j2,j3)*(-t(j3,j5,j6)+t(j4,j5,j6))))/
     .((s(j1,j2)**2-2._dp*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2._dp*s(j1,j2)*s(j5,
     .j6)-
     .2._dp*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*za(j1,j2)*za(j3,j4)*
     .(-(za(j3,j5)*zb(j4,j5))-za(j3,j6)*zb(j4,j6))**2))
      elseif(st=='q+qb-g+g+') then
      Fax=
     .(Lnrat(-t(j1,j2,j3),-s(j5,j6))*za(j2,j5)**2)/(za(j1,j2)*za(j3,j4)*
     .*2*za(j5,j6))-
     .(Lnrat(-t(j1,j2,j4),-s(j5,j6))*za(j2,j5)**2)/(za(j1,j2)*za(j3,j4)*
     .*2*za(j5,j6))-
     .(Lsm1_2me(t(j1,j2,j3),t(j1,j2,j4),s(j1,j2),s(j5,j6))*za(j2,j5)*
     .(za(j2,j4)*za(j3,j5)+za(j2,j3)*za(j4,j5)))/
     .(2._dp*za(j1,j2)*za(j3,j4)**3*za(j5,j6))-
     .(L0(-t(j1,j2,j3),-s(j1,j2))*za(j2,j5)*za(j3,j5)*zb(j1,j3))/
     .(s(j1,j2)*za(j3,j4)**2*za(j5,j6))+
     .(L0(-t(j1,j2,j4),-s(j1,j2))*za(j2,j5)*za(j4,j5)*zb(j1,j4))/
     .(s(j1,j2)*za(j3,j4)**2*za(j5,j6))+
     .(((L1(-t(j1,j2,j4),-s(j5,j6))*s(j3,j4))/s(j5,j6)**2+
     .L0(-t(j1,j2,j4),-s(j5,j6))/s(j5,j6))*za(j2,j3)*za(j2,j5)*zb(j3,j6)
     .)/
     .(za(j1,j2)*za(j3,j4)**2)-(L1(-t(j1,j2,j4),-s(j5,j6))*za(j2,j3)*za(
     .j2,j5)*
     .zb(j1,j3)*zb(j3,j6))/(s(j5,j6)**2*za(j2,j4)*za(j3,j4))-
     .(((L1(-t(j1,j2,j3),-s(j5,j6))*s(j3,j4))/s(j5,j6)**2+
     .L0(-t(j1,j2,j3),-s(j5,j6))/s(j5,j6))*za(j2,j4)*za(j2,j5)*zb(j4,j6)
     .)/
     .(za(j1,j2)*za(j3,j4)**2)+(L1(-t(j1,j2,j3),-s(j5,j6))*(s(j1,j4)+s(j
     .3,j4))*
     .za(j2,j5)*zb(j4,j6))/(s(j5,j6)**2*za(j1,j3)*za(j3,j4))-
     .(za(j2,j5)*(((-(za(j2,j3)*zb(j1,j3))-za(j2,j4)*zb(j1,j4))*zb(j3,j6
     .))/
     .za(j2,j4)+(zb(j4,j6)*t(j1,j3,j4))/za(j1,j3)))/
     .(12._dp*mtsq*s(j5,j6)*za(j3,j4))
      endif

      return
      end
