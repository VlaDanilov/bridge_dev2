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
 
      function a6treega(st,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: a6treega
      integer:: j1,j2,j3,j4,j5,j6
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: st
      real(dp):: t
c---- helicity stamps
c---- 'q+qb-g+ga+'=1
c---- 'q+qb-g+ga-'=2
c---- 'q+qb-g-ga+'=3
      if(st==3) then
      a6treega=
     .(za(j2,j3)**2*zb(j6,j1)**2)/(t(j2,j3,j4)*za(j2,j4)*
     .(za(j4,j5)*zb(j5,j1)+za(j4,j6)*zb(j6,j1))*zb(j6,j5)) +
     .(za(j2,j3)*zb(j2,j1)*(za(j5,j2)*zb(j2,j1)+za(j5,j3)*zb(j3,j1))**2)/
     .(s(j2,j3)*za(j5,j6)*(za(j4,j2)*zb(j2,j1)+za(j4,j3)*zb(j3,j1))*
     .(za(j4,j1)*zb(j1,j2)+za(j4,j3)*zb(j3,j2))*zb(j3,j1)) +
     .(zb(j4,j1)*(za(j3,j1)*zb(j1,j6)+za(j3,j4)*zb(j4,j6))**2)/
     .(t(j1,j3,j4)*s(j1,j4)*(za(j4,j1)*zb(j1,j2)+za(j4,j3)*zb(j3,j2))*zb(j6,j5))
      elseif(st==2) then
      a6treega=
     .(za(j2,j4)**2*zb(j6,j1)**2)/(t(j2,j3,j4)*za(j2,j3)*
     .(za(j3,j5)*zb(j5,j1)+za(j3,j6)*zb(j6,j1))*zb(j6,j5)) +
     .(za(j2,j4)*zb(j2,j1)*(za(j5,j2)*zb(j2,j1)+za(j5,j4)*zb(j4,j1))**2)/
     .(s(j2,j4)*za(j5,j6)*(za(j3,j2)*zb(j2,j1)+za(j3,j4)*zb(j4,j1))*
     .(za(j3,j1)*zb(j1,j2)+za(j3,j4)*zb(j4,j2))*zb(j4,j1)) +
     .(zb(j3,j1)*(za(j4,j1)*zb(j1,j6)+za(j4,j3)*zb(j3,j6))**2)/
     .(t(j1,j3,j4)*s(j1,j3)*(za(j3,j1)*zb(j1,j2)+za(j3,j4)*zb(j4,j2))*zb(j6,j5))
      elseif(st==1) then
      a6treega= 
     .(-za(j1,j2)*za(j2,j5)**2)/(za(j1,j3)*za(j1,j4)*za(j2,j3)*za(j2,j4)*za(j5,j6))
      else 
      write(6,*) 'unimplemented st'
      stop
      endif

      return
      end
