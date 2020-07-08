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
 
*
* ADJUSTED TO REMOVE ALL 1/eps**2 PIECES IF EPINV2 IN EPINV2.F = 0._dp
* ADJUSTED TO REMOVE ALL 1/eps**2 PIECES IF EPINV2 IN EPINV2.F = 0._dp
* ADJUSTED TO REMOVE ALL 1/eps**2 PIECES IF EPINV2 IN EPINV2.F = 0._dp
*
      function vvg(st,j1,j2,j3,j4,j5,j6)
      implicit none
      include 'types.f'
      complex(dp):: vvg
************************************************************************
*     Author: R.K. Ellis                                               *
*     July, 1998.                                                      *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'epinv.f'
      include 'scale.f'
      include 'epinv2.f'
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: Lnrat,xl12,xl34,xl23,xl56,Vcc,Vsc
      character*9 st


      xl12=Lnrat(musq,-s(j1,j2))
      xl34=Lnrat(musq,-s(j3,j4))
      xl23=Lnrat(musq,-s(j2,j3))
      xl56=Lnrat(musq,-s(j5,j6))

      if(st=='q+g-g+qb-') then
      Vcc=-four
     & -(epinv*epinv2+xl12*epinv+half*xl12**2)
     & -(epinv*epinv2+xl23*epinv+half*xl23**2)
     & -(epinv*epinv2+xl34*epinv+half*xl34**2)-two*(epinv+xl56)
      Vsc=half*(one+epinv+xl56)

      elseif(st=='q+g+g-qb-') then
      Vcc=-four
     &   -(epinv*epinv2+xl12*epinv+half*xl12**2)
     &   -(epinv*epinv2+xl23*epinv+half*xl23**2)
     &   -(epinv*epinv2+xl34*epinv+half*xl34**2)-two*(epinv+xl56)
      Vsc=half*(one+epinv+xl56)

      elseif(st=='q+g+g+qb-') then
      Vcc=-four
     &  -(epinv*epinv2+xl12*epinv+half*xl12**2)
     &  -(epinv*epinv2+xl23*epinv+half*xl23**2)
     &  -(epinv*epinv2+xl34*epinv+half*xl34**2)-two*(epinv+xl56)
      Vsc=half*(one+epinv+xl56)

      elseif(st=='q+g+qb-g-') then
      Vcc=-four
     & -(+(epinv*epinv2+xl12*epinv+half*xl12**2)
     &   +(epinv*epinv2+xl23*epinv+half*xl23**2))-two*(epinv+xl56)
      Vsc=half*(one+epinv+xl56)
 

      elseif(st=='q+g+qb-g+') then
      Vcc=-four
     & -(epinv*epinv2+xl12*epinv+half*xl12**2)
     & -(epinv*epinv2+xl23*epinv+half*xl23**2)-two*(epinv+xl56)
      Vsc=half*(one+epinv+xl56)

      elseif(st=='q+qb-g-g+') then 
      Vcc=-four-(epinv*epinv2+xl12*epinv+half*xl12**2)-two*(epinv+xl56)
      Vsc=half*(one+epinv+xl56)

      elseif(st=='q+qb-g+g-') then
      Vcc=-four-(epinv*epinv2+xl12*epinv+half*xl12**2)-two*(epinv+xl56)
      Vsc=half*(one+epinv+xl56)

      elseif(st=='q+qb-g+g+') then
      Vcc=-four-(epinv*epinv2+xl12*epinv+half*xl12**2)-two*(epinv+xl56)
      Vsc=half*(one+epinv+xl56)
      else 
      write(6,*) 'unimplemented st',st
      stop
      endif

      vvg=Vcc+Vsc
      end
