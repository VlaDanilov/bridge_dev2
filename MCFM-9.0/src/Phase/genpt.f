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
 
      subroutine genpt(xr,ptmin,exact,pt,xjac)
      implicit none
      include 'types.f'
c--- given a random number xr between 0 and 1, generate a transverse
c--- momentum pt according to the flag exact:
c---    exact=true:  generate according to dpt/pt down to ptmin
c---    exact=false: generate down to pt=0 with a shape determined by ptmin
c---
c--- returns: pt and xjac, the Jacobian of the transformation from pt dpt to dxr
      logical:: exact
      real(dp):: xr,ptmin,pt,xjac,hmin,hmax,h,delh,ptmax
      include 'energy.f'

      ptmax=sqrts/2._dp

      if (exact) then
        hmin=1._dp/ptmax
        hmax=1._dp/ptmin
        delh=hmax-hmin
        h=hmin+xr*delh
        pt=1._dp/h
        xjac=delh/h**3
      else
        pt=2._dp*ptmin*ptmax*xr/(2._dp*ptmin+ptmax*(1._dp-xr))
        xjac=pt*ptmax/2._dp/ptmin/(2._dp*ptmin+ptmax)*(2._dp*ptmin+pt)**2
      endif

      return
      end
      
