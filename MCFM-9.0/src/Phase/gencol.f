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
 
      subroutine gencol(x,xjac,xmin,emit,r)
      implicit none
      include 'types.f'
      
c---Generate an x value and store it for later retrieval
      integer:: emit,lemit
      real(dp):: x,xjac,xmin,xl,xljac,xlmin,r
      save xl,xljac,xlmin,lemit
      x=1._dp-(1._dp-xmin)*abs(1._dp-2._dp*r)
      xjac=2._dp*(1._dp-xmin)
      xl=x
      xljac=xjac
      xlmin=xmin
      lemit=emit
      return

      entry getcol(x,xjac,xmin,emit)
c---return the same values as last time
      x=xl
      xjac=xljac
      xmin=xlmin
      emit=lemit
      end
