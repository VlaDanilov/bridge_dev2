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
 
c--- A subroutine for locating x in xx (an ordered array where 
c--- xx(n) > xx(1) has a meaning) 

      subroutine locatemcfm(xx,n,x,j) 
      implicit none
      include 'types.f'
      
      integer:: n,inp,np,j,nav 
      real(dp):: xx(n),x
       

      inp=0
      np=n+1
 10   if((np-inp) > 1) then 
         nav=(np+inp)/2
         if((xx(n)>xx(1)).eqv.(x>xx(nav))) then 
            inp=nav
         else
            np=nav
         endif
         GOTO 10
         endif
         j=inp
         return 
         end
      
