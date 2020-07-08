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
 
      subroutine scaleset_shat(p,mu0)
      implicit none
      include 'types.f'
c--- subroutine to calculate dynamic scale equal to
c--- invariant mass of particles 3 and 4
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'kpart.f'
      real(dp):: p(mxpart,4),mu0

      if (kpart==klord) then
        mu0=(p(1,4)+p(2,4))**2-(p(1,1)+p(2,1))**2
     &     -(p(1,2)+p(2,2))**2-(p(1,3)+p(2,3))**2       
        mu0=sqrt(abs(mu0))
      else
        write(6,*) 'dynamicscale s-hat not supported beyond LO.'
        stop
      endif
      
      return
      end
      
