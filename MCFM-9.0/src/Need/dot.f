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
 
      function dot(p,i,j)
      implicit none
      include 'types.f'
      real(dp):: dot
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: i,j
      real(dp):: p(mxpart,4)
      dot=p(i,4)*p(j,4)-p(i,1)*p(j,1)-p(i,2)*p(j,2)-p(i,3)*p(j,3)
      return
      end

      function dotvec(pi,pj)
      implicit none
      include 'types.f'
      real(dp):: dotvec, pi(4), pj(4)
      
      dotvec=pi(4)*pj(4)-pi(1)*pj(1)-pi(2)*pj(2)-pi(3)*pj(3)
      return
      end

      function massvec(p)
        implicit none
        include 'types.f'

        real(dp) :: massvec
        real(dp), intent(in) :: p(4)

        massvec=p(4)*p(4)-p(1)*p(1)-p(2)*p(2)-p(3)*p(3)
      end

