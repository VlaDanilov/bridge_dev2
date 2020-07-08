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
 
      function etmiss(p,etvec)
      implicit none
      include 'types.f'
      real(dp):: etmiss
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: j,k
      logical:: is_neutrino,is_darkmatter
      real(dp):: etvec(4),p(mxpart,4)
      
      do k=1,4
        etvec(k)=0._dp
      enddo
      
      do j=1,mxpart
        if (is_neutrino(j) .or. is_darkmatter(j)) then
          do k=1,4
            etvec(k)=etvec(k)+p(j,k)
          enddo
        endif
      enddo
      
      etmiss=sqrt(etvec(1)**2+etvec(2)**2)
      
      return
      end
      
