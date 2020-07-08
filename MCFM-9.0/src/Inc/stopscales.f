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
 
      real(dp):: initfacscale_H,initfacscale_L,initrenscale_H,initrenscale_L
      common/stopscales1/initfacscale_H,initfacscale_L,initrenscale_H,initrenscale_L
!$omp threadprivate(/stopscales1/)
      real(dp):: facscale_H,facscale_L,renscale_H,renscale_L
      real(dp):: msqLH(-nf:nf,-nf:nf),msqHL(-nf:nf,-nf:nf)
      real(dp):: as_H,as_L
      common/stopscales2/facscale_H,facscale_L,renscale_H,renscale_L,as_H,as_L,msqLH,msqHL
!$omp threadprivate(/stopscales2/)
