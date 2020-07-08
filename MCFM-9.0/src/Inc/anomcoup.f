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
 
      real(dp) :: delg1_z,delg1_g,lambda_g,lambda_z
      real(dp) ::  h1Z,h2Z,h3Z,h4Z,h1gam,h2gam,h3gam,h4gam
      real(dp) :: delk_g,delk_z,tevscale
      real(dp) :: h1tZ,h2tZ,h3tZ,h4tZ,h1tgam,h2tgam,h3tgam,h4tgam
      logical :: anomtgc
      real(dp) :: hitZ(4), hitgam(4)
      common/anomcoup1/delg1_z,delg1_g,lambda_g,lambda_z,delk_g,delk_z
      common/anomcoup2/h1Z,h2Z,h3Z,h4Z,h1gam,h2gam,h3gam,h4gam,tevscale,hitZ,hitgam,anomtgc
!$omp threadprivate(/anomcoup1/,/anomcoup2/)

