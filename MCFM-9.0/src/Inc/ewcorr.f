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
 
      character*24 ewcorr
      integer kewcorr
      integer, parameter:: knone=1, ksudakov=2, kexact=3
      real(dp):: wt_noew
      common/kewcorr/kewcorr
      common/ewcorr/ewcorr
      common/wt_noew/wt_noew
!$omp threadprivate(/wt_noew/)
      
