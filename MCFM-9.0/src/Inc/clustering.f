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
 
      logical:: clustering,inclusive
      integer:: jetalgorithm
      character*4 algorithm
      common/clustering/clustering,inclusive
      common/algorithm/algorithm
      common/jetalgorithm/jetalgorithm
      integer, parameter ::
     & kt=1, antikt=2, Rsepcone=3, hqrk=4, noclustering=5
