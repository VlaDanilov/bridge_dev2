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
 
module SCET
    use types
    implicit none
    private

    ! whether to actually do multitau
    logical, public, save :: doMultitaucut = .false.

    real(dp), public, save, allocatable :: tcutarray(:)
    real(dp), public, save :: smallestTaucut

    real(dp), public, save, allocatable :: scetreweight(:)
!$omp threadprivate(scetreweight)

    ! this flag specifies for the above cut pieces the result of
    ! maketaucut for taucut; for each dipole contribution
    logical, public, save :: includeTaucutgrid(0:40)
!$omp threadprivate(includeTaucutgrid)

    ! used to transfor current dipole contribution information
    ! from nplotter to bookplot without changing the bookplot
    ! interface
    integer, public, save :: currentNd
!$omp threadprivate(currentNd)


end module
