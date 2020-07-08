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
 
module CPUTime
#ifdef HAVE_MPI
    use mpi
#endif
    use types

    public :: start_times, get_cputime, get_walltime

    private

    real(dp), save :: cputime_start
    integer, save :: walltime_start

    contains

    subroutine start_times
        implicit none

        call cpu_time(cputime_start)
        call system_clock(walltime_start)
    end subroutine

    ! returns cpu time in seconds
    function get_cputime()
        implicit none
        real(dp) :: get_cputime
        real(dp) :: cputime_end
#ifdef HAVE_MPI
        integer :: ierr
        include 'mpicommon.f'
#endif

        call cpu_time(cputime_end)

        get_cputime = cputime_end - cputime_start
#ifdef HAVE_MPI
        if (world_size > 1) then
            call mpi_allreduce(mpi_in_place, get_cputime, 1, mpi_double_precision, &
                mpi_sum, mpi_comm_world, ierr)
        endif
#endif
    end function

    ! returns walltime in seconds
    function get_walltime()
        implicit none
        real(dp) :: get_walltime
        integer :: walltime_end, count_rate

        call system_clock(walltime_end, count_rate)

        get_walltime = real(walltime_end - walltime_start,dp)/real(count_rate,dp)

    end function


end module
