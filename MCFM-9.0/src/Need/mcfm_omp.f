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
 
      program mcfm
          use omp_lib
          use types
#ifdef HAVE_MPI
          use mpi
#endif
          use CPUTime
          use MCFMBenchmark
          use m_config, only : cfg_var_configadded
          use parseinput, only : cfg
          implicit none
          real(dp) :: r0,er0
          integer :: threadmax
          integer :: ierr,mylen,support
#ifdef HAVE_MPI
          character (len=mpi_max_processor_name) :: procname
#endif
          include 'mpicommon.f'

          integer:: exitCode

          threadmax = omp_get_max_threads()
          call omp_set_num_threads(threadmax)

#ifdef HAVE_MPI
          call mpi_init_thread(mpi_thread_multiple,support,ierr)
          call mpi_comm_rank(mpi_comm_world,rank,ierr)
          call mpi_comm_size(mpi_comm_world,world_size,ierr)
          call mpi_get_processor_name(procname,mylen,ierr)
#else
          world_size = 1
          rank = 0
#endif

          call start_times()

          call mcfmmain(r0,er0)

          if (cfg_var_configadded(cfg, "extra%benchmark")) then
              exitCode = comparisonCode(2)
          else
              exitCode = 0
          endif

#ifdef HAVE_MPI
          ! this is to have a consistent exit code for all mpi ranks
          call mpi_bcast(exitCode, 1, mpi_integer, 0, mpi_comm_world, ierr)
          call mpi_finalize(ierr)
#endif

          if (exitCode /= 0) then
              error stop 1
          else
              stop 0
          endif

      end program
