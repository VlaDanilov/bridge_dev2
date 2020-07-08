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
 
      subroutine mcfm_init()
          use omp_lib
          use mod_qcdloop_c
          use parseinput
          use MCFMStorage
          use MCFMBenchmark, only : setupBenchmark
          use iso_fortran_env
      implicit none
************************************************************************
*                                                                      *
*  This routine should initialize any necessary variables and          *
*  perform the usual print-outs                                        *
*                                                                      *
************************************************************************
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'cutoff.f'
      include 'limits.f'
      include 'npart.f'
      include 'phasemin.f'
      include 'facscale.f'
      include 'scale.f'
      include 'verbose.f'
      include 'phot_dip.f'
      include 'includect.f'
      include 'TRtensorcontrol.f'
      include 'frag.f'
      include 'energy.f'
      include 'TRmaxindex.f'
      include 'masses.f'
      include 'nflav.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'mpicommon.f'
      include 'nproc.f'
      include 'runname.f'

      include 'rtsmin.f'
     
      include 'ptilde.f'
      include 'APPLinclude.f'
      real(dp):: p1ext(4),p2ext(4),p(mxpart,4),val
      integer:: j,k
      common/pext/p1ext,p2ext
      data p/mxpart*3._dp,mxpart*4._dp,mxpart*0._dp,mxpart*5._dp/

      integer :: nplotmax
      common/nplotmax/nplotmax

      character(len=:), allocatable :: cmd
      character(len=cfg_string_len) :: benchmark

!$omp threadprivate(/pext/)

      if (rank.eq.0) call banner

      if (rank == 0 .and. world_size > 1) then
          write (*,*) ""
          write (*,*) "Running MCFM with ", world_size, " MPI processes"
          write (*,*) ""
      endif
      if (rank == 0) then
          write (*,*) ""
          write (*,*) "Running MCFM with ", omp_get_max_threads(), " OMP threads" 
          write (*,*) ""
      endif

      if (rank == 0) then
          write (*,*) ""
          write (*,*) "MCFM compiled with ", compiler_version(), " using the options ",
     &          compiler_options()
          write (*,*) ""
      endif

      if (rank ==0) then
          allocate(character(len=1024) :: cmd)
          call get_command(cmd)
          write (*,*) ""
          write (*,*) "Running MCFM as ", trim(cmd)
          write (*,*) ""
          deallocate(cmd)
      endif


c initialize qcdloop cache for all omp threads
!$omp parallel
      call qlcachesize(1)
!$omp end parallel

*Initialise data in commonblock /pvmaxindex/
      maxcindex=3
      maxdindex=4
      maxeindex=5

      ! TODO: add this to configuration file
* Initialize parameter settings
      call mdata

      epinv=1.e1_dp
      epinv2=1.e1_dp

      if (rank.ge.1) verbose=.false.

      call default_config()

      call cfg_update_from_arguments(cfg)
      call cfg_sort(cfg)

      if (cfg_var_configadded(cfg, "extra%benchmark")) then
          call cfg_get(cfg, "extra%benchmark", benchmark)
          write(6,*)
          write(6,*) '****************************************'
          write(6,*) '*      Running in benchmark mode       *'
          write(6,*) '* bench = ', benchmark
          write(6,*) '****************************************'
          write(6,*)
          call setupBenchmark()
      endif

      call read_config()
      !call cfg_write(cfg, "stdout")

      if (verbose) then
      write(6,*)
      write(6,*) '****************************************'
      write(6,*) '*     Cross section in femtobarns      *'
      write(6,*) '****************************************'
      write(6,*)
      endif

* Counter-terms for radiation in top decay should be included
      includect=.true.
      
* Set-up incoming beams and PS integration cut-offs
c--- Note: version 6.4 onwards, scale cutoff with c.o.m. energy
c--- Note: since Sep. 2018 no longer scale since it is dimensionless
!      cutoff=cutoff*(sqrts/2000._dp)**2
      rtsmin=min(rtsmin,sqrt(wsqmin+cutoff))
      rtsmin=min(rtsmin,sqrt(bbsqmin+cutoff))
      taumin=(rtsmin/sqrts)**2
      xmin=1.e-8_dp

      p1ext(4)=-half*sqrts
      p1ext(1)=0._dp
      p1ext(2)=0._dp
      p1ext(3)=-half*sqrts

      p2ext(4)=-half*sqrts
      p2ext(1)=0._dp
      p2ext(2)=0._dp
      p2ext(3)=+half*sqrts

* Set-up run name
      call setrunname(scale,facscale,cfg%config_directory)

* npart=9 is a dummy value, to ensure that all histograms are included
      npart=9
      val=1.e-15_dp   
      nprocbelow = nproc

      ! This first call has two purposes:
      ! 1. determine nplotmax for histogram allocation.
      ! 2. must be run in parallel so that each thread initializes
      !    in "tagbook" mode and takes runs fully in "tagplot"
      !    for subsequent calls in the integration routines
!$omp parallel
      call nplotter(p,val,val**2,1)
!$omp end parallel

!$omp parallel
      call initHistogramStorage(nplotmax)
!$omp end parallel
      call initMasterStorage(nplotmax)

      p(:,:) = 0._dp
       
* Initialize flag for photon fragmentation dipoles
      phot_dip(:)=.false.
      fragint_mode=.false. 
* Initialize integer:: used in TensorReduction to zero
      TRtensorcontrol=0
      
      end
            
