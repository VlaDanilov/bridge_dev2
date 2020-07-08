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
 
      subroutine mcfm_vegas_adaptive(integ,integ_err)
          use types
          use Integration
          use mod_sobseq
          use MCFMStorage
          use MCFMPrint
          use PDFerrors
          use Scalevar
          use SCET, only : doMultitaucut, scetreweight, includeTaucutgrid
          use parseinput
#ifdef HAVE_MPI
          use mpi
#endif
! This is to avoid compiler version dependence of mpi.mod
!#ifdef HAVE_MPI
!          implicit none
!          include 'mpif.h'
!#endif
!          implicit none
          include 'ipsgen.f'
          include 'reset.f'
          include 'kpart.f'
          include 'kprocess.f'
          include 'vegas_common.f'
          include 'nproc.f'
          include 'taucut.f'
          include 'mpicommon.f'
          include 'frag.f'

          include 'vegas_damp.f'

          real(dp), intent(out) :: integ, integ_err

          real(dp) :: region(2*mxdim)
          real(dp) :: lowint, virtint, realint, scetint, fragint
          external lowint, virtint, realint, scetint, fragint
          character*1 getstr

          ! these hold values for all maxParts integration parts and a maximum of
          ! maxIps ipsgen contributions
          real(dp) :: totalsig, totalsd
          ! this array holds modified sd values to give
          ! easier integrations slight preference, see below how it's used
          real(dp) :: sdMod(maxParts,maxIps)
          integer :: initcalls(maxParts)
          integer :: maxipsgens(maxParts)
          logical :: computePart(maxParts)
          logical :: doFirstCall(maxParts,maxIps)
          logical :: warmupComplete(maxParts,maxIps)

          type(sobol_state) :: sstates(maxParts,maxIps,25)
          integer(kind=8) :: strides
          real(dp) :: skiprnd

          ! this sets the vegas grid dampening parameter for all parts
          ! in warmup(1) and final(2) runs
          ! setting it to 0 will disable grid adjustments
          real(dp) :: vegasDampPart(maxParts,2)

          ! saves whether initcalls has been specified through configuration
          logical :: callsConfigAdded(maxParts)

          real(dp) :: xreal,xreal2
          common/xreal/xreal,xreal2

          integer :: stage
          integer :: ierr, ierr2
          integer :: hist_readsnapshot, hist_writesnapshot
          integer :: nprocabove
          logical :: origCoeffonly

          ! number of iterations to compute in one batch
          integer :: iterBatchWarmup
          integer :: iterBatch, iterBatch1, iterBatch2
          ! multiplicator for number of calls if precision per batch is
          ! too low in warmup; also multiplier for full run
          real(dp) :: iterCallMult
          ! precision goal in percent of one batch for warmup
          real(dp) :: warmupPrecisionGoal
          real(dp) :: warmupChisqGoal
          ! precision goal of final result
          real(dp) :: resultPrecisionGoal

          ! factor to multiply number of calls after warmup for full run
          real(dp) :: callBoost

          logical :: usesobol
          logical :: writeIntermediate

          ! when this number of calls per iteration has been reached, 
          ! switch iterBatch to 1
          integer :: maxCallsPerIter

          integer :: nextLoc(2)
          integer :: i,j,k
          integer :: part,ips
          logical :: warmdone
          
          integer :: nprocextra

          logical :: bin
          common/bin/bin

          logical :: dryrun
          common/dryrun/dryrun

          logical :: readin

          character*255 runname
          common/runname/runname

          integer(kind=8), parameter, dimension(1:30)   :: s = [
     &      1,2,3,3,4,
     &      4,5,5,5,5,
     &      5,5,6,6,6,
     &      6,6,6,7,7,
     &      7,7,7,7,7,
     &      7,7,7,7,7
     &      ]
          integer(kind=8), parameter, dimension(1:30)   :: a = [
     &      0,1,1,2,1,
     &      4,2,4,7,11,
     &      13,14,1,13,16,
     &      19,22,25,1,4,
     &      7,8,14,19,21,
     &      28,31,32,37,41
     &      ]
          integer(kind=8), parameter, dimension(7,1:30) :: m = reshape([
     &                      1,0,0,0,0,0,0,
     &    					1,3,0,0,0,0,0,
     &    					1,3,1,0,0,0,0,
     &    					1,1,1,0,0,0,0,
     &    					1,1,3,3,0,0,0,   ! 5
     &    					1,3,5,13,0,0,0,
     &    					1,1,5,5,17,0,0,
     &    					1,1,5,5,5,0,0,
     &    					1,1,7,11,19,0,0,
     &    					1,1,5,1,1,0,0,   ! 10
     &    					1,1,1,3,11,0,0,
     &    					1,3,5,5,31,0,0,
     &                      1,3,3,9,7,49,0,
     &                      1,1,1,15,21,21,0,
     &                      1,3,1,13,27,49,0, ! 15
     &                      1,1,1,15,7,5,0,
     &                      1,3,1,15,13,25,0,
     &                      1,1,5,5,19,61,0,
     &                      1,3,7,11,23,15,103,
     &                      1,3,7,13,13,15,69, ! 20
     &                      1,1,3,13,7,35,63,
     &                      1,3,5,9,1,25,53,
     &                      1,3,1,13,9,35,107,
     &                      1,3,1,5,27,61,31,
     &                      1,1,5,11,19,41,61,   ! 25
     &                      1,3,5,3,3,13,69,
     &                      1,1,7,13,1,19,1,
     &                      1,3,7,5,13,19,59,
     &                      1,1,3,9,25,29,41,
     &                      1,3,5,13,23,1,55  ! 30
     &      ], [7,30])
          integer(kind=8) :: toskip

          interface
              function fxn(vector,wgt)
                  use types
                  implicit none
                  include 'mxdim.f'
                  real(dp) :: fxn
                  real(dp), intent(in) :: vector(mxdim)
                  real(dp), intent(in) :: wgt
              end function
          end interface

          procedure (fxn), pointer :: ifun => null ()

          call cfg_get(cfg, "integration%precisiongoal", resultPrecisionGoal)
          call cfg_get(cfg, "integration%warmupprecisiongoal", warmupPrecisionGoal)
          call cfg_get(cfg, "integration%warmupchisqgoal", warmupChisqGoal)
          call cfg_get(cfg, "integration%usesobol", usesobol)
          call cfg_get(cfg, "integration%writeintermediate", writeIntermediate)

          call cfg_get_add(cfg, "integration%ndmx", ndmx, 100, "Number of Vegas grid subdivisions.")
          call cfg_get_add(cfg, "integration%callboost", callBoost, 4._dp, "")
          call cfg_get_add(cfg, "integration%itercallmult", iterCallMult, 1.4_dp,
     &      "Multiply calls/it by this factor after every iterBatch number of iterations.")

          call cfg_get_add(cfg, "integration%iterbatch1", iterBatch1, 5,
     &      "Batch of iterations using same calls/it just after warmup.")

          call cfg_get_add(cfg, "integration%iterbatch2", iterBatch2, 1,
     &      "Batch of iterations using same calls/it after first iterbatch1.")

          call cfg_get_add(cfg, "integration%iterbatchwarmup", iterBatchWarmup, 5,
     &      "Number of iterations for warmup using same calls/it.")
          call cfg_get_add(cfg, "integration%maxcallsperiter", maxCallsPerIter, 500000000,
     &      "Switch to iterBatch=1 after this number of calls/it.")

          call cfg_get(cfg, "integration%readin", readin)

          ! no binning for warumup stage
          bin = .false.

          ! we don't need or want this here
          dryrun = .false.
          nprn = 0

          do part=1,maxParts; do ips=1,maxIps
              iterationStorage(part,ips)%used = .false.
              associate (info => iterationStorage(part,ips)%vinfo)
                  info%si = 0._dp
                  info%swgt = 0._dp
                  info%schi = 0._dp
                  info%lastIter = 0
                  info%doFirstCall = .false. ! ensure a first full run?
                  info%warmupComplete = .false.
                  info%damp = [1.5_dp, 0.8_dp]
                  info%useSobol = usesobol
                  info%sobolInitialized = .false.
                  info%sobolSkip = 1000
              end associate
          enddo; enddo

c---      setup inital number of calls for warmup phase
c---      these could also be tuned to specific processes
c---
c---      note that the MC integration error is only accurate if the 
c---      number of calls is "sufficiently" large, especially for 
c---      difficult integrands

          callsConfigAdded(:) = .false.

          call cfg_get_add(cfg, "integration%initcallslord", initcalls(lord),
     &                              100000, "initial calls LO")
          if (cfg_var_configadded(cfg, "integration%initcallslord")) then
              callsConfigAdded(lord) = .true.
          endif

          call cfg_get_add(cfg, "integration%initcallsnlovirt", initcalls(nloVirt),
     &                              100000, "initial calls NLO virt")
          if (cfg_var_configadded(cfg, "integration%initcallsnlovirt")) then
              callsConfigAdded(nloVirt) = .true.
          endif

          call cfg_get_add(cfg, "integration%initcallsnloreal", initcalls(nloReal),
     &                              500000, "initial calls NLO real")
          if (cfg_var_configadded(cfg,"integration%initcallsnloreal")) then
              callsConfigAdded(nloReal) = .true.
          endif

          call cfg_get_add(cfg, "integration%initcallsnlofrag", initcalls(nloFrag),
     &                              100000, "initial calls NLO frag")
          if (cfg_var_configadded(cfg,"integration%initcallsnlofrag")) then
              callsConfigAdded(nloFrag) = .true.
          endif

          call cfg_get_add(cfg, "integration%initcallsnlorealextra", initcalls(nloRealExtra),
     &                              500000, "initial calls NLO real extra")
          if (cfg_var_configadded(cfg,"integration%initcallsnlorealextra")) then
              callsConfigAdded(nloRealExtra) = .true.
          endif

          call cfg_get_add(cfg, "integration%initcallssnlobelow", initcalls(snloBelow),
     &                              200000, "initial calls sNLO below")
          if (cfg_var_configadded(cfg, "integration%initcallssnlobelow")) then
              callsConfigAdded(snloBelow) = .true.
          endif


          call cfg_get_add(cfg, "integration%initcallssnloabove", initcalls(snloAbove),
     &                              400000, "initial calls sNLO above")
          if (cfg_var_configadded(cfg, "integration%initcallssnloabove")) then
              callsConfigAdded(snloAbove) = .true.
          endif

          call cfg_get_add(cfg, "integration%initcallsnnlobelow", initcalls(nnloBelow),
     &                              200000, "initial calls NNLO below")
          if (cfg_var_configadded(cfg, "integration%initcallsnnlobelow")) then
              callsConfigAdded(nnloBelow) = .true.
          endif


          call cfg_get_add(cfg, "integration%initcallsnnlovirtabove", initcalls(nnloVirtAbove),
     &                              1000000, "initial calls nnlo virt above")
          if (cfg_var_configadded(cfg, "integration%initcallsnnlovirtabove")) then
              callsConfigAdded(nnloVirtAbove) = .true.
          endif

          call cfg_get_add(cfg, "integration%initcallsnnlorealabove", initcalls(nnloRealAbove),
     &                              10000000, "initial calls nnlo real above")
          if (cfg_var_configadded(cfg, "integration%initcallsnnlorealabove")) then
              callsConfigAdded(nnloRealAbove) = .true.
          endif

          ! The common block variable maxipsgen is set
          ! to the corresponding maxipsgens value for each part below.
          ! The phase space generation routines can adapt on that value.
          ipsgen = 1
          maxipsgens(:) = 1

          origKpart = kpart
          origCoeffonly = coeffonly

c--- ===================================== 
c--- process specific initializations here
c--- ===================================== 

#define SETUP(part,num) if (.not. callsConfigAdded(part)) initcalls(part) = num

          if (nproc == 1 .or. nproc == 6) then
              SETUP(lord,5d4)
              SETUP(nloReal,1d6)
              SETUP(nloVirt,2d5)
              SETUP(nnloBelow,2d5)
              SETUP(nnloVirtAbove,6d5)
              SETUP(nnloRealAbove,8d6)
          endif

          if (nproc == 300 .or. nproc == 3000) then
            maxipsgens(:) = 2
            SETUP(nloReal,1d6)
            SETUP(nloVirt,2d5)
          endif

          if (nproc == 302 .or. nproc == 3002) then
            maxipsgens(nloVirt) = 2
            maxipsgens(nloReal) = 2
            maxipsgens(lord) = 2
          endif

          if (nproc == 304 .or. nproc == 3004) then
            maxipsgens(lord) = 2
          endif

          if (nproc == 164 .or. nproc==169) then
            iterationStorage(lord,1)%vinfo%damp = [1.5, 0.8]
            SETUP(lord,250000)

            iterationStorage(nloVirt,1)%vinfo%damp = [1.5, 0.8]
            SETUP(nloVirt,100000)

            maxipsgens(nloReal) = 2
            iterationStorage(nloReal,1)%vinfo%damp = [1.5, 0.8]
            SETUP(nloReal,500000)
          endif

          if (nproc == 1281 .or. nproc == 1291 .or. nproc == 1301 .or. 
     &        nproc == 1311 .or. nproc == 1321 .or. nproc == 1282 .or.
     &        nproc == 1292 .or. nproc == 1302 .or. nproc == 1302 .or.
     &        nproc == 1312 .or. nproc == 1322) then
              maxipsgens(:) = 2
          endif

          if (nproc == 226 .or. nproc == 2261) then
              maxipsgens(:) = 2
          endif

c---
c--- setup calls/it
c---
          do j=1,maxParts
              iterationStorage(j,:)%vinfo%callsPerIt = initcalls(j)
          enddo

c---
c--- determine parts to calculate
c--- 
          computePart(:) = .false.

          if (origKpart == klord) then
              computePart(lord) = .true.
          endif

          if (origKpart == kvirt) then
!              origCoeffonly = .true.
              computePart(nloVirt) = .true.
          endif

          if (origKpart == kreal) then
!              origCoeffonly = .true.
              computePart(nloReal) = .true.
          endif

          if (origKpart == kfrag) then
              origCoeffonly = .true.
              computePart(nloFrag) = .true.
          endif

          if (origKpart == ktota) then
              computePart(nloVirt) = .true.
              computePart(nloReal) = .true.

              if (frag) then
                  computePart(nloFrag) = .true.
              endif
          endif

          if (origKpart == ktodk) then
              computePart(nloVirt) = .true.
              computePart(nloReal) = .true.
              computePart(nloRealExtra) = .true.
              call setuprealextra(nprocextra)
          endif


          if (origKpart ==  ksnlo) then
              computePart(snloBelow) = .true.
              computePart(snloAbove) = .true.
              if (onlypowcorr) computePart(snloAbove) = .false.
          endif

          if (origKpart == knnlo) then
            if (origCoeffonly .eqv. .false.) then
              ! alternative enable snloBelow and snloAbove
              computePart(nloVirt) = .true.
              computePart(nloReal) = .true.
            endif

            computePart(nnloBelow) = .true.
            computePart(nnloVirtAbove) = .true.
            computePart(nnloRealAbove) = .true.
            if (onlypowcorr) then
              computePart(nnloVirtAbove) = .false.
              computePart(nnloRealAbove) = .false.
            endif
          endif

          if (knnlopart == knnloRR) then
            computePart(:) = .false.
            computePart(nnloRealAbove) = .true.
          endif

          if (knnlopart == knnloRV) then
            computePart(:) = .false.
            computePart(nnloVirtAbove) = .true.
          endif

          if (knnlopart == knnloVV) then
            computePart(:) = .false.
            computePart(nnloBelow) = .true.
          endif

          if (origKpart == ksnlo .or. origKpart == knnlo) then
            call setupscet(nprocabove)
          endif

c---      setup sobol generator
          ! for now we stride across mpi processes and use omp critical
          if (usesobol) then
              if (world_size > 1) then
                  write (*,*) "Striding sobol with world_size = ", world_size

                  if (world_size <= 2) then
                      strides = 1
                  else if (world_size == 4) then
                      strides = 2
                  else if (world_size == 8) then
                      strides = 3
                  else if (world_size == 16) then
                      strides = 4
                  else if (world_size == 32) then
                      strides = 5
                  else if (world_size == 64) then
                      strides = 6
                  else if (world_size == 128) then
                      strides = 7
                  else if (world_size == 256) then
                      strides = 8
                  else
                      write (*,*) "WARNING: Sobol sequence only usable with 2^n mpi processes."
                      write (*,*) "Falling back to using pseudo-random numbers."
                      usesobol = .false.
                      do part=1,maxParts; do ips=1,maxIps
                          associate (info => iterationStorage(part,ips)%vinfo)
                              info%useSobol = usesobol
                          end associate
                      enddo; enddo
                  endif
              else
                  strides = 0
              endif

              if (ndim > 22) then
                  write (*,*) "sobol direction vectors for ndim > ", ndim, " not implemented"
                  error stop
              endif
          endif ! usesobol

          if (readin) then
              ! for now let's just read the file from all processes
              !if (rank == 0) then
                  call deserializeMCFM()
              !endif
              !call mpi_broadcast_iterationStorage()
          endif

c ================ 
c === warmup phase
c ================ 

          warmdone = .false.

          do part=1,maxParts
              if (computePart(part)) then
                  ipsLordLoop: do ipsgen=1,maxipsgens(part)
                      currentPart = part
                      currentIps = ipsgen

                      if (iterationStorage(currentPart,currentIps)%vinfo%warmupComplete) then
                          cycle ipsLordLoop
                      endif

                      warmdone = .true.

                      fragint_mode = .false.

                      select case (part)
                        case (lord)
                            nproc = nprocbelow
                            abovecut = .false.
                            usescet = .false.
                            kpart = klord
                            coeffonly = origCoeffonly
                            ifun => lowint
                        case (nloReal)
                            nproc = nprocbelow
                            abovecut = .false.
                            usescet = .false.
                            kpart = kreal
                            coeffonly = .false.
                            ifun => realint
                        case (nloVirt)
                            nproc = nprocbelow
                            abovecut = .false.
                            usescet = .false.
                            kpart = kvirt
                            coeffonly = origCoeffonly
                            ifun => virtint
                        case (nloFrag)
                            nproc = nprocbelow
                            abovecut = .false.
                            usescet = .false.
                            kpart = kfrag
                            coeffonly = origCoeffonly
                            fragint_mode = .true.
                            ifun => fragint
                        case (nloRealExtra)
                            nproc = nprocextra
                            abovecut = .false.
                            usescet = .false.
                            kpart = kreal
                            coeffonly = .false.
                            ifun => realint
                        case (snloBelow)
                            nproc = nprocbelow
                            abovecut = .false.
                            usescet = .true.
                            kpart = ksnlo
                            coeffonly = origCoeffonly
                            ifun => scetint
                        case (snloAbove)
                            nproc = nprocabove
                            abovecut = .true.
                            usescet = .true.
                            kpart = klord
                            coeffonly = origCoeffonly
                            ifun => lowint
                        case (nnloBelow)
                            nproc = nprocbelow
                            abovecut = .false.
                            usescet = .true.
                            kpart = knnlo
                            coeffonly = .true.
                            ifun => scetint
                        case (nnloVirtAbove)
                            nproc = nprocabove
                            abovecut = .true.
                            usescet = .true.
                            kpart = kvirt
                            coeffonly = .true.
                            ifun => virtint
                        case (nnloRealAbove)
                            nproc = nprocabove
                            abovecut = .true.
                            usescet = .true.
                            kpart = kreal
                            coeffonly = .true.
                            ifun => realint
                        case default
                            error stop "unknown part"
                        end select

                        reset = .true.
                        scalereset = .true.
                        maxipsgen = maxipsgens(part)

                        call chooser
                        ndim = ndim + ndim_incr(part)

                        associate (info => iterationStorage(currentPart,currentIps)%vinfo )
                            if (info%useSobol) then
                                if (.not. info%sobolInitialized) then
                                    do i=1,ndim+2
                                        call info%sstate(i)%initialize(s(i),a(i),m(:,i), stride=strides)
                                        skiprnd = info%sstate(i)%skip_ahead(info%sobolSkip)
                                        ! we use the sobol generator in strided mode
                                        ! for MPI and "critical" for OMP
                                        do k=0, rank-1
                                            skiprnd = info%sstate(i)%next()
                                        enddo
                                    enddo
                                    info%sobolInitialized = .true.
                                endif
                            endif
                        end associate

                        do
                          if (rank == 0) then
                              write (*,*) ""
                              write (*,'(A,A,I1)') trim(partStrings(part)), " warmup integration, contribution ", currentIps
                              write (*,*) ""
                          endif

                          iterationStorage(currentPart,currentIps)%used = .true.
                          call integrate(ifun, 0, iterBatchWarmup, ndim, iterationStorage(currentPart,currentIps)%vinfo)

                          associate (info => iterationStorage(currentPart,currentIps)%vinfo)
                              if ( info%si == 0._dp ) then
                                  if (rank == 0) then
                                      write (*,*) "Integral zero. Skipping this contribution."
                                  endif
                                  exit
                              else if ( (info%sd() / abs(info%sig()) > warmupPrecisionGoal) ) then
                                  info%callsPerIt = info%callsPerIt * iterCallMult
                                  if (rank == 0) then
                                      write (*,'(A,I3,A)') "Relative warmup precision goal of ",
     &                                      nint(warmupPrecisionGoal*100),
     &                                      " percent not reached"
                                      write (*,*) "Increasing calls to ", info%callsPerIt
                                  endif
                              else if ( info%chisq() > warmupChisqGoal ) then
                                  info%callsPerIt = info%callsPerIt * iterCallMult
                                  if (rank == 0) then
                                      write (*,'(A,F5.3,A)') "Chisq/it warmup precision goal of ",
     &                                      warmupChisqgoal, " not reached"
                                      write (*,*) "Increasing calls to ", info%callsPerIt
                                  endif
                              else
                                  if (rank == 0) then
                                      write (*,*) "Reached warmup precisionGoal with ",
     &                                      info%callsPerIt, " calls per iteration"
                                  endif
                                  ! warmup complete
                                  exit
                              endif
                          end associate
                        enddo

                        ndim = ndim - ndim_incr(part)

                        if ( iterationStorage(currentPart,currentIps)%vinfo%si == 0._dp ) then
                            ! skip main integration
                            iterationStorage(currentPart,currentIps)%vinfo%doFirstCall = .false.
                            iterationStorage(currentPart,currentIps)%used = .false.
                        else
                            iterationStorage(currentPart,currentIps)%vinfo%doFirstCall = .true.
                        endif

                        iterationStorage(currentPart,currentIps)%vinfo%warmupComplete = .true.
                        
                        if (rank == 0) then
                            call serializeMCFM()
                        endif


                  enddo ipsLordLoop

              endif
          enddo

          if (rank == 0) then
              totalsig = 0._dp
              totalsd = 0._dp
              do part=1,maxParts
                do ips=1,maxIps
                    if (iterationStorage(part,ips)%used) then
                        totalsig = totalsig + iterationStorage(part,ips)%vinfo%sig()
                        totalsd = totalsd + iterationStorage(part,ips)%vinfo%sd()**2
                    endif
                enddo
              enddo
              totalsd = sqrt(totalsd)
              if (warmdone) then
                write (*,*) ""
                write (*,*) "WARMUP phase complete"
                write (*,*) "Intermediate warmup result"
              else
                write (*,*) "Result loaded from snapshot"
              endif
                call printcross(totalsig, totalsd, chisqMax())
                write (*,*) ""
          endif

c =========================== 
c === final integration phase
c =========================== 


          ! enable binning of histograms
          bin = .true.

          integrationLoop: do
            nextLoc(:) = 0

            ! first check if there has been a contribution
            ! with no final run, but only warmup integration
            loopContrib: do i=1,maxParts
                loopIpsgen: do j=1,maxIps
                    if (iterationStorage(i,j)%vinfo%doFirstCall) then
                        nextLoc(1) = i
                        nextLoc(2) = j
                        ! do not load result from saved file
                        stage = 1
                        if (rank == 0) then
                            write (*,'(A,A,A,I1)') "first full integration for ", trim(partStrings(i)), " contribution ", j
                        endif
                        iterBatch = iterBatch1
                        exit loopContrib
                    endif
                enddo loopIpsgen
            enddo loopContrib

            ! otherwise pick contribution with biggest sd
            if (nextLoc(1) == 0) then
                sdMod(:,:) = 0._dp
                do i=1,maxParts
                    do j=1,maxIps
                        if (iterationStorage(i,j)%used) then
                            sdMod(i,j) = iterationStorage(i,j)%vinfo%sd()
                        endif
                    enddo
                enddo

                ! we fake the "real emission" contributions uncertainty
                ! for the prioritization system here
                ! to boost the easier to compute contributions
                sdMod(nloReal,:) = sdMod(nloReal,:) / 1.5_dp
                sdMod(snloAbove,:) = sdMod(snloAbove,:) / 1.5_dp
                sdMod(nnloRealAbove,:) = sdMod(nnloRealAbove,:) / 1.5_dp

                nextLoc = maxloc(sdMod(:,:))
                ! load previous result from saved file
                stage = 2

                if (rank == 0) then
                    write (*,'(A,A,I1)') trim(partStrings(nextLoc(1))), " full integration, contribution ", nextLoc(2)
                endif

                iterBatch = iterBatch2
            endif

            currentPart = nextLoc(1)
            currentIps = nextLoc(2)

            if (.not. iterationStorage(currentPart,currentIps)%used) then
                error stop "something went wrong, this part should not be selected for integration"
            endif

            ! boost number of calls w.r.t. to warmup phase
            associate (info => iterationStorage(currentPart,currentIps)%vinfo)
                if (info%doFirstCall) then
                    info%doFirstCall = .false.
                    !info%callsPerIt = info%callsPerIt * 4
                    info%callsPerIt = int(real(info%callsPerIt,dp) * callBoost)
                endif
            end associate

            ! Once we have reached maxCallsPerIter we just run one iteration
            ! per batch to have frequent snapshots and updated results.
            ! We set it here again for resumed runs.
            if (iterationStorage(currentPart,currentIps)%vinfo%callsPerIt*iterCallMult
     &              > maxCallsPerIter) then
              !iterationStorage(currentPart,currentIps)%vinfo%callsPerIt = maxCallsPerIter
              iterBatch = 1
            endif

            ! now perform the full integration, similar to the warmup above

            reset = .true.
            scalereset = .true.
            fragint_mode = .false.

            select case (currentPart)
              case (lord)
                  nproc = nprocbelow
                  abovecut = .false.
                  usescet = .false.
                  kpart = klord
                  coeffonly = origCoeffonly
                  ifun => lowint
              case (nloReal)
                  nproc = nprocbelow
                  abovecut = .false.
                  usescet = .false.
                  kpart = kreal
                  coeffonly = .false.
                  ifun => realint
                  ! bin this contribution to all tau bins
                  if (origKpart == knnlo) then
!$omp parallel
                      scetreweight(:) = 1._dp
!$omp end parallel
                  endif
              case (nloVirt)
                  nproc = nprocbelow
                  abovecut = .false.
                  usescet = .false.
                  kpart = kvirt
                  coeffonly = origCoeffonly
                  ifun => virtint
                  ! bin this contribution to all tau bins
                  if (origKpart == knnlo) then
!$omp parallel
                      scetreweight(:) = 1._dp
!$omp end parallel
                  endif
              case (nloFrag)
                  nproc = nprocbelow
                  abovecut = .false.
                  usescet = .false.
                  kpart = kfrag
                  coeffonly = origCoeffonly
                  fragint_mode = .true.
                  ifun => fragint
              case (nloRealExtra)
                  nproc = nprocextra
                  abovecut = .false.
                  usescet = .false.
                  kpart = kreal
                  coeffonly = .false.
                  ifun => realint
              case (snloBelow)
                  nproc = nprocbelow
                  abovecut = .false.
                  usescet = .true.
                  kpart = ksnlo
                  coeffonly = origCoeffonly
                  ifun => scetint
              case (snloAbove)
                  nproc = nprocabove
                  abovecut = .true.
                  usescet = .true.
                  kpart = klord
                  coeffonly = origCoeffonly
                  ifun => lowint
              case (nnloBelow)
                  nproc = nprocbelow
                  abovecut = .false.
                  usescet = .true.
                  kpart = knnlo
                  coeffonly = .true.
                  ifun => scetint
              case (nnloVirtAbove)
                  nproc = nprocabove
                  abovecut = .true.
                  usescet = .true.
                  kpart = kvirt
                  coeffonly = .true.
                  ifun => virtint
              case (nnloRealAbove)
                  nproc = nprocabove
                  abovecut = .true.
                  usescet = .true.
                  kpart = kreal
                  coeffonly = .true.
                  ifun => realint
              case default
                  error stop "unknown part"
              end select

              ! sanity setup. this includeTaucutgrid is initialized to .true.
              ! in parseinput.f. But when run in a full nnlo calculation
              ! after some scet runs the real integration might run again
              ! with includeTaucutgrid modified.
!              if (usescet) then
!$omp parallel
              includeTaucutgrid(:) = .true.
!$omp end parallel

!              endif

              reset = .true.
              scalereset = .true.

              maxipsgen = maxipsgens(currentPart)
              ipsgen = currentIps

              call chooser
              ndim = ndim + ndim_incr(currentPart)

              associate (info => iterationStorage(currentPart,currentIps)%vinfo )
                  if (info%useSobol) then
                      if (.not. info%sobolInitialized) then
                          do i=1,ndim+2
                              call info%sstate(i)%initialize(s(i),a(i),m(:,i), stride=strides)
                              skiprnd = info%sstate(i)%skip_ahead(info%sobolSkip)
                              ! we use the sobol generator in strided mode
                              ! for MPI and "critical" for OMP
                              do k=0, rank-1
                                  skiprnd = info%sstate(i)%next()
                              enddo
                          enddo
                      endif
                  endif
              end associate

              iterationStorage(currentPart,currentIps)%used = .true.
              call integrate(ifun, stage, iterBatch, ndim, iterationStorage(currentPart,currentIps)%vinfo)
              ndim = ndim - ndim_incr(currentPart)

              iterationStorage(currentPart,currentIps)%vinfo%callsPerIt = 
     &              iterationStorage(currentPart,currentIps)%vinfo%callsPerIt * iterCallMult

              if (rank == 0) then
                  call serializeMCFM()
              endif
             

              ! determine if we had at least one full integration for all parts
              ! then show final results, save histogram
              ! and exit if precision goal reached
              if (all(iterationStorage(:,:)%vinfo%doFirstCall .eqv. .false.)) then
                  if (rank ==0) then
                      call finalizeStorage
                      ! generate and save user readable output
                  endif

                  totalsig = 0._dp
                  totalsd = 0._dp
                  do part=1,maxParts
                    do ips=1,maxIps
                        if (iterationStorage(part,ips)%used) then
                            totalsig = totalsig + iterationStorage(part,ips)%vinfo%sig()
                            totalsd = totalsd + iterationStorage(part,ips)%vinfo%sd()**2
                        endif
                    enddo
                  enddo
                  totalsd = sqrt(totalsd)

                  if (abs(totalsd/totalsig) < resultPrecisionGoal) then
                      exit integrationLoop
                  else ! print intermediate results
                      if (rank == 0 .and. writeintermediate) then
                          write (*,*) "Intermediate full result"

#ifdef HAVE_LHAPDF
                          call printallcross()

                          if (doPDFerrors) then
                            call printPDFuncertainties()
                          endif
#else
                          call printcross(totalsig, totalsd, chisqMax() )
#endif

                          if (doScalevar) then
                              call printScaleuncertainties()
                          endif

                          if (doMultitaucut) then
                              call printTaucuts()
                          endif

                          call writeAllHistograms()
                      endif

                      if (rank == 0) then
                          ! intermediate benchmark result
                          benchmarkinter: block
                              use MCFMBenchmark
                              integer :: exitCode
                              if (cfg_var_configadded(cfg, "extra%benchmark")) then
                                  exitCode = comparisonCode(1)
                              endif
                          end block benchmarkinter
                      endif
                  endif
              endif

              ! save snapshots

          enddo integrationLoop

          integ = totalsig
          integ_err = totalsd

      end subroutine

