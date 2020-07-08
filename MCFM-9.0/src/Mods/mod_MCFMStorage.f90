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
 
! storage of cross section information
! with the purpose of vegas integration
module MCFMStorage
    use types
    use Superhisto
    use Integration
    use Scalevar
    use PDFerrors
    use SCET
    use omp_lib
    use mod_sobseq
    use iso_fortran_env
    implicit none

    type IntegrationInfo
        ! mandatory raw vegas data to resume integration
        real(dp) :: si ! sum I_i / sigma_i^2
        real(dp) :: swgt ! sum 1 / sigma_i^2
        real(dp) :: schi
        integer :: lastIter
        real(dp), allocatable :: grid(:,:)
        ! further parameters
        logical :: doFirstCall, warmupComplete
        integer(int64) :: callsPerIt
        real(dp) :: damp(2)
        logical :: useSobol
        logical :: sobolInitialized
        integer(int64) :: sobolSkip
        type(sobol_state) :: sstate(25)
        contains
            procedure :: sd => info_sd
            procedure :: sig => info_sig
            procedure :: chisq => info_chisq
            procedure :: serialize => serializeIntegrationInfo
            procedure :: deserialize => deserializeIntegrationInfo
    end type

    type HistogramStorage
        type(sh_histogram), allocatable :: histos(:)
        contains
            procedure :: serialize => serializeHistogramStorage
            procedure :: deserialize => deserializeHistogramStorage
            procedure :: reset => histogramstorage_reset
    end type

    type PartStorage
        logical :: used
        type(IntegrationInfo) :: vinfo
        type(HistogramStorage) :: histCentral
        type(HistogramStorage), allocatable :: histPDFerrors(:)
        type(HistogramStorage), allocatable :: histScalevar(:)
        type(HistogramStorage), allocatable :: histTaucut(:)
      contains
          procedure :: reset => partstorage_reset
          procedure :: printme => partstorage_print
          procedure :: serialize => serializePartStorage
          procedure :: deserialize => deserializePartStorage
    end type

    ! saves OMP threadlocal integration information
    type(PartStorage), public, save :: threadStorage(maxParts,maxIps)
!$omp threadprivate(threadStorage)

    ! accumulated integration information for each iteration
    ! from OMP threads and MPI processes
    type(PartStorage), public, save :: masterStorage(maxParts,maxIps)

    ! accumulated results over multiple iterations
    ! xx in histograms will accumulate I_i / sigma_i^2
    ! xxsq in histograms will accumulate 1 / sigma_i^2

    ! This is the object that gets dumped and loaded
    ! for resumable integrations. It has valid data on the
    ! mpi master process.
    type(PartStorage), public, save :: iterationStorage(maxParts,maxIps)
    
    ! contains final results for presentation
    type(PartStorage), public, save :: finalSum

    logical, public, save :: storageAllocated = .false.
!$omp threadprivate(storageAllocated)

    integer, public, save :: currentPart
    integer, public, save :: currentIps

    public :: IntegrationInfo
    public :: HistogramStorage
    public :: integrate
    public :: finalizeStorage
    public :: threadStorageOp
    public :: initMasterStorage
    public :: initHistogramStorage
    public :: mpi_broadcast_iterationStorage

    public :: serializeMCFM
    public :: deserializeMCFM

    ! fixed parameter for now, we allow no dynamic resizing
    integer, save, public :: ndmx = 100
    private

    contains

    subroutine serializeMCFM()
        implicit none
        integer :: i,j
        integer :: ierr
        character(len=255) :: imsg

        character(len=255) :: runname
        common/runname/runname

        open(unit=11, file=trim(runname)//"_snapshot.dat", status='replace', &
            form='unformatted', iostat=ierr, iomsg=imsg)
        if (ierr == 0) then
            do i=1,maxParts
                do j=1,maxIps
                    call iterationStorage(i,j)%serialize(11)
                enddo
            enddo
            close(unit=11)

            write (*,*) ""
            write (*,*) "Snapshot written to " // trim(runname) // "_snapshot.dat"
            write (*,*) ""
        else
            write (*,*) "Problem writing snapshot file " // trim(runname) // "_snapshot.dat"
            write (*,*) trim(imsg)
            write (*,*) "Error code = ", ierr
        endif

    end subroutine

    subroutine deserializeMCFM()
        implicit none
        integer :: i,j
        integer :: ierr
        character(len=255) :: imsg

        character(len=255) :: runname
        common/runname/runname

        open(unit=11, file=trim(runname)//"_snapshot.dat", status='old', &
            form='unformatted', iostat=ierr, iomsg=imsg)
        if (ierr == 0) then
            do i=1,maxParts
                do j=1,maxIps
                call iterationStorage(i,j)%deserialize(11)
                enddo
            enddo
            close(unit=11)

            write (*,*) "Snapshot read from " // trim(runname) // "_snapshot.dat"
        else
            write (*,*) "Could not read snapshot file " // trim(runname) // "_snapshot.dat"
            write (*,*) trim(imsg)
            write (*,*) "Error code = ", ierr
        endif

    end subroutine

    subroutine mpi_broadcast_iterationStorage()
        implicit none
        include 'mpicommon.f'

        integer :: k
        integer :: ierr

        error stop "to do"

    end subroutine

    ! this reduces all histograms from masterStorage
#ifdef HAVE_MPI
    subroutine mpi_reduce_masterStorage()
        use Scalevar
        use PDFerrors
        use SCET
        use mpi
        implicit none
        include 'mpicommon.f'

        integer :: j,k
        integer :: ierr

        if (masterStorage(currentPart,currentIps)%used .eqv. .true.) then
            do k=1,size(masterStorage(currentPart,currentIps)%histCentral%histos)
                if (masterStorage(currentPart,currentIps)%histCentral%histos(k)%initialized()) then
                    call mpi_allreduce(mpi_in_place, &
                        masterStorage(currentPart,currentIps)%histCentral%histos(k)%xx, &
                        size(masterStorage(currentPart,currentIps)%histCentral%histos(k)%xx), &
                        mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
                    call mpi_allreduce(mpi_in_place, &
                        masterStorage(currentPart,currentIps)%histCentral%histos(k)%xxsq, &
                        size(masterStorage(currentPart,currentIps)%histCentral%histos(k)%xxsq), &
                        mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
                endif
            enddo

            if (doScalevar) then
                do j=1,maxScalevar 
                    do k=1,size(masterStorage(currentPart,currentIps)%histScalevar(j)%histos)
                        if (masterStorage(currentPart,currentIps)%histScalevar(j)%histos(k)%initialized()) then
                            call mpi_allreduce(mpi_in_place, &
                                masterStorage(currentPart,currentIps)%histScalevar(j)%histos(k)%xx, &
                                size(masterStorage(currentPart,currentIps)%histScalevar(j)%histos(k)%xx), &
                                mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
                            call mpi_allreduce(mpi_in_place, &
                                masterStorage(currentPart,currentIps)%histScalevar(j)%histos(k)%xxsq, &
                                size(masterStorage(currentPart,currentIps)%histScalevar(j)%histos(k)%xxsq), &
                                mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
                        endif
                    enddo
                enddo
            endif

            if (maxPDFsets > 0) then
                do j=1,maxPDFsets
                    do k=1,size(masterStorage(currentPart,currentIps)%histPDFerrors(j)%histos)
                        if (masterStorage(currentPart,currentIps)%histPDFerrors(j)%histos(k)%initialized()) then
                            call mpi_allreduce(mpi_in_place, &
                                masterStorage(currentPart,currentIps)%histPDFerrors(j)%histos(k)%xx, &
                                size(masterStorage(currentPart,currentIps)%histPDFerrors(j)%histos(k)%xx), &
                                mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
                            call mpi_allreduce(mpi_in_place, &
                                masterStorage(currentPart,currentIps)%histPDFerrors(j)%histos(k)%xxsq, &
                                size(masterStorage(currentPart,currentIps)%histPDFerrors(j)%histos(k)%xxsq), &
                                mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
                        endif
                    enddo
                enddo
            endif

            do j=1,size(tcutarray)
                do k=1,size(masterStorage(currentPart,currentIps)%histTaucut(j)%histos)
                    if (masterStorage(currentPart,currentIps)%histTaucut(j)%histos(k)%initialized()) then
                        call mpi_allreduce(mpi_in_place, &
                            masterStorage(currentPart,currentIps)%histTaucut(j)%histos(k)%xx, &
                            size(masterStorage(currentPart,currentIps)%histTaucut(j)%histos(k)%xx), &
                            mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
                        call mpi_allreduce(mpi_in_place, &
                            masterStorage(currentPart,currentIps)%histTaucut(j)%histos(k)%xxsq, &
                            size(masterStorage(currentPart,currentIps)%histTaucut(j)%histos(k)%xxsq), &
                            mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
                    endif
                enddo
            enddo
        endif
    end subroutine
#endif

    ! this can be called any time to populate finalSum
    ! from iterationStorage
    subroutine finalizeStorage()
        use SCET
        implicit none

        ! process iterationStorage for presentation
        integer :: nplotmax
        common/nplotmax/nplotmax
        integer :: i,j,k,l,m

        integer :: iused, jused


        ! find a storage block that has been used to initialize finalSum
        ! and also count number of contributions
        do i=1,maxParts
            do j=1,maxIps
                if (iterationStorage(i,j)%used .eqv. .true.) then
                    iused = i
                    jused = j
                endif
            enddo
        enddo

        finalSum = iterationStorage(iused,jused)
        call finalSum%reset

        ! this computes iteration averaged, weighted, integral and standard deviation in xx, xxsq
        do i=1,maxParts
            do j=1,maxIps
                if (iterationStorage(i,j)%used .eqv. .true.) then
                    do k=1,nplotmax
                        associate ( histo => iterationStorage(i,j)%histCentral%histos(k), &
                                    fin => finalSum%histCentral%histos(k) )
                            ! some histograms might not have been initialized
                            ! because they never received a binning
                            ! copy over the first initialized one
                            if ((fin%initialized() .eqv. .false.) .and. (histo%initialized() .eqv. .true.)) then
                                ! sum not initialized, but contribution
                                fin = histo
                                call fin%reset()
                            endif
                            if (histo%initialized()) then
                                do l=0,fin%nbins+1
                                    if (histo%xx(l) /= 0._dp) then
                                        fin%xx(l) = fin%xx(l) + histo%xx(l) / histo%xxsq(l)
                                        fin%xxsq(l) = fin%xxsq(l) + 1._dp/histo%xxsq(l)
                                    endif
                                enddo
                            endif
                        end associate

                        if (doScalevar) then
                            do m=1,maxscalevar
                                associate (histo => iterationStorage(i,j)%histScalevar(m)%histos(k), &
                                        fin => finalSum%histScalevar(m)%histos(k) )
                                    if ((fin%initialized() .eqv. .false.) .and. (histo%initialized() .eqv. .true.)) then
                                        ! sum not initialized, but contribution
                                        fin = histo
                                    endif
                                    if (histo%initialized()) then
                                        do l=0,fin%nbins+1
                                            if (histo%xx(l) /= 0._dp) then
                                                fin%xx(l) = fin%xx(l) + histo%xx(l) / histo%xxsq(l)
                                                fin%xxsq(l) = fin%xxsq(l) + 1._dp/histo%xxsq(l)
                                            endif
                                        enddo
                                    endif
                                end associate
                            enddo
                        endif

                        if (maxPDFsets > 0) then
                            do m=1,maxPDFsets
                                associate (histo => iterationStorage(i,j)%histPDFerrors(m)%histos(k), &
                                        fin => finalSum%histPDFerrors(m)%histos(k) )
                                    if ((fin%initialized() .eqv. .false.) .and. (histo%initialized() .eqv. .true.)) then
                                        ! sum not initialized, but contribution
                                        fin = histo
                                    endif
                                    if (histo%initialized()) then
                                        do l=0,fin%nbins+1
                                            if (histo%xx(l) /= 0._dp) then
                                                fin%xx(l) = fin%xx(l) + histo%xx(l) / histo%xxsq(l)
                                                fin%xxsq(l) = fin%xxsq(l) + 1._dp/histo%xxsq(l)
                                            endif
                                        enddo
                                    endif
                                end associate
                            enddo
                        endif

                        do m=1,size(tcutarray)
                            associate (histo => iterationStorage(i,j)%histTaucut(m)%histos(k), &
                                    fin => finalSum%histTaucut(m)%histos(k) )
                               if ((fin%initialized() .eqv. .false.) .and. (histo%initialized() .eqv. .true.)) then
                                   ! sum not initialized, but contribution
                                   fin = histo
                               endif
                                if (histo%initialized()) then
                                    do l=0,fin%nbins+1
                                        if (histo%xx(l) /= 0._dp) then
                                            fin%xx(l) = fin%xx(l) + histo%xx(l) / histo%xxsq(l)
                                            fin%xxsq(l) = fin%xxsq(l) + 1._dp/histo%xxsq(l)
                                        endif
                                    enddo
                                endif
                            end associate
                        enddo
                    enddo
                endif
            enddo
        enddo

        ! variance to standard deviation
        do k=1,nplotmax
            associate ( fin => finalSum%histCentral%histos(k) )
                if (fin%initialized()) then
                    fin%xxsq = sqrt(fin%xxsq)
                endif
            end associate

            if (doScalevar) then
                do m=1,maxscalevar
                    associate ( fin => finalSum%histScalevar(m)%histos(k) )
                        if (fin%initialized()) then
                            fin%xxsq = sqrt(fin%xxsq)
                        endif
                    end associate
                enddo
            endif

            if (maxPDFsets > 0) then
                do m=1,maxPDFsets
                    associate ( fin => finalSum%histPDFerrors(m)%histos(k) )
                        if (fin%initialized()) then
                            fin%xxsq = sqrt(fin%xxsq)
                        endif
                    end associate
                enddo
            endif

            do m=1,size(tcutarray)
                associate ( fin => finalSum%histTaucut(m)%histos(k) )
                    if (fin%initialized()) then
                        fin%xxsq = sqrt(fin%xxsq)
                    endif
                end associate
            enddo
        enddo

    end subroutine

    subroutine partstorage_print(this)
        use SCET
        implicit none
        class(PartStorage), intent(in) :: this

        integer :: nplotmax
        !common/nplotmax/nplotmax
        integer :: k,j

        nplotmax = 1

        do k=1,nplotmax
            if (shinitialized(this%histCentral%histos(k))) then
                call this%histCentral%histos(k)%printme
            endif
        enddo

        if (doScalevar) then
            do j=1,maxscalevar
                do k=1,nplotmax
                    if (shinitialized(this%histScalevar(j)%histos(k))) then
                        call this%histScalevar(j)%histos(k)%printme
                    endif
                enddo
            enddo
        endif

        if (maxPDFsets > 0) then
            do j=1,maxPDFsets
                do k=1,nplotmax
                    if (shinitialized(this%histPDFerrors(j)%histos(k))) then
                        call this%histPDFerrors(j)%histos(k)%printme
                    endif
                enddo
            enddo
        endif

        do j=1,size(tcutarray)
            do k=1,nplotmax
                if (shinitialized(this%histTaucut(j)%histos(k))) then
                    call this%histTaucut(j)%histos(k)%printme
                endif
            enddo
        enddo

    end subroutine

    subroutine histogramstorage_reset(this)
        implicit none
        class(HistogramStorage), intent(inout) :: this

        integer :: j

        do j=1,size(this%histos)
            if (this%histos(j)%initialized()) then
                call this%histos(j)%reset()
            endif
        enddo

    end subroutine

    subroutine partstorage_reset(this)
        use SCET
        implicit none
        class(PartStorage), intent(inout) :: this
        integer :: nplotmax
        common/nplotmax/nplotmax
        integer :: k,j

        call this%histCentral%reset()

        if (doScalevar) then
            do j=1,maxscalevar
                call this%histScalevar(j)%reset()
            enddo
        endif 

        if (maxPDFsets > 0) then
            do j=1,maxPDFsets
                call this%histPDFerrors(j)%reset()
            enddo
        endif

        do j=1,size(tcutarray)
            call this%histTaucut(j)%reset()
        enddo

    end subroutine

    subroutine serializePartStorage(this, unit)
        use SCET
        implicit none
        class(PartStorage), intent(in) :: this
        integer, intent(in) :: unit
        integer :: j

        write (unit) this%used
        call this%vinfo%serialize(unit)
        call this%histCentral%serialize(unit)
        write (unit) doPDFerrors
        if (doPDFerrors) then
            write (unit) size(this%histPDFerrors)
            do j=1,size(this%histPDFerrors)
                call this%histPDFerrors(j)%serialize(unit)
            enddo
        endif
        write (unit) doScalevar
        if (doScalevar) then
            write (unit) size(this%histScalevar)
            do j=1,size(this%histScalevar)
                call this%histScalevar(j)%serialize(unit)
            enddo
        endif
        write (unit) size(tcutarray)
        do j=1,size(tcutarray)
            call this%histTaucut(j)%serialize(unit)
        enddo

    end subroutine

    subroutine deserializePartStorage(this, unit)
        use SCET
        implicit none
        class(PartStorage), intent(inout) :: this
        integer, intent(in) :: unit
        logical :: read_pdf, read_scale
        integer :: npdf, nscale
        integer :: j
        integer :: tcutarraysize

        read (unit) this%used
        call this%vinfo%deserialize(unit)
        call this%histCentral%deserialize(unit)
        read (unit) read_pdf
        if (read_pdf .neqv. doPDFerrors) then
            write (*,*) "WARNING: snapshot contains PDF errors, &
                &but input file does not enable them"
        endif
        if (read_pdf) then
            read (unit) npdf
            if (npdf /= maxPDFsets) then
                write (*,'(A,I1,A,I1)') "WARNING: number of PDF error sets read is ", &
                    npdf, "and does not match current value from input file of ", maxPDFsets
            endif
            do j=1,npdf
                call this%histPDFerrors(j)%deserialize(unit)
            enddo
        endif
        read (unit) read_scale
        if (read_scale .neqv. doScalevar) then
            write (*,*) "WARNING: snapshot contains scale uncertainties, &
                &but input file does not enable them"
        endif
        if (read_scale) then
            read (unit) nscale
            if (nscale /= maxscalevar) then
                write (*,'(A,I1,A,I1)') "WARNING: number of scale variation points read is ", &
                    nscale, "and does not match current value in input file of ", maxscalevar
            endif
            do j=1,nscale
                call this%histScalevar(j)%deserialize(unit)
            enddo
        endif

        read (unit) tcutarraysize
        if (tcutarraysize /= size(tcutarray)) then
            write (*,'(A,I1,A,I1)') "WARNING: number of taucut values in tcutarray is ", &
                tcutarraysize, "and does not match current value in input file of ", size(tcutarray)
        endif
        do j=1,tcutarraysize
            call this%histTaucut(j)%deserialize(unit)
        enddo
    end subroutine

    subroutine accumulateIteration(ncall)
        use SCET
        use iso_fortran_env
        implicit none

        integer(int64), intent(in) :: ncall
        integer :: nplotmax
        common/nplotmax/nplotmax
        integer :: k,j

        !write (*,*) "Accumulating iteration for thread ", omp_get_thread_num(), &
            !" with calls = ", ncall

        associate ( itstor => iterationStorage(currentPart,currentIps), &
                    mastor => masterStorage(currentpart,currentIps) )
            do k=1,nplotmax
                if (itstor%histCentral%histos(k)%initialized()) then
                    call itstor%histCentral%histos(k)%iterprocess( &
                        mastor%histCentral%histos(k), real(ncall,dp))
                endif
            enddo

            if (doScalevar) then
                do j=1,maxscalevar
                    !write (*,*)  "Scale thread", j, masterStorage(currentPart,currentIps)%histScalevar(j)%histos(1)%xx(1)
                    do k=1,nplotmax
                        if (itstor%histScalevar(j)%histos(k)%initialized()) then
                            call itstor%histScalevar(j)%histos(k)%iterprocess( &
                                mastor%histScalevar(j)%histos(k), real(ncall,dp))
                        endif
                    enddo
                enddo
            endif

            if (maxPDFsets > 0) then
                do j=1,maxPDFsets
                    !write (*,*)  "PDF thread", j, masterStorage(currentPart,currentIps)%histPDFerrors(j)%histos(1)%xx(1)
                    do k=1,nplotmax
                        if (itstor%histPDFerrors(j)%histos(k)%initialized()) then
                            call itstor%histPDFerrors(j)%histos(k)%iterprocess( &
                                mastor%histPDFerrors(j)%histos(k), real(ncall,dp))
                        endif
                    enddo
                enddo
            endif

            do j=1,size(tcutarray)
                do k=1,nplotmax
                    if (itstor%histTaucut(j)%histos(k)%initialized()) then
                        call itstor%histTaucut(j)%histos(k)%iterprocess( &
                            mastor%histTaucut(j)%histos(k), real(ncall,dp))
                    endif
                enddo
            enddo

            call mastor%reset
        end associate
    end subroutine

    subroutine threadStorageOp(oper)
        implicit none
        integer :: nplotmax
        common/nplotmax/nplotmax
        integer :: k,j

        interface
            subroutine oper(hist)
                use Superhisto
                implicit none
                class(sh_histogram), intent(inout) :: hist
            end subroutine
        end interface

        do k=1,nplotmax
            if ( shinitialized(threadStorage(currentPart,currentIps)%histCentral%histos(k)) ) then
                call oper(threadStorage(currentPart,currentIps)%histCentral%histos(k))
            endif
        enddo

        if (doScalevar) then
            do j=1,maxscalevar
                do k=1,nplotmax
                    if ( shinitialized(threadStorage(currentPart,currentIps)%histScalevar(j)%histos(k)) ) then
                        call oper(threadStorage(currentPart,currentIps)%histScalevar(j)%histos(k))
                    endif
                enddo
            enddo
        endif

        if (maxPDFsets > 0) then
            do j=1,maxPDFsets
                do k=1,nplotmax
                    if ( shinitialized(threadStorage(currentPart,currentIps)%histPDFerrors(j)%histos(k)) ) then
                        call oper(threadStorage(currentPart,currentIps)%histPDFerrors(j)%histos(k))
                    endif
                enddo
            enddo
        endif

        do j=1,size(tcutarray)
            do k=1,nplotmax
                if ( shinitialized(threadStorage(currentPart,currentIps)%histTaucut(j)%histos(k)) ) then
                    call oper(threadStorage(currentPart,currentIps)%histTaucut(j)%histos(k))
                endif
            enddo
        enddo

    end subroutine

    subroutine accumulateThreadStorage()
        implicit none
        integer :: nplotmax
        common/nplotmax/nplotmax
        integer :: k,j

        do k=1,nplotmax
            ! add to accumulator
            if (masterStorage(currentPart,currentIps)%histCentral%histos(k)%initialized()) then
                call masterStorage(currentPart,currentIps)%histCentral%histos(k)%add( &
                    threadStorage(currentPart,currentIps)%histCentral%histos(k))
            endif
        enddo


        if (doScalevar) then
            !write (*,*) "scale variation accumulate thread storage"
            do j=1,maxscalevar
                do k=1,nplotmax
                    if (masterStorage(currentPart,currentIps)%histScalevar(j)%histos(k)%initialized()) then
                        call masterStorage(currentPart,currentIps)%histScalevar(j)%histos(k)%add( &
                            threadStorage(currentPart,currentIps)%histScalevar(j)%histos(k))
                    endif
                enddo
            enddo
        endif

        if (maxPDFsets > 0) then
            do j=1,maxPDFsets
                do k=1,nplotmax
                    if (masterStorage(currentPart,currentIps)%histPDFerrors(j)%histos(k)%initialized()) then
                        call masterStorage(currentPart,currentIps)%histPDFerrors(j)%histos(k)%add( &
                            threadStorage(currentPart,currentIps)%histPDFerrors(j)%histos(k))
                    endif
                enddo
            enddo
        endif

        do j=1,size(tcutarray)
            do k=1,nplotmax
                if (masterStorage(currentPart,currentIps)%histTaucut(j)%histos(k)%initialized()) then
                    call masterStorage(currentPart,currentIps)%histTaucut(j)%histos(k)%add( &
                        threadStorage(currentPart,currentIps)%histTaucut(j)%histos(k))
                endif
            enddo
        enddo

        call threadStorage(currentPart,currentIps)%reset

    end subroutine

    subroutine initMasterStorage(nplotmax)
        implicit none
        include 'scalevar.f'
        integer :: j,k,m
        integer, intent(in) :: nplotmax

        !write (*,*) "Allocating masterStorage for ", nplotmax, " histograms"

        do j=1,maxParts
            do k=1,maxIps
                allocate(masterStorage(j,k)%histCentral%histos(nplotmax))
                allocate(iterationStorage(j,k)%histCentral%histos(nplotmax))

                if (doScalevar) then
                    allocate(iterationStorage(j,k)%histScalevar(maxscalevar))
                    do m=1,maxscalevar
                        allocate(iterationStorage(j,k)%histScalevar(m)%histos(nplotmax))
                    enddo

                    allocate(masterStorage(j,k)%histScalevar(maxscalevar))
                    do m=1,maxscalevar
                        allocate(masterStorage(j,k)%histScalevar(m)%histos(nplotmax))
                    enddo
                endif

                if (maxPDFsets > 0) then
                    allocate(iterationStorage(j,k)%histPDFerrors(maxPDFsets))
                    do m=1,maxPDFsets
                        allocate(iterationStorage(j,k)%histPDFerrors(m)%histos(nplotmax))
                    enddo

                    allocate(masterStorage(j,k)%histPDFerrors(maxPDFsets))
                    do m=1,maxPDFsets
                        allocate(masterStorage(j,k)%histPDFerrors(m)%histos(nplotmax))
                    enddo
                endif

                allocate(iterationStorage(j,k)%histTaucut(size(tcutarray)))
                do m=1,size(tcutarray)
                    allocate(iterationStorage(j,k)%histTaucut(m)%histos(nplotmax))
                enddo
                allocate(masterStorage(j,k)%histTaucut(size(tcutarray)))
                do m=1,size(tcutarray)
                    allocate(masterStorage(j,k)%histTaucut(m)%histos(nplotmax))
                enddo

            enddo
        enddo

    end subroutine

    subroutine initHistogramStorage(nplotmax)
        implicit none
        include 'scalevar.f'
        integer :: j,k,m
        integer, intent(in) :: nplotmax

        !write (*,*) "Allocating threadStorage for ", nplotmax, " histograms"

        do j=1,maxParts
            do k=1,maxIps
                allocate(threadStorage(j,k)%histCentral%histos(nplotmax))

                if (doScalevar) then
                    !write (*,*) "Allocating threadStorage for scale variation", maxscalevar
                    allocate(threadStorage(j,k)%histScalevar(maxscalevar))
                    do m=1,maxscalevar
                        allocate(threadStorage(j,k)%histScalevar(m)%histos(nplotmax))
                    enddo
                endif

                if (maxPDFsets > 0) then
                    !write (*,*) "Allocating threadStorage for PDF errors", maxPDFsets
                    allocate(threadStorage(j,k)%histPDFerrors(maxPDFsets))
                    do m=1,maxPDFsets
                        allocate(threadStorage(j,k)%histPDFerrors(m)%histos(nplotmax))
                    enddo
                endif

                allocate(threadStorage(j,k)%histTaucut(size(tcutarray)))
                do m=1,size(tcutarray)
                    allocate(threadStorage(j,k)%histTaucut(m)%histos(nplotmax))
                enddo
            enddo
        enddo


        storageAllocated = .true.

    end subroutine

    subroutine serializeHistogramStorage(stor, unit)
        implicit none
        class(HistogramStorage), intent(in) :: stor
        integer, intent(in) :: unit
        integer :: j

        write (unit) size(stor%histos)
        do j=1,size(stor%histos)
            call stor%histos(j)%serialize(unit)
        enddo

    end subroutine

    subroutine deserializeHistogramStorage(stor, unit)
        implicit none
        class(HistogramStorage), intent(out) :: stor
        integer, intent(in) :: unit
        integer :: hnum, j

        read (unit) hnum
        allocate(stor%histos(hnum))
        do j=1,hnum
            call stor%histos(j)%deserialize(unit)
        enddo

    end subroutine

!!! =============================
!!! VEGAS integration below
!!! =============================

    subroutine serializeIntegrationInfo(info, unit)
        implicit none
        class(IntegrationInfo), intent(in) :: info
        integer, intent(in) :: unit
        integer :: gridshape(2)

        write (unit) info%si, info%swgt, info%schi
        write (unit) info%lastIter
        ! assume that if lastIter > 0, the grid has been allocated
        if (info%lastIter > 0) then
            gridshape = shape(info%grid)
            write (unit) gridshape
            write (unit) info%grid
        endif
        write (unit) info%doFirstCall, info%warmupComplete
        write (unit) info%callsPerIt
        write (unit) info%damp
        write (unit) info%useSobol
        write (unit) info%sobolSkip
    end subroutine

    subroutine deserializeIntegrationInfo(info, unit)
        implicit none
        class(IntegrationInfo), intent(out) :: info
        integer, intent(in) :: unit
        integer :: gridshape(2)

        read (unit) info%si, info%swgt, info%schi
        read (unit) info%lastIter
        if (info%lastIter > 0) then
            read (unit) gridshape
            allocate (info%grid(gridshape(1), gridshape(2)))
            read (unit) info%grid
        endif
        read (unit) info%doFirstCall, info%warmupComplete
        read (unit) info%callsPerIt
        read (unit) info%damp
        read (unit) info%useSobol
        read (unit) info%sobolSkip

        info%sobolInitialized = .false.

    end subroutine

    function info_sd(this)
        use types
        implicit none
        real(dp) :: info_sd
        class(IntegrationInfo), intent(in) :: this

        info_sd =  1._dp/sqrt(this%swgt)
    end function

    elemental function info_chisq(info)
        use types
        implicit none
        real(dp) :: info_chisq
        class(IntegrationInfo), intent(in) :: info

        if (info%lastIter == 0) then
            info_chisq = 0._dp
        else
            info_chisq = (info%schi - info%si* info%si / info%swgt) / (info%lastIter - 1._dp)
        endif
    end function

    function info_sig(this)
        use types
        implicit none
        real(dp) :: info_sig
        class(IntegrationInfo), intent(in) :: this

        info_sig = this%si / this%swgt
    end function

    subroutine uniformgrid(grid)
        use types
        implicit none
        real(dp), intent(inout) :: grid(:,:)

        integer :: j

        do j=1,ndmx
            grid(j,:) = (1._dp / ndmx ) * j
        enddo
    end subroutine

    subroutine rebin(mi, grid)
        use types
        implicit none

        real(dp), intent(in) :: mi(ndmx)
        real(dp), intent(inout) :: grid(ndmx)
        real(dp) :: gridnew(ndmx)

        integer :: i,k
        real(dp) :: x0, xn, delmi

        delmi = 0._dp
        xn = 0._dp
        x0 = 0._dp
        k = 0

        do i=1,ndmx-1
            do while (delmi < 1._dp/ndmx)
                k = k + 1
                delmi = delmi + mi(k)
                x0 = xn
                xn = grid(k)
            enddo
            delmi = delmi - 1._dp/ndmx
            gridnew(i) = xn - (xn-x0)*delmi/mi(k)
        enddo

        gridnew(ndmx) = 1._dp
        grid(:) = gridnew(:)
    end subroutine

    subroutine integrate(fxn, stage, niter, ndim, info)
        use types
        use cxx11random
        use iso_fortran_env
#ifdef HAVE_MPI
        use mpi
#endif
        use CPUTime
        implicit none
        type(IntegrationInfo), intent(inout) :: info
        integer, intent(in) :: stage, ndim, niter
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

        include 'maxwt.f'
        include 'mpicommon.f'

        ! these are omp copyin
        include 'mxpart.f'
        include 'nf.f'
        include 'phasemin.f'
        include 'cutoff.f'
        include 'jetcuts.f'
        include 'breit.f'
        include 'zerowidth.f'
        include 'srdiags.f'
        include 'interference.f'
        include 'qcdcouple.f'
        include 'ewcouple.f'
        include 'masses.f'
        include 'facscale.f'
        include 'scale.f'
        include 'stopscales.f'
        real(dp) :: p1ext(4), p2ext(4)
        common/pext/p1ext,p2ext
!$omp threadprivate(/pext/)
        include 'lc.f'
        include 'ptilde.f'
        include 'bitflags.f'
        include 'flags.f'
        include 'lastphot.f'
        include 'lhcb.f'
        include 'b0.f'
        include 'swapxz.f'
        include 'nodecay.f'
        include 'heavyflav.f'
        include 'notag.f'
        include 'nflav.f'
        include 'reset.f'
        include 'TRmaxindex.f'
        include 'nplot.f'
        include 'epinv.f'
        include 'epinv2.f'
        include 'anomcoup.f'
        include 'hbbparams.f'
        include 'ipsgen.f'
        include 'taucut.f'
        ! end omp copyin

        include 'mxdim.f'

        integer :: it,j,l,m
        integer(int64) :: k
        real(dp) :: xxn
        real(dp) :: x(mxdim)

        integer :: ia(1:ndim)
        real(dp) :: delx(1:ndim), wgt
        real(dp) :: f, f2, fsum, f2sum

        real(dp) :: f2binned(ndmx,ndim)
        real(dp) :: f2binsum(ndim)
        real(dp) :: mi(ndmx,ndim), dx, gridnew(ndmx)

        real(dp) :: alpha
        real(dp) :: sigsq, chisq
        integer(int64) :: adjustedCalls

        logical :: bin
        common/bin/bin
        integer :: ierr

        ! cpu time
        real :: time_start, time_step, time_end
        real(dp) :: timeacc

        ! cpu time total
        real(dp) :: cputime_total


        ! wall time
        integer :: count_rate
        integer :: count_start, count_step, count_end
        real(dp) :: countacc


        if (stage == 0) then ! reset grid and results
            info%si = 0._dp
            info%swgt = 0._dp
            info%schi = 0._dp
            info%lastIter = 0
            alpha = info%damp(1)
            if (.not. allocated(info%grid)) then
                allocate(info%grid(1:ndmx, 1:ndim))
            endif
            call uniformgrid(info%grid)
        else if (stage == 1) then  ! only reset results, keep grid
            info%si = 0._dp
            info%swgt = 0._dp
            info%schi = 0._dp
            info%lastIter = 0
            alpha = info%damp(2)
            if (.not. allocated(info%grid)) then
                write (*,*) "WARNING: grid not allocated in integration stage 1"
                allocate(info%grid(1:ndmx, 1:ndim))
                call uniformgrid(info%grid)
            endif
        else ! keep grid and results (just resume)
            alpha = info%damp(2)
            if (.not. allocated(info%grid)) then
                write (*,*) info%sig(), info%sd()
                error stop "integration stage > 1 but grid not allocated"
            endif
        endif

        ! make callsPerIt evenly divisible
        if (world_size > 1) then
            adjustedCalls = nint(real(info%callsPerIt,dp) / real(world_size,dp), int64)
            if ((adjustedCalls*world_size) /= info%callsPerIt) then
                info%callsPerIt = adjustedCalls*world_size
                if (rank == 0) then
                    write (*,*) " Adjusted number of calls per iteration to ", adjustedCalls
                    write (*,*) " to account for ", world_size, " mpi processes"
                endif
            endif
        else
            adjustedCalls = info%callsPerIt
        endif

        if (rank == 0) then
            write (*,*) ""
            write (*,'(A)') " Vegas integration parameters:"
            write (*,'(A,I2,A,I12,A,I2,A,F3.1)') " ndim = ", ndim, &
                "  ncall = ", info%callsPerIt, "  iter = ", niter, &
                "  alpha = ", alpha
            write (*,*) ""
        endif

        call cpu_time(time_start)
        call system_clock(count_start, count_rate)

        do it=1,niter
            ! reset accumulators for each iteration
            fsum = 0._dp
            f2sum = 0._dp
            f2binned(:,:) = 0._dp
            x(:) = 0._dp

            call cpu_time(time_step)
            call system_clock(count_step)
            timeacc = 0._dp
            countacc = 0

!$omp  parallel do &
!$omp& schedule(guided) &
!$omp& default(private) &
!$omp& shared(ndim,info,bin,adjustedCalls) &
!$omp& shared(rank,world_size,ndmx) &
!$omp& copyin(/xmin/,/taumin/,/cutoff/,/jetcuts/,/breit/,/zerowidth/) &
!$omp& copyin(/srdiags/,/vsymfact/,/qcdcouple/,/ewcouple/,/masses/) &
!$omp& copyin(/interference/,/facscale/,/mcfmscale/,/stopscales1/) &
!$omp& copyin(/stopscales2/, /mcfmtaucut/) &
!$omp& copyin(/ipsgen/,/pext/,/ColC/,/ptildes/) &
!$omp& copyin(/bitflags/,/flags/,/lastphot/) &
!$omp& copyin(/lhcb1/,/lhcb2/,/lhcb3/,/lhcb4/,/lhcb5/) &
!$omp& copyin(/QCDb0/,/nodecay/,/swapxz/,/heavyflav/) &
!$omp& copyin(/notag/,/nflav/,/reset/,/pvmaxindex/) &
!$omp& copyin(/epinv/,/epinv2/,/plotindex/,/hbbparams/) &
!$omp& copyin(/anomcoup1/,/anomcoup2/) &
!$omp& reduction(+:fsum) reduction(+:f2sum) reduction(+:f2binned)
            do k=1,adjustedCalls
                if (info%useSobol) then
!$omp critical(sobol)
                    do l=1,ndim+2
                        x(l) = info%sstate(l)%next_strided()
                    enddo
!$omp end critical(sobol)
                else
                    do l=1,ndim + 2
                        x(l) = cxx11_random_number()
                    enddo

                endif

                wgt = 1._dp / info%callsPerIt

                do j=1,ndim
                    xxn = x(j) * ndmx + 1._dp
                    ia(j) = max(min(int(xxn),ndmx),1)

                    if (ia(j) == 1) then
                        delx(j) = info%grid(1,j)
                        x(j) = (xxn-1._dp) * delx(j)
                    else
                        delx(j) = info%grid(ia(j),j) - info%grid(ia(j)-1,j)
                        x(j) = info%grid(ia(j)-1,j) + (xxn-ia(j))*delx(j)
                    endif

                    wgt = wgt * delx(j) * ndmx
                enddo

                f = wgt * fxn(x,wgt)
                f2 = f**2


                fsum = fsum + f
                f2sum = f2sum + f2

                do j=1,ndim
                    f2binned(ia(j), j) = f2binned(ia(j), j) + f2
                enddo

            enddo ! calls
!$omp end parallel do

            !!!
            ! OMP thread local reduction
            ! threadStorage results are accumulated into masterStorage
            !!!
            if (bin) then
!$omp parallel
!$omp critical
            call accumulateThreadStorage()
!$omp end critical
!$omp end parallel
            endif


#ifdef HAVE_MPI
            !!!
            ! MPI synchronization here
            ! synchronize fsum, f2sum, f2binned
            ! as well as masterStorages from processes
            !!!
            if (bin .and. world_size > 1) then
                call mpi_reduce_masterStorage()
            endif

            if (world_size > 1) then
                call mpi_allreduce(mpi_in_place, fsum, 1, mpi_double_precision, &
                    mpi_sum, mpi_comm_world, ierr)
                call mpi_allreduce(mpi_in_place, f2sum, 1, mpi_double_precision, &
                    mpi_sum, mpi_comm_world, ierr)
                call mpi_allreduce(mpi_in_place, f2binned, ndim*ndmx, &
                    mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
            endif
#endif

            !!!
            ! accumulation from masterStorage into iterationStorage
            !!!
            if (bin) then
                call accumulateIteration(info%callsPerIt)
            endif

            do j=1,ndim
                f2binsum(j) = 0._dp
                do k=1,ndmx
                    if (k == 1) then
                        dx = info%grid(k,j)
                    else
                        dx = info%grid(k,j) - info%grid(k-1,j)
                    endif
                    f2binsum(j) = f2binsum(j) + sqrt(f2binned(k,j)) * dx
                enddo

                do k=1,ndmx
                    !mi(k,j) = sqrt(f2binned(k,j)) * dx
                    if (f2binned(k,j) /= 0._dp) then
                        mi(k,j) = ((sqrt(f2binned(k,j))/f2binsum(j) - 1._dp) / &
                            log( sqrt(f2binned(k,j))/f2binsum(j) ) )**alpha
                    else
                        mi(k,j) = 0._dp
                    endif
                enddo
            enddo

            ! normalize and rebin
            do j=1,ndim
                mi(:,j) = mi(:,j) / sum(mi(:,j))
                call rebin(mi(:,j), info%grid(:,j))
            enddo

            sigsq = (info%callsperit * f2sum - fsum**2) / (info%callsperit - 1)

            info%lastIter = info%lastIter + 1
            if (fsum /= 0._dp) then
                info%si = info%si +  fsum / sigsq
                info%swgt = info%swgt + 1._dp / sigsq
                info%schi = info%schi + fsum**2 / sigsq
            endif

            if (stage > 0) then
                info%sobolSkip = info%sobolSkip + info%callsPerIt * niter
            endif

            if (info%lastIter == 1) then
                chisq = 0._dp
            else
                chisq = (info%schi - info%si* info%si / info%swgt) / (info%lastIter - 1._dp)
            endif


#ifdef HAVE_MPI
            !!!
            ! MPI synchronization
            ! push back grid, etc. to avoid small differences in numerical calculation
            ! in different processes
            ! these could lead to different branches in the outer code
            !!!
            if (world_size > 1) then
                call mpi_bcast(info%grid(:,:), product(shape(info%grid)), &
                    mpi_double_precision, 0, mpi_comm_world, ierr)
                call mpi_bcast(info%si, 1, mpi_double_precision, 0, mpi_comm_world, ierr)
                call mpi_bcast(info%swgt, 1, mpi_double_precision, 0, mpi_comm_world, ierr)
                call mpi_bcast(info%schi, 1, mpi_double_precision, 0, mpi_comm_world, ierr)
            endif
            !!!
            ! we do not synchronize histograms etc. here, since they just need
            ! to be correct on the root rank and master thread for outputting
            ! to files and to the screen
            !!!
#endif

            call cpu_time(time_end)
            timeacc = time_end - time_step
#ifdef HAVE_MPI
            if (world_size > 1) then
                call mpi_allreduce(mpi_in_place, timeacc, 1, mpi_double_precision, &
                    mpi_sum, mpi_comm_world, ierr)
            endif
#endif

            call system_clock(count_end)
            countacc = real(count_end - count_step,dp)/real(count_rate,dp)

            if (rank == 0) then
                write (*, '(A,I3,A)') '************** Integration by Vegas (iteration ', &
                            info%lastIter, ') ***************'
                write (*,'(A,A65,A)') '*','','*'
                write (*,'(A,G15.8,A,G15.8,A)') '*  integral  = ', fsum, &
                            '   accum. integral = ', info%si / info%swgt, '*'
                write (*,'(A,G15.8,A,G15.8,A)') '*  std. dev =  ', sqrt(sigsq), &
                            '   accum. std. dev = ', 1._dp/sqrt(info%swgt), '*'
                write (*,'(A,G15.6,A36,A)') '*   max. wt. = ', wtmax, '', '*'
                write (*,'(A,A65,A)') '*','','*'
                write (*,'(A,G15.8,A)') '*  CPU time used: ', timeacc,' seconds *'
                write (*,'(A,G15.8,A)') '*  Wall time used: ', countacc,' seconds*'
                write (*,'(A,A65,A)') '*','','*'
                write (*,'(A,G12.4,A)') '***************   chi**2/iteration = ', chisq, &
                            '   ***************'
                write (*,*) ""
            endif

            if (fsum == 0._dp) exit

        enddo ! iteration

        call cpu_time(time_end)
        timeacc = time_end - time_start
#ifdef HAVE_MPI
        if (world_size > 1) then
            call mpi_allreduce(mpi_in_place, timeacc, 1, mpi_double_precision, &
                mpi_sum, mpi_comm_world, ierr)
        endif
#endif

        countacc = real(count_end - count_start,dp)/real(count_rate,dp)

        if (rank == 0) then
            write (*,*) ""
            write (*,'(A,I2,A,G15.8)') "CPU time for last ", niter, " iterations: ", timeacc
            write (*,'(A,I2,A,G15.8)') "Wall time for last ", niter, " iterations: ", countacc
            write (*,*) ""
        endif

        cputime_total = get_cputime()
        if (rank == 0) then
            write (*,*) ""
            write (*,'(A,G15.8)') "CPU time total: ", cputime_total
            write (*,'(A,G15.8)') "Wall time total: ", get_walltime()
            write (*,*) ""
        endif

    end subroutine


end module
