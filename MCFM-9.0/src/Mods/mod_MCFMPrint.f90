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
 
module MCFMPrint
    use types
    implicit none

    contains

    function formatCross(cross, error)
        implicit none
        real(dp), intent(in) :: cross
        real(dp), intent(in), optional :: error
        character(len=:), allocatable :: formatCross
        character(len=255) :: tmp

        if (cross < 1d3) then
            if (present(error)) then
                write (tmp,'(G14.6,A,G13.5,A)') cross, " ± ", error, " fb"
            else
                write (tmp,'(G14.6,A)') cross, " fb"
            endif
        else if (cross < 1d6) then
            if (present(error)) then
                write (tmp,'(G14.6,A,G13.5,A)') cross/1d3, " ± ", error/1d3, " pb"
            else
                write (tmp,'(G14.6,A)') cross/1d3, " pb"
            endif
        else
            if (present(error)) then
                write (tmp,'(G14.6,A,G13.5,A)') cross/1d6, " ± ", error/1d6, " nb"
            else
                write (tmp,'(G14.6,A)') cross/1d6, " nb"
            endif
        endif

        formatCross = trim(tmp)
    end function

    function spacereplace(string)
        implicit none
        character(len=*), intent(in) :: string
        character(len=len(string)) :: spacereplace
        integer :: j

        do j=1,len(string)
            if (string(j:j) == ' ') then
                spacereplace(j:j) = '_'
            else
                spacereplace(j:j) = string(j:j)
            endif
        enddo

    end function

    function chisqMax()
        use MCFMStorage
        implicit none

        integer :: i,j
        real(dp) :: chisqmax

        chisqMax = maxval( iterationStorage(:,:)%vinfo%chisq() , iterationStorage(:,:)%used )

    end function

#ifdef HAVE_LHAPDF
    subroutine pdfUncertaintyHistograms(nset, histCentral, histSymm, histPlus, histMinus)
        use Superhisto
        use MCFMStorage
        use LHAPDF
        use PDFerrors
        implicit none
        integer, intent(in) :: nset
        type(HistogramStorage), intent(out) :: histCentral, histSymm, histPlus, histMinus
        integer :: j,k,l
        real(dp), allocatable :: values(:)
        type(pdfuncertainty) :: errors
        integer :: numpdfmembers
        integer :: nstart
        real(dp) :: errPdferr

        write (*,*) "Computing PDF uncertainties for histograms"

        histCentral = finalSum%histCentral
        histSymm = finalSum%histCentral
        call histSymm%reset()
        histPlus = finalSum%histCentral
        call histPlus%reset()
        histMinus = finalSum%histCentral
        call histMinus%reset()


        numpdfmembers = lhapdf_number(trim(PDFnames(nset)))
        allocate(values(numpdfmembers))

        nstart = 0
        do j=1,nset-1
            nstart = nstart + lhapdf_number(trim(PDFnames(j)))
        enddo

        do k=1,size(histCentral%histos)
            if (histCentral%histos(k)%initialized()) then
                do l=0,histCentral%histos(k)%nbins+1
                    do j=1,numpdfmembers
                        if (j==1 .and. nset == 1) then
                            values(1) = finalSum%histCentral%histos(k)%xx(l)
                        else
                            values(j) = finalSum%histCentral%histos(k)%xx(l) + &
                                finalSum%histPDFerrors(nstart + j - 1)%histos(k)%xx(l)
                        endif
                    enddo


                    errors = computeUncertainty(values, nstart)
                    
                    histCentral%histos(k)%xx(l) = errors%central
                    histSymm%histos(k)%xx(l) = errors%errsymm
                    histPlus%histos(k)%xx(l) = errors%errplus
                    histMinus%histos(k)%xx(l) = errors%errminus

                    ! save error on pdf error in histSymm
                    errPdferr = 0._dp
                    do j=nstart+1,nstart+numpdfmembers-1
                        errPdferr = max(errPdferr,finalSum%histPDFerrors(j)%histos(k)%xxsq(l))
                    enddo
                    histSymm%histos(k)%xxsq(l) = errPdferr

                enddo
            endif
        enddo

        write (*,*) "Done computing PDF uncertainties for histograms"

    end subroutine
#endif

    subroutine writeAllHistograms()
        use parseinput
        use Superhisto
        use PDFerrors
        implicit none

        logical :: writetxt, writetop
        call cfg_get(cfg, "histogram%writetxt", writetxt)
        call cfg_get(cfg, "histogram%writetop", writetop)

        if (writetxt) then
            call writeHistograms(shwrite,"txt")
            if (doPDFerrors) then
                call writeHistogramsPDFerrors(shwritepdf,"txt")
            endif
        endif
        if (writetop) then
            call writeHistograms(shwritetop,"top")
        endif

    end subroutine

    subroutine writeHistogramsPDFerrors(writefun,ext)
        use Superhisto
        use MCFMStorage
        use PDFerrors
        implicit none

        character(len=3), intent(in) :: ext

        character(len=255) :: runname
        common/runname/runname
        integer :: ierr, k,l,j
        character(len=255) :: imsg
        character(len=1024) :: filename

        interface
            subroutine writefun(histCentralX,histSymmX,histPlusX,histMinusX,iounit)
                use Superhisto
                implicit none
                type(sh_histogram), intent(in) :: histCentralX, histSymmX, histPlusX, histMinusX
                integer, intent(in) :: iounit
            end subroutine
        end interface

        type(HistogramStorage) :: histCentral, histSymm, histPlus, histMinus

#ifdef HAVE_LHAPDF
        if (doPDFerrors) then
            do j=1,numPDFsets
                call pdfUncertaintyHistograms(j, histCentral, histSymm, histPlus, histMinus)
                ! calculate actual pdf error histograms
                do k=1,size(histCentral%histos)
                    if (histCentral%histos(k)%initialized()) then
                        write (filename,'(A)') trim(runname)//"_"//spacereplace(trim(histCentral%histos(k)%title)) &
                                                //"_pdferr_"//trim(pdfnames(j))//"."//ext
                        open(unit=11, file=filename, status='replace', form='formatted', iostat=ierr, iomsg=imsg)
                        if (ierr == 0) then
                            call writefun(histCentral%histos(k),histSymm%histos(k), &
                                    histPlus%histos(k),histMinus%histos(k),11)
                            close(unit=11)
                        else
                            write (*,*) "Problem writing histogram file "//filename
                            write (*,*) trim(imsg)
                            write (*,*) "Error code = ", ierr
                        endif
                    endif
                enddo
            enddo
        endif
#endif

    end subroutine

    subroutine writeHistograms(writefun,ext)
        use Superhisto
        use MCFMStorage
        use Scalevar
        use PDFerrors
        use parseinput
        use SCET
        use MCFMTaufit, only : autofit
        implicit none

        character(len=3), intent(in) :: ext

        character(len=255) :: runname
        common/runname/runname
        integer :: ierr, k,l
        character(len=255) :: imsg
        character(len=1024) :: filename
        character(len=255) :: centralname

        interface
            subroutine writefun(histo,iounit)
                use Superhisto
                implicit none
                type(sh_histogram), intent(in) :: histo
                integer, intent(in) :: iounit
            end subroutine
        end interface

        if (numPDFsets > 1) then
            centralname = '' ! "_"//trim(PDFnames(1))
        else
            centralname = ''
        endif

        do k=1,size(finalSum%histCentral%histos)
            associate (histo => finalSum%histCentral%histos(k))
                if (histo%initialized()) then
                    write (filename,'(A)') trim(runname)//"_"// &
                        spacereplace(trim(histo%title))//trim(centralname)//"."//ext
                    open(unit=11, file=filename, status='replace', form='formatted', iostat=ierr, iomsg=imsg)
                    if (ierr == 0) then
                        call writefun(histo,11)
                        close(unit=11)
                    else
                        write (*,*) "Problem writing histogram file "//filename
                        write (*,*) trim(imsg)
                        write (*,*) "Error code = ", ierr
                    endif
                endif
            end associate
        enddo

        if ((doPDFerrors .eqv. .false.) .and. numPDFsets > 1) then
            do l=2,numPDFsets
                do k=1,size(finalSum%histPDFerrors(l-1)%histos)
                    associate (histo => finalSum%histPDFerrors(l-1)%histos(k))
                        if (histo%initialized()) then
                            write (filename,'(A)') trim(runname)//"_"// &
                                spacereplace(trim(histo%title))// & 
                                "_"//trim(PDFnames(l))//"."//ext
                            open(unit=11, file=filename, status='replace', form='formatted', iostat=ierr, iomsg=imsg)
                            if (ierr == 0) then
                                call writefun(histo,11)
                                close(unit=11)
                            else
                                write (*,*) "Problem writing histogram file "//filename
                                write (*,*) trim(imsg)
                                write (*,*) "Error code = ", ierr
                            endif
                        endif
                    end associate
                enddo
            enddo
        endif

        if (doScalevar) then
            do k=1,maxScalevar
                do l=1,size(finalSum%histScalevar(k)%histos)
                    associate (histo => finalSum%histScalevar(k)%histos(l))
                        if (histo%initialized()) then
                            write (filename,'(A,I1,A)') trim(runname)//"_"// &
                                spacereplace(trim(histo%title))//"_scalevar_", k, ".txt"
                            open(unit=11, file=filename, status='replace', form='formatted', iostat=ierr, iomsg=imsg)
                            if (ierr == 0) then
                                call writefun(histo,11)
                                close(unit=11)
                            else
                                write (*,*) "Problem writing histogram file "//filename
                                write (*,*) trim(imsg)
                                write (*,*) "Error code = ", ierr
                            endif
                        endif
                    end associate
                enddo
            enddo
        endif

        if (size(tcutarray) > 1) then
            taucutfit: block
                type(HistogramStorage) :: histFitted
                real(dp), allocatable :: ydat(:), weights(:)
                real(dp), allocatable :: relerrs(:)
                real(dp) :: fitcross, fiterr, reserr
                integer :: i, info

                histFitted = finalSum%histCentral
                allocate(ydat(size(tcutarray)))
                allocate(weights(size(tcutarray)))
                allocate(relerrs(size(tcutarray)))

                ! fill histFitted with fits
                do l=1,size(histFitted%histos)
                    if (histFitted%histos(l)%initialized()) then
                        allocate(histFitted%histos(l)%extras(2, 0:histFitted%histos(l)%nbins+1))

                        do k=0,histFitted%histos(l)%nbins+1
                            do i=1,size(tcutarray)
                                ydat(i) = finalSum%histTaucut(i)%histos(l)%xx(k)
                                weights(i) = 1._dp / finalSum%histTaucut(i)%histos(l)%xxsq(k)
                                relerrs(i) = abs( finalSum%histTaucut(i)%histos(l)%xxsq(k) / &
                                    finalSum%histTaucut(i)%histos(l)%xx(k))
                            enddo

                            call autofit(ydat, weights, fitcross, fiterr, reserr, info)
                            histFitted%histos(l)%xx(k) = fitcross
                            histFitted%histos(l)%xxsq(k) = fiterr
                            histFitted%histos(l)%extras(1,k) = maxval(relerrs)
                            histFitted%histos(l)%extras(2,k) = reserr
                        enddo
                    endif
                enddo

                ! and write out
                do k=1,size(histFitted%histos)
                    associate (histo => histFitted%histos(k))
                        if (histo%initialized()) then
                            write (filename,'(A)') trim(runname)//"_"// &
                                spacereplace(trim(histo%title))//"_taucutfit"//"."//ext
                            open(unit=11, file=filename, status='replace', form='formatted', iostat=ierr, iomsg=imsg)
                            if (ierr == 0) then
                                call writefun(histo,11)
                                close(unit=11)
                            else
                                write (*,*) "Problem writing histogram file "//filename
                                write (*,*) trim(imsg)
                                write (*,*) "Error code = ", ierr
                            endif
                        endif
                    end associate
                enddo
            end block taucutfit
        endif

        do k=1,size(tcutarray)
            do l=1,size(finalSum%histTaucut(k)%histos)
                associate (thist => finalSum%histTaucut(k)%histos(l))
                    if (thist%initialized()) then
                        write (filename,'(A,F0.6,A)') trim(runname)//"_"// &
                            spacereplace(trim(thist%title))//"_taucut_", tcutarray(k), ".txt"
                        open(unit=11, file=filename, status='replace', form='formatted', iostat=ierr, iomsg=imsg)
                        if (ierr == 0) then
                            call writefun(thist,11)
                            close(unit=11)
                        else
                            write (*,*) "Problem writing histogram file "//filename
                            write (*,*) trim(imsg)
                            write (*,*) "Error code = ", ierr
                        endif
                    endif
                end associate
            enddo
        enddo

    end subroutine

#ifdef HAVE_LHAPDF
    ! assumes finalizeStorage has been filled
    subroutine printPDFuncertainties()
        use MCFMStorage
        use PDFerrors
        use LHAPDF
        implicit none
        real(dp), allocatable :: values(:)
        type(pdfuncertainty) :: errors
        integer :: j,k,numpdfmembers
        real(dp) :: errPdferr
        integer :: l

        integer :: accum

        if (doPDFerrors) then
            accum = 0
            do j=1,numPDFsets
                numpdfmembers = lhapdf_number(trim(PDFnames(j)))
                if (.not. allocated(values)) then
                    allocate(values(numpdfmembers))
                else
                    deallocate(values)
                    allocate(values(numpdfmembers))
                endif

                errPdferr = 0._dp
                do l=accum+1,accum+numpdfmembers-1
                    errPdferr = max(errPdferr,finalSum%histPDFerrors(l)%histos(1)%xxsq(1))
                enddo

                do k=1,numpdfmembers
                    if (j==1 .and. k==1 ) then
                        values(1) = finalSum%histCentral%histos(1)%xx(1)
                    else
                        values(k) = finalSum%histCentral%histos(1)%xx(1) + &
                            finalSum%histPDFerrors(accum)%histos(1)%xx(1)
                    endif
                    accum = accum + 1
                enddo

                errors = computeUncertainty(values, accum-2)

                write(*,*) '********* PDF uncertainty analysis *********'
                write(*,'(A,A)')  " for PDF set ", trim(pdfnames(j))
                write(*,*) ""
                write(*,*) 'Central value',errors%central
                write(*,*) ""
                write(*,*) '         Absolute PDF uncertainties'
                write(*,'(A,G11.4)') '   Symmetric +/-',errors%errsymm
                write(*,'(A,G11.4)') '   +ve direction',errors%errplus
                write(*,'(A,G11.4)') '   -ve direction',errors%errminus
                write(*,*) ""
                write(*,*) "Estimated numerical uncertainty for PDF uncertainties"
                write(*,'(A,G10.2,A)') " for symmetric: ", abs(errPdferr/errors%errsymm)*100._dp, " percent"
                write(*,'(A,G10.2,A)') " for +ve direction: ", abs(errPdferr/errors%errplus)*100._dp, " percent"
                write(*,'(A,G10.2,A)') " for -ve direction: ", abs(errPdferr/errors%errminus)*100._dp, " percent"
                write(*,*) ""
                write(*,*) '         Relative PDF uncertainties'
                write(*,'(A,G10.2,A)') '   Symmetric +/- ',errors%errsymm/errors%central*100._dp, " percent"
                write(*,'(A,G10.2,A)') '   +ve direction ',errors%errplus/errors%central*100._dp, " percent"
                write(*,'(A,G10.2,A)') '   -ve direction ',errors%errminus/errors%central*100._dp, " percent"
                write(*,*) ""
                write(*,*) '********************************************'
            enddo
        endif

    end subroutine
#endif

    subroutine printScaleUncertainties()
        use MCFMStorage
        use Scalevar
        implicit none

        if (doScalevar) then
            write (*,*) ""
            write (*,*) "Scale variation results: "
            write (*,*) "Central value: ", finalSum%histCentral%histos(1)%xx(1), "+/-", &
                finalSum%histCentral%histos(1)%xxsq(1)
            write (*,*) "Differences from central value:"
            if (maxScalevar >= 2) then
                write (*,*) " muR*2, muF*2: ", finalSum%histScalevar(1)%histos(1)%xx(1), "+/-", &
                    finalSum%histScalevar(1)%histos(1)%xxsq(1)
                write (*,*) " muR/2, muF/2: ", finalSum%histScalevar(2)%histos(1)%xx(1), "+/-", &
                    finalSum%histScalevar(2)%histos(1)%xxsq(1)
            endif
            if (maxScalevar >= 6) then
                write (*,*) " muR*2,   muF: ", finalSum%histScalevar(3)%histos(1)%xx(1), "+/-", &
                    finalSum%histScalevar(3)%histos(1)%xxsq(1)
                write (*,*) " muR/2,   muF: ", finalSum%histScalevar(4)%histos(1)%xx(1), "+/-", &
                    finalSum%histScalevar(4)%histos(1)%xxsq(1)
                write (*,*) " muR,   muF*2: ", finalSum%histScalevar(5)%histos(1)%xx(1), "+/-", &
                    finalSum%histScalevar(5)%histos(1)%xxsq(1)
                write (*,*) " muR,   muF/2: ", finalSum%histScalevar(6)%histos(1)%xx(1), "+/-", &
                    finalSum%histScalevar(6)%histos(1)%xxsq(1)
            endif
            write (*,*) ""
        endif
    end subroutine

    subroutine printTaucuts()
        use MCFMStorage
        use SCET
        use MCFMTaufit
        implicit none
        include 'taucut.f'
        include 'kpart.f'

        integer :: j
        if (smallestTaucut < taucut) then
            write (*,*) "Take the values below the nominal taucut value with a grain of salt"
            write (*,*) "since the statistics will be significantly lower and errors might be"
            write (*,*) "severely underestimated."
        endif
        write (*,*) ""
        write (*,*) " Inclusive cross section taucut dependence: "
        do j=1,size(tcutarray)
            write (*,'(A,F0.6,A,F0.6,A,A)') " sigma(tcut=", taucut ,") - sigma(tcut=", tcutarray(j), ") = ", &
                formatCross(finalSum%histTaucut(j)%histos(1)%xx(1), finalSum%histTaucut(j)%histos(1)%xxsq(1))
        enddo
        write (*,*) ""

        ! tau cut fit
        taucutfit: block
            real(dp), allocatable :: ydat(:)
            real(dp), allocatable :: weights(:)
            real(dp) :: fitcross, fiterr, reserr
            integer :: i, info

            allocate(ydat(size(tcutarray)))
            allocate(weights(size(tcutarray)))

            do i=1,size(tcutarray)
                ydat(i) = finalSum%histTaucut(i)%histos(1)%xx(1)
                weights(i) = 1._dp / finalSum%histTaucut(i)%histos(1)%xxsq(1)
            enddo
            
            call autofit(ydat, weights, fitcross, fiterr, reserr, info)

            write (*,*) "Fitted correction: ", formatCross(fitcross, fiterr)
            write (*,*) "Reduced chisquare: ", reserr
            write (*,*) "The reduced chisquare should be small or close to one"
            write (*,*) "for a good fit. Results can be trusted if the taucut "
            write (*,*) "dependence uncertainties are reliable (small)."
            write (*,*) ""

        end block taucutfit
    end subroutine

    ! to print central values of all sets
#ifdef HAVE_LHAPDF
    subroutine printallcross()
        use PDFerrors
        use LHAPDF
        use MCFMStorage, only : finalSum
        implicit none
        real(dp) :: centralCross, centralError, cross, error
        integer, allocatable :: centralBins(:)
        integer :: j, accum

        allocate(centralBins(numPDFsets))

            ! the first set is our MC central value, saved separately
        accum = 1
        do j=2,numPDFsets
            if (doPDFerrors) then
                accum = accum + lhapdf_number(trim(PDFnames(j-1)))
            else
                accum = accum + 1
            endif
            
            centralBins(j) = accum
        enddo

        centralBins(2:) = centralBins(2:) - 1

        if (numPDFsets > 1) then
            write (*,'(A)') "=== Printing central cross section values for all PDF sets ==="
        endif

        write (*,'(A,A,A,I3,A)') "=== Result for PDF set ", trim(PDFnames(1)), " member ", &
            PDFmembers(1), " ==="
        centralCross = finalSum%histCentral%histos(1)%xx(1)
        centralError = finalSum%histCentral%histos(1)%xxsq(1)
        call printcross(centralCross,centralError,chisqmax())

        do j=2,numPDFsets
            cross = finalSum%histPDFerrors(centralBins(j))%histos(1)%xx(1)
            error = finalSum%histPDFerrors(centralBins(j))%histos(1)%xxsq(1)

            write (*,*) ""
            write (6,'(A,A,A,A)') " Difference from ", trim(PDFnames(j)), &
                " to first PDF set is ", formatCross(cross,error)
        enddo

    end subroutine
#endif

    subroutine printcross(cross,error,chisqmax)
        implicit none
        real(dp), intent(in) :: cross, error, chisqmax
        include 'kprocess.f'
        real(dp) :: rescale

        write (6,'(A,A)') "Value of integral is ", formatCross(cross, error)

        write (*,'(A,G10.3)') " Maximum chisq/it over all contributions is ", chisqmax

! JC: removed August 2019, probably no longer interesting
        !--- for gg->H+X processes, also write out the cross section
        !---  normalized by sigma(gg->H, finite mt)/sigma(gg->H, mt-> infinity)
!        if ( (kcase == kggfus0) .or. (kcase == kggfus1) &
!            .or.(kcase == kggfus2) .or. (kcase == kggfus3) &
!            .or.(kcase == kHWWjet) .or. (kcase == kHWW2jt) &
!            .or.(kcase == kHWW3jt) .or. (kcase == kHWW_4l) &
!            .or.(kcase == kHWW2lq) .or. (kcase == kHWWdkW) &
!            .or.(kcase == kHZZjet) .or. (kcase == kHZZ2jt) &
!            .or.(kcase == kHZZ3jt) .or. (kcase == kHZZpjt) &
!            .or.(kcase == kHZZ_jj) .or. (kcase == kHZZ_4l) &
!            .or.(kcase == kHZZqgI) ) then
!
!            call finitemtcorr(rescale)
!            write(6,*)
!            write(6,*) 'Cross section normalized by the ratio'
!            write(6,*) 'sigma(gg->H, finite mt)/sigma(gg->H, mt-> infinity)'
!            write(6,*) '(i.e. exact for gg->H process, but &
!                         &approx. for gg->H+n jets, n=1,2,3)'
!            write(6,*)
!            write(6,'(A,A)') ' Rescaled integral is', formatCross(cross*rescale, error*rescale)
!            write(6,'(a25,f9.5,a2)') '   (Rescaling factor is ',rescale,')'  
!        endif

    end subroutine

end module
