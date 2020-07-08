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
 
      subroutine pdfwrap
          use MCFMStorage
          use LHAPDF
          use PDFerrors
      implicit none
      include 'masses.f'
      include 'lhapdf.f'
      include 'nlooprun.f'
      include 'pdlabel.f'
      include 'couple.f'
      include 'mpicommon.f'
      integer :: i,j,iorder

      include 'first.f'

      if(first .eqv. .false.) then
          return
      else
          first = .false.
      endif

      if (rank == 0) then
!$omp master
      write(6,*)
      write(6,*) '*******************************************'
      write(6,*) '*     MCFM is calling LHAPDF              *'
      write(6,*) '*                                         *'
      do j=1,numPDFsets
         write (6,98) "PDFname", PDFnames(j)(1:49)
         write (6,99) "PDFmember", PDFmembers(j)
         write (6,*) ""
      enddo
      write(6,*) '*******************************************'
      write(6,*)
!$omp end master
      endif

      if (doPDFerrors) then
          maxPDFsets = 0
          do i=1,numPDFsets
              maxPDFsets = maxPDFsets + lhapdf_number(trim(PDFnames(i))) 
          enddo
          maxPDFsets = maxPDFsets -1

          call initall(PDFnames, PDFmembers)

!$omp parallel
        if (.not. allocated(pdfreweight)) then
            allocate(pdfreweight(maxPDFsets))
        endif
!$omp end parallel
        if (rank == 0) then
            write(6,*)
            write(6,*) '****************************************'        
            write(6,*) '*        Calculating errors using      *'
            write(6,97) maxPDFsets
            write(6,*) '****************************************'
        endif

      else
          maxPDFsets = size(PDFmembers) - 1
          call initcentral(PDFnames, PDFmembers)
!$omp parallel
          allocate(pdfreweight(maxPDFsets))
!$omp end parallel
      endif

!$omp parallel
      ! this is used in getalphas and getorderqcd below
      currentPDF = 0
!$omp end parallel


      amz = getalphas(zmass)
      nlooprun = getorderqcd()

      pdlabel = ''
      i=1
      do j=1,numPDFsets
        pdlabel(i:i+len(trim(PDFnames(j)))+1) = trim(PDFnames(j))//"_"
        i = i + len(trim(PDFnames(j))) + 1
      enddo
      pdlabel = pdlabel(1:len(trim(pdlabel))-1)

      ! Determine if there are different values of alphas.
      ! There is no reliable way to get this information from LHAPDF:
      ! some sets have ErrorType "XX+as", but there are sets that
      ! only have members with different alphas, and their ErrorType
      ! is just "replicas".

      doPDFAlphas = .false.
      do j=1,maxPDFsets
          currentPDF = j
          if (abs(getalphas(zmass) - amz) > 1d-5) then
              doPDFAlphas = .true.
          endif
      enddo
      currentPDF = 0

      if (rank == 0 .and. doPDFAlphas) then
          write(6,*) '****************************************'        
          write(6,*) '* Alphas variation included in PDF set *'
          write(6,*) '****************************************'
      endif

      return
 
   97 format(' *        ',i4,' sets of error PDFs       *')
   98 format(' *   ',a7,' ',a29,' *')
   99 format(' *  ',a10,i3,'                          *')

      end
 

