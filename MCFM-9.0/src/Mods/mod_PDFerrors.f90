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
 
module PDFerrors
    use types
    implicit none

    public

    logical, public, save :: doPDFerrors

    ! this is set to .true. when there are different as(mZ) values
    ! among the PDF sets and members
    logical, public, save :: doPDFAlphas

    ! support for multiple pdf sets
    character(len=256), allocatable :: PDFnames(:)
    integer, allocatable :: PDFmembers(:)
    ! this is the size of PDFnames and PDFmembers
    integer :: numPDFsets

    integer, public, save :: currentPDF
!$omp threadprivate(currentPDF)


    ! this is just the size of pdfreweight
    integer, save :: maxPDFsets

    real(dp), public, save, allocatable :: pdfreweight(:)
!$omp threadprivate(pdfreweight)

end module
