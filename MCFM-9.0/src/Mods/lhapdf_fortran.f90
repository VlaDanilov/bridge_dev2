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
 
module LHAPDF
      use iso_c_binding
      use types
      implicit none

      type pdf
          private
          character(:), allocatable :: setname
          type(c_ptr) :: ptr
          contains
          procedure :: xfxq2 => pdf_xfxq2
          procedure :: nmem => lhapdf_number_pdf
          procedure :: alphas => pdf_alphas
          procedure :: orderqcd => pdf_orderqcd
      end type

      type, bind(c) :: pdfuncertainty
          real(c_double) :: central
          real(c_double) :: errplus
          real(c_double) :: errminus
          real(c_double) :: errsymm
          real(c_double) :: scale
          real(c_double) :: errplus_pdf
          real(c_double) :: errminus_pdf
          real(c_double) :: errsymm_pdf
          real(c_double) :: err_par
      end type

      public :: pdf
      public :: pdf_xfxq2
      public :: lhapdf_loadMember
      public :: lhapdf_info
      public :: lhapdf_number
      public :: pdfuncertainty

      ! these routines work on the global pdfs
      public :: fdist_lhapdf
      public :: initcentral, initall
      public :: getalphas, getorderqcd
      public :: computeUncertainty

      interface lhapdf_number
          module procedure lhapdf_number_pdf
          module procedure lhapdf_number_name
      end interface

      private

      type(pdf), allocatable, save :: pdfs(:)
!$omp threadprivate(pdfs)

      interface
          subroutine lhapdf_info() bind(C, name="lhapdf_info")
          end subroutine

          function c_lhapdf_loadMember(setname, member) bind(C, name="lhapdf_loadMember")
              use iso_c_binding
              implicit none
              type(c_ptr) :: c_lhapdf_loadMember 
              character(kind=c_char), dimension(*), intent(in) :: setname
              integer(c_int), value, intent(in) :: member
          end function

          function c_lhapdf_evolve(pdfptr, x, q2, flav) bind(C, name="lhapdf_evolve")
              use iso_c_binding
              implicit none
              real(c_double) :: c_lhapdf_evolve
              type(c_ptr), value, intent(in) :: pdfptr
              real(c_double), value, intent(in) :: x, q2
              integer(c_int), value, intent(in) :: flav
          end function

          function c_lhapdf_number(setname) bind(C, name="lhapdf_number")
              use iso_c_binding
              implicit none
              integer(c_int) :: c_lhapdf_number
              character(kind=c_char), dimension(*), intent(in) :: setname
          end function

          function c_lhapdf_alphas(pdfptr, q2) bind(C, name="lhapdf_alphas")
              use iso_c_binding
              implicit none
              real(c_double) :: c_lhapdf_alphas
              type(c_ptr), value, intent(in) :: pdfptr
              real(c_double), value, intent(in) :: q2
          end function

          function c_lhapdf_orderqcd(pdfptr) bind(C, name="lhapdf_orderqcd")
              use iso_c_binding
              implicit none
              integer(c_int) :: c_lhapdf_orderqcd
              type(c_ptr), value, intent(in) :: pdfptr
          end function

          function c_lhapdf_computeUncertainty(pdfptr, array, num) bind(C, name="lhapdf_computeUncertainty")
              use iso_c_binding
              import :: pdfuncertainty
              implicit none
              type(pdfuncertainty) :: c_lhapdf_computeUncertainty
              type(c_ptr), value, intent(in) :: pdfptr
              type(c_ptr), value, intent(in) :: array
              integer, value, intent(in) :: num
          end function

          subroutine c_lhapdf_getconfig(pdfptr, key, outvalue, olen) bind(C, name="lhapdf_getconfig")
              use iso_c_binding
              implicit none
              type(c_ptr), value, intent(in) :: pdfptr
              character(kind=c_char), dimension(*), intent(in) :: key
              character(kind=c_char), dimension(*), intent(inout) :: outvalue
              integer(c_int), value, intent(in) :: olen
          end subroutine
      end interface

      contains

      subroutine lhapdf_getconfig(pdf_in, key, outvalue)
          implicit none
          class(pdf), intent(in) :: pdf_in
          character(len=*), intent(in) :: key
          character(len=*), intent(inout) :: outvalue

          call c_lhapdf_getconfig(pdf_in%ptr, key//c_null_char, outvalue, len(outvalue))
      end subroutine

      function computeUncertainty(array, pdfnum)
          implicit none
          type(pdfuncertainty) :: computeUncertainty
          integer, intent(in) :: pdfnum
          real(dp), target, intent(in) :: array(:)

          computeUncertainty = c_lhapdf_computeUncertainty(pdfs(pdfnum)%ptr, c_loc(array), size(array))
      end function

      function pdf_orderqcd(pdf_in)
          implicit none
          integer :: pdf_orderqcd
          class(pdf), intent(in) :: pdf_in

          pdf_orderqcd = c_lhapdf_orderqcd(pdf_in%ptr)
      end function

      function lhapdf_number_pdf(pdf_in)
          implicit none
          integer :: lhapdf_number_pdf
          class(pdf), intent(in) :: pdf_in

          lhapdf_number_pdf = c_lhapdf_number(trim(pdf_in%setname)//c_null_char)
      end function

      function lhapdf_number_name(setname)
          implicit none
          integer :: lhapdf_number_name
          character(*), intent(in) :: setname

          lhapdf_number_name = c_lhapdf_number(trim(setname)//c_null_char)
      end function


      function lhapdf_loadMember(setname,member)
          implicit none
          type(pdf) :: lhapdf_loadMember

          character(*), intent(in) :: setname
          integer, intent(in) :: member

          lhapdf_loadMember%setname = trim(setname)
          lhapdf_loadMember%ptr = c_lhapdf_loadMember(trim(setname)//c_null_char, member)

      end function

      function pdf_xfxq2(pdf_in, x, q2, flav)
          implicit none
          real(dp) :: pdf_xfxq2
          class(pdf), intent(in) :: pdf_in
          real(dp), intent(in) :: x,q2
          integer, intent(in) :: flav

          pdf_xfxq2 = c_lhapdf_evolve(pdf_in%ptr, x, q2, flav)
      end function

      function pdf_alphas(pdf_in, q)
          implicit none
          real(dp) :: pdf_alphas
          class(pdf), intent(in) :: pdf_in
          real(dp), intent(in) :: q

          pdf_alphas = c_lhapdf_alphas(pdf_in%ptr, q**2)
      end function


!!!! routines affecting global pdfs state

      function getorderqcd()
          use PDFerrors, only: currentPDF
          implicit none
          integer :: getorderqcd

          getorderqcd = pdfs(currentPDF)%orderqcd()
      end function

      ! get alphas from central pdf
      function getalphas(q)
          use PDFerrors, only: currentPDF
          implicit none
          real(dp) :: getalphas
          real(dp), intent(in) :: q

          getalphas = pdfs(currentPDF)%alphas(q)
      end function

      subroutine initcentral(setnames, setmembers)
          use omp_lib
          implicit none
          character(len=*), intent(in) :: setnames(:)
          integer, intent(in) :: setmembers(:)
          integer :: j

!$omp parallel
          if (.not. allocated(pdfs)) then
              allocate(pdfs(0:size(setnames)-1))
              do j=1,size(setnames)
!$omp critical(loadMember)
                  pdfs(j-1) = lhapdf_loadMember(trim(setnames(j))//c_null_char,setmembers(j))
!$omp end critical(loadMember)
              enddo
          endif
!$omp end parallel
      end subroutine

      subroutine initall(setnames, setmembers)
          use omp_lib
          implicit none
          character(len=*), intent(in) :: setnames(:)
          integer, intent(in) :: setmembers(:)

          integer :: nmax
          integer :: j,k
          integer :: accum
          integer, allocatable :: memberCounts(:)

          allocate(memberCounts(size(setnames)))
          do j=1,size(setnames)
              memberCounts(j) = lhapdf_number(trim(setnames(j)))
          enddo
          nmax = sum(memberCounts(:))

!$omp parallel private(accum)
          if (.not. allocated(pdfs)) then
              allocate(pdfs(0:nmax-1))

              accum = 0
              do j=1,size(setnames)
                  do k=0, memberCounts(j) - 1
!$omp critical(loadMember)
                      pdfs(accum) = lhapdf_loadMember(trim(setnames(j))//c_null_char, k)
!$omp end critical(loadMember)
                      accum = accum + 1
                  enddo
              enddo
          endif
!$omp end parallel
      end subroutine

      subroutine fdist_lhapdf(ih, x, xmu, fx)
          use PDFerrors, only: currentPDF
          implicit none
          integer, intent(in) :: ih
          ! ih = +1 for proton
          ! ih = -1 for anti-proton
          real(dp), intent(in) :: x, xmu
          real(dp), intent(out) :: fx(-5:5)

          integer :: i

          if (x > 1d0) then
              fx(:) = 0._dp
              return
          endif

          if (ih ==1) then
              do i=-5,5
                  fx(i) = pdfs(currentPDF)%xfxq2(x, xmu**2, i) / x
              enddo
          else
              do i=-5,5
                  fx(i) = pdfs(currentPDF)%xfxq2(x, xmu**2, -i) / x
              enddo
          endif

      end subroutine

end module
