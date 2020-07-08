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
 
      module qsort_m
          use iso_c_binding
          implicit none

          interface
              subroutine qsort_c(array,elem_count,elem_size,compare) bind(C,name="qsort")
                  import
                  type (c_ptr), value :: array
                  integer (c_size_t), value :: elem_count
                  integer (c_size_t), value :: elem_size
                  type (c_funptr), value :: compare 
              end subroutine
          end interface

          type, bind(C) :: s_jetvalue
              real(c_double) :: jetvalue
              integer(c_int) :: jetindex
          end type

          public :: qsort_dp, qsort_dp_with_jetindex
          private

          contains

          subroutine qsort_dp_with_jetindex(array, jetindices)
              use types
              use iso_c_binding
              implicit none

              real(dp), intent(inout) :: array(:)
              integer, intent(inout) :: jetindices(:)

              ! assume array and jetindices have same size
              type (s_jetvalue), target :: jetvalues(size(array))

              integer (c_size_t) :: elem_size, elem_count

              jetvalues(:)%jetvalue = array(:)
              jetvalues(:)%jetindex = jetindices

              elem_size = c_sizeof(jetvalues(1))
              elem_count = size(jetvalues)

              call qsort_c(c_loc(jetvalues), elem_count, elem_size, c_funloc(compare_dp_byfirst))

              array(:) = jetvalues(:)%jetvalue
              jetindices(:) = jetvalues(:)%jetindex


          end subroutine

          subroutine qsort_dp(array)
              use types
              use iso_c_binding
              implicit none
              real(dp), intent(inout), target :: array(:)

              integer (c_size_t) :: dpsize
              integer (c_size_t) :: elem_count

              ! size of dp in bytes
              dpsize = c_sizeof(array(1))
              elem_count = size(array)

              call qsort_c(c_loc(array), elem_count, dpsize, c_funloc(compare_dp))
          end subroutine

          ! this sorts largest by smallest
          function compare_dp_byfirst(d1,d2) bind(C)
              use types
              use iso_c_binding
              implicit none

              type (s_jetvalue), intent(in) :: d1,d2
              integer (c_int) :: compare_dp_byfirst

              if (d1%jetvalue > d2%jetvalue) then
                  compare_dp_byfirst = -1
              elseif (d2%jetvalue > d1%jetvalue) then
                  compare_dp_byfirst = 1
              else
                  compare_dp_byfirst = 0
              endif

          end function

          function compare_dp(d1,d2) bind(C)
              use types
              use iso_c_binding
              implicit none

              real (c_double), intent(in) :: d1, d2
              integer (c_int) :: compare_dp 

              if (d1 > d2) then
                  compare_dp = 1
              elseif (d2 > d1) then
                  compare_dp = -1
              else
                  compare_dp = 0
              endif
          end function

      end module
