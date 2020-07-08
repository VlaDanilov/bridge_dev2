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
 
module Integration
    use types
    implicit none

    enum, bind(c)
          enumerator :: nloReal=1, nloVirt=2, nnloBelow=3, &
               nnloVirtAbove=4, nnloRealAbove=5, &
               snloBelow = 6, snloAbove = 7, lord=8, nloFrag=9, &
               nloRealExtra = 10
    endenum

    ! adjust maxParts to be the maximum number used above
    integer, parameter :: maxParts = 10
    integer, parameter :: maxIPS = 2

    integer, parameter :: ndim_incr(10) = [3,1,2,1,3,2,0,0,1,3]

    character(*), parameter :: partStrings(10) = [ &
        "NLO real           ", &
        "NLO virt           ", &
        "NNLO below cut     ", &
        "NNLO virt above cut", &
        "NNLO real above cut", &
        "SNLO below cut     ", &
        "SNLO above cut     ", &
        "LO                 ", &
        "NLO frag           ", &
        "NLO real extra     "]

end module
