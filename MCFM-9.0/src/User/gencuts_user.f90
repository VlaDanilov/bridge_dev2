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
 
submodule (m_gencuts) m_gencuts_user
    implicit none


    contains

    module function reweight_user(pjet)
        use types
        implicit none
        include 'mxpart.f'
        real(dp) :: reweight_user
        real(dp), intent(in) :: pjet(mxpart,4)

        real(dp) :: ptpure

        reweight_user = 1._dp

    end function

    module logical function gencuts_user(pjet, njets)
      use types
      implicit none
      include 'mxpart.f'
      real(dp), intent(in) :: pjet(mxpart,4)
      integer, intent(in) :: njets

      real(dp) :: ptpure

      gencuts_user = .false.

      ! implement your own cuts here

!     if (njets /= 2) then
!         gencuts = .true.
!         return
!     endif

    end function

end submodule
