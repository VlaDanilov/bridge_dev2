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
 
! This is an example for studying H+jet at high pT.

! The cuts limit the maximum pT to 1000 GeV.
! The reweighting reweights with exp(ptH/65).

! This flattens the sampling of pT distribution, resulting in
! similar relative integration uncertainties at low and high pT.

! By using the reweighting factor the output from the Vegas integration
! does not correspond to a physical result. But all results in the
! histograms (including the total cross section) are correct.
! All the reweight_user does is to manually tune the Vegas importance sampling.

submodule (m_gencuts) m_gencuts_user
    implicit none


    logical :: studyHjetmass = .false.

    contains

    module function reweight_user(pjet)
        use types
        implicit none
        include 'mxpart.f'
        include 'first.f'
        include 'runstring.f'
        real(dp) :: reweight_user
        real(dp), intent(in) :: pjet(mxpart,4)

        real(dp) :: ptpure

        if (first) then
            studyHjetmass = (index(runstring,"hjetmass") > 0)
            first = .false.
        endif

        if (studyHjetmass .eqv. .true.) then
            reweight_user = exp(ptpure(pjet(3,:) + pjet(4,:))/65._dp)
            return
        endif
        
        reweight_user = 1._dp

    end function

    module logical function gencuts_user(pjet, njets)
        use types
        implicit none
        include 'mxpart.f'
        include 'first.f'
        include 'runstring.f'
        real(dp), intent(in) :: pjet(mxpart,4)
        integer, intent(in) :: njets

        real(dp) :: ptpure

        gencuts_user = .false.

          if (first) then
              studyHjetmass = (index(runstring,"hjetmass") > 0)
              first = .false.
          endif

        if (studyHjetmass) then
          if (ptpure(pjet(3,:) + pjet(4,:)) > 1000._dp) then
              gencuts_user = .true.
              return
          endif
        endif

        gencuts_user = .false.

    end function

end submodule
