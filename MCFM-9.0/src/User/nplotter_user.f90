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
 
subroutine nplotter_user(pjet, wt,wt2, nd)
    use types
    implicit none
!************************************************************************
!*   Subroutine that is called to allow the user to bin their own
!*   histograms
!*
!*   Variables passed to this routine:
!*   
!*     pjet:  4-momenta of incoming partons(i=1,2), outgoing leptons and
!*            jets in the format p(i,4) with the particles numbered
!*            according to the input file and components labelled by
!*            (px,py,pz,E)
!*   
!*       wt:  weight of this event
!*   
!*      wt2:  weight^2 of this event
!*   
!*       nd:  an integer:: specifying the dipole number of this contribution
!*            (if applicable), otherwise equal to zero
!************************************************************************
    include 'mxpart.f'
    include 'ptilde.f'
    include 'nplot.f'
    integer :: nplotmax
    common/nplotmax/nplotmax
    include 'first.f'

    real(dp), intent(in) :: pjet(mxpart,4)
    real(dp), intent(in) :: wt,wt2
    integer, intent(in) ::  nd

    integer, parameter :: tagbook=1, tagplot=2
    integer :: tag, n


    if (first) then
       ! mode for initialization
       ! initialize observables with dummy values
       tag = tagbook
    else                
       ! mode for binning events
       tag = tagplot
    end if

    n = nextnplot-1 ! nextnplot denotes the number for the _next_ plot

    ! call bookplot routines here, increase n for each call

    ! n = n + 1
    ! call bootplot(n, ...)
    ! n = n + 1
    ! call bootplot(n, ...)


    if (first) then
        first = .false.
        nextnplot = n+1
        nplotmax = n
    endif
      
end subroutine

