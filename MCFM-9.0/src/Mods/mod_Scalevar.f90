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
 
module Scalevar
    use Integration
    use types
    implicit none

    integer, public, save :: maxscalevar = -1
    logical, public, save :: doScalevar = .false.

    logical, public, save :: foundpow
!$omp threadprivate(foundpow)
    integer, public, save :: alphaspow
!$omp threadprivate(alphaspow)
    real(dp), public, save :: scalereweight(6)
!$omp threadprivate(scalereweight)

    real(dp), public, parameter :: &
        scalevarmult(1:6) = [2d0, 0.5d0, 2d0, 0.5d0, 1d0, 1d0], &
        facscalevarmult(1:6) = [2d0, 0.5d0, 1d0, 1d0, 2d0, 0.5d0]

    include 'initialscales.f'

    public :: usescales, updateAlphas

    private

    contains

    ! forces an update of all variables related to alphas
    subroutine updateAlphas(renscale_in)
        use constants
        implicit none
        include 'qcdcouple.f'
        include 'couple.f'
        include 'nlooprun.f'
        real(dp), intent(in) :: renscale_in
        real(dp) :: alphas

        as = alphas(renscale_in,amz,nlooprun)
        ason2pi=as/(2._dp*pi)
        ason4pi=as/(4._dp*pi)
        gsq=4._dp*pi*as

    end subroutine

    subroutine usescales(renscale_in,facscale_in)
        use constants
        implicit none
        real(dp), intent(in) :: renscale_in, facscale_in
        include 'nf.f'
        include 'qcdcouple.f'
        include 'stopscales.f'
        include 'facscale.f'
        include 'scale.f'
        include 'couple.f'
        include 'nlooprun.f'

        real(dp) :: alphas
        logical :: mustUpdateAlphas

        if (scale == renscale_in .and. facscale == facscale_in) then
            return
        endif

        mustUpdateAlphas = .false.

        if (scale /= renscale_in) then
            mustUpdateAlphas = .true.
        endif
      
        scale=renscale_in
        facscale=facscale_in

        if  (scale > 60000._dp) scale=60000._dp
        if  (facscale > 60000._dp) facscale=60000._dp
        if  (scale < 1._dp) scale=1._dp
        if  (facscale < 1._dp) facscale=1._dp

        ! these are additional scales used in the t-channel single top + b routines
        facscale_H=facscale
        facscale_L=facscale
        renscale_H=scale
        renscale_L=scale

        musq=scale**2

        if (mustUpdateAlphas) then
            call updateAlphas(scale)
        endif

    end subroutine

end module
