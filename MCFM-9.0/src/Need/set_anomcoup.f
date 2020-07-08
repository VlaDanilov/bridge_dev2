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
 
      subroutine set_anomcoup(p)
        implicit none
        include 'types.f'
        include 'mxpart.f'
        include 'anomcoup.f'
        include 'masses.f'

        real(dp), intent(in) :: p(mxpart,4)
        real(dp) :: shat, dotvec, xfac

        ! by using s345 (the lepton pair + the photon) we don't have to
        ! recompute the form factor for different parts
        shat = dotvec(p(3,:)+p(4,:)+p(5,:),p(3,:)+p(4,:)+p(5,:))

        if (tevscale > 0._dp) then
          xfac = 1._dp/(1 + shat/(tevscale*1.e3_dp)**2)
        else
          xfac = 1._dp
        endif

        hitZ(1) = xfac**3*h1Z/zmass**2
        hitZ(2) = xfac**4*h2Z/zmass**4
        hitZ(3) = xfac**3*h3Z/zmass**2
        hitZ(4) = xfac**4*h4Z/zmass**4
        hitgam(1) = xfac**3*h1gam/zmass**2
        hitgam(2) = xfac**4*h2gam/zmass**4
        hitgam(3) = xfac**3*h3gam/zmass**2
        hitgam(4) = xfac**4*h4gam/zmass**4


      end subroutine
