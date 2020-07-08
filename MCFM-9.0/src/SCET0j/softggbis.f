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
 
!#################################################################
!##### O(as) soft contribution to tau_1, calculates soft1(-1:1)
!##### contributions to delta(tau1), L0(tau1) and L1(tau1)
!##### returns integrated contribution.
!##### a factor of alphas/2/pi has been extracted out
!##### see Eq. (2.22) of 1012.4480,
!#################################################################
! This file is created using the color factor swap result of
! Idilbi et al., hep-ph/0605068v2
      subroutine softggbis(order,soft1,soft2)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'scale.f'
      include 'nf.f'
      include 'scet_const.f'
      integer :: order
      real(dp) :: soft1(-1:1),soft2(-1:3)

      if (order < 1) return
c--- coefficient of (alphas/2/pi)
!---- contribution to delta(tau1)
      soft1(-1) = CA*pisqo6
!---- contribution to L0(tau1)
      soft1(0) = 0.0_dp
!---- contribution to L1(tau1)
      soft1(1) =-8*CA

!---- L taucut
c      Lt1cut = log(taucut/scale)
c      soft(1)=soft1(-1)
c     & +soft1(0)*Lt1cut
c     & +soft1(1)*Lt1cut**2/two

      if (order < 2) return
c--- coefficient of (alphas/2/pi)^2
c--- note overall factor of 4*CA from Eq. (2.22) compared with Eq. (2.23)

!---- contribution to delta(tau1)
      soft2(-1) =
     & +CA**2*( - 27/ten*zeta2**2)
     & +CA*(- 5/27._dp*be0 - 160/27._dp*CA - 37/12._dp*be0*zeta2
     &      + 2*CA*zeta2 + 22/five*CA*zeta2**2 + 58/12._dp*be0*zeta3)
!---- contribution to L0(tau1)
      soft2(0) = 
     & + CA**2*64*zeta3
     & + CA*(28/nine*be0 - 2*be0*zeta2 + 32/nine*CA - 14*CA*zeta3)
!---- contribution to L1(tau1)
      soft2(1) = 
     & - CA**2*72*zeta2
     & + CA*(- 20/three*be0 - 16/three*CA + 8*CA*zeta2)
!---- contribution to L2(tau1)
      soft2(2) = 4*CA*be0
!---- contribution to L3(tau1)
      soft2(3) = 32*CA**2

c      soft(2)=soft2(-1)
c     & +soft2(0)*Lt1cut
c     & +soft2(1)*Lt1cut**2/two
c     & +soft2(2)*Lt1cut**3/three
c     & +soft2(3)*Lt1cut**4/four
      return
      end


