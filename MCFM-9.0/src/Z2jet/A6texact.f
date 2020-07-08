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
 
      function A6texact(s23,mt2)
        use mod_qcdloop_c
      implicit none
      include 'types.f'
      include 'constants.f'
      complex(dp)::A6texact
      real(dp)::s23,mt2,musq

c--- musq is irrelevant, set it to some value
      musq=abs(s23)

c--- this is the form of the top loop vacuum polarization obtained
c--- after renormalization by subtraction at zero momentum transfer
      A6texact=-two/three*(
     &  (one+two*mt2/s23)*(qlI2(s23,mt2,mt2,musq,0)-qlI2(0._dp,mt2,mt2,musq,0))
     &   +one/three)

      return
      end
