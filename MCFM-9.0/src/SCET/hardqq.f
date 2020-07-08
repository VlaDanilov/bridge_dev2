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
 
      subroutine hardqq(Qsq,musq,hard)
      implicit none
!    Hard function for qqbar in units of as/2/pi
      include 'types.f'
      include 'constants.f'
      real(dp),intent(in)::Qsq,musq
      real(dp),intent(out)::hard(2)
      complex(dp)::coeff(2)
      call qqcoeff(Qsq,musq,coeff)
! factors of 1/2 adjust for as/4/pi -> as/2/pi
      hard(1)=real(coeff(1),kind=dp)
      hard(2)=half**2*real(coeff(1)*conjg(coeff(1)),kind=dp)
     & +half*real(coeff(2),kind=dp)
      return
      end
