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
 
      subroutine ggdilep(p,msq)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'constants.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'sprods_com.f'
      real(dp) :: p(mxpart,4),msq(fn:nf,fn:nf),t1,t2,fac

      call dotem(3,p,s)

      t1=-s(1,3)/s(1,2)
      t2=-s(2,3)/s(1,2)

      msq(:,:)=zip
!      esq=fourpi/137.03599911_dp
      fac=esq**2*q1**4/four !four is spin-ave
      msq(0,0) = fac*eight*(t1/t2+t2/t1)

      return
      end
