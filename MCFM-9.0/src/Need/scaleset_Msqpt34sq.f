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
 
      subroutine scaleset_Msqpt34sq(p,mu0)
      implicit none
      include 'types.f'
c--- subroutine to calculate dynamic scale equal to
c---  sqrt(M^2+pt34^2), where M is the mass of the particle (34)
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'kprocess.f'
      include 'breit.f'
      real(dp):: p(mxpart,4),mu0,pttwo

      if((kcase==kW_only) .or.
     &   (kcase==kZ_only) .or.
     &   (kcase==kW_1jet) .or.
     &   (kcase==kW_2jet) .or.
     &   (kcase==kW_3jet) .or.
     &   (kcase==kZ_1jet) .or.
     &   (kcase==kZ_2jet) .or.
     &   (kcase==kZ_3jet) .or.
     &   (kcase==kWbbbar) .or.
     &   (kcase==kWbbmas) .or.
     &   (kcase==kZbbbar) .or.
     &   (kcase==kggfus0) .or.
     &   (kcase==kggfus1) .or.
     &   (kcase==khjetma) .or.
     &   (kcase==kggfus2) .or.
     &   (kcase==kggfus3) .or.
     &   (kcase==kgagajj) .or.
     &   (kcase==kh2jmas) .or.
     &   (kcase==khttjet) .or.
     &   (kcase==kHigaga) .or.
     &   (kcase==kHgagaj)) then
        mu0=mass3**2+pttwo(3,4,p)**2
        mu0=sqrt(abs(mu0))
      else
        write(6,*) 'dynamicscale sqrt(M^2+pt34^2)'//
     &             ' not supported for this process.'
        stop
      endif
      
      return
      end
      
