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
 
      function GammaHbb0(Msq,mbsq)
      implicit none
      include 'types.f'
      real(dp):: GammaHbb0      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      real(dp):: Msq,mbsq,beta,besq
c      write(6,*) 'GammaHbb0:Msq',Msq
c      write(6,*) 'GammaHbb0:mbsq',mbsq
c      pause
      besq=1._dp-4._dp*mbsq/Msq
      beta=sqrt(besq)
      GammaHbb0=3._dp/4._dp/pi*mbsq*Gf/rt2*sqrt(Msq)*beta**3
      return
      end

