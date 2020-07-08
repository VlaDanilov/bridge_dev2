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
 
      subroutine gagaAGTYassemble1loop(ss,tt,uu,M1finM0,BoldC0,M0sq,M1rensq)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'scale.f' 
      real(dp), intent(in) :: ss,tt,uu
      complex(dp), intent(in) :: M1finM0,BoldC0,M0sq
      real(dp), intent(out) :: M1rensq
      real(dp) :: AGTYG1s,BigX,BigY

c--- AGTY function for |M1fin|^2 pieces
      BigX=log(-tt/ss)
      BigY=log(-uu/ss)
      M1rensq=AGTYG1s(tt,uu,BigX,BigY)
      
      M1rensq=M1rensq+two*real(BoldC0*conjg(M1finM0),dp)
     & +real(BoldC0*conjg(BoldC0)*M0sq,dp)

c--- restore overall color factor
      M1rensq=M1rensq*CF**2
      
      return
      end
      
      
      
