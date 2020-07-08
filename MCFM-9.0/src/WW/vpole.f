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
 
      function Vpole(sij)
      implicit none
      include 'types.f'
      complex(dp):: Vpole
      
c---  DKS Eq. 2.12
      include 'epinv.f'
      include 'epinv2.f'
      include 'scale.f'
      real(dp):: sij
      complex(dp):: Lnrat,xl12
        
      xl12=Lnrat(-sij,musq)

      Vpole=-epinv*epinv2+epinv*(-1.5_dp+xl12)
     &   -0.5_dp*xl12**2+1.5_dp*xl12-3.5_dp

      return
      end
