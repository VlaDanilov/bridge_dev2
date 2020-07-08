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
 
      subroutine boost(mass,p1,p_in,p_out)
      implicit none
      include 'types.f'
c     take momenta p_in in frame in which particle one is at rest with mass 
c     "mass" 
c     and convert to frame in which particle one has fourvector p1
      
      real(dp):: mass,p1(4),p_in(4),p_out(4)
      real(dp):: gam,beta(1:3),bdotp,one
      parameter(one=1._dp)
      integer:: j,k
      gam=p1(4)/mass
      bdotp=0._dp
      do j=1,3
      beta(j)=-p1(j)/p1(4)
      bdotp=bdotp+p_in(j)*beta(j)
      enddo
      p_out(4)=gam*(p_in(4)-bdotp)
      do k=1,3
      p_out(k)=p_in(k)+gam*beta(k)*(gam/(gam+one)*bdotp-p_in(4))
      enddo
      return
      end      

