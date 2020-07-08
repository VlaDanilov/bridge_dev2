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
 

      subroutine scalar_dm(i1,i2,za,zb,amp) 
      implicit none
      include 'types.f'
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      include 'dm_params.f' 
      include 'zprods_decl.f' 
      integer:: i1,i2
      complex(dp):: amp(2,2) 
      real(dp):: Bp,Beta,sab

      amp(1,2)=czip 
      amp(2,1)=czip
      
      sab=Dble(za(i1,i2)*zb(i2,i1))
!      write(6,*) 'in scalar dm ',sab 
!      pause
      Beta=one-4d0*xmass**2/sab 
      Beta=sqrt(beta) 
      Bp=one+Beta
      
      amp(1,1)=-za(i1,i2)*Bp+xmass**2/zb(i2,i1)/Bp     
      amp(2,2)=zb(i2,i1)*Bp-xmass**2/za(i1,i2)/Bp

      return 
      end 
