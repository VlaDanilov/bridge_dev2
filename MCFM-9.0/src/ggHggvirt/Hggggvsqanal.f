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
 
      function Hggggvsqanal(j1,j2,j3,j4)
      implicit none
      include 'types.f'
      real(dp):: Hggggvsqanal
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
C---  matrix element squared for 0 --> H + g(j1)+g(j2)+g(j3)+g(j4)
      integer:: j1,j2,j3,j4,j,h1,h2,h3,h4
      complex(dp):: Avirt(3,2,2,2,2),Alo(3,2,2,2,2)
      real(dp):: temp,ren,H4prenorm
      
      call   Amplo(j1,j2,j3,j4,Alo)
      call Ampvirt_gggg(j1,j2,j3,j4,Avirt)
c--- get renormalization factor
      ren=H4prenorm()

      temp=0._dp
      do j=1,3
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2   
      Avirt(j,h1,h2,h3,h4)=Avirt(j,h1,h2,h3,h4)+ren*Alo(j,h1,h2,h3,h4)      

      temp=temp+real(Avirt(j,h1,h2,h3,h4)*conjg(Alo(j,h1,h2,h3,h4)))
      enddo
      enddo
      enddo
      enddo
      enddo

      Hggggvsqanal=xn**3*V/2._dp*temp
                  
      return
      end


