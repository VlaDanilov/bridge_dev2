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
 
      subroutine Ampvirt_gggg(j1,j2,j3,j4,A)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_com.f'
      complex(dp):: A(3,2,2,2,2)
      complex(dp):: A1Hggggpppp,A1Hggggpmmm,A1Hggggmmpp,A1Hggggmpmp
      integer:: j1,j2,j3,j4,i1(3),i2(3),i3(3),i4(4),j

      i1(1)=j1      
      i2(1)=j2      
      i3(1)=j3      
      i4(1)=j4      

      i1(2)=j1      
      i2(2)=j2      
      i3(2)=j4      
      i4(2)=j3      

      i1(3)=j1      
      i2(3)=j4      
      i3(3)=j2      
      i4(3)=j3      

      A=0._dp

      do j=1,3
      A(j,2,2,2,2)=A1Hggggpppp(i1(j),i2(j),i3(j),i4(j),za,zb)
      A(j,1,1,1,1)=A1Hggggpppp(i1(j),i2(j),i3(j),i4(j),zb,za)

      A(j,2,1,1,1)=A1Hggggpmmm(i1(j),i2(j),i3(j),i4(j),za,zb)
      A(j,1,2,1,1)=A1Hggggpmmm(i2(j),i3(j),i4(j),i1(j),za,zb)
      A(j,1,1,2,1)=A1Hggggpmmm(i3(j),i4(j),i1(j),i2(j),za,zb)
      A(j,1,1,1,2)=A1Hggggpmmm(i4(j),i1(j),i2(j),i3(j),za,zb)

      A(j,1,2,2,2)=A1Hggggpmmm(i1(j),i2(j),i3(j),i4(j),zb,za)
      A(j,2,1,2,2)=A1Hggggpmmm(i2(j),i3(j),i4(j),i1(j),zb,za)
      A(j,2,2,1,2)=A1Hggggpmmm(i3(j),i4(j),i1(j),i2(j),zb,za)
      A(j,2,2,2,1)=A1Hggggpmmm(i4(j),i1(j),i2(j),i3(j),zb,za)

      A(j,1,1,2,2)=A1Hggggmmpp(i1(j),i2(j),i3(j),i4(j),za,zb)
      A(j,1,2,2,1)=A1Hggggmmpp(i4(j),i1(j),i2(j),i3(j),za,zb)
      A(j,2,2,1,1)=A1Hggggmmpp(i3(j),i4(j),i1(j),i2(j),za,zb) 
      A(j,2,1,1,2)=A1Hggggmmpp(i2(j),i3(j),i4(j),i1(j),za,zb)

      A(j,1,2,1,2)=A1Hggggmpmp(i1(j),i2(j),i3(j),i4(j),za,zb)
      A(j,2,1,2,1)=A1Hggggmpmp(i2(j),i3(j),i4(j),i1(j),za,zb)
      enddo

      return
      end
