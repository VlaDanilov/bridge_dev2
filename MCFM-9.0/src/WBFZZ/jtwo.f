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
 
      subroutine jtwo(n2,n3,n4,n5,n6,n1,za,zb,zab,zba,j2,jw2)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      complex(dp):: zab(mxpart,4,mxpart),zba(mxpart,4,mxpart),
     & j2(4,2,2,2,2),j2_34_56_1(4,2,2,2,2),j2_56_34_1(4,2,2,2,2),
     & jw2(4,2,2,2),jw2_34_56_1(4,2,2,2),jw2_56_34_1(4,2,2,2)
      integer:: n1,n2,n3,n4,n5,n6,jdu,h21,h34,h56
C---The two Z-current divided by (-i)
      call jtwo3456(n2,n3,n4,n5,n6,n1,
     & za,zb,zab,zba,j2_34_56_1,jw2_34_56_1)
      call jtwo3456(n2,n5,n6,n3,n4,n1,
     & za,zb,zab,zba,j2_56_34_1,jw2_56_34_1)

      do jdu=1,2
      do h56=1,2
      do h34=1,2
      jw2(1:4,jdu,h34,h56)=
     &  jw2_34_56_1(1:4,jdu,h34,h56)
     & +jw2_56_34_1(1:4,jdu,h56,h34)
      do h21=1,2
      j2(1:4,jdu,h21,h34,h56)=
     &  j2_34_56_1(1:4,jdu,h21,h34,h56)
     & +j2_56_34_1(1:4,jdu,h21,h56,h34)
      enddo
      enddo
      enddo
      enddo

      return
      end 


