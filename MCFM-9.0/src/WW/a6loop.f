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
 
      function a6loops(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: a6loops
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
c---  DKS Eq. 3.15
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: a6loopa
      a6loops=a6loopa(j1,j2,j3,j4,j5,j6,za,zb)
     &       +a6loopa(j1,j2,j6,j5,j4,j3,za,zb)
      return
      end
        
      function a6loopa(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: a6loopa
      
c---  DKS Eq. 2.10 for alpha = a
c---  note that (-i) included in A(alpha),A(tree,alpha)
c---  so no factor of (+i) in front of F(alpha)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: tree,Vpole,a6treea,fa
      tree=a6treea(j1,j2,j3,j4,j5,j6,za,zb)
      a6loopa=tree*Vpole(s(1,2))+fa(j1,j2,j3,j4,j5,j6,za,zb)
      return 
      end

c      function a6loopb(j1,j2,j3,j4,j5,j6,za,zb)
c      implicit none
c      include 'types.f'
c      complex(dp):: a6loopb
c      
c---  DKS Eq. 2.10 for alpha = b
c---  note that (-i) included in A(alpha),A(tree,alpha)
c---  so no factor of (+i) in front of F(alpha)
c      include 'constants.f'
c      include 'nf.f'
c      include 'mxpart.f'
c      include 'cplx.h'
c      include 'sprods_com.f'
c      include 'zprods_decl.f'
c      integer:: j1,j2,j3,j4,j5,j6
c      complex(dp):: tree,Vpole,a6treeb
        
c      tree=a6treeb(j1,j2,j3,j4,j5,j6,za,zb)    
c      a6loopb=tree*Vpole(s(j1,j2))
c
c      return 
c      end

              

        
