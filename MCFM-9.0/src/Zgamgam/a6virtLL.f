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
 
      function a6virtLL(st,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: a6virtLL
****************************************************
* virtual amplitude for
* 0->q(p1)+qb(p2)+gam(p3)+gam(p4)+lb*p5)+l(p6)
* where both photon are coming from lepton line
****************************************************
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      character*9 st
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: a6treeg,vll
c-----
      a6virtLL = a6treeg(st,j1,j2,j3,j4,j5,j6,za,zb)
     &              *vll(st,j5,j6,j3,j4,j1,j2)
c-----
      return
      end


      function vll(st,j1,j2,j3,j4,j5,j6)
      implicit none
      include 'types.f'
      complex(dp):: vll
c-----divergent part
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'epinv.f'
      include 'scale.f'
      include 'epinv2.f'
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: Lnrat,xl12
      character*9 st
c-----
      xl12=Lnrat(musq,-s(j1,j2))
c-----
      vll = -4._dp
     &      -(epinv*epinv2+epinv*xl12+half*xl12**2)
     &      -(3._dp/2._dp)*(epinv+xl12)
c-----
      return
      end
