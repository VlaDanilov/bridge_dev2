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
 
C--- Implementation of arXiv:0906.0008v1, Eqs. (A.8), (A.9) and (A.10)
c--- with factor of i removed

      function A0phiAgQgmppm(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A0phiAgQgmppm
      
c--- (A.8) first line
c----with factor of i removed
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4
      
      A0phiAgQgmppm=-za(j1,j4)**2/(za(j1,j2)*za(j2,j3))
      
      return
      end
      
      function A0phiAgQgmmpp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A0phiAgQgmmpp
      
c--- (A.8) second line
c----with factor of i removed
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4
      
      A0phiAgQgmmpp=-za(j1,j2)**2/(za(j3,j4)*za(j4,j1))
      
      return
      end
      
      function A0phidAgQgmppm(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A0phidAgQgmppm
      
c--- (A.8) third line
c----with factor of i removed
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4
      
      A0phidAgQgmppm=zb(j2,j3)**2/(zb(j3,j4)*zb(j4,j1))
      
      return
      end
      
      function A0phidAgQgmmpp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A0phidAgQgmmpp
      
c--- (A.8) fourth line
c----with factor of i removed
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4
      
      A0phidAgQgmmpp=zb(j3,j4)**2/(zb(j1,j2)*zb(j2,j3))
      
      return
      end
      
      function A0phiAgQgmmpm(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A0phiAgQgmmpm
      
c--- (A.9)
c----with factor of i removed
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: j1,j2,j3,j4
      complex(dp):: zab2
      real(dp):: s3
      s3(j1,j2,j3)=s(j1,j2)+s(j2,j3)+s(j3,j1)
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      
      A0phiAgQgmmpm=-zab2(j4,j1,j2,j3)**2
     &               /(zb(j1,j2)*zb(j2,j3)*s3(j1,j2,j3))
     &              -zab2(j2,j1,j4,j3)**2
     &               /(zb(j3,j4)*zb(j4,j1)*s3(j3,j4,j1))
      
      return
      end
      
      function A0phidAgQgmppp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A0phidAgQgmppp
      
c--- (A.9)
c----with factor of i removed
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: j1,j2,j3,j4
      complex(dp):: zab2
      real(dp):: s3
      s3(j1,j2,j3)=s(j1,j2)+s(j2,j3)+s(j3,j1)
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      
      A0phidAgQgmppp=+zab2(j1,j2,j3,j4)**2
     &                /(za(j1,j2)*za(j2,j3)*s3(j1,j2,j3))
     &               +zab2(j1,j3,j4,j2)**2
     &                /(za(j3,j4)*za(j4,j1)*s3(j3,j4,j1))
      
      return
      end
      
