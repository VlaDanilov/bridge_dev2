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
 
      function ALggpppp(e1,p2,p3,e4,za,zb,zab,zba)
      implicit none
      include 'types.f'
      complex(dp):: ALggpppp
      
C-----Authors: John Campbell and Keith Ellis, November 2011
C---- taken from arXiv:1101.5947 [hep-ph], Eq. (92)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'zabprods_decl.f'
      include 'epinv.f'
      include 'scale.f'
      real(dp):: mt2
      complex(dp):: BSYA0ggpppp,BSYALggpppp,VL,lnrat
      integer:: e1,p2,p3,e4
 
      mt2=mt**2
      VL=-2d0*cplx1(epinv**2)+cplx1(0.5d0*epinv)
     & -epinv*(lnrat(musq,-s(1,p2))+lnrat(mt2,-s(1,p2)))
     & -epinv*lnrat(musq,-s(p2,p3))
      ALggpppp=VL*BSYA0ggpppp(e1,p2,p3,e4,za,zb,zab,zba)
     &           +BSYALggpppp(e1,p2,p3,e4,za,zb,zab,zba)
      return
      end
