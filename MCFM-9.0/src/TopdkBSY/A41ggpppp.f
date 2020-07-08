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
 
      function A41ggpppp(e1,p2,p3,e4,za,zb,zab,zba)
      implicit none
      include 'types.f'
      complex(dp):: A41ggpppp
      
C-----Authors: John Campbell and Keith Ellis, November 2011
C---- taken from arXiv:1101.5947 [hep-ph], Eq. (4)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'zabprods_decl.f'
      complex(dp):: ALggpppp,ARggpppp,
     & BSYA1fggpppp,BSYA1Hggpppp
      integer:: e1,p2,p3,e4
      real(dp):: nlf,nhf
      parameter(nlf=5d0,nhf=1d0)
 
      A41ggpppp=ALggpppp(e1,p2,p3,e4,za,zb,zab,zba)
     & -Ncinv**2*ARggpppp(e1,p2,p3,e4,za,zb,zab,zba)
     & -BSYA1fggpppp(e1,p2,p3,e4,za,zb,zab,zba)*nlf/xn
     & -BSYA1Hggpppp(e1,p2,p3,e4,za,zb,zab,zba)*nhf/xn

      return
      end

