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
 
      function A41ggppmp(e1,p2,p3,e4,za,zb,zab,zba)
      implicit none
      include 'types.f'
      complex(dp):: A41ggppmp
      
C-----Authors: John Campbell and Keith Ellis, November 2011
C---- taken from arXiv:1101.5947 [hep-ph], Eq. (4)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'zabprods_decl.f'
      complex(dp):: ALggppmp,ARggppmp,
     & BSYA1fggppmp,BSYA1Hggppmp
      integer:: e1,p2,p3,e4
      real(dp):: nlf,nhf
      parameter(nlf=5d0,nhf=1d0)
 
      A41ggppmp=ALggppmp(e1,p2,p3,e4,za,zb,zab,zba)
     & -Ncinv**2*ARggppmp(e1,p2,p3,e4,za,zb,zab,zba)
     & -BSYA1fggppmp(e1,p2,p3,e4,za,zb,zab,zba)*nlf/xn
     & -BSYA1Hggppmp(e1,p2,p3,e4,za,zb,zab,zba)*nhf/xn

      return
      end

