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
 
      function A41qqppmp(e1,p2,p3,e4,za,zb,zab,zba)
      implicit none
      include 'types.f'
      complex(dp):: A41qqppmp
C-----Authors: John Campbell and Keith Ellis, November 2011
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'zabprods_decl.f'
      complex(dp):: A1lcqqppmp,A1slcqqppmp,
     & A1fqqppmp,A1Hqqppmp
      integer:: e1,p2,p3,e4,e1p,e4p
      real(dp):: nlf,nhf
      parameter(nlf=5d0,nhf=1d0)

c--- perform the swaps 5 <-> 7 and 6 <-> 8     
      e1p=12-e1
      e4p=14-e4

c--- changed the sign of the term
c----   -2d0*Ncinv**2*A1lcqqppmp(e1p,p3,p2,e4p,zb,za,zba,zab)  
c--- although the sign appears to be the same as written in Eq. (11) of
c---  arXiv:1101.5947, because of complex conjugation we would have
c---  expected the fourth term in this expression to have the opposite sign 
      A41qqppmp=
     & (1d0-2d0*Ncinv**2)*A1lcqqppmp(e1,p2,p3,e4,za,zb,zab,zba)
     & -A1fqqppmp(e1,p2,p3,e4,za,zb,zab,zba)*nlf/xn
     & -A1Hqqppmp(e1,p2,p3,e4,za,zb,zab,zba)*nhf/xn
     & -2d0*Ncinv**2*A1lcqqppmp(e1p,p3,p2,e4p,zb,za,zba,zab)
     & -Ncinv**2*A1slcqqppmp(e1,p2,p3,e4,za,zb,zab,zba)

      return
      end

