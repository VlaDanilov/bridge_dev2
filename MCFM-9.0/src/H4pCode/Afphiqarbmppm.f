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
 
      function Afphiqarbmppm(j1,j2,j3,j4,za,zb) 
      implicit none
      include 'types.f'
      complex(dp):: Afphiqarbmppm
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'epinv.f'
      include 'scale.f'
      real(dp):: s12,s34
      complex(dp):: lnrat,A0phiqarbmppm
      integer:: j1,j2,j3,j4
      s12=s(j1,j2)
      s34=s(j3,j4)
C---arXIv:09060008v1 Eq.(4.12)
c--- the function defined in this routine is in fact (-i*A_4),
c---   i.e. complete LHS
      Afphiqarbmppm=A0phiqarbmppm(j1,j2,j3,j4,za,zb)*
     & (-2._dp/3._dp*(2._dp*epinv+lnrat(musq,-s12)+lnrat(musq,-s34))
     & - 20._dp/9._dp)

      return
      end
