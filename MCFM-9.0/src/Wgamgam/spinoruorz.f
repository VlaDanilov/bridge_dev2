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
 
      subroutine spinoruorz(N,p,za,zb)
      implicit none
      include 'types.f'
c--- This routine just provides an easy way of switching between
c--- spinoru (normal MCFM running) and spinorz (checks of virtual)
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      real(dp):: p(mxpart,4)
      integer:: N,j,k

      call spinoru(N,p,za,zb)
c      call spinorz(N,p,za,zb)
      
      return
      end
