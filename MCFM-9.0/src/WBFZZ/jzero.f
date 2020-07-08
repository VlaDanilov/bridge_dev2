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
 
      subroutine jzero(n2,n1,zab,zba,j0)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      complex(dp):: zab(mxpart,4,mxpart),zba(mxpart,4,mxpart),
     & j0(4,2)
C---The one Z-current multiplied by i
C---order of indices Lorentz,quark-line helicity

      integer:: n1,n2,nu
      do nu=1,4
      j0(nu,1)=zab(n2,nu,n1)
      j0(nu,2)=zba(n2,nu,n1)
      enddo
      return
      end
