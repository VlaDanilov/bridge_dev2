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
 
      subroutine spinorextend(za,zb)
      implicit none
      include 'types.f'
      
C-----Authors: John Campbell and Keith Ellis, November 2011
C-----given spinor products calculated for heavy quark auxiliary vectors
C-----eta_1(e1) and eta_4(e4), return spinor products extended to yield also products
c-----for the opposite polarization of the auxiliary vector.
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'masses.f'
      include 'etadef.f'
      integer:: j

      do j=1,6
      za(e1p,j)=zb(e1,1)*za(1,j)/mt 
      za(j,e4p)=-za(j,4)*zb(4,e4)/mt 
      za(j,e1p)=-za(e1p,j)
      za(e4p,j)=-za(j,e4p)

c--- use c.c. to obtain zb from za      
      zb(j,e1p)=za(e1,1)*zb(1,j)/mt 
      zb(e4p,j)=-zb(j,4)*za(4,e4)/mt
      zb(e1p,j)=-zb(j,e1p)
      zb(j,e4p)=-zb(e4p,j)
      enddo
      za(e1p,e4p)=-zb(e1,1)*za(1,4)*zb(4,e4)/mt**2
      za(e4p,e1p)=-za(e1p,e4p)
c--- use c.c. to obtain zb from za      
      zb(e4p,e1p)=+za(e1,1)*zb(1,4)*za(4,e4)/mt**2
      zb(e1p,e4p)=-zb(e4p,e1p)
      return
      end
