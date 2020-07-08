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
 
      subroutine VV_Hgaga_gs(p,msqc)
      implicit none
      include 'types.f'
      
c--- Weak Bosion Fusion : sums up WW and VV contributions
************************************************************************
*     Author: J. M. Campbell                                           *
*     July, 2002.                                                      *
*                                                                      *
*     Weak Boson Fusion : sums up WW and ZZ contributions              *
*     This routine calculates the dipole subtraction terms             *
*     for the process:                                                 *
*     q(-p1) + q(-p2) --> H(p34) + q(p5) + q(p6)                       *
*                                                                      *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ptilde.f'
      real(dp):: p(mxpart,4),msqc(maxd,-nf:nf,-nf:nf),
     & msqc_ww(maxd,-nf:nf,-nf:nf),msqc_zz(maxd,-nf:nf,-nf:nf)
      integer:: j,k,nd
 
      call WW_Hgaga_gs(p,msqc_ww)
      call ZZ_Hgaga_gs(p,msqc_zz)
      
      do nd=1,ndmax
      do j=-nf,nf
      do k=-nf,nf
        msqc(nd,j,k)=msqc_ww(nd,j,k)+msqc_zz(nd,j,k)
      enddo
      enddo      
      enddo
      
      return
      end
      
