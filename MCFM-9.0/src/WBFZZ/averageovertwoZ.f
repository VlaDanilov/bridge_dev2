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
 
      subroutine averageovertwoZ(fxn,p,msq)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'masses.f'
      integer:: kcount,jcount
      real(dp):: msq(-nf:nf,-nf:nf),msq1(-nf:nf,-nf:nf),
     & p(mxpart,4),pout1(mxpart,4),pout(mxpart,4),statfac
      pout(:,:)=0d0
      msq(:,:)=0d0
C---sum over 36 values;
      do jcount=1,6
      call pgen(p,jcount,3,4,pout1)
      do kcount=1,6
      call pgen(pout1,kcount,5,6,pout)
      call fxn(pout,msq1)
      msq(:,:)=msq(:,:)+msq1(:,:)
      enddo
      enddo
      statfac=0.5d0
      msq(:,:)=msq(:,:)*(zwidth**2/(4d0*esq*(le**2+re**2)))**2*statfac
      return
      end
