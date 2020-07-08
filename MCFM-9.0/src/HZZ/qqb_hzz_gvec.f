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
 
      subroutine qqb_hzz_gvec(p,n,in,msq)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
C  ip is the label of the emitting parton
C  kp is the label of the spectator parton
      integer:: j,k,in
      real(dp):: msq(-nf:nf,-nf:nf),msqt(-nf:nf,-nf:nf)
      real(dp):: n(4),nDn,p(mxpart,4)

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

      nDn=n(4)**2-n(3)**2-n(2)**2-n(1)**2
      call qqb_hzz(p,msqt)

      msq(0,0)=-0.5_dp*nDn*msqt(0,0)      

      return
      end

