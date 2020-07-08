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
 
      subroutine qqb_z2jetx_swap(pin,msq,msqd1,msqx,msqd2)
      implicit none
      include 'types.f'
c--- this is just a wrapper routine to qqb_dirgam_g,
c--- that interchanges p4 and p5 before the call
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: nu
      real(dp):: pin(mxpart,4),msq(-nf:nf,-nf:nf),
     & pswap(mxpart,4)
      real(dp):: msqd1(0:2,-nf:nf,-nf:nf),msqd2(0:2,-nf:nf,-nf:nf)
     & ,msqx(0:2,-nf:nf,-nf:nf)
     
      pswap(:,:)=pin(:,:)
      pswap(5,:)=pin(6,:)
      pswap(6,:)=pin(5,:)
      
      call qqb_z2jetx(pswap,msq,msqd1,msqx,msqd2)
      
      return
      end
      
