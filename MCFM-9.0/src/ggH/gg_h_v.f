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
 
      subroutine gg_h_v(p,msqv)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
      include 'scale.f'
      real(dp):: msq(-nf:nf,-nf:nf),msqv(-nf:nf,-nf:nf),
     & p(mxpart,4),dot,xl12
      integer:: j,k
      
      do j=-nf,nf
      do k=-nf,nf
      msqv(j,k)=0._dp
      enddo
      enddo

      call gg_h(p,msq)
      xl12=log(two*dot(p,1,2)/musq)

c--sum of virtual diagram in DRED and UV counterterm including 
C--term required to bring into the MSbar scheme
      scheme='dred'

      msqv(0,0)=ason2pi*xn*2._dp*(
     &  -epinv*(epinv2-xl12)-0.5_dp*xl12**2+11._dp/6._dp+0.5_dp*pisq
     &  -((11._dp-two*real(nf,dp)/xn)*epinv-1._dp)/6._dp)*msq(0,0)

      return
      end
     
