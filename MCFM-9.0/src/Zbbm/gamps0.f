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
 
      subroutine gamps0(mq,p1,p2,p3,p4,t5,t6,sum,coupL,coupR)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: p1,p2,p3,p4,t5,t6
      integer:: hz,h1,h2,h5,h6
      complex(dp)::
     & gg_a(2,2,2,2,2),gg_b(2,2,2,2,2),gg_c(2,2,2,2,2),
     & gg_d(2,2,2,2,2),gg_e(2,2,2,2,2),gg_f(2,2,2,2,2),
     & gg_g(2,2,2,2,2),gg_h(2,2,2,2,2),
     & sum1(2),sum2(2),sum0(2),coupL,coupR
      real(dp):: sum,mq

      call gampsabc(mq,p1,p2,p3,p4,t5,t6,gg_a,gg_b,gg_c)
      call gampsdef(mq,p1,p2,p3,p4,t5,t6,gg_d,gg_e,gg_f)
      call gampsgh(mq,p1,p2,p3,p4,t5,t6,gg_g,gg_h)

      sum=0._dp

      do h1=1,2      
      do h2=1,2      
      do h5=1,2      
      do h6=1,2      

      do hz=1,2      
      sum1(hz)=
     & +gg_a(hz,h1,h2,h5,h6)
     & +gg_b(hz,h1,h2,h5,h6)
     & +gg_c(hz,h1,h2,h5,h6)
     & -gg_g(hz,h1,h2,h5,h6)
     & -gg_h(hz,h1,h2,h5,h6)
      sum2(hz)=
     & +gg_d(hz,h1,h2,h5,h6)
     & +gg_e(hz,h1,h2,h5,h6)
     & +gg_f(hz,h1,h2,h5,h6)
     & +gg_g(hz,h1,h2,h5,h6)
     & +gg_h(hz,h1,h2,h5,h6)
      sum0(hz)=
     & +gg_a(hz,h1,h2,h5,h6)
     & +gg_b(hz,h1,h2,h5,h6)
     & +gg_c(hz,h1,h2,h5,h6)
     & +gg_d(hz,h1,h2,h5,h6)
     & +gg_e(hz,h1,h2,h5,h6)
     & +gg_f(hz,h1,h2,h5,h6)
      enddo

      sum=sum
     & +abs(sum1(1)*coupL+sum1(2)*coupR)**2
     & +abs(sum2(1)*coupL+sum2(2)*coupR)**2
     & -abs(sum0(1)*coupL+sum0(2)*coupR)**2/xn**2
             
      enddo
      enddo
      enddo
      enddo

      return
      end
