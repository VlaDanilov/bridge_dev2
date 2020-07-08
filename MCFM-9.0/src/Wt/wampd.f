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
 
      subroutine wampd(mq,qwidth,p1,pn,pe,pb,t1,amp)
      implicit none
      include 'types.f'
c     t(t1/p1) --> Pn(pn)+Pe^+(pe)+b(pb)
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer:: p1,pe,pn,pb,t1
      complex(dp):: amp(2)
      integer:: j
      real(dp):: propd,propt,tsq,mq,qwidth

      propd=sqrt((s(pe,pn)-wmass**2)**2+(wmass*wwidth)**2)
      tsq  =s(pe,pn)+s(pe,pb)+s(pn,pb)
      propt=sqrt((tsq-mq**2)**2+(mq*qwidth)**2) 

c      write(6,*) propd,tsq,propt

      amp(:)=zip

c--- label on amplitudes represents heavy quark helicity
c---  amp(1)= negative helicity, amp(2)= positive helicity    
      amp(1) =
     &  - za(pn,pb)*zb(pe,t1)

      amp(2) =
     &  - 1/(zb(p1,t1))*za(pn,pb)*zb(p1,pe)*mq

      do j=1,2
      amp(j)=amp(j)/propd/propt
      enddo

      return
      end 
