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
 
      subroutine qqb_QQb_mix_sft(p,msq)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'constants.f'
      include 'sprods_com.f'
      include 'qcdcouple.f'
      include 'masses.f'
      integer j,k,ii,ff
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf),pp(mxpart,4),
     .     msqb(-nf:nf,-nf:nf),ss(mxpart,mxpart)

      pp(1:4,:)=p(1:4,:)
      call qqb_QQb_mix(pp,msqb)

      call dotem(5,p,s)

      ss(:,:) = s(:,:)/2._dp

      msq = 0._dp

      do j = -nf,nf
         k = -j
         if (j .ne. 0) then
            msq(j,k) = + msq(j,k)
     .           + msqb(j,k)*(
     .           + (s(1,3)/(s(1,5)+s(3,5)))/s(1,5)
     .           - (s(2,3)/(s(2,5)+s(3,5)))/s(2,5)
     .           - (s(1,4)/(s(1,5)+s(4,5)))/s(1,5)
     .           + (s(2,4)/(s(2,5)+s(4,5)))/s(2,5)
     .           + (s(3,1)/(s(3,5)+s(1,5)) - mt**2/2._dp/s(3,5))/s(3,5)
     .           - (s(4,1)/(s(4,5)+s(1,5)) - mt**2/2._dp/s(4,5))/s(4,5)
     .           - (s(3,2)/(s(3,5)+s(2,5)) - mt**2/2._dp/s(3,5))/s(3,5)
     .           + (s(4,2)/(s(4,5)+s(2,5)) - mt**2/2._dp/s(4,5))/s(4,5)
!     .           + s(1,3)/s(1,5)/s(3,5)
!     .           - s(2,3)/s(2,5)/s(3,5)
!     .           - s(1,4)/s(1,5)/s(4,5)
!     .           + s(2,4)/s(2,5)/s(4,5)
     .           )    

c--- overall factor here is canonical one for eikonal;
c--- combined with 1/16/pi^2 from phase-space it gives (as/2/pi)    
         msq(j,k)=2._dp*gsq*(
     .           + s(1,3)/s(1,5)/s(3,5)
     .           - s(2,3)/s(2,5)/s(3,5)
     .           - s(1,4)/s(1,5)/s(4,5)
     .           + s(2,4)/s(2,5)/s(4,5))
         msq(j,k)=msq(j,k)*msqb(j,k)
         end if
      end do

      msq = 4._dp*msq

c--- soft limit does not affect the born assuming sum_{i=1,4}p_i = 0

      end subroutine qqb_QQb_mix_sft
            

