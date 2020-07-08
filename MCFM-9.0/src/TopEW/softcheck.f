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
 
      subroutine softcheck(p)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'constants.f'
      integer j,k
      real(dp):: p(mxpart,4),msqg(-nf:nf,-nf:nf),
     .     msqs(-nf:nf,-nf:nf)

      call writeout(p)

      write(6,*) 'j, k, msqg, msqs, g/s: '

      call qqb_QQb_mix_g(p,msqg)
      call qqb_QQb_mix_sft(p,msqs)


      do j=-nf,nf
         k = -j
         if (j .ne. 0) then
            write(6,*) j, k, 
     .           msqg(j,k),msqs(j,k),msqg(j,k)/msqs(j,k)
         end if
      end do

      pause

      end subroutine softcheck
