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
 
      subroutine storedip(msq_dip,msq_dipv,dsub,dsubv,
     &                    sub_dip,sub_dipv,n)
      implicit none
      include 'types.f'
c--- this routine transfers the information on the colour
c--- structure from a common block into separate arrays for
c--- each parton configuration
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'msq_cs.f'
      include 'msqv_cs.f'
      integer:: i,j,k,n
      real(dp):: msq_dip(6,0:2,-nf:nf,-nf:nf),dsub(4),sub_dip(6,4),
     &           msq_dipv(6,0:2,-nf:nf,-nf:nf),dsubv,sub_dipv(6)
      
      do i=0,2
        do j=-nf,nf
        do k=-nf,nf
          msq_dip(n,i,j,k)=msq_cs(i,j,k)
          msq_dipv(n,i,j,k)=msqv_cs(i,j,k)
        enddo
        enddo
      enddo
      
      do i=1,4
        sub_dip(n,i)=dsub(i)
      enddo
      sub_dipv(n)=dsubv
      
      return
      end
