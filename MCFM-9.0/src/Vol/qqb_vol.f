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
 
      subroutine qqb_vol(P,msq)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      
      real(dp):: P(mxpart,4),msq(-nf:nf,-nf:nf)
      integer:: j,k
      call dotem(mxpart,p,s)
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

      msq(2,-1)=1._dp
      return
      end

c      subroutine qqb_vol_g(P,msq)
c      implicit none
c      include 'types.f'
c       
c      include 'constants.f'
c      include 'nf.f'
c      include 'mxpart.f'
c      include 'cplx.h'
c      include 'masses.f'
c      include 'sprods_com.f'
c      real(dp):: P(mxpart,4),msq(-nf:nf,-nf:nf)
c      integer:: j,k,N
c      N=7
c      call dotem(N,p,s)
c      do j=-nf,nf
c      do k=-nf,nf
c      msq(j,k)=0._dp
c      enddo
c      enddo

c      msq(2,-1)=1._dp
c      return
c      end

c      subroutine qqb_vol_gs(P,msq)
c      implicit none
c      include 'types.f'
c       
c      include 'constants.f'
c      include 'nf.f'
c      include 'mxpart.f'
c      include 'cplx.h'
c      include 'masses.f'
c      include 'ptilde.f'
c      include 'sprods_com.f'
c      real(dp):: P(mxpart,4),msq(maxd,-nf:nf,-nf:nf)
c      integer:: j,k,nd,N
c      N=7
c      call dotem(N,p,s)
c      do nd=1,maxd
c      do j=-nf,nf
c      do k=-nf,nf
c      msq(nd,j,k)=0._dp
c      enddo
c      enddo
c      enddo
c      ndmax=0
c      
c      return
c      end
