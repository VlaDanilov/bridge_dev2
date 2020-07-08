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
 
      subroutine spinork(N,p,zabikj,zbaikj,k)
      implicit none
      include 'types.f'
c---Calculate spinor products dotted in with a vector k
c---extended to deal with negative energies ie with all momenta outgoing
C   zabikj=<i-|k|j-> zbaikj=<i+|k|j+> 
c---Arbitrary conventions of Bern, Dixon, Kosower, Weinzierl, 
c---za(i,j)*zb(j,i)=s(i,j)
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: p(mxpart,4),rt(mxpart),k(4),kp,km,flip(mxpart)
      complex(dp):: pr(mxpart),pl(mxpart),f(mxpart),kr,kl,
     & zabikj(mxpart,mxpart),zbaikj(mxpart,mxpart)
      integer:: i,j,N
      
C--setup components for vector which is contracted in
      kp=+k(4)+k(1)
      km=+k(4)-k(1)
      kr=cplx2(+k(3),-k(2))
      kl=cplx2(+k(3),+k(2))

c---if one of the vectors happens to be zero this routine fails.
      do j=1,N
            zabikj(j,j)=czip
            zbaikj(j,j)=czip
            
C-----positive energy case
            if (p(j,4) > 0._dp) then
            flip(j)=1._dp
            f(j)=cone
            else
            flip(j)=-1._dp
            f(j)=im
            endif
            rt(j)=sqrt(flip(j)*(p(j,4)+p(j,1)))
            pr(j)=cplx2(flip(j)*p(j,3),-flip(j)*p(j,2))
            pl(j)=conjg(pr(j))
      enddo
      do i=1,N
         do j=1,i
         zabikj(i,j)=f(i)*f(j)
     & *(pr(i)*pl(j)*cplx1(kp/(rt(i)*rt(j)))
     &    -pr(i)*kl*cplx1(rt(j)/rt(i))
     &    -cplx1(rt(i)/rt(j))*kr*pl(j)+cplx1(rt(i)*rt(j)*km))
         zbaikj(j,i)=zabikj(i,j) 
         zabikj(j,i)=flip(i)*flip(j)*conjg(zabikj(i,j))
         zbaikj(i,j)=flip(i)*flip(j)*conjg(zbaikj(j,i))

      enddo
      enddo

      return
      end
