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
 
      subroutine spinorz(N,p,za,zb)
      implicit none
      include 'types.f'
c---Calculate spinor products
c---extended to deal with negative energies ie with all momenta outgoing
c---Arbitrary conventions of Bern, Dixon, Kosower, Weinzierl, 
c---za(i,j)*zb(j,i)=s(i,j)
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      real(dp):: p(mxpart,4)
      complex(dp):: c12(mxpart),f(mxpart),rt(mxpart)
      integer:: i,j,N
      
c---if one of the vectors happens to be zero this routine fails.
      do j=1,N
         za(j,j)=czip
         zb(j,j)=za(j,j)

C-----positive energy case
            rt(j)=cplx1(p(j,4)+p(j,3))
            if (real(rt(j)) == 0._dp) then
            write(6,*) 'spinorz:j',j
            write(6,*) 'spinorz fails for momenta directed along z axis'
            stop
            endif
            rt(j)=sqrt(rt(j))
            c12(j)=cplx2(p(j,1),p(j,2))
            f(j)=cone
      enddo
      do i=2,N
         do j=1,i-1
         s(i,j)=two*(p(i,4)*p(j,4)-p(i,1)*p(j,1)
     &              -p(i,2)*p(j,2)-p(i,3)*p(j,3))
         za(i,j)=f(i)*f(j)
     &   *(c12(i)*(rt(j)/rt(i))-c12(j)*(rt(i)/rt(j)))

! Disabled this check August 2018 since it can cause poor cancellation
!         if (abs(s(i,j))<1.e-9_dp) then
!         zb(i,j)=-(f(i)*f(j))**2*conjg(za(i,j))
!         else
         zb(i,j)=-s(i,j)/za(i,j)
!         endif
         za(j,i)=-za(i,j)
         zb(j,i)=-zb(i,j)
         s(j,i)=s(i,j)
         enddo
      enddo

      return
      end
