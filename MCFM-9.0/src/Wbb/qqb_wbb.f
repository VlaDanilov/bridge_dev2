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
 
      subroutine qqb_wbb(p,msq)
      implicit none
      include 'types.f'
c---  Matrix elements squared 
c     q(-p1)+qb(-p2) --> nu(p3)+e^+(p4)+b(p5)+bb(p6)
c---  averaged(summed) over initial(final) colours and spins
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ckm.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'masses.f'
      include 'first.f'
      integer:: j,k
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf),msqwbb
      real(dp):: qqb,qbq
      
      if (first) then
       write(6,*)
       write(6,*) '****************** Process info ********************'
       write(6,*) '*                                                  *'
       write(6,*) '* mb=0 for this process, although cuts are applied *'
       write(6,*) '* to simulate the effect of the b-mass:            *'
       write(6,*) '*                                                  *'
       write(6,99) ' *                pt(b) > ',sqrt(mbsq),
     &  '                *'
       write(6,99) ' *                m(bb) > ',two*sqrt(mbsq),
     &  '                *'
       write(6,*) '****************************************************'
       first=.false.
      endif
      
C---Initialize to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

C---Fill spinor products
      call spinoru(6,p,za,zb)

c ensure that we have a hard process
      if (
     &      (s(5,6) < four*mbsq) 
     & .or. (s(1,5)*s(2,5)/s(1,2) < mbsq) 
     & .or. (s(1,6)*s(2,6)/s(1,2) < mbsq) ) return
      
C--calculate matrix element squared
      qqb=msqwbb(1,2,5,6)
      qbq=msqwbb(2,1,5,6)

      do j=-(nf-1),(nf-1)
      do k=-(nf-1),(nf-1)
      if     ((j > 0) .and. (k < 0)) then
               msq(j,k)=Vsq(j,k)*qqb
      elseif ((j < 0) .and. (k > 0)) then
               msq(j,k)=Vsq(j,k)*qbq
      endif
      enddo
      enddo

      return
   
   99 format(a26,f6.3,a21)   
      
      end







