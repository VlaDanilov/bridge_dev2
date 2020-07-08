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
 
      subroutine phase5h(r,p1,p2,p3,p4,p5,p6,p7,wt)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'mxdim.f'
      include 'breit.f'
      include 'kprocess.f'
c********* generate phase space for 2-->5 process
c********* r(mxdim),p1(4),p2(4) are inputs 
c--------- incoming p1 and p2 reversed in sign from physical values 
c---- i.e. phase space for -p1-p2 --> p3+p4+p5+p6+p7
c---- with all 2 pi's (ie 1/(2*pi)^11)

      real(dp):: r(mxdim),wt,wt34567,wt34,wt567,wt67,tmp
      real(dp):: p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),p7(4)
      real(dp):: p567(4),p12(4),p67(4),p34(4),p56(4),s3min
      integer:: oldn3
      real(dp), parameter:: wt0=1._dp/twopi**3, tiny=1e-6_dp

      p12(:)=-p1(:)-p2(:)

      oldn3=n3
      n3=1
      call phi1_2(r(1),r(2),r(3),r(4),p12,p567,p34,wt34567,*99)
      n3=0
      call phi3m0(r(5),r(6),p34,p3,p4,wt34,*99)
c--- handle case where 5 and 6 are massive b-quarks
      if (((kcase == kWHbbdk) .and. (mb > tiny))
     &   .or.((kcase==kZHbbdk).and.(mb > tiny))) then
        s3min=four*mb**2
        call phi1_2m(0._dp,r(7),r(8),r(11),s3min,p567,p7,p56,wt567,*99)
        call phi3m(r(12),r(13),p56,p5,p6,mb,mb,wt67,*99)
c--- otherwise, 5, 6 and 7 are all massless
      else      
        call phi1_2m(0._dp,r(7),r(8),r(11),0._dp,p567,p5,p67,wt567,*99)
        call phi3m0(r(12),r(13),p67,p6,p7,wt67,*99)
      endif
      
      wt=wt0*wt34567*wt567*wt67*wt34
  
      n3=oldn3
  
      return
 99   continue
      wt=0._dp
      n3=oldn3
      return
      end

