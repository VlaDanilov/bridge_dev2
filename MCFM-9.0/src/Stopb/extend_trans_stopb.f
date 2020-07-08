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
 
      subroutine extend_trans_stopb(pold,p,ptrans,pext)
      implicit none
      include 'types.f'
      
c--- take vector 
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: pold(mxpart,4),p(mxpart,4),ptrans(mxpart,4),
     & pext(mxpart,4),p3(4),p4(4),p5(4),pt(4),ptt(4),
     & p3out(4),p4out(4),p5out(4)
      integer:: j,nu
      
      do nu=1,4
        pt(nu)=p(3,nu)
        ptt(nu)=ptrans(3,nu)
        p3(nu)=pold(3,nu)
        p4(nu)=pold(4,nu)
        p5(nu)=pold(5,nu)
      enddo
      
      call boostx(p3,pt,ptt,p3out)
      call boostx(p4,pt,ptt,p4out)
      call boostx(p5,pt,ptt,p5out)

      do nu=1,4
        pext(1,nu)=ptrans(1,nu)
        pext(2,nu)=ptrans(2,nu)
        pext(3,nu)=p3out(nu)
        pext(4,nu)=p4out(nu)
        pext(5,nu)=p5out(nu)
        pext(6,nu)=ptrans(4,nu)
        pext(7,nu)=ptrans(5,nu)
        do j=8,mxpart
        pext(j,nu)=0._dp
        enddo
      enddo
            
      return
      end
      
