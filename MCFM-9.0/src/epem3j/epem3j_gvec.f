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
 
      subroutine epem3j_gvec(p,n,in,msq)
      implicit none
      include 'types.f'
      
c--- simple modification of qqb_w1jet_gs.f: permuted 1 and 4, 2 and 3
c--- to switch leptons with quarks and added a factor of Nc      
      
c----Matrix element for W production
C----averaged over initial colours and spins
c    contracted with the vector v(mu)
C For nwz=+1
c     u(-p1)+dbar(-p2)--> g(p5)+ W^+(n(p3)+e^+(p4))
C For nwz=-1
c     d(-p1)+ubar(-p2)--> g(p5)+ W^-(e^-(p3)+nbar(p4))
c---
c---ip emitter
c---kp spectator
c---in label of gluon which is contracted with n
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'sprods_com.f'
      integer:: j,k,in
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4)
      real(dp):: w1jetn,p1p2,n(4),fac

      fac=two*gsq*V*xn*gwsq**2
      call dotem(5,p,s)

      if (in == 5) then     
        p1p2=+aveqq*fac*w1jetn(4,3,2,1,5,p,n)
      else
        write(6,*) 'Error in epem3j_gvec.f: in = ',in
        stop
      endif

      msq(0,0)=p1p2

      return
      end
 
