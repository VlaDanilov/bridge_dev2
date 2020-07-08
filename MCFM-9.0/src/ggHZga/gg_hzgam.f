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
 
      subroutine gg_hzgam(p,msq)
      implicit none
      include 'types.f'
      
C----Author: R.K.Ellis July 2012
c----Lowest order matrix element for H production
c----in the heavy quark (mt=Infinity) limit.
C----averaged over initial colours and spins
c    g(-p1)+g(-p2)-->H --> Z/gamma^*--(l(p3)+a(p4)) + gamma(p5) 
c---
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'sprods_com.f'
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4)
      real(dp):: decay,gg,Asq,HZgamMSQ

c---set msq=0 to initialize
      msq(:,:)=0._dp

      call dotem(5,p,s)
      decay=HZgamMSQ(3,4,5)
     & /((s(1,2)-hmass**2)**2+(hmass*hwidth)**2)
      
      Asq=(as/(3._dp*pi))**2/vevsq
      gg=0.5_dp*V*Asq*s(1,2)**2

c---calculate propagators
      msq(0,0)=avegg*gg*decay

      return
      end


