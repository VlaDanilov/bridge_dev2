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
 
      subroutine qqb_Hg_v(p,msq)
      implicit none
      include 'types.f'
      
c---Matrix element squared averaged over initial colors and spins
c     parton(-p1)+parton(-p2) --> H(p)+parton(p5)
c                                  |
c                                   --> b(p3)+bb(p4)
c                            
c--all momenta incoming
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'msbarmasses.f'
      include 'ewcouple.f'
      include 'susycoup.f'
      include 'scheme.f'
      include 'scale.f'
      include 'couple.f'
      integer:: j,k
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),s
      real(dp):: coupsq_eff,ghbb_eff
      real(dp):: mb_eff,massfrun
      real(dp):: fac,propsq,hdecay,bbghvirt

c--susycoup is the deviation of Higgs coupling 
c-- from the standard model value

c--statement function 
      s(j,k)=2._dp*(p(j,4)*p(k,4)-p(j,1)*p(k,1)
     &           -p(j,2)*p(k,2)-p(j,3)*p(k,3))
C---ur-amplitude is b(p1)+bbar(p2)+g(p3)+H(q)
c      amp(j1,j2,j3)=4._dp*(s(j1,j2)**2+s(3,4)**2)/(s(j1,j3)*s(j2,j3))

      scheme='dred'

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

      if (s(3,4) < 4._dp*mbsq) return

c--- run mb to appropriate scale
      mb_eff=massfrun(mb_msbar,scale,amz,2)
c      mb_eff=mb_msbar
      
      call hbbdecay(p,3,4,hdecay)
      hdecay=hdecay*susycoup**2
      propsq=1._dp/((s(3,4)-hmass**2)**2+(hmass*hwidth)**2)
c--- The _eff couplings include the running mass
c--- We need to separate these from the factors associated with the
c--- Higgs decay, because the Br. Ratio does not include running mb
      ghbb_eff=sqrt(esq/xw)*mb_eff/2._dp/wmass
      coupsq_eff=susycoup**2*ghbb_eff**2

      fac=coupsq_eff*propsq*hdecay
c      fac=CF*xn*gsq*coupsq*propsq*hdecay

      msq(0,+5)=+fac*aveqg*bbghvirt(2,5,1)
      msq(0,-5)=+fac*aveqg*bbghvirt(5,2,1)
      msq(+5,0)=+fac*aveqg*bbghvirt(1,5,2)
      msq(-5,0)=+fac*aveqg*bbghvirt(5,1,2)

      return
      end
          
    

      
     
