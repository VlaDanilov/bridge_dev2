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
 
      subroutine dk2qqb_QQb_g(p,msq)
      implicit none
      include 'types.f'
      
************************************************************************
*     Author: R.K. Ellis                                               *
*     January, 2012.                                                   *
*     calculate the element squared and subtraction terms              *
*     for the process                                                  *
*                                                                      *
*     q(-p1) +qbar(-p2)=nu(p3)+e+(p4)+b(p5)+bbar(p6)+e-(p7)+nubar(p8)  *
*     +g(p9) radiated from antitop quark in decay                      *
*                                                                      *
*     Top is kept strictly on-shell although all spin correlations     *
*     are retained.                                                    *
*                                                                      *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'plabel.f'
      integer:: j,k,hb,hc,h12,j1,j2,h1,h2,hg,j1max
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4)
      real(dp):: fac,qqb,gg
      complex(dp)::  prop,
     & mtop(2,2),manti(2,2,2),mprod(2,2,2),mtot(2,2,2,2),
     & mabtot(2,2,2,2,2),mbatot(2,2,2,2,2),mqed(2,2,2,2,2),
     & mab(2,2,2,2),mba(2,2,2,2)
       parameter(j1max=1)

C----set all elements to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

      call toppaironshell(p,-1,mprod,mab,mba)
      call tdecay(p,3,4,5,mtop)
      call adecayg(p,7,8,6,9,manti)

c--- q-qbar amplitudes
      do hb=1,2
      do hg=1,2
      do hc=1,2
      do h12=1,2
      mtot(hb,hg,hc,h12)=czip
      do j1=1,j1max
      do j2=1,2
      mtot(hb,hg,hc,h12)=mtot(hb,hg,hc,h12)+
     & mtop(hb,j1)*mprod(j1,h12,j2)*manti(j2,hg,hc)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      prop=cplx2(zip,mt*twidth)**2
      fac=V*gwsq**4*gsq**2/abs(prop)**2*gsq*V/xn
c--- include factor for W hadronic decays
      if (plabel(3) == 'pp') fac=2._dp*xn*fac
      if (plabel(7) == 'pp') fac=2._dp*xn*fac
      qqb=0._dp
      do hb=1,2
      do hg=1,2
      do hc=1,2
      do h12=1,2
      qqb=qqb+fac*aveqq*abs(mtot(hb,hg,hc,h12))**2
      enddo
      enddo
      enddo
      enddo
      
c--- gg amplitudes
      do hb=1,2
      do hg=1,2
      do hc=1,2
      do h1=1,2
      do h2=1,2
      mabtot(hb,hg,h1,h2,hc)=czip
      mbatot(hb,hg,h1,h2,hc)=czip

      do j1=1,j1max
      do j2=1,2
      mabtot(hb,hg,h1,h2,hc)=mabtot(hb,hg,h1,h2,hc)+
     & mtop(hb,j1)*mab(j1,h1,h2,j2)*manti(j2,hg,hc)
      mbatot(hb,hg,h1,h2,hc)=mbatot(hb,hg,h1,h2,hc)+
     & mtop(hb,j1)*mba(j1,h1,h2,j2)*manti(j2,hg,hc)
      mqed(hb,hg,h1,h2,hc)=mabtot(hb,hg,h1,h2,hc)+mbatot(hb,hg,h1,h2,hc)
      enddo
      enddo

      enddo
      enddo
      enddo
      enddo
      enddo

      gg=0._dp
      do hb=1,2
      do hg=1,2
      do hc=1,2
      do h1=1,2
      do h2=1,2
      gg=gg+fac*avegg*xn
     & *(abs(mabtot(hb,hg,h1,h2,hc))**2+abs(mbatot(hb,hg,h1,h2,hc))**2
     & -abs(mqed(hb,hg,h1,h2,hc))**2/xnsq)
      enddo
      enddo
      enddo
      enddo
      enddo
       
C---fill qb-q, gg and q-qb elements
      do j=-nf,nf
      if ((j < 0) .or. (j > 0)) then
          msq(j,-j)=qqb
      elseif (j == 0) then
          msq(0,0)=gg
      endif
      enddo

      return
      end
