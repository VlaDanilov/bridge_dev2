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
 
      subroutine schantoponshellv(q1,q2,p,m,mv)
      implicit none
      include 'types.f'
      
************************************************************************
*     Author: R.K. Ellis, January 2012                                 *
*                                                                      *
*     Virtual corrections to:-                                         *
*     u(-p1)+d~(-p2)--> t(nu,eb,pb)+pc~(p6)                            *
*                                                                      *
*     keeping polarization information for t                           *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'zprods_decl.f'
      real(dp):: p(mxpart,4),q(mxpart,4),
     & alt,alc,dot,s12,mt2,ct
      complex(dp):: m(2,2),mv(2,2),iza,izb,wprop
      complex(dp):: X0L,C0L,C0R,C1L,C1R,C1Lon2,C1Ron2
      integer:: p1,p2,t,eb,c,si,i,j,q1,q2
      logical:: oldincludect
      parameter(p1=1,p2=2,eb=4,c=5,t=3)
C-----matrix element for .e+_dpu~ -> t+b~ where both t and b~ are on shell
C-----t rendered massless wrt eb, and b rendered massless wrt p1
c--- statement functions
      iza(i,j)=cone/za(i,j)
      izb(i,j)=cone/zb(i,j)
c--- end statement functions
     
c--- corrections in production do not need CT to be included     

C---zero all arrays
      do i=1,2
      do j=1,2
      m(i,j)=czip
      mv(i,j)=czip
      enddo
      enddo
      do si=1,4
      q(p1,si)=p(q1,si)
      q(p2,si)=p(q2,si)
      q(eb,si)=p(4,si)
      q(c,si)=p(6,si)
      q(t,si)=p(3,si)+p(4,si)+p(5,si)
      enddo
      mt2=mt**2
      s12=2._dp*dot(q,p1,p2)
c--- corrections in production do not need CT to be included     
      call coefsdkmass(.false.,s12,0._dp,0._dp,ct,X0L,c0R,c1L,C1R)
      call coefsdkmass(.false.,s12,mt,mb,ct,C0L,c0R,c1L,C1R)
      C1Lon2=0.5_dp*c1L
      C1Ron2=0.5_dp*c1R
      wprop=cplx2(s12-wmass**2,wmass*wwidth)
     
C---- now render "t" massless wrt to vector eb
C---- now render "pc" massless wrt to vector p1
      alt=mt2/(2._dp*dot(q,t,eb))
      alc=mb**2/(2._dp*dot(q,c,p1))
      do si=1,4
      q(t,si)=q(t,si)-alt*q(eb,si)
      q(c,si)=q(c,si)-alc*q(p1,si)
      enddo
      call spinoru(5,q,za,zb)
      
C----order of indices is polt,polc
      m(1,2)= - za(t,p2)*zb(c,p1)*wprop**(-1)

      m(1,1)=czip

      m(2,2)=za(eb,p2)*zb(c,p1)*iza(t,eb)*wprop**(-1)*mt

      m(2,1)=czip

      mv(1,2)= - za(c,p2)*za(t,p1)*zb(c,p1)*iza(c,p1)*c1Lon2*
     & wprop**(-1)*mt**(-1)*mb + za(c,p2)*zb(c,eb)*zb(c,p1)*izb(t,eb)*
     & c1Ron2*wprop**(-1) - za(t,p2)*zb(c,p1)*x0L*wprop**(-1) - za(t,p2
     & )*zb(c,p1)*c0L*wprop**(-1) - za(p1,p2)*zb(eb,p1)*iza(c,p1)*izb(t
     & ,eb)*c0R*wprop**(-1)*mt*mb

      mv(1,1)= - za(c,t)*za(c,p2)*zb(c,p1)*c1Lon2*wprop**(-1)*mt**(-1)
     &  + za(c,p2)*zb(eb,p1)*izb(t,eb)*c1Ron2*wprop**(-1)*mb + za(c,p2)
     & *zb(eb,p1)*izb(t,eb)*c0R*wprop**(-1)*mt

      mv(2,2)=za(c,p2)*za(eb,p1)*zb(c,p1)*iza(c,p1)*iza(t,eb)*c1Lon2*
     & wprop**(-1)*mb - za(c,p2)*zb(c,t)*zb(c,p1)*c1Ron2*wprop**(-1)*
     & mt**(-1) + za(eb,p2)*zb(c,p1)*iza(t,eb)*x0L*wprop**(-1)*mt + za(
     & eb,p2)*zb(c,p1)*iza(t,eb)*c0L*wprop**(-1)*mt + za(p1,p2)*zb(t,p1
     & )*iza(c,p1)*c0R*wprop**(-1)*mb

      mv(2,1)=za(c,eb)*za(c,p2)*zb(c,p1)*iza(t,eb)*c1Lon2*wprop**(-1)
     &  - za(c,p2)*zb(t,p1)*c1Ron2*wprop**(-1)*mt**(-1)*mb - za(c,p2)*
     & zb(t,p1)*c0R*wprop**(-1)

      return
      end
