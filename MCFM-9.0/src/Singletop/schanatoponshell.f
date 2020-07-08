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
 
      subroutine schanatoponshell(q1,q2,p,iswitch,m)
      implicit none
      include 'types.f'
      
************************************************************************
*     Author: R.K. Ellis, January 2012                                 *
*                                                                      *
*     d(-q1)+u~(-q2)--> t~(ee,nb,pc)+pb(p6)                            *
*                                                                      *
*     keeping polarization information for t~ named "a"                *
*     iswitch= 0 for no gluon emission                                 *
*     iswitch=+1 for gluon emission in top decay                       *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'zprods_decl.f'
      real(dp):: p(mxpart,4),q(mxpart,4),
     & ala,alb,dot,s12,mt2
      complex(dp):: m(2,2),izb,cprop
      integer:: p1,p2,a,ee,b,si,i,j,iswitch,q1,q2
      parameter(p1=1,p2=2,ee=3,a=4,b=5)
C-----matrix element for .e+_dpu~ -> t~+b where t~ and b are on shell
C-----a rendered massless wrt ee, and pb rendered massless wrt p2
c--- statement functions
      izb(i,j)=cone/zb(i,j)
c--- end statement functions
C---zero all arrays
      do i=1,2
      do j=1,2
      m(i,j)=czip
      enddo
      enddo
      do si=1,4
      q(p1,si)=p(q1,si)
      q(p2,si)=p(q2,si)
      if (iswitch == 0) then
      q(a,si)=p(3,si)+p(4,si)+p(5,si)
      elseif (iswitch == 1) then
      q(a,si)=p(3,si)+p(4,si)+p(5,si)+p(7,si)
      endif
      q(ee,si)=p(3,si)
      q(b,si)=p(6,si)
      enddo
      mt2=mt**2
      s12=2._dp*dot(q,p1,p2)
      cprop=cplx2(s12-wmass**2,wmass*wwidth)
      
C---- now render "a" massless wrt to vector ee
C---- now render "pb" massless wrt to vector p2
      ala=mt2/(2._dp*dot(q,a,ee))
      alb=mb**2/(2._dp*dot(q,b,p2))
      do si=1,4
      q(a,si)=q(a,si)-ala*q(ee,si)
      q(b,si)=q(b,si)-alb*q(p2,si)
      enddo
      call spinoru(5,q,za,zb)
      
C----order of indices is polb,pola
      m(1,2)= - za(b,p2)*zb(a,p1)*cprop**(-1)

      m(2,2)=czip

      
      m(1,1)=za(b,p2)*zb(ee,p1)*izb(a,ee)*cprop**(-1)*mt

      m(2,1)=czip

      return
      end
