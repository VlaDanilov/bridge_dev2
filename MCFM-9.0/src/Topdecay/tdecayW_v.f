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
 
c--- File written by FORM program tdecayW_v.frm on Thu Mar  1 12:29:11 CST 2012
      subroutine tdecayW_v(p,pqq,pqb,pb,m)
      implicit none
      include 'types.f'
      
************************************************************************
*     Author: R.K. Ellis, January 2012                                 *
*     Virtual corrections to                                           *
*     top decay  t --> q(pqq)+qb(pqb)+b(pb)                            *
*     with bottom and top masses (and no radiation)                    *
*     and with overall factor of (4 pi)^e/Gamma(1-e)                   *
*     in massless spinor notation                                      *
*     pe,pnb,pc are integer::s that point to                             *
*     the appropriate four-momenta in p                                *
*     pqq=quark                                                        *
*     pqb=antiquark                                                    *
*     pb=bottom quark                                                  *
*     q(t) is rendered massless wrt to pqb                             *
*     q(b) is rendered massless wrt to pqq                             *
*     returned m(bpol,tpol)                                            *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      real(dp):: p(mxpart,4),q(mxpart,4),dot,sw,alt,alb,c0,ct
      complex(dp):: m(2,2),cprop,c0L
      integer:: qq,qb,b,t,si,pb,pqq,pqb
      parameter(t=1,qq=3,qb=4,b=2)
      do si=1,4
      q(t,si)=p(pqq,si)+p(pqb,si)+p(pb,si)
      q(qq,si)=p(pqq,si)
      q(qb,si)=p(pqb,si)
      q(b,si)=p(pb,si)
      enddo
      alt=mt**2/(2._dp*dot(q,t,qb))
      alb=mb**2/(2._dp*dot(q,b,qq))
      do si=1,4
      q(t,si)=q(t,si)-alt*q(qb,si)
      q(b,si)=q(b,si)-alb*q(qq,si)
      enddo
      call spinoru(4,q,za,zb)
      sw=s(qq,qb)
      call coefswdk(sw,ct,c0)
      c0L=cplx1(c0+ct)
      cprop=cplx2(sw-wmass**2,wmass*wwidth)
C---order of polarizations is m(bpol,tpol)
      m(1,1)= + cprop**(-1) * (  - za(qq,b)*zb(qb,t)*c0L )

      m(1,2)= czip

      m(2,1)= czip

      m(2,2)= czip

      return
      end
