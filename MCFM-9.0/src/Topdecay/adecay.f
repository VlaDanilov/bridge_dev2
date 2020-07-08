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
 
      subroutine adecay(p,pe,pnb,pc,m)
      implicit none
      include 'types.f'
      
************************************************************************
*     Author: R.K. Ellis, January 2012                                 *
*     antitop decay  a --> e(pe)+nb(pnb)+bbar(pc)                      *
*     with bottom and top masses (and no radiation)                    *
*     in massless spinor notation                                      *
*     pe,pnb,pc are integer::s that point to                             *
*     the appropriate four-momenta in p                                *
*     pe=electron                                                      *
*     pnb=antineutrino                                                 *
*     pc=anti-bottom quark                                             *
*     q(c) is rendered massless wrt to pnb                             *
*     q(a) is rendered massless wrt to pe                              *
*     returned m(apol,cpol)                                            *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      real(dp):: p(mxpart,4),q(mxpart,4),dot,sw,ala,alc
      complex(dp):: m(2,2),cprop
      integer:: nb,a,c,e,si,pc,pe,pnb
      parameter(a=1,e=3,nb=4,c=2)
      do si=1,4
      q(a,si)=p(pe,si)+p(pnb,si)+p(pc,si)
      q(e,si)=p(pe,si)
      q(nb,si)=p(pnb,si)
      q(c,si)=p(pc,si)
      enddo
      ala=mt**2/(2._dp*dot(q,a,e))
      alc=mb**2/(2._dp*dot(q,c,nb))
      do si=1,4
      q(a,si)=q(a,si)-ala*q(e,si)
      q(c,si)=q(c,si)-alc*q(nb,si)
      enddo

      call spinoru(4,q,za,zb)      
      sw=s(e,nb)
      cprop=cplx2(sw-wmass**2,wmass*wwidth)
C---order of polarizations is the m(apol,cpol)
      m(1,1)=czip
      m(1,2)=czip
      m(2,1)=czip
      m(2,2)=-za(e,a)*zb(c,nb)/cprop
      return
      end
