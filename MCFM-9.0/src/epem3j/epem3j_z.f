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
 
      subroutine epem3j_z(p,z)
      implicit none
      include 'types.f'
************************************************************************
*     John M. Campbell                                                 *
*     November, 2008.                                                  *
************************************************************************
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'PR_new.f'
      include 'agq.f'
      include 'colstruc.f'
      integer:: is
      real(dp):: z,xl34,xl35,xl45,p(mxpart,4),dot
      real(dp):: ff_qq,ff_gg

      xl34=log(+two*dot(p,3,4)/musq)
      xl35=log(+two*dot(p,3,5)/musq)
      xl45=log(+two*dot(p,4,5)/musq)

      do is=1,3
c--- leading colour (Ca) goes in (0,0)
c--- subleading colour (2*Cf-Ca) goes in (0,1)
c--- Tr terms go in (1,0)
      caonly=.true.
      nfonly=.false.
      Q1(g,g,g,is)=ason4pi*ca*(ff_qq(z,xl45,is)+half*ff_gg(z,xl45,is)
     &                        +ff_qq(z,xl35,is)+half*ff_gg(z,xl35,is))

      Q1(g,g,q,is)=ason4pi*ca*(-ff_qq(z,xl34,is)-ff_qq(z,xl34,is))/xnsq

      caonly=.false.
      nfonly=.true.
      Q1(q,q,g,is)=ason4pi*ca*(half*ff_gg(z,xl45,is)
     &                        +half*ff_gg(z,xl35,is))
      enddo
      
      return
      end
