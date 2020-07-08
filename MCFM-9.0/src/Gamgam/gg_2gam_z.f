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
 
      subroutine gg_2gam_z(p,z)
      implicit none
      include 'types.f'
************************************************************************
*     John M. Campbell                                                 *
*     February, 2016.                                                  *
************************************************************************
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      include 'qcdcouple.f'
      include 'facscale.f'
      include 'scale.f'
      include 'PR_new.f'
      include 'agq.f'
      integer:: is
      real(dp):: z,xl12,p(mxpart,4),dot,ii_gg,ii_gq,tempgg,tempgq

      xl12=log(two*dot(p,1,2)/musq)
      
      do is=1,3
      tempgg=ason2pi*xn*ii_gg(z,xl12,is)
      Q1(g,g,g,is)=tempgg
      Q2(g,g,g,is)=tempgg

      tempgq=ason4pi*two*cf*ii_gq(z,xl12,is)
      Q1(g,q,g,is)=tempgq
      Q1(g,a,g,is)=tempgq
      Q2(g,q,g,is)=tempgq
      Q2(g,a,g,is)=tempgq
      enddo

      return
      end
