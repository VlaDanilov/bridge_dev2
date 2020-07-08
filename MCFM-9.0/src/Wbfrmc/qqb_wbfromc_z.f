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
 
      subroutine qqb_wbfromc_z(p,z)
      implicit none
      include 'types.f'
************************************************************************
*     Author: J.M. Campbell                                            *
*     November, 2011.                                                  *
************************************************************************
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'PR_new.f'
      include 'agq.f'
      include 'breit.f'
      integer:: is
      real(dp):: z,xl12,xl15,xl25,p(mxpart,4),dot
      real(dp):: mbar12,mbar15,mbar25
      real(dp):: ii_mqq,ii_mqg,ii_mgq,ii_mgg,
     &                 if_mqq,if_mgg,
     &                 fi_mqq

      xl12=log(+two*dot(p,1,2)/musq)
      xl15=log(-two*dot(p,1,5)/musq)
      xl25=log(-two*dot(p,2,5)/musq)
      mbar12=0._dp
CDTS (5.45,5.77)
      mbar15=mass2/sqrt(-2._dp*dot(p,1,5))
      mbar25=mass2/sqrt(-2._dp*dot(p,2,5))
c--- sum over regular and plus terms
      do is=1,3
c--- (q,g)
      Q2(g,g,q,is)=ason4pi*xn
     & *(ii_mgg(z,xl12,mbar12,is)+if_mgg(z,xl25,mbar25,is)
     &  +fi_mqq(z,xl25,mbar25,is))
      Q1(q,q,g,is)=ason4pi*xn
     & *(ii_mqq(z,xl12,mbar12,is)
     & -(if_mqq(z,xl15,mbar15,is)+fi_mqq(z,xl15,mbar15,is))/xnsq)
c--- (qb,g)
      Q2(g,g,a,is)=Q2(g,g,q,is)
      Q1(a,a,g,is)=Q1(q,q,g,is)

c--- (g,q)
      Q1(g,g,q,is)=ason4pi*xn
     & *(ii_mgg(z,xl12,mbar12,is)+if_mgg(z,xl15,mbar15,is)
     &  +fi_mqq(z,xl15,mbar15,is))
      Q2(q,q,g,is)=ason4pi*xn
     & *(ii_mqq(z,xl12,mbar12,is)
     & -(if_mqq(z,xl25,mbar25,is)+fi_mqq(z,xl25,mbar25,is))/xnsq)
c--- (g,qb)
      Q1(g,g,a,is)=Q1(g,g,q,is)
      Q2(a,a,g,is)=Q2(q,q,g,is)

c--- (g,g)
      Q1(q,g,g,is)=ason4pi*2._dp*tr*ii_mqg(z,xl12,mbar12,is)
      Q2(q,g,g,is)=Q1(q,g,g,is)
 
c--- these are terms coming from the diagrams which are
c--- a 4-quark contribution, with (q,q~)->g in the initial state
      Q2(g,q,a,is)=ason4pi*2._dp*cf*ii_mgq(z,xl12,mbar12,is)
      Q2(g,a,q,is)=Q2(g,q,a,is)
      Q2(g,a,a,is)=Q2(g,q,a,is)
      Q2(g,q,q,is)=Q2(g,q,a,is)

      Q1(g,q,a,is)=Q2(g,q,a,is)
      Q1(g,a,q,is)=Q2(g,q,a,is)
      Q1(g,a,a,is)=Q2(g,q,a,is)
      Q1(g,q,q,is)=Q2(g,q,a,is)
      
      enddo
      
      return
      end
