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
 
      subroutine gg_2gam_g(p,msq)
      implicit none
      include 'types.f'
      
************************************************************************
*     Authors: R.K. Ellis and John M. Campbell                         *
*     December, 2010.                                                  *
************************************************************************
*                                                                      *
*     Matrix element for gamma + gamma + 1 parton production           *
*     in order alpha_s, averaged over initial colours and spins        *
*                                                                      *
*     g(-p1)+g(-p2) --> gamma(p3)+gamma(p4)+gluon(p5)                  *
*                                                                      *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_com.f'
      include 'noglue.f'
      integer j
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),ggtogagag,qqbgnfsq,
     & qa_gagag,aq_gagag,ag_gagaa,qg_gagaq,ga_gagaa,gq_gagaq
      
c--set msq=0 to initalize
      msq(:,:)=0._dp
      
      call spinoru(5,p,za,zb)

      qg_gagaq=aveqg*qqbgnfsq(1,5,2,3,4)
      gq_gagaq=aveqg*qqbgnfsq(2,5,1,3,4)
      qa_gagag=aveqq*qqbgnfsq(1,2,5,3,4)

      ag_gagaa=qg_gagaq
      ga_gagaa=gq_gagaq
      aq_gagag=qa_gagag
!      ag_gagaa=-aveqg*qqbgnfsq(5,1,2,3,4)
!      ga_gagaa=-aveqg*qqbgnfsq(5,2,1,3,4)
!      aq_gagag=+aveqq*qqbgnfsq(2,1,5,3,4)

      do j=1,nf
      msq(+j,0)=qg_gagaq
      msq(-j,0)=ag_gagaa

      msq(0,+j)=gq_gagaq
      msq(0,-j)=ga_gagaa

      msq(j,-j)=qa_gagag
      msq(-j,j)=aq_gagag
      enddo
      
c--- averaging factor already included in ggtogagag 
      msq(0,0)=ggtogagag()
      
      return
      end
