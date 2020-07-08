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
 
      subroutine gg_2gam_v(p,msq)
      implicit none
      include 'types.f'      
************************************************************************
*     Authors: R.K. Ellis and John M. Campbell                         *
*     December, 2010.                                                  *
************************************************************************
*                                                                      *
*     Matrix element for gamma + gamma production,                     *
*     averaged over initial colours and spins                          *
*                                                                      *
*     g(-p1)+g(-p2) --> gamma(p3)+gamma(p4)                            *
*                                                                      *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'scheme.f'
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),facgg,Qsum,
     & virtgamgam
      real(dp),parameter::statfac=0.5_dp

c--set msq=0 to initalize
      msq(:,:)=0._dp

!      scheme='dred'
! This bug pointed out by Xiaran Zhao, Aug 21 2018
      scheme='tH-V'

      call dotem(3,p,s)

c--- initialize gg 2-loop matrix elements
      Qsum=Q(1)**2+Q(2)**2+Q(3)**2+Q(4)**2+Q(5)**2 
      facgg=4._dp*esq*gsq/(16._dp*pisq)*Qsum

      msq(0,0)=avegg*V*facgg**2*statfac*virtgamgam(s(1,2),s(1,3),s(2,3))

      return
      end


