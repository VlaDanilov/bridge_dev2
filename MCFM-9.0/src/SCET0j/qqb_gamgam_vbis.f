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
 
      subroutine qqb_gamgam_vbis(p,msq,order) 
!------ routine for qqb_gamgam using hard coeff
      implicit none 
      include 'types.f' 
      include 'constants.f' 
      include 'mxpart.f'
      include 'nf.f' 
      include 'qcdcouple.f'
      include 'ewcharge.f'
      include 'ewcouple.f' 
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf)
      integer:: order
      real(dp) :: hard(2),qqb,qbq
      real(dp)::statfac,fac
      integer:: i,j,k
      parameter(statfac=0.5_dp)

      fac=8._dp*xn*esq**2*statfac

!======= compute Matrix elements 
      call gamgamampsq(order,p,1,2,3,4,qqb,hard)

c---- apply overall factor 
      qqb=fac*aveqq*qqb
      qbq=qqb

!----- higher order corrections
      if(order==1) then 
         qqb=qqb*ason2pi*hard(1)
         qbq=qbq*ason2pi*hard(1)
      endif

      msq(:,:)=0._dp
      do j=-nf,nf
         k=-j
         if (j*k >= 0) cycle    ! skip gluons, qq, aa
      
         if ((j > 0) .and. (k < 0)) then
            msq(j,k)=qqb*Q(j)**4
         elseif ((j < 0) .and. (k > 0)) then
            msq(j,k)=qbq*Q(k)**4
         endif
         
         
      enddo
      return
      end
