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
 
      subroutine qqb_QQb_mix(p,msq)
c--- Leading order matrix element for mixed QCD/Z matrix element
c--- for top pair production, q(-p1)+qbar(-p2) -->  Q(p3)+Q~(P4) 
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'sprods_com.f'
      include 'breit.f'
      include 'first.f'
      integer j,k
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),gaq(nf),gvq(nf),
     .     gat,gvt,ss,beta,z,qqb(nf),qbq(nf),cs,sw2,cw2,mz,T3(nf)

      if(first) then
         first=.false.
         write(6,*) 'Heavy Quark mass:', mass2
      endif

      sw2 = xw
      cw2 = 1._dp - sw2
      T3 = (/-0.5_dp,0.5_dp,-0.5_dp,0.5_dp,-0.5_dp/)
      gvq(1:5) = (T3(1:5)-2._dp*sw2*Q(1:5))/2._dp/sqrt(sw2*cw2)
      gaq(1:5) = T3(1:5)/2._dp/sqrt(sw2*cw2)
      gw = 0.5_dp/sqrt(2._dp*sw2)

      gvt = gvq(2)
      gat = gaq(2)

      mz = zmass

      call dotem(4,p,s)

      ss = s(1,2)
      beta = sqrt(1._dp - 4._dp*mass2**2/ss)
      z = (s(1,3)-s(2,3))/ss/beta

      qqb = 4._dp*ss/(ss - mz**2)*(
     .     + (2._dp - beta**2*(1._dp - z**2))*gvq*gvt
     .     + 2._dp*beta*z*gaq*gat )

      qbq = 4._dp*ss/(ss - mz**2)*(
     .     + (2._dp - beta**2*(1._dp - z**2))*gvq*gvt
     .     - 2._dp*beta*z*gaq*gat )

      msq=zip
      do j=-nf,nf
         k=-j
         if((j == 0) .and. (k == 0)) then
            msq(j,k)=0._dp
         else if((j > 0) .and. (k < 0)) then
            msq(j,k) = qqb(j)
         else if((j < 0 ) .and. (k > 0)) then
            msq(j,k) = qbq(k)
         end if
      end do
 
! factor of two for the interference, 1/Nc**2/4 is the color and spin average
      msq = msq*2._dp*esq*gsq/Nc**2/4._dp
!      msq = 0._dp

      end subroutine qqb_QQb_mix
