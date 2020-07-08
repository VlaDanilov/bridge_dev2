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
 
      subroutine qqb_QQb_mix_v(p,msq)
c--- Exact electroweak corrections to heavy quark pair production,
c--- including mixed QCD-weak LO process
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'constants.f'
      include 'breit.f'
      include 'sprods_com.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'scheme.f'
!      include 'first.f'
      include 'noglue.f'
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),ss,beta,z,
     &     qqbew(5),qbqew(5),ggew1,ggew2,m2(-nf:nf,-nf:nf)
      integer j,k
      
      scheme = 'dred'

!      if(first) then
!         first=.false.
!         write(6,*) 'Heavy Quark mass:', mass2
!      endif

      call dotem(4,p,s)

      ss = s(1,2)
      beta = sqrt(1._dp-4._dp*mt**2/ss)
      z = (s(1,3)-s(2,3))/ss/beta

c--- avoid calculating unnecessarily by checking input file flags
      if (ggonly) then
        qqbew(:)=zip
      else
        call qqbQQb_ew_oneloop(qqbew,ss,beta,z)
        call qqbQQb_ew_oneloop(qbqew,ss,beta,-z)
      endif
      
      if (noglue .or. omitgg) then
        ggew1=zip
        ggew2=zip
      else
        call ggQQb_ew_oneloop(ggew1,ss,beta,z)
        call ggQQb_ew_oneloop(ggew2,ss,beta,-z)
      endif
      
      call qqb_QQb(p,m2)

      msq = 0._dp

      do j=-nf,nf
         k=-j
         if((j == 0) .and. (k == 0)) then
            msq(j,k) = ggew1 + ggew2
         else if((j > 0) .and. (k < 0)) then
            msq(j,k) = qqbew(j)
         else if((j < 0 ) .and. (k > 0)) then
            msq(j,k) = qbqew(k)
         end if
      end do

      msq = msq*m2!/gsq**2*16._dp*pi**2!*0.01_dp

      end subroutine qqb_QQb_mix_v
