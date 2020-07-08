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
 
!-----q(i1) + qb(i2)  +  + g(i3)+gamma(i4)+gamma(i5) nf loops squared
!=====extracted from gmgmjt routines
!=====CW Feb 2016
      function qqbgnfsq(i1,i2,i3,i4,i5)
      implicit none
      include 'types.f'
      real(dp):: qqbgnfsq
      integer i1,i2,i3,i4,i5
      include 'constants.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'nf.f'
      include 'ewcharge.f'
      include 'zprods_com.f'
      complex(dp) qqbg_nf(2,2,2,2)
      real(dp) nf_fac,fac,faclo,Qsum
      real(dp),parameter::statfac=0.5_dp
      integer h1,h2,h3,h4

      Qsum=Q(1)**2+Q(2)**2+Q(3)**2+Q(4)**2+Q(5)**2
      faclo=8._dp*cf*xn*gsq*esq**2*statfac
      fac=faclo*(xn*ason2pi/2._dp)**2
      nf_fac=fac/xn**2

      call amp_virt_nf_gmgmjt(i1,i2,i3,i4,i5,za,zb,qqbg_nf)

      qqbgnfsq=zip
      do h1=1,2
        do h2=1,2
          do h3=1,2
            do h4=1,2
              qqbgnfsq=qqbgnfsq+real(qqbg_nf(h1,h2,h3,h4)
     &                        *conjg(qqbg_nf(h1,h2,h3,h4)),dp)
            enddo
          enddo
        enddo
      enddo
      
      qqbgnfsq=qqbgnfsq*nf_fac*Qsum**2
      
      return
      end
