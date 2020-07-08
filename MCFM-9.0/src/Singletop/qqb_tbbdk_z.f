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
 
      subroutine qqb_tbbdk_z(p,z)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'PR_new.f'
      include 'agq.f'
      integer:: is
      real(dp):: z,xl12,p(mxpart,4),dot,ii_qq,ii_qg,tempqq,tempqg,
     & mbar12,mbar12a,mbar12b,ff_1mqq,ff_2mqq

      xl12=log(two*dot(p,1,2)/musq)

c---- note: the q-qbar integrated counterterms corresponding to final-final
c----       dipoles depend on whether the b~ quark is massless or not
      do is=1,3
      
      if (abs(mb) < 1.e-6_dp) then
c--- integrated dipole for mb=0
        mbar12=mt/sqrt(two*dot(p,1,2))
        tempqq=+ason2pi*cf*(ii_qq(z,xl12,is)+ff_1mqq(z,xl12,mbar12,is))
      else
c--- integrated dipole for mb>0
        mbar12a=mt/sqrt(two*dot(p,1,2))
        mbar12b=mb/sqrt(two*dot(p,1,2))
        tempqq=+ason2pi*cf*(ii_qq(z,xl12,is)
     &                     +ff_2mqq(z,xl12,mbar12a,mbar12b,is)/2._dp
     &                     +ff_2mqq(z,xl12,mbar12b,mbar12a,is)/2._dp)
      endif
      tempqg=+ason2pi*tr*ii_qg(z,xl12,is)

      Q1(q,q,a,is)=tempqq
      Q2(a,a,q,is)=tempqq
      Q1(a,a,q,is)=tempqq
      Q2(q,q,a,is)=tempqq

      Q2(a,g,q,is)=tempqg
      Q2(q,g,a,is)=tempqg
      Q1(a,g,q,is)=tempqg
      Q1(q,g,a,is)=tempqg
      enddo
      
      return
      end
