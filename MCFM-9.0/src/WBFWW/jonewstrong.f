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
 
      subroutine jonewstrong(p7,p3,p4,p1,za,zb,zab,jqcd)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'cmplxmass.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: p1,p2,p3,p4,p7,ro
      complex(dp):: zab(mxpart,4,mxpart),jqcd(2,4),propw34
      real(dp):: t3,s34,s134,s347
C-----Begin statement functions
      t3(p1,p2,p3)=s(p1,p2)+s(p2,p3)+s(p3,p1)
C-----end statement functions

      s34=s(p3,p4)
      s134=t3(p1,p3,p4)
      s347=t3(p3,p4,p7)

      propw34=s34-cwmass2

      do ro=1,4
      jqcd(1,ro)= + propw34**(-1) * ( za(p7,p3)*zb(p7,p4)*zab(p7,ro,p1)
     &    *s347**(-1) + za(p7,p3)*zb(p3,p4)*zab(p3,ro,p1)*s347**(-1) - 
     &    za(p1,p3)*zb(p1,p4)*zab(p7,ro,p1)*s134**(-1) + za(p3,p4)*zb(
     &    p1,p4)*zab(p7,ro,p4)*s134**(-1) )
      jqcd(2,ro)=jqcd(1,ro)
      enddo
      jqcd(:,:)=jqcd(:,:)/cxw
      return
      end
