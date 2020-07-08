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
 
      subroutine jonewWpWp(p7,p3,p4,p1,za,zb,zab,jZ,jg)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cmplxmass.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'zcouple.f'
      include 'nwz.f'
      integer:: p1,p2,p3,p4,p7,ro
      complex(dp):: zab(mxpart,4,mxpart),zab2,jZ(2,4),jg(2,4),
     & rxw,propw34,propw17
      real(dp):: qn,p17(4),p34(4),t3,s134,s347,
     & s34,s137,s147,s17,q3,q4,l3,l4
C     amp(jdu1,jdu2,h17,h28)
      parameter(qn=0._dp)
C-----Begin statement functions
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      t3(p1,p2,p3)=s(p1,p2)+s(p2,p3)+s(p3,p1)
C-----end statement functions

      rxw=sqrt((cone-cxw)/cxw)
      s17=s(p1,p7)
      s34=s(p3,p4)
      s134=t3(p1,p3,p4)
      s137=t3(p1,p3,p7)
      s147=t3(p1,p4,p7)
      s347=t3(p3,p4,p7)
      p17(:)=0.5_dp*real(zab(p1,:,p1)+zab(p7,:,p7))
      p34(:)=0.5_dp*real(zab(p3,:,p3)+zab(p4,:,p4))

      propw17=s17-cwmass2
      propw34=s34-cwmass2

c--- determine correct couplings from nwz
      if (nwz == +1) then
        q3=qn
        l3=ln
        q4=qe
        l4=le
      elseif (nwz == -1) then
        q3=qe
        l3=le
        q4=qn
        l4=ln
      else
        write(6,*) 'Unexpected nwz in jonewWpWp.f: ',nwz
        stop
      endif

      do ro=1,4
      jZ(1,ro)= + propw34**(-1) * ( L(1)*za(p7,p3)*zb(p7,p4)*zab(p7,ro,
     &    p1)*s347**(-1) + L(1)*za(p7,p3)*zb(p3,p4)*zab(p3,ro,p1)*
     &    s347**(-1) - L(2)*za(p1,p3)*zb(p1,p4)*zab(p7,ro,p1)*
     &    s134**(-1) + L(2)*za(p3,p4)*zb(p1,p4)*zab(p7,ro,p4)*
     &    s134**(-1) + za(p7,p3)*zb(p1,p4)*p17(ro)*propw17**(-1)*rxw - 
     &    za(p7,p3)*zb(p1,p4)*p34(ro)*propw17**(-1)*rxw - zab2(p7,p3,p4
     &    ,p1)*zab(p3,ro,p4)*propw17**(-1)*rxw + zab2(p3,p1,p7,p4)*zab(
     &    p7,ro,p1)*propw17**(-1)*rxw )
      jZ(1,ro) = jZ(1,ro) - za(p7,p1)*zb(p1,p4)*zab(p3,ro,p1)*
     & propw17**(-1)*l3*s147**(-1) - za(p7,p3)*zb(p7,p1)*zab(p7,ro,p4)*
     &    propw17**(-1)*l4*s137**(-1) + za(p7,p3)*zb(p1,p3)*zab(p3,ro,
     &    p4)*propw17**(-1)*l4*s137**(-1) - za(p7,p4)*zb(p1,p4)*zab(p3,
     &    ro,p4)*propw17**(-1)*l3*s147**(-1)
      jZ(2,ro)= + propw34**(-1) * (  - L(1)*za(p1,p3)*zb(p1,p4)*zab(p7,
     &    ro,p1)*s134**(-1) + L(1)*za(p3,p4)*zb(p1,p4)*zab(p7,ro,p4)*
     &    s134**(-1) + L(2)*za(p7,p3)*zb(p7,p4)*zab(p7,ro,p1)*
     &    s347**(-1) + L(2)*za(p7,p3)*zb(p3,p4)*zab(p3,ro,p1)*
     &    s347**(-1) - za(p7,p3)*zb(p1,p4)*p17(ro)*propw17**(-1)*rxw + 
     &    za(p7,p3)*zb(p1,p4)*p34(ro)*propw17**(-1)*rxw + zab2(p7,p3,p4
     &    ,p1)*zab(p3,ro,p4)*propw17**(-1)*rxw - zab2(p3,p1,p7,p4)*zab(
     &    p7,ro,p1)*propw17**(-1)*rxw )
      jZ(2,ro) = jZ(2,ro) - za(p7,p1)*zb(p1,p4)*zab(p3,ro,p1)*
     & propw17**(-1)*l3*s147**(-1) - za(p7,p3)*zb(p7,p1)*zab(p7,ro,p4)*
     &    propw17**(-1)*l4*s137**(-1) + za(p7,p3)*zb(p1,p3)*zab(p3,ro,
     &    p4)*propw17**(-1)*l4*s137**(-1) - za(p7,p4)*zb(p1,p4)*zab(p3,
     &    ro,p4)*propw17**(-1)*l3*s147**(-1)
      jg(1,ro)= + propw34**(-1) * ( za(p7,p3)*zb(p7,p4)*zab(p7,ro,p1)*
     &    Qd*s347**(-1) + za(p7,p3)*zb(p1,p4)*p17(ro)*propw17**(-1) - 
     &    za(p7,p3)*zb(p1,p4)*p34(ro)*propw17**(-1) + za(p7,p3)*zb(p3,
     &    p4)*zab(p3,ro,p1)*Qd*s347**(-1) - za(p1,p3)*zb(p1,p4)*zab(p7,
     &    ro,p1)*Qu*s134**(-1) + za(p3,p4)*zb(p1,p4)*zab(p7,ro,p4)*Qu*
     &    s134**(-1) - zab2(p7,p3,p4,p1)*zab(p3,ro,p4)*propw17**(-1) + 
     &    zab2(p3,p1,p7,p4)*zab(p7,ro,p1)*propw17**(-1) )
      jg(1,ro) = jg(1,ro) - za(p7,p1)*zb(p1,p4)*zab(p3,ro,p1)*
     & propw17**(-1)*q3*s147**(-1) - za(p7,p3)*zb(p7,p1)*zab(p7,ro,p4)*
     &    propw17**(-1)*q4*s137**(-1) + za(p7,p3)*zb(p1,p3)*zab(p3,ro,
     &    p4)*propw17**(-1)*q4*s137**(-1) - za(p7,p4)*zb(p1,p4)*zab(p3,
     &    ro,p4)*propw17**(-1)*q3*s147**(-1)
      jg(2,ro)= + propw34**(-1) * ( za(p7,p3)*zb(p7,p4)*zab(p7,ro,p1)*
     &    Qu*s347**(-1) - za(p7,p3)*zb(p1,p4)*p17(ro)*propw17**(-1) + 
     &    za(p7,p3)*zb(p1,p4)*p34(ro)*propw17**(-1) + za(p7,p3)*zb(p3,
     &    p4)*zab(p3,ro,p1)*Qu*s347**(-1) - za(p1,p3)*zb(p1,p4)*zab(p7,
     &    ro,p1)*Qd*s134**(-1) + za(p3,p4)*zb(p1,p4)*zab(p7,ro,p4)*Qd*
     &    s134**(-1) + zab2(p7,p3,p4,p1)*zab(p3,ro,p4)*propw17**(-1) - 
     &    zab2(p3,p1,p7,p4)*zab(p7,ro,p1)*propw17**(-1) )
      jg(2,ro) = jg(2,ro) - za(p7,p1)*zb(p1,p4)*zab(p3,ro,p1)*
     & propw17**(-1)*q3*s147**(-1) - za(p7,p3)*zb(p7,p1)*zab(p7,ro,p4)*
     &    propw17**(-1)*q4*s137**(-1) + za(p7,p3)*zb(p1,p3)*zab(p3,ro,
     &    p4)*propw17**(-1)*q4*s137**(-1) - za(p7,p4)*zb(p1,p4)*zab(p3,
     &    ro,p4)*propw17**(-1)*q3*s147**(-1)
      enddo
      jZ(:,:)=jZ(:,:)/cxw
      jg(:,:)=jg(:,:)/cxw
      return
      end
