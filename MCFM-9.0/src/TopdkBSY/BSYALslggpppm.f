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
 
      function BSYALslggpppm(e1,p2,e4,p3,za,zb,zab,zba)
      implicit none
      include 'types.f'
      complex(dp):: BSYALslggpppm
      
C-----Authors: John Campbell and Keith Ellis, November 2011
C---- arXiv:1101.5947 [hep-ph], Eq. (97)
C---- (These are twiddle functions, c.f.arXiv:1101.5947[hep-ph],Eq.(91))
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'zabprods_decl.f'
      include 'qdef.f'
      include 'massiveintegrals.f'
      real(dp):: s12,s23,mt2,mt3
      integer:: e1,p2,p3,e4,j
      j=p2-1
      mt2=mt**2
      mt3=mt**3
      s12=mt2+s(1,p2)
      s23=s(p2,p3)
      BSYALslggpppm =  + I3m12x3x4(j) * (  - 1.D0/2.D0/(zab(p2,q1,p3))*
     &    za(e1,p3)*za(e4,p3)*zb(p3,p2)*zab(p2,q1,p2)*mt + 1.D0/2.D0/(
     &    zab(p2,q1,p3))*za(e4,e1)*zab(p2,q1,p2)*zab(p3,q1,p2)*mt - 1.D0
     &    /2.D0/(zab(p2,q1,p3))/(za(p3,p2))*za(e1,p3)*za(e4,p3)*zab(p2,
     &    q1,p2)*mt3 + 1/(zab(p2,q1,p3))/(zab(p2,q1,p3))*za(e4,p3)*zab(
     &    p2,q1,p2)*zab(e1,q1,p3)*mt3 + 1/(zab(p2,q1,p3))/(zab(p2,q1,p3
     &    ))*za(e4,e1)*zab(p2,q1,p2)*s12*mt3 - 1/(zab(p2,q1,p3))/(zab(
     &    p2,q1,p3))/(za(p3,p2))*za(e1,p2)*za(e4,p3)*zab(p2,q1,p2)**2*
     &    mt3 + 1.D0/2.D0/(zab(p2,q1,p3))/(zab(p2,q1,p3))/(zab(p2,q1,p3
     &    ))/(za(p3,p2))*za(e1,p2)*za(e4,p2)*zab(p2,q1,p2)**2*zab(p3,q1
     &    ,p3)*mt3 )
      BSYALslggpppm = BSYALslggpppm + I312x3x4(j) * ( 1.D0/2.D0/(zab(p2
     &    ,q1,p3))*za(e1,p3)*za(e4,p3)*zb(p3,p2)*zab(p2,q1,p2)*mt - 1.D0
     &    /2.D0/(zab(p2,q1,p3))*za(e4,e1)*zab(p2,q1,p2)*zab(p3,q1,p2)*
     &    mt + 1.D0/2.D0/(zab(p2,q1,p3))/(za(p3,p2))*za(e1,p3)*za(e4,p3
     &    )*zab(p2,q1,p2)*mt3 - 1/(zab(p2,q1,p3))/(zab(p2,q1,p3))*za(e4
     &    ,p3)*zab(p2,q1,p2)*zab(e1,q1,p3)*mt3 - 1/(zab(p2,q1,p3))/(
     &    zab(p2,q1,p3))*za(e4,e1)*zab(p2,q1,p2)*s12*mt3 + 1/(zab(p2,q1
     &    ,p3))/(zab(p2,q1,p3))/(za(p3,p2))*za(e1,p2)*za(e4,p3)*zab(p2,
     &    q1,p2)**2*mt3 - 1.D0/2.D0/(zab(p2,q1,p3))/(zab(p2,q1,p3))/(
     &    zab(p2,q1,p3))/(za(p3,p2))*za(e1,p2)*za(e4,p2)*zab(p2,q1,p2)
     &    **2*zab(p3,q1,p3)*mt3 )
      BSYALslggpppm = BSYALslggpppm + I3m13x2x4(j) * (  - 1.D0/2.D0/(
     &    zab(p2,q1,p3))*za(e1,p3)*za(e4,p3)*zb(p3,p2)*zab(p3,q1,p3)*mt
     &     + 1.D0/2.D0/(zab(p2,q1,p3))*za(e4,e1)*zab(p3,q1,p2)*zab(p3,
     &    q1,p3)*mt - 3.D0/2.D0/(zab(p2,q1,p3))/(za(p3,p2))*za(e1,p3)*
     &    za(e4,p3)*zab(p3,q1,p3)*mt3 + 1/(zab(p2,q1,p3))/(zab(p2,q1,p3
     &    ))*za(e1,p2)*za(e4,p3)*zb(p3,p2)*zab(p3,q1,p3)*mt3 + 1/(zab(
     &    p2,q1,p3))/(zab(p2,q1,p3))*za(e4,e1)*zab(p3,q1,p3)*s12*mt3 - 
     &    2.D0/(zab(p2,q1,p3))/(zab(p2,q1,p3))/(za(p3,p2))*za(e1,p2)*
     &    za(e4,p3)*zab(p2,q1,p2)*zab(p3,q1,p3)*mt3 + 1.D0/2.D0/(zab(p2
     &    ,q1,p3))/(zab(p2,q1,p3))/(zab(p2,q1,p3))/(za(p3,p2))*za(e1,p2
     &    )*za(e4,p2)*zab(p2,q1,p2)*zab(p3,q1,p3)**2*mt3 )
      BSYALslggpppm = BSYALslggpppm + I313x2x4(j) * ( 1.D0/2.D0/(zab(p2
     &    ,q1,p3))*za(e1,p3)*za(e4,p3)*zb(p3,p2)*zab(p3,q1,p3)*mt - 1.D0
     &    /2.D0/(zab(p2,q1,p3))*za(e4,e1)*zab(p3,q1,p2)*zab(p3,q1,p3)*
     &    mt + 3.D0/2.D0/(zab(p2,q1,p3))/(za(p3,p2))*za(e1,p3)*za(e4,p3
     &    )*zab(p3,q1,p3)*mt3 - 1/(zab(p2,q1,p3))/(zab(p2,q1,p3))*za(e1
     &    ,p2)*za(e4,p3)*zb(p3,p2)*zab(p3,q1,p3)*mt3 - 1/(zab(p2,q1,p3)
     &    )/(zab(p2,q1,p3))*za(e4,e1)*zab(p3,q1,p3)*s12*mt3 + 2.D0/(
     &    zab(p2,q1,p3))/(zab(p2,q1,p3))/(za(p3,p2))*za(e1,p2)*za(e4,p3
     &    )*zab(p2,q1,p2)*zab(p3,q1,p3)*mt3 - 1.D0/2.D0/(zab(p2,q1,p3))
     &    /(zab(p2,q1,p3))/(zab(p2,q1,p3))/(za(p3,p2))*za(e1,p2)*za(e4,
     &    p2)*zab(p2,q1,p2)*zab(p3,q1,p3)**2*mt3 )
      BSYALslggpppm = BSYALslggpppm + I41x2x4x3 * (  - 1.D0/2.D0/(zab(
     &    p2,q1,p3))*za(e1,p3)*za(e4,p3)*zb(p3,p2)*s23*mt3 + 1.D0/2.D0
     &    /(zab(p2,q1,p3))*za(e1,p3)*za(e4,p3)*zb(p3,p2)*zab(p2,q1,p3)*
     &    zab(p3,q1,p2)*mt - 1.D0/2.D0/(zab(p2,q1,p3))*za(e4,e1)*zab(p2
     &    ,q1,p3)*zab(p3,q1,p2)**2*mt + 1.D0/2.D0/(zab(p2,q1,p3))*za(e4
     &    ,e1)*zab(p3,q1,p2)*s23*mt3 - 3.D0/2.D0/(zab(p2,q1,p3))/(za(p3
     &    ,p2))*za(e1,p3)*za(e4,p3)*zab(p2,q1,p2)*zab(p3,q1,p3)*mt3 + 
     &    1/(zab(p2,q1,p3))/(zab(p2,q1,p3))*za(e1,p2)*za(e4,p3)*zb(p3,
     &    p2)*zab(p2,q1,p2)*zab(p3,q1,p3)*mt3 + 1/(zab(p2,q1,p3))/(zab(
     &    p2,q1,p3))*za(e4,e1)*zab(p2,q1,p2)*zab(p3,q1,p3)*s12*mt3 - 2.D
     &    0/(zab(p2,q1,p3))/(zab(p2,q1,p3))/(za(p3,p2))*za(e1,p2)*za(e4
     &    ,p3)*zab(p2,q1,p2)**2*zab(p3,q1,p3)*mt3 + 1.D0/2.D0/(zab(p2,
     &    q1,p3))/(zab(p2,q1,p3))/(zab(p2,q1,p3))/(za(p3,p2))*za(e1,p2)
     &    *za(e4,p2)*zab(p2,q1,p2)**2*zab(p3,q1,p3)**2*mt3 )

c---- SB overall sign change to put in agreement with eq. (97)
      BSYALslggpppm = -BSYALslggpppm

      return
      end
