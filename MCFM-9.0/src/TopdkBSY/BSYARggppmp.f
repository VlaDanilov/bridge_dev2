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
 
      function BSYARggppmp(e1,p2,p3,e4,za,zb,zab,zba)
      implicit none
      include 'types.f'
      complex(dp):: BSYARggppmp
      
C-----Authors: John Campbell and Keith Ellis, March 2012
C---- arXiv:1101.5947 [hep-ph], Eq. (95)
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
      real(dp):: xbeta2,s12,s23,mt2,mt3
      complex(dp):: BSYA0ggppmp,AT0ggppmp
      integer:: e1,p2,p3,e4,j
      j=p2-1
      mt2=mt**2
      mt3=mt**3
      s12=mt2+s(1,p2)
      s23=s(p2,p3)
      xbeta2=1d0-4d0*mt**2/s23

c---- SB sign changes equivalent to those in BSYALggppmp.f
c---- a) tree level sign, b) overall sign change
      AT0ggppmp=-BSYA0ggppmp(e1,p2,p3,e4,za,zb,zab,zba)
      BSYARggppmp =  + F212(j) * (  - 2.D0/(zab(p2,q1,p2))*s12*
     &    AT0ggppmp + 1.D0/2.D0/(zab(p2,q1,p2))/(za(p3,p2))*za(e1,p2)*
     &    za(e4,p3)*zab(p3,q1,p2)**2*s12**(-1)*mt + 1.D0/2.D0/(zab(p2,
     &    q1,p2))/(za(p3,p2))*za(e1,p3)*za(e4,p3)*zab(p3,q1,p2)*
     &    s12**(-1)*mt3 + 1.D0/2.D0/(zab(p2,q1,p2))/(zab(p2,q1,p3))/(
     &    za(p3,p2))*za(e1,p2)*za(e4,p2)*zab(p3,q1,p2)**2*mt - 1/(zab(
     &    p2,q1,p2))/(zab(p2,q1,p3))/(zab(p2,q1,p3))/(za(p3,p2))*za(e1,
     &    p2)*za(e4,p2)*zab(p2,q1,p2)**2*zab(p3,q1,p2)*mt - 2.D0/(zab(
     &    p2,q1,p2))/(zab(p2,q1,p3))/(zab(p2,q1,p3))/(za(p3,p2))*za(e1,
     &    p2)*za(e4,p2)*zab(p2,q1,p2)*zab(p3,q1,p2)*mt3 - 2.D0/(zab(p2,
     &    q1,p2))/(zab(p2,q1,p3))/(zab(p2,q1,p3))/(za(p3,p2))*za(e1,p2)
     &    *za(e4,p3)*zab(p2,q1,p2)*zab(p2,q1,p3)*zab(p3,q1,p2)*mt - 4.D0
     &    /(zab(p2,q1,p2))/(zab(p2,q1,p3))/(zab(p2,q1,p3))/(za(p3,p2))*
     &    za(e1,p2)*za(e4,p3)*zab(p2,q1,p3)*zab(p3,q1,p2)*mt3 - 2.D0/(
     &    zab(p2,q1,p3))*za(e4,e1)*zab(p3,q1,p2)*s12*s23**(-1)*mt - 1/(
     &    zab(p2,q1,p3))/(za(p3,p2))*za(e1,p2)*za(e4,p3)*zab(p3,q1,p2)*
     &    mt )
      BSYARggppmp = BSYARggppmp + I2m * ( 1.D0/2.D0*AT0ggppmp )
      BSYARggppmp = BSYARggppmp + F4m1x2x3x4(j) * (  - 2.D0/(za(p3,p2))
     &    *za(e1,p2)*za(e4,p3)*zab(p3,q1,p2)**2*mt - 3.D0/(za(p3,p2))*
     &    za(e1,p3)*za(e4,p3)*zab(p3,q1,p2)*mt3 )
      BSYARggppmp = BSYARggppmp + I3m1x23x4 * ( 2.D0*mt2*AT0ggppmp - 
     &    s23*AT0ggppmp )
      BSYARggppmp = BSYARggppmp + I3m12x3x4(j) * ( 4.D0/(zab(p2,q1,p2))
     &    /(za(p3,p2))*za(e1,p2)*za(e4,p3)*zab(p2,q1,p2)*zab(p3,q1,p2)
     &    **2*s23**(-1)*mt*xbeta2**(-1) + 8.D0/(zab(p2,q1,p2))/(za(p3,
     &    p2))*za(e1,p2)*za(e4,p3)*zab(p3,q1,p2)**2*s23**(-1)*
     &    xbeta2**(-1)*mt3 + 8.D0/(zab(p2,q1,p2))/(za(p3,p2))*za(e1,p3)
     &    *za(e4,p3)*zab(p2,q1,p2)*zab(p3,q1,p2)*s23**(-1)*xbeta2**(-1)
     &    *mt3 + 1/(zab(p2,q1,p2))/(za(p3,p2))*za(e1,p3)*za(e4,p3)*zab(
     &    p3,q1,p2)*xbeta2**(-1)*mt3 - 5.D0/(zab(p2,q1,p2))/(za(p3,p2))
     &    *za(e1,p3)*za(e4,p3)*zab(p3,q1,p2)*mt3 - 2.D0/(zab(p2,q1,p2))
     &    /(zab(p2,q1,p3))*za(e4,e1)*zab(p2,q1,p2)**2*zab(p3,q1,p2)*
     &    s23**(-1)*xbeta2**(-1)*mt3 - 4.D0/(zab(p2,q1,p2))/(zab(p2,q1,
     &    p3))*za(e4,e1)*zab(p2,q1,p2)*zab(p3,q1,p2)*s23**(-1)*
     &    xbeta2**(-1)*mt2*mt3 - 2.D0/(zab(p2,q1,p2))/(zab(p2,q1,p3))*
     &    za(e4,e1)*zab(p2,q1,p3)*zab(p3,q1,p2)**2*s23**(-1)*
     &    xbeta2**(-1)*mt3 + 2.D0/(zab(p2,q1,p2))/(zab(p2,q1,p3))*za(e4
     &    ,e1)*zab(p3,q1,p2)*s12*mt3 )
      BSYARggppmp = BSYARggppmp + I3m12x3x4(j) * (  - 2.D0/(zab(p2,q1,
     &    p2))/(zab(p2,q1,p3))/(zab(p2,q1,p3))/(za(p3,p2))*za(e1,p2)*
     &    za(e4,p2)*zab(p2,q1,p2)**2*zab(p3,q1,p2)*mt3 - 4.D0/(zab(p2,
     &    q1,p2))/(zab(p2,q1,p3))/(zab(p2,q1,p3))/(za(p3,p2))*za(e1,p2)
     &    *za(e4,p3)*zab(p2,q1,p2)*zab(p2,q1,p3)*zab(p3,q1,p2)*mt3 - 2.D
     &    0/(zab(p2,q1,p3))/(za(p3,p2))*za(e1,p2)*za(e4,p2)*zab(p3,q1,
     &    p2)**2*s23**(-1)*xbeta2**(-1)*mt3 + 4.D0/(zab(p2,q1,p3))/(za(
     &    p3,p2))*za(e1,p2)*za(e4,p3)*zab(p2,q1,p2)*zab(p3,q1,p2)*
     &    s23**(-1)*xbeta2**(-1)*mt3 + 2.D0/(zab(p2,q1,p3))/(za(p3,p2))
     &    *za(e1,p2)*za(e4,p3)*zab(p3,q1,p2)*xbeta2**(-1)*mt3 )
      BSYARggppmp = BSYARggppmp + I3m2x3x41 * (  - 2.D0/(zab(p2,q1,p2))
     &    *s23*xbeta2**(-1)*mt2*AT0ggppmp + 4.D0/(zab(p2,q1,p2))*s12*
     &    xbeta2**(-1)*mt2*AT0ggppmp - 2.D0/(zab(p2,q1,p2))*za(e4,e1)*
     &    zab(p3,q1,p2)**2*s23**(-1)*xbeta2**(-1)*mt3 - 6.D0/(zab(p2,q1
     &    ,p2))/(zab(p2,q1,p2))*za(e1,p3)*za(e4,p3)*zb(p3,p2)*zab(p3,q1
     &    ,p2)*s23**(-1)*xbeta2**(-1)*mt2*mt3 - 2.D0/(zab(p2,q1,p2))/(
     &    zab(p2,q1,p2))/(zab(p2,q1,p3))/(za(p3,p2))*za(e1,p2)*za(e4,p3
     &    )*zab(p2,q1,p2)**2*zab(p3,q1,p2)*xbeta2**(-1)*mt3 + 1/(zab(p2
     &    ,q1,p2))/(zab(p2,q1,p2))/(zab(p2,q1,p3))/(za(p3,p2))*za(e1,p2
     &    )*za(e4,p3)*zab(p2,q1,p2)*zab(p2,q1,p3)*zab(p3,q1,p2)**2*mt*
     &    xbeta2**(-1) + 4.D0/(zab(p2,q1,p2))/(zab(p2,q1,p2))/(zab(p2,
     &    q1,p3))/(za(p3,p2))*za(e1,p2)*za(e4,p3)*zab(p2,q1,p3)*zab(p3,
     &    q1,p2)**2*xbeta2**(-1)*mt3 + 1/(zab(p2,q1,p2))/(zab(p2,q1,p3)
     &    )/(za(p3,p2))*za(p3,p2)*za(e1,p2)*zab(p2,q1,p2)*zab(p3,q1,p2)
     &    *zab(e4,q1,p2)*mt*xbeta2**(-1) - 1/(zab(p2,q1,p2))/(zab(p2,q1
     &    ,p3))/(za(p3,p2))*za(e1,p2)*za(e4,p2)*zab(p3,q1,p2)**2*s12*mt
     &    *xbeta2**(-1) )
      BSYARggppmp = BSYARggppmp + I3m2x3x41 * ( 1/(zab(p2,q1,p2))/(zab(
     &    p2,q1,p3))/(zab(p2,q1,p3))*za(e1,p2)*za(e4,p2)*zb(p3,p2)*zab(
     &    p2,q1,p2)*zab(p3,q1,p2)*mt3 + 2.D0/(zab(p2,q1,p2))/(zab(p2,q1
     &    ,p3))/(zab(p2,q1,p3))*za(e1,p2)*za(e4,p3)*zb(p3,p2)*zab(p2,q1
     &    ,p3)*zab(p3,q1,p2)*mt3 + 1/(zab(p2,q1,p3))*za(e1,p2)*za(e4,p3
     &    )*zb(p3,p2)*zab(p3,q1,p2)*mt + 2.D0/(zab(p2,q1,p3))*za(e4,e1)
     &    *zab(p2,q1,p2)*zab(p3,q1,p2)*s23**(-1)*xbeta2**(-1)*mt3 + 4.D0
     &    /(zab(p2,q1,p3))*za(e4,e1)*zab(p3,q1,p2)*s23**(-1)*
     &    xbeta2**(-1)*mt2*mt3 )
      BSYARggppmp = BSYARggppmp + I46m1x2x3x4(j) * ( 1/(zab(p2,q1,p2))
     &    /(za(p3,p2))*za(e1,p2)*za(e4,p2)*zab(p3,q1,p2)**3*s23**(-1)*
     &    mt*xbeta2**(-1) + 1/(zab(p2,q1,p2))/(za(p3,p2))*za(e1,p2)*za(
     &    e4,p3)*zab(p2,q1,p2)*zab(p3,q1,p2)**2*s23**(-1)*mt*
     &    xbeta2**(-1) + 2.D0/(zab(p2,q1,p2))/(za(p3,p2))*za(e1,p3)*za(
     &    e4,p3)*zab(p2,q1,p2)*zab(p3,q1,p2)*mt*xbeta2**(-1) + 4.D0/(
     &    zab(p2,q1,p2))/(za(p3,p2))*za(e1,p3)*za(e4,p3)*zab(p3,q1,p2)*
     &    xbeta2**(-1)*mt3 + 1/(zab(p2,q1,p2))/(zab(p2,q1,p2))*za(e1,p3
     &    )*za(e4,p3)*zb(p3,p2)*zab(p2,q1,p2)*zab(p3,q1,p2)*mt + 5.D0/2.
     &    D0/(zab(p2,q1,p2))/(zab(p2,q1,p2))*za(e1,p3)*za(e4,p3)*zb(p3,
     &    p2)*zab(p3,q1,p2)*mt3 + 3.D0/(zab(p2,q1,p2))/(zab(p2,q1,p2))
     &    /(za(p3,p2))*za(e1,p2)*za(e4,p3)*zab(p2,q1,p2)**2*zab(p3,q1,
     &    p2)**2*s23**(-1)*mt*xbeta2**(-1) + 3.D0/(zab(p2,q1,p2))/(zab(
     &    p2,q1,p2))/(za(p3,p2))*za(e1,p2)*za(e4,p3)*zab(p3,q1,p2)**2*
     &    s12*mt*xbeta2**(-1) + 1/(zab(p2,q1,p2))/(zab(p2,q1,p2))/(za(
     &    p3,p2))*za(e1,p3)*za(e4,p2)*zab(p3,q1,p2)**2*xbeta2**(-1)*mt3
     &     )
      BSYARggppmp = BSYARggppmp + I46m1x2x3x4(j) * ( 1/(zab(p2,q1,p2))
     &    /(zab(p2,q1,p2))/(za(p3,p2))*za(e1,p3)*za(e4,p3)*zab(p2,q1,p2
     &    )*zab(p3,q1,p2)*xbeta2**(-1)*mt3 + 1.D0/2.D0/(zab(p2,q1,p2))
     &    /(zab(p2,q1,p2))/(za(p3,p2))*za(e1,p3)*za(e4,p3)*zab(p3,q1,p2
     &    )*s23*xbeta2**(-1)*mt3 - 1/(zab(p2,q1,p2))/(zab(p2,q1,p2))/(
     &    zab(p2,q1,p3))*za(e4,e1)*zab(p2,q1,p2)**2*zab(p3,q1,p2)*s12*
     &    mt + 1/(zab(p2,q1,p2))/(zab(p2,q1,p2))/(zab(p2,q1,p3))*za(e4,
     &    e1)*zab(p3,q1,p2)*s12**2*s23*mt + 2.D0/(zab(p2,q1,p2))/(zab(
     &    p2,q1,p3))/(za(p3,p2))*za(e1,p2)*za(e4,p3)*zab(p2,q1,p2)**2*
     &    zab(p3,q1,p2)*mt + 1/(zab(p2,q1,p2))/(zab(p2,q1,p3))/(za(p3,
     &    p2))*za(e1,p2)*za(e4,p3)*zab(p2,q1,p2)*zab(p3,q1,p2)*s23*mt
     &     + 2.D0/(zab(p2,q1,p2))/(zab(p2,q1,p3))/(za(p3,p2))*za(e1,p2)
     &    *za(e4,p3)*zab(p3,q1,p2)*s23*mt3 + 1/(zab(p2,q1,p2))/(zab(p2,
     &    q1,p3))/(zab(p2,q1,p3))*za(e1,p2)*za(e4,p2)*zb(p3,p2)*zab(p2,
     &    q1,p2)**2*zab(p3,q1,p2)*mt + 2.D0/(zab(p2,q1,p2))/(zab(p2,q1,
     &    p3))/(zab(p2,q1,p3))*za(e1,p2)*za(e4,p2)*zb(p3,p2)*zab(p2,q1,
     &    p2)*zab(p3,q1,p2)*mt3 )
      BSYARggppmp = BSYARggppmp + I46m1x2x3x4(j) * ( 2.D0/(zab(p2,q1,p2
     &    ))/(zab(p2,q1,p3))/(zab(p2,q1,p3))*za(e1,p2)*za(e4,p3)*zb(p3,
     &    p2)*zab(p2,q1,p2)*zab(p2,q1,p3)*zab(p3,q1,p2)*mt + 4.D0/(zab(
     &    p2,q1,p2))/(zab(p2,q1,p3))/(zab(p2,q1,p3))*za(e1,p2)*za(e4,p3
     &    )*zb(p3,p2)*zab(p2,q1,p3)*zab(p3,q1,p2)*mt3 + 1/(zab(p2,q1,p3
     &    ))/(za(p3,p2))*za(e1,p2)*za(e4,p2)*zab(p2,q1,p2)*zab(p3,q1,p2
     &    )**2*s23**(-1)*mt*xbeta2**(-1) + 2.D0/(zab(p2,q1,p3))/(za(p3,
     &    p2))*za(e1,p2)*za(e4,p2)*zab(p3,q1,p2)**2*s23**(-1)*
     &    xbeta2**(-1)*mt3 )
      BSYARggppmp = BSYARggppmp + F2m23 * ( 1/(za(p3,p2))*za(e1,p2)*za(
     &    e4,p3)*zab(p3,q1,p2)**2*s23**(-2)*mt*xbeta2**(-1) + 1/(za(p3,
     &    p2))*za(e1,p3)*za(e4,p2)*zab(p3,q1,p2)**2*s23**(-2)*mt*
     &    xbeta2**(-1) + 1/(zab(p2,q1,p2))/(zab(p2,q1,p3))/(zab(p2,q1,
     &    p3))/(za(p3,p2))*za(e1,p2)*za(e4,p2)*zab(p2,q1,p2)**2*zab(p3,
     &    q1,p2)*mt + 2.D0/(zab(p2,q1,p2))/(zab(p2,q1,p3))/(zab(p2,q1,
     &    p3))/(za(p3,p2))*za(e1,p2)*za(e4,p3)*zab(p2,q1,p2)*zab(p2,q1,
     &    p3)*zab(p3,q1,p2)*mt + 1/(zab(p2,q1,p3))*za(e4,e1)*zab(p2,q1,
     &    p2)*zab(p3,q1,p2)*s23**(-1)*mt*xbeta2**(-1) + 1/(zab(p2,q1,p3
     &    ))*za(e4,e1)*zab(p2,q1,p2)*zab(p3,q1,p2)*s23**(-1)*mt + 2.D0
     &    /(zab(p2,q1,p3))*za(e4,e1)*zab(p3,q1,p2)*s23**(-1)*
     &    xbeta2**(-1)*mt3 + 1/(zab(p2,q1,p3))/(za(p3,p2))*za(e1,p2)*
     &    za(e4,p2)*zab(p2,q1,p2)*zab(p3,q1,p2)**2*s23**(-2)*mt*
     &    xbeta2**(-1) + 3.D0/2.D0/(zab(p2,q1,p3))/(za(p3,p2))*za(e1,p2
     &    )*za(e4,p2)*zab(p3,q1,p2)**2*s23**(-1)*mt*xbeta2**(-1) - 2.D0
     &    /(zab(p2,q1,p3))/(za(p3,p2))*za(e1,p2)*za(e4,p3)*zab(p2,q1,p2
     &    )*zab(p3,q1,p2)*s23**(-1)*mt*xbeta2**(-1) )
      BSYARggppmp = BSYARggppmp + F2m23 * (  - 1/(zab(p2,q1,p3))/(za(p3
     &    ,p2))*za(e1,p2)*za(e4,p3)*zab(p3,q1,p2)*mt*xbeta2**(-1) - 1/(
     &    zab(p2,q1,p3))/(za(p3,p2))*za(e1,p3)*za(e4,p3)*zab(p2,q1,p2)*
     &    zab(p2,q1,p3)*zab(p3,q1,p2)*s23**(-2)*mt*xbeta2**(-1) - 3.D0/
     &    2.D0/(zab(p2,q1,p3))/(za(p3,p2))*za(e1,p3)*za(e4,p3)*zab(p2,
     &    q1,p3)*zab(p3,q1,p2)*s23**(-1)*mt*xbeta2**(-1) )
      BSYARggppmp = BSYARggppmp - 2.D0/(za(p3,p2))*za(e1,p2)*za(e4,p3)*
     & zab(p3,q1,p2)**2*s23**(-2)*mt*xbeta2**(-1) - 2.D0/(za(p3,p2))*
     &    za(e1,p3)*za(e4,p2)*zab(p3,q1,p2)**2*s23**(-2)*mt*
     &    xbeta2**(-1) + 1.D0/2.D0/(zab(p2,q1,p2))/(za(p3,p2))*za(e1,p2
     &    )*za(e4,p3)*zab(p3,q1,p2)**2*s12**(-1)*mt - 1.D0/2.D0/(zab(p2
     &    ,q1,p2))/(za(p3,p2))*za(e1,p3)*za(e4,p3)*zab(p2,q1,p2)*zab(p3
     &    ,q1,p2)*s12**(-1)*mt - 1.D0/2.D0/(zab(p2,q1,p2))/(zab(p2,q1,
     &    p3))*za(e1,p2)*za(e4,p2)*zb(p3,p2)*zab(p3,q1,p2)**2*s23**(-1)
     &    *mt - 1.D0/2.D0/(zab(p2,q1,p2))/(zab(p2,q1,p3))*za(e4,e1)*
     &    zab(p2,q1,p3)*zab(p3,q1,p2)**2*s23**(-1)*mt - 2.D0/(zab(p2,q1
     &    ,p2))/(zab(p2,q1,p3))/(zab(p2,q1,p3))/(za(p3,p2))*za(e1,p2)*
     &    za(e4,p2)*zab(p2,q1,p2)**2*zab(p3,q1,p2)*mt - 4.D0/(zab(p2,q1
     &    ,p2))/(zab(p2,q1,p3))/(zab(p2,q1,p3))/(za(p3,p2))*za(e1,p2)*
     &    za(e4,p3)*zab(p2,q1,p2)*zab(p2,q1,p3)*zab(p3,q1,p2)*mt - 2.D0
     &    /(zab(p2,q1,p3))*za(e4,e1)*zab(p2,q1,p2)*zab(p3,q1,p2)*
     &    s23**(-1)*mt*xbeta2**(-1)
      BSYARggppmp = BSYARggppmp - 2.D0/(zab(p2,q1,p3))*za(e4,e1)*zab(p2
     & ,q1,p2)*zab(p3,q1,p2)*s23**(-1)*mt - 4.D0/(zab(p2,q1,p3))*za(e4,
     &    e1)*zab(p3,q1,p2)*s23**(-1)*xbeta2**(-1)*mt3 - 2.D0/(zab(p2,
     &    q1,p3))/(za(p3,p2))*za(e1,p2)*za(e4,p2)*zab(p2,q1,p2)*zab(p3,
     &    q1,p2)**2*s23**(-2)*mt*xbeta2**(-1) - 3.D0/(zab(p2,q1,p3))/(
     &    za(p3,p2))*za(e1,p2)*za(e4,p2)*zab(p3,q1,p2)**2*s23**(-1)*mt*
     &    xbeta2**(-1) + 4.D0/(zab(p2,q1,p3))/(za(p3,p2))*za(e1,p2)*za(
     &    e4,p3)*zab(p2,q1,p2)*zab(p3,q1,p2)*s23**(-1)*mt*xbeta2**(-1)
     &     + 2.D0/(zab(p2,q1,p3))/(za(p3,p2))*za(e1,p2)*za(e4,p3)*zab(
     &    p3,q1,p2)*mt*xbeta2**(-1) + 2.D0/(zab(p2,q1,p3))/(za(p3,p2))*
     &    za(e1,p3)*za(e4,p3)*zab(p2,q1,p2)*zab(p2,q1,p3)*zab(p3,q1,p2)
     &    *s23**(-2)*mt*xbeta2**(-1) + 3.D0/(zab(p2,q1,p3))/(za(p3,p2))
     &    *za(e1,p3)*za(e4,p3)*zab(p2,q1,p3)*zab(p3,q1,p2)*s23**(-1)*mt
     &    *xbeta2**(-1)
      
c---- SB sign change, see point b) above
      BSYARggppmp = -BSYARggppmp
      return
      end
