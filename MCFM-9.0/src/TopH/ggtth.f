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
 
      subroutine ggtth(denr,denb,denq1,denq2,dena1,dena2,wtgg)
      implicit none
      include 'types.f'
      
C***********************************************************************
C     Author: R.K. Ellis                                               *
C     April, 2010.                                                     *
C     calculate the element squared and subtraction terms              *
C     for the process                                                  *
c     My notation                                                      *
C     g(-p1) + g(-p2)=bbar(p6)+e-(p7)+nubar(p8)+b(p5)+nu(p3)+e+(p4)    *
C     +b(p9)+bbar(p10)                                                 *
C                                                                      *
C h=9+10                                                               *
C q=3+4+5                                                              *
C r=3+4+5+h                                                            *
C a=-6-7-8                                                             *
C b=-6-7-8-h                                                           *
C since momenta 3,5,6,8,9,10 are not needed they are reused to         *
C represent the de-massified vectors.                                  *
C Thus q4 = q de-massified wrt p4 etc;                                 *
C this vector has been stored in wrapper routine (qqb_tth) in p(3,mu)  *
C q4(mu)=q(mu)-qsq/2/p4Dq*p4(mu)-->p(3,mu)                             *
C a7(mu)=a(mu)-asq/2/p7Da*p7(mu)-->p(5,mu)                             *
C r1(mu)=r(mu)-rsq/2/p1Dr*p1(mu)-->p(6,mu)                             *
C r2(mu)=r(mu)-rsq/2/p2Dr*p2(mu)-->p(8,mu)                             *
C b1(mu)=b(mu)-bsq/2/p1Db*p1(mu)-->p(9,mu)                             *
C b2(mu)=b(mu)-bsq/2/p2Db*p2(mu)-->p(10,mu)                            *
C Formula generated by the program ggtth.frm                           *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_com.f'
      include 'masses.f'
      complex(dp):: mab(2,2),mba(2,2)
      real(dp):: denr,denb,denq1,denq2,dena1,dena2,wtgg,s12
      integer:: h1,h2
      integer:: q4,a7,r1,r2,b1,b2
      parameter(q4=3,a7=5,r1=6,r2=8,b1=9,b2=10)
      integer:: p1,p2,p4,p7
      parameter(p1=1,p2=2,p4=4,p7=7)
      s12=real(za(p1,p2)*zb(p2,p1))
      mab(1,1)= + mt*denr**(-1)*dena2**(-1) * (  - za(p1,p2)**3*za(p7,
     &    a7)*za(q4,r2)*zb(p1,a7)*zb(p2,r2)*zb(p4,q4)*s12**(-2) + za(p1
     &    ,p2)**2*za(p1,q4)*za(p2,r2)*za(p7,a7)*zb(p1,a7)*zb(p2,r2)*zb(
     &    p4,q4)*s12**(-2) - za(p1,p2)**2*za(p1,q4)*za(p7,a7)*zb(p1,a7)
     &    *zb(p4,q4)*s12**(-1) - za(p1,p2)**2*za(p1,r1)*za(p2,p7)*za(q4
     &    ,r2)*zb(p1,r1)*zb(p2,r2)*zb(p4,q4)*s12**(-2) + za(p1,p2)**2*
     &    za(p1,r1)*za(p2,r2)*za(p7,a7)*zb(p1,a7)*zb(p2,r2)*zb(p4,r1)*
     &    s12**(-2) - za(p1,p2)**2*za(p1,r1)*za(p7,a7)*zb(p1,a7)*zb(p4,
     &    r1)*s12**(-1) )
      mab(1,1) = mab(1,1) + mt*denr**(-1) * ( za(p1,p2)**2*za(p1,p7)*
     &    za(q4,r1)*zb(p1,r1)*zb(p4,q4)*s12**(-2) - za(p1,p2)**2*za(p1,
     &    q4)*za(p7,a7)*zb(p1,a7)*zb(p4,q4)*s12**(-2) - za(p1,p2)**2*
     &    za(p1,r1)*za(p7,a7)*zb(p1,a7)*zb(p4,r1)*s12**(-2) )
      mab(1,1) = mab(1,1) + mt*denb**(-1)*denq1**(-1) * ( za(p1,p2)**2*
     &    za(p1,q4)*za(p2,b2)*za(p7,a7)*zb(p1,a7)*zb(p2,b2)*zb(p4,q4)*
     &    s12**(-2) + za(p1,p2)**2*za(p1,q4)*za(p2,b2)*za(p7,b1)*zb(p1,
     &    b1)*zb(p2,b2)*zb(p4,q4)*s12**(-2) - za(p1,p2)**2*za(p1,b1)*
     &    za(p2,b2)*za(p7,a7)*zb(p1,b1)*zb(p2,p4)*zb(a7,b2)*s12**(-2)
     &     - za(p1,p2)**2*za(p2,b2)*za(p7,a7)*zb(p2,p4)*zb(a7,b2)*
     &    s12**(-1) + za(p1,p2)*za(p1,q4)*za(p2,b2)*za(p7,a7)*zb(p4,q4)
     &    *zb(a7,b2)*s12**(-1) )
      mab(1,1) = mab(1,1) + mt*denb**(-1) * (  - za(p1,p2)**2*za(p1,q4)
     &    *za(p7,a7)*zb(p1,a7)*zb(p4,q4)*s12**(-2) - za(p1,p2)**2*za(p1
     &    ,q4)*za(p7,b1)*zb(p1,b1)*zb(p4,q4)*s12**(-2) + za(p1,p2)**2*
     &    za(p1,b1)*za(p7,a7)*zb(p1,p4)*zb(a7,b1)*s12**(-2) )
      mab(1,1) = mab(1,1) + mt*denq1**(-1)*dena2**(-1) * ( za(p1,p2)**3
     &    *za(p1,b1)*za(p7,a7)*zb(p1,a7)*zb(p1,b1)*zb(p2,p4)*s12**(-2)
     &     - za(p1,p2)**3*za(p2,r2)*za(p7,a7)*zb(p1,a7)*zb(p2,p4)*zb(p2
     &    ,r2)*s12**(-2) + za(p1,p2)**3*za(p7,a7)*zb(p1,a7)*zb(p2,p4)*
     &    s12**(-1) + za(p1,p2)**2*za(p1,q4)*za(p2,p7)*za(r1,b2)*zb(p1,
     &    r1)*zb(p2,b2)*zb(p4,q4)*s12**(-2) + za(p1,p2)**2*za(p1,q4)*
     &    za(p2,r2)*za(p7,a7)*zb(p1,a7)*zb(p2,r2)*zb(p4,q4)*s12**(-2)
     &     + za(p1,p2)**2*za(p1,q4)*za(p2,b2)*za(p7,a7)*zb(p1,a7)*zb(p2
     &    ,b2)*zb(p4,q4)*s12**(-2) - za(p1,p2)**2*za(p1,q4)*za(p7,a7)*
     &    zb(p1,a7)*zb(p4,q4)*s12**(-1) - za(p1,p2)**2*za(p1,b1)*za(p2,
     &    r2)*za(p7,a7)*zb(p1,a7)*zb(p2,p4)*zb(r2,b1)*s12**(-2) )
      mab(1,1) = mab(1,1) + mt**3*denr**(-1)*dena2**(-1) * (  - za(p1,
     &    p2)**3*za(p7,a7)*zb(p1,a7)*zb(p2,p4)*s12**(-2) - za(p1,p2)**2
     &    *za(p1,r1)*za(p2,p7)*zb(p1,r1)*zb(p2,p4)*s12**(-2) + za(p1,p2
     &    )*za(p1,q4)*za(p2,p7)*zb(p4,q4)*s12**(-1) + za(p1,p2)*za(p1,
     &    r1)*za(p2,p7)*zb(p4,r1)*s12**(-1) )
      mab(1,1) = mab(1,1) + mt**3*denr**(-1) * ( za(p1,p2)**2*za(p1,p7)
     &    *zb(p1,p4)*s12**(-2) )
      mab(1,1) = mab(1,1) + mt**3*denb**(-1)*denq1**(-1) * (  - za(p1,
     &    p2)**3*za(p7,a7)*zb(p1,a7)*zb(p2,p4)*s12**(-2) - za(p1,p2)**3
     &    *za(p7,b1)*zb(p1,b1)*zb(p2,p4)*s12**(-2) - za(p1,p2)**2*za(p1
     &    ,b1)*za(p2,p7)*zb(p1,b1)*zb(p2,p4)*s12**(-2) - za(p1,p2)**2*
     &    za(p2,p7)*zb(p2,p4)*s12**(-1) + za(p1,p2)*za(p1,q4)*za(p2,p7)
     &    *zb(p4,q4)*s12**(-1) )
      mab(1,1) = mab(1,1) + mt**3*denb**(-1) * ( za(p1,p2)**2*za(p1,p7)
     &    *zb(p1,p4)*s12**(-2) )
      mab(1,1) = mab(1,1) + mt**3*denq1**(-1)*dena2**(-1) * (  - za(p1,
     &    p2)**3*za(p7,a7)*zb(p1,a7)*zb(p2,p4)*s12**(-2) - za(p1,p2)**2
     &    *za(p1,r1)*za(p2,p7)*zb(p1,r1)*zb(p2,p4)*s12**(-2) - za(p1,p2
     &    )**2*za(p1,b1)*za(p2,p7)*zb(p1,b1)*zb(p2,p4)*s12**(-2) - za(
     &    p1,p2)**2*za(p2,p7)*zb(p2,p4)*s12**(-1) + za(p1,p2)*za(p1,q4)
     &    *za(p2,p7)*zb(p4,q4)*s12**(-1) )

      mab(2,2)= + mt*denr**(-1)*dena2**(-1) * (  - za(p1,p7)*za(p2,r2)*
     &    za(q4,r1)*zb(p1,p2)**2*zb(p1,r1)*zb(p2,r2)*zb(p4,q4)*
     &    s12**(-2) + za(p1,p7)*za(q4,r1)*zb(p1,p2)**2*zb(p1,r1)*zb(p4,
     &    q4)*s12**(-1) + za(p1,r1)*za(p2,q4)*za(p7,a7)*zb(p1,p2)**2*
     &    zb(p1,r1)*zb(p2,a7)*zb(p4,q4)*s12**(-2) + za(p1,r1)*za(p2,r2)
     &    *za(p7,a7)*zb(p1,p2)**2*zb(p1,r1)*zb(p2,a7)*zb(p4,r2)*
     &    s12**(-2) - za(p7,a7)*za(q4,r1)*zb(p1,p2)*zb(p1,r1)*zb(p2,a7)
     &    *zb(p4,q4)*s12**(-1) )
      mab(2,2) = mab(2,2) + mt*denr**(-1) * ( za(p1,p7)*za(q4,r1)*zb(p1
     &    ,p2)**2*zb(p1,r1)*zb(p4,q4)*s12**(-2) - za(p1,q4)*za(p7,a7)*
     &    zb(p1,p2)**2*zb(p1,a7)*zb(p4,q4)*s12**(-2) - za(p1,r1)*za(p7,
     &    a7)*zb(p1,p2)**2*zb(p1,a7)*zb(p4,r1)*s12**(-2) )
      mab(2,2) = mab(2,2) + mt*denb**(-1)*denq1**(-1) * ( za(p1,b1)*za(
     &    p2,q4)*za(p7,a7)*zb(p1,p2)**3*zb(p4,q4)*zb(a7,b1)*s12**(-2)
     &     + za(p1,b1)*za(p2,q4)*za(p7,a7)*zb(p1,p2)**2*zb(p1,b1)*zb(p2
     &    ,a7)*zb(p4,q4)*s12**(-2) + za(p1,b1)*za(p2,q4)*za(p7,b2)*zb(
     &    p1,p2)**2*zb(p1,b1)*zb(p2,b2)*zb(p4,q4)*s12**(-2) - za(p1,b1)
     &    *za(p2,b2)*za(p7,a7)*zb(p1,p2)**2*zb(p1,p4)*zb(p2,b2)*zb(a7,
     &    b1)*s12**(-2) + za(p2,q4)*za(p7,a7)*zb(p1,p2)**2*zb(p2,a7)*
     &    zb(p4,q4)*s12**(-1) + za(p2,q4)*za(p7,b2)*zb(p1,p2)**2*zb(p2,
     &    b2)*zb(p4,q4)*s12**(-1) )
      mab(2,2) = mab(2,2) + mt*denb**(-1) * (  - za(p1,q4)*za(p7,a7)*
     &    zb(p1,p2)**2*zb(p1,a7)*zb(p4,q4)*s12**(-2) - za(p1,q4)*za(p7,
     &    b1)*zb(p1,p2)**2*zb(p1,b1)*zb(p4,q4)*s12**(-2) + za(p1,b1)*
     &    za(p7,a7)*zb(p1,p2)**2*zb(p1,p4)*zb(a7,b1)*s12**(-2) )
      mab(2,2) = mab(2,2) + mt*denq1**(-1)*dena2**(-1) * (  - za(p1,p7)
     &    *za(p1,b1)*za(p2,q4)*zb(p1,p2)**3*zb(p1,b1)*zb(p4,q4)*
     &    s12**(-2) + za(p1,p7)*za(p2,q4)*za(p2,r2)*zb(p1,p2)**3*zb(p2,
     &    r2)*zb(p4,q4)*s12**(-2) + za(p1,p7)*za(p2,q4)*za(r2,b1)*zb(p1
     &    ,p2)**2*zb(p1,b1)*zb(p2,r2)*zb(p4,q4)*s12**(-2) - za(p1,p7)*
     &    za(p2,q4)*zb(p1,p2)**3*zb(p4,q4)*s12**(-1) + za(p1,r1)*za(p2,
     &    q4)*za(p7,a7)*zb(p1,p2)**2*zb(p1,r1)*zb(p2,a7)*zb(p4,q4)*
     &    s12**(-2) - za(p1,r1)*za(p2,b2)*za(p7,a7)*zb(p1,p2)**2*zb(p1,
     &    p4)*zb(p2,a7)*zb(r1,b2)*s12**(-2) + za(p1,b1)*za(p2,q4)*za(p7
     &    ,a7)*zb(p1,p2)**2*zb(p1,b1)*zb(p2,a7)*zb(p4,q4)*s12**(-2) + 
     &    za(p2,q4)*za(p7,a7)*zb(p1,p2)**2*zb(p2,a7)*zb(p4,q4)*
     &    s12**(-1) )
      mab(2,2) = mab(2,2) + mt**3*denr**(-1)*dena2**(-1) * ( za(p1,p7)*
     &    za(p2,q4)*zb(p1,p2)**3*zb(p4,q4)*s12**(-2) + za(p1,p7)*za(p2,
     &    r2)*zb(p1,p2)**3*zb(p4,r2)*s12**(-2) - za(p1,p7)*za(p2,r2)*
     &    zb(p1,p2)**2*zb(p1,p4)*zb(p2,r2)*s12**(-2) + za(p1,p7)*zb(p1,
     &    p2)**2*zb(p1,p4)*s12**(-1) - za(p7,a7)*zb(p1,p2)*zb(p1,p4)*
     &    zb(p2,a7)*s12**(-1) )
      mab(2,2) = mab(2,2) + mt**3*denr**(-1) * ( za(p1,p7)*zb(p1,p2)**2
     &    *zb(p1,p4)*s12**(-2) )
      mab(2,2) = mab(2,2) + mt**3*denb**(-1)*denq1**(-1) * ( za(p1,p7)*
     &    za(p2,q4)*zb(p1,p2)**3*zb(p4,q4)*s12**(-2) - za(p1,p7)*za(p2,
     &    b2)*zb(p1,p2)**2*zb(p1,p4)*zb(p2,b2)*s12**(-2) - za(p7,a7)*
     &    zb(p1,p2)*zb(p1,p4)*zb(p2,a7)*s12**(-1) - za(p7,b2)*zb(p1,p2)
     &    *zb(p1,p4)*zb(p2,b2)*s12**(-1) )
      mab(2,2) = mab(2,2) + mt**3*denb**(-1) * ( za(p1,p7)*zb(p1,p2)**2
     &    *zb(p1,p4)*s12**(-2) )
      mab(2,2) = mab(2,2) + mt**3*denq1**(-1)*dena2**(-1) * ( za(p1,p7)
     &    *za(p2,q4)*zb(p1,p2)**3*zb(p4,q4)*s12**(-2) - za(p1,p7)*za(p2
     &    ,r2)*zb(p1,p2)**2*zb(p1,p4)*zb(p2,r2)*s12**(-2) - za(p1,p7)*
     &    za(p2,b2)*zb(p1,p2)**2*zb(p1,p4)*zb(p2,b2)*s12**(-2) + za(p1,
     &    p7)*zb(p1,p2)**2*zb(p1,p4)*s12**(-1) - za(p7,a7)*zb(p1,p2)*
     &    zb(p1,p4)*zb(p2,a7)*s12**(-1) )

      mab(1,2)= + mt*denr**(-1)*dena2**(-1) * (  - za(p1,p7)*za(p1,r1)*
     &    za(q4,r2)*zb(p2,r1)*zb(p2,r2)*zb(p4,q4)*s12**(-1) + za(p1,q4)
     &    *za(p1,r2)*za(p7,a7)*zb(p2,a7)*zb(p2,r2)*zb(p4,q4)*s12**(-1)
     &     + za(p1,r1)*za(p1,r2)*za(p7,a7)*zb(p2,a7)*zb(p2,r2)*zb(p4,r1
     &    )*s12**(-1) )
      mab(1,2) = mab(1,2) + mt*denb**(-1)*denq1**(-1) * ( za(p1,q4)*za(
     &    p1,b2)*za(p7,a7)*zb(p2,a7)*zb(p2,b2)*zb(p4,q4)*s12**(-1) + 
     &    za(p1,q4)*za(p1,b2)*za(p7,b2)*zb(p2,b2)**2*zb(p4,q4)*
     &    s12**(-1) - za(p1,b1)**2*za(p7,a7)*zb(p2,p4)*zb(p2,b1)*zb(a7,
     &    b1)*s12**(-1) )
      mab(1,2) = mab(1,2) + mt*denq1**(-1)*dena2**(-1) * (  - za(p1,p2)
     &    *za(p1,r2)*za(p7,a7)*zb(p2,p4)*zb(p2,a7)*zb(p2,r2)*s12**(-1)
     &     - za(p1,p7)*za(p1,q4)*za(p1,b2)*zb(p1,p2)*zb(p2,b2)*zb(p4,q4
     &    )*s12**(-1) + za(p1,p7)*za(p1,q4)*za(r2,b2)*zb(p2,r2)*zb(p2,
     &    b2)*zb(p4,q4)*s12**(-1) + za(p1,q4)*za(p1,r2)*za(p7,a7)*zb(p2
     &    ,a7)*zb(p2,r2)*zb(p4,q4)*s12**(-1) + za(p1,q4)*za(p1,b2)*za(
     &    p7,a7)*zb(p2,a7)*zb(p2,b2)*zb(p4,q4)*s12**(-1) - za(p1,r1)*
     &    za(p1,b1)*za(p7,a7)*zb(p2,p4)*zb(p2,a7)*zb(r1,b1)*s12**(-1) )
      mab(1,2) = mab(1,2) + mt**3*denr**(-1)*dena2**(-1) * (  - za(p1,
     &    p7)*za(p1,r1)*zb(p2,p4)*zb(p2,r1)*s12**(-1) )
      mab(1,2) = mab(1,2) + mt**3*denb**(-1)*denq1**(-1) * (  - za(p1,
     &    p7)*za(p1,b1)*zb(p2,p4)*zb(p2,b1)*s12**(-1) )
      mab(1,2) = mab(1,2) + mt**3*denq1**(-1)*dena2**(-1) * (  - za(p1,
     &    p7)*za(p1,r1)*zb(p2,p4)*zb(p2,r1)*s12**(-1) - za(p1,p7)*za(p1
     &    ,b1)*zb(p2,p4)*zb(p2,b1)*s12**(-1) )

      mab(2,1)= + mt*denr**(-1)*dena2**(-1) * (  - za(p2,p7)*za(p2,r2)*
     &    za(q4,r1)*zb(p1,r1)*zb(p1,r2)*zb(p4,q4)*s12**(-1) + za(p2,q4)
     &    *za(p2,r1)*za(p7,a7)*zb(p1,a7)*zb(p1,r1)*zb(p4,q4)*s12**(-1)
     &     + za(p2,r1)*za(p2,r2)*za(p7,a7)*zb(p1,a7)*zb(p1,r1)*zb(p4,r2
     &    )*s12**(-1) )
      mab(2,1) = mab(2,1) + mt*denb**(-1)*denq1**(-1) * ( za(p2,q4)*za(
     &    p2,b1)*za(p7,a7)*zb(p1,a7)*zb(p1,b1)*zb(p4,q4)*s12**(-1) + 
     &    za(p2,q4)*za(p2,b1)*za(p7,b1)*zb(p1,b1)**2*zb(p4,q4)*
     &    s12**(-1) - za(p2,b2)**2*za(p7,a7)*zb(p1,p4)*zb(p1,b2)*zb(a7,
     &    b2)*s12**(-1) )
      mab(2,1) = mab(2,1) + mt*denq1**(-1)*dena2**(-1) * ( za(p1,p2)*
     &    za(p2,b2)*za(p7,a7)*zb(p1,p4)*zb(p1,a7)*zb(p1,b2)*s12**(-1)
     &     + za(p2,p7)*za(p2,q4)*za(p2,r2)*zb(p1,p2)*zb(p1,r2)*zb(p4,q4
     &    )*s12**(-1) + za(p2,p7)*za(p2,q4)*za(r1,b1)*zb(p1,r1)*zb(p1,
     &    b1)*zb(p4,q4)*s12**(-1) + za(p2,q4)*za(p2,r1)*za(p7,a7)*zb(p1
     &    ,a7)*zb(p1,r1)*zb(p4,q4)*s12**(-1) + za(p2,q4)*za(p2,b1)*za(
     &    p7,a7)*zb(p1,a7)*zb(p1,b1)*zb(p4,q4)*s12**(-1) - za(p2,r2)*
     &    za(p2,b2)*za(p7,a7)*zb(p1,p4)*zb(p1,a7)*zb(r2,b2)*s12**(-1) )
      mab(2,1) = mab(2,1) + mt**3*denr**(-1)*dena2**(-1) * (  - za(p2,
     &    p7)*za(p2,r2)*zb(p1,p4)*zb(p1,r2)*s12**(-1) )
      mab(2,1) = mab(2,1) + mt**3*denb**(-1)*denq1**(-1) * (  - za(p2,
     &    p7)*za(p2,b2)*zb(p1,p4)*zb(p1,b2)*s12**(-1) )
      mab(2,1) = mab(2,1) + mt**3*denq1**(-1)*dena2**(-1) * (  - za(p2,
     &    p7)*za(p2,r2)*zb(p1,p4)*zb(p1,r2)*s12**(-1) - za(p2,p7)*za(p2
     &    ,b2)*zb(p1,p4)*zb(p1,b2)*s12**(-1) )

      mba(1,1)= + mt*denr**(-1)*dena1**(-1) * ( za(p1,p2)**3*za(p7,a7)*
     &    za(q4,r1)*zb(p1,r1)*zb(p2,a7)*zb(p4,q4)*s12**(-2) - za(p1,p2)
     &    **2*za(p1,p7)*za(p2,r2)*za(q4,r1)*zb(p1,r1)*zb(p2,r2)*zb(p4,
     &    q4)*s12**(-2) + za(p1,p2)**2*za(p1,r1)*za(p2,q4)*za(p7,a7)*
     &    zb(p1,r1)*zb(p2,a7)*zb(p4,q4)*s12**(-2) + za(p1,p2)**2*za(p1,
     &    r1)*za(p2,r2)*za(p7,a7)*zb(p1,r1)*zb(p2,a7)*zb(p4,r2)*
     &    s12**(-2) - za(p1,p2)**2*za(p2,q4)*za(p7,a7)*zb(p2,a7)*zb(p4,
     &    q4)*s12**(-1) - za(p1,p2)**2*za(p2,r2)*za(p7,a7)*zb(p2,a7)*
     &    zb(p4,r2)*s12**(-1) )
      mba(1,1) = mba(1,1) + mt*denr**(-1) * (  - za(p1,p2)**2*za(p1,p7)
     &    *za(q4,r1)*zb(p1,r1)*zb(p4,q4)*s12**(-2) + za(p1,p2)**2*za(p1
     &    ,q4)*za(p7,a7)*zb(p1,a7)*zb(p4,q4)*s12**(-2) + za(p1,p2)**2*
     &    za(p1,r1)*za(p7,a7)*zb(p1,a7)*zb(p4,r1)*s12**(-2) )
      mba(1,1) = mba(1,1) + mt*denb**(-1)*denq2**(-1) * ( za(p1,p2)**2*
     &    za(p1,b1)*za(p2,q4)*za(p7,a7)*zb(p1,b1)*zb(p2,a7)*zb(p4,q4)*
     &    s12**(-2) + za(p1,p2)**2*za(p1,b1)*za(p2,q4)*za(p7,b2)*zb(p1,
     &    b1)*zb(p2,b2)*zb(p4,q4)*s12**(-2) - za(p1,p2)**2*za(p1,b1)*
     &    za(p2,b2)*za(p7,a7)*zb(p1,p4)*zb(p2,b2)*zb(a7,b1)*s12**(-2)
     &     - za(p1,p2)**2*za(p1,b1)*za(p7,a7)*zb(p1,p4)*zb(a7,b1)*
     &    s12**(-1) - za(p1,p2)*za(p1,b1)*za(p2,q4)*za(p7,a7)*zb(p4,q4)
     &    *zb(a7,b1)*s12**(-1) )
      mba(1,1) = mba(1,1) + mt*denb**(-1) * ( za(p1,p2)**2*za(p1,q4)*
     &    za(p7,a7)*zb(p1,a7)*zb(p4,q4)*s12**(-2) + za(p1,p2)**2*za(p1,
     &    q4)*za(p7,b1)*zb(p1,b1)*zb(p4,q4)*s12**(-2) - za(p1,p2)**2*
     &    za(p1,b1)*za(p7,a7)*zb(p1,p4)*zb(a7,b1)*s12**(-2) )
      mba(1,1) = mba(1,1) + mt*denq2**(-1)*dena1**(-1) * ( za(p1,p2)**3
     &    *za(p1,r1)*za(p7,a7)*zb(p1,p4)*zb(p1,r1)*zb(p2,a7)*s12**(-2)
     &     - za(p1,p2)**3*za(p2,b2)*za(p7,a7)*zb(p1,p4)*zb(p2,a7)*zb(p2
     &    ,b2)*s12**(-2) - za(p1,p2)**3*za(p7,a7)*zb(p1,p4)*zb(p2,a7)*
     &    s12**(-1) + za(p1,p2)**2*za(p1,p7)*za(p2,q4)*za(r2,b1)*zb(p1,
     &    b1)*zb(p2,r2)*zb(p4,q4)*s12**(-2) + za(p1,p2)**2*za(p1,r1)*
     &    za(p2,q4)*za(p7,a7)*zb(p1,r1)*zb(p2,a7)*zb(p4,q4)*s12**(-2)
     &     - za(p1,p2)**2*za(p1,r1)*za(p2,b2)*za(p7,a7)*zb(p1,p4)*zb(p2
     &    ,a7)*zb(r1,b2)*s12**(-2) + za(p1,p2)**2*za(p1,b1)*za(p2,q4)*
     &    za(p7,a7)*zb(p1,b1)*zb(p2,a7)*zb(p4,q4)*s12**(-2) - za(p1,p2)
     &    **2*za(p2,q4)*za(p7,a7)*zb(p2,a7)*zb(p4,q4)*s12**(-1) )
      mba(1,1) = mba(1,1) + mt**3*denr**(-1)*dena1**(-1) * ( za(p1,p2)
     &    **3*za(p7,a7)*zb(p1,p4)*zb(p2,a7)*s12**(-2) - za(p1,p2)**2*
     &    za(p1,p7)*za(p2,r2)*zb(p1,p4)*zb(p2,r2)*s12**(-2) - za(p1,p2)
     &    *za(p1,p7)*za(p2,q4)*zb(p4,q4)*s12**(-1) - za(p1,p2)*za(p1,p7
     &    )*za(p2,r2)*zb(p4,r2)*s12**(-1) )
      mba(1,1) = mba(1,1) + mt**3*denr**(-1) * (  - za(p1,p2)**2*za(p1,
     &    p7)*zb(p1,p4)*s12**(-2) )
      mba(1,1) = mba(1,1) + mt**3*denb**(-1)*denq2**(-1) * ( za(p1,p2)
     &    **3*za(p7,a7)*zb(p1,p4)*zb(p2,a7)*s12**(-2) + za(p1,p2)**3*
     &    za(p7,b2)*zb(p1,p4)*zb(p2,b2)*s12**(-2) - za(p1,p2)**2*za(p1,
     &    p7)*za(p2,b2)*zb(p1,p4)*zb(p2,b2)*s12**(-2) - za(p1,p2)**2*
     &    za(p1,p7)*zb(p1,p4)*s12**(-1) - za(p1,p2)*za(p1,p7)*za(p2,q4)
     &    *zb(p4,q4)*s12**(-1) )
      mba(1,1) = mba(1,1) + mt**3*denb**(-1) * (  - za(p1,p2)**2*za(p1,
     &    p7)*zb(p1,p4)*s12**(-2) )
      mba(1,1) = mba(1,1) + mt**3*denq2**(-1)*dena1**(-1) * ( za(p1,p2)
     &    **3*za(p7,a7)*zb(p1,p4)*zb(p2,a7)*s12**(-2) - za(p1,p2)**2*
     &    za(p1,p7)*za(p2,r2)*zb(p1,p4)*zb(p2,r2)*s12**(-2) - za(p1,p2)
     &    **2*za(p1,p7)*za(p2,b2)*zb(p1,p4)*zb(p2,b2)*s12**(-2) - za(p1
     &    ,p2)**2*za(p1,p7)*zb(p1,p4)*s12**(-1) - za(p1,p2)*za(p1,p7)*
     &    za(p2,q4)*zb(p4,q4)*s12**(-1) )

      mba(2,2)= + mt*denr**(-1)*dena1**(-1) * ( za(p1,q4)*za(p2,r2)*za(
     &    p7,a7)*zb(p1,p2)**2*zb(p1,a7)*zb(p2,r2)*zb(p4,q4)*s12**(-2)
     &     - za(p1,r1)*za(p2,p7)*za(q4,r2)*zb(p1,p2)**2*zb(p1,r1)*zb(p2
     &    ,r2)*zb(p4,q4)*s12**(-2) + za(p1,r1)*za(p2,r2)*za(p7,a7)*zb(
     &    p1,p2)**2*zb(p1,a7)*zb(p2,r2)*zb(p4,r1)*s12**(-2) + za(p2,p7)
     &    *za(q4,r2)*zb(p1,p2)**2*zb(p2,r2)*zb(p4,q4)*s12**(-1) + za(p7
     &    ,a7)*za(q4,r2)*zb(p1,p2)*zb(p1,a7)*zb(p2,r2)*zb(p4,q4)*
     &    s12**(-1) )
      mba(2,2) = mba(2,2) + mt*denr**(-1) * (  - za(p1,p7)*za(q4,r1)*
     &    zb(p1,p2)**2*zb(p1,r1)*zb(p4,q4)*s12**(-2) + za(p1,q4)*za(p7,
     &    a7)*zb(p1,p2)**2*zb(p1,a7)*zb(p4,q4)*s12**(-2) + za(p1,r1)*
     &    za(p7,a7)*zb(p1,p2)**2*zb(p1,a7)*zb(p4,r1)*s12**(-2) )
      mba(2,2) = mba(2,2) + mt*denb**(-1)*denq2**(-1) * (  - za(p1,q4)*
     &    za(p2,b2)*za(p7,a7)*zb(p1,p2)**3*zb(p4,q4)*zb(a7,b2)*
     &    s12**(-2) + za(p1,q4)*za(p2,b2)*za(p7,a7)*zb(p1,p2)**2*zb(p1,
     &    a7)*zb(p2,b2)*zb(p4,q4)*s12**(-2) + za(p1,q4)*za(p2,b2)*za(p7
     &    ,b1)*zb(p1,p2)**2*zb(p1,b1)*zb(p2,b2)*zb(p4,q4)*s12**(-2) + 
     &    za(p1,q4)*za(p7,a7)*zb(p1,p2)**2*zb(p1,a7)*zb(p4,q4)*
     &    s12**(-1) + za(p1,q4)*za(p7,b1)*zb(p1,p2)**2*zb(p1,b1)*zb(p4,
     &    q4)*s12**(-1) - za(p1,b1)*za(p2,b2)*za(p7,a7)*zb(p1,p2)**2*
     &    zb(p1,b1)*zb(p2,p4)*zb(a7,b2)*s12**(-2) )
      mba(2,2) = mba(2,2) + mt*denb**(-1) * ( za(p1,q4)*za(p7,a7)*zb(p1
     &    ,p2)**2*zb(p1,a7)*zb(p4,q4)*s12**(-2) + za(p1,q4)*za(p7,b1)*
     &    zb(p1,p2)**2*zb(p1,b1)*zb(p4,q4)*s12**(-2) - za(p1,b1)*za(p7,
     &    a7)*zb(p1,p2)**2*zb(p1,p4)*zb(a7,b1)*s12**(-2) )
      mba(2,2) = mba(2,2) + mt*denq2**(-1)*dena1**(-1) * (  - za(p1,q4)
     &    *za(p1,r1)*za(p2,p7)*zb(p1,p2)**3*zb(p1,r1)*zb(p4,q4)*
     &    s12**(-2) + za(p1,q4)*za(p2,p7)*za(p2,b2)*zb(p1,p2)**3*zb(p2,
     &    b2)*zb(p4,q4)*s12**(-2) + za(p1,q4)*za(p2,p7)*za(r1,b2)*zb(p1
     &    ,p2)**2*zb(p1,r1)*zb(p2,b2)*zb(p4,q4)*s12**(-2) + za(p1,q4)*
     &    za(p2,p7)*zb(p1,p2)**3*zb(p4,q4)*s12**(-1) + za(p1,q4)*za(p2,
     &    r2)*za(p7,a7)*zb(p1,p2)**2*zb(p1,a7)*zb(p2,r2)*zb(p4,q4)*
     &    s12**(-2) + za(p1,q4)*za(p2,b2)*za(p7,a7)*zb(p1,p2)**2*zb(p1,
     &    a7)*zb(p2,b2)*zb(p4,q4)*s12**(-2) + za(p1,q4)*za(p7,a7)*zb(p1
     &    ,p2)**2*zb(p1,a7)*zb(p4,q4)*s12**(-1) - za(p1,b1)*za(p2,r2)*
     &    za(p7,a7)*zb(p1,p2)**2*zb(p1,a7)*zb(p2,p4)*zb(r2,b1)*
     &    s12**(-2) )
      mba(2,2) = mba(2,2) + mt**3*denr**(-1)*dena1**(-1) * (  - za(p1,
     &    q4)*za(p2,p7)*zb(p1,p2)**3*zb(p4,q4)*s12**(-2) - za(p1,r1)*
     &    za(p2,p7)*zb(p1,p2)**3*zb(p4,r1)*s12**(-2) - za(p1,r1)*za(p2,
     &    p7)*zb(p1,p2)**2*zb(p1,r1)*zb(p2,p4)*s12**(-2) + za(p2,p7)*
     &    zb(p1,p2)**2*zb(p2,p4)*s12**(-1) + za(p7,a7)*zb(p1,p2)*zb(p1,
     &    a7)*zb(p2,p4)*s12**(-1) )
      mba(2,2) = mba(2,2) + mt**3*denr**(-1) * (  - za(p1,p7)*zb(p1,p2)
     &    **2*zb(p1,p4)*s12**(-2) )
      mba(2,2) = mba(2,2) + mt**3*denb**(-1)*denq2**(-1) * (  - za(p1,
     &    q4)*za(p2,p7)*zb(p1,p2)**3*zb(p4,q4)*s12**(-2) - za(p1,b1)*
     &    za(p2,p7)*zb(p1,p2)**2*zb(p1,b1)*zb(p2,p4)*s12**(-2) + za(p7,
     &    a7)*zb(p1,p2)*zb(p1,a7)*zb(p2,p4)*s12**(-1) + za(p7,b1)*zb(p1
     &    ,p2)*zb(p1,b1)*zb(p2,p4)*s12**(-1) )
      mba(2,2) = mba(2,2) + mt**3*denb**(-1) * (  - za(p1,p7)*zb(p1,p2)
     &    **2*zb(p1,p4)*s12**(-2) )
      mba(2,2) = mba(2,2) + mt**3*denq2**(-1)*dena1**(-1) * (  - za(p1,
     &    q4)*za(p2,p7)*zb(p1,p2)**3*zb(p4,q4)*s12**(-2) - za(p1,r1)*
     &    za(p2,p7)*zb(p1,p2)**2*zb(p1,r1)*zb(p2,p4)*s12**(-2) - za(p1,
     &    b1)*za(p2,p7)*zb(p1,p2)**2*zb(p1,b1)*zb(p2,p4)*s12**(-2) + 
     &    za(p2,p7)*zb(p1,p2)**2*zb(p2,p4)*s12**(-1) + za(p7,a7)*zb(p1,
     &    p2)*zb(p1,a7)*zb(p2,p4)*s12**(-1) )

      mba(1,2)= + mt*denr**(-1)*dena1**(-1) * (  - za(p1,p7)*za(p1,r1)*
     &    za(q4,r2)*zb(p2,r1)*zb(p2,r2)*zb(p4,q4)*s12**(-1) + za(p1,q4)
     &    *za(p1,r2)*za(p7,a7)*zb(p2,a7)*zb(p2,r2)*zb(p4,q4)*s12**(-1)
     &     + za(p1,r1)*za(p1,r2)*za(p7,a7)*zb(p2,a7)*zb(p2,r2)*zb(p4,r1
     &    )*s12**(-1) )
      mba(1,2) = mba(1,2) + mt*denb**(-1)*denq2**(-1) * ( za(p1,q4)*za(
     &    p1,b2)*za(p7,a7)*zb(p2,a7)*zb(p2,b2)*zb(p4,q4)*s12**(-1) + 
     &    za(p1,q4)*za(p1,b2)*za(p7,b2)*zb(p2,b2)**2*zb(p4,q4)*
     &    s12**(-1) - za(p1,b1)**2*za(p7,a7)*zb(p2,p4)*zb(p2,b1)*zb(a7,
     &    b1)*s12**(-1) )
      mba(1,2) = mba(1,2) + mt*denq2**(-1)*dena1**(-1) * (  - za(p1,p2)
     &    *za(p1,b1)*za(p7,a7)*zb(p2,p4)*zb(p2,a7)*zb(p2,b1)*s12**(-1)
     &     - za(p1,p7)*za(p1,q4)*za(p1,r1)*zb(p1,p2)*zb(p2,r1)*zb(p4,q4
     &    )*s12**(-1) + za(p1,p7)*za(p1,q4)*za(r2,b2)*zb(p2,r2)*zb(p2,
     &    b2)*zb(p4,q4)*s12**(-1) + za(p1,q4)*za(p1,r2)*za(p7,a7)*zb(p2
     &    ,a7)*zb(p2,r2)*zb(p4,q4)*s12**(-1) + za(p1,q4)*za(p1,b2)*za(
     &    p7,a7)*zb(p2,a7)*zb(p2,b2)*zb(p4,q4)*s12**(-1) - za(p1,r1)*
     &    za(p1,b1)*za(p7,a7)*zb(p2,p4)*zb(p2,a7)*zb(r1,b1)*s12**(-1) )
      mba(1,2) = mba(1,2) + mt**3*denr**(-1)*dena1**(-1) * (  - za(p1,
     &    p7)*za(p1,r1)*zb(p2,p4)*zb(p2,r1)*s12**(-1) )
      mba(1,2) = mba(1,2) + mt**3*denb**(-1)*denq2**(-1) * (  - za(p1,
     &    p7)*za(p1,b1)*zb(p2,p4)*zb(p2,b1)*s12**(-1) )
      mba(1,2) = mba(1,2) + mt**3*denq2**(-1)*dena1**(-1) * (  - za(p1,
     &    p7)*za(p1,r1)*zb(p2,p4)*zb(p2,r1)*s12**(-1) - za(p1,p7)*za(p1
     &    ,b1)*zb(p2,p4)*zb(p2,b1)*s12**(-1) )

      mba(2,1)= + mt*denr**(-1)*dena1**(-1) * (  - za(p2,p7)*za(p2,r2)*
     &    za(q4,r1)*zb(p1,r1)*zb(p1,r2)*zb(p4,q4)*s12**(-1) + za(p2,q4)
     &    *za(p2,r1)*za(p7,a7)*zb(p1,a7)*zb(p1,r1)*zb(p4,q4)*s12**(-1)
     &     + za(p2,r1)*za(p2,r2)*za(p7,a7)*zb(p1,a7)*zb(p1,r1)*zb(p4,r2
     &    )*s12**(-1) )
      mba(2,1) = mba(2,1) + mt*denb**(-1)*denq2**(-1) * ( za(p2,q4)*za(
     &    p2,b1)*za(p7,a7)*zb(p1,a7)*zb(p1,b1)*zb(p4,q4)*s12**(-1) + 
     &    za(p2,q4)*za(p2,b1)*za(p7,b1)*zb(p1,b1)**2*zb(p4,q4)*
     &    s12**(-1) - za(p2,b2)**2*za(p7,a7)*zb(p1,p4)*zb(p1,b2)*zb(a7,
     &    b2)*s12**(-1) )
      mba(2,1) = mba(2,1) + mt*denq2**(-1)*dena1**(-1) * ( za(p1,p2)*
     &    za(p2,r1)*za(p7,a7)*zb(p1,p4)*zb(p1,a7)*zb(p1,r1)*s12**(-1)
     &     + za(p2,p7)*za(p2,q4)*za(p2,b1)*zb(p1,p2)*zb(p1,b1)*zb(p4,q4
     &    )*s12**(-1) + za(p2,p7)*za(p2,q4)*za(r1,b1)*zb(p1,r1)*zb(p1,
     &    b1)*zb(p4,q4)*s12**(-1) + za(p2,q4)*za(p2,r1)*za(p7,a7)*zb(p1
     &    ,a7)*zb(p1,r1)*zb(p4,q4)*s12**(-1) + za(p2,q4)*za(p2,b1)*za(
     &    p7,a7)*zb(p1,a7)*zb(p1,b1)*zb(p4,q4)*s12**(-1) - za(p2,r2)*
     &    za(p2,b2)*za(p7,a7)*zb(p1,p4)*zb(p1,a7)*zb(r2,b2)*s12**(-1) )
      mba(2,1) = mba(2,1) + mt**3*denr**(-1)*dena1**(-1) * (  - za(p2,
     &    p7)*za(p2,r2)*zb(p1,p4)*zb(p1,r2)*s12**(-1) )
      mba(2,1) = mba(2,1) + mt**3*denb**(-1)*denq2**(-1) * (  - za(p2,
     &    p7)*za(p2,b2)*zb(p1,p4)*zb(p1,b2)*s12**(-1) )
      mba(2,1) = mba(2,1) + mt**3*denq2**(-1)*dena1**(-1) * (  - za(p2,
     &    p7)*za(p2,r2)*zb(p1,p4)*zb(p1,r2)*s12**(-1) - za(p2,p7)*za(p2
     &    ,b2)*zb(p1,p4)*zb(p1,b2)*s12**(-1) )

      wtgg=0d0
      do h1=1,2
      do h2=1,2
      wtgg=wtgg+abs(mab(h1,h2))**2+abs(mba(h1,h2))**2
     & -1d0/xn**2*abs(mab(h1,h2)+mba(h1,h2))**2
      enddo
      enddo
      return
      end
