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
 
      subroutine qqbZtt1(propz,wtqqb,wtqbq)
      implicit none
      include 'types.f'
      
C***********************************************************************
C     Author: R.K. Ellis                                               *
C     December, 2010.                                                  *
C     calculate the Born matrix element squared                        *
C     for the process                                                  *
c     My notation                                                      *
C     qb(-p1) + q(-p2)=bbar(p6)+e-(p7)+nubar(p8)+b(p5)+nu(p3)+e+(p4)   *
C                                                                      *
C q=3+4+5                                                              *
C a=-6-7-8                                                             *
C since momenta 3,5,6,8 are not needed they are reused to              *
C represent the de-massified vectors.                                  *
C Thus q4 = q de-massified wrt p4 ;                                    *
C Thus a7 = a de-massified wrt p7 etc;                                 *
C this vector has been stored in wrapper routine (qqbZtt) in p(3,mu)   *
C q4(mu)=q(mu)-qsq/2/p4Dq*p4(mu)-->p(3,mu)                             *
C a7(mu)=a(mu)-asq/2/p7Da*p7(mu)-->p(5,mu)                             *
C Formula generated by the program qqbZtt1.frm                         *
C First index of m is the qqb polarization                             *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'zprods_com.f'
      include 'masses.f'
      include 'topzlabels.f'
      complex(dp):: mqqb(2),mqbq(2),propz
      real(dp):: wtqqb(2),wtqbq(2),mtsq
      integer:: h1,j,j12
C For charge and couplings of top quark, use up-charge
      j=2
      mtsq=mt**2
      do j12=1,2
      wtqbq(j12)=0._dp
      wtqqb(j12)=0._dp
      enddo
      do j12=1,2
      mqbq(1)=(Q(j)*Q(j12)+L(j12)*R(j)*propz)*za(p1,p7)*zb(p2,p4)*mtsq
     & -(Q(j)*Q(j12)+L(j)*L(j12)*propz)
     & *za(p1,q4)*za(p7,a7)*zb(p2,a7)*zb(p4,q4)

      mqbq(2)=(Q(j)*Q(j12)+R(j)*R(j12)*propz)*za(p2,p7)*zb(p1,p4)*mtsq
     & -(Q(j)*Q(j12)+L(j)*R(j12)*propz)
     & *za(p2,q4)*za(p7,a7)*zb(p1,a7)*zb(p4,q4)

      mqqb(1)=(Q(j)*Q(j12)+L(j12)*R(j)*propz)*za(p2,p7)*zb(p1,p4)*mtsq 
     & -(Q(j)*Q(j12)+L(j)*L(j12)*propz)
     & *za(p2,q4)*za(p7,a7)*zb(p1,a7)*zb(p4,q4)

      mqqb(2)=(Q(j)*Q(j12)+R(j)*R(j12)*propz)*za(p1,p7)*zb(p2,p4)*mtsq 
     & -(Q(j)*Q(j12)+L(j)*R(j12)*propz)
     & *za(p1,q4)*za(p7,a7)*zb(p2,a7)*zb(p4,q4)

      do h1=1,2
      wtqqb(j12)=wtqqb(j12)+abs(mqqb(h1))**2
      wtqbq(j12)=wtqbq(j12)+abs(mqbq(h1))**2
      enddo
      enddo
      return
      end
