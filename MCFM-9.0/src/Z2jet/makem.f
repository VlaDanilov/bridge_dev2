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
 
      subroutine makem(i1,i2,i3,i4,i5,i6,i7,
     & m1_1234,m2_1234,m3_3412,m4_3412)
      implicit none
      include 'types.f'
C     Author: R.K. Ellis, March 2001
C     A subroutine calculating Nagy and Trocsnayi, PRD59 014020 (1999) 
C     Eq. A.23 with a factor of 2*i*e^2*g^3/s removed
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'sprods_com.f'
      integer:: fa,fb,ha,hb,hg,lh,i1,i2,i3,i4,i5,i6,i7,f(7),j
      real(dp):: vl(2),gq(2,2)
      complex(dp):: m1_1234(2,2,2,2,2,2),m2_1234(2,2,2,2,2,2),
     &               m3_3412(2,2,2,2,2,2),m4_3412(2,2,2,2,2,2),
     &               a1(2,2,2,2),a2(2,2,2,2),a3(2,2,2,2),a4(2,2,2,2),
     &               b1(2,2,2,2),b2(2,2,2,2),b3(2,2,2,2),b4(2,2,2,2),
     & prop
      include 'cplx.h'
      vl(1)=l1
      vl(2)=r1
      do j=1,2
      gq(j,1)=L(j)
      gq(j,2)=R(j)
      enddo

      call nagyqqQQg(i1,i2,i3,i4,i5,i6,i7,a1,a2,a3,a4)
      call nagyqqQQg(i3,i4,i1,i2,i5,i6,i7,b1,b2,b3,b4)
      prop=s(i6,i7)/cplx2((s(i6,i7)-zmass**2),zmass*zwidth)


      do fa=1,2
         f(i1)=fa
         f(i2)=fa
         do fb=1,2
            f(i3)=fb
            f(i4)=fb
            do ha=1,2
               do hb=1,2
                  do hg=1,2
                     do lh=1,2

      m1_1234(f(i1),f(i3),ha,hb,hg,lh)=
     & (Q(f(i1))*q1+gq(f(i1),ha)*vl(lh)*prop)*a1(ha,hb,hg,lh)
      m2_1234(f(i1),f(i3),ha,hb,hg,lh)=
     & (Q(f(i1))*q1+gq(f(i1),ha)*vl(lh)*prop)*a2(ha,hb,hg,lh)
      m3_3412(f(i3),f(i1),hb,ha,hg,lh)=
     & (Q(f(i3))*q1+gq(f(i3),hb)*vl(lh)*prop)*b3(hb,ha,hg,lh)
      m4_3412(f(i3),f(i1),hb,ha,hg,lh)=
     & (Q(f(i3))*q1+gq(f(i3),hb)*vl(lh)*prop)*b4(hb,ha,hg,lh)

                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      return
      end

