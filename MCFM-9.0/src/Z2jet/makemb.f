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
 
      subroutine makemb(i1,i2,i3,i4,i5,i6,i7,mb1,mb2)
      implicit none
      include 'types.f'
      
C     Author: R.K. Ellis, March 2001
C     A subroutine calculating Nagy and Trocsanyi, PRD59 014020 (1999) 
C     Eq. B.56 with a factor of 2*i*e^2*g^3/s removed
      integer:: f1,f3,hq,Qh,hg,lh,i1,i2,i3,i4,i5,i6,i7
      complex(dp):: mb1(2,2,2,2,2,2),mb2(2,2,2,2,2,2),
     &               m1_1234(2,2,2,2,2,2),m2_1234(2,2,2,2,2,2),
     &               m3_3412(2,2,2,2,2,2),m4_3412(2,2,2,2,2,2)

      call makem(i1,i2,i3,i4,i5,i6,i7,m1_1234,m2_1234,m3_3412,m4_3412)
      do f1=1,2
      do f3=1,2
      do hq=1,2
      do Qh=1,2
      do hg=1,2
      do lh=1,2
      mb1(f1,f3,hq,Qh,hg,lh)=
     & m1_1234(f1,f3,hq,Qh,hg,lh)+m3_3412(f3,f1,Qh,hq,hg,lh)
      mb2(f1,f3,hq,Qh,hg,lh)=
     & m2_1234(f1,f3,hq,Qh,hg,lh)+m4_3412(f3,f1,Qh,hq,hg,lh)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
