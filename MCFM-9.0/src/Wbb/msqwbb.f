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
 
      function msqwbb(i1,i2,i5,i6)
      implicit none
      include 'types.f'
      real(dp):: msqwbb
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'ewcouple.f'
C---This is the matrix element squared for 
c    |q_L(p1)+q_L(p5) --> q_L(p2)+q_L(pfive)+W(l(p6)+antilepton(p7))|^2
c   +|q_L(p1)+q_R(p5) --> q_L(p2)+q_R(pfive)+W(l(p6)+antilepton(p7))|^2
c---with couplings for W included
      integer:: i1,i2,i5,i6
      complex(dp):: aqqb_wbb
      real(dp):: faclo
      faclo=V*gsq**2*gwsq**2*aveqq
      msqwbb=faclo
     &*(abs(aqqb_wbb(i1,i2,i5,i6,3,4))**2
     & +abs(aqqb_wbb(i1,i2,i6,i5,3,4))**2)
      return 
      end
