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
 
      function F42mhF(psq,qsq,s,t)
      implicit none
      include 'types.f'
      complex(dp):: F42mhF
      
c      include 'scale.f'
c      include 'epinv.f'
      real(dp):: s,t,psq,qsq
c      real(dp):: den
c      complex(dp):: qlI4,qlI3
      complex(dp):: Lsm1_2mht
c      integer:: ep

c      den=s*t
c--- note: added a factor of 1/2 here
c      F42mhF=
c     & +den*qlI4(0._dp,0._dp,psq,qsq,s,t,0._dp,0._dp,0._dp,0._dp,musq,0)/2._dp
c     & -  s*qlI3(0._dp,0._dp,s,0._dp,0._dp,0._dp,musq,0)/2._dp
c     & -  t*qlI3(0._dp,0._dp,t,0._dp,0._dp,0._dp,musq,0)
c     & +psq*qlI3(0._dp,0._dp,psq,0._dp,0._dp,0._dp,musq,0)/2._dp
c     & +qsq*qlI3(0._dp,0._dp,qsq,0._dp,0._dp,0._dp,musq,0)/2._dp

c--- NOTE: checked on 8/30/09 that Lsm1_2mht == (expression above)
      F42mhF=Lsm1_2mht(s,t,psq,qsq)
      
      return
      end

