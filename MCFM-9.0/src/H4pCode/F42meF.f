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
 
      function F42meF(psq,qsq,s,t)
      implicit none
      include 'types.f'
      complex(dp):: F42meF
c--- note: ordering of arguments to function is taken from, e.g.
c---       arXiV:0804.4149v3 (App. B) and not hep-ph/0607139 (Eq. 22)
      
c      include 'scale.f'
      real(dp):: s,t,psq,qsq
c      real(dp):: den
c      complex(dp):: qlI4,qlI3
      complex(dp):: Lsm1_2me
c      den=s*t-psq*qsq
c--- note: added a factor of 1/2 here
c      F42meF=
c     & +den*qlI4(0._dp,psq,0._dp,qsq,s,t,0._dp,0._dp,0._dp,0._dp,musq,0)/2._dp
c     & -  s*qlI3(0._dp,0._dp,s,0._dp,0._dp,0._dp,musq,0)
c     & -  t*qlI3(0._dp,0._dp,t,0._dp,0._dp,0._dp,musq,0)
c     & +psq*qlI3(0._dp,0._dp,psq,0._dp,0._dp,0._dp,musq,0)
c     & +qsq*qlI3(0._dp,0._dp,qsq,0._dp,0._dp,0._dp,musq,0)
      
c--- NOTE: checked on 8/30/09 that Lsm1_2me == (expression above)
      F42meF=Lsm1_2me(s,t,psq,qsq)
     
      return
      end

