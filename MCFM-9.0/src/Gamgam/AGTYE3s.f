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
 
      function AGTYE3s(t,u,Lx,Ly,Ls)
      implicit none
!     Results taken from hep-ph/0201274
!  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
!     Eq. A.9
      include 'types.f'
      real(dp):: AGTYE3s,E3sx
      real(dp):: t,u,Lx,Ly,Ls
      AGTYE3s=E3sx(t,u,Lx,Ly,Ls)
     &       +E3sx(u,t,Ly,Lx,Ls)
      return
      end


      function E3sx(t,u,Lx,Ly,Ls)
      implicit none
      include 'types.f'
      real(dp):: E3sx
      include 'constants.f'
      include 'zeta.f'
      real(dp)::t,u,Lx,Ly,Ls

      E3sx=(16._dp/9*Lx**3+ (-76._dp/9+8._dp/3._dp*Ls)*Lx**2
     & +16._dp/9*pisq*Lx+8._dp/9*Ly**3
     & +(4._dp/3._dp*Ls-2._dp/9._dp)*Ly**2
     & +(8._dp/9*pisq+4*Ls-10._dp)*Ly  
     & -1._dp/3._dp*pisq*Ls-202._dp/27*Ls+19._dp/9._dp*pisq
     & -2._dp/9*zeta3+3401._dp/162._dp)*(t/u)  
     & + (16._dp/9*pisq*Lx+16._dp/9*Ly**3+
     & (8._dp/3._dp*Ls-52._dp/9._dp)*Ly**2  
     & + (-76._dp/9._dp+8._dp/3._dp*Ls )*Ly+8._dp/9._dp*pisq )
!    + \t \leftrightarrow u \ 
      return
      end

