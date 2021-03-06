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
 
      function AGTYX4s(t,u,Lx,Ly,Ls)
      implicit none
!     Results taken from hep-ph/0201274
!  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
!     Eq. B.7
      include 'types.f'
      real(dp):: AGTYX4s,X4sx,t,u,Lx,Ly,Ls
      AGTYX4s=+X4sx(t,u,Lx,Ly,Ls)+X4sx(u,t,Ly,Lx,Ls)
      return
      end

      function X4sx(t,u,Lx,Ly,Ls)
      implicit none
      include 'types.f'
      include 'zeta.f'
      real(dp):: X4sx
      real(dp)::t,u,Lx,Ly,Ls

      X4sx= (-32._dp/3._dp*Ly*pisq+32._dp/
     & 3._dp*Ls*Lx**2+16._dp*Ls*Ly-64._dp/
     & 3._dp*Lx*pisq-16._dp*pisq+16._dp/3._dp*Ls*Ly**2-112._dp/
     & 3._dp*Ls)*t/u 
     & +(32._dp/3._dp*Ls*Lx**2+32._dp/3._dp*Ls*Ly-64._dp/
     & 3._dp*Lx*pisq-32._dp/3._dp*pisq)
!      + t leftrightarrow u 

      return
      end
