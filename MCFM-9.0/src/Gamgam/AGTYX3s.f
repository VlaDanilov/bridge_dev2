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
 
      function AGTYX3s(t,u,Lx,Ly,Ls)
      implicit none
!     Results taken from hep-ph/0201274
!  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
!     Eq. B.6
      include 'types.f'
      real(dp):: AGTYX3s,X3sx,t,u,Lx,Ly,Ls
      AGTYX3s=+X3sx(t,u,Lx,Ly,Ls)+X3sx(u,t,Ly,Lx,Ls)
      return
      end

      function X3sx(t,u,Lx,Ly,Ls)
      implicit none
      include 'types.f'
      real(dp):: X3sx
      real(dp)::t,u,Lx,Ly,Ls

      X3sx=1._dp/18._dp*(2*Ls+Ly+Lx)**2*(t/u) 
! + \t \leftrightarrow u \ 
      return
      end
