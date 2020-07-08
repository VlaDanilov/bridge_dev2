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
 
      function AGTYX1s(t,u,Lx,Ly,Ls)
      implicit none
!     Results taken from hep-ph/0201274
!  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
!     Eq. B.4
      include 'types.f'
      real(dp):: AGTYX1s,X1sx,t,u,Lx,Ly,Ls
      AGTYX1s=X1sx(t,u,Lx,Ly,Ls)+X1sx(u,t,Ly,Lx,Ls)
      return
      end

      function X1sx(t,u,Lx,Ly,Ls)
      implicit none
      include 'types.f'
      real(dp):: X1sx
      real(dp)::t,u,Lx,Ly,Ls

      X1sx=2._dp/3._dp*(3*Ly+Ly**2-7._dp+2*Lx**2)*(2*Ls+Ly+Lx )*(t/u)  
     & +(4._dp/3._dp*Lx**2*Ly+8._dp/3._dp*Lx*Ls+4._dp/3*Lx*Ly
     & +4._dp/3._dp*Ly**3+4._dp/3._dp*Ly**2+8._dp/3._dp*Ly**2*Ls)

! + \t \leftrightarrow u \ 
      return
      end
