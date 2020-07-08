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
 
      function AGTYX2s(t,u,Lx,Ly,Ls)
      implicit none
!     Results taken from hep-ph/0201274
!  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
!     Eq. B.5
      include 'types.f'
      real(dp):: AGTYX2s,X2sx,t,u,Lx,Ly,Ls
      AGTYX2s=+X2sx(t,u,Lx,Ly,Ls)+X2sx(u,t,Ly,Lx,Ls)
      return
      end

      function X2sx(t,u,Lx,Ly,Ls)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'zeta.f'
      real(dp):: X2sx
      real(dp)::t,u,Lx,Ly,Ls

      X2sx=-1._dp/9._dp*(2*Ls+Ly+Lx)
     & *(11*Ls+3*Lx**2+6*Lx*Ly+10*Ly+10*Lx-3*pisq)*(t/u)  
     & + (-2._dp/3._dp*Lx*Ly-2._dp/3._dp*Ly**2-2._dp/
     & 3*Lx**2*Ly-4._dp/3._dp*Ly**2*Ls-4._dp/
     & 3*Lx*Ls-2._dp/3._dp*Ly**3) 
! + \t \leftrightarrow u \ 
      return
      end
