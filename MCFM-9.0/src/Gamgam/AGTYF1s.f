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
 
      function AGTYF1s(t,u,Lx,Ly,Ls)
      implicit none
!     Results taken from hep-ph/0201274
!  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
!     Eq. A.7
      include 'types.f'
      real(dp):: AGTYF1s,F1sx,t,u,Lx,Ly,Ls

      AGTYF1s=F1sx(t,u,Lx,Ly,Ls)
     &       +F1sx(u,t,Ly,Lx,Ls)

      return
      end


      function F1sx(t,u,Lx,Ly,Ls)
      implicit none
      include 'types.f'
      real(dp):: F1sx
      include 'constants.f'
      include 'zeta.f'
      real(dp):: t,u,Lx,Ly,Ls

      F1sx= (5._dp/36._dp*Lx**2 
     & +(-10._dp/27._dp+1._dp/3._dp*Ls+1._dp/18._dp*Ly)*Lx
     & +5._dp/36._dp*Ly**2
     & +(-10._dp/27._dp+1._dp/3._dp*Ls)*Ly  
     & +1._dp/3._dp*Ls**2+1._dp/54._dp*pisq-20._dp/27._dp*Ls)*(t/u) 

!     + \t \leftrightarrow u \ 
      return
      end
