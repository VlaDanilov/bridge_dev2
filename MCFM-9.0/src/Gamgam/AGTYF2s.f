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
 
      function AGTYF2s(t,u,Ls)
      implicit none
!     Results taken from hep-ph/0201274
!  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
!     Eq. A.11
      include 'types.f'
      real(dp):: AGTYF2s,F2sx,t,u,Ls

      AGTYF2s=F2sx(t,u,Ls)
     &       +F2sx(u,t,Ls)

      return
      end


      function F2sx(t,u,Ls)
      implicit none
      include 'types.f'
      real(dp):: F2sx
      include 'zeta.f'
      real(dp):: t,u,Ls

      F2sx=(((-160*Ls)/27. + (32*Ls**2)/9. - (92*pisq)/27.)*t)/u 

!     + \t \leftrightarrow u \ 
      return
      end
