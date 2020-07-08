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
 
      function AGTYG1s(t,u,Lx,Ly)
      include 'types.f'
!     Results taken from hep-ph/0201274
!  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
!     Eq. B.1
      real(dp):: AGTYG1s,G1sx,t,u,Lx,Ly
      AGTYG1s=+G1sx(t,u,Lx,Ly)+G1sx(u,t,Ly,Lx)
      return
      end

      function G1sx(t,u,Lx,Ly)
      implicit none
      include 'types.f'
      real(dp):: G1sx
      include 'constants.f'
      include 'zeta.f'
      real(dp)::t,u,Lx,Ly

      G1sx=
     & (14*Lx**4+28*Lx**3+8*Lx**2*Ly**2+56*Lx**2*pisq-48*Lx**2
     & +12*Lx**2*Ly+32*Lx*Ly*pisq+80*pisq*Lx
     & +2*Ly**4+12*Ly**3-10*Ly**2+8*Ly**2*pisq+26*pisq+24*pisq*Ly
     * -84*Ly+102._dp)*(t/u)  
     & +8*Lx*(Lx**3+Lx**2+4*pisq*Lx+2*pisq)*(t/u)**2
     & +2*Lx**2*(Lx**2+4*pisq)*(t/u)**3  
     & +(32*Lx**3+8*Lx**2*Ly**2+80*pisq*Lx+32*Lx*Ly*pisq+8*Ly**2*Lx
     & +8*Ly**4+32*Ly**2*pisq-32*Ly**2-4._dp-56*Ly+24*pisq)
!     + \t \leftrightarrow u \ 
      return
      end
