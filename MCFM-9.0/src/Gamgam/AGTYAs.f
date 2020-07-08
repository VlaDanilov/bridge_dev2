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
 
      function AGTYAs(s,t,u,Lx,Ly,Li2x,Li2y,Li3x,Li3y,Li4x,Li4y,
     & Li4z,Li4zinv)
!     Results taken from hep-ph/0201274
!  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
!     Eq. A.1
      implicit none
      include 'types.f'
      real(dp):: AGTYAs,Asx
      real(dp)::s,t,u,Lx,Ly,Li2x,Li2y,Li3x,Li3y,Li4x,Li4y,Li4z,Li4zinv
      AGTYAs=Asx(s,t,u,Lx,Ly,Li2x,Li3x,Li3y,Li4x,Li4y,Li4z)
     &      +Asx(s,u,t,Ly,Lx,Li2y,Li3y,Li3x,Li4y,Li4x,Li4zinv)
      return
      end

      function Asx(s,t,u,Lx,Ly,Li2x,Li3x,Li3y,Li4x,Li4y,Li4z)
      implicit none
      include 'types.f'
      real(dp):: Asx
      include 'constants.f'
      include 'zeta.f'
      real(dp)::s,t,u,Lx,Ly,Li4x,Li4y,Li4z,Li3x,Li3y,Li2x

       Asx=(128*Li4z-128*Li4x+128*Li4y+
     & (-64._dp/3._dp+128*Ly )*Li3x  
     & +(64._dp/3._dp*Lx-64._dp/3._dp*pisq )*Li2x+16._dp/3._dp*Lx**4
     & -64._dp/3._dp*Lx**3*Ly+
     & (-16._dp+32._dp/3._dp*pisq+32*Ly**2 )*Lx**2  
     & +(-64._dp/3._dp*pisq*Ly+48._dp+160._dp/9._dp*pisq)*Lx
     & +64._dp/3._dp*zeta3+224._dp/45._dp*pi**4-128*Ly*zeta3)*(t/u)  
     & +(32._dp/3._dp*Li3x-32._dp/3._dp*Li3y+
     & (-32._dp/3._dp*Lx-32._dp/3._dp*Ly )*Li2x  
     & +(-32._dp/3._dp*Ly**2-80._dp/9._dp*pisq-64._dp/3._dp)*Lx
     & +(64._dp/3._dp+32._dp/3._dp*pisq )*Ly)*(t/s)**2
     & +24*Lx**2*(t/u)**2  +(416._dp/3._dp*Li3x
     & +64*Li3y*Lx-416._dp/3*Li2x*Lx
     & +(8*Ly**2+16._dp)*Lx**2  
     & +(-8._dp/3._dp*Ly+80._dp/3._dp
     & +112._dp/9._dp*pisq-64*zeta3-64*Ly**2 )*Lx
     & -416._dp/3._dp*zeta3-148._dp/9._dp*pisq+44._dp/45._dp*pi**4 ) 
!     & & & + \t \leftrightarrow u \ 
      return
      end
