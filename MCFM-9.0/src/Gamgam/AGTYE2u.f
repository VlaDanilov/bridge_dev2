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
 
      function AGTYE2u(s,t,Lx,Ly,Lu,Li2x,Li3x)
      implicit none
!     Results taken from hep-ph/0201274
!  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
!     Eq. A.17
      include 'types.f'
      real(dp):: AGTYE2u
      include 'constants.f'
      include 'zeta.f'
      real(dp)::s,t,u,Lx,Ly,Lu,Li2x,Li3x

      u=-s-t

      AGTYE2u= (4._dp/3._dp*Li3x-4._dp/
     & 3*Li2x*Lx-11._dp/36._dp*Lx**3+ (4._dp/
     & 3*Ly-13._dp/18-1._dp/2*Lu )*Lx**2  
     & + (-7._dp/2*Ly**2+ (1._dp/3._dp+3*Lu
     & )*Ly+1._dp/36._dp*pisq-31._dp/6*Lu+40._dp/9
     & )*Lx  
     & +7._dp/3._dp*Ly**3+ (-1._dp/3._dp-3*Lu
     & )*Ly**2+ (31._dp/3._dp*Lu-5._dp/9*pisq-80._dp/
     & 9 )*Ly  
     & -25._dp/9*zeta3-13._dp/9*pisq*Lu-11._dp/
     & 3*Lu**2+206._dp/27*Lu+65._dp/81+487._dp/108*pisq
     & )*(t**2+s**2)/(s*t)  
     & + (-11._dp/36._dp*Lx**3+ (14._dp/9-1._dp/2*Lu+Ly )*Lx**2 
     & +(-Ly**2+(-28._dp/9+Lu )*Ly-11._dp/36._dp*pisq-1._dp/3._dp)*Lx  
     & +pisq*Ly-1._dp/18*pisq
     & *(-28._dp+9*Lu ) )*(t**2-s**2)/(s*t)  
     & -Lx**2*(t**2+s**2)/u**2+2*Lx*(t**2-s**2)/u**2  
     & + (-11._dp/9*Lx**3+ (14._dp/9._dp-2*Lu+4*Ly )*Lx**2  
     & + (-14._dp/3._dp*Ly**2+ (-40._dp/9._dp+4*Lu )*Ly
     & +38._dp/9-2*Lu-11._dp/9*pisq )*Lx
     & +(-76._dp/9+4*Lu-8._dp/9*pisq )*Ly  
     & +28._dp/9*Ly**3+ (40._dp/9-4*Lu)*Ly**2-pisq*(-5._dp+2*Lu))
      return
      end

