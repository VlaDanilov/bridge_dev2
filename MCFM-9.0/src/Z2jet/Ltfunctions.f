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
 
      function Ltm1(bx,by,x,y,xInt)
      implicit none
      include 'types.f'
      integer bx,by
      real(dp)::x,y
      complex(dp)::Ltm1,xInt(15)
      
      Ltm1=xInt(by)-xInt(bx)
      
      return
      end
      

      function Lt0(bx,by,x,y,xInt)
      implicit none
      include 'types.f'
      integer bx,by
      real(dp)::x,y
      complex(dp)::Lt0,Ltm1,xInt(15)
      
      Lt0=y/(y-x)*Ltm1(bx,by,x,y,xInt)
      
      return
      end
      

      function Lt1(bx,by,x,y,xInt)
      implicit none
      include 'types.f'
      integer bx,by
      real(dp)::x,y
      complex(dp)::Lt1,Lt0,xInt(15)
      
      Lt1=y/(y-x)*(Lt0(bx,by,x,y,xInt)+1._dp)
      
      return
      end
