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
 
c--- suite of dummy routines that replace OneLOop routines if
c--- that library is not linked
      subroutine error_olo_notlinked
      implicit none
      include 'types.f'
      
      write(6,*) 'OneLOop library has not been linked in'
      write(6,*) 'the MCFM makefile. Please recompile with'
      write(6,*) 'appropriate flag set.'  
      stop
          
      return
      end
      

      subroutine olo_unit(i1,aa)
      implicit none
      include 'types.f'
      
      integer:: i1
      character*(*) aa
      
      call error_olo_notlinked()
      
      return
      end
      
      subroutine olo_onshell(x1)
      implicit none
      include 'types.f'
      
      real(dp):: x1
      
      call error_olo_notlinked()
      
      return
      end
      
      subroutine olo_a0(c0,x1,x2)
      implicit none
      include 'types.f'
      
      real(dp):: x1,x2
      complex(dp):: c0(0:2)
      
      call error_olo_notlinked()
      
      return
      end
      
      subroutine olo_b0(c0,x1,x2,x3,x4)
      implicit none
      include 'types.f'
      
      real(dp):: x1,x2,x3,x4
      complex(dp):: c0(0:2)
      
      call error_olo_notlinked()
      
      return
      end
      
      subroutine olo_c0(c0,x1,x2,x3,x4,x5,x6,x7)
      implicit none
      include 'types.f'
      
      real(dp):: x1,x2,x3,x4,x5,x6,x7
      complex(dp):: c0(0:2)
      
      call error_olo_notlinked()
      
      return
      end
      
      subroutine olo_d0(c0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11)
      implicit none
      include 'types.f'
      
      real(dp):: x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11
      complex(dp):: c0(0:2)
      
      call error_olo_notlinked()
      
      return
      end
      
      
