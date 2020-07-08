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
 
      subroutine qqbgen(e1,p2,p3,e4,za,zb,zab,zba,ppmp,amp)
      implicit none
      include 'types.f'
      
C-----Authors: John Campbell and Keith Ellis, November 2011
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'zabprods_decl.f'
      complex(dp):: amp(2),ppmp
      integer:: e1,p2,p3,e4,e1p,e4p
C---we are calculating the amplitudes mmpp and mpmp
C---from our base amplitudes which are ppmp

c--- perform the swaps 5 <-> 7 and 6 <-> 8     
      e1p=12-e1
      e4p=14-e4
            
c      write(6,*) 'e1',e1
c      write(6,*) 'e4',e4
C---mpmp
      amp(1)=+ppmp(e1p,p2,p3,e4,za,zb,zab,zba)
c      amp(2)=ppmp(e1p,p3,p2,e4,za,zb,zab,zba)
c      write(6,*) 'amp(1)',amp(1),abs(amp(1))
c      write(6,*) 'amp(2)',amp(2),abs(amp(2))
C---mmpp
c--- JC: added minus sign here based on form of tree level amplitude
c---     (complex conjugation = (-1) * [za <-> zb])
      amp(2)=-ppmp(e1,p2,p3,e4p,zb,za,zba,zab)
c      write(6,*) 'amp(2)',amp(2),abs(amp(2))
c      write(6,*) 'rat',abs(amp(1))/abs(amp(2))
c      pause    
      return
      end
