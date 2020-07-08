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
 
      integer :: lbjscheme
      common/lbjscheme_com/lbjscheme
      logical:: jetsopphem
      common /cjetsopphem/jetsopphem
      real(dp):: leptptmin,leptptmax,leptrapmin,leptrapmax,
     & misspt,Rjlmin,Rllmin,delyjjmin,
     & leptpt2min,leptpt2max,leptrap2min,leptrap2max,
     & gammptmin,gammptmax,gammrapmin,gammrapmax,
     & Rgalmin,mtrans34cut,
     & gammpt2,Rgagamin,gammpt3,Rgajetmin,
     & leptveto1min,leptveto1max,leptveto2min,leptveto2max 
      common/leptcuts/leptptmin,leptptmax,leptrapmin,leptrapmax,
     & misspt,Rjlmin,Rllmin,delyjjmin,
     & leptpt2min,leptpt2max,leptrap2min,leptrap2max,
     & gammptmin,gammptmax,gammrapmin,gammrapmax,
     & Rgalmin,mtrans34cut,
     & gammpt2,Rgagamin,gammpt3,Rgajetmin,
     & leptveto1min,leptveto1max,leptveto2min,leptveto2max 
