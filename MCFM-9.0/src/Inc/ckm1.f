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
 
      real(dp):: VV(-nf:nf,-nf:nf),gl(-nf:nf,-nf:nf),
     & gr(-nf:nf,-nf:nf),e(-nf:nf,-nf:nf),fl,fr
      real(dp):: glsq(-nf:nf,-nf:nf),grsq(-nf:nf,-nf:nf),
     & flsq,frsq
      common/CKM1/VV,gl,gr,fl,fr,e,glsq,grsq,flsq,frsq
