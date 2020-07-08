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
 
      real(dp):: xmass,medmass,dm_lam
      real(dp):: medwidth
      real(dp):: gdm,g_dmx,g_dmq
      logical:: effective_th
      character*6 dm_mediator 
      real(dp):: dmL(5),dmR(5)
      logical:: yukawa_scal
      common/yuk_scal/yukawa_scal
      common/dm_params/xmass,medmass,dm_lam,medwidth
      common/dm_coup/dmL,dmR
      common/dm_med/dm_mediator
      common/effec_dm/effective_th
      common/dm_g/gdm,g_dmx,g_dmq
