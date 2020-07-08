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
 
      integer, parameter :: klord=1, kvirt=2, kreal=3, ktota=4, kfrag=5, ktodk=6
      integer, parameter :: ksnlo=7, knnlo=8
      integer, parameter :: knnloVV=1, knnloRV=2, knnloRR=3
      integer kpart,knnlopart
      integer origKpart
      logical coeffonly
      common/kpart/kpart
      common/origKpart/origKpart
      common/knnlopart/knnlopart
      common/coeffonly/coeffonly
      
