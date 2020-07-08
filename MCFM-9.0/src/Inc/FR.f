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
 
!---------------------------------------------
! FR contains regular and plus pieces associated 
! with radiation from of a photon by a quark qu,or antiquark
! qub 
      integer:: quf,qubf
      parameter(quf=1,qubf=-1) 

      real(dp):: FR(-1:1,-1:1,1:3) 

      common/FR/FR
