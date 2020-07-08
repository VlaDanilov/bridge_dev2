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
 
c--- ancillary array that is used for fragmentation
c--- contributions in gmgmjt process
      real(dp):: msqbits(12)
      common/msqbits/msqbits
      integer::
     & ddb_ddb,ddb_ssb,ddb_uub,
     & uub_uub,uub_ccb,uub_ddb,
     & dbd_ddb,dbd_ssb,dbd_uub,
     & ubu_uub,ubu_ccb,ubu_ddb
      parameter(
     & ddb_ddb=1,ddb_ssb=2,ddb_uub=3,
     & uub_uub=4,uub_ccb=5,uub_ddb=6,
     & dbd_ddb=7,dbd_ssb=8,dbd_uub=9,
     & ubu_uub=10,ubu_ccb=11,ubu_ddb=12)
     
!$omp threadprivate(/msqbits/)     
     
