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
 
      function Fcc(st,j1,j2,j3,j4,j5,j6,za,zb) 
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: Fcc
      integer:: j1,j2,j3,j4,j5,j6
      character*9 st
      complex(dp):: Fcc_qpgmgpqm,Fcc_qpgpgmqm,Fcc_qpgpgpqm,
     & Fcc_qpgpqmgm,Fcc_qpgpqmgp,Fcc_qpqmgmgp,Fcc_qpqmgpgm,Fcc_qpqmgpgp

      if     (st=='q+g-g+qb-') then
        Fcc=Fcc_qpgmgpqm(j1,j2,j3,j4,j5,j6,za,zb)
      elseif (st=='q+g+g-qb-') then 
        Fcc=Fcc_qpgpgmqm(j1,j2,j3,j4,j5,j6,za,zb)
      elseif (st=='q+g+g+qb-') then
        Fcc=Fcc_qpgpgpqm(j1,j2,j3,j4,j5,j6,za,zb)
      elseif (st=='q+g+qb-g-') then
        Fcc=Fcc_qpgpqmgm(j1,j2,j3,j4,j5,j6,za,zb)
      elseif (st=='q+g+qb-g+') then
        Fcc=Fcc_qpgpqmgp(j1,j2,j3,j4,j5,j6,za,zb)
      elseif (st=='q+qb-g-g+') then
        Fcc=Fcc_qpqmgmgp(j1,j2,j3,j4,j5,j6,za,zb)
      elseif (st=='q+qb-g+g-') then
        Fcc=Fcc_qpqmgpgm(j1,j2,j3,j4,j5,j6,za,zb)
      elseif (st=='q+qb-g+g+') then
        Fcc=Fcc_qpqmgpgp(j1,j2,j3,j4,j5,j6,za,zb)
      else
        write(6,*) 'Error in Fcc: argument st is ',st
        stop
      endif

      return
      end
      
