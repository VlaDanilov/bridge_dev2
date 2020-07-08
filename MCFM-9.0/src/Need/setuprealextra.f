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
 
      subroutine setuprealextra(nprocextra)
      implicit none
! Routines to setup process number corresponding to additional real radiation
! that must also run when using current value of nproc
      include 'nproc.f'
      integer nprocextra
      
      if     (nproc == 92) then
        nprocextra = 920
      elseif (nproc == 97) then
        nprocextra = 970
      elseif (nproc == 101) then
        nprocextra = 1010
      elseif (nproc == 114) then
        nprocextra = 115
      elseif (nproc == 141) then
        nprocextra = 142
      elseif (nproc == 146) then
        nprocextra = 147
      elseif (nproc == 161) then
        nprocextra = 162
      elseif (nproc == 166) then
        nprocextra = 167
      elseif (nproc == 171) then
        nprocextra = 172
      elseif (nproc == 176) then
        nprocextra = 177
      elseif (nproc == 181) then
        nprocextra = 182
      elseif (nproc == 186) then
        nprocextra = 187
      elseif (nproc == 233) then
        nprocextra = 234
      elseif (nproc == 238) then
        nprocextra = 239
      elseif (nproc == 501) then
        nprocextra = 502
      elseif (nproc == 511) then
        nprocextra = 512
      else
        write(6,*) 'Unexpected process in setuprealextra: ',nproc
        stop
      endif
      
      return
      end
