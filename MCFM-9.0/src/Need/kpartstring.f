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
 
c--- converts integer kpart to corresponding 4-character string
      function kpartstring(k)
      implicit none
      include 'kpart.f'
      character*15 kpartstring
      integer k

      if     (k == klord) then
        kpartstring='lo'
      elseif (k == kvirt) then
        kpartstring='virt'
      elseif (k == kreal) then
        kpartstring='real'
      elseif (k == ktota) then
        kpartstring='nlo'
      elseif (k == kfrag) then
        kpartstring='frag'
      elseif (k == ktodk) then
        kpartstring='todk'
      elseif (k == ksnlo) then
        kpartstring='snlo'
      elseif (k == knnlo) then
        kpartstring='nnlo'
      else
        write(6,*) 'Unexpected kpart in kpartstring: ',k
        stop
      endif

      if (coeffonly) then
        kpartstring=trim(kpartstring)//'coeff'
      endif

      return
      end

      
