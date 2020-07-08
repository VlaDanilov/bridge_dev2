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
 
      subroutine checkorder(order)
      implicit none
      include 'types.f'
c--- checks the value of nproc and part against a list of processes that
c--- are calculable only at LO. If the calculation is not possible,
c--- writes an error message and aborts
      
      include 'frag.f'
      include 'kpart.f'
      include 'nproc.f'
      character*1 order

c--- special cases where there is no LO calculation
      if ((kpart.ne.kreal) .and. (order == 'R')) then 
        write(6,*)
        write(6,*)'This process can only be calculated with part=real,'
        write(6,*)'it is a subset of a NLO calculation only.'
      stop
      endif
      
c--- if we're calculating LO only, there's no problem      
      if (kpart==klord) return
      
c--- otherwise, we must be performing a NLO calculation, and this list of
c--- process numbers can't be calculated beyond LO 
      if (order == 'L') then
        write(6,*)
        write(6,*)'This process cannot be calculated beyond LO - please'
        write(6,*)'check the values of nproc and part then try again'
        stop
      endif

c--- check that fragmentation is not turned on for a non-fragmentation process
      if ((frag) .and. (order .ne. 'F')) then
        write(6,*)
        write(6,*) 'This process does not include photon fragmentation.'
        write(6,*) 'Please set frag=.false. in the input file.'
        stop
      endif
       
      return
      end
  
