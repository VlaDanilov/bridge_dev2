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
 
      subroutine mcfmmain(r,er)
          use Integration
          use parseinput
      implicit none
************************************************************************
*                                                                      *
*  This is the main program for MCFM                                   *
*                                                                      *
*  The sequence of calls should always be:                             *
*   call mcfm_init          : basic variable initialization, print-out *
*   call mcfm_vegas(warmup) : warm-up the Vegas grid                   *
*   call mcfm_vegas(accum)  : accumulate results                       *
*   call mcfm_exit          : final processing and print-out           *
*                                                                      *
************************************************************************
      
      real(dp):: integ,integ_err,r,er
      logical:: dryrun
      common/dryrun/dryrun

      logical :: forceAdaptive, forceOld

      integer*4 old_cw

      call f_fpu_fix_start(old_cw)

* basic variable initialization, print-out
      call mcfm_init()

      call mcfm_vegas_adaptive(integ,integ_err)

* final processing and print-out
      r=integ
      er=integ_err
      call mcfm_exit(integ,integ_err)
      end
       
