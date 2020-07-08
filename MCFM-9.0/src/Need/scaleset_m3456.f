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
 
      subroutine scaleset_m3456(p,mu0)
      implicit none
      include 'types.f'
c--- subroutine to calculate dynamic scale equal to
c--- invariant mass of particles 3, 4, 5 and 6
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'kprocess.f'
      real(dp):: p(mxpart,4),mu0

      if((kcase==kWWqqbr) .or.
     &   (kcase==kWWnpol) .or.
     &   (kcase==kWW_jet) .or.
     &   (kcase==kWZbbar) .or.
     &   (kcase==kZZlept) .or.
     &   (kcase==kWHbbar) .or.
     &   (kcase==kWHbbdk) .or.
     &   (kcase==kWHgaga) .or.
     &   (kcase==kWH1jet) .or.
     &   (kcase==kZH1jet) .or.
     &   (kcase==kZHbbar) .or.
     &   (kcase==kZHgaga) .or.
     &   (kcase==kHWW_4l) .or.
     &   (kcase==kHWW_tb) .or.
     &   (kcase==kHWWint) .or.
     &   (kcase==kHWWHpi) .or.
     &   (kcase==kggWW4l) .or.
     &   (kcase==kggWWbx) .or.
     &   (kcase==kHZZ_4l) .or.
     &   (kcase==kHZZ_tb) .or.
     &   (kcase==kHZZint) .or.
     &   (kcase==kHZZHpi) .or.
     &   (kcase==kHZZqgI) .or.
     &   (kcase==kHZZpjt) .or.
     &   (kcase==kHWWjet) .or.
     &   (kcase==kHZZjet) .or.
     &   (kcase==kHWW2jt) .or.
     &   (kcase==kHZZ2jt) .or.
     &   (kcase==kqq_HZZ) .or.
     &   (kcase==kfourga) .or.
     &   (kcase==kW_2gam) .or.
     &   (kcase==kZ_2gam) .or.
     &   (kcase==kHmZZ4l) .or.
     &   (kcase==kHVV_tb) .or.
     &   (kcase==kggVV4l) .or.
     &   (kcase==kggVVbx) .or.
     &   (kcase==kggZZ4l) .or.
     &   (kcase==kggZZbx) .or.
     &   (kcase==kHmZZ4l) .or. 
     &   (kcase==kHHpair)) then 
        mu0=(p(3,4)+p(4,4)+p(5,4)+p(6,4))**2
     &     -(p(3,1)+p(4,1)+p(5,1)+p(6,1))**2
     &     -(p(3,2)+p(4,2)+p(5,2)+p(6,2))**2
     &     -(p(3,3)+p(4,3)+p(5,3)+p(6,3))**2       
        mu0=sqrt(abs(mu0))
      else
        write(6,*)'dynamicscale m(3456) not supported for this process.'
        stop
      endif
      
      return
      end
      
