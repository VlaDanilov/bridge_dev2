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
 
      function getdynamictau(p)
      implicit none
      include 'types.f'
      include 'kprocess.f'
      include 'mxpart.f'
      include 'taucut.f'
      real(dp):: getdynamictau,p(mxpart,4),pt,pttwo

! Return normal value of tau if not(dynamictau)
      if (dynamictau .eqv. .false.) then
        getdynamictau=taucut
        return
      endif
      
      if     ( (kcase == kW_only) .or. (kcase == kW_1jet)
     &     .or.(kcase == kZ_only) .or. (kcase == kZ_1jet)
     &     .or.(kcase == kggfus0) .or. (kcase == kggfus1)
     &     .or.(kcase == kgamgam) .or. (kcase == kgmgmjt) ) then
        getdynamictau=sqrt((p(3,4)+p(4,4))**2-(p(3,1)+p(4,1))**2
     &                    -(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2)
      elseif ( (kcase == kZgamma) .or. (kcase == kZgajet) ) then
        getdynamictau=sqrt((p(3,4)+p(4,4)+p(5,4))**2-(p(3,1)+p(4,1)+p(5,1))**2
     &                    -(p(3,2)+p(4,2)+p(5,2))**2-(p(3,3)+p(4,3)+p(5,3))**2)
      elseif ( (kcase == kWHbbar) .or. (kcase == kWH1jet) 
     &     .or.(kcase == kZHbbar) .or. (kcase == kZH1jet) ) then
        getdynamictau=sqrt((p(3,4)+p(4,4)+p(5,4)+p(6,4))**2-(p(3,1)+p(4,1)+p(5,1)+p(6,1))**2
     &                    -(p(3,2)+p(4,2)+p(5,2)+p(6,2))**2-(p(3,3)+p(4,3)+p(5,3)+p(6,3))**2)
      else
        write(6,*) 'Invalid process for dynamic taucut!'
        stop
      endif
      
      getdynamictau=getdynamictau*taucut
      
      return
      end
      
