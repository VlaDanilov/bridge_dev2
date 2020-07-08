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
 
      subroutine scaleset_ptj1(p,mu0)
c--- subroutine to calculate dynamic scale equal to
c---  ptj1, which is the pt of the leading jet in the event
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'constants.f'
      include 'kprocess.f'
      include 'jetlabel.f'
      include 'npart.f'
      include 'breit.f'
      integer isub,oldjets
      real(dp):: p(mxpart,4),pjet(mxpart,4),mu0,pt,rcut,ptj1,ptj2
      common/rcut/rcut

      if((kcase==ktwojet) .or.
     &   (kcase==ktwo_ew)) then

c--- first work out whether this point is real radiation or not
        if (abs(p(npart+2,4)) .gt. 1.e-8_dp) then
          isub=0  ! real term
        else
          isub=1  ! subtraction term
        endif
      
c-- cluster jets but make sure recorded number of jets is not changed
        oldjets=jets     
        call genclust2(p,rcut,pjet,isub)
      
c        write(6,*) 'partons:'
c        if (p(5,4) .ge. 1.e-8_dp) write(6,*) 'pt5',pt(5,p)
c        if (p(6,4) .ge. 1.e-8_dp) write(6,*) 'pt6',pt(6,p)  
c      write(6,*) 'jets:'
c        if (pjet(5,4) .ge. 1.e-8_dp) write(6,*) 'pt5',pt(5,pjet)
c        if (pjet(6,4) .ge. 1.e-8_dp) write(6,*) 'pt6',pt(6,pjet)
c        write(6,*)
      
c--- restore old value of jets
        jets=oldjets

c--- asssign scale
        mu0=max(pt(3,pjet),pt(4,pjet),pt(5,pjet))
c        write(6,*) mu0,pt(3,pjet),pt(4,pjet),pt(5,pjet)

      else
        write(6,*) 'dynamicscale ptj1'//
     &             ' not supported for this process.'
        stop
      endif
      
      return
      end
      
