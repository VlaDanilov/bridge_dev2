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
 
      subroutine scaleset_Msqptj1sq(p,mu0)
      implicit none
      include 'types.f'
c--- subroutine to calculate dynamic scale equal to
c---  sqrt(M^2+ptj1^2), where M is the mass of the particle (34)
c--- and ptj1 is the pt of the leading jet in the event
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'kprocess.f'
      include 'jetlabel.f'
      include 'npart.f'
      include 'breit.f'
      integer:: isub,oldjets
      real(dp):: p(mxpart,4),pjet(mxpart,4),mu0,pt,rcut,ptj1,ptj2
      common/rcut/rcut

      if((kcase==kW_1jet) .or.
     &   (kcase==kZ_1jet) .or.
     &   (kcase==kggfus1) .or.
     &   (kcase==khjetma)) then

c--- first work out whether this point is real radiatio or not
        if (abs(p(npart+2,4)) > 1.e-8_dp) then
          isub=0  ! real term
        else
          isub=1  ! subtraction term
        endif
      
c-- cluster jets but make sure recorded number of jets is not changed
        oldjets=jets     
        call genclust2(p,rcut,pjet,isub)
      
c        write(6,*) 'partons:'
c        if (p(5,4) >= 1.e-8_dp) write(6,*) 'pt5',pt(5,p)
c        if (p(6,4) >= 1.e-8_dp) write(6,*) 'pt6',pt(6,p)  
c      write(6,*) 'jets:'
c        if (pjet(5,4) >= 1.e-8_dp) write(6,*) 'pt5',pt(5,pjet)
c        if (pjet(6,4) >= 1.e-8_dp) write(6,*) 'pt6',pt(6,pjet)
c        write(6,*)
      
c--- restore old value of jets
        jets=oldjets

c--- order according to pt      
      ptj1=pt(5,pjet)
      ptj2=pt(6,pjet)
      if (ptj2 > ptj1) ptj1=ptj2      

c--- assign scale
        mu0=mass3**2+ptj1**2
        mu0=sqrt(abs(mu0))
      
      else
        write(6,*) 'dynamicscale sqrt(M^2+ptj1^2)'//
     &             ' not supported for this process.'
        stop
      endif
      
      return
      end
      
