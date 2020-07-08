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
 
      subroutine phase2(r,p1,p2,p3,p4,wt,*)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'mxdim.f'
      include 'x1x2.f'
c---- generate phase space for 2-->2 process
c---- r(mxdim),p1(4),p2(4) are inputs reversed in sign 
c---- from physical values 
c---- phase space for -p1-p2 --> p3+p4
c---- with all 2 pi's (ie 1/(2*pi)^2)
      real(dp):: r(mxdim)
      real(dp):: p1(4),p2(4),p3(4),p4(4)
      real(dp):: cosphi,sinphi,u,phi,rtshat,costh,sinth
      real(dp):: wt,wt0
      include 'energy.f'
      parameter(wt0=1._dp/8._dp/pi)
      rtshat=sqrt(xx(1)*xx(2))*sqrts
C write out vectors in +,-,T,T notation
      u=r(3)
      costh=2._dp*u-1._dp
      sinth=sqrt(1._dp-costh**2)
      phi=two*pi*r(4)
      sinphi=sin(phi)
      cosphi=cos(phi)

      p3(4)=+half*sqrts*(u*xx(1)+(1._dp-u)*xx(2))
      p3(1)=+half*sinth*sinphi*rtshat
      p3(2)=+half*sinth*cosphi*rtshat
      p3(3)=+half*sqrts*(u*xx(1)-(1._dp-u)*xx(2))


      p4(4)=+half*sqrts*((1._dp-u)*xx(1)+u*xx(2))
      p4(1)=-half*sinth*sinphi*rtshat
      p4(2)=-half*sinth*cosphi*rtshat
      p4(3)=+half*sqrts*((1._dp-u)*xx(1)-u*xx(2))

c---debug
      write(6,*) 'p3',p3(4),p3(3),p3(2),p3(1)
      write(6,*) 'p4',p4(4),p4(3),p4(2),p4(1)

      p3(4)=+u*p1(4)+(one-u)*p2(4)
      p3(1)=+half*sinth*sinphi*rtshat
      p3(2)=+half*sinth*cosphi*rtshat
      p3(4)=+u*p1(3)+(one-u)*p2(3)


      p4(4)=+(one-u)*p1(4)+u*p2(4)
      p4(1)=-half*sinth*sinphi*rtshat
      p4(2)=-half*sinth*cosphi*rtshat
      p4(3)=+(one-u)*p1(3)+u*p2(3)

      write(6,*) 'p3',p3(4),p3(3),p3(2),p3(1)
      write(6,*) 'p4',p4(4),p4(3),p4(2),p4(1)
c      pause
      wt=wt0

      end

