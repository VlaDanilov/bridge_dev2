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
 
      function ehsva4_dm(s,t,u,s34)
      implicit none
      include 'types.f'
      complex(dp):: ehsva4_dm
C     ehsv:EqnA.8
      
      real(dp):: s,t,u,s34
      complex(dp):: ehsvb4_dm
      ehsva4_dm=ehsvb4_dm(s,t,u,s34)+ehsvb4_dm(u,s,t,s34)
     &     +ehsvb4_dm(t,u,s,s34)
      return 
      end

      function ehsva2_dm(s,t,u,s34)
      implicit none
      include 'types.f'
      complex(dp):: ehsva2_dm
C     ehsv:EqnA.9
      
      real(dp):: s,t,u,s34
      complex(dp):: ehsvb2_dm
      ehsva2_dm=ehsvb2_dm(s,t,u,s34)+ehsvb2_dm(s,u,t,s34)
      return 
      end

      function ehsvb4_dm(s,t,u,s34)
      implicit none
      include 'types.f'
      complex(dp):: ehsvb4_dm
      
C     ehsv:EqnA.10
      include 'masses.f'
      real(dp):: hmass2,s,t,u,s34
      complex(dp):: w2,w3

      hmass2=s34
c--- The Fermilab preprint has w2(s), but it makes no difference due
c---  to symmetrization in ehsva4 above 
      ehsvb4_dm=mbsq/hmass2*(-2d0/3d0
     & +(mbsq/hmass2-0.25d0)*(w2(t)-w2(hmass2)+w3(s,t,u,hmass2)))
      return 
      end

      function ehsvb2_dm(s,t,u,s34)
      implicit none
      include 'types.f'
      complex(dp):: ehsvb2_dm
C     ehsv:EqnA.11
      
      include 'masses.f'
      real(dp):: hmass2,s,t,u,s34
      complex(dp):: w1,w2,w3
      hmass2=s34
      ehsvb2_dm=mbsq/hmass2**2*(s*(u-s)/(s+u)
     & +2d0*u*t*(u+2d0*s)/(s+u)**2*(w1(t)-w1(hmass2))
     & +(mbsq-0.25d0*s)
     & *(0.5d0*w2(s)+0.5d0*w2(hmass2)-w2(t)+w3(s,t,u,hmass2))
     & +s**2*(2d0*mbsq/(s+u)**2-0.5d0/(s+u))*(w2(t)-w2(hmass2))
     & +0.5d0*u*t/s*(w2(hmass2)-2d0*w2(t))
     & +0.125d0*(s-12d0*mbsq-4d0*u*t/s)*w3(t,s,u,hmass2))
      return 
      end

      function ehsva5_dm(s,t,u,s34)
      implicit none
      include 'types.f'
      complex(dp):: ehsva5_dm
C     ehsv:EqnA.14
      
      include 'masses.f'
      real(dp):: hmass2,s,t,u,s34
      complex(dp):: w1,w2
      hmass2=s34
      ehsva5_dm=mbsq/hmass2*(4d0+4d0*s/(u+t)*(w1(s)-w1(hmass2))
     & +(1d0-4d0*mbsq/(u+t))*(w2(s)-w2(hmass2)))
      return 
      end




