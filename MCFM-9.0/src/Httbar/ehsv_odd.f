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
 
c--- These are the functions similar to equations (A.8)-(A.14)
c--- but which describe the case of a pseudoscalar (CP odd) Higgs
c--- These results are adapted from Spira et al., hep-ph/9504378.

      function ehsva4_odd(s,t,u)
      implicit none
      include 'types.f'
      complex(dp):: ehsva4_odd
C     ehsv:EqnA.8
      
      real(dp):: s,t,u
      complex(dp):: ehsvb4_odd
      ehsva4_odd=ehsvb4_odd(s,t,u)+ehsvb4_odd(u,s,t)+ehsvb4_odd(t,u,s)
      return 
      end

      function ehsva2_odd(s,t,u)
      implicit none
      include 'types.f'
      complex(dp):: ehsva2_odd
C     ehsv:EqnA.9
      
      real(dp):: s,t,u
      complex(dp):: ehsvb2_odd
      ehsva2_odd=ehsvb2_odd(s,t,u)+ehsvb2_odd(s,u,t)
      return 
      end

      function ehsvb4_odd(s,t,u)
      implicit none
      include 'types.f'
      complex(dp):: ehsvb4_odd
      
C     ehsv:EqnA.10
      include 'masses.f'
      real(dp):: hmass2,s,t,u
      complex(dp):: w2,w3
      hmass2=s+t+u
      ehsvb4_odd=mbsq/hmass2*(w2(hmass2)-w2(s)-w3(s,t,u,hmass2))/6._dp
      return 
      end

      function ehsvb2_odd(s,t,u)
      implicit none
      include 'types.f'
      complex(dp):: ehsvb2_odd
C     ehsv:EqnA.11
      
      include 'masses.f'
      real(dp):: hmass2,s,t,u
      complex(dp):: w2,w3
      hmass2=s+t+u
      ehsvb2_odd=mbsq/hmass2**2*(
     &                       -s*w3(s,t,u,hmass2)+s/2._dp*w3(t,s,u,hmass2)
     &                       +2._dp*s*(0.75_dp-u/(s+u))*w2(hmass2)
     &                       -s/2._dp*w2(s)
     &                       -s*(1._dp-2._dp*u/(s+u))*w2(t))/6._dp
      return 
      end

      function ehsva5_odd(s,t,u)
      implicit none
      include 'types.f'
      complex(dp):: ehsva5_odd
C     ehsv:EqnA.14
      
      include 'masses.f'
      real(dp):: hmass2,s,t,u
      complex(dp):: w2
      hmass2=s+t+u
      ehsva5_odd=mbsq/hmass2*2._dp/3._dp*(w2(s)-w2(hmass2))
      return 
      end

