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
 
c--- Combinations of bubble integrals, 0704.3914v3 Eq. (3.21) 
c--- The hatted versions of these functions are further defined
c--- by Eqs. (3.25), (3.26) and (3.27)
      function BGRL1(s,t)
      implicit none
      include 'types.f'
      include 'cplx.h'
      complex(dp):: BGRL1
      
      complex(dp):: lnrat
      real(dp):: s,t
      BGRL1=lnrat(-s,-t)/cplx1(s-t)
      return
      end

      function BGRL2(s,t)
      implicit none
      include 'types.f'
      include 'cplx.h'
      complex(dp):: BGRL2
      
      complex(dp):: lnrat
      real(dp):: s,t
      BGRL2=lnrat(-s,-t)/cplx1(s-t)**2
      return
      end

      function BGRL2hat(s,t)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'cplx.h'
      complex(dp):: BGRL2hat
      
      complex(dp):: lnrat
      real(dp):: s,t
      BGRL2hat=lnrat(-s,-t)/cplx1(s-t)**2
     & -cplx1(half*(s+t)/((s-t)*s*t))
      return
      end

      function BGRL3(s,t)
      implicit none
      include 'types.f'
      include 'cplx.h'
      complex(dp):: BGRL3
      complex(dp):: lnrat
      real(dp):: s,t
      BGRL3=lnrat(-s,-t)/cplx1(s-t)**3
      return
      end

      function BGRL3hat(s,t)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'cplx.h'
      complex(dp):: BGRL3hat
      complex(dp):: lnrat
      real(dp):: s,t
      BGRL3hat=lnrat(-s,-t)/cplx1(s-t)**3
     & -cplx1(half*(s+t)/((s-t)**2*s*t))
      return
      end
