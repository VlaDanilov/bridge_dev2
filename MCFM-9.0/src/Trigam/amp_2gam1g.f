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
 
!====== C.Williams March 2013 
!====== Amplitude for q(i1)^-qb(i2)^+gluon(i3)^-gamma(i4)^+gamma(i5)^+ 
      
      function amp_2gam1g(p1,p2,p3,p4,p5,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: amp_2gam1g
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      include 'zprods_decl.f' 
      integer:: p1,p2,p3,p4,p5 

      amp_2gam1g=za(p2,p1)*za(p1,p3)**2/
     &     (za(p1,p5)*za(p1,p4)*za(p2,p5)*za(p2,p4))
      
      return 
      end
