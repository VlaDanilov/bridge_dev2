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
 
      subroutine qqbggAxtri123x4x56(p1,p2,p3,p4,p5,p6,za,zb,
     & coeff0,coeff2)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer p1,p2,p3,p4,p5,p6
      real(dp):: s12,s123,s56
      complex(dp):: coeff0(2,2),coeff2(2,2),zab2

C--- begin statement function
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
C--- end statement functions
      s12=s(p1,p2)
      s123=s(p1,p2)+s(p1,p3)+s(p2,p3)
      s56=s(p5,p6)

c--- Note: coeff0 will be filled later using infrared relations
      coeff0(:,:) = czip

      include 'tri123x4x56coeffs.f'

      return
      end

