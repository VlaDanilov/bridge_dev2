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
 
      subroutine BLC_mm(k1,k2,k3,k4,k5,k6,coeff,ampBLC)
      implicit none
      include 'types.f'
c--- This is the leading colour amplitude in the notation
c---  (q1, Qb2, Q3, qb4)

c--- This routine makes use of the identity B(1,2,3,4)=-A(1,3,2,4) 
      
      include 'Wbbmlabels.f'
      integer:: k1,k2,k3,k4,k5,k6      
      complex(dp):: ampBLC(2,2),coefftmp(0:4,20)
      
      call ALC_mm(k1,k3,k2,k4,k5,k6,coefftmp,ampBLC)

c--- minus sign applied in the identities below
      coeff(4,d2x3x4)=-coefftmp(4,d2x3x4)
      coeff(4,d1x2x3)=-coefftmp(4,d1x2x3)
      coeff(4,d1x23x4)=-coefftmp(4,d1x23x4)
      coeff(4,d1x2x34)=-coefftmp(4,d1x2x34)
      coeff(4,d12x3x4)=-coefftmp(4,d12x3x4)
      
      coeff(3,c23x4)=-coefftmp(3,c23x4)
      coeff(3,c1x23)=-coefftmp(3,c1x23)
      coeff(3,c2x3)=-coefftmp(3,c2x3)
      coeff(3,c2x34)=-coefftmp(3,c2x34)
      coeff(3,c3x4)=-coefftmp(3,c3x4)
      coeff(3,c12x3)=-coefftmp(3,c12x3)
      coeff(3,c1x2)=-coefftmp(3,c1x2)
      coeff(3,c1x234)=-coefftmp(3,c1x234)
      coeff(3,c123x4)=-coefftmp(3,c123x4)
      coeff(3,c12x34)=-coefftmp(3,c12x34)

      coeff(2,b123)=-coefftmp(2,b123)
      coeff(2,b234)=-coefftmp(2,b234)
      coeff(2,b23)=-coefftmp(2,b23)
      coeff(2,b1234)=-coefftmp(2,b1234)
      coeff(2,b12)=-coefftmp(2,b12)
      coeff(2,b34)=-coefftmp(2,b34)
      coeff(2,b2x1m)=-coefftmp(2,b2x1m)
      
      coeff(1,a0m)=-coefftmp(1,a0m)
      
      coeff(0,irat)=-coefftmp(0,irat)
      
      
      return
      end
      
