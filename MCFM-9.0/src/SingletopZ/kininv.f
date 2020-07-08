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
 
      subroutine kininv(p)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'decl_kininv.f'
      real(dp):: p(mxpart,4)
      
      p12(:)=p(1,:)+p(2,:)
      p26(:)=p(2,:)+p(6,:)
      p25(:)=p(2,:)+p(5,:)
      p34(:)=p(3,:)+p(4,:)
      p346(:)=p(3,:)+p(4,:)+p(6,:)
      p134(:)=p(1,:)+p(3,:)+p(4,:)
      p1346(:)=p(1,:)+p(3,:)+p(4,:)+p(6,:)
      p16(:)=p(1,:)+p(6,:)
      p234(:)=p(2,:)+p(3,:)+p(4,:)
      p345(:)=p(3,:)+p(4,:)+p(5,:)
      p235(:)=p(2,:)+p(3,:)+p(5,:)

      s34=p34(4)**2-p34(1)**2-p34(2)**2-p34(3)**2
      s25=p25(4)**2-p25(1)**2-p25(2)**2-p25(3)**2
      s12=p12(4)**2-p12(1)**2-p12(2)**2-p12(3)**2
      s26=p26(4)**2-p26(1)**2-p26(2)**2-p26(3)**2
      s346=p346(4)**2-p346(1)**2-p346(2)**2-p346(3)**2
      s134=p134(4)**2-p134(1)**2-p134(2)**2-p134(3)**2
      s16=p16(4)**2-p16(1)**2-p16(2)**2-p16(3)**2
      s234=p234(4)**2-p234(1)**2-p234(2)**2-p234(3)**2
      s235=p235(4)**2-p235(1)**2-p235(2)**2-p235(3)**2
      s345=p345(4)**2-p345(1)**2-p345(2)**2-p345(3)**2

      return
      end

      
