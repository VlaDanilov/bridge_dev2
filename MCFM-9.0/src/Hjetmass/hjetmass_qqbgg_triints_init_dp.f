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
 
      subroutine hjetmass_qqbgg_triints_init_dp
     &    (i1,i2,i3,i4,triints,ord,s)
        use mod_qcdloop_c
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'constants.f'
      include 'sprods_decl.f'
      include 'scale.f'
      include 'masses.f'
      integer i1,i2,i3,i4
      integer ord
      include 'qqbggnames.f'
      complex*16 triints(Ntriints)
      real*8 t
      real*8 mt2, mhsq

        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        mt2 = mt**2
        mhsq = s(i1,i2)+s(i1,i3)+s(i1,i4)+s(i2,i3)+s(i2,i4)+s(i3,i4)

        triints(:) = (0d0,0d0)

c   "(s123,mhsq,   0;mtsq,mtsq,mtsq)",
        triints(c4x123) = 
     &  qlI3(t(i1,i2,i3), mhsq, 0d0, mt2, mt2, mt2, musq, ord)

c   "(s124,mhsq,   0;mtsq,mtsq,mtsq)",
        triints(c3x124) = 
     &  qlI3(t(i1,i2,i4), mhsq, 0d0, mt2, mt2, mt2, musq, ord)

c   "( s12,s123,   0;mtsq,mtsq,mtsq)",
        triints(c3x12) = 
     &  qlI3(s(i1,i2), t(i1,i2,i3), 0d0, mt2, mt2, mt2, musq, ord)

c   "( s12,s124,   0;mtsq,mtsq,mtsq)",
        triints(c4x12) = 
     &  qlI3(s(i1,i2), t(i1,i2,i4), 0d0, mt2, mt2, mt2, musq, ord)

c   "( s34,   0,   0;mtsq,mtsq,mtsq)",
        triints(c3x4) = 
     &  qlI3(s(i3,i4), 0d0, 0d0, mt2, mt2, mt2, musq, ord)

c   "( s34,mhsq, s12;mtsq,mtsq,mtsq)"
        triints(c12x34) = 
     &  qlI3(s(i3,i4), mhsq, s(i1,i2), mt2, mt2, mt2, musq, ord)

      end subroutine
