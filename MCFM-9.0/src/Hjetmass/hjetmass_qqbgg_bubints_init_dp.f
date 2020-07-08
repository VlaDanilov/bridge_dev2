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
 
      subroutine hjetmass_qqbgg_bubints_init_dp
     &   (i1,i2,i3,i4,bubints,ord,s)
        use mod_qcdloop_c
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'sprods_decl.f'
      include 'scale.f'
      include 'masses.f'
      integer i1,i2,i3,i4
      include 'qqbggnames.f'
      integer ord
      complex(dp)::bubints(Nbubints)
      real*8 t
      real*8 mt2, mhsq,s1234

      t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

      mt2 = mt**2
      s1234=s(i1,i2)+s(i1,i3)+s(i1,i4)+s(i2,i3)+s(i2,i4)+s(i3,i4)
      bubints(:) = (0d0,0d0)

c   "(mhsq;mtsq,mtsq)",
        bubints(b1234) = qlI2(s1234, mt2, mt2, musq, ord)

c   "(s123;mtsq,mtsq)",
        bubints(b123) = qlI2(t(i1,i2,i3), mt2, mt2, musq, ord)

c   "(s124;mtsq,mtsq)",
        bubints(b124) = qlI2(t(i1,i2,i4), mt2, mt2, musq, ord)

c   "( s12;mtsq,mtsq)",
        bubints(b12) = qlI2(s(i1,i2), mt2, mt2, musq, ord)

c   "( s34;mtsq,mtsq)"]
        bubints(b34) = qlI2(s(i3,i4), mt2, mt2, musq, ord)

      end subroutine
