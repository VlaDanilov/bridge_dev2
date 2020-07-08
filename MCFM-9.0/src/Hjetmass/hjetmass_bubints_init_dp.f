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
 
      subroutine hjetmass_bubints_init_dp(i1,i2,i3,i4,bubints,ord,s)
        use mod_qcdloop_c
          implicit none
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'scale.f'
          include 'masses.f'
          integer i1,i2,i3,i4
          integer ord
          complex*16 bubints(9)
          real*8 mhsq_q, qmt2
          real*8 t, s(mxpart,mxpart)
          real*8 mt2, mhsq

        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        mt2 = mt**2
        mhsq = s(i1,i2)+s(i1,i3)+s(i1,i4)+s(i2,i3)+s(i2,i4)+s(i3,i4)

        bubints(:) = (0d0,0d0)

!           (mhsq;mtsq,mtsq)
        bubints(1) = qlI2(mhsq, mt2, mt2, musq, ord)

!           (s123;mtsq,mtsq)
        bubints(2) = qlI2(t(i1,i2,i3), mt2, mt2, musq, ord)

!           (s124;mtsq,mtsq)
        bubints(3) = qlI2(t(i1,i2,i4), mt2, mt2, musq, ord)

!           ( s12;mtsq,mtsq)
        bubints(4) = qlI2(s(i1,i2), mt2, mt2, musq, ord)

!           (s134;mtsq,mtsq)
        bubints(5) = qlI2(t(i1,i3,i4), mt2, mt2, musq, ord)

!           ( s14;mtsq,mtsq)
        bubints(6) = qlI2(s(i1,i4), mt2, mt2, musq, ord)

!           (s234;mtsq,mtsq)
        bubints(7) = qlI2(t(i2,i3,i4), mt2, mt2, musq, ord)

!           ( s23;mtsq,mtsq)
        bubints(8) = qlI2(s(i2,i3), mt2, mt2, musq, ord)

!           ( s34;mtsq,mtsq)
        bubints(9) = qlI2(s(i3,i4), mt2, mt2, musq, ord)

      end subroutine
