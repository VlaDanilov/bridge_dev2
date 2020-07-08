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
 
      subroutine hjetmass_bubints_init(i1,i2,i3,i4,bubints,ord,s)
          use mod_qcdloop_c
          implicit none
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'scale.f'
          include 'masses.f'
          integer i1,i2,i3,i4
          integer ord
          complex*32 bubints(9)
          real*16 mhsq_q, qmt2
          real*16 t, s(mxpart,mxpart)
          real*16 mt2, mhsq, qmusq

        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        mt2 = mt**2
        mhsq = s(i1,i2)+s(i1,i3)+s(i1,i4)+s(i2,i3)+s(i2,i4)+s(i3,i4)
        qmusq = musq

        bubints(:) = (0q0,0q0)

!           (mhsq;mtsq,mtsq)
        bubints(1) = qlI2q(mhsq, mt2, mt2, qmusq, ord)

!           (s123;mtsq,mtsq)
        bubints(2) = qlI2q(t(i1,i2,i3), mt2, mt2, qmusq, ord)

!           (s124;mtsq,mtsq)
        bubints(3) = qlI2q(t(i1,i2,i4), mt2, mt2, qmusq, ord)

!           ( s12;mtsq,mtsq)
        bubints(4) = qlI2q(s(i1,i2), mt2, mt2, qmusq, ord)

!           (s134;mtsq,mtsq)
        bubints(5) = qlI2q(t(i1,i3,i4), mt2, mt2, qmusq, ord)

!           ( s14;mtsq,mtsq)
        bubints(6) = qlI2q(s(i1,i4), mt2, mt2, qmusq, ord)

!           (s234;mtsq,mtsq)
        bubints(7) = qlI2q(t(i2,i3,i4), mt2, mt2, qmusq, ord)

!           ( s23;mtsq,mtsq)
        bubints(8) = qlI2q(s(i2,i3), mt2, mt2, qmusq, ord)

!           ( s34;mtsq,mtsq)
        bubints(9) = qlI2q(s(i3,i4), mt2, mt2, qmusq, ord)

      end subroutine
