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
 
      subroutine hjetmass_boxints_init(i1,i2,i3,i4,boxints,ord,s)
        use mod_qcdloop_c
          implicit none
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'scale.f'
          include 'masses.f'

          integer i1,i2,i3,i4
          integer ord
          complex*32 boxints(16)
          real*16 t, s(mxpart,mxpart)
          real*16 mt2, mhsq, qmusq

        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        mt2 = mt**2
        mhsq = s(i1,i2)+s(i1,i3)+s(i1,i4)+s(i2,i3)+s(i2,i4)+s(i3,i4)
        qmusq = musq

        boxints(:) = (0q0,0q0)

!      (   0,   0,   0,s123, s12, s23;mts,mts,mts,mts)
        boxints(1) = qlI4q(0q0, 0q0, 0q0, t(i1,i2,i3), s(i1,i2), s(i2,i3),
     &          mt2, mt2, mt2, mt2, qmusq, ord)

!      (   0,   0,   0,s124, s14, s12;mts,mts,mts,mts)
        boxints(2) = qlI4q(0q0, 0q0, 0q0, t(i1,i2,i4), s(i1,i4), s(i1,i2),
     &          mt2, mt2, mt2, mt2, qmusq, ord)

!      (   0,   0,   0,s134, s34, s14;mts,mts,mts,mts)
        boxints(3) = qlI4q(0q0, 0q0, 0q0, t(i1,i3,i4), s(i3,i4), s(i1,i4),
     &          mt2, mt2, mt2, mt2, qmusq, ord)

!      (   0,   0,   0,s234, s23, s34;mts,mts,mts,mts)
        boxints(4) = qlI4q(0q0, 0q0, 0q0, t(i2,i3,i4), s(i2,i3), s(i3,i4),
     &          mt2, mt2, mt2, mt2, qmusq, ord)

!      (   0,   0, s12,mhsq, s34,s123;mts,mts,mts,mts)
        boxints(5) = qlI4q(0q0, 0q0, s(i1,i2), mhsq, s(i3,i4), t(i1,i2,i3),
     &          mt2, mt2, mt2, mt2, qmusq, ord)

!      (   0,   0, s12,mhsq, s34,s124;mts,mts,mts,mts)
        boxints(6) = qlI4q(0q0, 0q0, s(i1,i2), mhsq, s(i3,i4), t(i1,i2,i4),
     &          mt2, mt2, mt2, mt2, qmusq, ord)

!      (   0,   0, s14,mhsq, s23,s124;mts,mts,mts,mts)
        boxints(7) = qlI4q(0q0, 0q0, s(i1,i4), mhsq, s(i2,i3), t(i1,i2,i4),
     &          mt2, mt2, mt2, mt2, qmusq, ord)

!      (   0,   0, s14,mhsq, s23,s134;mts,mts,mts,mts)
        boxints(8) = qlI4q(0q0, 0q0, s(i1,i4), mhsq, s(i2,i3), t(i1,i3,i4),
     &          mt2, mt2, mt2, mt2, qmusq, ord)

!      (   0,   0, s23,mhsq, s14,s123;mts,mts,mts,mts)
        boxints(9) = qlI4q(0q0, 0q0, s(i2,i3), mhsq, s(i1,i4), t(i1,i2,i3),
     &          mt2, mt2, mt2, mt2, qmusq, ord)

!      (   0,   0, s23,mhsq, s14,s234;mts,mts,mts,mts)
        boxints(10) = qlI4q(0q0, 0q0, s(i2,i3), mhsq, s(i1,i4), t(i2,i3,i4),
     &          mt2, mt2, mt2, mt2, qmusq, ord)

!      (   0,   0, s34,mhsq, s12,s134;mts,mts,mts,mts)
        boxints(11) = qlI4q(0q0, 0q0, s(i3,i4), mhsq, s(i1,i2), t(i1,i3,i4),
     &          mt2, mt2, mt2, mt2, qmusq, ord)

!      (   0,   0, s34,mhsq, s12,s234;mts,mts,mts,mts)
        boxints(12) = qlI4q(0q0, 0q0, s(i3,i4), mhsq, s(i1,i2), t(i2,i3,i4),
     &          mt2, mt2, mt2, mt2, qmusq, ord)

!      (   0, s12,   0,mhsq,s123,s124;mts,mts,mts,mts)
        boxints(13) = qlI4q(0q0, s(i1,i2), 0q0, mhsq, t(i1,i2,i3), t(i1,i2,i4),
     &          mt2, mt2, mt2, mt2, qmusq, ord)

!      (   0, s14,   0,mhsq,s124,s134;mts,mts,mts,mts)
        boxints(14) = qlI4q(0q0, s(i1,i4), 0q0, mhsq, t(i1,i2,i4), t(i1,i3,i4),
     &          mt2, mt2, mt2, mt2, qmusq, ord)

!      (   0, s23,   0,mhsq,s123,s234;mts,mts,mts,mts)
        boxints(15) = qlI4q(0q0, s(i2,i3), 0q0, mhsq, t(i1,i2,i3), t(i2,i3,i4),
     &          mt2, mt2, mt2, mt2, qmusq, ord)

!      (   0, s34,   0,mhsq,s134,s234;mts,mts,mts,mts)
        boxints(16) = qlI4q(0q0, s(i3,i4), 0q0, mhsq, t(i1,i3,i4), t(i2,i3,i4),
     &          mt2, mt2, mt2, mt2, qmusq, ord)

      end subroutine
