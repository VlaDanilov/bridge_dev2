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
 

      double complex function hjetmass_qqbgg_box_pmmp_0_0_s12_mhsq_s34_s124_dp
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double complex ret
          double precision mt

      t1 = za(i1, i2)
      t2 = za(i3, i4)
      t3 = zb(i2, i1)
      t4 = zb(i4, i3)
      t5 = za(i1, i4)
      t6 = zb(i4, i1)
      t7 = za(i2, i4)
      t8 = zb(i4, i2)
      t9 = za(i2, i3)
      t10 = t6 * za(i1, i3) + t8 * t9
      t11 = mt ** 2
      t12 = zb(i3, i1)
      t13 = t5 * t12
      t14 = t7 * zb(i3, i2)
      t15 = t14 + t13
      t16 = t2 * t4
      t17 = t16 * t15
      t18 = t1 * t3
      t5 = t5 * t6
      t8 = t7 * t8
      t19 = t16 * (t8 + t18 + t5)
      t20 = 16 * t11 * t10 * t17 + 4 * t19 ** 2
      t20 = cdsqrt(t20)
      t19 = 2 * t19
      t21 = -t20 - t19
      t19 = t20 - t19
      t17 = 0.1D1 / t17
      t22 = t8 + t5
      t13 = t16 * t17 * (t13 * t22 + t18 * (t14 + t13) + t14 * t22)
      t14 = t20 * t17 * t15 / 2
      t5 = 2 * t8 + 2 * t18 + 2 * t5
      t1 = 0.1D1 / t1
      t8 = 0.1D1 / t10
      t3 = 0.1D1 / t3
      t10 = t21 ** 2
      t15 = t19 ** 2
      t18 = t21 + t19
      ret = t16 * t3 * t17 ** 2 * t1 * ((t10 * (-t13 + t5 - t14) + t15 *
     # (-t13 + t5 + t14)) * t8 * t9 * t6 - t17 * t12 * t7 * (t21 * t10 +
     # t19 * t15)) / 4 - t17 * (-t2 * t6 ** 2 * t3 * t18 + (-t18 * t9 **
     # 2 + (t15 + t10) * t3 * t17 * t7 * t6 * t2) * t1 * t4) + 2 * t11 *
     # t17 * t9 * t6 * t1 * t3 * (t21 + t19)

      hjetmass_qqbgg_box_pmmp_0_0_s12_mhsq_s34_s124_dp
     & = ret/32d0/(0,1d0)
      return

      end function
