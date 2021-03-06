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
 

      double complex function hjetmass_triangle_pmpm_s34_0_0_dp 
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
      t2 = za(i1, i3)
      t3 = za(i2, i3)
      t4 = zb(i2, i1)
      t5 = zb(i4, i1)
      t6 = zb(i4, i2)
      t7 = t3 * t6
      t8 = t2 * t5
      t9 = t7 + t8
      t10 = za(i1, i4)
      t11 = zb(i3, i2)
      t12 = zb(i3, i1)
      t13 = t12 * t6
      t14 = t11 * t5
      t15 = -t13 + t14
      t16 = za(i2, i4)
      t17 = t16 * t6
      t18 = 3 * t10 * t5
      t19 = t2 ** 2
      t20 = t2 * t19
      t21 = t3 ** 2
      t22 = t5 ** 2
      t23 = t6 ** 2
      t24 = t17 * t5
      t25 = -3 * t3 * t15 + t24
      t26 = t10 * t22
      t27 = 0.1D1 / t9
      t28 = 0.1D1 / t4
      t29 = 0.1D1 / t1
      t30 = t27 ** 2
      t31 = 0.1D1 / t2 ** 2
      t32 = 0.1D1 / t6 ** 2
      ret = 16 * za(i3, i4) * (t10 * t21 ** 2 * t6 * t23 * t15 - t7 * t1
     #9 * (t26 * (t25 + t26) + t7 * (t12 * t3 + 3 * t16 * t5) * t15) + t
     #20 * t22 * (-t8 * t16 * t15 + t17 * t25 + t10 * t5 * (t3 * t15 + t
     #24)) + t1 * t3 * t4 * t5 * (t6 * (t19 * (t3 * (-2 * t13 + 3 * t14)
     # - 2 * t17 * t5) + t21 * t10 * t23) + t20 * t11 * t22 + t18 * t2 *
     # t3 * t23) - t2 * t3 * t21 * t23 * (t11 * t3 + t17 - t18) * t15 + 
     #t8 * t7 * t1 ** 2 * t4 ** 2 * t9) * zb(i4, i3) * t29 * t31 * t28 *
     # t32 * t27 * t30

      hjetmass_triangle_pmpm_s34_0_0_dp = ret/32d0/(0,1d0)
      return

      end function
