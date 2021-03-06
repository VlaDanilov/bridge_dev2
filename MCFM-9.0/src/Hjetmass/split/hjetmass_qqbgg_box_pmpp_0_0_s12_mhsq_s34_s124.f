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
 

      complex*32 function hjetmass_qqbgg_box_pmpp_0_0_s12_mhsq_s34_s124
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          double precision mt

      t1 = za(i1, i2)
      t2 = zb(i2, i1)
      t3 = za(i1, i4)
      t4 = zb(i4, i1)
      t5 = za(i2, i4)
      t6 = zb(i4, i2)
      t7 = za(i3, i4)
      t8 = zb(i3, i1)
      t9 = zb(i4, i3)
      t10 = t3 * t8
      t11 = t5 * zb(i3, i2)
      t12 = t11 + t10
      t13 = t7 * t9
      t14 = t13 * t12
      t15 = za(i1, i3)
      t16 = za(i2, i3)
      t17 = t15 * t4
      t18 = t16 * t6 + t17
      t19 = mt ** 2
      t20 = t1 * t2
      t21 = t3 * t4
      t6 = t5 * t6
      t22 = t13 * (t6 + t20 + t21)
      t23 = 16 * t19 * t18 * t14 + 4 * t22 ** 2
      t23 = sqrt(t23)
      t14 = 0.1q1 / t14
      t24 = t6 + t21
      t11 = t13 * t14 * (t10 * t24 + t11 * t24 + t20 * (t11 + t10))
      t13 = 2
      t12 = t23 * t14 * t12 / 2
      t6 = t13 * (t6 + t20 + t21)
      t20 = -t12 + t6 - t11
      t22 = t13 * t22
      t24 = -t22 - t23
      t22 = -t22 + t23
      t6 = t12 + t6 - t11
      t11 = 0.1q1 / t18
      t10 = -t10 / 4
      t12 = t17 / 2
      t18 = t10 * t24 * t14 + t12 * t20 * t11
      t23 = -t15 * t5 + t16 * t3
      t25 = t20 * t24
      t26 = t25 * t14 * t11 * t9 * t23
      t27 = t22 * t6
      t23 = t27 * t14 * t11 * t9 * t23
      t10 = t10 * t22 * t14 + t12 * t6 * t11
      t12 = 0.1q1 / t1
      t28 = 0.1q1 / t15
      t2 = 0.1q1 / t2
      t3 = 0.1q1 / t3
      t29 = t27 + t25
      t30 = t8 * t2
      t31 = t2 * t28
      t32 = t20 + t6
      t17 = t17 * t3
      t33 = t3 * t4
      t34 = t2 * t11 * t9
      t35 = t2 * t9
      t36 = t31 * t14 ** 2 * t9 ** 2 * t11 * t7 * t8 * (t20 * t24 ** 2 +
     # t6 * t22 ** 2)
      t37 = t5 * t9 * t12
      t1 = t9 * ((-t1 * t28 + t30) * t3 * t19 - t4 * t2 * (t37 + t8))
      t38 = t3 * (t18 + t10)
      t37 = t35 * (t28 * (t3 * (t10 ** 2 + t10 * t19 + t18 ** 2 + t19 * 
     #t18) + t4 * (-t18 - t10)) + t38 * t8 + t38 * t37)
      ret = t13 * t34 * ((t16 * (t3 * (-t10 * t6 - t18 * t20) + t32 * t4
     #) - t17 * t32 * t5) * t12 * t9 + t33 * (t20 * (t19 - t18) + t6 * (
     #t19 - t10) + (-t20 - t6) * t15 * t8)) + t36 / 8 + 8 * t1 - t33 * t
     #34 * t12 * (t20 * t26 + t23 * t6) / 4 - t9 * ((t28 * (t16 * t29 + 
     #(-t27 - t25) * t2 * t7 * t4) + (t30 * t29 * t12 * t15 - t25 - t27)
     # * t3 * t5) * t14 * t11 * t9 + t31 * t12 * (t3 * (-t10 * t23 - t18
     # * t26) + t4 * (t26 + t23))) / 2 + t35 * (t9 * (t14 * (t25 * (t5 *
     # (-t18 * t3 + t4) + (t21 - t18) * t28 * t16) + t27 * (t5 * (-t10 *
     # t3 + t4) + (t21 - t10) * t28 * t16)) * t12 * t11 + t17 * t16 * (t
     #20 ** 2 + t6 ** 2) * t12 * t11 ** 2) + t28 * t14 * t8 * (t22 * (t1
     #9 + t10) + t24 * (t19 + t18) + t21 * (-t24 - t22))) + 4 * t37

      hjetmass_qqbgg_box_pmpp_0_0_s12_mhsq_s34_s124 = ret/32q0/(0,1q0)
      return

      end function
