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
 

      complex*32 function hjetmass_qqbgg_triangle_pmmp_s124_mhsq_0_rat
     &     (i1,i2,i3,i4,za,zb)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          real*16 cg

      t1 = za(i2, i3)
      t2 = zb(i2, i1)
      t3 = za(i3, i4)
      t4 = zb(i4, i1)
      t5 = t1 * t2 - t3 * t4
      t6 = za(i1, i3)
      t7 = zb(i3, i1)
      t8 = zb(i3, i2)
      t9 = zb(i4, i3)
      t10 = t6 * t7
      t11 = t1 * t8
      t12 = t3 * t9
      t13 = t10 + t11 + t12
      if ( real(t13) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t11 = cg * sqrt(t13 ** 2) + t10 + t11 + t12
      t12 = za(i1, i2)
      t13 = za(i2, i4)
      t14 = za(i1, i4)
      t15 = zb(i4, i2)
      t16 = t12 * t2
      t17 = t14 * t4
      t18 = t13 * t15 + t16 + t17
      t11 = 0.1q1 / t11
      t19 = -0.2q1 * t1 * t7 * t18 * t11 + t13 * t4
      t20 = 0.2q1 * t3 * t7 * t18 * t11 - t13 * t2
      t10 = -0.2q1 * t10 * t18 * t11 + t16 + t17
      t16 = t6 * (t10 * t7 + t19 * t8)
      t17 = t1 * t15 + t4 * t6
      t21 = 0.2q1 * t6
      t8 = t7 * (t1 * (-t21 * t8 * t18 * t11 + t14 * t15) + t10 * t6)
      t21 = t21 * t18 * t9 * t11 - t12 * t15
      t22 = 0.1q1 / t7
      t18 = 0.1q1 / t18
      t23 = 0.1q1 / t12
      t24 = 0.1q1 / t2
      t14 = 0.1q1 / t14
      t25 = t4 * t23 * t24
      t8 = 0.1q1 / t8
      t21 = 0.1q1 / t21
      t26 = 0.1q1 / t3
      t13 = 0.1q1 / t13
      t16 = 0.1q1 / t16
      t27 = 0.1q1 / t9
      t20 = 0.1q1 / t20
      t28 = t6 * t16
      t21 = t26 * t22 * t21
      t29 = (t20 * t27 + t28) * t19
      t30 = t1 * t10
      t31 = t5 * t1
      t32 = t24 * t23 * t13 * t11
      t33 = t32 * t1
      t34 = t9 * t22
      ret = -0.64q2 * t1 * t4 * t22 * t18 * (t25 + t14) - 0.16q2 * t33 *
     # (t31 * ((-t21 - t8) * t10 * t9 + t26) + (-t30 * t8 - t29) * t17 *
     # t7) - 0.128q3 * t31 * t11 * t18 * (t14 * (t12 * t13 + t34) + (t34
     # * t23 + t13) * t24 * t4) + 0.8q1 * t25 * t13 * t1 * (t30 * (t21 +
     # t8) + t29) + 0.32q2 * t32 * t19 ** 2 * t7 * (t3 * t17 * t20 ** 2 
     #* t27 + t6 ** 2 * t9 * t16 ** 2 * (t1 * (t15 * t3 + t2 * t6) - t5 
     #* t6)) + 0.48q2 * t33 * t5 * t19 * (t28 * t9 + t20)

      hjetmass_qqbgg_triangle_pmmp_s124_mhsq_0_rat = ret/32q0*(0,1q0)
      return

      end function
