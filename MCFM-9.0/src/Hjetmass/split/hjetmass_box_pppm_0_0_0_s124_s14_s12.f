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
 

      complex*32 function hjetmass_box_pppm_0_0_0_s124_s14_s12
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
      t7 = mt ** 2
      t8 = t1 * t2
      t9 = t3 * t4
      t10 = 0.1q1 / t4
      t11 = 0.1q1 / t3
      t12 = (0.1q1 / 0.2q1)
      t13 = t12 *sqrt(-t9 * t8 * (-2 * t8 * t3 * t4 - 8 * t7 * t5 * t6
     #)) * sqrt(0.2q1) * t11 * t10
      t14 = -t13 + t8
      t15 = zb(i3, i2)
      t16 = za(i3, i4)
      t9 = t9 * t8
      t7 = t9 * (4 * t7 * t5 * t6 + t9)
      t7 = sqrt(t7)
      t17 = zb(i3, i1)
      t18 = za(i1, i3)
      t19 = zb(i4, i3)
      t20 = 0.1q1 / t5
      t1 = 0.1q1 / t1
      t6 = 0.1q1 / t6
      t21 = 0.1q1 / t2
      t22 = t12 * ((t7 + t9) * t11 * t20 * t21 * t10 * t16 * t17 + t14 *
     # t1 * t6 * t18 * t19)
      t23 = t16 * t19
      t24 = t22 + t23
      t8 = t13 + t8
      t7 = t12 * (-(t7 - t9) * t11 * t20 * t21 * t10 * t16 * t17 + t8 * 
     #t1 * t6 * t18 * t19)
      t9 = t7 + t23
      t10 = t15 * za(i2, i3) + t17 * t18
      t11 = -t10 + t22
      t7 = -t10 + t7
      t10 = 0.1q1 / t24
      t11 = 0.1q1 / t11
      t9 = 0.1q1 / t9
      t7 = 0.1q1 / t7
      t13 = 0.1q1 / t16
      t16 = t7 - t9
      t18 = -t10 + t11
      t21 = t14 ** 2
      t22 = t14 * t21
      t23 = t8 ** 2
      t24 = t8 * t23
      t25 = t16 * t23 + t18 * t21
      t26 = t11 * t21 + t23 * t7
      t27 = t10 - t11
      t28 = t9 - t7
      t29 = t6 ** 2
      t30 = t1 ** 2
      t31 = t1 * t30
      t32 = t3 ** 2
      t33 = t3 * t32
      t34 = t21 * t10
      t35 = t23 * t9
      t36 = t4 * t15
      t37 = t20 * t3
      t38 = t13 * t30
      t39 = t22 * t27 + t24 * t28
      t40 = t22 * t11
      t41 = t24 * t7
      t42 = t41 + t40
      t9 = t13 * (t10 * t22 + t24 * t9)
      t10 = t9 * t2
      t43 = t20 * t2
      t40 = t41 + t40
      t41 = t14 * t18 + t16 * t8
      t7 = t11 * t14 + t7 * t8
      t8 = t43 * t31 * t6 * t29 * t33 * t4 * t19 * t13 * t40
      ret = -t12 * t6 * t31 * t32 * (-t37 * t13 * t29 * t42 * t15 * t4 *
     #* 2 + (t10 + (t16 * t24 + t18 * t22) * t20 * t17) * t6 * t19 - t9 
     #* t17 + t36 * t6 * (t39 * t6 * t20 * t19 + t9) + t43 * t39 * t29 *
     # t19 ** 2) + (0.3q1 / 0.2q1) * t31 * t29 * t33 * t17 * t4 * t20 * 
     #t13 * t40 + (0.5q1 / 0.2q1) * t8 + t6 * t30 * t3 * (t17 * (t15 * t
     #25 + (-t35 - t34) * t13 * t3 * t2) + t2 * ((t26 * t19 * t2 + t36 *
     # (t35 + t34)) * t13 * t3 + t19 * t15 * (t21 * t27 + t23 * t28)) * 
     #t6 + t37 * t25 * t17 ** 2 + t38 * t20 * t19 * t4 * t32 * (t21 ** 2
     # * t27 + t23 ** 2 * t28) * t29) - 2 * t6 * t1 * t3 * (t2 * (t17 * 
     #(t37 * t41 * t17 + t15 * t41) + (t5 * t15 * t7 + t7 * t17 * t3) * 
     #t13 * t2) + (t10 * t32 * t19 * t20 * t29 - t3 * t15 * t13 * t42 * 
     #t6) * t30 * t4) - 4 * t38 * t29 * t4 * t32 * t2 * (t37 * t26 * t17
     # + t15 * t26)

      hjetmass_box_pppm_0_0_0_s124_s14_s12 = ret/32q0*(0,1q0)
      return

      end function
