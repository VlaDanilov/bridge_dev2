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
 

      complex*32 function hjetmass_triangle_ppmm_s23_mhsq_s14_rat
     &     (i1,i2,i3,i4,za,zb)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          
          complex*32 ret

          parameter (cg = 1q0)

      t1 = za(i1, i2)
      t2 = zb(i2, i1)
      t3 = za(i1, i3)
      t4 = zb(i3, i1)
      t5 = za(i2, i4)
      t6 = zb(i4, i2)
      t7 = za(i3, i4)
      t8 = zb(i4, i3)
      t9 = za(i1, i4)
      t10 = za(i2, i3)
      t11 = zb(i3, i2)
      t12 = zb(i4, i1)
      t13 = t1 * t2
      t14 = t3 * t4
      t15 = t5 * t6
      t16 = t7 * t8
      t17 = t14 + t15 + t16 + t13
      t17 = -4 * t9 * t10 * t11 * t12 + t17 ** 2
      t15 = cg * sqrt(t17) + t13 + t14 + t15 + t16
      t16 = t2 * t5 + t4 * t7
      t17 = t15 / 2
      t18 = t9 * t10
      t19 = t18 * t2
      t20 = t17 * t7 + t19
      t21 = -t11 * t20
      t22 = t17 * t5 - t18 * t4
      t23 = t12 * t22
      t24 = t7 * t11 * t12
      t25 = t17 * t2 + t24
      t26 = -t10 * t25
      t27 = (0.1q1 / 0.4q1)
      t28 = t15 ** 2
      t29 = t18 * t11 * t12
      t30 = t27 * t28 - t29
      t13 = t14 + t13
      t14 = -t15 * t13 / 2 + t29
      t4 = -t5 * t11 * t12 + t17 * t4
      t5 = t10 * t4
      t29 = t9 * t12
      t13 = -t29 * t13 + t29 * t17
      t31 = -t3 * t11 * t12 + t17 * t6
      t32 = t9 * t31
      t33 = t2 * t3 + t6 * t7
      t34 = t1 * t6 + t3 * t8
      t35 = t9 * (t1 * t11 * t12 + t17 * t8)
      t3 = t17 * t3 - t18 * t6
      t6 = t3 * t11
      t1 = -t11 * (t17 * t1 + t18 * t8)
      t8 = t12 * t20
      t17 = t9 * t25
      t18 = t10 * t11 * t33
      t20 = t29 * t16
      t25 = t15 * t27
      t24 = t25 * t2 + t24 / 2
      t27 = t15 * t9 * t24
      t29 = 0.1q1 / t11
      t36 = 0.1q1 / t10
      t37 = 0.1q1 / t12
      t38 = 0.1q1 / t15
      t34 = 0.1q1 / t34
      t1 = 0.1q1 / t1
      t39 = 0.1q1 / t9
      t30 = 0.1q1 / t30 ** 2
      t40 = t14 ** 2
      t2 = t2 * t7
      t35 = 0.1q1 / t35
      t41 = 0.1q1 / t14
      t42 = 0.1q1 / t5
      t43 = 0.1q1 / t23
      t44 = 0.1q1 / t13
      t45 = t26 ** 2
      t46 = t8 ** 2
      t47 = t13 ** 2
      t48 = t21 * t8
      t49 = t39 * t30 * t37
      t50 = t49 * t36 * t29
      t51 = -0.1q1 / t5
      t52 = t26 * t27 * t30 + t2
      t53 = t6 * t34 * t30
      t4 = t39 * t37 * t36 * t29 * (t40 * (-t53 * t18 * t9 * t4 * t1 ** 
     #2 + t30 * t34 * (t17 * t18 + t33 * t15 * t11 * (t25 * t7 + t19 / 2
     #)) * t1) + t52 * t44 * t43 * t20 * t8 + (-t32 * t27 * t10 * t31 * 
     #t34 ** 2 * t30 * t39 ** 2 * t37 ** 2 - t49 * t18 * t27 * t34) * t3
     #5 * t47 - t46 * t20 * t5 * t27 * t44 * t43 ** 2 * t30 + t32 * t52 
     #* t39 * t34 * t37 * t35 * t13 - t53 * (t18 * t20 - t48) * t1 * t14
     # + t2 * t16 * t15 * t10 * t24 * t41 * t51)
      t5 = t20 ** 2
      t7 = t21 * t20 * t15
      ret = -32 * t6 * t38 * t39 * t36 * t29 * t37 * t34 * t1 * (t21 * t
     #40 * t12 * t3 * t30 * t34 - t2 * t14) - 2 * t50 * t15 * t28 * t16 
     #* t32 * (t13 * t33 * t39 * t37 * t34 * t35 + t16 * t45 * t41 ** 2 
     #* t42) + 16 * t4 + 8 * t50 * (-t7 * t45 * t41 * t42 + (-t15 * t6 *
     # t5 * t44 ** 2 + t7 * t44) * t43 * t46 - t15 * t8 * t5 * t18 * t44
     # * t43) + 4 * t50 * t28 * (t16 * (t48 * t26 * t41 * t42 + (-t21 * 
     #t23 * t42 ** 2 + t17 * t42) * t41 * t45 - t46 * t17 * t44 * t43) +
     # (-t11 * t22 * t32 * t37 * t34 * t39 * t35 ** 2 + t21 * t37 * t34 
     #* t39 * t35) * t33 * t47) - t50 * t28 ** 2 * t16 ** 2 * t26 * t33 
     #* t41 * t42

      hjetmass_triangle_ppmm_s23_mhsq_s14_rat = ret/32q0/(0,1q0)
      return

      end function
