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
 

      complex*32 function hjetmass_box_pppm_0_s23_0_mhsq_s123_s234
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          
          complex*32 ret
          double precision mt

      t1 = za(i2, i4)
      t2 = zb(i2, i1)
      t3 = zb(i3, i1)
      t4 = zb(i4, i2)
      t5 = za(i3, i4)
      t6 = zb(i4, i3)
      t7 = t1 * t4
      t8 = t5 * t6
      t9 = t8 + t7
      t10 = za(i1, i4)
      t11 = zb(i4, i1)
      t12 = za(i1, i3)
      t13 = za(i2, i3)
      t14 = zb(i3, i2)
      t15 = za(i1, i2)
      t16 = t15 * t4
      t17 = t12 * t6
      t18 = t17 + t16
      t19 = t11 ** 2
      t20 = t10 ** 2 * t19
      t21 = t20 * t18
      t22 = t8 + t7
      t23 = t10 * t13 * t14 * t11
      t24 = t15 * t2
      t25 = t10 * t11
      t22 = t25 * (t22 * t12 * t3 + t24 * t22 - t23)
      t26 = t5 * t3
      t27 = t1 * t2 + t26
      t28 = mt ** 2
      t29 = t25 * t28
      t30 = 16 * t29 * t27 * t21 + 4 * t22 ** 2
      t30 = sqrt(t30)
      t21 = 0.1q1 / t21
      t16 = t17 + t16
      t7 = -t7 * t16 - t8 * t16
      t8 = t9 * t12
      t31 = -t24 * t9 - t8 * t3 + t23
      t18 = t25 * t30 * t21 * t18 / 4
      t7 = t20 * t21 * (t12 * t7 * t3 + t23 * t16 + t24 * t7) / 2
      t16 = t18 - t31 + t7
      t23 = t25 * t27
      t22 = 2 * t22
      t27 = t22 - t30
      t23 = 0.1q1 / t23
      t32 = 0.1q1 / t10
      t33 = 0.1q1 / t11
      t8 = t8 * t32 * t33
      t34 = t23 * t16
      t17 = t17 / 4
      t35 = t3 * (-t34 * t5 + t8) - t17 * t27 * t21
      t7 = -t18 - t31 + t7
      t18 = t22 + t30
      t22 = t7 * t23
      t8 = t3 * (-t22 * t5 + t8) - t17 * t18 * t21
      t9 = t9 * t33
      t17 = t9 * t3
      t30 = t27 * t21 / 4
      t31 = t30 * t10
      t33 = -t31 * t6 + t17
      t36 = t9 * t2
      t37 = t18 * t21 / 4
      t38 = t37 * t10
      t39 = -t38 * t4 + t36
      t17 = -t38 * t6 + t17
      t9 = t9 * t32 * t15
      t32 = -t34 * t1 + t9
      t30 = t30 * t15
      t38 = t3 * t32 - t30 * t6
      t9 = -t22 * t1 + t9
      t37 = t37 * t15
      t40 = t9 * t3 - t37 * t6
      t31 = -t31 * t4 + t36
      t9 = t9 * t2 - t37 * t4
      t4 = t2 * t32 - t30 * t4
      t30 = 0.1q1 / t15
      t13 = 0.1q1 / t13
      t6 = 0.1q1 / t6
      t32 = t30 * (t28 + t8)
      t36 = t28 * t30
      t37 = t2 * t40
      t41 = t2 * t38
      t42 = t1 * t14
      t43 = t6 * t13
      t44 = 0.1q1 / t1
      t45 = t17 * t44
      t46 = t45 - t14
      t47 = t38 ** 2
      t48 = t40 ** 2
      t49 = t30 * t10
      t50 = t38 * t27
      t51 = t8 * t3
      t52 = t5 * t30
      t53 = t13 * t5
      t54 = 0.1q1 / t4
      t55 = 0.1q1 / t9
      t56 = t38 * t54
      t57 = t40 * t55
      t58 = t57 + t56
      t45 = -t45 + t14
      t59 = -t33 * t44 + t14
      t60 = t38 + t40
      t61 = t1 ** 2
      t62 = t2 ** 2
      t63 = t3 * t33
      t64 = t14 * t30
      t65 = t13 * t3
      t66 = t28 * t3
      t67 = t40 * t18
      t68 = t2 * t6
      t69 = t28 * t1
      t70 = t41 * t54 * t6
      t45 = t13 * (t3 * (t69 * (-t68 + t52) * t13 + t6 * (t2 * (t56 * t5
     # - t1) - t1 * (t54 + t55) * t30 * t28 ** 2) * t14) + t70 * (t36 - 
     #t2) * t33 - t57 * t62 * t17 * t6 - t43 * t5 ** 2 * t44 * t60 * t3 
     #** 2) + t13 * (t37 * t36 * t17 * t55 * t6 + t6 * (t2 * (t57 * t3 *
     # t45 + t56 * (t59 * t30 * t38 - t63 * t44) + t30 * t45 * t55 * t48
     #) + t3 * (t14 * (-t30 * t60 - t3) + t30 * (t17 * t40 + t33 * t38) 
     #* t44) + t66 * (-t64 * t58 - t65)) * t5 + t42 * t62 * t6 * t58 + t
     #36 * t61 * t2 * t13 + t25 * t6 * t30 * (t67 * t45 + t50 * t59) * t
     #21)
      t58 = -t35 + t4
      t59 = t27 + t18
      t60 = t14 * t59
      t62 = t27 * t31
      t71 = t18 * t39
      t72 = t17 * t18
      t73 = t27 * t33
      t74 = t27 * t4
      t75 = t21 * t13 * t11
      t76 = t72 + t73
      t77 = t12 * t13
      t78 = t12 * t17
      t79 = t43 * t42
      t24 = t75 * (t13 * (t1 * (-t71 - t62) + t6 * (t18 * (t10 * (-t37 +
     # t51) - t78 * t3) + t10 * (t3 * t35 - t41) * t27)) + t30 * (t1 * (
     #t60 * t10 + t77 * t76) + t53 * t10 * (t50 + t67) + t10 * (t3 * t7 
     #* t18 * (t1 * t9 * t13 + t17) + (t41 * t1 * t13 + t63) * t27 * t16
     #) * t6 * t23 - t77 * t60 * t61) + t79 * t59 * t28) + t75 * (t6 * (
     #-t72 * t14 - t73 * (t65 * t12 + t14) + t42 * (t65 * t59 * t12 + t6
     #0 + t24 * (-t27 - t18) * t13) + t65 * (t71 + t62) * t15) + t49 * (
     #(t18 * t9 + t74 + (t3 * t16 * t27 * t58 + (t37 - t51) * t18 * t7) 
     #* t6 * t23) * t13 * t1 - t73 - t72) + t53 * t42 * t59)
      t59 = t14 ** 2
      t61 = t67 * t7
      t62 = t50 * t16
      t63 = t27 * t54
      t65 = t18 * t55
      t71 = t26 * t13
      t77 = t43 * t21
      t80 = t40 * t30
      t81 = t67 * t55
      t5 = t75 * (-t72 * t53 + t73 * (-t70 * t10 * t44 + t53 * (t3 * t15
     # * t44 * t6 - 1)) + t6 * (t1 * (t13 * (t18 * (t80 * t12 - t9) - t7
     #4) + t49 * (-t63 * t16 * (t41 + t66) + t65 * (-t37 - t66) * t7) * 
     #t23) + t10 * (t27 * (-t30 * t54 * t47 - t56 * t3) - t81 * (t80 + t
     #3)) * t44 * t5) * t14 + (-t72 * t37 * t55 + t71 * (-t50 - t67)) * 
     #t6 * t44 * t10 + t65 * t69 * t59 * t6) + t77 * t11 * (t13 * (t50 *
     # (t30 * (t10 * t58 + t12 * (t42 - t33)) + t31) + t67 * (t30 * (t10
     # * (-t8 + t9) - t78) + t39)) + t10 * (t44 * t76 + (-t65 - t63) * t
     #14 * t28) * t3 + t71 * (-t14 * t27 + t18 * t46) * t15 + t63 * t69 
     #* t59 + t49 * (t61 * t17 * t55 + t62 * t33 * t54) * t23 * t2)
      t12 = t21 ** 2
      t31 = t18 ** 2
      t27 = t27 ** 2
      t25 = t77 * t25 * t3 * (t52 * (t62 + t61) * t13 * t23 - t60)
      t39 = t13 ** 2
      t6 = t12 * t6 * t39 * t14 * t19 * t10 * (-t15 * (t27 + t31) + (t16
     # * t27 + t31 * t7) * t23 * t1)
      t7 = t68 * t26 * t39 * (t38 + t40)
      ret = 2 * t43 * (-t36 * t3 * (t33 + t17) + (-t51 * t1 - t52 * (t48
     # + t47)) * t13 * t2 + t53 * (t18 * (t49 * t44 * t48 + t40 * t46) +
     # t50 * (t44 * (t10 * t38 * t30 + t33) - t14)) * t21 * t11) + 2 * t
     #43 * (t13 * (t1 * (t3 * (t35 * (t30 * (t28 - t4) - t2) - t9 * t32 
     #+ t36 * (t8 - t4)) + t37 * (t32 + t2) + t41 * (t30 * (t28 + t35) +
     # t2)) + t26 * (t3 * (-t35 - t8) + t30 * (t38 * (-t28 + t4) + t40 *
     # (-t28 + t9)))) + t42 * t30 * (t2 * (-t38 - t40) + t3 * (t9 + t4))
     #) + t24 / 2 + 8 * t79 * t36 * t2 * (t57 + t56) + (0.3q1 / 0.2q1) *
     # t25 - t6 / 8 - 3 * t77 * t64 * t29 * (t50 * t54 + t81) + 6 * t7 -
     # 4 * t45 + t5 + t43 * t20 * t12 * t14 * (t56 * (-t34 * t30 + t44) 
     #* t27 + t57 * t31 * (-t22 * t30 + t44)) / 4

      hjetmass_box_pppm_0_s23_0_mhsq_s123_s234 = ret/32q0/(0,1q0)
      return

      end function
