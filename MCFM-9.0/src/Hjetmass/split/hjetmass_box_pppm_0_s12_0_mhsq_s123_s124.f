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
 

      complex*32 function hjetmass_box_pppm_0_s12_0_mhsq_s123_s124
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          
          complex*32 ret
          double precision mt

      t1 = za(i1, i3)
      t2 = za(i1, i4)
      t3 = za(i3, i4)
      t4 = zb(i3, i1)
      t5 = zb(i4, i1)
      t6 = zb(i4, i3)
      t7 = za(i2, i3)
      t8 = zb(i3, i2)
      t9 = za(i2, i4)
      t10 = zb(i4, i2)
      t11 = za(i1, i2)
      t12 = zb(i2, i1)
      t13 = t2 * t5
      t14 = t9 * t10
      t15 = t14 + t13
      t16 = t11 * t3 * t12 * t6
      t17 = t7 * t8
      t18 = t1 * t4
      t19 = t3 * t6
      t15 = t19 * (t17 * t15 + t18 * t15 - t16)
      t20 = t9 * t8
      t21 = t2 * t4
      t22 = t21 + t20
      t23 = mt ** 2
      t24 = t19 * t23
      t25 = t1 * t5
      t26 = t7 * t10
      t27 = t25 + t26
      t28 = t3 ** 2
      t29 = t28 * t6 ** 2
      t30 = t29 * t27
      t31 = 16 * t24 * t22 * t30 + 4 * t15 ** 2
      t31 = sqrt(t31)
      t15 = 2 * t15
      t32 = -t31 + t15
      t33 = t14 + t13
      t34 = 0.1q1 / t6
      t30 = 0.1q1 / t30
      t35 = (0.1q1 / 0.4q1)
      t36 = t33 * t34
      t37 = t36 * t8
      t38 = t35 * t32 * t30
      t39 = t38 * t3
      t40 = t39 * t10 - t37
      t15 = t31 + t15
      t41 = t35 * t15 * t30
      t42 = t41 * t3
      t10 = t42 * t10 - t37
      t37 = t18 + t17
      t37 = t25 * t37 + t26 * t37
      t33 = t17 * t33 + t18 * t33 - t16
      t43 = (0.1q1 / 0.2q1)
      t27 = t35 * t19 * t31 * t30 * t27
      t13 = t43 * t29 * t30 * (t13 * t37 + t14 * t37 + t16 * (-t26 - t25
     #))
      t14 = t13 - t27 - t33
      t16 = t19 * t22
      t25 = 0.1q1 / t3
      t16 = 0.1q1 / t16
      t25 = t36 * t25 * t7
      t31 = t14 * t16
      t37 = t31 * t9 + t25
      t44 = -t38 * t7 * t5 + t4 * t37
      t13 = t13 + t27 - t33
      t27 = t36 * t4
      t33 = t42 * t5 - t27
      t36 = t13 * t16
      t25 = t36 * t9 + t25
      t42 = t8 * t25 - t26 * t41
      t25 = -t41 * t7 * t5 + t4 * t25
      t26 = -t26 * t38 + t8 * t37
      t5 = t39 * t5 - t27
      t11 = 0.1q1 / t11
      t27 = 0.1q1 / t26
      t37 = 0.1q1 / t7
      t38 = 0.1q1 / t42
      t39 = 0.1q1 / t9
      t41 = t15 * t25
      t45 = t32 * t44
      t46 = t45 * t27 + t41 * t38
      t47 = t32 * t5
      t48 = t15 * t33
      t49 = t48 + t47
      t50 = t23 * t11
      t51 = t32 + t15
      t52 = t11 ** 2
      t53 = t4 ** 2
      t54 = (t50 - t12) * t9
      t55 = t9 * t11
      t56 = t32 * t14
      t57 = t15 * t13
      t58 = t32 * t27
      t59 = t15 * t38
      t60 = t2 * t3 * t39
      t61 = t47 * t27
      t62 = t2 * t52
      t63 = t9 * t12
      t64 = t63 + t33
      t65 = t63 + t5
      t66 = t32 ** 2
      t67 = t8 ** 2
      t68 = t15 ** 2
      t69 = t3 * t4
      t70 = t33 * t39
      t71 = t3 * t44
      t72 = t69 * t39
      t73 = t39 * (t41 + t45)
      t74 = t68 * t33
      t75 = t19 * t11
      t76 = t75 * t39
      t18 = t18 * t37 * t11
      t77 = t8 * t25
      t78 = t39 * t5
      t79 = t4 * t7
      t80 = -t72 * t61 * t8 + (t58 * t8 * (t79 * t12 + t78 * t44) + t59 
     #* (t77 * t70 + (t17 - t23) * t12 * t4) + t69 * ((t57 + t56) * t16 
     #* t4 - t73) * t37) * t11 * t2
      t46 = t30 * t80 + t30 * (t8 * (t15 * (t38 * (-t70 * t69 + t50 * (t
     #69 - t33)) - t3 * t25 * t11 * t37) + t11 * t46 * t12 * t2 + t11 * 
     #(-t71 * t37 + (t69 - t5) * t27 * t23) * t32) + t21 * (t73 * t62 + 
     #t58 * (t72 - t50) * t12) + t76 * (t5 * t66 + t74) * t30 + t18 * (t
     #63 * (-t32 - t15) - t47 - t48) + t11 * t7 * (t58 * t65 + t59 * t64
     #) * t67) + t30 * (t4 * (t37 * (t12 * (t60 * t46 + t55 * (-t59 - t5
     #8) * t1 * t23) + t3 * (t57 * (t38 * (t54 - t33) + t55) + t56 * (t2
     #7 * (t54 - t5) + t55)) * t16 * t8) - t62 * t23 * t51 + t17 * t2 * 
     #t11 * t39 * (t48 * t38 + t61)) + t62 * t49 * t6 + t20 * t52 * (t15
     # * (-t23 + t42) + t32 * (-t23 + t26)) + t59 * t60 * t12 * t53)
      t54 = t30 ** 2
      t56 = t69 * t12 * t37
      t57 = t54 * t6
      t60 = t68 + t66
      t61 = t66 * t27
      t63 = t68 * t38
      t47 = -t48 - t47
      t48 = t9 ** 2
      t73 = t32 * t40
      t80 = t15 * t10
      t51 = t12 * t51
      t81 = t51 * t37
      t82 = t4 * t26
      t83 = t30 * t11
      t84 = t2 * t25
      t85 = t8 * t44
      t6 = t83 * (t37 * (t55 * t1 * (t32 * (-t85 + t82) - t41 * t8) + t1
     #9 * (t15 * (t11 * (t42 * t9 + t84) + t33) + t32 * (t55 * t26 + t5)
     #)) - t84 * t76 * t68 * t30 + t51 * t55 * t2 * t6) + t83 * (t2 * (t
     #11 * (t32 * (t19 * t44 * t37 + t82) + t4 * t15 * t42) - t51 * t4) 
     #+ t4 * (t80 + t73) + t47 * t8 + t9 * (t6 * (t11 * (t47 * t37 * t1 
     #+ t73 + t80) + t81 * t3) - t51 * t8 + t18 * t15 * t42) + t83 * t6 
     #* t2 * (t66 * (t5 * t7 - t71) + t74 * t7) * t39 - t81 * t1 * t6 * 
     #t11 * t48)
      t18 = -t32 - t15
      t19 = t4 * (-t42 - t26) + t8 * (t44 + t25)
      t47 = t23 * t37
      t19 = t37 * t34 * t11 * (t20 * t19 + t21 * t19) + t69 * (t11 * (-t
     #21 * t18 * t39 - t18 * t8) + t47 * (t11 * t18 + t12 * (t59 + t58))
     #) * t30
      t51 = -t8 + t47
      t69 = t23 * t4
      t71 = t8 * t33
      t73 = t34 * t8
      t74 = t34 * t4
      t76 = t55 * t23 * t8
      t18 = t74 * (t38 * (t2 * (t4 * t8 * (t50 - t70) + (t69 - t77) * t3
     #7 * t12) + t71 * t47 + t76 * (-t47 + t8)) + t8 * (t37 * (-t12 * t2
     # * t44 - t55 * t23 ** 2) - t78 * t21) * t27) + t37 * (t2 * (t69 * 
     #t9 * t52 + t34 * t39 * (t33 + t5) * t53) + t23 * t48 * t8 * t52) +
     # t73 * t4 * (-t12 * t22 + t84 * (t50 - t70) * t37 - t71) * t38 + t
     #74 * (t2 * (t4 * (t12 * t51 + t50 * t8) + t85 * t37 * (t50 - t78))
     # + t8 * t5 * t51 + t67 * (t50 - t12) * t9) * t27 + t56 * t18 * t30
      t20 = t74 * (t55 * t67 + t21 * (t37 * (-t50 + t12) + t21 * t11 * t
     #39) + t47 * t20 * t12 * (t27 + t38))
      t21 = t7 * t39
      t15 = t29 * t30 * t54 * t11 * t12 * ((t31 + t21) * t27 * t32 * t66
     # + t15 * t68 * t38 * (t36 + t21))
      t21 = t30 * (-t72 * t37 * t49 + t76 * (t59 + t58) * t12)
      t22 = t24 * t54 * t12 * t11 * (t63 + t61)
      ret = -t35 * (t57 * (t52 * t12 * (t60 * t9 * t1 - t2 * t7 * t60) +
     # (t39 * (t11 * t60 + t12 * (t63 + t61)) + t11 * (t13 * t68 + t14 *
     # t66) * t16 * t37) * t4 * t28) + t57 * (t66 * (t11 * (t11 * (t1 * 
     #t5 + t3 * (-t1 * t44 * t37 + t26) - t40 * t7) + t3 * (t17 * t5 + (
     #-t79 - t44) * t12 * t2) * t39 * t27) + t31 * t27 * t3 * (t65 * t11
     # * t8 + t56)) + t68 * (t11 * (t11 * (t1 * t33 - t7 * t10 + t3 * (-
     #t1 * t25 * t37 + t42)) + t3 * (t17 * t33 + (-t79 - t25) * t12 * t2
     #) * t39 * t38) + t36 * t38 * t3 * (t64 * t11 * t8 + t56)))) + t43 
     #* t6 - 8 * t20 + (0.3q1 / 0.2q1) * t62 * t30 * t8 * (t41 + t45) + 
     #(0.5q1 / 0.4q1) * t75 * t54 * t12 * (t68 + t66) - t15 / 16 - 3 * t
     #21 + (0.3q1 / 0.4q1) * t22 + t46 + 2 * t19 - 4 * t18 - 16 * t73 * 
     #t2 * t53 * t11

      hjetmass_box_pppm_0_s12_0_mhsq_s123_s124 = ret/32q0/(0,1q0)
      return

      end function
