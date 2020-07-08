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
 

      complex*32 function hjetmass_box_pppm_0_s14_0_mhsq_s124_s134
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
      t2 = za(i1, i3)
      t3 = za(i2, i3)
      t4 = zb(i2, i1)
      t5 = zb(i3, i1)
      t6 = zb(i3, i2)
      t7 = za(i1, i4)
      t8 = zb(i4, i1)
      t9 = za(i2, i4)
      t10 = zb(i4, i2)
      t11 = za(i3, i4)
      t12 = zb(i4, i3)
      t13 = t10 * t11 + t2 * t4
      t14 = mt ** 2
      t15 = t1 * t5
      t16 = t9 * t12
      t17 = t16 + t15
      t18 = t6 ** 2
      t19 = t3 ** 2 * t18
      t20 = t19 * t17
      t21 = t2 * t5
      t22 = t11 * t12
      t23 = t22 + t21
      t24 = t7 * t3 * t6
      t25 = t24 * t8
      t26 = t9 * t10
      t27 = t1 * t4
      t28 = t3 * t6
      t23 = t28 * (t26 * t23 + t27 * t23 - t25)
      t29 = 16 * t14 * t3 * t6 * t13 * t20 + 4 * t23 ** 2
      t29 = sqrt(t29)
      t23 = 2 * t23
      t30 = -t29 + t23
      t20 = 0.1q1 / t20
      t31 = t16 + t15
      t31 = t26 * t31 + t27 * t31
      t32 = t22 + t21
      t33 = t26 * t32 + t27 * t32 - t25
      t34 = (0.1q1 / 0.2q1)
      t17 = t28 * t29 * t20 * t17 / 4
      t19 = t34 * t19 * t20 * (t21 * t31 + t22 * t31 + t25 * (-t16 - t15
     #))
      t21 = -t17 + t19 - t33
      t13 = t28 * t13
      t22 = 0.1q1 / t3
      t13 = 0.1q1 / t13
      t25 = 0.1q1 / t6
      t31 = t32 * t22
      t25 = t31 * t25
      t32 = t25 * t9
      t35 = t21 * t13
      t36 = t35 * t11
      t37 = t32 + t36
      t16 = t16 / 4
      t38 = -t16 * t30 * t20 + t10 * t37
      t23 = t29 + t23
      t17 = t17 + t19 - t33
      t19 = t17 * t13
      t11 = t19 * t11
      t29 = t32 + t11
      t16 = -t16 * t23 * t20 + t10 * t29
      t32 = t31 * t1
      t33 = t35 * t2
      t39 = t33 * t6 + t32
      t2 = t19 * t2
      t32 = t2 * t6 + t32
      t25 = t25 * t1
      t40 = t15 / 4
      t33 = -t40 * t30 * t20 + t4 * (t33 + t25)
      t2 = -t40 * t23 * t20 + t4 * (t2 + t25)
      t25 = t5 * t9
      t40 = t25 / 4
      t37 = t40 * t30 * t20 - t4 * t37
      t29 = t40 * t23 * t20 - t4 * t29
      t31 = t31 * t9
      t36 = t36 * t6 + t31
      t11 = t11 * t6 + t31
      t12 = 0.1q1 / t12
      t2 = 0.1q1 / t2
      t31 = 0.1q1 / t7
      t40 = 0.1q1 / t1
      t33 = 0.1q1 / t33
      t8 = 0.1q1 / t8
      t41 = t23 * t29
      t42 = t41 * t2
      t43 = t30 * t37
      t44 = t43 * t33
      t45 = t44 + t42
      t46 = t5 * t7
      t47 = t30 + t23
      t48 = t17 ** 2
      t49 = t21 ** 2
      t50 = t9 ** 2
      t51 = t5 ** 2
      t52 = t21 * t30
      t53 = t23 * t17
      t54 = t8 * t5
      t55 = t14 * t2
      t56 = t31 * t8
      t57 = t14 * t30
      t58 = t6 * t31
      t59 = t23 * t2
      t60 = t59 * t46
      t61 = (t46 - t36) * t33
      t62 = t28 * t13 ** 2
      t63 = t30 * t33
      t64 = t30 * t36
      t65 = t23 * t11
      t66 = t5 * t12
      t67 = t51 * t12
      t68 = t20 * t6
      t69 = t53 + t52
      t70 = t53 * t11
      t71 = t14 * t6
      t72 = t5 * t4
      t73 = t6 * t29
      t15 = t68 * (t12 * (t25 * (t54 * t69 * t10 - t71 * (t59 * t17 + t5
     #2 * t33)) * t40 + t72 * (t70 * (t56 - t2) + t52 * (t56 * t36 + t61
     #)) + t18 * t31 * t69 * t40 * t50) * t13 + t62 * t40 * t12 * (t44 *
     # t6 * t49 + t59 * (-t11 * t4 + t73) * t48) + t5 * (t45 * t9 + t56 
     #* (-t50 * t30 * t38 * t40 + t15 * (-t43 - t41) * t12)) * t22) + t6
     #8 * (t40 * (t9 * (t13 * (t6 * (t42 * t17 + t44 * t21) + t54 * (-t5
     #3 - t52) * t4) + t58 * t47 * t22 * t50 + t25 * (t23 * (-t56 * t16 
     #- t55) - t57 * t33) * t22) + t62 * t4 * (t61 * t49 * t30 + t60 * t
     #48) * t12) + t66 * (t60 * t19 * t4 + (t31 * (t65 + t64) + (-t63 - 
     #t59) * t5 * t14) * t22 * t9) + t67 * t22 * t45 * t1)
      t41 = t17 * t23 ** 2 + t21 * t30 ** 2
      t43 = 0.1q1 / t9
      t44 = t7 * t2
      t45 = t11 * t2
      t49 = t10 * t29
      t59 = t21 * t36
      t60 = t59 * t33
      t62 = t17 * t11
      t74 = t62 * t2
      t75 = t4 * t38
      t76 = t10 * t37
      t77 = t4 * t36
      t78 = t28 * t13 * t4
      t79 = t36 * t33
      t80 = t25 * t14
      t81 = t46 * t4
      t82 = t54 * t14
      t83 = t4 * t9
      t84 = t74 + t60
      t85 = t33 * t7 - t8
      t86 = t2 + t33
      t87 = t11 ** 2
      t88 = t4 ** 2
      t89 = t36 ** 2
      t90 = t6 * t40
      t91 = t4 * t22
      t92 = t40 * t8
      t93 = t92 * t91
      t94 = t4 * t12
      t95 = t86 * t40
      t96 = t95 * t22
      t97 = t6 * t4
      t98 = t11 * t22
      t99 = t10 * t11
      t100 = t22 * t12
      t101 = t100 * (-t4 * (t89 + t87) + t99 * t90 * t50 - t26 * t87 * t
     #40) * t31
      t102 = t26 + t38
      t103 = -t14 + t38
      t104 = t29 * t17
      t105 = t21 * t33
      t106 = t17 * t2
      t107 = t105 * t37
      t108 = t104 * t2
      t109 = t26 * t5
      t61 = t22 * (t9 * (-t45 * t72 * t10 + t58 * (t102 * t36 + t11 * t1
     #6) * t40) + t82 * (t109 * t40 + t77 * t31)) + t90 * (t10 * (-t105 
     #* t77 * t9 + t46 * (t108 + t107)) + t4 * (t105 * (-t3 * t89 * t43 
     #+ t46 * t103) - t106 * (t3 * t87 * t43 + t46 * t14)) + t6 * (t14 *
     # (-t108 - t107) + t3 * (t43 * ((-t54 + t61) * t37 * t21 + t104 * (
     #t2 * (t46 - t11) - t54)) + t81 * (t106 + t105)) - t70 * t9 * t31 *
     # t20)) * t13
      t24 = t12 * t61 + t5 * (t9 * (t96 * t6 * t12 * t14 ** 2 + t93 * t1
     #4 + t94 * (t44 * t90 * t19 - t79 * t22) * t10) + t94 * (-t97 * t19
     # * t8 + t98) + t35 * t88 * t6 * t12 * t85 - t93 * t50 * t10) + t12
     # * (t40 * (t6 * (t45 * t14 * t5 + t77 * t35 * (t14 - t38) * t33 + 
     #t19 * (t4 * (t16 * (t44 - t8) * t5 + t45 * (t14 - t16)) + t49 * (-
     #t45 - t54)) + t78 * (t54 * (-t21 - t17) - t60 - t74) - t54 * t35 *
     # (t76 + t75)) + t80 * (t10 * (t46 * (-t2 - t33) + t79) - t6) * t22
     #) + t5 * (t22 * (t77 + t73) + t45 * (t81 - t73) * t43) + t83 * t58
     # * t22 * (t82 + t11)) + t67 * (t1 * t88 * t7 * t22 * t86 + t92 * t
     #50 * t10 ** 2 * t22 + t26 * (t90 - t91) * t8) + t97 * (-t94 * t13 
     #* t84 + t96 * t50 * t14) + t66 * t24 * t4 * t13 * t40 * t84 * t43 
     #+ t101
      t60 = t6 * t9
      t61 = t29 * t2
      t70 = t66 * t6
      t72 = t76 + t75
      t74 = t5 * t36
      t75 = t5 * t10
      t55 = t100 * (-t75 * t71 * t56 * t50 * t40 + t5 * (t11 * (t4 * (-t
     #27 - t16) - t49) + t46 * (t16 * t4 + t49)) * t2 + (t10 * (-t31 * t
     #89 + t74 + t5 * (1 + t55) * t11) - t5 * (t38 + t16) * t6) * t40 * 
     #t9 - t58 * t29 * t11 + t5 * (-t36 * t72 + t46 * t72) * t33)
      t72 = t36 + t11
      t51 = t4 * (t5 * (-t89 * t33 * t43 - t45 * t6 - t79 * t6) + t51 * 
     #(-t8 * t43 * t72 + (t2 + t33) * t7 * t6)) + t58 * t36 * (-t22 * t3
     #7 + t90 * t9)
      t1 = t51 * t12 + t4 * (t40 * (-t80 * t7 * t22 * t86 + (-t54 * t3 *
     # t43 * (t62 + t59) + t26 * (t17 * (-t45 - t54) + t5 * t85 * t21)) 
     #* t13 * t12 * t6) + t100 * (t82 * t11 + t60 * t36) * t31 + t67 * t
     #7 * (t26 * t86 * t22 + t79 * t43)) + t24 + t5 * (t22 * (-t79 * t1 
     #* t12 - t8 * t9) + t44 * t19 * t6 * t12) * t88 + t70 * (t22 * (-t1
     #4 * t33 + 1) + (t46 - t36) * t33 * t43) * t37 + t70 * (t14 * (-t61
     # * t22 + t79 * t40) + t46 * (-t95 * t14 + t61 * t43)) + t66 * (-t5
     #4 * t6 * (t37 + t29) - t4 * t87 * t2) * t43 + t90 * (t58 * t9 * t1
     #1 * t12 + t67 * t14 * t8 - t25 * t4 * t8) + t55 - t90 * t35 * t36 
     #* t12 * (t68 * t9 * t30 * t31 + t76 * t33)
      t2 = t4 * (t39 + t32) - t99
      t3 = t9 * t40
      t7 = t39 * t37
      t24 = t32 * t29
      t33 = t10 * t40
      t43 = t9 * t22
      t44 = -t38 - t16
      t42 = t66 * (t56 * t26 * (-t90 * t36 + t91 * (-t36 - t11)) + t54 *
     # (t40 * (t10 * t72 + t44 * t6) + t68 * t69 * t13 * t4) - t42 * t19
     # * t20 * t18)
      t2 = t42 + t5 * (t56 * t12 * (t4 * (t11 * (t16 * t22 + t32 * t40) 
     #+ t36 * (t22 * t38 + t39 * t40) + t3 * t22 * (t103 * t39 + t32 * (
     #-t14 + t16))) + t90 * (t16 * (t43 * t14 + t11) + t36 * t38 + t24) 
     #+ t98 * t49 + t43 * t88 * t39) - t6 * (t63 * t68 * t35 * t12 + t92
     #) * t37 + t92 * (t4 * (-t43 * (t38 + t16) - t36 - t11) - t73)) + t
     #54 * (-t26 * t40 * t22 * (t37 + t29) + (t6 * (t40 * (t2 * t9 + t7)
     # + (t3 * t38 + t29 + t37) * t22 * t14) + (t32 * t88 + (-t11 * t14 
     #+ t24 + t7) * t40 * t10) * t22 * t9 + t33 * t2 * t22 * t50 - t89 *
     # t33 - t33 * t87 + t22 * t10 * (-t3 * (t26 + t14) + t37) * t36) * 
     #t31 * t12)
      t7 = t68 * t5
      t24 = t22 * t50
      t28 = t68 * t66 * t40 * t13 * (t23 * (t78 * t8 * t48 + t9 * (t99 *
     # t56 + t6) * t17) + t52 * (t8 * (t26 * t36 * t31 + t35 * t28 * t4)
     # + t60))
      t33 = t9 * (-t66 * t18 * t40 + t91 * (-t70 + t14 * t40 * (t45 + t7
     #9))) + t100 * t90 * (t58 * t14 - t75) * t50 + t91 * t66 * t14 * (-
     #t46 * t86 + t45 + t79) + t67 * t4 * (t22 * (-t27 + t14) - t6) * t8
      ret = t34 * (t7 * t56 * (t12 * (t23 * (-t5 * (t60 + t11) + t19 * t
     #6 * (-t83 + t29)) + t30 * (t9 * (-t109 * t22 - t97 * t35) - t74)) 
     #+ t24 * (-t14 * t47 + t26 * t47) * t40) + t7 * (t43 * (t3 * t47 - 
     #t66 * t47) + t56 * (t6 * t30 * t12 * (t35 * t37 - t25) + t43 * (t8
     #3 * t47 + t66 * (t23 * (-t26 + t14) + t57 + t27 * (-t30 - t23))) +
     # t3 * (t60 * t47 + t64 + t65 + (t52 * (t102 * t6 - t39 * t4) + t53
     # * (-t32 * t4 + t6 * (t26 + t16))) * t13 * t12)))) - t56 * t25 * t
     #13 * t20 ** 2 * t18 * (t41 * t40 * t9 - t66 * t41) / 8 + (0.3q1 / 
     #0.2q1) * t28 + 6 * t66 * t8 * ((t10 * (-t37 - t29) + t4 * t44) * t
     #22 * t5 + t71 * t13 * t4 * t40 * (t21 + t17)) + 8 * t33 - 3 * t24 
     #* t68 * t40 * t31 * (t65 + t64) - t15 + 4 * t1 + 2 * t2

      hjetmass_box_pppm_0_s14_0_mhsq_s124_s134 = ret/32q0/(0,1q0)
      return

      end function
