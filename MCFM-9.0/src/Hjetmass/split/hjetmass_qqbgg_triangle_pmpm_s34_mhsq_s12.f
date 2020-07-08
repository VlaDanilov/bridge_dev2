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
 

      complex*32 function hjetmass_qqbgg_triangle_pmpm_s34_mhsq_s12
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
      t3 = za(i1, i3)
      t4 = zb(i3, i1)
      t5 = za(i2, i3)
      t6 = zb(i3, i2)
      t7 = za(i1, i4)
      t8 = zb(i4, i1)
      t9 = za(i2, i4)
      t10 = zb(i4, i2)
      t11 = t3 * t4
      t12 = t5 * t6
      t13 = t7 * t8
      t14 = t9 * t10
      t15 = t11 + t12 + t13 + t14
      t16 = za(i3, i4)
      t17 = zb(i4, i3)
      t15 = -4 * t1 * t16 * t2 * t17 + t15 ** 2
      t15 = t11 - sqrt(t15) + t12 + t13 + t14
      t18 = 2 * t1 * t2
      t19 = t15 + t18
      t20 = (0.1q1 / 0.2q1)
      t18 = t20 * t15 ** 2 - t18 * t16 * t17
      t21 = 2 * t16 * t17 + t15
      t22 = (0.1q1 / 0.4q1)
      t23 = t15 ** 2
      t24 = t1 * t16
      t25 = t24 * t2 * t17
      t26 = t23 * t22 - t25
      t26 = 0.1q1 / t26
      t27 = t15 * t22
      t24 = t24 * t20
      t28 = t24 * t4 - t27 * t9
      t29 = t15 * t2 * t26
      t30 = t29 * t28
      t31 = -t20 * t9 * t2 * t17 + t27 * t4
      t32 = t15 * t1 * t26
      t33 = t32 * t31
      t34 = t11 + t13
      t23 = t23 * t26
      t25 = t25 * t20 * t15 * t26
      t35 = t22 * t23 * t34 - t25
      t36 = t4 * t7 + t6 * t9
      t37 = t15 * t16
      t38 = t37 * t17 * t26
      t39 = t38 * t36
      t40 = t4 * t5 + t8 * t9
      t41 = t29 * t1
      t42 = t41 * t40
      t43 = t10 * t7 + t3 * t6
      t44 = t23 * t22
      t45 = t44 * t1 * t2
      t34 = -t20 * t41 * t34 + t45
      t46 = t15 * t17 * t26
      t47 = t34 * t35
      t48 = -t42 * t23 * t43 / 8
      t49 = t48 + t30 * t46 * (t24 * t10 - t27 * t3) + t47
      t50 = t24 * t6 + t27 * t7
      t51 = t46 * t50
      t36 = t23 * t36
      t26 = t37 * t26
      t31 = -t26 * t31
      t28 = -t46 * t28
      t37 = t23 * t40
      t6 = t20 * t7 * t2 * t17 + t27 * t6
      t7 = -t32 * t6
      t40 = -t37 * t41 * t43 / 8
      t43 = t7 * t26 * (t20 * t5 * t2 * t17 + t27 * t8) + t47 + t40
      t3 = t40 + t31 * t32 * (t20 * t3 * t2 * t17 - t27 * t10) + t47
      t10 = t12 + t14
      t32 = -t20 * t41 * t10 + t45
      t13 = t13 + t14
      t14 = t22 * t23 * t13 - t25
      t5 = t48 - t51 * t29 * (t24 * t8 + t27 * t5) + t47
      t8 = t11 + t12
      t11 = t22 * t23 * t8 - t25
      t10 = t22 * t23 * t10 - t25
      t12 = t44 * t16 * t17
      t13 = -t20 * t38 * t13 + t12
      t8 = -t20 * t38 * t8 + t12
      t12 = t29 * t50
      t3 = 0.1q1 / t3
      t18 = 0.1q1 / t18
      t22 = 0.1q1 / t2
      t15 = 0.1q1 / t15
      t5 = 0.1q1 / t5
      t23 = 0.1q1 / t17
      t24 = 0.1q1 / t43
      t25 = 0.1q1 / t16
      t27 = 0.1q1 / t49
      t29 = 0.1q1 / t1
      t38 = -t34 - t32 - t11
      t40 = t34 + t32 + t14
      t41 = -t35 - t10 - t13
      t43 = t35 + t10 + t8
      t44 = t32 + t14
      t45 = t30 ** 2
      t46 = t27 ** 2
      t48 = t27 * t46
      t49 = t24 ** 2
      t50 = t24 * t49
      t52 = t18 ** 2
      t53 = t51 ** 2
      t54 = t21 ** 2
      t55 = mt ** 2
      t56 = t5 ** 2
      t57 = t5 * t56
      t58 = t3 ** 2
      t59 = t3 * t58
      t60 = t34 * t48
      t61 = t40 * t48
      t62 = t10 * t21
      t63 = t55 * t22
      t64 = t63 * t15
      t65 = t64 * t23 * t25
      t66 = t65 * t46 * t29
      t16 = t16 * t17
      t17 = t16 * t19 ** 2
      t67 = t8 * t56
      t68 = t17 * t33 * t22
      t69 = t42 * t35
      t70 = t41 * t58
      t71 = t43 * t49
      t72 = t7 * t34
      t73 = t72 * t31
      t74 = t4 * t9
      t75 = t39 * t34
      t76 = t75 * t33
      t77 = t76 * t3
      t63 = t63 * t29
      t78 = t19 * t33
      t79 = t38 * t56
      t80 = -t40 * t35 * t46 - t79 * t35 - t5
      t81 = -t10 - t8
      t82 = t32 + t11
      t83 = t35 ** 2
      t84 = t35 * t83
      t85 = t31 ** 2
      t86 = t34 ** 2
      t87 = t81 * t57
      t88 = t28 * t56
      t89 = t42 * t53
      t90 = t72 * t85
      t91 = t83 * t28
      t92 = t13 * t48
      t93 = t51 * t30
      t94 = t93 * t31
      t95 = t7 * t31
      t96 = t86 * t31
      t97 = t96 * t58
      t98 = t30 * t24
      t99 = t82 * t49
      t100 = t31 * t21
      t101 = t78 * t16
      t102 = t19 * t52
      t76 = t102 * (t101 * (t80 * t51 * t45 + t99 * t90) * t29 * t22 + t
     #100 * (-t76 * t24 + (t98 + t97) * t7 * t28))
      t103 = -t35 - t10
      t104 = t103 * t49 - t70
      t105 = t69 * t51
      t106 = t1 * t2
      t107 = t106 * t23 * t25
      t108 = t44 * t58
      t109 = t93 * t35
      t110 = t109 * (t107 * (-t28 * t13 * t46 + (-t35 * t48 + t57 * t8) 
     #* t10 * t51 * t42) * t52 * t30 * t54 - t66 * t12 * t42 * t28 + (t1
     #07 * t46 * (t105 * t27 + t28) * t18 + t102 * (t46 * (t13 * t33 + t
     #14 * t28) - t67 * t33 + t61 * t105)) * t30 * t21)
      t101 = t102 * t90 * (t21 * (t104 * t33 + t28 * (-t99 + t108)) + t1
     #01 * (t34 * t49 - t108) * t29 * t22)
      t108 = -t56 + t46
      t111 = t34 + t32
      t112 = t3 - t24
      t113 = t5 - t27
      t114 = t108 * t83
      t115 = t30 * t3
      t116 = t96 * t49
      t117 = t35 * t10
      t118 = t69 * t19 * t45 * t53 * t57 * t111 + t19 * (t45 * (t114 * t
     #33 * t51 + t69 * (t56 * (t5 * (-t35 * t82 + t38 * t8) + 1) - t46) 
     #* t53) + t77 * t31 - t95 * (t115 + t116) * t28 + t94 * t33 * t113)
     # * t18 * t21 + t117 * t107 * t45 * t51 * (-t92 * t42 * t51 + t88) 
     #* t18 * t54
      t119 = t17 * t32
      t55 = (t33 * (t31 * (-t75 * t55 * t15 * t23 * t24 * t25 + t17 * t1
     #12 * t52 * t7 * t30) - t17 * t86 * t7 * t58 * t52 * t85) + t69 * t
     #45 * t53 * t56 * (t119 * t11 * t52 * t5 + t55 * t15 * t23 * t25)) 
     #* t29 * t22
      t1 = t18 * t118 + t31 * (-t78 * t73 * t52 * t21 * t8 * t49 + (t63 
     #* (-t74 * t7 * t24 + t77 * t15) + (-t73 * t28 * t18 * t49 * t21 + 
     #t72 * (t71 + t70) * t52 * t28 * t31 * t54) * t2 * t1) * t25 * t23)
     # + t45 * (t69 * (t19 * (t18 * (t11 * t57 - t61) + t62 * (t38 * t57
     # + t60) * t52) - t66 - t60 * t17 * t22 * t29 * t44 * t52) * t53 + 
     #t52 * ((-t10 * t46 + t67) * t25 * t23 * t28 * t35 * t54 * t2 * t1 
     #+ t68 * t27 * t29) * t51) + (t63 * t74 * (t95 * t3 - t93 * t5) + (
     #(t45 * (-t91 * t46 * t51 - t92 * t89 * t83) + t94 * t28 * t27 + t8
     #5 * t28 * t7 * t3) * t52 * t54 + (t45 * (t35 * (-t88 * t51 + t42 *
     # (t10 * t48 + t87) * t53) - t89 * t57 * t83) + t90 * t28 * t58) * 
     #t18 * t21) * t2 * t1) * t25 * t23 + t76 + t110 + t101 + t55
      t2 = -t58 + t49
      t38 = t107 * t21
      t55 = t18 * t19
      t60 = t31 * t5
      t61 = t107 * t54
      t76 = t61 * t18
      t77 = t34 + t32
      t89 = t32 ** 2
      t90 = t10 ** 2
      t94 = t8 ** 2
      t101 = t11 ** 2
      t110 = t7 ** 2
      t118 = t13 ** 2
      t120 = t30 * t29
      t6 = -t33 * t26 * t6
      t26 = t6 * t58
      t121 = t59 - t50
      t9 = t9 ** 2
      t122 = t93 * t69
      t123 = t33 * t34 * t36
      t124 = t17 * t29 * t22
      t125 = t124 * t31
      t126 = t102 * t21
      t127 = t61 * t52
      t128 = t127 * (t42 * t39 * t51 * t27 * (-t30 * t83 * t27 + t31) + 
     #t85 * ((-t90 - t83) * t50 * t34 + t70) * t110 * t37)
      t4 = t4 ** 2
      t129 = t14 ** 2
      t130 = t36 * t34 * t31
      t131 = t34 * t31
      t132 = t9 * t29 * t25
      t133 = t95 * (-t112 * t23 * t22 * t4 - t132 * t3)
      t4 = t93 * (t113 * t23 * t22 * t4 + t132 * t113)
      t132 = t107 * t34 * t37 * t121 * t110 * t85
      t134 = t126 * (-t130 * t28 * t3 + (-t30 * t8 - t31 * t82) * t31 * 
     #t49 * t110 * t37 + t122 * (t36 * t8 - t75) * t56)
      t111 = t125 * t52 * t3 * (t123 + (-t115 * t111 + t131 * (t129 + t8
     #9) * t58) * t110 * t37)
      t4 = t18 * (t37 * (t55 * t73 * t21 * (t39 * t11 * t49 + t26) + t31
     # * (t107 * t100 * ((t35 * t49 + (t118 + t83 + t90) * t34 * t59 - t
     #34 * t94 * t50) * t18 * t21 - t49 + t58) + t55 * (t100 * t58 * (t3
     #2 + t14) + t16 * (t120 * t49 - t31 * (t101 + t89) * t29 * t50) * t
     #22 * t34 * t19)) * t110) + (t68 * t36 * t18 * t29 * (-t5 + t27) + 
     #t76 * t56 * (t10 + t8) * t51 * t42 * t39 + (t55 * t46 * t77 - t107
     # * t56) * t51 * t42 * t39 * t21) * t35 * t30 + t55 * t100 * t36 * 
     #t34 * t28 * t24) + t126 * (-t122 * t39 * t56 * t82 + (t30 * t104 *
     # t110 - t6 * t72 * t49) * t37 * t31) + t128 + t95 * t9 * t29 * t25
     # * t24 + t125 * (-t123 * t24 + (t30 * (-t14 * t58 + t99) + t121 * 
     #t31 * t34 * t86) * t110 * t37) * t52 + t111 + t132 + t134 - t127 *
     # t85 * t37 * t110 * t49 * t81 + t133 + t4
      t9 = t10 + t13
      t68 = t46 * t9
      t72 = t30 * t35
      t99 = t72 * t39
      t4 = t4 + (t95 * t64 * t75 * t37 * t11 * t49 * t29 + (t120 * t67 *
     # t64 * t35 * t36 + t106 * ((t30 * (-t68 * t35 + t56 * t83) - t60) 
     #* t52 * t39 * t54 + t99 * t18 * t46 * t21)) * t51 * t42) * t25 * t
     #23 + t55 * t30 * t37 * t110 * t31 * t2
      t67 = t127 * t31 * t39
      t83 = t21 * t13
      t104 = t65 * t35
      t111 = t36 * t42
      t122 = t120 * t55 * t16
      t123 = t7 * t37
      t67 = (t21 * (-t107 * t75 * t31 * t2 * t18 + t102 * (-t116 * t39 +
     # t98 * t39 + t131 * t58 * (-t13 * t36 + t39 * t44))) + t65 * t36 *
     # t31 * t29 * t112 + t67 * (t70 * t34 + t71 * t34 + t3)) * t7 * t37
     # + t111 * (-t104 * t39 * t29 * t113 + t124 * t93 * t52 * t27 + t93
     # * (t35 * t83 * t46 * t52 + t35 * t56 * t18) * t19) + t123 * (-t13
     #0 * t65 * t13 * t29 * t58 - t67 * t24 + t126 * (t39 * (-t131 * t32
     # * t49 - t115 + t97) + t31 * t112 * t36)) + t111 * t55 * ((t113 * 
     #t31 + t30 * (t117 * t108 + t114)) * t18 * t21 + t122 * t80 * t22 -
     # t72 * t46) * t51
      t17 = t120 * t17 * t22
      t65 = t65 * t29
      t70 = t126 * t27
      t43 = t43 * t50
      t71 = t64 * t49 * t29
      t80 = t16 * t19 * t22
      t97 = t102 * t85
      t98 = t18 * t48
      t111 = t52 * t21
      t26 = -t65 * t26 * t95 * t37 + t111 * (-t19 * t58 + t19 * (-t10 * 
     #t11 - t8 * t82) * t50 - t107 * t62 * t13 * t59) * t110 * t37 * t85
      t83 = t97 * t59 * (-t80 * t32 * t29 + t83) * t110 * t37 * t86
      t47 = t30 * (t107 * (t41 * t46 * t52 * t31 * t54 + t72 * (t48 - t5
     #7)) + t55 * (t30 * t56 + (t31 * (t34 * t46 + t79) + t68 * t30) * t
     #18 * t21 + t122 * t34 * t56 * (-t47 * t5 + 1) * t22)) * t53 + t126
     # * t109 * t12 * t28 * t46
      t68 = t36 * t28
      t79 = t110 * t37
      t62 = t34 * (-t119 * t37 * t110 * t14 * t59 * t22 * t29 * t85 - t6
     #1 * t39 * t28 * t24 * t31) + t79 * t19 * (t59 * (-t80 * t14 * t29 
     #+ t62) - t43 * t21) * t85 * t86 + t72 * t19 * t21 * (t33 * t39 * t
     #5 + t68 * t27)
      t17 = t42 * t47 + t52 * t62 + t34 * t26 + t42 * (t93 * t39 * (-t65
     # * t5 + t70) + t18 * t30 * (t38 * (t108 * t31 + (t98 * t118 * t35 
     #+ t98 * t84) * t30 * t21) + t17 * t35 * t18 * (t48 * (t129 + t89 +
     # t86) + t57 * (-t101 - t89)) + t46 * t19 * (t100 * t18 * t44 - t30
     #)) * t53) + t30 * (t42 * (t39 * t27 * (-t70 * t35 * t14 + t65) * t
     #51 + t52 * (t17 * t56 * t82 + t61 * (-t103 * t56 * t31 + t72 * t90
     # * t48)) * t53) - t104 * t36 * t28 * t29 * t113) + t85 * ((t106 * 
     #(-t41 * t59 - t43) * t18 * t21 + t71) * t25 * t23 + t55 * (t11 * t
     #50 + (t50 * (t103 * t32 - t11 * t35) + t35 * t32 * t59) * t18 * t2
     #1)) * t110 * t37 * t34 + t97 * (t80 * t50 * t29 * t82 + t21 * t35 
     #* t59) * t110 * t37 * t86 + t99 * t52 * t21 * t27 * (t38 * t28 - t
     #78) + t83
      t26 = t79 * t85
      t6 = t42 * (-t109 * t66 * t14 * t39 + t52 * t45 * (-t124 * t46 * t
     #40 - t61 * t84 * t57) * t53) + t17 + t26 * ((-t64 * t58 * t29 + t1
     #06 * (t50 * (-t103 * t8 + t117) - t9 * t35 * t59) * t52 * t54) * t
     #25 * t23 + t55 * (t59 * ((-t14 * t41 + t32 * t9) * t18 * t21 - t32
     # - t14) + t32 * (t55 * t16 * t11 * t22 * t29 + 1) * t50 + t21 * t1
     #8 * t49)) * t34 + t111 * t5 * t30 * (t42 * (t5 * (t19 * t30 * t81 
     #+ t38 * (t31 * t8 + t72 * (-t90 - t94) * t5)) * t53 - t19 * (t28 *
     # t35 * t12 * t5 + t39) * t51) - t68 * t19 * t35) - t26 * t55 * t12
     #1 * t86 - t99 * t127 * t28 * t5 + t131 * t23 * t25 * (t106 * t54 *
     # t39 * t28 * t3 * t52 + t71 * t123 * t6)
      t9 = -t58 + t49
      t3 = t3 - t24
      t17 = t73 * t37
      t3 = t36 * (t19 * (t17 * t9 * t18 + t111 * (t17 * (t49 * (-t35 - t
     #10 - t8) + t58 * (t35 + t10)) + t69 * t39 * (t5 - t27))) + t65 * t
     #75 * t37 * t3 + t123 * t124 * (t30 * t3 + t131 * (t49 * (t32 + t11
     #) + t58 * (-t32 - t14)) + t96 * t9) * t52)
      t5 = 16
      ret = t20 * t126 * t75 * t36 * t37 * t112 + t5 * (t102 * t51 * t35
     # * t45 * t21 * (-t105 * t34 * t57 + (-t56 + t46) * t10 * t33 - t88
     # * (t34 + t32 + t11) + t28 * t77 * t46) + t45 * (t69 * t18 * (t38 
     #* (-t87 * t21 * t35 * t18 + t92) + t55 * (t48 * t21 * (t10 * t44 +
     # t13 * t40) + t16 * (-t32 * t14 * t48 + t82 * t57 * t34) * t29 * t
     #22 * t19)) * t53 + t18 * (t107 * t91 * t54 * t18 * t56 + t55 * t11
     #3 * t28 * t21 - t78 * t35 * t108) * t51) + t1 + t18 * t7 * t85 * (
     #t78 * (-t112 * t18 * t21 + t2 * t34) - t76 * t28 * t24) + t93 * t2
     #3 * t25 * (-t60 * t106 * t54 * t28 * t52 + t63 * (t88 * t69 * t12 
     #* t15 + t74 * t27))) - 8 * t6 + 4 * t4 - 2 * t67 + t3

      hjetmass_qqbgg_triangle_pmpm_s34_mhsq_s12 = ret/32q0/(0,1q0)
      return

      end function
