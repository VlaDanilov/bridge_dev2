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
 

      complex*32 function hjetmass_qqbgg_triangle_pmpp_s34_mhsq_s12
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          double precision mt
          real*16 cg

          cg = 1q0

      t1 = zb(i4, i3)
      t2 = za(i1, i3)
      t3 = zb(i3, i1)
      t4 = za(i2, i3)
      t5 = zb(i3, i2)
      t6 = za(i1, i4)
      t7 = zb(i4, i1)
      t8 = za(i2, i4)
      t9 = zb(i4, i2)
      t10 = t2 * t3
      t11 = t4 * t5
      t12 = t6 * t7
      t13 = t8 * t9
      t14 = t10 + t11 + t12 + t13
      t15 = za(i1, i2)
      t16 = za(i3, i4)
      t17 = zb(i2, i1)
      t14 = -4 * t15 * t16 * t17 * t1 + t14 ** 2
      t14 = cg * sqrt(t14) + t10 + t11 + t12 + t13
      t18 = (0.1q1 / 0.4q1)
      t19 = t14 ** 2
      t20 = t15 * t16
      t21 = t20 * t17 * t1
      t22 = t19 * t18 - t21
      t22 = 0.1q1 / t22
      t23 = t3 * t4 + t7 * t8
      t19 = t19 * t22
      t24 = t19 * t23
      t12 = t10 + t12
      t25 = t14 * t15
      t26 = t25 * t17 * t22
      t27 = (0.1q1 / 0.2q1)
      t28 = t19 * t18
      t29 = t28 * t15 * t17
      t30 = -t12 * t26 * t27 + t29
      t31 = t14 * t18
      t25 = t25 * t22
      t32 = t25 * (t27 * t2 * t17 * t1 - t31 * t9)
      t33 = t2 * t5 + t6 * t9
      t34 = t21 * t27 * t14 * t22
      t12 = t18 * t19 * t12 - t34
      t35 = -t25 * (t27 * t6 * t17 * t1 + t31 * t5)
      t36 = t27 * t4 * t17 * t1 + t31 * t7
      t37 = t14 * t16 * t22
      t38 = t30 * t12
      t39 = t35 * t37 * t36 + t38 - t24 * t26 * t33 / 8
      t40 = 2 * t15 * t17 + t14
      t21 = t27 * t14 ** 2 - 2 * t21
      t41 = 2 * t16 * t1 + t14
      t23 = t26 * t23
      t10 = t10 + t11
      t42 = t18 * t19 * t10 - t34
      t20 = t20 * t27
      t43 = t14 * t1 * t22
      t5 = t43 * (t20 * t5 + t31 * t6)
      t6 = t37 * t1
      t10 = t28 * t16 * t1 - t27 * t6 * t10
      t28 = t43 * (-t31 * t2 + t20 * t9)
      t37 = t20 * t7 + t31 * t4
      t22 = -t5 * t14 * t17 * t22 * t37 + t38 - t23 * t19 * t33 / 8
      t11 = t11 + t13
      t13 = -t27 * t26 * t11 + t29
      t11 = t18 * t19 * t11 - t34
      t18 = t25 * (-t27 * t8 * t17 * t1 + t31 * t3)
      t20 = t43 * (-t20 * t3 + t31 * t8)
      t26 = -t43 * t37
      t25 = t25 * t36
      t2 = t2 * t7 + t4 * t9
      t9 = t19 * t2
      t2 = t6 * t2
      t6 = 0.1q1 / t8
      t19 = 0.1q1 / t15
      t22 = 0.1q1 / t22
      t29 = 0.1q1 / t39
      t31 = 0.1q1 / t4
      t21 = 0.1q1 / t21
      t33 = 0.1q1 / t17
      t34 = t32 * t20
      t36 = t35 * t26
      t37 = t36 + t34
      t39 = t1 * t19
      t43 = t31 * t7
      t44 = t39 + t43
      t45 = -t20 * t7 + t26 * t3
      t46 = t30 + t42
      t47 = t10 + t11
      t48 = t41 ** 2
      t49 = t12 ** 2
      t50 = t12 * t49
      t51 = t1 ** 2
      t52 = t22 ** 2
      t53 = t22 * t52
      t54 = t23 ** 2
      t55 = t23 * t54
      t56 = t3 ** 2
      t57 = t21 ** 2
      t58 = t10 * t25
      t59 = t30 * t31
      t60 = t1 * t35
      t61 = t59 * t29
      t62 = t61 * t6
      t63 = t3 * t32
      t64 = t1 * t7
      t65 = t15 * t17
      t66 = t65 * t24
      t67 = t66 * t6
      t68 = t67 * t31
      t69 = t39 * t16 * t40
      t70 = t7 * t5
      t71 = t1 * t12
      t72 = t3 * t28
      t73 = t51 * t16
      t74 = t73 * t40
      t75 = t24 * t29
      t76 = t6 * t23
      t77 = -t12 - t10 - t11
      t78 = -t30 - t42 - t13
      t79 = t30 + t42 + t13
      t80 = t42 + t13
      t81 = t30 + t13
      t82 = t29 ** 2
      t83 = t29 * t82
      t84 = t24 ** 2
      t85 = t28 * t6
      t86 = t85 * t5
      t87 = t35 * t25
      t88 = t32 * t18
      t89 = t30 * t6
      t90 = t89 * t84
      t91 = t52 * t28 * t5
      t92 = t91 * t54
      t93 = t90 * t35
      t94 = t12 * t31
      t95 = t74 * t33 * t19
      t96 = t68 * t21 * t48 * (t92 * t77 + (-t34 * t10 + t36 * t77) * t8
     #2 * t30 * t24) + t31 * ((t33 * (t70 * t28 * t19 * t80 + (t70 * t78
     # + t72 * t81) * t6 * t12) + t86 * t79 * t21 * t41 * t24) * t52 * t
     #54 + t90 * t21 * t82 * t41 * (t12 * (t88 + t87) + t34 * t13)) * t4
     #0 * t16 * t1 + t95 * (-t93 * t82 * t80 + t94 * (-t30 - t13) * t52 
     #* t28 * t54)
      t97 = t12 + t10 + t11
      t98 = t3 * t12
      t99 = t1 * t28
      t100 = t16 * t1 * t40
      t101 = t100 * t6 * t31
      t102 = t101 * t21
      t103 = t21 * t31
      t104 = t41 * t11 * t21
      t105 = t100 * t21
      t106 = t29 * t6 * t84
      t107 = t106 * (-t88 * t100 * t57 * t41 * t31 + (t105 * t31 * (-t87
     # + t88 * (-1 + t104)) - t60 + (t65 * t37 * t31 + t60 * t97) * t21 
     #* t41) * t29 * t30)
      t108 = t88 + t87
      t109 = -t30 - t42 - t13
      t110 = t30 ** 2
      t111 = t28 ** 2
      t112 = t16 ** 2
      t113 = t40 ** 2
      t114 = t108 * t19
      t115 = t25 * t19
      t116 = t30 * t1
      t117 = t1 * t49
      t118 = t31 * t33
      t119 = (-t98 * t85 * t15 * t54 * t10 * t52 + t75 * (-t116 * t26 + 
     #t34 * t7)) * t21 * t41 - t3 * t54 * t111 * t52 - t90 * t65 * t34 *
     # t82 * (t12 + t11) * t57 * t48
      t120 = t114 * t51 * t84 * t110 * t6 * t82 * t31 * t57 * t33 * t113
     # * t112
      t90 = t105 * (t118 * t75 * (t7 * (t89 * t18 - t114) - t89 * t3 * t
     #25 + t116 * (t115 + t76)) + t90 * (-t60 * t30 * t19 * t33 + t103 *
     # (t88 * t10 + t36 * (t42 + t13)) * t41) * t82 + t118 * t52 * t54 *
     # (t19 * t111 * t3 * t109 - t117 * t109 * t6))
      t111 = t12 + t10 + t11
      t114 = t89 * t75
      t36 = (t28 * (t72 * t111 + t117) * t21 * t41 + t15 * t12 * t6 * (-
     #t70 + t72)) * t52 * t54 + t75 * t21 * (t36 * t41 * t7 + t114 * (-t
     #87 * t80 + t88 * (-t42 - t13)) * t33 * t21 * t19 * t113 * t112 * t
     #51)
      t34 = t31 * t36 + t31 * t119 + t54 * (t102 * t33 * (t70 - t72) * t
     #22 + t103 * (t41 * (t28 * (t71 * t47 - t70 * t97) + (-t1 * t50 + t
     #28 * (t17 * t24 * t5 - t98 * t11) + t70 * t12 * t47) * t6 * t15) +
     # t99 * (t70 * t30 * t19 + t98 * t42 * t6) * t33 * t40 * t16) * t52
     #) + t21 * t96 + t54 * (t74 * t28 * t19 * t31 * t21 * t33 * t22 + t
     #31 * (t28 * (-t71 * (t69 * t42 * t21 * t33 + 1) + t70) + t49 * (t1
     # + (-t1 * t47 + t70 - t72) * t21 * t41) * t6 * t15) * t52) + t75 *
     # (t33 * (t56 * t32 * t6 + t63 * t44 + t64 * (-t19 * t35 + t59)) + 
     #t41 * (t62 * t40 * t24 * t16 * t1 * (t35 * (t26 * t30 + t58) + t34
     # * t46) * t57 + t6 * (t59 * t15 * t45 - t60 * t24) * t21) + t68 * 
     #t37 * t57 * t48) - t76 * t71 * t15 * t41 * t21 * t24 * t31 * t22 +
     # t107 - t120 + t90
      t36 = 0.1q1 / t16
      t14 = 0.1q1 / t14
      t68 = t100 * t41 * t57
      t14 = mt ** 2 * t14
      t74 = t14 * t19
      t88 = t74 * t33
      t90 = t88 + t68
      t96 = t18 * t26 + t20 * t25
      t107 = t103 * t65
      t117 = t95 * t18
      t119 = t41 * t20
      t120 = t20 * t6
      t121 = t26 * t31
      t122 = t5 * t26
      t123 = t122 * t6
      t124 = t10 + t12
      t125 = t10 + t11
      t126 = t5 * t25
      t127 = t74 * t118
      t128 = t28 * t20
      t129 = t65 * t31
      t130 = t18 * t28
      t131 = t12 * t10
      t132 = t94 * t65
      t133 = t65 * t41
      t134 = t94 * t21
      t135 = t134 * t6
      t65 = -t88 * t94 * t86 * t52 + t134 * t86 * (t65 * t97 * t41 + t10
     #0 * t109) * t53
      t97 = t114 * t21 * (t107 * t48 * t20 * t26 + t117)
      t46 = t54 * (t101 * t57 * (t126 * t69 * t33 - t119 * t28) * t22 + 
     #t135 * (t133 * (t122 * (-t12 * t41 * t21 + 1) + t128) + t100 * ((t
     #28 * (t11 * t18 + t20 * (t30 + t13)) + t5 * (t26 * t46 + t58)) * t
     #21 * t41 - t130 - t126)) * t52) + t55 * t65 + t6 * (t54 * (t95 * t
     #21 * (t130 * t103 * t16 * t40 + t5) * t22 + t12 * (-t5 * (t127 * t
     #18 * t2 + t1) + t105 * (t39 * t33 * t5 * t78 + t103 * (t28 * (t124
     # * t18 + t20 * t42) + t126 * t12) * t41) + t129 * (-t122 * t125 + 
     #t128 * t77) * t57 * t48 + t1 * t5 * t125 * t21 * t41) * t52) + t55
     # * (-t91 * t71 * t16 * t40 * t41 * t31 * t57 + t132 * t57 * t28 * 
     #t5 * t48 * (t11 * (-t10 - t12) - t131) * t53) + t115 * t61 * t112 
     #* t51 * t113 * t57 * t24 * t18 * t33) - t103 * t72 * t69 * t23 * t
     #25 * t33 * t22 + t97
      t58 = t28 * t5
      t61 = t58 * t103 * t40 * t16
      t65 = t5 * t52 * t54
      t69 = t30 * t42
      t95 = t11 ** 2
      t97 = t10 ** 2
      t101 = t21 * t41
      t134 = t57 * t48
      t136 = t42 ** 2
      t137 = t13 ** 2
      t138 = t7 * t33
      t139 = t101 * t23
      t140 = t94 * t86 * t33 * t57 * t53 * t19 * t55 * t51 * (t137 + t11
     #0 + t136) * t113 * t112
      t43 = t101 * (t121 * t75 * t63 + (t1 * (-t24 * t5 * t6 - t94 * t26
     #) + t43 * (t122 + t128)) * t22 * t23)
      t2 = (t33 * (t72 * t44 + t85 * t56 + t64 * (-t19 * t5 + t94)) + t1
     #35 * t45 * t41 * t15) * t22 * t23 + t102 * (-t87 * t104 * t84 * t3
     #0 * t82 + t71 * t54 * t22 * t33 - t24 * t41 * t21 * (t130 + t126) 
     #* t22 * t23 + t91 * (t101 * t47 - 1) * t55) + t127 * t87 * t106 + 
     #t118 * t76 * t57 * t19 * t51 * (t75 * t108 + t92 * t109) * t113 * 
     #t112 + t132 * t86 * (1 + t134 * (t95 + t97)) * t53 * t55 + t43 - t
     #88 * t87 * t59 * t84 * t11 * t6 * t82 + t105 * (t75 * (-t115 * t11
     #8 * t63 + t76 * (-t103 * t37 * t41 + t60 * t19 * t33)) + t87 * t10
     #3 * t106 * t41 + t31 * t22 * t23 * (t18 * (-t138 * t28 * t19 + (t1
     #39 * t5 * t2 * t22 + t138) * t6 * t12) + t33 * (t19 * (t71 - t70) 
     #- t98 * t6) * t25)) + t140 + t134 * t129 * t76 * t22 * (t28 * (t20
     # * t24 + t65 * t50) + t122 * t24)
      t5 = t6 * t35
      t28 = t5 * t32
      t37 = t7 * t35
      t39 = t39 * t32
      t43 = t23 * t33
      t44 = t105 * t29
      t45 = t31 * t24
      t14 = t45 * (t139 * t99 * t22 + t5 * t33 * t19 * (t79 * t57 * t23 
     #* t113 * t32 * t112 * t51 + t14 * t30 * t20 * t9) * t82 * t24 + t9
     #3 * t21 * t32 * (t100 * (t101 * (t13 * t77 + t42 * t77 - t38) + t4
     #2 + t13) + t33 * t21 * t19 * (t13 * (t30 + t42) + t69) * t113 * t1
     #12 * t51 - t133 * t125) * t83) + t45 * (t6 * (t15 * (t101 * (-t38 
     #* t17 * t84 * t32 * t35 * t83 + t116 * t75 + t23 * t22 * (t70 - t7
     #2)) + t134 * t83 * t35 * t32 * t30 * t84 * t17 * (t11 * t124 + t13
     #1)) + t88 * t84 * t30 * t32 * t35 * t82) + t44 * (t84 * (t101 * t8
     #9 * t32 * t35 * t29 + t28 * t110 * (-t101 * t125 + 1) * t82) + t76
     # * t75 * t35 * t32 * (t101 * t77 + 1) + t43 * (t6 * (-t37 + t63) -
     # t39)))
      t38 = t1 * t110
      t45 = t32 * t19
      t17 = t75 * t17 * t32
      t47 = t32 ** 2
      t50 = t82 * t31 * t84
      t59 = t89 * t15
      t59 = t50 * (t3 * (t101 * t111 * t47 + t59 * t32) + t105 * (t89 * 
     #t63 * t42 + t37 * (-t89 * t13 + t45 * t81)) * t33 + t45 * t114 * t
     #35 * t57 * (t137 + t136) * t33 * t113 * t112 * t51 - t59 * t37)
      t5 = t59 + t50 * (-t3 * t47 + t134 * t67 * t35 * t32 * ((t97 + t49
     # + t95) * t29 * t30 - t12 - t10 - t11) + t105 * (t110 * (t1 * t80 
     #- t37 + t63) * t6 + t45 * (t63 * t78 - t38)) * t33 + t101 * (t32 *
     # (t116 * t10 + t37 * (-t10 - t11)) + t89 * (t37 * t111 + t63 * t77
     #) * t15)) + t31 * t29 * t84 * (t29 * (t32 * (-t116 + t37) + t89 * 
     #(t17 * t35 + t116) * t15) + t28 * t112 * t51 * t113 * t57 * t24 * 
     #t30 * t110 * t19 * t33 * t82 + t44 * (t33 * (t30 * (t6 * (t63 * t1
     #3 + t38) - t39 * t80) + t37 * t42 * (t45 - t89)) + t5 * t101 * (t3
     #0 * t20 * t9 + t24 * t32 * t80)) + t101 * (t32 * (t29 * (-t37 * t1
     #2 + t116 * (t12 + t11)) - t1) + (t35 * (t17 - t7) + t63 - t38 * t1
     #11 * t29) * t6 * t15))
      t9 = t18 * t6 + t25 * t31
      t15 = t94 * t76 * t96 * t22
      t17 = t24 * t9
      t3 = t33 * (-t3 * t7 * t36 + t74 * (t17 * t36 + t15)) + t41 * (t15
     # * t16 + t17) * t57 * t40 * t1 + t43 * t73 * t19 * t9 * t57 * t113
      t9 = t121 + t120
      t15 = -16
      t17 = 8
      ret = t27 * t5 + t15 * (t36 * (t33 * (t8 * t7 ** 2 * t31 + t4 * t5
     #6 * t6) + t66 * (t121 + t120) * t57 * t48) + t54 * (-t123 * t94 * 
     #t13 * t90 * t52 + t123 * t31 * t90 * t22) + t76 * t21 * t12 * (t11
     #5 * t103 * t51 * t112 * t113 * t18 * t33 + t117 + t119 * (t107 * t
     #41 * t26 - t1)) * t22 + t62 * t24 * (t68 * t96 + t88 * t96)) + t17
     # * (t21 * (t121 * t72 * t41 * t23 * t22 + t135 * t33 * t52 * t19 *
     # t54 * t113 * t112 * (t126 * t78 + t58 * (t13 * (-t30 - t42) - t69
     #) * t22 * t23 + t130 * t78) * t51 + t41 * (t12 * (t126 * t103 * t1
     #6 * t40 * t11 * t52 * t54 + t61 * (t10 * t79 + t11 * t79) * t53 * 
     #t55) + t49 * (t61 * t79 * t53 * t55 + t65) - t75 * t30 * t20) * t6
     # * t1) + t46) + 32 * t3 - 64 * t23 * (t40 * t1 * t57 * t41 * t9 + 
     #t88 * t9 * t36) + 2 * t34 - 4 * t2 + t14

      hjetmass_qqbgg_triangle_pmpp_s34_mhsq_s12 = ret/32q0/(0,1q0)
      return

      end function
