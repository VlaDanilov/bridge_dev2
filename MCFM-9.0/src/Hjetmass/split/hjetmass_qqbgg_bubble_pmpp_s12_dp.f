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
 

      double complex function hjetmass_qqbgg_bubble_pmpp_s12_dp 
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
      t2 = za(i2, i3)
      t3 = zb(i3, i1)
      t4 = za(i2, i4)
      t5 = zb(i4, i1)
      t6 = t2 * t3
      t7 = t4 * t5
      t8 = t6 + t7
      t9 = za(i1, i3)
      t10 = zb(i3, i2)
      t11 = za(i1, i4)
      t12 = zb(i4, i2)
      t13 = t11 * t12
      t14 = t9 * t10
      t15 = t13 + t14
      t16 = t9 * t3
      t17 = t2 * t10
      t18 = t11 * t5
      t19 = t4 * t12
      t20 = t16 - t17 + t18 - t19
      t20 = t20 ** 2
      t21 = 4 * t8 * t15 + t20
      t21 = cdsqrt(t21)
      t22 = t21 - t16 + t17 - t18 + t19
      t23 = 0.1D1 / t15
      t24 = (0.1D1 / 0.2D1)
      t25 = t24 * t12
      t26 = t25 * t22 * t23 + t5
      t27 = za(i3, i4)
      t28 = zb(i4, i3)
      t29 = 0.1D1 / t15
      t30 = t9 * t22
      t31 = t3 * t12
      t32 = t5 * t10
      t33 = t32 + t31
      t34 = t27 * t28
      t35 = t34 * t10
      t36 = t34 * t3
      t37 = t33 * t4
      t33 = (t11 * t33 + t35) * t29
      t38 = 2
      t15 = 4 * t8 * t15 + t20
      t15 = cdsqrt(t15)
      t20 = t15 - t16 + t17 - t18 + t19
      t39 = 0.1D1 / t9
      t40 = t2 * t39
      t41 = t24 * t22 * t23
      t42 = t41 - t40
      t43 = 0.1D1 / t10
      t44 = t3 * t43
      t45 = t41 + t44
      t20 = 0.1D1 / t20
      t46 = t38 * t8
      t47 = t46 * t20 + t41
      t21 = t21 + t16 - t17 + t18 - t19
      t48 = t23 * (t22 + t21)
      t49 = -t18 - t16 + t17
      t50 = t29 * t22
      t51 = t32 + t31
      t52 = t29 ** 2
      t53 = t29 * t52
      t54 = (t4 * t51 - t36) * t23
      t35 = t11 * t51 + t35
      t51 = t22 ** 2
      t55 = t52 * t51
      t56 = (0.1D1 / 0.4D1)
      t57 = t38 * t3 * t8
      t58 = -t50 * t3 * t49 - t24 * t22 * (-t30 * t3 * t10 * t52 + t54) 
     #+ t56 * t55 * t35 - t57
      t59 = zb(i2, i1)
      t15 = t15 + t16 - t17 + t18 - t19
      t60 = 0.1D1 / t11
      t61 = t4 * t60
      t62 = t24 * t21 * t23
      t63 = -t62 - t61
      t64 = 0.1D1 / t12
      t65 = t5 * t64
      t66 = -t62 + t65
      t15 = 0.1D1 / t15
      t46 = -t46 * t15 - t62
      t67 = t41 - t61
      t68 = t41 + t65
      t25 = -t25 * t21 * t23 + t5
      t69 = t9 * t21
      t70 = -t62 - t40
      t71 = -t62 + t44
      t62 = -t62 * t10 + t3
      t41 = t41 * t10 + t3
      t72 = -t40 + t61
      t73 = t40 + t44
      t74 = t40 + t65
      t75 = t61 + t44
      t76 = t61 + t65
      t77 = t4 * t10
      t78 = t11 * t3
      t79 = -t77 + t78
      t80 = t77 * t60 + t3
      t81 = t4 * t80
      t82 = t29 * t21
      t83 = t21 ** 2
      t84 = t52 * t83
      t49 = t82 * t3 * t49 + t24 * t21 * (t69 * t3 * t10 * t52 + t54) + 
     #t56 * t84 * t35 - t57
      t54 = t34 + t16 + t18
      t57 = t34 + t16 + t18
      t85 = t22 * t23
      t86 = t24 * t85 * t57 - t8
      t87 = t34 + t17 + t19
      t85 = t24 * t85 * t87 + t8
      t88 = t34 + t17 + t19
      t89 = t21 * t23
      t57 = -t24 * t89 * t57 - t8
      t87 = -t24 * t89 * t87 + t8
      t89 = t40 * t12 + t5
      t90 = t31 * t43 - t5
      t91 = t32 * t64 - t3
      t9 = t77 * t9
      t77 = t11 * (t34 + t17) + t9
      t92 = t78 * t2
      t93 = (-t16 + t34) * t4
      t23 = (t4 * (-t34 + t16) + t92) * t23
      t94 = t13 * t4
      t95 = t19 + t17 - t18
      t9 = t11 * (t34 + t17) + t9
      t96 = t38 * t4 * t8
      t97 = -t24 * t21 * (-t94 * t21 * t52 + t23) + t56 * t84 * t9 + t82
     # * t4 * t95 - t96
      t9 = -t50 * t4 * t95 + t24 * t22 * (t94 * t22 * t52 + t23) + t56 *
     # t55 * t9 - t96
      t23 = t4 ** 2
      t24 = t60 ** 2
      t56 = t60 * t24
      t6 = t4 * (-t14 * t23 * t24 + t6 + t61 * (t17 - t16))
      t94 = t3 * t5
      t95 = t3 ** 2
      t96 = t43 ** 2
      t98 = t43 * t96
      t99 = 0.1D1 / t63
      t100 = 0.1D1 / t2
      t48 = 0.1D1 / t48
      t27 = 0.1D1 / t27
      t46 = 0.1D1 / t46
      t101 = 0.1D1 / t67
      t47 = 0.1D1 / t47
      t102 = 0.1D1 / t4
      t103 = t54 + t88
      t104 = t48 ** 2
      t105 = t47 ** 2
      t106 = t101 ** 2
      t107 = t99 ** 2
      t108 = t62 * t46
      t109 = t12 * t100
      t110 = t25 * t100
      t111 = t41 * t47
      t112 = t46 * t15
      t113 = t112 * t21
      t114 = t53 * t104 * t8 * t1
      t115 = t46 ** 2
      t116 = t62 * t102
      t117 = t47 * t20
      t54 = t114 * t27 * (t21 * (-t87 * (t116 + t110) * t15 * t115 + t11
     #2 * (t102 * (t10 * t87 + t54 * t62) + t109 * t87)) + t117 * t86 * 
     #t22 * (t102 * (-t111 + t10) + t109)) + t114 * (t113 * (t102 * (t27
     # * (t57 * (-t108 + t10) + t62 * t88) + t82 * t1 * t62 * t97 * t99 
     #* t107 * t100 * t24) + t109 * t57 * t27 + t110 * t27 * (-t46 * t57
     # + t54 + t88)) + (t27 * (t100 * (t103 * t26 + t12 * t85) + t102 * 
     #(t10 * t85 + t103 * t41)) * t47 + t27 * (-t41 * t85 * t102 + (-t86
     # - t85) * t100 * t26) * t105) * t20 * t22 - t111 * t1 * t51 * t9 *
     # t24 * t100 * t102 * t29 * t20 * t101 * t106)
      t85 = t86 + t85
      t86 = t41 * t102
      t85 = t85 * t100 * t26 + t86 * t85
      t57 = -t116 * (t57 + t87) + t110 * (-t57 - t87)
      t87 = 0.1D1 / t45
      t88 = 0.1D1 / t42
      t103 = 0.1D1 / t70
      t109 = 0.1D1 / t71
      t114 = t41 * t9
      t118 = t114 * t24 * t106
      t119 = t26 * t58 * t39 * t87 * t88 * t43
      t120 = t62 * t97 * t107 * t24
      t121 = t25 * t49
      t122 = t121 * t39 * t103 * t43 * t109
      t123 = 0.1D1 / t76
      t63 = -0.1D1 / t63
      t124 = 0.1D1 / t75
      t125 = 0.1D1 / t72
      t67 = -0.1D1 / t67
      t126 = t100 * t5
      t127 = t102 * t3 + t126
      t128 = t63 ** 2
      t129 = t67 ** 2
      t130 = t121 * t15 * t109
      t131 = t117 * t88
      t132 = t131 * t87 * t51
      t133 = t117 * t22
      t134 = t102 * t52 ** 2 * t100
      t135 = t81 * t39 * t124 * t125 * t43
      t136 = t80 * t6 * t52 * t128 * t129
      t137 = t65 * t56 * t123 * t100
      t138 = t137 * t1 * (t136 * (t61 * (t63 + t67) - 1) - t135)
      t9 = t1 * (t8 * (t52 * t27 * (t113 * t127 - t133 * t127) * t48 + t
     #134 * t1 * (t24 * ((t47 * (t9 * t10 + t41 * (-t38 * t4 * (t12 * (-
     #t50 * t11 + t4) + t17 - t18) + t50 * t77 + t92 - t93)) - t114 * t1
     #05) * t106 * t20 * t51 + t112 * t107 * t83 * (-t62 * (-t38 * t4 * 
     #(t12 * (t82 * t11 + t4) + t17 - t18) - t82 * t77 + t92 - t93) + (t
     #108 - t10) * t97)) + (t83 * (t46 * (-t130 * t103 ** 2 + (t109 * (t
     #12 * t49 + t25 * (t38 * t3 * (t10 * (-t69 * t29 - t2) + t16 + t18)
     # - t33 * t21 + t36 - t37)) - t121 * t109 ** 2) * t15 * t103) - t13
     #0 * t103 * t115) + t132 * (-t26 * (t38 * t3 * (t10 * (t30 * t29 - 
     #t2) + t16 + t18) + t33 * t22 + t36 - t37) + (t26 * (t87 + t88 + t4
     #7) - t12) * t58)) * t43 * t39) * t104) + t138)
      t17 = 0.1D1 / t68
      t30 = 0.1D1 / t73
      t33 = 0.1D1 / t66
      t37 = 0.1D1 / t59
      t46 = 0.1D1 / t72
      t47 = t30 ** 2
      t49 = t123 ** 2
      t58 = t64 ** 2
      t69 = t5 ** 2
      t72 = t5 * t69
      t77 = t39 ** 2
      t93 = t39 * t77
      t97 = t24 * t4
      t86 = t1 * ((t133 * (-t44 * t26 ** 2 * t100 * t39 * t87 * t88 + t8
     #6 * t60 * t101 * (t65 * t41 * t17 - t3) + t126 * t26 * t39 * t88) 
     #+ t113 * (t116 * t3 * t99 * t60 + t110 * t39 * t103 * (t44 * t25 *
     # t109 - t5) - t65 * t62 ** 2 * t99 * t33 * t60 * t102)) * t37 * t5
     #2 * t48 * t8 - t1 * t39 * t60 * (t2 * t95 * t77 * t102 * t46 * t47
     # * t96 + t97 * t69 * t100 * t125 * t49 * t58))
      t71 = -0.1D1 / t71
      t45 = -0.1D1 / t45
      t70 = -0.1D1 / t70
      t42 = -0.1D1 / t42
      t75 = 0.1D1 / t75
      t73 = 0.1D1 / t73 ** 2
      t81 = t81 * t60 * t43
      t106 = t1 * t100
      t75 = t60 * t75
      t107 = t73 * t3
      t108 = t2 ** 2 * t28
      t88 = t39 * t88
      t110 = t100 * t60
      t114 = t39 * t103
      t116 = t112 * t83
      t62 = t62 * t99
      t121 = t44 * t102
      t11 = t121 * t77 * (t1 * (t107 * (t90 * t3 * (-t13 * t95 * t96 + t
     #7 + t44 * (-t19 + t18)) * t71 ** 2 * t45 ** 2 * t52 + t75 * (t75 *
     # t3 * (t78 * t43 + t4) - t3)) * t96 + t40 * t47 * (-t89 * (t38 * t
     #94 * (t40 * t11 - t4) + t40 * (t4 * (-t32 - t31) + t36 + t40 * t35
     #)) * t52 * t70 ** 2 * t42 ** 2 + (-t38 * t4 * t3 + t40 * t79) * t2
     #4 * t46 ** 2) * t43) - t108 * t39 * t30 * t29 * t70 * t42)
      t13 = t137 * t23 * t28
      t6 = t1 * (t106 * t102 * t53 * (-t51 * t20 ** 2 * t105 * (t41 * t6
     #0 * t101 + t88 * t26) + (t25 * t39 * t103 + t62 * t60) * t15 ** 2 
     #* t115 * t83) * t48 * t8 ** 2 + t7 * (t106 * t56 * (t39 * (t124 * 
     #(t81 * t125 ** 2 + (t60 * t79 - t3) * t43 * t125) + t81 * t125 * t
     #124 ** 2) + t52 * t128 * t60 * t129 * (-t10 * t6 + t80 * (t38 * (t
     #4 * (t14 * t60 + t12) + t18) * t4 + t4 * (t16 + t34) + t92))) * t6
     #4 * t123 + t106 * t24 ** 2 * (t136 + t135) * t64 * t49) + t53 * t2
     #8 * (t117 * (-t110 * t101 - t88 * t102) * t51 - t116 * (t114 * t10
     #2 + t110 * t99)) * t48 * t8 + t11 - t13 * t29 * t63 * t67)
      t11 = t1 ** 2
      t14 = t116 * t103
      t15 = t29 * t11
      t16 = 0.1D1 / t74
      t18 = -0.1D1 / t66
      t19 = -0.1D1 / t68
      t31 = 0.1D1 / t74
      t32 = 0.1D1 / t76 ** 2
      t34 = t107 * t75
      t35 = t28 * t64
      t36 = t106 * t59
      t41 = t40 * t28
      t66 = t24 * t58
      t65 = t65 * t24 * t123
      t25 = t14 * t25 * t109 + t132 * t26
      t26 = t126 * t1
      t68 = t95 * t90
      t10 = t95 * t90 ** 2 * t77 * t37 * t73 * t96 * t29 * t71 * t45 + t
     #77 * t3 * (t3 * (t29 * t70 * t37 * t42 * t89 ** 2 + t60 * t37 * t4
     #6 * t89) - t41 * t60 * t46) * t96 * t47 - t66 * t69 * t91 ** 2 * t
     #37 * t32 * t29 * t18 * t19 + t13 * t39 * t125 * t124 * t43 - t94 *
     # t77 * t37 * t60 * t46 * (t2 * (t40 * t10 + t3) * t16 * t64 * t102
     # + 1) * t43 * t30 - t106 * t72 * t91 * t24 * t32 * t64 * t58 * t29
     # * t18 * t19 + t53 * (-t28 * t100 * t39 * t43 * t25 + (-t116 * t62
     # * t33 * (t26 + t28) + t111 * (-t26 - t28) * t101 * t20 * t17 * t5
     #1) * t102 * t64 * t60) * t48 * t8 - t97 * t94 * (t61 * t12 + t5) *
     # t39 * t100 * t125 * t37 * t124 * t43 * t123 * t64 + t65 * t37 * (
     #t5 * (-t123 * t29 * t63 * t64 * t67 * t80 ** 2 + t39 * t123 * t64 
     #* t125 * t80) - t3 * t39 * t125) - t68 * t28 * t77 * t73 * t98 * t
     #29 * t71 * t45
      t3 = t1 * t10 + t1 * (t3 * (t77 * (-t1 * t59 * t95 * t28 * t102 * 
     #t73 * t96 ** 2 * t29 * t71 * t45 + t121 * t2 * t60 * t46 * t37 * t
     #30 - t34 * t90 * t37 * t96 + t34 * t28 * t98) + t110 * t1 * t8 * t
     #102 * t53 * t48 * (t117 * t51 * t101 + t116 * t99) - t108 * t93 * 
     #t60 * t102 * t46 * t30 * t43 * t16 * t64) + t4 * (t69 * t39 * t24 
     #* t100 * t125 * t37 * t123 * t64 + t5 * t28 * t39 * t56 * t125 * t
     #49 * t58) + t66 * (-t35 * t91 * t29 * t18 * t19 + (t37 * t91 + t35
     #) * t31 * t39) * t32 * t69 - t36 * t72 * t28 * t24 * t32 * t58 ** 
     #2 * t29 * t18 * t19 + t44 * t29 * t30 * t77 * (-t108 * t1 * t59 * 
     #t77 * t30 * t43 * t102 - t41 * t89 * t30 * t43 - t5 * t89 * t37) *
     # t42 * t70 + t65 * t29 * (-t35 * t36 * t23 * t123 * t24 + t61 * (-
     #t35 * t80 * t123 + t106 * t3) + t3 * t80 * t37) * t67 * t63)
      t4 = t60 * t64
      t4 = t134 * t48 * t59 * t28 * t8 * t11 * (-t117 * t22 * t51 * (t4 
     #* t101 * t17 + t88 * t43 * t87) + t112 * (t114 * t43 * t109 + t4 *
     # t99 * t33) * t21 * t83)
      ret = -t38 * t4 - 64 * t54 + 256 * t53 * t27 * t48 * t104 * t8 * t
     #1 * (t117 * t85 * t22 + t113 * t57) - 128 * t52 * t104 * t8 * t1 *
     # (t117 * (t1 * (t55 * t100 * t48 * t102 * (t119 - t118) + t50 * t1
     #00 * t102 * (-t119 + t118)) + t27 * t85) + t112 * (t1 * (t84 * t10
     #0 * t48 * t102 * (t122 - t120) + t82 * t100 * t102 * (-t122 + t120
     #)) + t27 * t57)) - 32 * t9 + 8 * t86 - 16 * t6 - 12 * t15 * (t2 * 
     #(t89 * t93 * t47 * t70 * t42 * t96 * t102 * t95 + t94 * t93 * t30 
     #* t70 * t42 * t43 * t102) + t126 * (-t7 * t80 * t56 * t49 * t58 * 
     #t63 * t67 + t8 * t39 * t102 * t52 * t48 * (t131 * t51 + t14))) - 2
     #0 * t121 * t15 * t39 * (t68 * t39 * t73 * t96 * t71 * t45 + t8 * t
     #100 * t52 * t48 * t25) - 4 * t3

      hjetmass_qqbgg_bubble_pmpp_s12_dp = -ret/16d0*(0,1d0)
      return

      end function
