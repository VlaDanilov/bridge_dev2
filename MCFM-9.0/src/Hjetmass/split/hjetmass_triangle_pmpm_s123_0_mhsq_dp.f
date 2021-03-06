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
 

      double complex function hjetmass_triangle_pmpm_s123_0_mhsq_dp 
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double complex ret
          double precision mt
          double precision cg

      t1 = za(i1, i3)
      t2 = za(i1, i4)
      t3 = zb(i2, i1)
      t4 = zb(i3, i1)
      t5 = za(i2, i4)
      t6 = za(i3, i4)
      t7 = t6 * t4
      t8 = t3 * t5 + t7
      t9 = za(i2, i3)
      t10 = za(i1, i2)
      t11 = zb(i3, i2)
      t12 = t10 * t3
      t13 = t1 * t4
      t14 = t9 * t11
      t15 = t12 + t13 + t14
      t16 = zb(i4, i1)
      t17 = zb(i4, i2)
      t18 = zb(i4, i3)
      t19 = t2 * t16
      t20 = t5 * t17
      t21 = t6 * t18
      t22 = t20 + t21 + t19
      if ( dreal(t22) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t22 = cg * cdsqrt(t22 ** 2) + t19 + t20 + t21
      t23 = 0.1D1 / t22
      t24 = -2 * t5 * t15 * t16 * t23 + t4 * t9
      t25 = -2 * t19 * t15 * t23 + t12 + t13
      t26 = t24 * t17
      t27 = t16 * t25
      t28 = t2 * (t27 + t26)
      t29 = -t11 * t6 + t2 * t3
      t30 = t2 * t4
      t31 = t11 * t5 + t30
      t32 = -2 * t2 * t15
      t33 = t1 * t11
      t34 = t5 * (t32 * t17 * t23 + t33)
      t35 = t16 * (t2 * t25 + t34)
      t36 = -2 * t6 * t15 * t16 * t23 - t3 * t9
      t37 = t36 * t18
      t38 = t2 * (t26 + t37)
      t32 = t32 * t18 * t23 - t10 * t11
      t34 = t16 * (t32 * t6 + t34)
      t39 = 0.1D1 / t11
      t40 = 0.1D1 / t18
      t41 = 0.1D1 / t16
      t38 = 0.1D1 / t38
      t42 = 0.1D1 / t3
      t43 = 0.1D1 / t25
      t44 = 0.1D1 / t9
      t45 = 0.1D1 / t36
      t28 = 0.1D1 / t28
      t46 = 0.1D1 / t10
      t22 = 0.1D1 / t22
      t47 = t21 + t20
      t48 = t20 * t41
      t49 = t48 + t2
      t50 = t24 ** 2
      t51 = t24 * t50
      t52 = t45 ** 2
      t53 = t45 * t52
      t54 = t43 ** 2
      t55 = t43 * t54
      t56 = t28 ** 2
      t57 = t28 * t56
      t58 = t20 * t40
      t59 = t4 * t39
      t60 = t4 * t5
      t61 = t2 * t24
      t62 = t6 * t16
      t63 = t62 * t50
      t64 = t4 * t42
      t65 = t60 * t42
      t66 = t61 * t39
      t67 = t44 * t42
      t68 = t67 * t38
      t69 = t68 * t16
      t70 = t18 * t39
      t71 = 0.1D1 / t32
      t34 = 0.1D1 / t34
      t72 = 0.1D1 / t6
      t73 = 0.1D1 / t2
      t35 = 0.1D1 / t35
      t74 = t32 ** 2
      t75 = t18 ** 2
      t76 = t38 ** 2
      t77 = t38 * t76
      t78 = t16 ** 2
      t79 = t25 ** 2
      t80 = mt ** 2
      t81 = t2 ** 2
      t82 = t81 ** 2
      t83 = t2 * t81
      t84 = t5 ** 2
      t85 = t5 * t25
      t86 = t5 * t36
      t87 = t6 * t24
      t88 = t5 * t18
      t89 = t5 * t16
      t90 = t41 * t43
      t91 = t49 * t28
      t92 = t2 * t28
      t93 = t84 * t17
      t94 = t2 * t78 * t50
      t95 = t5 * t39
      t96 = t40 * t45
      t97 = t96 * t16
      t98 = t59 * t42
      t99 = t5 * t80
      t100 = t46 * t44
      t101 = t100 * t4
      t102 = t98 * t21
      t103 = t6 * t78
      t104 = t75 * t36
      t105 = t103 * t15 * t40 * t38
      t106 = t105 * t42
      t107 = t58 * t78
      t108 = t98 * t20
      t109 = t70 * t65
      t110 = t25 * t34
      t27 = t101 * (t81 * (t106 * t52 * t22 * t51 - t98 * t93 * (t27 * t
     #76 + t37 * t56) * t22 * t15 * t24 + t5 * (-t103 * t25 * t76 * t42 
     #* t45 - t104 * t39 * t41 * t54 * t28 + t102 * t17 * t56) * t22 * t
     #15 * t50) + t83 * (t107 * t15 * t76 * t42 * t45 * t22 * t51 - t109
     # * t16 * t36 * t15 * t56 * t22 * t24 + t108 * t16 * t15 * t76 * t2
     #2 * t50) - t110 * t98 * t80 * t84)
      t27 = t27 + t101 * (t2 * (t39 * (t92 * t43 * t50 * t75 * (-t91 * t
     #86 + t87 * (t91 + t90)) + t64 * (t89 * (t85 - t61) * t38 + t21 * t
     #19 * t24 * (-t85 + t61) * t76 + t63 * t81 * t18 * t56 + t88 * (-t8
     #7 + t86) * t28)) + t94 * t76 * t42 * (-t93 * t25 * t40 + t61 * t6)
     # * t45 - t85 * t6 * t78 * t50 * t40 * t38 * t42 * t52) * t22 * t15
     # + t99 * (t2 * (t50 * (-t90 * t70 * t28 - t97 * t38 * t42) - t98 *
     # (t38 + t28) * t24) + t5 * (t16 * (-t95 * t40 * t43 * t73 * t35 * 
     #t74 - t59 * t40 * t42 * t35 * t32) - t5 * t79 * t72 * t42 * t71 * 
     #t34)))
      t36 = t8 * t16
      t37 = t4 ** 2
      t91 = t4 * t37
      t93 = t1 ** 2 * t37
      t98 = t93 * t44
      t111 = -t57 + t77
      t112 = t40 ** 2
      t113 = t16 * t32 * t40
      t114 = t59 * t18
      t115 = t110 * t42
      t116 = t42 * t46
      t117 = t116 * (t9 * t11 ** 2 + t98)
      t118 = t39 * t44
      t93 = t118 * (-t10 * t3 ** 2 - t93 * t46)
      t119 = t52 * t76
      t120 = t119 * t40
      t121 = t118 * t12
      t122 = t4 * t18
      t123 = t122 * t36
      t124 = t100 * t23
      t125 = 0.1D1 / t15
      t126 = t116 + t118
      t127 = t41 ** 2
      t128 = t8 * t44
      t129 = t124 * t66 * t65 * t80 * t31 * t38 + t122 * t83 * t8 * t50 
     #* (t116 * t14 + t121) * t77
      t130 = t46 * t41
      t131 = -t130 * t13 * t81 * t8 * t75 * t56 * t54 + t93 * t92 * t8 *
     # t75 * t127 * t55 + t93 * t83 * t8 * t75 * t57 * t43
      t132 = t2 * t46
      t133 = t6 * t44
      t134 = t67 * t80
      t135 = t42 * t39
      t136 = -t134 * t84 * t8 * t79 * t18 * t23 * t46 * t72 * t71 * t34 
     #- t5 * t126 * t125 * t37 - t135 * (t133 + t132) * t125 * t91
      t137 = t95 * t100 * t92 * t90 * t75 * (t7 * t15 * t22 + t80 * t8 *
     # t23) * t50
      t12 = t16 * t129 + t5 * t136 + t51 * t131 + t78 * (t45 * t51 * t8 
     #* t83 * (t12 * t44 + t117) * t77 + t96 * t68 * t61 * t5 * t46 * (-
     #t60 * t25 * t15 * t22 + t80 * t24 * t31 * t23) + t128 * t13 * t120
     # * t81 * t51) + t2 * (t109 * t100 * t80 * t8 * t24 * t23 * t28 + t
     #117 * t8 * t78 * t112 * t38 * t53 * t51) + t81 * (t100 * t96 * t65
     # * t15 * t22 * t78 * t50 * t38 + t8 * ((-t3 + t93) * t56 * t54 * t
     #41 * t75 + t120 * t78 * (t117 + t11)) * t51) + t83 * (-t14 * t8 * 
     #t51 * t75 * t46 * t43 * t57 + t123 * (-t121 * t57 + t116 * (t98 * 
     #t111 * t39 - t14 * t57)) * t50) + t124 * t84 * t80 * (-t115 * (t89
     # * t72 + t114) * t8 + t39 * t16 * (-t5 * t32 * t73 * t35 * (t113 *
     # t43 + 1) + t64 * (t35 * (-t25 - t113) + t110)) * t31) - t104 * t1
     #01 * t90 * t66 * t15 * t22 * t84 * t28 + t137
      t14 = t118 * t3
      t25 = t46 + t14
      t93 = t20 * t15
      t98 = t6 * t31
      t101 = t5 * t29
      t104 = t134 * t23 * t46
      t109 = t80 * t23 * t46
      t113 = t118 * t28
      t33 = t33 * t64
      t117 = t28 * t55
      t120 = t44 * t40
      t121 = t135 * t99
      t129 = t13 * t44
      t131 = t129 + t11
      t134 = t119 * t116
      t136 = t4 * t16
      t137 = t8 * t24
      t138 = t50 * t81
      t48 = t138 * (t137 * (-t134 * t107 * t11 + (t14 * t48 + t132) * t5
     #6 * t54 * t75) * t22 * t15 + t104 * t114 * t16 * (-t98 * t76 + t10
     #1 * (t56 - t76)) + (t16 * (t131 * t45 * t24 * t16 + t122) * t77 - 
     #t18 * (t3 * t18 * t24 * t43 + t136) * t57) * t8 * t2)
      t29 = t48 + t83 * (t123 * (-t104 * t39 * t56 + t13 * (t116 * t111 
     #+ t118 * t111)) * t50 + t8 * (t33 * t78 * t77 * t45 * t46 + (-t124
     # * t80 * t39 * t56 - t13 * t25 * t57) * t43 * t75) * t51) + t121 *
     # t91 * t125 * (t130 + t120) + (t78 * (-t98 * t68 * t80 * t23 * t52
     # * t46 * t40 + t33 * t8 * t38 * t53 * t46 * t112) - t117 * t128 * 
     #t59 * t1 * t3 * t75 * t127) * t51 * t2 + ((t8 * (-t133 * t15 * t22
     # + t64 * (-t93 * t44 * t22 + t11) * t46 * t40 * t1) * t52 - t104 *
     # (t101 + t98) * t45) * t76 * t78 + t113 * t43 * t75 * (t99 * t29 *
     # t23 * t28 * t46 + t90 * (t13 * (t93 * t46 * t22 - t3) * t28 - t10
     #9) * t8)) * t51 * t81
      t33 = t64 * t1
      t48 = -t33 * t46 - 1
      t80 = t118 * t48 - t116
      t98 = t100 * t84
      t101 = t138 * t4
      t104 = t13 * t46
      t107 = -t104 - t3
      t123 = t6 * t37
      t124 = t11 * t42
      t132 = t85 * t67 * t1
      t68 = t24 * (t132 * t81 * t91 * t46 * t39 * t76 + t95 * t68 * t81 
     #* t37 * t46) + t51 * (-t64 * t81 * t11 * t16 * t46 * t52 * t112 * 
     #t38 - t97 * t64 * t83 * t11 * t46 * t76) + t101 * (t30 * t80 + t85
     # * t97 * (t116 * t11 + t44)) * t76 + t98 * (t124 * t84 * t79 * t72
     # * t71 * t34 + t115 * t2 * t37 * t39 + t32 * t39 * t35 * (t84 * t3
     # * t32 * t73 * t43 + t123 * t42))
      t97 = t5 * t46
      t7 = t68 * t16 + t18 * (t86 * t67 * t1 * t81 * t91 * t24 * t46 * t
     #39 * t56 + t123 * t97 * t67 * t66 * t28) + t75 * (t118 * t92 * t90
     # * t50 * (t97 * t8 + t90 * t4 * (t87 * t107 + t86 * t3)) + t101 * 
     #t90 * (-t100 * t87 * t59 * t1 + t86 * t25) * t56) + t4 * (t18 * (t
     #138 * t7 * t80 * t56 + t92 * t39 * (t98 * t24 + (-t116 * t5 + t90 
     #* t24) * t8 * t4)) + t75 * (-t90 * t81 * t6 * t51 * t25 * t56 - t9
     #0 * t66 * t5 * t8 * t46 * t28) + t98 * t16 * (t42 * (t5 * (t2 * t7
     #1 * t34 * t72 * t79 + t110) + t61 * t38) + t95 * t32 * t35 + t95 *
     # t6 * t74 * t73 * t43 * t35)) + t94 * t64 * t46 * t38 * (t85 * t11
     # + t129 * (t85 - t61)) * t52 * t112 + t96 * t19 * t38 * t24 * (-t3
     #7 * t8 * t42 + (t116 * t84 * t11 * t24 + t101 * t48 * t38 + t65 * 
     #(t61 * t46 + t8)) * t44 * t16)
      t32 = t116 + t118
      t34 = -t20 - t21
      t35 = t15 ** 2
      t48 = t6 ** 2
      t68 = t22 ** 2
      t71 = t59 * t1
      t72 = t71 * t44
      t73 = t20 * t77
      t74 = t24 * t45
      t79 = t3 * t39
      t80 = t84 * t17 ** 2
      t91 = t36 * t67 * t46
      t94 = t8 * t50
      t95 = t119 * t67
      t97 = t121 * t37 * t23 * t125
      t98 = t20 * t24
      t101 = t100 * t35 * t68
      t62 = t94 * t18 * (t118 * (-t122 * t62 * t77 + t20 * (-t107 * t43 
     #* t24 * t18 + t136) * t57) + t18 * t46 * (t98 * t43 * t57 - t62 * 
     #t64 * t77)) * t22 * t15 + t101 * t94 * t70 * (t80 * t64 * t16 * t7
     #7 + (-t20 * t56 * t54 - t80 * t57 * t43 - t117) * t24 * t18)
      t106 = -t106 * t22 * t8 * t11 * t51 * t46 * t53 + t91 * (t48 * t16
     # * t51 * t38 * t53 - t114 * t84 * t28) * t68 * t35
      t107 = t75 * t36 * t51
      t110 = t101 * t2 * t82 * t8 * t78 * t51 * t75 * t39 * t43 * t57
      t14 = t137 * (t50 * (t117 * t14 * t75 * t41 - t134 * t103 * t11) +
     # t60 * t18 * t16 * (-t116 * t56 + t118 * (-t56 + t76))) * t22 * t1
     #5 * t81
      t115 = t1 * t50
      t13 = t137 * t22 * t46 * t15 * t81 * (t122 * t28 * (t98 * t19 * t4
     #2 * t56 + t118 * (t18 * (t138 * t1 * t16 * t56 * t43 + t115 * t41 
     #* t55 + t115 * t92 * t54) + t28 * t42 * t16 * (-t83 * t78 * t24 * 
     #t15 * t28 * t22 + t13 * t81 * t16 * t24 * t28 - t60 * t1))) + t95 
     #* t103 * t50 * (t93 * t22 - t13) + t67 * t59 * t19 * t15 * t22 * t
     #48 * t24 * t18 * t75 * t77)
      t13 = t13 + t2 * t106 + t82 * (t107 * t43 * t57 * t25 * t22 * t15 
     #- t107 * t101 * t39 * t54 * t56) + t83 * t62 + t2 * (-t33 * t105 *
     # t100 * t22 * t8 * t51 * t53 + t114 * t69 * t35 * t68 * t84 * t8 *
     # t46) + t81 * (t95 * t35 * t68 * t48 * t8 * t78 * t51 * t18 * t46 
     #+ t36 * t65 * t18 * t24 * t76 * t46 * (1 + t72) * t22 * t15) + t83
     # * (t94 * (t44 * (t79 * t75 * t24 * t54 * t56 + t74 * t34 * t77 * 
     #t78 - t73 * t114 * t16) + t116 * t16 * (t74 * t77 * (t11 * t34 + t
     #129 * t34) * t16 + t122 * (t72 * (-t20 * t111 - t21 * t77) - t73))
     #) * t22 * t15 + t91 * t50 * (-t80 * t114 * t57 + t74 * (t48 * t75 
     #+ t80) * t77 * t16) * t68 * t35) + t122 * t8 * t78 * t50 * t57 * t
     #126 * t22 * t15 * t82 - t97 * (t120 * t31 * t16 + t130 * t8 * t18)
     # - t110 + t97 * (t31 * t46 + t128) + t14
      t14 = t2 * t38
      t21 = t14 + t96
      t25 = t104 + t3
      t3 = t94 * t22 * t15 * t5 * t2 * (t75 * (-t113 * t41 * t25 * t54 +
     # t2 * (t118 * (-t104 - t3) - t46) * t56 * t43) + t45 * t38 * t78 *
     # (t14 * t44 + t116 * (t11 * t21 + t129 * t21)))
      t11 = t18 * t24
      t14 = t101 * t36 * t88 * t24 * t81 * (t102 * t76 + t42 * (-t59 * t
     #16 * t56 + (-t77 * t16 * t45 * t50 - t114 * t77 * t24) * t17 * t6)
     # * t2 + t26 * t39 * t57 * (t11 * t43 + t64 * t16) * t81 + t108 * (
     #-t56 + t76))
      t9 = t71 + t9
      t10 = -t33 - t10
      t17 = t61 * t8
      t9 = t17 * t4 * (t78 * (-t24 * t112 * t42 * t52 * t131 * t38 + t96
     # * t61 * (t10 * t44 - t124) * t76) + t28 * t18 * (t30 * t28 * (t11
     #6 * t9 + t39) + t11 * (t39 * t25 * t54 * t127 + t92 * t90 * (t46 *
     # t9 + t79))) + t19 * t4 * (t118 * t10 - t42) * t76)
      t10 = t61 * t60 * (t90 * t15 * t22 * t8 * t75 * t46 * t39 * t28 + 
     #t67 * (-t96 * t8 * t78 * t15 * t38 * t22 + t109 * (-t31 * t16 * t2
     #8 - t8 * t18 * t38) * t39))
      t5 = t17 * t101 * t5 * (-t61 * t78 * t47 * t45 * t42 * t76 + (t5 *
     # t78 * t45 - t87 * t78 * t52) * t42 * t38 + t28 * t43 * t39 * t75 
     #* (t81 * t16 * t24 * t28 - t5 + t61 * (t20 * t28 + t43)))
      t11 = -4
      ret = -24 * t30 * t22 * t15 * t8 * (t70 * (t2 * t18 * t50 * t41 * 
     #t54 - t65) * t46 * t28 + t66 * t18 * (t49 * t43 * t24 * t18 + t64 
     #* (t20 + t19)) * t46 * t56 + t69 * (-t63 * t52 * t40 + t60 * t39 +
     # t61 * (-t59 * t47 + (-t58 - t6) * t45 * t24 * t16) * t38)) - 8 * 
     #t27 - 80 * t100 * t99 * t2 * t50 * t23 * (t70 * t31 * t43 * t28 + 
     #t36 * t42 * t45 * t38) + 16 * t12 + 32 * t29 + t11 * (t37 * t2 * (
     #t2 * (t132 * t96 * t78 * t76 * t46 * t50 + t85 * t16 * t76 * t32 *
     # t24) + t100 * t86 * t90 * t1 * t50 * t39 * t28 * (t92 + t90) * t7
     #5 + t86 * t61 * t56 * t32 * t18 + t135 * t8 * (t28 * t4 + t38 * (t
     #89 * t44 - t4))) + t7) + 64 * t13 + 96 * t3 - 128 * t14 + 12 * t9 
     #+ 48 * t10 + 192 * t5


      hjetmass_triangle_pmpm_s123_0_mhsq_dp = ret/32d0/(0,1d0)
      return

      end function
