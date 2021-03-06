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
 

      double complex function hjetmass_triangle_ppmm_s234_0_mhsq_dp 
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

      t1 = zb(i2, i1)
      t2 = za(i3, i4)
      t3 = zb(i4, i2)
      t4 = za(i1, i3)
      t5 = za(i2, i3)
      t6 = zb(i3, i2)
      t7 = za(i2, i4)
      t8 = zb(i4, i3)
      t9 = t5 * t6
      t10 = t7 * t3
      t11 = t2 * t8
      t12 = t9 + t10 + t11
      t13 = za(i1, i2)
      t14 = zb(i3, i1)
      t15 = za(i1, i4)
      t16 = zb(i4, i1)
      t17 = t13 * t1
      t18 = t4 * t14
      t19 = t15 * t16
      t20 = t17 + t18 + t19
      if ( dreal(t20) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t20 = cg * cdsqrt(t20 ** 2) + t17 + t18 + t19
      t21 = 0.1D1 / t20
      t22 = -2 * t4 * t1 * t12 * t21 + t2 * t3
      t23 = -2 * t17 * t12 * t21 + t10 + t9
      t24 = t2 * t6
      t25 = -2 * t15 * t1 * t12 * t21 - t24
      t26 = t7 * t16
      t27 = t14 * t5 + t26
      t28 = -2 * t13
      t29 = t7 * t8
      t30 = t28 * t14 * t12 * t21 + t29
      t31 = t4 * t30
      t32 = t1 * (t13 * t23 + t31)
      t28 = t28 * t16 * t12 * t21 - t5 * t8
      t33 = t1 * t23
      t34 = t14 * t22
      t35 = t13 * (t33 + t34)
      t36 = t15 * t28
      t37 = t1 * (t36 + t31)
      t38 = t16 * t25
      t39 = t13 * (t38 + t34)
      t20 = 0.1D1 / t20
      t40 = 0.1D1 / t8
      t37 = 0.1D1 / t37
      t32 = 0.1D1 / t32
      t41 = 0.1D1 / t5
      t42 = 0.1D1 / t6
      t39 = 0.1D1 / t39
      t43 = 0.1D1 / t13
      t44 = 0.1D1 / t3
      t45 = 0.1D1 / t15
      t46 = 0.1D1 / t28
      t47 = 0.1D1 / t23
      t48 = t38 + t34
      t49 = t2 * t41
      t50 = t6 * t40
      t51 = -t50 - t49
      t52 = t49 * t42
      t53 = t52 + t40
      t54 = t7 * t41
      t55 = t54 * t42
      t56 = t55 + t44
      t57 = t11 * t44
      t58 = t57 + t7
      t59 = mt ** 2
      t60 = t22 ** 2
      t61 = t1 ** 2
      t62 = t1 * t61
      t63 = t43 ** 2
      t64 = t43 * t63
      t65 = t32 ** 2
      t66 = t32 * t65
      t67 = t4 ** 2
      t68 = t4 * t67
      t69 = t46 ** 2
      t70 = t46 * t69
      t71 = t37 ** 2
      t72 = t37 * t71
      t73 = t23 ** 2
      t74 = t23 * t73
      t75 = t3 * t15
      t76 = t4 * t15
      t77 = t63 * t32
      t78 = t33 * t43 * t65
      t79 = t64 * t47
      t80 = t4 * t44
      t81 = t80 * t61
      t82 = t17 * t39
      t83 = t60 * t44
      t84 = t83 * t47
      t85 = t76 * t42
      t86 = t40 * t71
      t87 = t41 * t1
      t88 = t61 * t23
      t89 = t88 * t12 * t20
      t90 = t77 * t1
      t91 = t45 * t46
      t92 = t1 * t37
      t93 = t74 * t37
      t94 = t93 * t45
      t95 = t23 * t71
      t96 = t86 * t15
      t64 = t64 * t32
      t97 = t16 * t73
      t98 = t97 * t71
      t99 = t27 * t4
      t34 = t99 * (t81 * t23 * (t76 * t14 * t65 * t63 + t98 * t69) * t20
     # * t12 + (t15 * (-t88 * t63 * t65 - t64 * t1) - t62 * t74 * t72 * 
     #t46 + t50 * t79) * t7 * t4 + t87 * (t33 * (-t96 * t48 + t80 * (t15
     # * ((-t34 - t33) * t43 * t65 - t77) - t95 * t46 * t48)) * t21 * t5
     #9 + t7 * t4 * (t15 * (t65 * (t18 * t33 * t12 * t20 * t63 - t89 * t
     #43) - t11 * t61 * t73 * t43 * t66 - t90 * t12 * t20) + t94 * t69 *
     # (t92 * t18 * t12 * t20 - t91 * t11))) * t42 - t15 * t2 * t62 * t7
     #3 * t72)
      t35 = 0.1D1 / t35
      t100 = 0.1D1 / t4
      t101 = t19 + t18
      t102 = t1 * t30
      t103 = t14 * t23
      t104 = t103 - t102
      t105 = t33 * t37
      t106 = t16 * t44
      t107 = t92 * t14 * t40
      t108 = t43 * t32
      t109 = t108 * t44
      t110 = t23 * t37
      t111 = t110 * t101
      t112 = t28 * t44
      t113 = t40 * t37
      t114 = t110 * t91
      t115 = t40 * t47
      t116 = t115 * t22
      t117 = t84 * t17
      t118 = t42 * t59
      t119 = 0.1D1 / t25
      t120 = 0.1D1 / t16
      t121 = t1 * t28
      t122 = t23 * t16
      t123 = -t122 + t121
      t124 = t15 * t23
      t125 = t43 * t40
      t126 = t52 * t1
      t107 = t126 * (t59 * (-t83 * t13 * t61 * t100 * t120 * t119 * t39 
     #- t82 * t22 * t100 * t40 + t125) + t61 * (t124 * t123 * t44 * t65 
     #+ t113 * t30 * (-t110 * t19 + 1) - t36 * t109) * t20 * t12 * t4 + 
     #t78 * t14 * t15 * t44 * t123 * t20 * t12 * t67) + t49 * t1 * (t4 *
     # (t105 * t106 * t40 * (-t111 + 1) + t23 * (t91 * t33 * t14 * t71 *
     # t44 * t104 * t67 + t19 * t1 * (t103 * t86 + t109) + t37 * (t107 *
     # t104 + t106 * (t23 * t45 * t104 * t69 + t105 * t104 * t46)) * t4 
     #- t107) * t42 + t113 * t112 * t61 * (t111 - 1)) * t20 * t12 + t118
     # * (-t117 * t100 * t35 + t33 * (-t113 + t80 * (-t108 - t114)) + t1
     #16 * t100))
      t111 = t92 + t91
      t123 = t92 * t44
      t127 = t46 * t37
      t111 = t127 * t20 * t12 * t27 * t73 * t67 * t61 * (t123 + (t57 * t
     #111 + t111 * t7) * t42 * t41)
      t128 = 0.1D1 / t12
      t129 = t52 * t8
      t130 = t6 * t44
      t131 = t54 + t130
      t132 = t2 ** 2
      t133 = t16 ** 2
      t134 = t14 ** 2
      t135 = t15 ** 2
      t136 = t44 * (t129 + 1)
      t137 = t20 * t12
      t138 = t137 * t41
      t139 = t138 * t44
      t140 = t12 ** 2 * t20 ** 2
      t141 = t140 * t40
      t142 = t22 * t44
      t143 = t49 * t44
      t144 = t137 * t73
      t145 = t18 * t43 + t1
      t146 = t80 * t16
      t147 = t97 * t4
      t148 = t55 * t1
      t149 = t135 * t133
      t150 = t67 * t134
      t151 = t40 * t15
      t152 = t80 * t42
      t153 = t73 * t72
      t154 = t15 * t42
      t155 = t154 * t68
      t156 = t140 * t41
      t157 = t116 * t59 * t21
      t50 = t99 * t1 * (-t157 * t43 * t41 * t44 + t137 * t23 * (-t129 * 
     #t76 * t61 * t65 * t43 * t44 + t88 * t85 * t41 * (t57 * t1 + t145 *
     # t7) * t66 + t88 * (t15 * (t146 * t50 + t49 * (t101 * t42 + t146))
     # + t55 * t4 * t101 * t46 * t23) * t72 + t148 * (-t75 * t1 * t40 + 
     #t147 * t69) * t71 + t147 * t55 * t37 * t45 * t70) + t156 * t61 * (
     #-t109 * t85 + t153 * (-t152 * (t150 + t149) * t46 * t23 + t151 * (
     #-t76 * t133 * t44 - t149 * t42 - t150 * t42)) - t155 * t134 * t73 
     #* t66 * t43 * t44))
      t101 = t54 * t40
      t127 = t127 * t67
      t147 = t62 * t13
      t149 = t1 * t71
      t95 = t89 * t99 * (t152 * t138 * t23 * (-t95 * t15 * t133 * t69 - 
     #t95 * t18 * t16 * t69 - t147 * t15 * t66) + t149 * (t4 * (t44 * t5
     #1 - t101) + t127 * t14 * t73 * t44 + t151 * t110 * (t16 * (t54 * t
     #4 + t15) + t18)))
      t50 = t95 + t50 + t1 * (-t140 * t67 ** 2 * t61 * t134 * t27 * t73 
     #* t41 * t44 * t40 * t72 + t92 * t99 * (t1 * (t15 * (-t110 * t53 * 
     #t20 * t12 - t141 * t41 * t42) + t97 * t86 * t55 * t3 * t12 * t20 *
     # t135) + t142 * t59 * t40 * t21 * t41) + t137 * (t16 * (t126 * t8 
     #* t74 * t71 * t44 * t69 + t136 * t61 * t15 * t74 * t72 * t46 + t12
     #9 * t94 * t44 * t70) + t61 * (t40 * (t75 * t55 * t14 * t73 * t72 -
     # t139 * t37) + t15 * t1 * t73 * t44 * t66) - t139 * t93 * t133 * t
     #42 * t70) * t27 * t67 - t118 * t132 * t3 * t100 * t41 * t40 * t128
     # + t144 * t14 * t61 * (t72 * (t131 * t40 + t143 * (t8 * t23 * t42 
     #* t46 + 1)) + t136 * t43 * t66 * t15) * t27 * t68)
      t94 = t153 * t62
      t95 = -t79 + t94
      t129 = t8 ** 2
      t133 = t45 ** 2
      t134 = t7 ** 2
      t136 = t153 * t61
      t70 = t93 * t42 * t133 * t70
      t93 = t73 * t71
      t138 = t80 * t1
      t139 = t27 * t67
      t150 = t25 * t42
      t152 = t3 * t134
      t153 = t43 * t47
      t158 = t142 * t35
      t159 = t40 * t39
      t35 = t47 * t25 * t35
      t9 = t139 * (t61 * (t124 * t55 * t65 * t63 * (t11 + t10) + t152 * 
     #t71 * t41 * t42 * t45 * t69 * t74) + t62 * (t152 * t72 * t41 * t42
     # * t46 * t74 + t43 * t66 * t15 * ((t132 * t129 * t44 + t152) * t42
     # * t41 + t9 * t44) * t73) - t152 * t79 * t40 * t41 + t70 * t87 * t
     #132 * t129 * t44)
      t152 = t2 * t1
      t25 = t118 * t87 * t27 * (t25 * t4 * (t40 * (t153 - t92) - t138 * 
     #(t108 + t114)) - t117 * (t1 * t120 * t39 + t35) + t60 * t47 * (t17
     # * (-t158 * t47 - t159) + t115) * t100 * t15) * t21
      t5 = t5 * (t139 * t40 * t44 * t95 * t6 ** 2 + t99 * (t151 * t95 + 
     #t138 * (t15 * (t33 * t63 * t65 + t64) + t61 * t74 * t72 * t46)) * 
     #t6) + t59 * (t150 * t116 * t87 * t27 * (t47 - t82) * t21 - t147 * 
     #t154 * t27 * t22 * t60 * t21 * t100 * t41 * t120 * t44 * t47 * t11
     #9 * t39) + t85 * t3 ** 2 * t134 * t40 * t41 * t95 * t27 - t126 * t
     #3 * t40 * t128 * (t1 * t7 + t14 * t2) + t139 * t1 * (t132 * (t136 
     #* t41 * t44 * t8 + t87 * t74 * t71 * t42 * t46 * (t92 + t91) * t44
     # * t129) + t7 * (t149 * t74 * t45 * t69 + t54 * (t136 * t40 + t154
     # * t64 + t70) * t3) + t33 * t57 * (t15 * t65 * t63 + t93 * t45 * t
     #69)) + t91 * t89 * t52 * t67 * t44 * t37 * (t103 - t102) + t25 + t
     #9 + t94 * t85 * t132 * t27 * t8 * t41 + t152 * t40 * t128 * (t49 *
     # t16 - t1)
      t9 = t75 * t42
      t25 = t9 + t4
      t57 = t91 * t67
      t6 = t27 * t2 * (t43 * (t126 * t75 * t32 + t115 * (t54 * t3 * t25 
     #+ t4 * t6 + t75) * t43) + t73 * (t87 * t67 * t58 * t69 * t133 * t3
     #7 + t57 * t61 * (t44 * (t49 * t8 + t6) + t54) * t71))
      t60 = t14 * t42
      t64 = t41 * (t106 - t60)
      t70 = t19 * t65 * t43
      t89 = t154 * t43
      t25 = t2 * (t40 * (t41 * (t16 * (t43 * (-t154 - t80) + t47 * (-t15
     #0 - t142)) + t153 * t25 * t27) - t98 * t4 * t131 * t1) + t90 * t76
     # * t56 * t28) + t87 * (t4 * (t93 * (-t106 + t60) + t33 * (-t30 * t
     #42 + t112) * t71) - t36 * t108 * t42) * t132 + t132 * t67 * t8 * t
     #73 * t37 * t41 * t42 * t44 * t104 * t69 * t133 + t57 * t143 * t105
     # * (t110 * t11 * t42 * t104 + t27) - t61 * t22 * t42 * (t101 * t39
     # + t35 * t44) * t13 - t88 * t4 * t42 * (t109 * t15 + t113 * t54)
      t17 = t4 * (t61 * (t73 * (t86 * t14 - t70 * t55) + t109 * t19 * t4
     #1 * t42 * t23) + t16 * t40 * t63 * (t54 + t130) + (-t86 * t30 + t5
     #5 * (-t86 * t3 * t30 + t36 * t65 * t43)) * t23 * t62 - t19 * t77 *
     # t33 * t55) + t1 * (t36 * (-t55 * t3 - 1) * t63 * t40 + t158 * t38
     # * t17 * t41 * t42) * t47 + t19 * t63 * t40
      t30 = t108 * t87 * t4 * (t154 * t33 * t8 * t44 * t32 * (-t122 + t1
     #21) - t27) * t132
      t3 = t2 * t17 + t61 * (t1 * (-t57 * t55 * t37 * t44 * t73 - t110 *
     # t80 * t40) + t40 * (t43 * (t154 + t80) + t47 * (t142 + t150))) + 
     #t1 * t25 + t2 * (t13 * (-t60 * t83 * t62 * t119 * t120 * t41 * t39
     # + t64 * t159 * t61 * t22) + t33 * (t15 * (t32 * (t52 * t16 * t43 
     #- t146 * t63) + t81 * t28 * t43 * t65) + t92 * t40 * t4 * (t92 * t
     #131 * t28 + t64)) + t89 * t41 * (t26 * t125 * t3 - t99 * t61 * t32
     #) + t92 * t55 * t67 * t73 * t104 * t69 * t133 + t4 * t61 * (t86 * 
     #t55 * t3 * t14 - t70 * t44) * t73 - t121 * t115 * t4 * t63 * t131 
     #+ t105 * t91 * t4 * (-t49 * t27 + (-t103 * t41 * t42 * t44 + t14 *
     # t37 * t56 * t73) * t4 * t1 - t31 * t88 * t37 * t56)) - t147 * t14
     #2 * (t148 * t22 * t119 * t120 + t40) * t39 + t30
      t13 = t87 * t27
      t9 = t13 * (t61 * (t67 * (-t60 * t141 * t97 * t135 * t72 + t140 * 
     #t96 * (t106 + t60) * t23) + t68 * (-t60 * t140 * t19 * t72 * t44 *
     # t46 * t74 + t140 * t103 * t44 * (t89 * t65 + t86) - t141 * t97 * 
     #t14 * t15 * t72 * t44) + t122 * t140 * t86 * t135 * t42 * t4) + t6
     #2 * (t140 * t124 * t65 * t42 * t44 * t67 - t155 * t140 * t14 * t73
     # * t66 * t44) + t157 * t2 * t128 * (t9 * t100 + 1))
      t1 = t127 * t156 * t88 * t44 * t42 * t27 * (t23 * (-t16 * t46 + t9
     #2 * (-t19 - t18)) + t1)
      t14 = t139 * t144 * t123 * t49 * t46 * (t91 * t16 + t92 * (t18 * t
     #45 + t16))
      t2 = t13 * (t57 * t137 * t105 * t2 * t44 + t151 * t118 * (t153 - t
     #92) * t21 * t22)
      t10 = t152 * t108 * t76 * t27 * (t33 * t32 + t43 + (t10 * t43 + t3
     #3 * (t11 + t10) * t32) * t42 * t41)
      t13 = t137 * t76 * t52 * t32 * t27 * t61 * (t145 * t32 * t23 - t43
     #)
      t16 = -32
      ret = t16 * (t34 + t27 * (t4 * (t62 * (t73 * (t72 * (t4 * (-t24 * 
     #t44 + t51 * t7) - t75 * t7 * t53) - t76 * t58 * t43 * t66) - t11 *
     # t4 * t72 * t46 * t56 * t74) - t29 * t52 * t4 * t61 * t74 * t45 * 
     #t69 * t71 + t79 * t75 * t7 * t40 + t81 * (t15 * (-t78 - t77) + t18
     # * t52 * t8 * t74 * t71 * t45 * t69) * t20 * t12) + t87 * (t40 * (
     #t85 * t63 + t84 * (t82 - t47) + t67 * t63 * t44) - t38 * t67 * t73
     # * t37 * t42 * t44 * t45 * t69 - t86 * t33 * t67 * t44 * t48) * t2
     #1 * t59)) - 8 * t107 + 96 * t111 - 64 * t50 + 16 * t5 + 12 * t6 + 
     #4 * t3 - 128 * t9 + 192 * t1 - 24 * t14 - 80 * t118 * t81 * t27 * 
     #t22 * t21 * t41 * (t108 * t15 + t110 * t46) + 48 * t2 - 28 * t10 +
     # 56 * t13

      hjetmass_triangle_ppmm_s234_0_mhsq_dp = ret/32d0/(0,1d0)
      return

      end function
