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
 

      complex*32 function hjetmass_triangle_ppmm_s124_0_mhsq
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

      t1 = za(i1, i3)
      t2 = za(i2, i4)
      t3 = za(i3, i4)
      t4 = zb(i2, i1)
      t5 = zb(i3, i1)
      t6 = zb(i4, i2)
      t7 = za(i1, i2)
      t8 = za(i1, i4)
      t9 = zb(i4, i1)
      t10 = t7 * t4
      t11 = t8 * t9
      t12 = t2 * t6
      t13 = t10 + t11 + t12
      t14 = za(i2, i3)
      t15 = zb(i3, i2)
      t16 = zb(i4, i3)
      t17 = t1 * t5
      t18 = t14 * t15
      t19 = t3 * t16
      t20 = t19 + t17 + t18
      if ( real(t20) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t19 = cg * sqrt(t20 ** 2) + t17 + t18 + t19
      t20 = t10 + t11
      t21 = 0.1q1 / t19
      t22 = -2 * t17 * t13 * t21 + t20
      t23 = 2 * t3
      t24 = t2 * t4
      t25 = t23 * t5 * t13 * t21 - t24
      t26 = -2 * t14 * t5 * t13 * t21 + t2 * t9
      t27 = t5 * t22
      t28 = t1 * (t15 * t26 + t27)
      t29 = t14 * t4 - t3 * t9
      t30 = t14 * t6
      t31 = t1 * (-t16 * t25 + t27)
      t32 = 2 * t1 * t13
      t33 = t32 * t16 * t21 - t6 * t7
      t32 = -t32 * t15 * t21 + t6 * t8
      t34 = t1 * t22
      t35 = t14 * t32
      t36 = t5 * (t34 + t35)
      t34 = t5 * (-t3 * t33 + t34)
      t37 = t23 * t15 * t13 * t21 + t4 * t8
      t23 = -t23 * t13 * t16 * t21 + t11 + t12
      t38 = t3 * t6
      t39 = t1 * t4
      t40 = t39 + t38
      t41 = 0.1q1 / t2
      t28 = 0.1q1 / t28
      t19 = 0.1q1 / t19
      t42 = 0.1q1 / t7
      t43 = 0.1q1 / t8
      t44 = 0.1q1 / t9
      t45 = 0.1q1 / t15
      t46 = 0.1q1 / t26
      t31 = 0.1q1 / t31
      t47 = 0.1q1 / t16
      t48 = t2 * t15
      t49 = t48 * t43
      t50 = t49 + t5
      t51 = t31 ** 2
      t52 = t31 * t51
      t53 = mt ** 2
      t54 = t28 ** 2
      t55 = t28 * t54
      t56 = t3 ** 2
      t57 = t3 * t56
      t58 = t25 ** 2
      t59 = t25 * t58
      t60 = t19 ** 2
      t61 = t15 ** 2
      t62 = t46 ** 2
      t63 = t46 * t62
      t64 = t5 ** 2
      t65 = t5 * t64
      t66 = t1 ** 2
      t67 = t66 ** 2
      t68 = t1 * t67
      t69 = t1 * t66
      t70 = t14 ** 2
      t71 = t13 ** 2
      t72 = t41 * t42
      t73 = t55 * t42
      t74 = t73 * t43
      t75 = t5 * t41
      t76 = t75 * t56
      t77 = t6 * t44
      t78 = t58 * t15
      t79 = t15 * t43
      t80 = t6 * t42
      t81 = t15 * t42
      t82 = t75 * t43
      t83 = t82 * t44
      t84 = t70 * t59
      t85 = t3 * t54
      t86 = t85 * t25
      t87 = t47 * t31
      t88 = t7 * t25
      t89 = t88 * t43 * t46
      t90 = t44 * t4
      t91 = t8 * t42
      t92 = t80 * t44
      t93 = t3 * t47
      t94 = t65 * t58
      t95 = t94 * t41
      t96 = t55 * t58
      t97 = t96 * t5
      t98 = t29 * t5
      t99 = 0.1q1 / t13
      t100 = 0.1q1 / t5
      t101 = t90 * t7 * t43 + 1
      t102 = t4 ** 2
      t103 = t3 * t52
      t104 = t80 * t2
      t105 = t44 * t29
      t106 = t105 * t13
      t107 = t41 * t19 * t13
      t108 = t107 * t78 * (-t73 * t8 * t14 - t103 + t90 * ((t89 - 1) * t
     #55 * t14 - t103 * t7 * t43)) * t29 * t64 - t106 * t97 * t14 * t61 
     #* t19 * t43 * (t4 + t104)
      t109 = t19 * t13
      t110 = t53 * t25 * t41
      t111 = t110 * t21
      t112 = t79 * t3
      t113 = t28 * t42
      t106 = t113 * t98 * t3 * t44 * (t112 * t71 * t60 + t111) + t106 * 
     #t41 * t28 * t19 * (t10 * t14 * t59 * t43 * t45 * t63 + t109 * t56 
     #* t42) * t64
      t114 = t25 * t55
      t115 = t114 * t46
      t116 = t15 * t52 * t47
      t117 = t81 * t96
      t118 = t95 * t19 * t13
      t119 = t14 * t58
      t120 = t119 * t54
      t121 = t93 * t15 * t51
      t122 = t107 * t25
      t123 = t122 * (-t91 * t85 + t90 * (-t85 + (t120 * t62 - t121) * t4
     #3 * t7)) * t29 * t64 * t66
      t94 = t105 * t74 * t94 * t71 * t60 * t68 * t15
      t124 = t6 * t28
      t125 = t107 * t14
      t126 = t125 * t28
      t127 = t6 * t46
      t128 = t1 * t28
      t125 = t125 * t46
      t129 = t116 * t66 * t6
      t130 = t1 * t54
      t131 = t54 * t42
      t114 = t109 * t98 * t25 * t66 * (t18 * t96 * t17 * t41 * t46 + (t5
     # * (-t114 * t66 * t4 * t15 + t120 * t46 * (t1 * (t124 * t15 - t126
     # * t61) + t127)) + t64 * (-t129 * t25 + t130 * t46 * (t128 * t6 - 
     #t125) * t58) - t118 * t69 * t55 * t46 - t131 * t38 * t48) * t44 * 
     #t43)
      t13 = t114 + t1 * t106 + t67 * (-t117 * t13 * t19 * t64 * t29 + t1
     #18 * (-t116 * t101 + t115) * t29) + t69 * t108 + t98 * (t1 * (t30 
     #* t13 * t19 * t5 * t59 * t43 * t45 * t44 * t63 * t28 + t83 * (t87 
     #* t56 * t15 - t84 * t28 * t63) * t60 * t71) + t66 * (-t84 * t83 * 
     #t71 * t60 * t15 * t62 * t54 + t86 * (t44 * (-t79 * t4 - t80 * t5) 
     #- t81) * t19 * t13) + t69 * (t78 * (-t38 * t5 * t43 * t44 * t52 + 
     #t73 * (-t77 * t5 - t15) * t14) * t19 * t13 + t58 * t44 * t15 * (t7
     #0 * (t72 * t5 * t55 * t15 + t74 * t61) + t76 * t16 * t43 * t52) * 
     #t60 * t71) + t72 * t93 * t53 * t21 * t44 + t95 * t44 * (t79 * t52 
     #* t47 + t73) * t60 * t71 * t68 + t97 * (-t92 * t50 + t75 * (-t91 +
     # t90 * (t89 - 1))) * t19 * t13 * t67) + t53 * t2 * t3 * t102 * t42
     # * t43 * t100 * t44 * t99 + t123 + t94
      t16 = t14 * t41
      t19 = t3 * t43
      t68 = t86 * t71 * t60
      t74 = t79 * t51
      t83 = t3 * t41
      t84 = t53 * t3
      t36 = 0.1q1 / t36
      t34 = 0.1q1 / t34
      t89 = 0.1q1 / t32
      t91 = 0.1q1 / t14
      t33 = 0.1q1 / t33
      t94 = t3 * t26
      t95 = t14 * t25
      t106 = t19 * t22
      t108 = t41 * (t95 + t94) - t106
      t114 = t22 ** 2
      t118 = t26 * t41
      t123 = t45 * t46
      t132 = t123 * t25
      t133 = t132 * t28
      t134 = t18 * t25 * t41
      t135 = t118 * t74
      t136 = t109 * t25
      t137 = t34 * t33
      t138 = t137 * t83
      t139 = t123 * t107
      t140 = t131 * t109
      t141 = t84 * t43
      t142 = t107 * t64
      t143 = t22 * t43
      t144 = t1 * t25
      t86 = t90 * (t66 * (t142 * t143 * t119 * t85 * t46 + t109 * t81 * 
     #t86 * t5 * t14 * (-t143 + t118)) + t69 * (t125 * t64 * t59 * t43 *
     # t54 + t142 * t43 * (t27 * t123 * t85 + t18 * t51 * t47) * t58 + t
     #74 * t122 * t93 * t64 * t26) + t143 * t53 * t56 * t42 * t36 + t144
     # * t19 * (t18 * t109 * t87 * t75 - t113 * t53))
      t70 = t86 + t90 * (t1 * (t106 * t119 * t107 * t64 * t45 * t62 * t2
     #8 + t3 * (t109 * (t118 * t79 * t93 * t31 + t113 * t108) + t110 * t
     #43 * (t133 - t87)) * t5) + t66 * (t126 * t64 * t59 * t43 * t45 * t
     #62 + t136 * (t72 * t70 * t15 * t25 * t54 + t135 * t56 + t19 * (t13
     #4 * t51 - t113)) * t5) + t67 * (t139 * t65 * t59 * t43 * t54 - t14
     #0 * t64 * t58 * t43) + t69 * (-t109 * t120 * t79 * t5 * t42 + t140
     # * t25 * t108 * t64) + t141 * (t42 * (t22 * t100 * t33 - t47) + t5
     #6 * t114 * t91 * t41 * t89 * t36 - t138 * t114))
      t86 = 0.1q1 / t25
      t108 = t95 + t94
      t110 = t8 * t41
      t118 = -t77 - t110
      t119 = t45 ** 2
      t120 = t47 ** 2
      t125 = t47 * t120
      t126 = t95 * t90
      t140 = t76 * t25
      t145 = t77 * t43
      t146 = t2 * t29
      t147 = t75 * t14
      t148 = t15 * t25
      t149 = t148 * t5
      t150 = t90 * t42 * t3
      t151 = t42 * t86
      t152 = t32 * t43
      t153 = t22 * t41
      t154 = t3 * t22
      t155 = t128 * t90
      t156 = t155 * t64 * t58 * t43
      t157 = t79 + t75
      t158 = t29 * t4
      t159 = t105 * t17
      t160 = t1 * t43
      t161 = t72 * t56
      t162 = t66 * t102 * t25
      t163 = t90 * t17
      t164 = t163 * t57 * t114 * t43 * t91 * t41 * t89 * t36
      t165 = t4 * (t42 * (t15 * (-t94 * t86 - t14) - t147 * t8) + t77 * 
     #(-t14 * t5 * t42 + t79 * (-t17 * t108 * t31 - t2 * t14 * t42))) * 
     #t120
      t133 = t133 * t17 * (t140 * t145 + t90 * (t124 * t66 * t5 * t58 * 
     #t43 - t158) + t75 * t4 * (t128 * t22 * t25 + t105) * t3)
      t166 = t17 * t28
      t34 = t47 * (t74 * t30 * t90 * t66 * t5 * t58 + t159 * t102 * t31 
     #- t56 * t42 * t157) + t47 * (t79 * (t108 * t44 * t102 + t126 * t75
     # * t3 - t140) * t31 * t1 + t149 * t4 * (t16 * t25 * t101 + t94 * (
     #t145 + t41)) * t51 * t66 + t150 * (t147 + t79 * (t146 * t86 + t14)
     #)) + t5 * t4 * (-t134 * t1 * t31 + t94 * (-t1 * t15 * t41 * t31 + 
     #t151 * t118)) * t120 + t126 * t72 * t17 * t3 * t28 + t156 * (t83 *
     # t10 * t22 + t144 * t6) * t62 * t119 + t56 * (t152 * t42 + t153 * 
     #(t35 * t90 * t5 * t43 * t34 + t42)) * t33 + t123 * t90 * t43 * t28
     # * t41 * t58 * t64 * t66 * (t3 + t10 * (t154 + t144) * t28) + t150
     # * t27 * ((t160 - t16) * t36 * t3 - t12 * t66 * t25 * t54 * t43) +
     # t156 * (t39 * t88 * t41 + t38 * t22) * t62 * t119 + t164 + t165 +
     # t133 - t166 * t25 * (t162 * t43 * t44 * t28 + t161)
      t35 = t123 * t66
      t26 = t4 * (t25 * (t85 * t5 * t42 * t66 * (-t118 * t26 - t22) - t1
     #9 * t113 * t66 * t5 * t44) + t131 * t5 * t66 * (t16 * t8 - t1 + t7
     #7 * (-t160 * t2 + t14)) * t58 + t44 * t3 * (t42 * (-t49 * t6 * t26
     # * t86 * t120 + (-t153 - t152) * t33 * t14) + t98 * t47 * (t79 * t
     #1 * t31 + t151))) + t34 + t35 * t4 * t64 * t58 * t54 * (t143 * t38
     # * t44 + t144 * t41) - t56 * t5 * t43 * (t128 * t92 * t25 + t138 *
     # t22 * t32) + t162 * t44 * t5 * (t54 * (t108 * t41 - t106) + t135 
     #* t93 * t7) + t27 * t57 * (t72 + t145 * (t83 * t22 * t89 * t91 + t
     #42)) * t36
      t27 = t15 * t51
      t32 = t25 * t43 * t46
      t34 = t79 * t42
      t38 = t3 * t42
      t39 = t66 * t5
      t56 = t39 * t25
      t57 = t29 * t1
      t85 = t10 * t41
      t88 = t85 + t6
      t92 = t5 * t8
      t94 = t17 * t15
      t106 = t58 * t54
      t72 = t5 * (-t94 * t6 * t125 * t31 - t80 * t125 * (t92 + t48) * t8
     #6 + t58 * (t55 * (t5 * (t8 * (-t4 * t41 - t80) - t90 * t6) - t4 * 
     #(t145 * t2 + 1) * t15) - t116 * t5 * t88 + t115 * t5 * t88) * t69 
     #+ t149 * t66 * t6 * t120 * t51) + t72 * t44 * (-t106 * t69 * t64 +
     # t154 * t33) * t21 * t53 + t136 * t64 * t66 * (t27 * (-t17 * t41 *
     # t120 + t93 * (t145 + t41)) - t106 * t62 * (t145 * t17 * t45 + t16
     #))
      t23 = t29 * t72 + t5 * (t57 * (-t15 * (t109 * t75 * t3 * t120 * t3
     #1 + t104 * t96 * t66) + t5 * (t10 * t58 * (t25 * (t124 * t119 * t6
     #3 + t127 * t66 * t55 + t130 * t45 * (-t107 * t17 + t6) * t62) - t1
     #29) - t109 * t6 * t15 * t31 * t120 * (t56 * t31 + t3)) * t44 * t43
     #) + t44 * (t29 * (t121 * t82 * t66 * t25 * t23 + t58 * (-t81 * t54
     # + t75 * (t25 * t54 * t46 - t27 * t47)) * t43 * t69 - t112 * t17 *
     # t23 * t41 * t31 * t120 + t161 * t22 * t36) + t38 * t120 * (t79 + 
     #t75) * (t1 * t9 + t30) + t95 * t128 * t29 * (t82 * t25 * t45 * t62
     # + t128 * (-t34 + t75 * (t32 - t42))) * t37) * t21 * t53)
      t27 = t6 ** 2
      t30 = t102 * t7 ** 2 * t41
      t72 = t2 * t27
      t80 = t30 + t72
      t82 = t2 ** 2
      t8 = t8 ** 2
      t95 = t84 * t37
      t104 = t3 * t36
      t106 = t11 * t125
      t107 = t11 * t66
      t72 = t72 * t44
      t108 = t106 * t98 * t15 * t86 + t3 * t4 * t99 * (t90 * t14 - t3) +
     # (t8 * t9 * t41 + t72) * t125 * t86 * t29 * t64
      t7 = t97 * t29 * (t157 * t44 * t7 * t102 - t11 * t75 * t25 * t46 +
     # t34 * t82 * t27 * t44 + t72 * t5 * t42) * t69
      t11 = t166 * t44 * t43 * (t5 * (-t146 * t27 * t119 * t63 * t59 + t
     #132 * t83 * (-t53 * t29 * t37 * t21 + t154 * t109 * t4)) + t95 * t
     #29 * t21 * t42)
      t33 = t141 * t44 * (t40 * (t42 * (t22 * (t104 * t5 + t33) - t47 * 
     #t5) - t137 * t75 * t3 * t114) - t22 * (t81 * t100 * t33 + t76 * t3
     #6 * t91) * t29) * t21
      t7 = t42 * t108 + t17 * t29 * (t107 * t117 + t75 * (t15 * (t107 * 
     #t58 * t47 * t52 - t144 * t20 * t120 * t51 + t106 * t31) + t96 * t6
     #6 * t8 * t9 * t42)) + (t28 * (-t30 * t57 * t64 * t59 * t119 * t63 
     #+ t139 * t66 * t3 * t4 * t64 * t58) + t3 * (-t24 * t42 * t99 * t40
     # + (-t38 * t29 * t36 * t21 * t22 + t83 * (-t104 * t89 * t91 + t137
     #) * t21 * t29 * t114) * t15 * t53) + t64 * (t57 * t48 * t27 * t31 
     #* t125 - t148 * t66 * t6 * t29 * t51 * (t12 + t10) * t120 + t57 * 
     #t87 * (t78 * t80 * t51 * t66 + t95 * t41 * t21)) + t81 * t98 * t82
     # * t27 * t86 * t125 - t146 * t66 * t64 * t27 * t59 * t45 * t62 * t
     #54 - t69 * t64 * t59 * t29 * t80 * t46 * t55) * t44 * t43 + t45 * 
     #t62 * t29 * t59 * t64 * (-t85 * t101 - t6) * t54 * t66 + t7 + t11 
     #+ t33
      t8 = t128 * t64
      t1 = t158 * (t47 * (-t90 * t49 * t1 * t31 + t151 * (t77 * t2 * t50
     # + t48 + t92) * t47) + t58 * (-t8 * t44 * t88 * t62 * t119 + t35 *
     # t64 * (t44 * (-t85 - t6) - t110) * t54))
      t2 = t158 * t94 * t87 * (t144 * t31 - t47 + (-t12 * t47 + t144 * (
     #t12 + t10) * t31) * t44 * t43)
      t9 = t123 + t128
      t6 = t8 * t109 * t46 * t29 * t58 * t3 * (t128 * t41 + (t6 * t9 + t
     #85 * t9) * t44 * t43)
      t9 = t39 * t28
      t10 = t8 * t32 * t83 * t105 * t60 * t71 * (t25 * (t128 * t18 + t14
     # * t46 + t9) + t3)
      t8 = t105 * t3 * (t8 * t123 * t122 * t4 + t34 * (-t128 * t25 - t47
     #) * t21 * t53)
      t9 = t155 * t142 * t46 * t29 * t58 * (t123 * t14 + t128 * t14 + t9
     # * t45)
      ret = 64 * t13 - 128 * t105 * (t66 * (-t68 * t14 * t5 * t61 * t42 
     #* t43 - t3 * t15 * t41 * (t131 * t14 + t19 * t51) * t60 * t71 * t2
     #5 * t64) + t67 * (t65 * (t16 * t79 * t71 * t55 * t60 * t46 * t59 -
     # t78 * t41 * (t103 * t43 + t73 * t14) * t60 * t71) - t96 * t71 * t
     #60 * t14 * t64 * t61 * t42 * t43) + t69 * (-t68 * t79 * t64 * t42 
     #- t83 * (t74 * t47 + t131) * t60 * t71 * t25 * t65) + t84 * t4 * t
     #21 * t42 * t99 * (t49 * t100 + 1)) + 8 * t70 - 4 * t26 - 32 * t23 
     #+ 16 * t7 - 12 * t1 - 80 * t159 * t19 * t111 * (t25 * t46 * t28 - 
     #t87 * t15) - 56 * t163 * t109 * t79 * t31 * t29 * (t144 * t3 * t31
     # + t87 * t56 + t93) + 28 * t2 + 96 * t6 - 192 * t10 - 48 * t8 - 24
     # * t9

      hjetmass_triangle_ppmm_s124_0_mhsq = ret/32q0/(0,1q0)
      return

      end function
