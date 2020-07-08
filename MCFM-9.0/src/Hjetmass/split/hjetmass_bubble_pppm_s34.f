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
 

      complex*32 function hjetmass_bubble_pppm_s34
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          
          complex*32 ret
          double precision mt

      t1 = za(i1, i4)
      t2 = za(i3, i4)
      t3 = zb(i3, i1)
      t4 = zb(i3, i2)
      t5 = zb(i4, i2)
      t6 = za(i1, i3)
      t7 = 0.1q1 / t6
      t8 = t1 * t7
      t9 = t8 * t5 + t4
      t10 = zb(i4, i1)
      t11 = 0.1q1 / t10
      t12 = t3 * t11
      t13 = t12 + t8
      t14 = za(i2, i3)
      t15 = t14 * t5
      t16 = t6 * t10
      t17 = t15 + t16
      t18 = za(i2, i4)
      t19 = t1 * t3
      t20 = t18 * t4
      t21 = t19 + t20
      t22 = t6 * t3
      t23 = t14 * t4
      t24 = t1 * t10
      t25 = t18 * t5
      t26 = -t24 - t25 + t22 + t23
      t26 = t26 ** 2
      t27 = 0.4q1 * t21 * t17 + t26
      t27 = sqrt(t27)
      t28 = -t24 - t25 + t22 + t23 + t27
      t29 = 0.1q1 / t17
      t30 = 0.1q1 / 0.2q1
      t31 = t30 * t28 * t29
      t32 = t8 + t31
      t27 = -t24 - t25 + t22 + t23 - t27
      t33 = t30 * t27 * t29
      t34 = t8 + t33
      t35 = zb(i4, i3)
      t36 = 0.1q1 / t14
      t37 = t18 * t10
      t38 = t37 * t36 + t3
      t39 = 0.1q1 / t5
      t40 = t4 * t39
      t41 = t18 * t36
      t42 = t41 + t40
      t43 = t41 + t31
      t44 = t41 + t33
      t45 = zb(i2, i1)
      t46 = za(i1, i2)
      t47 = t1 ** 2
      t48 = t1 * t47
      t49 = t7 ** 2
      t50 = t49 ** 2
      t51 = t7 * t49
      t52 = t46 * t45
      t53 = t52 * t3
      t54 = t8 * t14
      t55 = t3 * t4
      t56 = t3 * t5
      t57 = t4 * t10
      t58 = t57 + t56
      t59 = t58 * t18
      t60 = 0.2q1
      t61 = t18 ** 2
      t62 = t18 * t61
      t63 = t36 ** 2
      t64 = t63 ** 2
      t65 = t36 * t64
      t66 = t36 * t63
      t67 = t18 * (t16 * t61 * t63 - t19 + t41 * (-t24 + t22))
      t68 = t1 * t14
      t69 = t68 * t3
      t26 = 0.4q1 * t21 * t17 + t26
      t26 = sqrt(t26)
      t70 = -t24 - t25 + t22 + t23 + t26
      t70 = 0.1q1 / t70
      t71 = t60 * t21
      t72 = -t71 * t70 - t31
      t73 = t29 * (t28 - t27)
      t74 = -t31 * t10 + t3
      t75 = -t31 * t5 + t4
      t26 = -t24 - t25 + t22 + t23 - t26
      t26 = 0.1q1 / t26
      t71 = -t71 * t26 - t33
      t76 = -t33 * t10 + t3
      t77 = -t33 * t5 + t4
      t78 = -t56 * t11 + t4
      t79 = t3 ** 2
      t80 = t79 ** 2
      t81 = t3 * t79
      t82 = t11 ** 2
      t83 = t82 ** 2
      t84 = t11 * t82
      t85 = -t12 + t31
      t86 = -t12 + t33
      t31 = -t40 + t31
      t33 = -t40 + t33
      t87 = t12 + t41
      t88 = t40 + t8
      t89 = 0.1q1 / t17
      t90 = t89 ** 2
      t91 = t90 ** 2
      t92 = t89 * t91
      t93 = t89 * t90
      t94 = t53 * t29
      t95 = (-t57 - t56) * t90
      t96 = t22 + t23
      t97 = t3 * t96
      t98 = t28 * t89
      t99 = t52 - t24 - t25
      t100 = t28 ** 2
      t101 = t28 * t100
      t102 = t10 * t100 * t90
      t103 = -0.3q1 / 0.2q1
      t104 = 0.1q1 / 0.4q1
      t105 = -0.1q1 / 0.8q1
      t106 = 0.3q1 * t3 * t21
      t107 = -t30 * t28 * (t95 * t28 * t14 + t94) + t103 * t97 * t28 * t
     #29 + t104 * t102 * t99 + t105 * t10 * t101 * t93 * t17 + t98 * (t9
     #8 * t16 * t3 + t59) - t106 + t19 * t98 * t60 * t10
      t108 = t52 - t25
      t109 = t1 * t10 ** 2
      t110 = t3 * t10
      t56 = t57 + t56
      t111 = t56 * t18
      t112 = t27 ** 2
      t113 = t27 * t112
      t114 = t10 * t112 * t90
      t115 = 0.3q1 / 0.4q1
      t116 = 0.3q1 * t97
      t117 = t27 * t89
      t16 = t103 * t97 * t27 * t29 + t104 * t114 * t99 + t105 * t10 * t1
     #13 * t93 * t17 - t30 * t27 * (t95 * t27 * t14 + t94) - t106 + t117
     # * (t117 * t16 * t3 + t59) + t117 * t19 * t60 * t10
      t19 = -t41 + t8
      t94 = t52 - t22 - t23
      t95 = t28 * t29
      t97 = t24 + t25
      t106 = t100 * t90
      t118 = t106 * t17
      t119 = t104 * t118
      t120 = t60 * t21
      t121 = t30 * t95 * t94 + t98 * t97 + t119 - t120
      t119 = -t30 * t95 * t99 - t98 * t96 + t119 - t120
      t122 = t98 * t17
      t123 = t60 * t97
      t124 = t27 * t29
      t125 = t112 * t90
      t126 = t125 * t17
      t127 = t104 * t126
      t99 = -t30 * t124 * t99 - t117 * t96 - t120 + t127
      t120 = t30 * t124 * t94 + t117 * t97 - t120 + t127
      t96 = t60 * t96
      t127 = t18 * t38
      t128 = t117 * t17
      t129 = t6 * t18
      t130 = t52 * t18
      t29 = t130 * t29
      t131 = (t129 + t68) * t90
      t97 = t18 * t97
      t132 = 0.3q1 * t18 * t21
      t95 = t103 * t97 * t95 + t104 * t106 * t14 * (-t52 + t22 + t23) + 
     #t105 * t14 * t101 * t93 * t17 - t30 * t28 * (t131 * t28 * t10 + t2
     #9) - t98 * (t18 * (t98 * t15 - t22) - t69) + t132 + t23 * t98 * t6
     #0 * t18
      t106 = -t52 + t22
      t133 = t14 ** 2
      t68 = t129 + t68
      t129 = t68 * t3
      t134 = t14 * t18
      t135 = 0.3q1 * t97
      t29 = -t104 * t125 * t14 * t94 + t103 * t97 * t124 + t105 * t14 * 
     #t113 * t93 * t17 - t30 * t27 * (t131 * t27 * t10 + t29) + t132 - t
     #117 * (t18 * (t117 * t15 - t22) - t69) + t23 * t117 * t60 * t18
      t30 = -t57 * t39 + t3
      t57 = -0.1q1 / t31
      t94 = 0.1q1 / t1
      t72 = 0.1q1 / t72
      t97 = 0.1q1 / t18
      t103 = -0.1q1 / t32
      t104 = -0.1q1 / t33
      t73 = 0.1q1 / t73
      t105 = -0.1q1 / t85
      t124 = -0.1q1 / t34
      t131 = -0.1q1 / t44
      t132 = -0.1q1 / t43
      t71 = 0.1q1 / t71
      t136 = -0.1q1 / t86
      t137 = t101 * t72 * t70
      t138 = t137 * t132
      t139 = t76 * t71
      t140 = t139 * t113
      t141 = t45 * t39
      t142 = t77 * t71
      t143 = t142 * t113 * t26 * t124 * t136
      t144 = t137 * t75 * t103 * t105
      t145 = t45 * t7
      t146 = t91 * t73
      t147 = t146 * t21
      t31 = 0.1q1 / t31
      t32 = 0.1q1 / t32
      t85 = 0.1q1 / t85
      t148 = 0.1q1 / t19
      t149 = 0.1q1 / t35
      t33 = 0.1q1 / t33
      t86 = 0.1q1 / t86
      t150 = 0.1q1 / t87
      t151 = 0.1q1 / t13
      t44 = 0.1q1 / t44
      t19 = 0.1q1 / t19
      t152 = 0.1q1 / t42
      t34 = 0.1q1 / t34
      t87 = 0.1q1 / t87
      t153 = 0.1q1 / t88
      t88 = 0.1q1 / t88
      t43 = 0.1q1 / t43
      t154 = t152 ** 2
      t155 = t151 ** 2
      t13 = 0.1q1 / t13 ** 2
      t156 = t4 ** 2
      t157 = t4 * t156
      t42 = 0.1q1 / t42 ** 2
      t158 = t39 ** 2
      t159 = t158 ** 2
      t160 = t39 * t158
      t161 = t75 * t94
      t162 = t3 * t97
      t163 = t112 * t26 * t71
      t164 = t51 * t151
      t165 = t148 * t7
      t166 = t36 * t2
      t167 = t55 * t66
      t168 = t40 * t2
      t88 = t7 * t88
      t169 = t2 * t45
      t170 = t93 * t149
      t171 = t4 * t94
      t172 = t45 * t81 * t49
      t173 = t172 * t36 * t13 * t87 * t83
      t174 = t1 * t79 * t9 ** 2 * t51 * t155 * t82 * t89 * t32 * t34 * t
     #149
      t175 = t45 * t157 * t63 * t42 * t159 * (t168 * t35 * t94 * t89 * t
     #31 * t33 - t88)
      t153 = t2 * (t3 * (-t47 * t45 * t9 * t50 * t155 * t82 * t89 * t32 
     #* t34 + t164 * t1 * (t141 * t47 * t36 * t97 * t19 * t7 * t153 - t4
     # * t9 * t32 * t149 * t34 * t89 + t40 * t1 * (t8 * t10 + t3) * t36 
     #* t97 * t149 * t19 * t153) * t11) + t61 * (-t45 * t4 * t38 * t64 *
     # t154 * t158 * t89 * t43 * t44 + t167 * t94 * t39 * (t166 * t44 * 
     #t89 * t43 + t165 * (t41 * t5 + t4) * t11 * t149 * t150) * t152) + 
     #t169 * t80 * t35 * t49 * t97 * t13 * t11 * t83 * t89 * t85 * t86 +
     # t170 * ((-t40 * t74 ** 2 * t57 * t36 * t97 * t132 + t161 * t103 *
     # t7 * (t12 * t75 * t105 - t4) + t162 * t74 * t36 * t132) * t70 * t
     #72 * t100 + t163 * (t40 * t76 ** 2 * t36 * t97 * t104 * t131 - t12
     # * t77 ** 2 * t94 * t124 * t7 * t136)) * t73 * t21 - t165 * t171 *
     # t141 * t62 * t64 * t150 * t11 * t152 - t18 * t156 * t38 ** 2 * t6
     #6 * t154 * t158 * t89 * t43 * t44 * t149 - t173 + t174 + t175)
      t173 = t38 * t44
      t174 = t165 * t156
      t31 = t157 * t63 * t89 * t31 * t33
      t33 = t79 * t9
      t175 = t81 * t78 ** 2 * t49 * t13 * t84 * t89 * t85 * t86 * t149
      t30 = t2 * (t1 * (-t33 * t51 * t36 * t19 * t155 * t82 * t149 + t16
     #4 * t55 * t36 * t19 * t11 * t149) + t18 * (-t174 * t38 * t66 * t15
     #4 * t158 * t149 + t167 * t149 * t39 * (t173 * t89 * t43 + t165) * 
     #t152) + t30 * (t88 * t157 * t63 * t42 * t160 * t149 - t31 * (t171 
     #* t2 + t45) * t42 * t159) + t47 * (-t164 * t79 * t36 * t97 * t19 *
     # t11 * t149 + t45 * t3 * t50 * t36 * t19 * t155 * t82) + t61 * (-t
     #145 * t4 * t64 * t148 * t154 * t158 - t174 * t94 * t66 * t152 * t3
     #9 * t149) + (-t81 * t36 * t149 * t87 * t49 * t84 - t172 * t85 * t8
     #6 * t89 * t83) * t13 * t78 - t171 * t169 * t62 * t35 * t65 * t154 
     #* t158 * t89 * t43 * t44 - t162 * t169 * t48 * t35 * t7 * t50 * t1
     #55 * t82 * t89 * t32 * t34 + t163 * t170 * t21 * t73 * (t171 * t77
     # * t7 * t124 - t162 * t76 * t36 * t131) + t31 * t30 ** 2 * t42 * t
     #160 * t149 - t175) + t153
      t31 = 0.1q1 / t46
      t42 = t131 ** 2
      t46 = t77 * t16
      t88 = t46 * t7
      t149 = t88 * t11 * t124 * t136 + t76 * t29 * t63 * t42
      t153 = t162 + t171
      t157 = t73 ** 2
      t159 = t132 ** 2
      t160 = t44 ** 2
      t164 = t71 ** 2
      t167 = t43 ** 2
      t169 = t149 * t26
      t170 = t92 * t97
      t172 = t170 * t94
      t41 = -t141 * t62 * t94 * t64 * t152 * t89 * t43 * t44 + (-t2 * t3
     #8 * t67 * t65 * t160 * t90 * t39 * t167 * t154 + t64 * t2 * (t11 *
     # (t127 * t36 * t150 * t7 * t148 ** 2 - t165 * t3 * t150) + t36 * t
     #160 * t90 * t167 * (t10 * t67 - t38 * (-t60 * (t18 * (t25 + t23) +
     # t69) - t18 * (t10 * (t41 * t6 + t1) + t52)))) * t39 * t152) * t94
     # * t61
      t9 = t12 * t49 * t97 * (t2 * (-t47 * t9 * (0.3q1 * t55 * (t54 - t1
     #8) + t60 * t8 * (t54 * t58 - t59) + t8 * (t15 * t47 * t10 * t49 + 
     #t53 + t8 * (t52 - t25) * t10)) * t32 ** 2 * t155 * t34 ** 2 * t49 
     #* t90 * t11 + t79 * t13 * (t3 * t36 * t87 + t78 * t3 * (-t15 * t79
     # * t82 + t20 + t12 * (-t25 + t23)) * t85 ** 2 * t86 ** 2 * t90) * 
     #t84) - t48 * t45 * t49 * t151 * t89 * t32 * t34)
      t15 = t72 ** 2
      t20 = t75 * t107
      t48 = t20 * t157 * t105
      t8 = t2 * (t162 * t47 * (t60 * t8 * (t14 * t3 - t37) - 0.3q1 * t18
     # * t3 + t47 * t14 * t10 * t49) * t50 * t63 * t19 ** 2 * t155 * t82
     # + t94 * (-t40 * t61 * t127 * t65 * t148 * t150 * t152 * (t150 + t
     #152) + (t48 * t170 * t70 * t103 * t15 + (t48 * t92 * t103 ** 2 + (
     #t20 * t92 * t105 ** 2 + (-t107 * t5 - t75 * (-t98 * t108 * t10 + t
     #109 * t98 + t115 * t102 * t17 - t60 * (t98 * t56 * t14 + t111) + t
     #116 + t53 - 0.4q1 * t110 * (t98 * t6 + t1))) * t92 * t105) * t157 
     #* t103) * t70 * t97 * t72) * t21 * t101) * t7 * t11 - t81 * t3 * (
     #t12 * t14 + t18) * t49 * t63 * t97 * t13 * t87 ** 2 * t83) + t163 
     #* t21 * t31 * t93 * t73 * t153
      t1 = t2 * t8 + t2 * (t21 * (-t100 * t31 * t93 * t70 * t72 * t153 *
     # t73 + t172 * t2 * (t113 * (t169 * t164 + (t11 * (t88 * t136 * t12
     #4 ** 2 + (t136 * (-t5 * t16 - t77 * (-t108 * t89 * t27 * t10 + t10
     #9 * t27 * t89 - t60 * (t56 * t89 * t27 * t14 + t111) + t115 * t114
     # * t17 + t53 - 0.4q1 * t110 * (t6 * t27 * t89 + t1) + t116)) + t46
     # * t136 ** 2) * t7 * t124) - t63 * t42 * (t10 * t29 + t76 * (-t117
     # * t106 * t14 - t117 * t133 * t4 + t115 * t126 * t14 + t60 * (t117
     # * t68 * t10 - t129) + t130 + 0.4q1 * t134 * (t117 * t5 - t4) + t1
     #35))) * t26 * t71) + t137 * t159 * t63 * (-t74 * (-t98 * t106 * t1
     #4 - t98 * t133 * t4 + t115 * t118 * t14 + t60 * (t98 * t68 * t10 -
     # t129) + t130 + t135 + 0.4q1 * t134 * (t98 * t5 - t4)) + t95 * (t7
     #2 * t74 - t10))) * t157) + t4 * t41 + t9)
      t4 = t36 * t94
      t6 = t113 * t71
      t8 = t6 * t26
      t4 = t2 * (-t166 * t7 * (t66 * t61 * t156 * t94 * t148 * t154 * t1
     #58 + t47 * t79 * t51 * t97 * t19 * t155 * t82) + t147 * t45 * (t13
     #7 * (t97 * t103 * t7 + t4 * t132) - t8 * (t97 * t124 * t7 + t4 * t
     #131)) + t146 * t97 * t94 * t2 * ((t77 * t124 * t7 + t76 * t36 * t1
     #31) * t26 ** 2 * t164 * t113 - t101 * t70 ** 2 * t15 * (t75 * t7 *
     # t103 + t74 * t36 * t132)) * t21 ** 2)
      t9 = t2 ** 2
      t14 = -t121 + t119
      t14 = t14 * t97 * t74 + t161 * t14
      t16 = -t123 - 0.2q1 * t52 + t22 + t23 - t96 + t24 + t25
      t17 = t121 - t119
      t19 = t77 * t94
      t37 = t76 * t97
      t41 = t91 * t157 * t21
      t22 = -t123 - 0.2q1 * t52 + t22 + t23 - t96 + t24 + t25
      t5 = t90 * t2 * (t125 * t31 * (t94 * (t120 * t5 + t22 * t77) + t97
     # * (t10 * t120 + t22 * t76)) * t26 * t71 * t157 * t21 - t168 * t61
     # * t38 * t67 * t94 * t65 * t152 * t167 * t160 * (t44 + t43) + t140
     # * t2 * t21 * t29 * t94 * t63 * t97 * t93 * t26 * t131 * t42 * t15
     #7) + t41 * t2 * (t100 * (t31 * t14 * t70 * t15 + t31 * (t94 * (t16
     # * t75 + t17 * t5) + t97 * (t10 * t17 + t16 * t74)) * t70 * t72) +
     # t163 * t31 * (t99 * (t94 * (t142 - t5) + t97 * (t139 - t10)) - t7
     #1 * (t37 + t19) * t120) + t137 * t2 * t74 * t95 * t94 * t63 * t97 
     #* t89 * t132 * t159)
      t10 = t20 * t11 * t103 * t105 * t7 + t74 * t95 * t63 * t159
      t15 = -t99 + t120
      t16 = t37 * t15
      t15 = t19 * t15
      t14 = t93 * t157 * t31 * t21 * t2 * (t89 * t73 * t14 * t70 * t72 *
     # t100 + t71 * t26 * t27 * (t117 * (t15 + t16) * t73 + t16 + t15) +
     # (t17 * t97 * t74 + t161 * t17) * t70 * t72 * t28)
      t15 = t89 * t9
      t16 = t15 * (t47 * (t55 * t32 * t11 * t97 * t151 * t34 * t50 + t33
     # * t32 * t82 * t97 * t155 * t34 * t50) - t173 * t61 * t156 * t94 *
     # t64 * t154 * t158 * t43)
      t6 = t172 * t73 * t157 * t21 * t9 * (-t137 * t10 + t6 * t169)
      t17 = t40 * t18 * t9 * t94 * t64 * t152 * (t38 * t67 * t90 * t167 
     #* t160 + t165 * t127 * t150 * t11)
      ret = -t60 * t147 * t2 * (t145 * t94 * t11 * (t144 - t143) + ((t13
     #8 * (t40 * t74 * t57 - t3) + (-t40 * t76 * t104 + t3) * t26 * t71 
     #* t131 * t113) * t94 * t2 + t141 * (-t140 * t104 * t131 * t26 + t1
     #38 * t74 * t57)) * t97 * t36) - 0.16q2 * t1 - 0.8q1 * t4 - 0.6q1 *
     # t171 * t147 * t9 * t7 * t97 * (-t137 * t103 + t8 * t124) - 0.32q2
     # * t5 - 0.96q2 * t41 * t97 * t94 * t9 * (t10 * t70 * t72 * t100 + 
     #t163 * t149) + 0.10q2 * t147 * t12 * t9 * t7 * t94 * t97 * (t144 -
     # t143) + 0.128q3 * t14 - 0.12q2 * t16 - 0.64q2 * t6 - 0.48q2 * t17
     # - 0.4q1 * t30 + t172 * t73 * t21 * t45 * t35 * t9 * ((t11 * t103 
     #* t105 * t7 + t57 * t36 * t132 * t39) * t70 * t72 * t100 ** 2 - t1
     #12 ** 2 * t26 * t71 * (t36 * t39 * t131 * t104 + t11 * t7 * t124 *
     # t136)) - 0.20q2 * t15 * t80 * t78 * t49 * t97 * t13 * t83 * t85 *
     # t86

      hjetmass_bubble_pppm_s34 = ret/16q0*(0,1q0)
      return

      end function
