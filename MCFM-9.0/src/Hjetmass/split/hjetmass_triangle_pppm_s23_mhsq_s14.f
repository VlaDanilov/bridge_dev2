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
 

      complex*32 function hjetmass_triangle_pppm_s23_mhsq_s14
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          
          complex*32 ret
          double precision mt

          parameter (cg = 1q0)

      t1 = za(i1, i4)
      t2 = zb(i4, i1)
      t3 = za(i1, i2)
      t4 = zb(i2, i1)
      t5 = za(i1, i3)
      t6 = zb(i3, i1)
      t7 = za(i2, i3)
      t8 = zb(i3, i2)
      t9 = za(i2, i4)
      t10 = zb(i4, i2)
      t11 = za(i3, i4)
      t12 = zb(i4, i3)
      t13 = t3 * t4
      t14 = t5 * t6
      t15 = t9 * t10
      t16 = t11 * t12
      t17 = t14 + t15 + t16 + t13
      t17 = -4 * t1 * t7 * t8 * t2 + t17 ** 2
      t17 = cg * sqrt(t17) + t13 + t14 + t15 + t16
      t18 = 2 * t7 * t8
      t19 = t18 + t17
      t18 = t17 ** 2 / 2 - t18 * t1 * t2
      t20 = (0.1q1 / 0.4q1)
      t21 = t17 ** 2
      t22 = t1 * t7
      t23 = t22 * t8 * t2
      t24 = t21 * t20 - t23
      t24 = 0.1q1 / t24
      t25 = t11 * t6 + t4 * t9
      t21 = t21 * t24
      t26 = t21 * t25
      t27 = t17 * t20
      t22 = t22 / 2
      t28 = -t22 * t6 + t27 * t9
      t29 = t17 * t8 * t24
      t30 = t29 * t28
      t31 = t17 * t1
      t32 = t31 * t2 * t24
      t25 = t32 * t25
      t33 = t14 + t13
      t34 = t21 * t20
      t35 = t34 * t1 * t2
      t36 = -t32 * t33 / 2 + t35
      t37 = t27 * t11 + t22 * t4
      t38 = -t29 * t37
      t23 = t23 * t17 * t24 / 2
      t33 = t20 * t21 * t33 - t23
      t39 = t15 + t16
      t40 = t20 * t21 * t39 - t23
      t41 = t10 * t3 + t12 * t5
      t42 = t21 * t41
      t43 = t29 * (-t22 * t10 + t27 * t5)
      t22 = -t29 * (t22 * t12 + t27 * t3)
      t44 = 2 * t1 * t2 + t17
      t45 = t27 * t6 - t9 * t8 * t2 / 2
      t31 = t31 * t24
      t46 = t31 * t45
      t47 = t27 * t4 + t11 * t8 * t2 / 2
      t48 = t31 * t47
      t49 = t31 * (t27 * t12 + t3 * t8 * t2 / 2)
      t50 = t17 * t7 * t24
      t45 = t50 * t45
      t27 = t31 * (t27 * t10 - t5 * t8 * t2 / 2)
      t24 = t17 * t2 * t24
      t28 = t24 * t28
      t31 = t32 * t41
      t14 = t14 + t16
      t16 = t20 * t21 * t14 - t23
      t41 = -t50 * t47
      t32 = -t32 * t39 / 2 + t35
      t24 = t24 * t37
      t13 = t15 + t13
      t15 = t20 * t21 * t13 - t23
      t23 = t29 * t7
      t29 = t34 * t7 * t8
      t14 = -t23 * t14 / 2 + t29
      t13 = -t23 * t13 / 2 + t29
      t5 = t10 * t11 + t4 * t5
      t10 = t23 * t5
      t3 = t12 * t9 + t3 * t6
      t12 = t23 * t3
      t5 = t21 * t5
      t23 = 0.1q1 / t2
      t29 = 0.1q1 / t11
      t34 = 0.1q1 / t45
      t35 = 0.1q1 / t33
      t37 = 0.1q1 / t1
      t39 = 0.1q1 / t9
      t45 = 0.1q1 / t27
      t18 = 0.1q1 / t18
      t47 = t36 * t29
      t50 = t49 * t37
      t51 = t8 * t36
      t52 = t4 * t49
      t53 = t52 + t51
      t54 = t36 ** 2
      t55 = t13 ** 2
      t56 = t35 ** 2
      t57 = t44 ** 2
      t58 = t8 ** 2
      t59 = t45 ** 2
      t60 = t26 ** 2
      t61 = t26 * t60
      t62 = t34 ** 2
      t63 = t7 ** 2
      t64 = t8 * t49
      t65 = t8 * t54
      t66 = t65 * t29
      t67 = (-t50 - t47) * t6
      t68 = t36 * t38
      t3 = t68 * t21 * t3
      t21 = t26 * t49
      t69 = t21 * t15
      t70 = t21 * t35
      t71 = t52 * t45
      t72 = t45 * t19
      t73 = t4 * t29
      t74 = t73 * t1
      t75 = t50 * t47 * t63 * t58 * t57 * t18 * t26 * t23
      t76 = t18 * t39
      t77 = t76 * t35 * t34
      t78 = t77 * t61
      t79 = t15 ** 2
      t80 = t1 * t19
      t81 = t80 * t2
      t82 = t81 * t26
      t83 = t49 * t19
      t84 = t83 * t13
      t85 = t47 * t7
      t86 = t80 * t6
      t53 = t78 * (t34 * (t45 * (-t75 * t15 * t35 + (t70 * t47 * t19 * t
     #13 * t18 + t67 * t23 * t15) * t44 * t8 * t7 + t84 * (t82 * t29 * t
     #18 + t6)) + t52 * (t85 * t8 * t44 * t15 * t23 - t84) * t59) - t86 
     #* t29 - t75 * t79 * t62 * t59) + t78 * (t72 * (t36 * (-t72 * t21 *
     # t1 * t2 * t29 * t62 * t18 * t55 + (-t64 * t45 + (-t71 + t6) * t29
     # * t1) * t34 * t13) + t49 * (t74 + t8) - t66 * t1 * t13 * t34 * t4
     #5) + t8 * (t23 * (t35 * (t45 * (t47 * t53 + t50 * t53) + t67) + t5
     #9 * (t4 * t49 ** 2 * t37 + t64 * t36 * t37 + t66) * t34 * t15) + t
     #29 * ((-t69 + t3) * t45 * t34 - t70) * t18 * t19) * t44 * t7 - t75
     # * t56)
      t75 = 0.1q1 / t36
      t28 = 0.1q1 / t28
      t78 = 0.1q1 / t43
      t84 = t49 * t39
      t87 = t27 * t29
      t88 = t87 + t84
      t89 = t25 * t49
      t90 = t36 * t46
      t91 = t90 + t89
      t92 = -t90 - t89
      t93 = t25 ** 2
      t94 = t93 ** 2
      t95 = t25 * t93
      t96 = t19 ** 2
      t97 = t18 ** 2
      t98 = t34 * t15
      t99 = t39 * t45
      t100 = t46 * t37
      t101 = t25 * t29
      t102 = t8 * t37
      t103 = t6 * t39
      t104 = t30 * t36
      t105 = t104 * t35
      t106 = t76 * t26
      t107 = t106 * t29
      t108 = t18 * t34
      t109 = t108 * t26
      t110 = t60 * t58
      t111 = t47 * t8
      t112 = t8 * t29
      t113 = t6 * t38
      t114 = t4 * t30
      t115 = t38 * t49
      t116 = t93 * t22
      t117 = t116 * t75
      t118 = t117 * t28
      t119 = t34 * t26
      t120 = t36 * t48
      t121 = t120 * t49
      t122 = t97 * t39
      t123 = t122 * t37
      t124 = t122 * t26
      t125 = t124 * t29 * t2 * t1
      t126 = t35 * t26
      t127 = t4 ** 2
      t128 = t49 * t35
      t129 = t38 * t29
      t130 = t129 * t122
      t131 = t123 * t29 * t23
      t132 = t76 * t8
      t133 = t19 * t6
      t134 = t38 * t15
      t135 = t48 * t13
      t136 = t76 * t93
      t137 = t18 * t26
      t138 = t137 * t19
      t21 = t126 * (t138 * (t126 * t88 * t18 * t44 * t8 + t103 * t30 * t
     #34) + t119 * (t23 * (t127 * t49 * t29 + t51 * t103 + t52 * (-t102 
     #+ t103)) + t130 * t21 * t1 * t2 * t96 + t131 * t26 * t58 * (-t121 
     #* t35 + t98 * t92) * t57 * t63 + t132 * (t23 * (t67 * t48 + t73 * 
     #t91) + t29 * (t34 * (t104 * t15 + t13 * t91) + t128 * t68) * t18 *
     # t26 * t19) * t44 * t7 + t113 * t83 * t76) * t45 - t136 * (t109 * 
     #t50 * t63 * t58 * t57 * t29 * t23 + t133 * t22 * t28) * t75 + t122
     # * t85 * t64 * t19 * t60 * t44 * t62 * (t135 + t134) * t59)
      t21 = t21 + t126 * (t7 * (t110 * t36 * t23 * t56 * t37 * t97 * t88
     # * t57 + t109 * t8 * (t23 * (-t103 * (t101 + t100) + t102 * (t99 *
     # t91 + t47 * (t98 * t45 + t35) * t26)) + t107 * (-t48 * t49 * t45 
     #+ t105) * t19) * t44) + (t60 * (-t111 * t13 * t45 * t62 + t112 * t
     #34) + t99 * t36 * (-t30 * t8 + (-t114 + t113) * t29 * t1) * t34 * 
     #t26 + t4 * t93 * t22 ** 2 * t75 * t39 * t28 * t78) * t18 * t19 + t
     #123 * t34 * t29 * t23 * t60 * t58 * (-t121 * t15 * t34 * t59 + t35
     # * t92) * t57 * t63 + t125 * (t119 * ((-t115 * t59 - t30 * t45) * 
     #t34 * t13 * t36 + t30) - t118) * t96 + t102 * t119 * t6 * t23)
      t17 = 0.1q1 / t17
      t42 = 0.1q1 / t42
      t24 = 0.1q1 / t24
      t41 = 0.1q1 / t41
      t31 = 0.1q1 / t31
      t52 = 0.1q1 / t22
      t59 = 0.1q1 / t7
      t67 = t80 * t40
      t83 = t7 * t8 * t44
      t88 = t83 * t32
      t91 = t88 * t23
      t92 = t91 - t67
      t98 = t30 * t39
      t109 = t129 + t98
      t121 = t46 * t39
      t123 = t48 * t29
      t139 = t123 + t121
      t140 = t6 * t43
      t141 = t140 * t52 - t4
      t142 = mt ** 2
      t143 = t40 ** 2
      t144 = t33 ** 2
      t145 = t31 ** 2
      t146 = t31 * t145
      t147 = t43 ** 2
      t148 = t42 ** 2
      t149 = t40 * t48
      t150 = t87 * t32
      t151 = t23 * t37
      t152 = t44 * t8
      t153 = t122 * t1 * t2 * t96
      t154 = t153 * t40
      t155 = t122 * t47 * t1 * t2 * t96
      t156 = t129 * t142 * t17
      t157 = t139 * t97 * t57 * t58
      t158 = t157 * t7
      t159 = t151 * t25
      t160 = t8 * t43 * t52
      t161 = t141 * t39
      t162 = t91 * t4 * t39
      t163 = t19 * t40
      t164 = t129 * t1
      t165 = t164 * t2 * t96 * t30 * t18
      t166 = t100 * t7 * t58 * t44 * t23
      t167 = t76 * t25
      t168 = t167 * t24 * t52
      t169 = t101 * t18
      t170 = t169 * t24
      t171 = t109 * t59
      t172 = t123 * t151
      t173 = t31 * t36
      t174 = t32 ** 2
      t175 = t151 * t136
      t176 = t19 * t30
      t177 = t176 * t112
      t178 = t4 * t6
      t179 = t142 * t26 * t17
      t131 = t131 * t63
      t180 = t131 * t58
      t181 = t180 * t93 * t43 * t57
      t182 = t59 * t23
      t183 = t26 * t46
      t184 = t123 * t26
      t185 = t97 * t26
      t186 = t25 * t43
      t187 = t186 * t112
      t188 = t132 * t25 * t24 * t52 * (t100 * t112 * t63 * t57 * t48 * t
     #23 * t18 - t176 * t123 * t7 * t44 * t18 - t176) * t42 + t187 * t7 
     #* t44 * t32 * t23 * t24 * t18 * t52 * (-t103 + t102) * t148
      t189 = t173 * (t185 * t19 * (t87 * t81 * t40 * t59 * t31 + t152 * 
     #t139) + t151 * (t124 * t47 * t63 * t58 * t57 * t27 * t174 * t41 * 
     #t145 + (t39 * (-t25 * t30 + t183) + t184) * t59 * t17 * t142))
      t168 = t144 * t188 + t33 * (t42 * (t7 * (-t177 * t122 * t93 * t43 
     #* t44 * t24 * t52 + t175 * t44 * t24 * t58) + t181 * t46 * t24 * t
     #52 + t182 * (-t179 * t37 * t139 + t178)) + t170 * t43 * (-t102 * t
     #141 * t23 * t32 * t44 * t7 - t163 * t4) * t148) + t144 * (t168 * (
     #t165 + t166) * t42 + t170 * (t162 + t163 * (t161 * t1 - t160)) * t
     #148) + t33 * (t101 * t133 * t18 * t147 * t40 * t24 * t148 * t52 + 
     #t97 * (t172 * t63 * t58 * t93 * t57 * t39 * t24 + t171 * t26 * t96
     # * t2 * t1) * t42) + t173 * (t26 * (t31 * (t154 * t49 * t59 + t149
     # * t151 * t47 * t142 * t17 * t39 * t41 + t152 * (t39 * (t149 * t85
     # * t41 - t32 * t49) - t150) * t97 * t19) + t155 * t143 * t27 * t41
     # * t145) + t159 * (-t156 * t59 + t158)) + t168 * t112 * t92 * t148
     # * t33 * t144 + t189
      t170 = 0.1q1 / t49
      t188 = t102 * t7 * t44
      t189 = t188 * t48 * t23
      t190 = t19 * t38
      t191 = -t190 + t189
      t192 = t25 * t75
      t193 = t192 * t33
      t194 = t6 ** 2
      t195 = t75 ** 2
      t196 = t78 ** 2
      t197 = t41 ** 2
      t198 = t41 * t197
      t199 = t28 ** 2
      t200 = t28 * t199
      t201 = t7 * t46
      t202 = t201 * t29
      t203 = t159 * t7
      t204 = t83 * t25
      t205 = t81 * (-t193 - t26)
      t206 = t33 * t39
      t207 = t206 * t1
      t208 = t78 * t39
      t209 = t208 * t75
      t210 = t209 * t29
      t211 = t1 * t2 * t96
      t212 = t211 * t26
      t213 = t206 * t7
      t214 = t100 * t142
      t215 = t22 * t4
      t216 = t93 * t23
      t55 = t170 * (-t86 * t76 * t47 * t60 * t27 * t41 * t31 - t155 * t6
     #1 * t27 * t14 * t197 * t31) + t18 * (t213 * t58 * t44 * t95 * t29 
     #* t23 * t75 * t28 * t78 + t206 * (t176 * (-t74 - t8) + t166) * t78
     # * t28 * t75 * t93 + t111 * t60 * t35 * t34 * t45 * t191) + t23 * 
     #(-t214 * t61 * t17 * t39 * t29 * t35 * t34 - t126 * t9 * t127 * t5
     #9 * t29 - t102 * t118 * t4 * t78) + t97 * (t22 * (-t207 * t2 * t96
     # * t94 * t55 * t29 * t75 * t200 * t196 + t210 * t13 * t95 * t19 * 
     #(-t204 + t205) * t199 - t212 * t95 * t39 * t29 * t195 * t28) + t12
     #6 * (-t119 * t99 * t47 * t1 * t2 * t96 * t38 * t30 + t203 * (t192 
     #* t87 + (-t202 * t34 + t128) * t39 * t26) * t57 * t58 - t190 * t12
     #6 * t111 * t44)) + t216 * t75 * t28 * (t215 * t208 + t102) * t6 - 
     #t182 * t126 * t11 * t194 * t39
      t86 = t4 * t46
      t119 = t6 * t48
      t128 = t86 - t119
      t155 = t22 * t48 * t78
      t217 = t155 + t46
      t218 = t170 ** 2
      t219 = t29 * t33
      t220 = t219 * t78
      t221 = t6 * t37
      t222 = t33 * t48
      t223 = t26 * t22
      t224 = t7 * t19
      t225 = t25 * t27
      t226 = t151 * t97
      t227 = t226 * t7
      t228 = t211 * t171
      t229 = t228 * t97
      t230 = t19 * t26
      t231 = t19 * t14
      t232 = t231 * t122
      t233 = t98 * t151 * t142
      t234 = t27 * t6
      t235 = t76 * t19
      t236 = t235 * t60
      t237 = t83 * t27 * t170
      t238 = t8 * t33
      t180 = t180 * t57
      t79 = -t180 * t94 * t33 * t79 * t22 * t200 * t196 + t93 * (t127 * 
     #t22 * t29 * t23 * t78 + t103 * (t238 * t23 * t78 + (t164 * t33 * t
     #78 + t30) * t18 * t19)) * t28
      t164 = t219 * t153 * t94 * t22 * t75 * t195 * t28
      t239 = t224 * t112 * t122 * t44 * t94 * t22 * t195 * t28
      t240 = t137 * t29
      t77 = t240 * (t173 * t80 * t4 * t26 * t39 * t41 - t151 * t77 * t26
     # * t48 * t45 * (t90 + t89) * t57 * t58 * t63 + t237 * t41 * t31 * 
     #(t236 * t36 * t16 * t41 - t234 * t159))
      t55 = t75 * t79 + t58 * (t227 * (t210 * t7 * t94 * t15 * t22 * t19
     #9 + t70 * t93 * t39 * t75 + (t29 * (t225 + t120) + t90 * t39) * t5
     #6 * t60) * t57 + t203 * t108 * t44 * t60 * t29 * t35) + t55 - t232
     # * t66 * t7 * t44 * t61 * t27 * t16 * t198 * t218 * t31 - t233 * t
     #60 * t54 * t5 * t17 * t29 * t197 * t170 * t31 + t137 * t112 * t27 
     #* (t230 * t36 * t170 + t203 * t4 * t44) * t31 * t41 + t132 * (t75 
     #* (t224 * t95 * t29 * t18 * t78 * (t222 * t12 + t223 * t15) * t199
     # + t216 * t7 * (t220 * t128 - t221 * t217) * t28) + t220 * t7 * t1
     #8 * t19 * t94 * t15 * t22 * t195 * t199 + t18 * t35 * t60 * t19 * 
     #(-t105 + (t115 * t25 + t120 * t30) * t45 * t34 * t29 * t7)) * t44 
     #+ (t98 * t159 * t142 * t29 * t34 * t17 - t229) * t35 * t60 + t133 
     #* t118 * t76 * t38 * t78 - t164 - t239 + t77
      t9 = t11 * t194 * t39 + t9 * t127 * t29
      t11 = -t114 + t113
      t77 = t24 ** 2
      t79 = t24 * t77
      t105 = t52 ** 2
      t113 = t25 * t39
      t115 = t52 * t33
      t120 = t112 * t18
      t127 = t43 * t29
      t164 = t151 * t142
      t210 = t83 * t19
      t216 = t210 * t97
      t224 = t164 * t17
      t239 = t182 * t9
      t241 = t101 * t39
      t242 = t154 * t87 * t60
      t243 = t106 * t41
      t244 = t243 * t170
      t245 = t222 + t186
      t246 = t238 * t52 + t4
      t247 = t16 ** 2
      t248 = t219 * t52
      t249 = t43 * t37
      t250 = t33 * t46
      t251 = t26 * t43
      t252 = t33 * t25
      t253 = t180 * t95
      t82 = t252 * (t23 * (-t194 * t43 * t39 * t52 + (t248 + t37) * t8 *
     # t4) - t129 * t246 * t18 * t19 + t120 * t52 * t23 * (t86 * (t206 +
     # t249) - t103 * t245) * t44 * t7) * t24 + t253 * t144 * t43 * t247
     # * t79 * t105 + t219 * t122 * t52 * t93 * t19 * (t82 * t43 * t14 +
     # t83 * (-t250 * t10 - t251 * t16)) * t77
      t254 = t201 * t129
      t255 = t87 * t25
      t256 = t44 * t19
      t257 = t256 * ((t39 * (-t254 * t54 * t170 * t41 + t89) + t255) * t
     #31 * t26 - t85 * t60 * t27 * t32 * t39 * t41 * t145 + t213 * t123 
     #* t30 * t93 * t75 * t28 * t78) * t97 - t244 * t176 * t54 * t31
      t82 = t42 * t82 + t58 * (t203 * t31 * t57 * (t39 * (t85 * t26 * t4
     #8 * t41 + t89) + t255) * t97 + t243 * t203 * t173 * t44) + t8 * t2
     #57 + t36 * (t31 * (t241 * (t100 * t63 * t58 * t57 * t27 * t23 * t1
     #70 * t97 - t164 * t38 * t17 - t237 * t176 * t97) * t41 * t26 + t12
     #3 * t39 * (t224 + t216) * t41 * t60 - t239) + t242 * t41 * t145) +
     # t42 * (t127 * t153 * t95 * t144 * t14 ** 2 * t79 * t105 + t25 * (
     #t181 + t115 * (-t140 * t23 * (t73 + t102) + t29 * (t207 * t11 - t1
     #14 * t43) * t18 * t19) + t120 * (t186 * t106 * t19 + (t37 * (t160 
     #* t25 + t4 * t48) + t113 * t4) * t23 * t33 + t102 * t144 * t48 * t
     #23 * t52) * t44 * t7) * t24) + t206 * t172 * t142 * t30 * t93 * t1
     #7 * t75 * t28 * t78 + t244 * (t123 * t100 * t63 * t58 * t57 * t23 
     #* t18 + t165 + t166) * t31 * t54 + t233 * t95 * t17 * t29 * t75 * 
     #t28
      t89 = t22 * t39
      t140 = t89 + t127
      t160 = t30 * t43 * t52
      t164 = t160 + t38
      t165 = t26 * t32
      t166 = t93 * t144 * t16 * t42 * t77
      t172 = t102 * t23
      t181 = t226 * t58
      t186 = t224 * t48
      t207 = t211 * t97
      t226 = -t14 * t46 - t16 * t30
      t243 = t39 * t41
      t244 = t152 * t19
      t255 = t244 * t25 * t42
      t257 = t163 * t65
      t88 = t97 * (t36 * (t150 * t159 * t63 * t58 * t57 * t39 * t41 * t1
     #45 + t228 * t31) + t54 * (t177 * t7 * t44 * t48 * t170 * t39 * t41
     # * t31 + t243 * t29 * t19 * (-t88 * t38 + (t81 * t30 - t83 * t46) 
     #* t170 * t27 * t40) * t145) + t255 * t140) - t257 * t76 * t41 * t1
     #45
      t150 = t39 * t52
      t177 = t150 * t42
      t228 = t177 * t29
      t10 = t232 * t112 * t7 * t44 * t95 * t43 * t16 * t79 * t42 * t105 
     #+ t228 * t93 * (t181 * t63 * t57 * t48 * t16 + t214 * t10 * t23 * 
     #t17 + t216 * (t226 * t52 * t43 - t16 * t38)) * t77
      t79 = t230 * t39
      t216 = t8 * t25
      t232 = t216 * t18
      t10 = t144 * t10 + t26 * t88 + t33 * (t253 * t43 * t16 * t77 * t42
     # * t52 + t232 * t42 * (-t79 + t127 * (-t106 * t19 * t46 + t119 * t
     #151) * t52 * t44 * t7) * t24) + t7 * (t181 * t93 * t42 * t140 * t5
     #7 + t132 * ((t26 * (-t30 * t27 * t32 * t170 * t41 * t145 * t54 - t
     #225 * t40 * t41 * t145 * t36) - t115 * t93 * t14 * t42 * t77 * t24
     #5) * t18 * t29 * t19 + t172 * (t165 * t54 * t41 * t145 + t166 * t5
     #2)) * t44) + t173 * (t178 * t182 + t47 * (t207 * t38 * t40 * t31 +
     # t186 * t30 * t170) * t41 * t39 * t26) + (t101 * ((t207 * t164 - t
     #186) * t39 * t26 - t133 * t38 * t43 * t18 * t52) * t24 - t239) * t
     #42 * t33 + t180 * (t165 * (t46 * t27 * t170 + t48) * t145 * t41 * 
     #t54 + t166 * t43 * t46 * t105) + t231 * t136 * t52 * t77 * (t81 * 
     #t29 * t18 * t164 - t8) * t42 * t144 - t224 * t171 * t93 * t75
      t88 = t188 * t16 * t23 - t231
      t119 = t150 * t7
      t133 = t18 * t29
      t136 = t153 * t30
      t150 = t103 * t1
      t164 = t112 * t122 * t7 * t60 * t44
      t166 = t107 * t8
      t171 = t27 ** 2
      t181 = t256 * t33
      t141 = t8 * (t93 * (t42 * (-t151 * t73 * t7 * t43 * t44 * t24 + t2
     #48 * ((-t249 * t141 - t161 * t33) * t23 * t16 * t44 * t7 - t115 * 
     #t231 * t43) * t77) + t181 * t195 * t18 * t139) - t115 * t230 * t10
     #1 * t43 * t42 * t24) + t211 * t18 * t60 * t59 * t42 * t140 - t203 
     #* t126 * t57 * t18 * t139 * t58
      t182 = t188 * t32 * t23
      t91 = (-t91 * t107 * t6 * t170 * t145 * t41 - t164 * t163 * t16 * 
     #t170 * t145 * t197) * t27 * t54
      t186 = t224 * t59
      t203 = t186 * (t129 * t126 * t25 - t139 * t35 * t60 + t89 * t95 * 
     #t195)
      t91 = t18 * t141 + t36 * (t188 * t137 * t73 * t27 * t32 * t23 * t4
     #1 * t145 + t240 * t6 * t170 * (-t182 + t163) * t145 * t41 * t171) 
     #+ t25 * (t42 * (t133 * (t188 * t6 * t25 * t147 * t23 * t52 + t161 
     #* t80 * t33 * t26) * t24 + t115 * t169 * (t19 * (t6 * t14 * t52 * 
     #t147 + (t115 * t103 * t1 - t4) * t14 * t43) + t119 * t58 * t44 * t
     #144 * t16 * t23) * t77) + (t252 * (-t8 * t191 * t78 * t18 + t136 *
     # t38 * t78) * t28 + t179 * t151 * t48 * t59) * t75 * t29 + t233 * 
     #t126 * t17 * t59) + t54 * (t164 * t27 * t32 * t170 * t88 * t145 * 
     #t197 + t240 * (t162 + t163 * (t150 - t8) * t170 * t27) * t145 * t4
     #1) + t166 * t170 * t92 * t145 * t41 * t36 * t54 + (-t138 * t73 * t
     #40 * t27 * t145 + t87 * t137 * t170 * (t211 * t106 * t30 + t189 * 
     #t6) * t31) * t41 * t36 + t91 + t203
      t67 = t25 * (t7 * (t252 * t120 * (t172 * t33 * t43 * t16 * t42 * t
     #77 * t105 - t190 * t76 * t46 * t75 * t28 * t78) * t44 + t159 * t15
     #7 * t75) - t167 * t80 * t144 * t14 * t29 * t77 * t42 * t52 * t246 
     #+ t25 * (t127 * t159 * t142 * t17 - t207 * t109 * t33) * t59 * t19
     #5 + t180 * t193 * t46 * t48 * t28 * t78 + (t23 * (t179 * t100 * t3
     #9 - t178) - t211 * t185 * t109) * t59 * t75) + t54 * (t242 * t14 *
     # t197 * t170 * t145 + t240 * (t151 * t7 * t58 * t44 * t27 * t32 * 
     #t170 - t67 * t4 * t39) * t145 * t41) + t91 + t173 * t138 * (-t234 
     #* t129 * t170 + (t133 * (-t237 * t46 + t81 * t38) - t8) * t39 * t2
     #6) * t41
      t91 = t173 * t41
      t92 = t23 * t26
      t105 = t93 * t75 * t28
      t141 = t60 * t43
      t157 = t92 * (t91 * t170 * (-t194 * t27 * t39 + t111 * t4) + t214 
     #* t29 * t39 * (t126 * t68 * t34 * t45 - t105) * t17)
      t117 = t207 * t29 * (t113 * (t117 * t38 * t33 * t13 * t199 * t196 
     #+ t141 * t42 * t24) + t60 * (t98 * t54 * t14 * t218 * t197 + t59) 
     #* t31 * t27)
      t91 = t120 * (t60 * (t235 * t27 * t218 * t226 * t31 * t197 * t54 -
     # t231 * t173 * t167 * t27 * t170 * t197) - t250 * t235 * t95 * t19
     #5 * t28 + t92 * t91 * (t4 * (t37 * t48 + t113) + (t37 * (t216 + t8
     #6) - t103 * t25) * t170 * t27)) * t44 * t7
      t161 = t118 * t26
      t162 = t47 * t27
      t164 = t59 * t35
      t92 = t92 * t178 * (t162 * t41 * t170 * t31 + t164)
      t104 = (-t114 * t208 * t118 + (t25 * (-t6 * t42 * t24 * t52 * t147
     # + t4 * t42 * t24 * t43) - t104 * t4 * t27 * t170 * t41 * t31) * t
     #29 * t26) * t18 * t19
      t114 = t132 * (t23 * (t128 * t31 * t41 * t170 * t29 * t54 * t26 + 
     #t100 * t118 * t4 * t78) + t133 * (t95 * ((-t250 * t13 * t78 - t33 
     #* t22 * (t135 + t134) * t196) * t199 * t75 - t155 * t33 * t28 * t1
     #95) + t94 * t33 * t15 * t13 * t22 * t75 * t200 * t196 + t68 * t60 
     #* (-t173 * t16 * t197 * t170 + t46 * t35 * t34 * t45) - t161 * t48
     # * t78) * t19) * t44 * t7
      t110 = t131 * t110 * t57 * t54 * t48 * t16 * t197 * t170 * t31
      t49 = t153 * (t93 * (t26 * t30 * t29 * t28 * t75 - t223 * t59 * t1
     #95) + t60 * t49 * t31 * t59 + t219 * t30 * t95 * t195 * t28)
      t9 = t192 * (t199 * (t180 * t222 * t116 * t15 * t196 + t220 * t93 
     #* (t136 * t13 + t151 * (t142 * t48 * t12 * t17 * t39 - t7 * t58 * 
     #t44 * t15 * t18))) + t59 * (t23 * t9 + t211 * (-t33 * t140 * t97 *
     # t195 * t93 - t185 * t101 * t43 * t75)) + t130 * t22 * t25 * t19 *
     # (t204 - t205) * t78 * t28) + t240 * (t7 * (t151 * t54 * t44 * t48
     # * t170 * t58 + t244 * t225 * t106) + t175 * t63 * t58 * t57 * t27
     # + t80 * t54 * t170 * t39 * t11) * t31 * t41 + t173 * t76 * t170 *
     # t60 * (t51 * t88 + t133 * (t211 * t68 * t14 + t151 * (t90 * t170 
     #+ t25) * t16 * t27 * t57 * t58 * t63)) * t197 + t157 - t253 * t75 
     #* t28 * t217 + t133 * (t26 * (-t65 * t38 * t170 * t41 * t31 - t68 
     #* t4 * t41 * t31) + t238 * t95 * t13 * t75 * t199 * t78) * t19 + t
     #91 + t117 + t104 - t92 + t49 + t110 + t114
      t11 = t18 * t19
      t9 = t95 * (t120 * t28 * t44 * t7 * (t76 * (t188 * t46 * t23 - t17
     #6) * t78 * t28 * t15 * t33 + t172) * t75 + t11 * t8 * (-t18 * t140
     # * t44 + t219 * t28) * t195) + t9 - t89 * t224 * t94 * t29 * t195 
     #* t28 + t105 * t138 * t112 + t51 * t31 * t41 * t26 * (-t231 * t124
     # * t85 * t44 * t48 * t41 * t170 + t151 * (-t234 * t170 + t4))
      t12 = t123 + t121
      t49 = t129 + t98
      t68 = t151 * t32
      t76 = t113 * t24
      t85 = -t87 - t84
      t88 = t40 * t46
      t90 = t54 * t40
      t91 = t90 * t145
      t38 = t18 * (t33 * (-t255 * t18 * t49 + t18 * (t79 * t152 * t22 * 
     #t32 + t127 * (-t68 * t63 * t58 * t93 * t57 * t39 * t24 + t163 * t8
     #3 * t93 * t39 * t24 - t212 * t40 * t59)) * t148) + t76 * (t8 * (-t
     #182 + t163) + t133 * (t210 * (t160 * t32 + t149) - t68 * (t43 * t4
     #6 * t52 + t48) * t57 * t58 * t63 - t211 * t40 * (t160 + t38))) * t
     #148 * t144 + t65 * t145 * t18 * t44 * (-t163 * t12 + t188 * (-t87 
     #- t84) * t31 * t23 * t174)) + t207 * t59 * t40 * (t90 * t85 * t146
     # - t223 * t206 * t148) + t181 * t112 * t43 * t148 * (t88 * t213 * 
     #t25 * t24 * t52 + t165) * t97 + t186 * (-t252 * t109 * t42 - t91 *
     # t139)
      t48 = (t238 + t215) * t78
      t79 = t48 - t6
      t47 = t31 * t47 * t197 * t170 * t60
      t84 = t87 * t26 * t41 * t31
      t74 = t11 * (t27 * (t73 * t60 * t31 * t41 + t47 * t14 * (t150 * t3
     #6 * t170 - t4)) + t39 * ((t74 * t33 * t78 - t6) * t28 * t195 * t22
     # * t95 + t161 * t8 * t78 - t71 * t60 * t30 * t34 * t35) - t6 * t60
     # * t171 * t29 * t170 * t41 * t31)
      t5 = t83 * t18 * (t60 * (t86 * t99 * t50 * t23 * t35 * t34 + (t172
     # * t27 * t16 * t197 * t31 * t218 + t39 * (t4 * t16 * t23 - t176 * 
     #t5 * t18) * t31 * t197 * t170) * t29 * t54 + t162 * t151 * t4 * t1
     #6 * t197 * t170 * t31) + (-t103 * t23 * t29 * t28 + t208 * t23 * t
     #15 * (t6 * (t22 * t37 + t219) + t22 * (t37 * (-t238 - t215) - t73 
     #* t33) * t78) * t199) * t75 * t95)
      t2 = t18 * (t95 * (t209 * t28 * (t28 * (t219 * t79 * t1 + t22 * t7
     #9) * t13 * t19 + t83 * t23 * (t29 * (-t8 * t144 * t15 * t28 * t78 
     #+ t215) + t102 * t22)) + t28 * t39 * t19 * (t22 * t48 + t219 * (t2
     #38 * t78 - t6) * t1) * t195) + t47 * (t231 * (t27 * t170 * (t234 -
     # t51) + t36 * (t170 * (t138 * t2 * t27 * t14 * t41 - t51) - t4) * 
     #t39 * t1) + t83 * t170 * (-t103 * t36 * t27 - t221 * t171 + t65 * 
     #t39) * t23 * t16) - t164 * t141 * t211 * t169 * t75 + t105 * t80 *
     # t26 * t29 * t39 * (t215 * t78 - t6)) + t153 * t60 * (-t164 * t192
     # * t22 + t84) + t74 + t5 - t186 * t61 * t56 * t85 + t180 * t61 * t
     #54 * t27 * t247 * t198 * t218 * t31
      t4 = t129 * t7
      t5 = t32 * t25
      t8 = t244 * t97
      t14 = t25 * (t125 * t96 * t40 * t43 * t148 * t24 * t33 + t177 * t1
     #00 * t156 * t23 * t24 * t144) + t93 * (t127 * t154 * t14 * t148 * 
     #t77 * t52 * t144 + t206 * t151 * t156 * t42 * t24) + t227 * (-t5 *
     # t85 * t145 * t36 + t32 * t139 * t145 * t54 + t252 * t42 * (t39 * 
     #(t115 * t101 * t7 * t43 * t16 * t32 * t42 * t77 + t46) + t123)) * 
     #t57 * t58 + t91 * t229 + t8 * (t144 * (t76 * t254 * t52 * t42 - t1
     #19 * t127 * t93 * (t14 * t32 + t16 * t40) * t77 * t148) + t33 * (-
     #t251 * t241 * t7 * t32 * t24 * t148 + t4 * t93 * t39 * t24 * t42) 
     #+ t145 * t36 * (t25 * t40 * t85 + (t39 * (t84 * t7 * t40 - t30) - 
     #t129) * t32 * t36))
      t15 = t152 * t137 * (t113 * t66 * t7 * t23 * t41 * t170 * t31 + t1
     #1 * (-t105 * t202 * t39 + t126 * t12))
      t16 = t7 * t30
      t22 = t232 * t19 * (t228 * t1 * t26 * t144 * t24 + (t75 * (t39 * (
     #t16 * t93 * t29 * t28 - t183) - t184) - t126 * t49) * t18 * t44)
      t27 = t148 * t144
      t17 = t27 * (t97 * t19 * (t241 * t81 * t43 * t143 * t24 * t42 + t1
     #52 * t49 * t32) + t68 * (t101 * t122 * t63 * t58 * t43 * t57 * t32
     # * t42 * t24 + t49 * t59 * t17 * t142))
      t28 = t89 + t127
      t30 = t97 * t42 * t148 * t144
      t28 = t30 * (t7 * (t187 * t163 * t44 * t32 * t39 * t24 + t151 * t2
     #8 * t174 * t57 * t58) + t211 * t143 * t59 * t28)
      t36 = t163 * t140
      t4 = -t257 * t44 * t97 * t32 * t146 * t85 + t252 * t152 * t97 * (t
     #182 * t140 - t36) * t148 + t27 * (t159 * t156 * t32 * t39 * t24 + 
     #t8 * (t39 * (t4 * t5 * t24 - t88) - t149 * t29) + t229 * t40 + t68
     # * t158)
      t5 = t166 * t19 * (t201 * t108 * t44 * t60 * t35 - t105 * t1 * t33
     # * t78)
      t8 = t8 * (t12 * t42 * t33 * t26 + t173 * t25 * t49)
      ret = -2048 * t30 * t36 * t152 * t32 - 512 * t17 - 6 * t22 + 4 * t
     #67 - 3 * t15 - t53 / 16 + (0.3q1 / 0.4q1) * t5 + t204 * t107 * (-(
     #0.3q1 / 0.8q1) * t51 * t126 * t23 * t34 * t45 - 24 * t11 * t222 * 
     #t24 * t42) - t21 * t20 - t39 * t35 * t34 * t60 * (t137 * t111 * t7
     #2 * (t69 * t7 * t44 * t18 * t13 * t62 * t45 + t1) + t151 * (t192 *
     # t64 * (t71 - t6) * t18 * t44 * t7 + t179 * t29 * (t3 * t34 * t45 
     #- t70))) / 8 + 48 * t8 + 64 * t14 + 128 * t38 - 32 * t168 - 8 * t8
     #2 - t55 - t2 / 2 + 256 * t4 + 1024 * t28 - (0.3q1 / 0.2q1) * t236 
     #* t112 * (t1 * t54 * t41 * t170 * t31 + t16 * t108 * t44 * t25 * t
     #35) - 2 * t9 + 12 * t152 * t18 * t25 * (-t113 * t112 * t7 * t144 *
     # t23 * t24 * t42 * t52 + (-t98 * t192 + t129 * (t243 * t173 * t7 *
     # t26 - t192)) * t18 * t19) + 16 * t10

      hjetmass_triangle_pppm_s23_mhsq_s14 = ret/32q0/(0,1q0)
      return

      end function
