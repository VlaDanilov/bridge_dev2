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
 

      complex*32 function hjetmass_triangle_pppm_s34_mhsq_s12
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
      t11 = za(i3, i4)
      t12 = zb(i4, i3)
      t13 = t3 * t4
      t14 = t5 * t6
      t15 = t7 * t8
      t16 = t9 * t10
      t17 = t14 + t15 + t16 + t13
      t17 = -4 * t1 * t11 * t2 * t12 + t17 ** 2
      t17 = cg * sqrt(t17) + t13 + t14 + t15 + t16
      t18 = 2 * t11 * t12
      t19 = t18 + t17
      t18 = t17 ** 2 / 2 - t18 * t1 * t2
      t20 = (0.1q1 / 0.4q1)
      t21 = t17 ** 2
      t22 = t1 * t11
      t23 = t22 * t2 * t12
      t24 = t20 * t21
      t25 = -t23 + t24
      t25 = 0.1q1 / t25
      t26 = t15 + t13
      t27 = t17 * t1
      t28 = t27 * t2 * t25
      t29 = t24 * t1 * t2 * t25
      t30 = -t28 * t26 / 2 + t29
      t31 = t17 * t20
      t32 = -t31 * t4 + t9 * t2 * t12 / 2
      t33 = t17 * t11 * t25
      t34 = t33 * t32
      t13 = t14 + t13
      t35 = t33 * t12
      t24 = t24 * t11 * t12 * t25
      t36 = -t35 * t13 / 2 + t24
      t37 = t7 * t4
      t38 = t6 * t9 + t37
      t39 = t35 * t38
      t40 = t4 * t5 + t8 * t9
      t21 = t21 * t25
      t41 = t21 * t40
      t42 = t31 * t6 + t7 * t2 * t12 / 2
      t27 = t27 * t25
      t43 = -t27 * t42
      t44 = t33 * (t31 * t8 + t5 * t2 * t12 / 2)
      t22 = t22 / 2
      t45 = t22 * t4 - t31 * t9
      t46 = t17 * t2 * t25
      t47 = t46 * t45
      t23 = t23 * t17 * t25 / 2
      t26 = t20 * t21 * t26 - t23
      t40 = t28 * t40
      t5 = -t46 * (t22 * t8 + t31 * t5)
      t8 = t22 * t6 + t31 * t7
      t25 = t17 * t12 * t25
      t48 = t25 * t8
      t14 = t14 + t16
      t49 = t20 * t21 * t14 - t23
      t50 = t10 * t7 + t3 * t6
      t51 = t21 * t50
      t42 = -t33 * t42
      t50 = t28 * t50
      t15 = t15 + t16
      t16 = -t35 * t15 / 2 + t24
      t24 = -t31 * t10 + t3 * t2 * t12 / 2
      t35 = t27 * t24
      t3 = t22 * t10 - t31 * t3
      t10 = t25 * t3
      t22 = 2 * t1 * t2 + t17
      t31 = t21 * t38
      t8 = t46 * t8
      t14 = -t28 * t14 / 2 + t29
      t15 = t20 * t21 * t15 - t23
      t13 = t20 * t21 * t13 - t23
      t21 = -t25 * t45
      t23 = 0.1q1 / t5
      t25 = 0.1q1 / t43
      t28 = 0.1q1 / t40
      t29 = 0.1q1 / t9
      t38 = 0.1q1 / t7
      t40 = 0.1q1 / t51
      t45 = 0.1q1 / t12
      t17 = 0.1q1 / t17
      t52 = 0.1q1 / t11
      t41 = 0.1q1 / t41
      t53 = 0.1q1 / t50
      t54 = 0.1q1 / t44
      t55 = 0.1q1 / t48
      t18 = 0.1q1 / t18
      t56 = t6 * t38
      t57 = t2 * t52
      t58 = t56 + t57
      t59 = t34 * t48
      t60 = t39 * t26
      t61 = -t60 + t59
      t62 = t2 * t48
      t63 = t6 * t26
      t64 = t63 + t62
      t65 = t6 * t47
      t66 = t4 * t8
      t67 = t66 - t65
      t68 = mt ** 2
      t69 = t23 ** 2
      t70 = t23 * t69
      t71 = t2 ** 2
      t72 = t47 ** 2
      t73 = t72 ** 2
      t74 = t47 * t72
      t75 = t13 ** 2
      t76 = t34 ** 2
      t77 = t76 ** 2
      t78 = t34 * t76
      t79 = t22 ** 2
      t80 = t26 ** 2
      t81 = t54 ** 2
      t82 = t54 * t81
      t83 = t19 ** 2
      t84 = t18 ** 2
      t85 = t1 ** 2
      t86 = t30 ** 2
      t87 = t22 * t13
      t88 = t56 * t1 * t19
      t89 = t11 * t38
      t90 = t89 * t12
      t91 = t90 * t79
      t92 = t91 * t31
      t93 = t87 * t2
      t94 = t93 * t19
      t95 = t34 * t54
      t96 = t57 * t1
      t97 = t96 * t19
      t98 = t22 * t47
      t99 = t1 * t19
      t100 = t99 * t34
      t101 = t100 * t45
      t102 = t6 * t34
      t103 = t2 * t43
      t104 = t103 * t53 * t41
      t105 = t38 * t29
      t106 = t105 * t40
      t107 = t40 * t29
      t108 = t107 * t45 * t52
      t109 = t57 * t45
      t110 = t98 * t29
      t111 = t110 * t18
      t112 = t6 * t80 * t55
      t113 = t39 * t47
      t114 = t34 * t31
      t115 = t68 * t38 * t17
      t116 = t85 * t71 * t83
      t117 = t116 * t34
      t118 = t117 * t39
      t119 = t118 * t38 * t84
      t120 = t89 * t48
      t121 = t91 * t84
      t122 = t29 * t72
      t123 = t18 * t74
      t124 = t18 * t54
      t125 = t124 * t29 * t76
      t126 = t43 * t38
      t127 = t30 * t52
      t128 = -t126 - t127
      t129 = t6 ** 2
      t130 = t36 ** 2
      t131 = t56 * t39
      t132 = t4 * t42
      t133 = t132 * t52
      t134 = t47 * t23
      t135 = t134 * t87
      t136 = t99 * t45
      t137 = t1 * t39
      t138 = t137 * t28
      t139 = t28 * t40
      t140 = t129 * t38
      t141 = t116 * t105
      t142 = t141 * t84
      t143 = t89 * t67
      t144 = t6 * t30
      t145 = t105 * t1
      t146 = t145 * t72 * t31 * t28 * t23
      t147 = t8 * t38
      t148 = t113 * t29
      t149 = t11 * t43
      t150 = t30 * t25
      t151 = t68 * t31 * t17
      t152 = t151 * t52
      t153 = t152 * t45
      t154 = t23 * t13
      t155 = t107 * t48
      t156 = t123 * t22
      t157 = t105 * t54 * t41
      t158 = t19 * t2
      t159 = t158 * t22
      t160 = t159 * t84
      t161 = t160 * (t145 * t74 * t13 * t28 * t69 + t147 * t39 * t54 + t
     #148 * t54 - t146) * t34
      t162 = t54 * t41
      t163 = t162 * (t29 * (-t92 * t47 * t84 + t53 * (t22 * (t143 * t43 
     #+ t66 * t30) + t136 * (t132 * t128 - t131 * t30) * t2) * t18) + t1
     #44 * t45 * t53 * (-t56 + t57)) * t76
      t164 = t72 * t23
      t131 = t164 * (-t140 * t45 * t28 * t40 + (t139 * (-t136 * (t131 + 
     #t133) + t135) * t18 - t135 * t19 * t55 * (t138 * t38 + t23 * t36) 
     #* t84) * t29 * t2 + t142 * t72 * t130 * t52 * t45 * t28 * t69 * t5
     #5) * t26
      t67 = t121 * t73 * t75 * t26 * t29 * t28 * t70 * t55 + t122 * (t12
     #0 * t22 * t40 * t18 * t67 + (t40 * t61 * t18 * t19 * t71 * t1 - t1
     #19 + t115 * (-t114 + t113)) * t52 * t45) * t28 * t23 + t123 * (t60
     # * t85 * t71 * t83 * t18 * t36 * t38 * t29 * t52 * t45 * t55 + t2 
     #* (t38 * (-t109 + t111) + t6 * (-t106 * t45 * t26 - t108 * t55 * t
     #80)) * t36 * t19 * t1 + t87 * t107 * (t89 * t64 + t112)) * t28 * t
     #69 + t125 * (t30 * (-t94 * t34 * t36 * t25 * t18 * t81 + t95 * t41
     # * (t92 * t13 * t25 * t18 + (-t88 * t36 * t45 + t87) * t53 * t2)) 
     #+ t102 * t25 * t41 * t53 * (-t97 * t36 * t45 + t87) * t54 * t86 + 
     #t104 * (t101 * t58 - t98)) + t156 * t2 * t28 * t23 * (t154 * t38 -
     # t155) + t161 + t131 + t157 * (t153 + (-t150 * t99 * t36 * t31 * t
     #84 + t149 * t13 * t53 * t18) * t54 * t22 * t2) * t78 + t163
      t80 = 0.1q1 / t1
      t10 = 0.1q1 / t10
      t35 = 0.1q1 / t35
      t131 = t47 * t29
      t135 = t147 + t131
      t161 = t26 * t52
      t163 = t48 * t38
      t165 = t163 + t161
      t166 = t43 ** 2
      t167 = t43 * t166
      t168 = t10 ** 2
      t169 = t10 * t168
      t170 = t48 ** 2
      t171 = t48 * t170
      t172 = t35 ** 2
      t173 = t35 * t172
      t174 = t22 * t8
      t175 = t174 * t26 * t18
      t176 = t4 * t29
      t177 = t136 * t18
      t178 = t177 * t2
      t179 = t99 * t71
      t180 = t22 * t2
      t181 = t180 * t31
      t182 = t36 * t45
      t183 = t30 * t166 * t173
      t184 = t26 * t170 * t169
      t185 = t123 * t29
      t64 = t67 + t160 * t15 * t16 * t29 * (t184 + t183) + t18 * t41 * t
     #29 * (t1 * (t182 * t53 * t128 * t19 * t71 + t159 * t39 * t5 * t38 
     #* t18) + t144 * t87 * t89 * t53) * t81 * t78 + t164 * (t40 * (t63 
     #* t109 + t176 * (-t45 * t64 + t175)) + t178 * (-t57 * t39 + t155 *
     # (t102 - t132)) * t38) * t28 + t95 * t68 * t39 * t45 * t17 * t52 *
     # t80 * t135 + t185 * (-t117 * t36 * t38 * t45 * t52 * t18 - t136 *
     # t40 * t165 * t36 * t71 - t91 * t47 * t13 * t18) * t28 * t69 + t16
     #2 * (t181 * t38 * t18 + (t30 * (-t179 * t39 * t52 * t18 - t4 * t6)
     # - t103 * t4) * t53 * t45 * t29) * t76
      t67 = 0.1q1 / t26
      t87 = 0.1q1 / t30
      t91 = t97 * t16 * t45
      t155 = t22 * t15
      t186 = -t155 + t91
      t187 = t47 * t67
      t188 = -t10 * t15 - t187
      t189 = t134 * t55
      t190 = t189 - t10
      t191 = t39 ** 2
      t192 = t87 ** 2
      t193 = t31 ** 2
      t194 = t50 ** 2
      t195 = t34 * t67
      t196 = t16 * t10
      t197 = t155 * t6
      t198 = t4 ** 2 * t29
      t199 = t51 * t10
      t200 = t47 * t87
      t201 = t200 * t199
      t202 = t150 * t95
      t203 = t80 * t29
      t204 = t203 * t12 * t11
      t205 = t176 * t11
      t206 = t6 * t45
      t207 = t206 * t97
      t208 = t195 * t18
      t209 = t208 * t38
      t210 = t67 ** 2
      t211 = t99 * t39
      t212 = t48 * t67
      t213 = t170 * t16
      t179 = t179 * t18
      t214 = t97 * t45
      t215 = t22 * t87 * t210
      t216 = t204 * t79
      t217 = t45 * (t212 * t57 * t4 * t87 * t10 + t212 * t56 * t4 * t87 
     #* t10 + t179 * t38 * (t87 * (-t213 * t67 * t168 - t39 * (t211 * t2
     #9 * t18 + t212) * t10) + t19 * t191 * t23 * t55 * t18) * t52) * t5
     #1 * t47
      t218 = t35 * t18
      t188 = t170 * (t187 * t199 * t18 * t38 * (-t180 * t188 + t176 * (t
     #188 * t22 * t11 + t196 * t136 * t2)) * t87 - t199 * t177 * t71 * t
     #72 * t38 * t52 * t192 * t67) + t48 * (-t199 * t109 * t88 * t18 * t
     #72 * t192 + t201 * ((t197 * t10 + t136 * (t176 * t39 * t67 + (-t19
     #6 + t195) * t52 * t6) * t2) * t18 * t38 + t198 * t67 * t45)) + t84
     # * (-t199 * t1 * t71 * t83 * t191 * t38 * t52 * t45 + t158 * t98 *
     # t31 * t51 * t39 * t38 * t23 * t55 + t204 * (t190 * t26 + t202) * 
     #t193 * t79) + t209 * ((t2 - t205) * t87 * t31 * t22 * t43 + t207 *
     # t39) * t35 * t50 + t209 * t4 * t43 * t87 * t186 * t172 * t194 + t
     #199 * t178 * t105 * t4 * t72 * t170 * t192 * t67 + t218 * (-t216 *
     # t18 * t193 * t30 + t126 * (t214 * t192 * t67 - t215) * t194 * t76
     # * t4) + t217
      t204 = t114 * t22
      t217 = t204 * t18
      t219 = t34 * t29
      t201 = t201 * t48
      t220 = t176 * t34
      t221 = t22 * t39
      t222 = t43 * t35
      t223 = t176 * t89
      t224 = t223 * t50
      t225 = t26 * t31
      t226 = t225 * t28
      t227 = t226 * t29
      t228 = t200 * t29
      t229 = t164 * t55
      t230 = t199 * t115
      t231 = t158 * t18
      t232 = t231 * t38
      t233 = t71 * t1
      t234 = t233 * t83
      t235 = t234 * t84 * t191
      t236 = t136 * t71
      t237 = t236 * t43
      t238 = t155 * t103
      t57 = t35 * (t76 * (t180 * t166 * t38 * t18 * t87 * t210 + t178 * 
     #t166 * t38 * t192 * (t176 - t57) * t67) - t235 * t38 * t52 * t45 +
     # t195 * (-t142 * t191 * t45 * t52 + (-t98 * t56 * t18 + t45 * t58 
     #* t4 + t198 * t45) * t87 * t43)) + t209 * t43 * (-t237 * t16 * t52
     # * t87 + t238 * t87 + t197) * t172
      t58 = t48 * t10
      t197 = t58 * t72
      t239 = t200 * t10
      t240 = t36 * t72 * t69
      t241 = t240 * t55
      t242 = t38 * t51
      t196 = t242 * (t22 * (-t58 * t187 * t66 * t87 * t18 + t158 * (t31 
     #* (t48 * t16 * t168 - t241) + (t197 * t31 * t192 + t239 * (t31 * (
     #t196 * t48 + t39) + t212 * t113)) * t29 * t1) * t84) + t153 * t148
     # * t87 * t10)
      t57 = t50 * t57 + t71 * (t95 * t1 * t83 * t84 * t50 * t191 * t38 *
     # t52 * t45 * t25 - t218 * t126 * t101 * t50 * t39 * t52 * t87 * t6
     #7) + t92 * (t229 * (t154 * t51 * t80 - t227) - t51 * (t228 + t80) 
     #* t168 * t15 * t48) * t84 + t230 * t31 * t39 * t80 * t52 * t45 + (
     #t34 * (-t224 * t15 * t67 * t87 * t172 * t166 + t222 * t89 * t6 * t
     #31 * t29 * t67) + t76 * (-t224 * t210 * t87 * t35 * t166 + t222 * 
     #t56 * t50 * t210) + t131 * t31 * (t120 * t6 * t87 * t10 - t189 * t
     #4 * t26 * t28)) * t18 * t22 + t232 * ((t201 * (t133 * t45 - t217 *
     # t29) + t220 * t166 * t50 * t16 * t45 * t87 * t172 + (-t102 * t16 
     #* t45 * t52 * t172 + t219 * t87 * (t4 * t39 * t45 + t217) * t35) *
     # t50 * t43) * t67 * t1 + t221 * t31 * t18 * (t95 * t50 * t25 + t19
     #9)) + t196
      t196 = t14 * t42 + t49 * t8
      t212 = t53 ** 2
      t217 = t41 ** 2
      t224 = t28 ** 2
      t243 = t28 * t224
      t244 = t40 ** 2
      t245 = t74 * t26
      t246 = t245 * t48 * t23 * t224 * t244
      t247 = t30 * t78
      t248 = t247 * t43 * t217 * t212 * t54
      t249 = t180 * t99
      t250 = t249 * t84
      t251 = t116 * t49
      t27 = -t42 * t27 * t32
      t32 = t11 * t12
      t252 = t32 * t79
      t253 = t252 * t14 ** 2
      t254 = t116 * t49 ** 2 * t52 * t45
      t255 = t105 * t84
      t256 = t100 * t29 * t18
      t257 = -1 - t256
      t258 = -t195 + t200
      t259 = t15 ** 2
      t260 = t16 ** 2
      t261 = t39 * t50
      t262 = t30 * t36
      t263 = t182 * t137
      t264 = t76 * t81
      t265 = t264 * t84
      t266 = t121 * t50
      t267 = t234 * t260 * t45 * t52
      t268 = t200 * t99 * t2
      t269 = t255 * t195 * t22
      t270 = t50 * t38
      t271 = t16 * t29
      t272 = t271 * t45
      t273 = t225 * t122
      t274 = t242 * t74
      t275 = t39 * t31 * t29
      t276 = t72 * t26
      t277 = t276 * t13
      t278 = t30 * t43
      t279 = t95 * t80
      t280 = t78 * t130
      t281 = t39 * t43
      t282 = t252 * t31
      t283 = t159 * (t36 * (-t274 * t13 * t55 * t70 + t273 * t55 * t69) 
     #- t275 * (t10 * t26 + t30 * t35) + t187 * t145 * t76 * t166 * t50 
     #* t192 * t35 - t261 * t126 * t15 * t172) + t270 * t234 * t45 * t52
     # * (t281 * t16 * t172 - t280 * t25 * t82) + t282 * (t29 * (t80 * (
     #t278 * t15 * t172 - t277 * t69 * t55) - t58 * t200 * t147) + t279 
     #* t135)
      t284 = t162 * t4 * t76
      t285 = t60 * t80
      t286 = t206 * t80 * t4 * (t222 - t95 + t58 - t134)
      t287 = t29 * (-t22 * (t245 * t4 * t13 * t28 * t69 * t55 + t195 * t
     #103 * t31 * t35) - t284 * (t214 * t150 * t39 + t89 * t22 * t31)) *
     # t18
      t119 = t84 * t283 + t166 * (-t121 * t76 * t50 * t29 * t210 * t258 
     #* t35 + t269 * t50 * (t32 * t155 * t195 + t268 * t16) * t172 + t27
     #0 * t84 * (t252 * t195 * t259 * t29 + t267) * t173) + t25 * (-t266
     # * t75 * t78 * t80 * t82 - t125 * t22 * t6 * t31 * t86 * t41 * t53
     # + t265 * (-t263 * t71 * t83 * t50 * t38 * t52 - t216 * t30 * t31 
     #* t13 + t159 * (t261 * t13 * t38 + t262 * t31 * t29))) + t43 * (t2
     #69 * (t211 * t2 * t50 * t258 - t32 * t174 * t31) * t35 + t272 * t5
     #2 * (t119 * t50 * t67 - t151 * t30 * t80) * t172) + t185 * t99 * t
     #2 * (t18 * t44 * t22 * t31 * t38 + t161 * t182 * t4 * t55) * t28 *
     # t69 + t225 * t216 * t84 * t15 * t48 * t168 + t164 * (-t112 * t107
     # * t22 * t31 * t18 + t109 * t257 * t4) * t28 + t153 * t29 * (t126 
     #* t195 * t42 * t35 - t285 * t10) + t286 + t287
      t24 = -t24 * t33 * t38 + t29 * t44
      t33 = t228 * t1
      t125 = t42 * t38
      t126 = t164 * t4
      t159 = t29 * t35
      t185 = t58 * t135
      t211 = t252 * t84
      t216 = t195 * t166
      t3 = -t3 * t38 * t46 + t29 * t5
      t46 = t261 * t255 * t238 * t100 * t172 * t67
      t238 = (-t250 * t105 * t76 * t16 * t210 * t172 + t121 * t259 * t80
     # * t173) * t50 * t166
      t46 = t18 * (t2 * (t136 * t105 * t4 * t74 * t36 * t28 * t69 + t242
     # * t213 * t155 * t18 * t19 * t169 - t58 * t110 * t31 * t87) + t282
     # * t134 * t18 * t135 * t80) + t76 * (-t160 * t39 * t3 * t81 - t203
     # * t153 * t54) + t119 + t231 * (t126 * t138 * t105 * t45 + (t48 * 
     #(-t271 * t225 * t168 + t31 * (t125 * (1 + t33) + t219) * t10) + t1
     #46 * t60 * t55 + t213 * t200 * t145 * t51 * t15 * t169) * t18 * t2
     #2) + ((t125 * (t222 - t95) + t159 * (-t30 * t39 + t34 * t43) + t24
     # * t69 * t72) * t52 * t17 * t45 * t68 + t211 * (-t222 * t135 - t18
     #5)) * t80 * t31 + t142 * t45 * t52 * (t216 * t50 * t260 * t173 + t
     #229 * t191 * t26 * t28) + t131 * (-t230 * t72 * t170 * t52 * t45 *
     # t192 * t67 + t126 * t28 * (-t154 * t89 + 1) * t18 * t22 - t266 * 
     #t195 * t166 * t15 * t87 * t172) - t46 + t238
      t110 = -t222 + t134
      t119 = t222 * t16 + t39
      t126 = t6 * t39
      t138 = t99 * t2
      t146 = t68 * t15 * t17
      t229 = t100 * t2
      t7 = t170 * (t98 * t105 * t15 * t18 * t87 * (t229 * t51 * t67 * t1
     #8 - t11 * t6) * t168 - t255 * t199 * t117 * t72 * t52 * t45 * t192
     # * t67) + t171 * (-t180 * t122 * t89 * t18 * t87 * t210 * t10 + t1
     #87 * t177 * t105 * t71 * t16 * t87 * t168) + t2 * (t221 * t19 * (t
     #242 * t72 * t13 * t69 * t55 - t195 * t145 * t50 * t31 * t35 - t225
     # * t189 * t29) * t84 + t222 * t209 * (-t159 * t155 * t11 * t166 * 
     #t87 + t136 * (-t133 * t50 * t87 + t126 * t29))) + t45 * (t80 * (t1
     #40 * t110 * t9 + t198 * t7 * t110) - t195 * t152 * t105 * t35 * t1
     #19 * t50) + t48 * (t239 * t38 * t45 * t52 * (t148 * t68 * t51 * t6
     #7 * t17 - t138 * t6 * t42 * t18) + t242 * t39 * t45 * t52 * (t228 
     #* t116 * t16 * t84 - t146 * t80) * t168) + t71 * (t101 * t105 * t1
     #8 * t167 * t16 * t87 * t172 * t67 - t263 * t242 * t83 * t84 * t72 
     #* t52 * t69 * t55)
      t101 = t90 * t22
      t110 = t159 * t195
      t133 = t218 * t90
      t140 = t252 * t259
      t117 = t255 * t187 * t117
      t115 = t115 * t31
      t152 = t177 * t105 * t71
      t198 = t22 * t18
      t209 = (t1 * (-t231 * t122 * t60 * t4 * t23 * t28 * t55 + t71 * (-
     #t274 * t130 * t70 * t55 + t30 * t191 * t29 * (-t95 * t25 + t35)) *
     # t84 * t83) + t222 * (-t115 * t50 * t16 * t80 * t35 + t102 * t2 * 
     #t67)) * t52 * t45
      t221 = t2 * t18
      t230 = t174 * t18
      t238 = t68 * t45 * t17 * t52
      t263 = t198 * t12
      t266 = t187 * t87
      t9 = t35 * (t34 * (t230 * (t4 * t50 * t87 + t6) * t67 * t43 + t221
     # * t87 * (t136 * t132 * t29 + t174) * t67 * t166) - t180 * t11 * t
     #76 * t167 * t29 * t210 * t18 * t87 - t261 * t31 * (t238 * t80 + t1
     #60)) + t155 * t149 * t18 * (-t102 * t43 * t29 * t67 + t263 * t50 *
     # t31 * t80) * t172 + t266 * t2 * (t256 * (t98 * t51 * t67 * t18 - 
     #t206) + t206 + t230) * t10 * t170 - t282 * t265 * t50 * t13 * t80 
     #* t25 + t279 * t9 * t129 * t45 + t45 * t58 * (t148 * t116 * t51 * 
     #t52 * t84 * t87 * t258 - t9 * t129 * t80)
      t129 = t150 * t124
      t7 = t38 * t9 + t170 * (-t117 * t51 * t16 * t52 * t45 * t87 * t168
     # + t242 * t84 * (t140 * t80 + t267) * t169 - t187 * t38 * t18 * t8
     #7 * (t66 * t11 * t22 * t29 + t236 * t42 * t52) * t10) + t171 * (-t
     #111 * t89 * t2 * t15 * t87 * t67 * t168 + t152 * t72 * t192 * t67 
     #* t10) + t48 * (t242 * t137 * t71 * t83 * t84 * t16 * t52 * t45 * 
     #t168 + t239 * t109 * t6) + t7 + t198 * (t166 * (-t159 * t89 * t6 *
     # t76 * t210 - t110 * t143 * t87) + t18 * (t76 * (t270 * t158 * t36
     # * t25 * t81 * t31 + t101 * t162 * t150 * t29 * t193) - t101 * t74
     # * t51 * t75 * t80 * t70 * t55 - t275 * t202 * t158) + t133 * t22 
     #* t31 * t50 * t43 * t76 * t29 * t210) + t209 - t161 * t235 * t29 *
     # t45 * t190 + t176 * (t129 * t22 * t31 * t76 * t41 + t37 * (t95 - 
     #t58) * t80 * t45)
      t9 = t2 * t38
      t37 = t136 * t2
      t7 = t170 * (t239 * t178 * t105 * (t132 * t67 + t200 * (t200 * t97
     # * t51 * t18 + t6)) + t255 * t200 * t51 * (t116 * t260 * t45 * t52
     # + t140) * t169 + t268 * t105 * t18 * (t200 * t186 * t18 * t51 + t
     #206 * t16) * t168) + t34 * ((t221 * t272 * t88 * t172 + t9 * t45 *
     # t87 * (-t97 * t42 * t18 + t6) * t35) * t67 * t166 + (t92 * t50 * 
     #t15 * t29 * t84 * t172 - t218 * t38 * (t282 * t228 * t50 * t18 + t
     #207 * t42)) * t67 * t43) + t48 * (-t146 * t113 * t105 * t51 * t52 
     #* t45 * t87 * t168 + t239 * t56 * t18 * (t37 * t39 * t29 + t174)) 
     #+ t76 * (t222 * t153 * t105 * t50 * t87 * t67 + t152 * t167 * t192
     # * t35 * t67) + t7 + t198 * t164 * t31 * (-t231 * t24 * t23 - t223
     # * t28)
      t24 = t134 * t36
      t66 = t222 * t39
      t88 = t139 * t164
      t97 = t231 * t29
      t101 = t17 * t29 * t68
      t111 = t222 * t195
      t132 = t198 * t32
      t136 = t198 * t29
      t140 = t47 * t168
      t143 = t47 * t13
      t146 = t170 * t15 * t168
      t148 = t147 * t1
      t152 = t162 * t30
      t153 = t39 * t29
      t158 = t153 * t150 * t36
      t159 = t78 * t36
      t167 = t122 * t60
      t130 = (t17 * t80 * t39 * t68 * (-t147 * t134 + t185) + t233 * (t2
     #08 * t29 * (-t239 * t170 + t66) * t19 + (t39 * (t278 * t271 * t172
     # - t125 * t58) + t76 * (t153 * t54 - t158 * t81) + t159 * t29 * t8
     #1 - t167 * t36 * t69 * t55 - t245 * t130 * t70 * t29 * t55) * t84 
     #* t83)) * t52 * t45
      t70 = t198 * (t29 * (-t4 * t13 * (t89 + t150) * t81 * t41 * t78 + 
     #t231 * (t148 * t72 * t170 * t192 * t10 - t281 * t30 * t15 * t172) 
     #+ t152 * t31 * t53 * (t89 * t6 + t2) * t76) + t132 * (t147 * (t80 
     #* (t264 * t13 - t146) - t216 * t15 * t29 * t172) + t203 * (t170 * 
     #(t26 * t169 * t259 - t140 * t15) + t143 * (-t277 * t70 * t55 + t26
     #4))))
      t153 = t125 + t219
      t169 = t147 + t131
      t171 = t36 * t47
      t176 = t213 * t168
      t178 = t31 * t28
      t32 = t166 * (t180 * t100 * t105 * t18 * t50 * t15 * t16 * t173 * 
     #t67 + t110 * t180 * (t200 * (-t125 * t99 * t18 + 1) - t195) + t2 *
     # (-t96 * t45 * t153 * t18 * t16 * t83 + t198 * t169 * t16 * t19 - 
     #t155 * t195 * t29) * t172) + t22 * (t47 * (t29 * (-t146 * t2 * t87
     # + t284) + t178 * t9 * t134) + t32 * t156 * t23 * t29 * (t154 * t8
     #0 - t178 * t38) + t97 * (t47 * (t47 * (-t34 * t13 * t69 - t171 * t
     #69) + t114 * t23 + t176) - t13 * t78 * t81)) + t234 * t18 * (t48 *
     # (-t219 * t39 * t10 + t271 * t60 * t168) - t213 * t153 * t168 - t2
     #80 * t150 * t29 * t82) * t52 * t45 + t129 * t141 * t76 * t191 * t5
     #2 * t45 * t41
      t96 = t15 * t42 + t16 * t8
      t100 = -t13 * t42 - t36 * t8
      t110 = t96 * t38
      t114 = (t142 * t42 * (t47 * (-t58 * t39 * t87 - t176 * t87) - t72 
     #* t170 * t192 * t10 - t111 * t119) + t151 * t134 * (t125 + t219) *
     # t80 - t284 * t2 + t58 * t228 * t179 * (t200 * t48 + t39)) * t52 *
     # t45
      t119 = t198 * (t170 * (-t266 * t229 * t105 * t18 * t8 * t10 + t231
     # * (t219 * t15 + t110) * t168) + t2 * (t273 * t139 * t23 + (t76 * 
     #(t81 * (t100 * t38 - t131 * t36) + t150 * t157 * t137 * t31) - t15
     #9 * t270 * t13 * t25 * t82 - t240 * t38 * (t227 * t1 * t47 * t55 +
     # t8)) * t18 * t19) - t133 * t174 * t166 * t76 * t29 * t210 + t218 
     #* t103 * t19 * t39 * t135)
      t129 = t95 * t36 + t39
      t3 = t238 * t39 * (t264 * t80 * t3 + t131 * (t58 * t147 * t87 - t1
     #34 * t80))
      t58 = t271 * t177 * t52 * (t200 * t170 * t168 + t216 * t172) * t71
      t133 = t211 * (t166 * (-t15 * t80 * t135 * t172 + t203 * t30 * t25
     #9 * t173) + t147 * t47 * (t143 * t69 * t80 - t146 * t29 * t87))
      t96 = t231 * (t284 * t145 * t45 * t129 + t198 * (t166 * (t172 * (t
     #125 * t15 + t219 * (t110 * t67 * t1 + t15)) + t145 * t76 * t42 * t
     #210 * t35 + t270 * t16 * t15 * t173) + t38 * t47 * (t96 * t168 * t
     #87 * t29 * t170 * t1 + t42 * t23 * (-t134 * t13 + t31))))
      t3 = t18 * t32 + (t1 * (t97 * (t126 * t162 * t86 * t76 * t25 * t53
     # + t88 * t112 * t39 + (t150 * t4 * t36 * t81 - t4 * t54) * t41 * t
     #78) + (t29 * (t134 * t34 * (t24 + t39) + t184 * t260 + t183 * t260
     # - t66 * t34) + t125 * (t134 * t39 + t264 * t36 + t95 * t39 + t240
     # - t66)) * t84 * t83 * t71) + t101 * (-t38 * (t78 * t166 * t50 * t
     #210 * t87 * t35 + t44 * t74 * t31 * t69 * t28) - t285 * t15 * t48 
     #* t168) + t117 * t42 * t170 * t87 * t10) * t52 * t45 + t136 * (t23
     #1 * (t13 * (t276 * t69 * t55 + t264 * t150) + t111 * t147 * t1) * 
     #t39 + t132 * (t226 * t74 * t13 * t69 * t38 * t55 + t266 * t147 * t
     #34 * t166 * t35 - t150 * t78 * t75 * t82 * t80) + t88 * t63 * t89 
     #* t31) + t130 + t70 + t114 + t119 + t96 + t58 + t133 + t3
      t32 = t144 + t103
      t44 = t198 * t38
      t55 = t36 * t14
      t58 = t53 * t54
      t66 = t58 * t217
      t70 = t73 * t23
      t96 = t78 * t43
      t97 = t22 * t14
      t110 = t29 * t18
      t111 = t63 + t62
      t94 = -t94 * t106 * t18 * t1 * t73 * t26 * t49 * t224 * t69 + (t29
     # * (-t236 * t170 * t49 * t38 + t97 * (t120 * t111 + t111 * t26)) *
     # t244 + t9 * (-t214 * t48 * t49 + t97 * (-t110 * t99 * t60 + t48))
     # * t40) * t224 * t23 * t74 - t127 * t237 * t49 * t78 * t29 * t217 
     #* t212 * t54
      t9 = t18 * t94 + t110 * (t251 * t123 * t52 * t45 * t40 * t38 * t22
     #4 * t23 * (t26 * (t24 + t39) - t59) + t138 * (t73 * (-t55 * t198 *
     # t26 * t38 * t224 * t69 * t40 + t163 * t198 * t49 * t224 * t23 * t
     #40) + t66 * t49 * t78 * (-t44 * t31 * t30 + (-t6 * t86 * t52 - t9 
     #* t166 - t278 * t56) * t53 * t45) + (t44 * (t26 * t8 * t21 + t59 *
     # t14) * t40 + t49 * t45 * t26 * (t52 * (-t63 - t62) - t56 * t48) *
     # t244) * t224 * t23 * t74) + t97 * (t247 * t217 * t212 * t54 * t32
     # + t89 * (t96 * t217 * t212 * t54 * t32 + t263 * (t66 * t247 * t31
     # + t70 * t224 * t40 * (t154 * t26 - t48)))))
      t24 = t43 * t47
      t32 = -t95 * t30 * t13 + t24
      t44 = t252 * t8
      t56 = t106 * t84
      t59 = t182 * t116
      t86 = t162 * t76
      t94 = t88 * t48
      t106 = t86 * t43 * t53 + t94
      t111 = t150 * t13
      t99 = t99 * t29
      t32 = t38 * (t2 * (t94 * t230 + t99 * t22 * ((-t262 * t8 * t53 * t
     #81 + t58 * t43 * t8) * t41 * t78 + t86 * t39 * (-t30 * t8 * t53 + 
     #t47) + t225 * t74 * t49 * t224 * t23 * t40 - t111 * t77 * t36 * t4
     #1 * t82) * t84) + t71 * (-t96 * t162 * t85 * t83 * t84 * t42 * t29
     # * t52 * t45 * t53 - t177 * t42 * t52 * t106) - t167 * t238 * t139
     # * t8 * t23) + t23 * (-t245 * t107 * t92 * t84 * t14 * t224 + t56 
     #* t72 * (-t116 * t42 * t45 * t52 * t61 - t252 * t47 * t48 * t8 + t
     #249 * t48 * (t34 * t8 + t42 * t47)) * t28) + t56 * t245 * (t59 * t
     #42 * t52 + t249 * t100 + t44 * t13) * t28 * t69 + t86 * (t53 * (t1
     #8 * (t147 * t103 * t22 + t220 * (t37 * t128 * t49 + t97 * (t89 * t
     #43 + t30)) * t41) + t255 * (t127 * t116 * t45 * t129 * t42 + t249 
     #* t32 * t42 - t44 * t32)) + t238 * t113 * t105)
      t44 = -t214 * t49 + t97
      t56 = t30 * t29
      t21 = t38 * ((t99 * t103 * t44 * t84 * t54 + t56 * (t252 * t14 * t
     #13 - t249 * (t13 * t49 + t55) + t251 * t182 * t52) * t84 * t81) * 
     #t53 * t217 * t77 + t66 * (t103 * t44 * t18 + t29 * (-t149 * t12 * 
     #t47 * t79 * t14 + t127 * t251 * t39 * t45 + t249 * (t30 * (-t14 * 
     #t39 + t27) + t24 * t49)) * t84) * t78 - t245 * t108 * t68 * t8 * t
     #21 * t17 * t224 * t23)
      t24 = t51 ** 2
      t44 = t22 * t31
      t55 = t18 * t38
      t6 = t55 * (t4 * (-t197 * t215 * t24 + t87 * (t197 * t214 * t24 * 
     #t87 + t140 * (-t155 + t91) * t24 * t48 + t204 * t194 * t35) * t67)
     # + t44 * t51 * (-t239 * t6 + t132 * (t10 * (-t228 - t80) + t189 * 
     #t80) * t31))
      t58 = t19 * t39
      t61 = t45 * t52
      t66 = t4 * t47
      t69 = t31 * t42
      t5 = t29 * (t86 * t115 * t127 * t42 * t45 * t53 + t250 * t38 * (t1
     #64 * t26 * t28 * (t241 * t13 + t69 * t40) + (-t111 * t39 + t171) *
     # t81 * t41 * t78) - t156 * t4 * t14 * t26 * t224 * t23 * t40 + t59
     # * t77 * t38 * t52 * t41 * t81 * t84 * (t150 * t36 * t54 - 1)) + t
     #41 * (t150 * t121 * t75 * t77 * t29 * t82 + t38 * t78 * (-t211 * t
     #131 * t13 + t61 * (-t101 * t5 * t39 + t158 * t116 * t84 - t179 * t
     #36) - t93 * t18 * t257) * t81 + t124 * t76 * (t1 * (-t56 * t109 * 
     #t102 * t19 * t53 - t58 * t38 * t45 * t52 * t71) + t56 * t22 * t53 
     #* (-t230 * t90 * t31 + t65) - t118 * t105 * t18 * t52 * t45)) + t1
     #39 * t122 * t18 * t23 * (t22 * (t63 * t47 + t89 * (-t66 * t48 * t1
     #4 * t28 - t175 * t12 * t31)) + t37 * (t66 * t165 * t28 * t49 - t16
     #1 * t102))
      t13 = t47 * t170 * t10
      t19 = t34 * t166 * t35
      t16 = t232 * (t198 * ((-t33 - 1) * t168 * t15 * t51 * t48 * t39 - 
     #t50 * t43 * t16 * t31 * t172 * (t195 * t1 * t29 + 1)) + t137 * t2 
     #* t29 * t45 * t87 * t67 * (t19 + t13))
      t2 = t55 * (t44 * t67 * (t102 * t50 * t35 + t201 * (t205 - t2)) + 
     #t214 * t39 * t87 * (t195 * t4 * t194 * t35 - t65 * t199) + t252 * 
     #t18 * t50 * (t35 * (t195 * t29 + t80) - t279 * t25) * t193)
      t15 = t180 * t18 * (t89 * t31 * t29 * t106 + t58 * (t48 * (t10 * (
     #t131 * (t148 * t87 + 1) + t147) - t26 * t15 * t29 * t168) - t134 *
     # t169) * t18)
      t1 = t160 * t31 * (t43 * (t35 * (t219 * (-t125 * t1 * t67 - 1) - t
     #125) + t271 * t30 * t172) + t95 * t153)
      t18 = t138 * t105 * t18 * (-t104 * t76 * t39 * t45 * t54 + t44 * t
     #124 * t78 * t41 + t164 * t39 * t28 * (-t62 * t45 * t40 + t98 * t18
     #))
      ret = t55 * (t138 * (-t61 * t113 * t4 * t24 * t87 * t67 * t10 / 8 
     #+ t136 * (-192 * t88 * t60 * t8 + t152 * t76 * t53 * (-8192 * t14 
     #* t43 * t49 * t76 * t217 * t53 + 96 * t69))) - t98 * t4 * t31 * t2
     #4 * t87 * t67 * t10 / 16) - 2048 * t70 * t255 * t26 * t48 * t243 *
     # t244 * (t254 + t253) - 512 * t21 - 48 * t18 - t20 * t6 - t2 / 2 +
     # 256 * t9 + 1024 * t105 * (-t27 * t127 * t68 * t78 * t17 * t45 * t
     #217 * t53 * t54 + t250 * (t246 * t196 + t248 * t196) + (-t248 - t2
     #46) * t84 * t14 * t8 * t79 * t12 * t11 - t251 * t42 * t45 * t52 * 
     #t84 * (t248 + t246)) + 4096 * t255 * (t249 * t73 * t14 * t26 * t49
     # * t48 * t243 * t23 * t244 + (t254 + t253) * t54 * t212 * t41 * t2
     #17 * t43 * t77 * t30) - t188 + 64 * t5 + 128 * t32 + 12 * t1 + 16 
     #* t3 + 24 * t15 + 2 * t57 + 3 * t110 * t181 * t89 * t87 * t67 * (t
     #19 + t13) + 6 * t16 - 32 * t64 - 8 * t46 + 4 * t7

      hjetmass_triangle_pppm_s34_mhsq_s12 = ret/32q0/(0,1q0)
      return

      end function
