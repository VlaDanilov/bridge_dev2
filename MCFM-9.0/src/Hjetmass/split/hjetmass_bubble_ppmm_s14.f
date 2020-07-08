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
 

      complex*32 function hjetmass_bubble_ppmm_s14
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          
          complex*32 ret
          double precision mt

      t1 = za(i1, i3)
      t2 = zb(i2, i1)
      t3 = za(i2, i4)
      t4 = za(i3, i4)
      t5 = zb(i3, i1)
      t6 = t3 * t2
      t7 = t4 * t5
      t8 = t7 + t6
      t9 = zb(i3, i2)
      t10 = zb(i4, i3)
      t11 = za(i1, i2)
      t12 = zb(i4, i2)
      t13 = t11 * t12
      t14 = t1 * t10
      t15 = t14 + t13
      t16 = t11 * t2
      t17 = t1 * t5
      t18 = t3 * t12
      t19 = t4 * t10
      t20 = -t18 - t19 + t16 + t17
      t21 = 0.4q1
      t20 = t20 ** 2
      t22 = t21 * t8
      t23 = t22 * t15 + t20
      t23 = sqrt(t23)
      t24 = -t18 - t19 + t23 + t16 + t17
      t23 = t18 + t19 + t23 - t16 - t17
      t25 = 0.1q1 / t24
      t26 = 0.1q1 / t10
      t27 = 0.1q1 / t23
      t28 = -t25 + t27
      t29 = t8 ** 2
      t30 = t29 ** 2
      t31 = t8 * t30
      t32 = t8 * t29
      t33 = t5 * t32 * t26
      t34 = 0.8q1
      t35 = t5 ** 2
      t36 = t26 ** 2
      t37 = t21 * t35
      t38 = 0.16q2
      t39 = t25 ** 2
      t40 = t25 * t39
      t41 = t37 * t36
      t42 = t38 * t29 * t39
      t43 = t5 * t8 * t26 * t25
      t44 = 0.2q1
      t22 = t22 * t25
      t45 = 0.32q2
      t46 = t8 * t27
      t47 = 0.1q1 / t8
      t48 = 0.1q1 / t5
      t49 = 0.1q1 / t9
      t50 = t2 ** 2
      t51 = t12 ** 2
      t52 = t21 * t2
      t53 = t52 * t12
      t54 = t46 * t34
      t55 = za(i2, i3)
      t56 = t1 * t15
      t9 = t55 * t9
      t57 = (-t9 + t16 + t17) * t1
      t58 = -t44 * t1 * (t19 + t18) + t57
      t59 = t1 * t8
      t60 = t47 ** 2
      t61 = t60 ** 2
      t62 = t47 * t60
      t63 = t44 * t2
      t64 = t12 * t60
      t65 = t27 ** 2
      t66 = t27 * t65
      t67 = t5 * t29 * t26 * t39
      t68 = t8 * t25
      t69 = t68 * t38
      t70 = t65 * t29
      t71 = 0.1q1 / t2
      t72 = 0.1q1 / t3
      t73 = t11 * t72
      t74 = t12 * t71 + t73
      t75 = t10 * t48
      t76 = t73 + t75
      t77 = t11 ** 2
      t78 = t72 ** 2
      t79 = t1 * t12
      t80 = t4 * t2
      t81 = 0.1q1 / t4
      t82 = t1 * t81
      t83 = t82 + t75
      t84 = t10 ** 2
      t85 = t48 ** 2
      t86 = t73 * t4
      t87 = t1 - t86
      t88 = 0.1q1 / t8
      t89 = 0.1q1 / 0.2q1
      t90 = t89 * t24 * t88
      t91 = t90 - t73
      t92 = t89 * t23 * t88
      t93 = -t92 - t73
      t94 = t79 * t15
      t95 = t1 * t3
      t96 = t4 * t11
      t97 = t9 * t1
      t98 = t2 * t10
      t99 = t5 * t12
      t100 = t99 + t98
      t101 = t1 ** 2
      t102 = -t44 * t101 * t100 - t12 * (t12 * (t96 - t95) + t97) - t53 
     #* t11 * t1
      t103 = t29 * t39
      t104 = t3 * t51
      t105 = t9 * t12
      t106 = t4 ** 2
      t107 = (-t105 - t104) * t4
      t108 = t97 * t2
      t109 = t106 * t12
      t110 = t109 * t10
      t111 = -t110 + t108 + t107
      t112 = t1 * t4
      t113 = t1 * t2
      t114 = t113 * (t16 + t17)
      t115 = t34 * t2
      t116 = -t25 + t27
      t117 = t72 * t49
      t118 = t99 * t26
      t119 = t118 - t2
      t120 = t35 * t49
      t121 = t5 * t49
      t122 = t121 * t26
      t123 = t17 * t26
      t124 = 0.1q1 / t11
      t125 = t124 ** 2
      t126 = t124 * t125
      t127 = t3 * t125
      t128 = t73 * t2 + t12
      t129 = t19 * t48 + t1
      t130 = -t98 * t48 + t12
      t131 = t4 * t12
      t132 = t101 * t12
      t133 = -t99 - t98
      t98 = t99 + t98
      t99 = t1 * t98
      t134 = 0.3q1
      t135 = (t9 + t18) * t4
      t136 = t90 + t75
      t75 = -t92 + t75
      t137 = -t90 * t2 - t12
      t138 = -t7 - t6
      t139 = (-t9 + t16 + t17) * t88
      t140 = t24 * t1
      t141 = (t19 + t18) * t47
      t142 = t9 + t18 + t19
      t143 = t24 ** 2
      t144 = t24 * t143
      t145 = t4 * t143 * t60
      t146 = -0.1q1 / 0.4q1
      t147 = 0.1q1 / 0.8q1
      t148 = -t89 * t140 * (t138 * t60 * t24 + t139) + t146 * t145 * t14
     #2 - t147 * t4 * t144 * t62 * t8 + t1 * (t141 * t24 - t13 - t14)
      t20 = t21 * t8 * t15 + t20
      t20 = sqrt(t20)
      t149 = -t18 - t19 + t16 + t17 + t20
      t150 = t88 * (t24 + t23)
      t149 = 0.1q1 / t149
      t151 = t44 * t15
      t152 = -t151 * t149 - t90
      t90 = -t90 * t4 + t1
      t153 = t12 * (t14 + t13)
      t154 = t6 * t64
      t88 = t88 * t12 * (t12 * (-t96 + t95) - t97)
      t155 = t17 * t80
      t156 = t112 * t98 * t60
      t104 = t4 * (t105 + t104) - t108 + t110
      t105 = t143 * t60
      t157 = t9 + t16
      t144 = t80 * t144 * t62
      t158 = 0.3q1 / 0.4q1
      t159 = t13 * t44
      t160 = t113 * t47
      t161 = t160 * t24
      t140 = t146 * t105 * t104 - t147 * t144 * t157 - t89 * t24 * (t155
     # * t143 * t62 + t156 * t24 + t88) + t158 * t114 * t143 * t60 + t80
     # * t143 ** 2 * t61 * t8 / 0.16q2 + t1 * (t140 * t98 * t47 - t154 *
     # t143 + t153) + t161 * (t159 - 0.3q1 / 0.8q1 * t105 * t6)
      t143 = t135 * t47
      t162 = t106 * t10
      t7 = t7 + t6
      t163 = t79 * t2
      t96 = t12 * (t12 * (t96 - t95) + t97)
      t97 = t92 * t4 + t1
      t164 = t92 * t2 - t12
      t17 = -t18 - t19 + t16 + t17 - t20
      t17 = 0.1q1 / t17
      t20 = -t151 * t17 + t92
      t92 = t23 * t1
      t151 = t23 ** 2
      t165 = t23 * t151
      t166 = t4 * t151 * t60
      t167 = t151 * t60
      t168 = t80 * t165 * t62
      t160 = t160 * t23
      t6 = t146 * t167 * t104 + t147 * t168 * t157 + t89 * t23 * (t155 *
     # t151 * t62 - t156 * t23 + t88) + t158 * t114 * t151 * t60 + t1 * 
     #(-t92 * t98 * t47 - t154 * t151 + t153) + t80 * t151 ** 2 * t61 * 
     #t8 / 0.16q2 + t160 * (-t159 + 0.3q1 / 0.8q1 * t167 * t6)
      t55 = 0.1q1 / t55
      t61 = 0.1q1 / t76
      t76 = 0.1q1 / t76
      t88 = 0.1q1 / t93
      t104 = 0.1q1 / t75
      t151 = 0.1q1 / t91
      t83 = 0.1q1 / t83
      t74 = 0.1q1 / t74
      t153 = 0.1q1 / t136
      t154 = t130 ** 2
      t155 = t81 ** 2
      t156 = t3 ** 2
      t159 = t123 * t129
      t61 = t61 * t48
      t169 = t12 * t87
      t170 = t49 * t12
      t123 = t125 * (t34 * t79 * t46 * t72 * t49 * t25 * (t59 * t116 + t
     #4) - t22 * t101 * t2 * t72 * t49 * t27 + t117 * t26 * (t12 * (-t12
     #3 - t4) + t113)) + t170 * (t1 * t26 - t169 * t128 * t71 ** 2 * t74
     # ** 2 * (t87 * t47 * t88 * t151 + t61)) * t126
      t171 = t4 * t71
      t100 = t124 * (t72 * t60 * t49 * (t1 * t70 * (t21 * t103 * (t21 * 
     #(-t21 * t112 * t100 - t115 * t79 * t3 + t111 * t44 + 0.6q1 * t114)
     # + t68 * (-t102 * t45 + 0.96q2 * t68 * t94)) + t46 * (t45 * t103 *
     # (t102 * t21 - t69 * t94) + 0.384q3 * t32 * t94 * t39 * t27)) + t4
     #4 * t4 * (t38 * t30 * t102 * t39 * t65 + 0.64q2 * t31 * t94 * t39 
     #* t65 * t28)) + t72 * t48 * (t1 * t120 * t36 * t119 + t4 * t122 * 
     #t119)) + t61 * t74 * t49 * t12 * (-t171 * (t12 + t128) - t1 - t87)
     # * t125
      t114 = t121 * t36
      t30 = t55 * (t3 * t100 - t4 * (-t34 * t46 * t26 * t49 * t25 * t12 
     #* (t8 * t12 * (t25 - t27) + t2) + t22 * t5 * t51 * t49 * t36 * t27
     #) + t48 * (-t49 * (t42 * t5 * (t12 * (t12 * (t21 * t59 - t44 * (t9
     # + t18 + t19) * t4) - t52 * t58) - t44 * t50 * t56) * t26 * t65 + 
     #t51 * t56 * t70 * (t21 * t67 * (-t34 * t35 * t36 - 0.96q2 * t29 * 
     #t39 + t43 * t45) + t46 * (t45 * t67 * (-t21 * t5 * t26 + t69) - 0.
     #384q3 * t33 * t39 * t27)) * t60) + t44 * t64 * t49 * (t12 * t58 + 
     #t63 * t56) * (-0.64q2 * t5 * t31 * t26 * t39 * t65 * t28 - t38 * t
     #35 * t30 * t36 * t39 * t65)) + t123 * t156 + t114 * t83 * t1 * (-t
     #159 * t154 * t155 * t47 * t83 * t153 * t104 + (t130 * (-t159 * t83
     # * t155 + 0.1q1) + t12 + t81 * (t1 + t129) * t2) * t72 * t76))
      t31 = t113 + t131
      t9 = t55 * (t124 * (t156 * ((t38 * t117 * t70 * t39 * (t1 * t102 +
     # t4 * t94) + 0.64q2 * t117 * t1 * t32 * t94 * t39 * t65 * t28) * t
     #124 + t170 * t71 * t74 * (t87 * (t44 * t73 * t1 * (t86 * t133 + t9
     #9) + t132 * t10 + (-t110 + t9 * (t113 - t131)) * t78 * t77 - t9 * 
     #t80 * t11 * t77 * t72 * t78 + t9 * t79 * t73 + t77 ** 2 * t106 * t
     #2 * t5 * t78 ** 2 + t113 * t77 * t5 * t78 * (t1 * t134 - t86 * t21
     #)) * t60 * t88 ** 2 * t151 ** 2 + t61 * (-t73 * t63 * t1 + t80 * t
     #77 * t78 - t79)) * t125) - t49 * t26 * t31 * t124 * t3 + t114 * t3
     #1) + t82 * t120 * t26 * t36 * (-t154 * (-t44 * t14 * t3 * t48 * t1
     #30 - t79 * t11 + t14 * (-t9 + t16) * t48 - t135 * t85 * t84 + t80 
     #* t3 * t10 * t84 * t48 * t85) * t60 * t104 ** 2 * t153 ** 2 + (t14
     # * t63 * t48 + t80 * t84 * t85 - t79) * t76 * t72) * t83)
      t10 = t46 * t49 * t25
      t56 = t169 * t127 * t55 * t47 * t49 * t74 * t151 * t88
      t22 = t55 * (t1 * t49 * (t51 * t46 * (-t44 * t43 * (t34 * t5 * t8 
     #* t26 * t25 - t41 - t42) + t46 * (-t34 * t43 * (-t44 * t5 * t26 + 
     #t22) + t45 * t5 * t29 * t26 * t25 * t27)) * t48 * t47 - t53 * (t25
     # * t27 * t28 * t33 * t34 + t37 * t29 * t36 * t25 * t27) * t48 * t4
     #7 + t54 * t50 * t26 * t25) + t3 * t117 * (t112 * t21 * (t34 * t32 
     #* t12 * t25 * t27 * t28 - t52 * t29 * t25 * t27) * t47 + t101 * t4
     #6 * (t44 * t103 * (t69 * t12 + t115) + t46 * (-t34 * t68 * (t22 * 
     #t12 + t63) + t45 * t29 * t12 * t25 * t27)) * t47 + t109 * t54 * t2
     #5) * t124)
      t20 = 0.1q1 / t20
      t28 = -0.1q1 / t75
      t32 = -0.1q1 / t136
      t33 = 0.1q1 / t150
      t37 = -0.1q1 / t91
      t42 = -0.1q1 / t93
      t43 = 0.1q1 / t152
      t46 = t137 ** 2
      t51 = t43 ** 2
      t52 = t20 ** 2
      t53 = t149 ** 2
      t54 = t97 * t72
      t58 = t54 * t42
      t59 = t97 * t164
      t61 = t164 ** 2
      t13 = (t146 * t166 * t142 + t147 * t4 * t165 * t62 * t8 + t89 * t9
     #2 * (t7 * t60 * t23 + t139) - t1 * (t141 * t23 + t13 + t14)) * t28
      t14 = t13 * t61 * t48 + t58 * t6
      t60 = t33 ** 2
      t62 = t148 * t48
      t64 = t72 * t37
      t67 = t43 * t149
      t68 = t67 * t40
      t7 = t68 * (t46 * (-t62 * t32 ** 2 + (-t148 * t43 + t158 * t145 * 
     #t8 + t162 * t24 * t47 + t24 * t143 - t44 * t1 * (t7 * t47 * t24 + 
     #t18 + t19) + t57) * t48 * t32) + t64 * (t140 * (-t90 * (t37 + t43)
     # + t4) + t90 * (-t134 * t161 * (t5 * (-t24 * t4 * t47 + t1) + t16)
     # + t158 * t105 * t80 * t157 - t21 * t163 * (-t3 * t24 * t47 + t11)
     # - t44 * t1 * (t133 * t47 * t24 * t4 + t99) - t89 * t144 * t8 - t1
     #11 * t47 * t24 - t96 + 0.9q1 / 0.4q1 * t95 * t105 * t50)))
      t24 = t14 * t66 * t17
      t43 = t49 * t60 * t55 * t15
      t69 = t43 * t8
      t5 = t69 * (t24 * t52 + (t61 * t28 * t48 * (t13 - t158 * t166 * t8
     # + t162 * t23 * t47 + t143 * t23 + t44 * t1 * (t138 * t47 * t23 + 
     #t18 + t19) - t57) + (-t4 * t6 - t97 * (t134 * t160 * (t5 * (t23 * 
     #t4 * t47 + t1) + t16) + t158 * t167 * t80 * t157 - t21 * t163 * (t
     #3 * t23 * t47 + t11) - t44 * t1 * (t98 * t47 * t23 * t4 + t99) + t
     #89 * t168 * t8 - t96 - (t110 - t108 - t107) * t47 * t23 + 0.9q1 / 
     #0.4q1 * t167 * t95 * t50)) * t72 * t42 + t54 * t6 * t42 ** 2) * t6
     #6 * t17 * t20 + t7)
      t6 = t67 * t90 * t39 * t37
      t7 = t117 * t33 * t55 * t15 * t8
      t11 = t64 * t90 * t140 + t62 * t46 * t32
      t13 = t69 * t2 * t48 * (t13 * t164 * t66 * t17 * t20 - t68 * t137 
     #* t148 * t32)
      t16 = t27 * t49 * t25 * t124 * t55 * t29 * t1 * (t156 * t94 * t125
     # * t25 * t27 + t113 * t116 + t131 * t116)
      t3 = -t132 * t10 * t45 * t156 * t126 + t38 * t127 * (t48 * t122 * 
     #t1 * (t118 - t63) - t79 * t3 * t124 * t49 * t26) - t34 * (t44 * t1
     #22 * (t79 * t35 * t36 - t80) - t113 * t41 * t49) * t124 * t48
      ret = -t21 * t22 + t30 * t34 - t38 * t9 + t55 * t3 - 0.128q3 * t10
     # * t124 * t55 * (t95 * t31 * t124 - t31 * t4) - 0.12q2 * t56 * (t1
     #71 * t128 + t87) - 0.512q3 * t49 * t33 * t55 * t15 ** 2 * t29 * (t
     #90 * t46 * t48 * t40 * t53 * t32 * t51 + t59 * t66 * t17 ** 2 * t5
     #2 * (t164 * t48 * t28 + t58) + t90 ** 2 * t137 * t72 * t40 * t53 *
     # t37 * t51) + 0.1024q4 * t5 + 0.192q3 * t7 * (t2 * t97 ** 2 * t65 
     #* t17 * t42 * t20 + t6 * (t137 * t4 - t90 * t2) - t59 * t4 * t65 *
     # t17 * t42 * t20) + 0.4096q4 * t49 * t33 * t60 * t55 * t15 * t8 * 
     #(t68 * t11 + t24 * t20) + 0.6144q4 * t43 * t29 * (t65 ** 2 * t17 *
     # t20 * t14 + t67 * t11 * t39 ** 2) + 0.320q3 * t7 * (t31 * t65 * t
     #42 * t17 * t20 * t97 - t6 * t31) - 0.2048q4 * t13 + 0.256q3 * t16 
     #- 0.20q2 * t56 * (t171 * t12 + t1)

      hjetmass_bubble_ppmm_s14 = ret/16q0*(0,1q0)
      return

      end function
