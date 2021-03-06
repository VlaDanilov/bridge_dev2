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
 

      complex*32 function hjetmass_bubble_ppmm_s23
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          
          complex*32 ret
          double precision mt

      t1 = za(i3, i4)
      t2 = zb(i3, i1)
      t3 = zb(i4, i2)
      t4 = zb(i2, i1)
      t5 = zb(i4, i3)
      t6 = 0.1q1 / t5
      t7 = t2 * t3
      t8 = t7 * t6
      t9 = t8 - t4
      t10 = za(i1, i3)
      t11 = za(i1, i2)
      t12 = za(i2, i4)
      t13 = t3 ** 2
      t14 = t6 ** 2
      t15 = t6 * t14
      t16 = t1 * t2
      t17 = t4 * t12
      t18 = t1 * t11
      t19 = t18 * t4
      t20 = t11 * t4
      t21 = t10 * t2
      t22 = t18 * t2
      t23 = t3 * t6
      t24 = t1 * t10
      t25 = t11 * t2
      t26 = t25 * t12
      t27 = t3 * (-t26 * t3 * t13 * t15 + (t12 * (-t21 + t20) - t22) * t
     #14 * t13 + t23 * (t10 * (t17 - t16) + t19) + t24 * t4)
      t28 = za(i1, i4)
      t29 = zb(i4, i1)
      t30 = 0.1q1 / t12
      t31 = t1 * t30
      t32 = t31 + t23
      t33 = t12 * t5
      t34 = t33 + t25
      t35 = t10 * t4
      t36 = t1 * t3
      t37 = t35 + t36
      t38 = t12 * t3
      t39 = t1 * t5
      t40 = t20 - t21 + t38 - t39
      t40 = t40 ** 2
      t41 = 0.4q1 * t37 * t34 + t40
      t41 = sqrt(t41)
      t42 = t20 - t21 + t38 - t39 + t41
      t43 = 0.1q1 / t34
      t44 = 0.1q1 / 0.2q1
      t45 = t44 * t42 * t43
      t46 = t45 - t23
      t41 = t20 - t21 + t38 - t39 - t41
      t47 = t44 * t41 * t43
      t48 = t47 - t23
      t49 = -t47 * t2 + t4
      t40 = 0.4q1 * t37 * t34 + t40
      t40 = sqrt(t40)
      t50 = t20 - t21 + t38 - t39 - t40
      t50 = 0.1q1 / t50
      t51 = 0.2q1
      t52 = t51 * t37
      t53 = -t52 * t50 - t47
      t54 = t43 * (t42 - t41)
      t55 = -t47 * t12 - t1
      t56 = 0.1q1 / t34
      t57 = t28 * t29
      t58 = t12 * (-t39 - t57 + t20 + t38) + t22
      t59 = t41 ** 2
      t60 = t59 ** 2
      t61 = t41 * t59
      t62 = t56 ** 2
      t63 = t62 ** 2
      t64 = t56 * t62
      t65 = t62 * t59
      t66 = t65 * t11
      t67 = t11 * t12
      t68 = (-t39 - t21) * t56
      t69 = t10 * t12
      t70 = t4 * t11 ** 2
      t71 = t24 * (t57 + t21 + t39)
      t72 = (t1 * (t11 * (-t38 + t57) + t70) + t57 * t69) * t56
      t73 = 0.3q1 / 0.2q1
      t74 = 0.3q1 / 0.4q1
      t75 = t3 * t51
      t76 = t5 * t41 * t56
      t77 = 0.4q1 * t20
      t78 = t69 * t41
      t79 = t78 * t56
      t80 = 0.1q1 / t11
      t81 = t10 * t80
      t47 = -t47 - t81
      t70 = t1 * (t11 * (-t38 + t57) + t70)
      t82 = t39 + t21
      t83 = -t21 - t57
      t84 = t1 ** 2
      t85 = t84 * t5
      t86 = t85 * t43
      t87 = t12 * t41
      t88 = t61 * t64
      t89 = t66 * t12
      t90 = t84 * t3
      t91 = 0.1q1 / 0.4q1
      t92 = 0.1q1 / 0.8q1
      t93 = 0.1q1 / 0.16q2
      t94 = -t44 * t41 * t10 * (t83 * t43 * t1 - t86 + t87 * t62 * (-t25
     # * t41 * t56 + t38)) - t73 * t78 * t43 * t37 + t74 * t69 * t65 * t
     #82 + t91 * t65 * (t69 * (t76 * t12 + t57) + t70) - t92 * t88 * t11
     # * t58 + t93 * t67 * t60 * t63 * t34 - t10 * (t4 * (t89 + t24) + t
     #90)
      t95 = t4 * t5
      t96 = t7 + t95
      t97 = t12 ** 2
      t98 = t36 * t43 * (t20 + t57)
      t99 = t97 * t3
      t100 = t99 * t5
      t101 = t33 * t65
      t102 = t38 + t20
      t103 = t21 + t57
      t104 = t38 * t57
      t84 = t84 * t5 ** 2
      t105 = t67 * t96 * t56
      t106 = t12 * (t21 - t57) - t22
      t107 = t35 * t73 * t43 + t36 * t51 * t56
      t60 = -t44 * t41 * (t100 * t59 * t64 + t78 * t96 * t62 - t98) + t7
     #4 * t38 * t65 * t102 - t91 * t65 * (t39 * t103 + t105 * t41 - t104
     # + t84) + t92 * t88 * t5 * t106 + t93 * t33 * t60 * t63 * t34 + t3
     #6 * (-t101 + t35 + t36) + t87 * t3 * t107
      t78 = t23 * t12
      t87 = t78 + t1
      t108 = -t45 * t12 - t1
      t109 = -t45 * t2 + t4
      t110 = t42 ** 2
      t111 = t110 ** 2
      t112 = t62 * t110
      t22 = t12 * (t39 + t57 - t20 - t38) - t22
      t113 = t64 * t42 * t110
      t114 = t69 * t42
      t115 = t112 * t67
      t43 = t44 * t42 * t10 * (t103 * t43 * t1 + t86 + t12 * t42 * t62 *
     # (t25 * t42 * t56 - t38)) - t73 * t114 * t43 * t37 + t74 * t112 * 
     #t69 * t82 + t91 * t112 * (t69 * (t33 * t42 * t56 + t57) + t70) + t
     #92 * t113 * t11 * t22 + t93 * t67 * t111 * t63 * t34 - t10 * (t4 *
     # (t115 + t24) + t90)
      t40 = t20 - t21 + t38 - t39 + t40
      t70 = -t45 - t81
      t40 = 0.1q1 / t40
      t45 = -t52 * t40 - t45
      t52 = t5 * t42 * t56
      t82 = t114 * t56
      t86 = t112 * t33
      t103 = t69 * t80
      t116 = t103 - t1
      t117 = t21 * t80
      t118 = t117 + t4
      t119 = t10 ** 2
      t120 = t80 ** 2
      t121 = t80 * t120
      t122 = t81 * t5
      t123 = t1 * t12
      t90 = t51 * t123 * t119 * t80 * (t122 + t3) - t10 * (t10 * t119 * 
     #t97 * t5 * t121 + t99 * t119 * t120 + t81 * t85 + t90)
      t99 = 0.1q1 / t2
      t124 = t4 * t99 + t81
      t125 = t23 + t81
      t126 = t10 * (-t33 * t119 * t120 + t36 + t81 * (t39 - t38))
      t127 = t36 * t12
      t128 = (t38 + t20) * t56
      t7 = t7 + t95
      t95 = t36 * (t20 + t57)
      t129 = (t39 * (-t21 - t57) - t84 + t104) * t56
      t130 = t1 * t13
      t131 = t3 * (t67 * t13 * t14 + t24 + t23 * (t69 + t18))
      t32 = 0.1q1 / t32
      t132 = 0.1q1 / t48
      t48 = -0.1q1 / t48
      t133 = -0.1q1 / t46
      t53 = 0.1q1 / t53
      t28 = 0.1q1 / t28
      t54 = 0.1q1 / t54
      t45 = 0.1q1 / t45
      t29 = 0.1q1 / t29
      t46 = 0.1q1 / t46
      t134 = t17 + t16
      t135 = t49 * t53
      t16 = t16 * t30
      t136 = t9 * t32
      t28 = t28 * t29
      t29 = t28 * t56
      t124 = 0.1q1 / t124
      t137 = -0.1q1 / t70
      t138 = 0.1q1 / t125
      t139 = -0.1q1 / t47
      t125 = 0.1q1 / t125
      t140 = t124 ** 2
      t141 = t9 ** 2
      t142 = t87 * t32
      t87 = t4 * t87
      t143 = t1 * t9
      t144 = t6 * t125
      t145 = t144 * t4
      t146 = t36 * t14
      t147 = t125 ** 2
      t148 = t46 ** 2
      t149 = t132 ** 2
      t150 = t137 ** 2
      t151 = t138 ** 2
      t152 = t139 ** 2
      t153 = t116 * t90
      t154 = t153 * t62 * t150 * t152
      t155 = t126 * t14
      t156 = t155 * t147
      t157 = t27 * t32
      t158 = t141 * t149 * t62 * t148
      t159 = t131 * t120
      t8 = t28 * (t4 * (t10 * t118 * t121 * t99 * (-t156 + t154) * t140 
     #+ t120 * (t81 * (t156 - t154) + ((-t153 + t81 * (-t116 * (-t51 * t
     #10 * (t1 * (-t21 - t20) + t69 * (t117 + t4) + t119 * t97 * t5 * t1
     #20) - t10 * (t3 * (t81 * t97 - t123) + t85 + t57 * (t103 - t1)) + 
     #0.3q1 * t33 * t119 * t1 * t80) - t12 * t90)) * t152 * t150 * t62 +
     # t14 * t147 * (t81 * (t69 * (-t122 + t3) + t24 * t51 * t5) + t126)
     #) * t99 * t118) * t124) + t31 * t14 * t32 * (-t9 * (t9 * t27 * t62
     # * t148 * t149 + t159 * t151) + t23 * (t158 * (t157 + t51 * t3 * (
     #t1 * (t39 + t21) + t13 * (t26 * t14 + t6 * t97)) + t3 * (t69 * (t8
     # + t4) - t19 + t57 * (t78 + t1)) + t130 * (0.3q1 * t25 * t6 + 0.4D
     #1 * t12)) + (t131 * (t136 - t2) + t9 * (-t75 * t18 + t38 * (-t23 *
     # t11 + t10))) * t151 * t120)))
      t9 = 0.1q1 / t47
      t13 = 0.1q1 / t70
      t18 = t109 ** 2
      t19 = t13 ** 2
      t20 = t133 ** 2
      t21 = t48 ** 2
      t24 = t54 ** 2
      t26 = t9 ** 2
      t31 = t56 * t41
      t47 = t55 * t94 * t120 * t26
      t57 = t56 * t42
      t39 = (t44 * t42 * (-t100 * t110 * t64 - t114 * t96 * t62 + t98) +
     # t74 * t112 * t38 * t102 + t91 * t112 * (-t105 * t42 + t39 * t83 +
     # t104 - t84) + t92 * t113 * t5 * t106 + t93 * t33 * t111 * t63 * t
     #34 - t36 * (t86 - t35 - t36) + t38 * t42 * t107) * t14 * t20
      t63 = t135 * t50
      t70 = t28 * t24
      t31 = t70 * t62 * t37 * ((-t39 * t57 * t2 + (t57 * t13 * t19 + t19
     #) * t120 * t43 * t108) * t45 * t40 * t109 + t39 * (t57 * t133 + 0.
     #1q1) * t45 * t40 * t18 + t63 * (t47 * (t31 * t9 + 0.1q1) + (t21 * 
     #(-t2 * t41 * t56 + t49) + t31 * t49 * t48 * t21) * t14 * t60))
      t57 = t40 ** 2
      t78 = t45 ** 2
      t9 = t28 * (t62 * (t42 * (t109 * t80 * t13 * t57 * t78 * t108 ** 2
     # + t18 * t6 * t133 * t57 * t78 * t108) - t41 * t49 * t55 * t50 ** 
     #2 * t53 ** 2 * (t49 * t6 * t48 + t55 * t80 * t9)) * t54 * t37 ** 2
     # + t36 * (t136 * t15 * (t2 * t27 * t149 * t62 * t148 - t159 * t138
     # * t151) * t30 + t157 * t158 * t15 * (t132 + t46) * t30) + t35 * t
     #99 * t121 * t124 * t118 * (t155 * t125 * t147 - t154 * (t137 + t13
     #9)))
      t13 = t108 * t109
      t15 = t53 * t50 * t41 * (t49 ** 2 * t14 * t21 * (t60 * t53 + t44 *
     # t33 * t88 * t34 - t51 * t79 * t7 - t73 * t89 * t7 + t74 * t65 * t
     #5 * t106 + 0.4q1 * t127 * (-t76 + t3) + 0.3q1 * t38 * (t128 * t41 
     #- t101 + t35) + t129 * t41 + t95) + (t94 * (-t12 * t49 + t55 * (t1
     #35 - t2)) - t49 * t55 * (-t44 * t67 * t61 * t64 * t34 + t74 * t66 
     #* t58 + 0.3q1 * t69 * (-t65 * t25 + t68 * t41 + t35 + t36) - t72 *
     # t41 - t71 + t79 * (t12 * (-t76 * t73 + t75) + t77))) * t26 * t120
     #)
      t2 = t70 * t64 * t37 * (t15 + (t45 * (-t18 * (-t44 * t113 * t33 * 
     #t34 + t51 * t82 * t7 + t73 * t115 * t7 - t74 * t112 * t5 * t106 - 
     #0.4q1 * t127 * (-t52 + t3) - 0.3q1 * t38 * (t128 * t42 + t35 - t86
     #) - t129 * t42 - t95) * t14 * t20 + (t43 * (-t108 * t2 - t109 * t1
     #2) - t13 * (-t44 * t113 * t67 * t34 - t74 * t112 * t11 * t22 + 0.3
     #q1 * t69 * (-t112 * t25 + t68 * t42 + t35 + t36) - t72 * t42 - t71
     # + t82 * (t12 * (-t52 * t73 + t75) + t77))) * t19 * t120) + t109 *
     # (t108 * t43 * t120 * t19 + t39 * t109) * t78) * t40 * t42)
      t3 = t63 * t41
      t5 = t28 * t54 * t62 * t6 * t37 * (t42 * (t1 * t133 * t40 * t45 * 
     #t18 + t13 * t4 * t133 * t40 * t45) - t3 * t48 * (t1 * t49 + t4 * t
     #55))
      t3 = t28 * t54 * t24 * t64 * t37 * (t42 * (t13 * t43 * t120 * t19 
     #* t40 * t45 + t39 * t40 * t45 * t18) - t3 * (t49 * t60 * t14 * t21
     # + t47))
      ret = -0.20q2 * t29 * t6 * ((-t110 * t109 * t133 * t40 * t45 * t13
     #4 + t135 * t134 * t48 * t50 * t59) * t54 * t62 * t37 + t136 * t130
     # * t14 * t46 * t132 * (t16 + t4)) - 0.8q1 * t28 * (t10 * (t4 ** 2 
     #* t116 * t118 * t120 * t99 ** 2 * (t116 * t56 * t137 * t139 - t144
     #) * t140 + t145 * t120 * t99 * (t1 * t118 + t116 * t4) * t124) + t
     #145 * t121 * (t17 * t99 + t1) * t124 * t119 + t146 * t32 * (t142 *
     # t1 * t141 * t30 ** 2 * t56 * t46 * t132 + (t30 * (t143 * (-t142 *
     # t30 - 0.1q1) - t87) + t23 * (-t16 - t4)) * t138 * t80)) + 0.16q2 
     #* t8 - 0.12q2 * t146 * t29 * t136 * t30 * t46 * t132 * (t143 + t87
     #) + 0.128q3 * t31 - 0.32q2 * t9 + 0.64q2 * t2 - 0.24q2 * t5 - 0.25
     #6q3 * t3

      hjetmass_bubble_ppmm_s23 = ret/16q0*(0,1q0)
      return

      end function
