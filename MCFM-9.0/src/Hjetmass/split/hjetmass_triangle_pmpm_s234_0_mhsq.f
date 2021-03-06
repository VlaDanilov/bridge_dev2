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
 

      complex*32 function hjetmass_triangle_pmpm_s234_0_mhsq
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

      t1 = za(i2, i4)
      t2 = zb(i3, i1)
      t3 = za(i2, i3)
      t4 = zb(i3, i2)
      t5 = zb(i4, i2)
      t6 = za(i3, i4)
      t7 = zb(i4, i3)
      t8 = t3 * t4
      t9 = t1 * t5
      t10 = t6 * t7
      t11 = t9 + t10 + t8
      t12 = za(i1, i2)
      t13 = zb(i2, i1)
      t14 = za(i1, i3)
      t15 = za(i1, i4)
      t16 = zb(i4, i1)
      t17 = t12 * t13
      t18 = t14 * t2
      t19 = t15 * t16
      t20 = t19 + t17 + t18
      if ( real(t20) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t20 = cg * sqrt(t20 ** 2) + t17 + t18 + t19
      t21 = 0.1q1 / t20
      t22 = -2 * t12
      t23 = t22 * t2 * t11 * t21 + t1 * t7
      t24 = t22 * t16 * t11 * t21 - t3 * t7
      t22 = t22 * t13 * t11 * t21 + t8 + t9
      t25 = t14 * t23
      t26 = t12 * t22
      t27 = t13 * (t25 + t26)
      t28 = t1 * t16
      t29 = t2 * t3 + t28
      t30 = -2 * t18 * t11 * t21 + t10 + t8
      t31 = t15 * t24
      t32 = t13 * (t25 + t31)
      t33 = -2 * t15
      t34 = t33 * t2 * t11 * t21 + t1 * t4
      t33 = t33 * t13 * t11 * t21 - t4 * t6
      t35 = t2 * (-2 * t14 * t13 * t11 * t21 + t5 * t6)
      t36 = t12 * (t22 * t13 + t35)
      t37 = t16 * t33
      t35 = t12 * (t35 + t37)
      t38 = 0.1q1 / t12
      t39 = 0.1q1 / t11
      t40 = 0.1q1 / t3
      t32 = 0.1q1 / t32
      t41 = 0.1q1 / t4
      t42 = 0.1q1 / t22
      t27 = 0.1q1 / t27
      t43 = 0.1q1 / t15
      t20 = 0.1q1 / t20
      t44 = 0.1q1 / t24
      t45 = 0.1q1 / t6
      t46 = 0.1q1 / t7
      t47 = t27 ** 2
      t48 = t27 * t47
      t49 = t32 ** 2
      t50 = t32 * t49
      t51 = t50 - t48
      t52 = t12 ** 2
      t53 = mt ** 2
      t54 = t13 ** 2
      t55 = t54 ** 2
      t56 = t13 * t54
      t57 = t20 ** 2
      t58 = t23 ** 2
      t59 = t23 * t58
      t60 = t2 ** 2
      t61 = t60 ** 2
      t62 = t16 ** 2
      t63 = t44 ** 2
      t64 = t44 * t63
      t65 = t15 ** 2
      t66 = t42 ** 2
      t67 = t1 ** 2
      t68 = t1 * t67
      t69 = t11 ** 2
      t70 = t1 * t12
      t71 = t52 * t23
      t72 = t71 * t50 * t44
      t73 = t65 * t23
      t74 = t73 * t42
      t75 = t74 * t48 * t45
      t76 = t1 * t15
      t77 = t12 * t6
      t78 = t77 * t23 * t44
      t79 = t40 * t41
      t80 = t79 * t12
      t81 = t45 * t46
      t82 = t81 * t67
      t83 = t82 * t5
      t84 = t83 * t15
      t85 = t9 * t81
      t86 = t13 * t32
      t87 = t86 * t79
      t88 = t79 * t70
      t89 = t88 * t2
      t90 = t70 * t40
      t91 = t42 * t66 * t27 * t45
      t92 = t79 * t57
      t93 = t46 * t41
      t94 = t93 * t58
      t95 = t43 * t46
      t96 = t38 * t41
      t97 = t45 * t21
      t98 = t19 + t18
      t99 = t5 * t15
      t100 = t65 * t16
      t101 = t15 * t2
      t102 = t101 * t45
      t103 = t23 * t15
      t104 = t46 * t23
      t105 = t95 * t79 * t52
      t106 = t105 * t16
      t107 = t54 * t12
      t108 = t107 * t23
      t109 = t93 * t23
      t110 = t81 * t54
      t111 = t23 * t13
      t112 = t111 * t20 * t29 * t11 * (t1 * (-t91 * t96 * t65 * t13 * t5
     # * t58 * t46 + t103 * t56 * t12 * (t81 * (-t99 * t23 * t41 * t42 -
     # t12) - t80) * t48 + t104 * t54 * t12 * (t80 * t98 * t44 * t23 * t
     #5 + t100 * t45) * t50 + t106 * t5 * t58 * t64 * t32 + t17 * t101 *
     # (t81 + t79) * t47 + t17 * t46 * (t80 * t5 * t16 * t58 * t63 - t10
     #2) * t49) + t109 * (-t71 * t18 * t13 * t16 * t40 * t49 * t63 + (-t
     #70 * t54 * t15 * t65 * t45 * t50 - t71 * t64 * t32) * t40 * t62 + 
     #t74 * t54 * t47 * t45 * (t52 * t54 * t27 + t18 * t42)) * t20 * t11
     # - t110 * t65 * t58 * t42 * t47 * (t17 * t27 + t42) * t3 - t54 * t
     #65 * t58 * t41 * t42 * t48 * (t18 + t17) + t80 * t84 * t13 * (t2 *
     # (t25 * t13 * t50 - t49) - t108 * t48))
      t113 = t18 * t27
      t114 = t113 + t42
      t115 = t12 * t16
      t116 = t16 * t46
      t117 = t79 * t6
      t118 = t117 * t52
      t119 = t91 * t3 * t65
      t120 = t58 * t15
      t121 = t93 * t57
      t83 = t11 * (t120 * (t116 * t72 + t76 * (-t81 * t5 * t23 * t42 * t
     #47 * t114 + t115 * t40 * t50) * t41 - t80 * t83 * t18 * t48) * t20
     # * t56 + t59 * (t118 * t16 * t49 * t63 - t119 * t38 * t46) * t20 *
     # t54) + t69 * (t121 * t12 * t65 * t55 * t59 * t45 * t66 * t47 + t8
     #1 * t92 * t70 * t15 * t13 * t60 * t27) - t82 * t53 * t2 * t21 * t4
     #0 * t39
      t51 = t29 * t83 + t112 + t29 * (t11 * (t89 * t15 * t23 * (t85 * t4
     #7 - t49) * t20 * t54 + t58 * (t19 * t80 * t50 * (t84 + t78) + t18 
     #* (t46 * (t70 * t51 * t45 * t15 - t75 * t3 + t72) + t80 * (t78 * t
     #50 + t76 * t51))) * t20 * t56 + t87 * t20 * t52 * t6 * t16 * t59 *
     # t43 * t64) + t69 * (-t92 * t52 * t15 * t54 * t62 * t59 * t46 * t6
     #3 * t49 + t94 * (t73 * (-t52 * t62 * t40 * t50 * t44 + t91) + (-t7
     #2 * t40 + t75 + t90 * (-t50 + t48) * t45 * t15) * t14 ** 2 * t60) 
     #* t57 * t56 + t81 * t76 * t92 * t12 * t52 * t13 * t55 * t58 * t48 
     #- t86 * t81 * t92 * t70 * t15 * t60) + t97 * t40 * t39 * t2 * t67 
     #* t53 * (t96 * t15 + (-t95 * t12 + t41) * t42 * t33))
      t55 = 0.1q1 / t13
      t35 = 0.1q1 / t35
      t57 = 0.1q1 / t16
      t36 = 0.1q1 / t36
      t62 = 0.1q1 / t33
      t72 = t16 * t23
      t75 = t2 * t24
      t83 = -t72 + t75
      t84 = t22 ** 2
      t92 = t43 ** 2
      t112 = t43 * t44
      t122 = t112 * t12
      t123 = t22 * t40
      t124 = t123 * t35
      t125 = t33 * t45
      t126 = t125 * t36
      t127 = t22 * t2
      t128 = t112 * t86
      t129 = t81 * t42
      t130 = t129 * t38 * t58
      t131 = t81 + t79
      t132 = t38 ** 2
      t133 = t33 ** 2
      t134 = t23 * t38
      t135 = t29 * t32
      t136 = t12 * t2
      t137 = t136 * t40
      t138 = t31 * t47
      t139 = t23 * t2
      t140 = t45 * (-t100 * t93 * t5 * t59 * t132 * t66 * t27 + t137 * t
     #46 * (t127 * t41 * t35 + t135) + t29 * (t134 * t42 - t79 * t2) * t
     #27 * t15) * t13 + t139 * (t26 * t131 * t49 + t138 * t131) * t54
      t141 = t104 * t2
      t74 = -t100 * t96 * t54 * t59 * t42 * t47 + t141 * (t102 * t41 * t
     #27 + t80 * t2 * t32 + t74 * t38 * (t3 * t24 * t38 * t42 + t16 * t4
     #1) * t45 * t27) * t13
      t102 = t67 * t29
      t142 = t29 * t45
      t141 = t128 * t23 * t12 * (t12 * (-t117 * t1 * t54 * t32 * t58 + t
     #141 * (t86 * t1 * t22 + t117 * t2)) - t102 * t40)
      t80 = t37 * t80 * t82 * t60 * t36
      t101 = t101 * t27
      t143 = t111 * t1
      t144 = t143 * (t12 * (t87 * t81 * t1 * t2 + t23 * t54 * t1 * (-t81
     # + t79 * (-t85 - 1)) * t49) + t52 * (t112 * t127 * t117 * t111 * t
     #49 + t117 * t23 * t92 * t63 * (t127 - t111) * t32) + t101 * t41 * 
     #(t28 * t81 * t40 + (t111 * t24 * t27 - t142) * t42 * t38 * t15))
      t37 = t144 + t1 * t74 + t67 * t140 + t65 * (t110 * t42 * t38 * t58
     # * t1 * (t9 * t83 * t41 + t3 * t83) * t47 + t130 * t13 * (-t72 * t
     #1 * t3 * t38 * t42 + t75 * t96 * t67 * t5 * t42 + t60 * t3 * t41) 
     #* t27) + t103 * t54 * t67 * (t72 * (-t81 + t79 * (-t85 - 1)) + t75
     # * t85 * t79) * t47 + t70 * t46 * (t23 * (t127 * t79 * t67 * t5 * 
     #t54 * t45 * t49 + t128 * t2 * t12 * t29 * t40) + t58 * (t112 * t89
     # * t5 * t54 * t22 * t49 + t122 * t87 * t2 * (t112 * t9 * t22 + t13
     #)) + t59 * (-t88 * t5 * t54 * t92 * t63 * t32 - t122 * t56 * (t9 *
     # t79 + 1) * t49) + t41 * t2 * t60 * (t13 * t84 * t40 * t57 * t62 *
     # t35 + t124 + t126 * (t37 * t55 * t42 + 1))) + t142 * t40 * (-t32 
     #+ t27) * t13 * t68 + t77 * t79 * t61 * t84 * t57 * t46 * t62 * t35
     # + t141 + t93 * t12 * t3 * t61 * t133 * t45 * t55 * t42 * t36 + t8
     #0
      t61 = -t50 + t48
      t62 = t5 ** 2
      t74 = t67 * t62
      t77 = t6 ** 2 * t7
      t80 = t3 ** 2 * t4
      t83 = t135 * t53
      t84 = t117 * t7
      t85 = t58 * t44
      t87 = t52 * t40
      t89 = t65 * t42
      t110 = t89 * t48
      t67 = t39 * t67
      t122 = t71 * t13
      t135 = t54 * t59
      t5 = t2 * (t82 * t16 * t40 * t39 + (t11 * (t127 * t95 * t71 * t40 
     #* t32 * t20 * t44 - t130 * t100 * t27 * t20) + t67 * t40 * t45) * 
     #t41 * t13) + t94 * t13 * t1 * (-t122 * t40 * t49 * t43 * t63 + t90
     # * t15 * t45 * t61 * t54 + t73 * t91 * t132) * t29 * t62 + t135 * 
     #(t96 * t65 * t66 * t47 - t95 * t52 * t49 * t63) * t29 * t5
      t73 = t12 * t33
      t82 = t1 * t40
      t94 = t90 * t13
      t100 = t45 * t2
      t22 = t93 * t29 * t2 * (-t94 * t34 * t27 * t45 + t100 * (-t136 * t
     #133 * t55 * t42 + t82 * (-t15 * t22 - t73)) * t36 + t137 * t33 * (
     #t136 * t57 * t43 + (t73 * t42 * t43 + 1) * t45 * t1) * t35) * t21 
     #* t53
      t5 = t1 * t5 + t13 * (t109 * t76 * t2 * t45 * (t101 * t11 * t24 * 
     #t38 * t42 * t20 + t83 * t40 * t21) + t29 * (t80 * t91 * t65 * t132
     # * t46 + t79 * (-t74 * t46 - t77) * t64 * t92 * t32 * t52) * t59) 
     #+ t54 * (-t85 * t105 * t11 * t20 * t1 * t2 * t32 + t29 * ((t3 + t8
     #1 * (t74 * t41 + t80)) * t47 * t66 * t38 * t65 - t52 * t6 * t49 * 
     #t43 * t63 * (1 + t84)) * t59) + t56 * (t102 * t93 * t62 * (-t87 * 
     #t50 * t44 + t110 * t45) * t59 + t120 * t70 * t29 * (t8 * t81 * t61
     # + t84 * t61)) + t60 * (t67 * t131 + t93 * (t124 * t12 * t57 - t12
     #6 * t15 * t55) * t21 * t29 * t2 * t53) - t8 * t52 * t56 * t29 * t5
     #9 * t46 * t44 * t50 + t10 * t110 * t56 * t29 * t59 * t41 + t22 - t
     #77 * t79 * t52 * t56 * t29 * t59 * t44 * t50 + t80 * t129 * t65 * 
     #t56 * t29 * t59 * t48
      t8 = t103 * t45
      t10 = t87 * t95
      t11 = t20 * t11
      t20 = t11 * t1
      t22 = t29 * t13
      t33 = t22 * t2
      t62 = t18 * t43
      t67 = t1 * t45
      t73 = t45 * t27
      t74 = t12 * t40
      t77 = t73 * t89
      t80 = t18 * t49
      t84 = t93 * t1
      t91 = t89 * t47
      t101 = t91 * t45
      t105 = t87 * t49 * t44
      t24 = t84 * (t40 * t2 * t53 * (t128 * t12 * t58 + t67 * (t127 * t3
     #6 + t111 * (t32 + t27))) + t11 * t13 * (t13 * (-t101 * t60 * t14 *
     # t24 * t38 * t58 - t139 * t90 * t32 * t45) + t82 * t60 * t45 * (t2
     #6 * t32 + t31 * t27) + t135 * t16 * (t101 + t105)))
      t19 = t24 + t84 * (t13 * (-t11 * t100 * t72 * t76 * t40 * t27 + t2
     # * (-t11 * t123 * t52 * t16 * t32 * t43 * t63 + t73 * t53 * t15 * 
     #t38 * t42) * t58) + t54 * (t11 * t16 * (t87 * t32 * t43 * t63 + t7
     #7 * t38 * t114) * t59 + t11 * t82 * t139 * t45 * (-t138 * t18 + t2
     #6 * (-t19 - t18) * t49) + t11 * t2 * (t73 * t15 * (t28 * t14 * t40
     # * t27 - t31 * t38 * t66) - t123 * (t62 + t16) * t44 * t49 * t52) 
     #* t58) + t56 * (-t11 * t8 * t75 * t90 * t47 + t11 * t45 * (-t89 * 
     #t75 * t47 + t90 * (t19 * (t49 + t47) + t80)) * t58 + t80 * t11 * t
     #87 * t112 * t59) + t60 * t53 * (t126 * t2 * t55 + t74 * (t125 * t1
     # + t127 * t57) * t43 * t35))
      t24 = t121 * t111 * t29 * t69 * t2 * (t40 * (t32 * (-t71 * t16 * t
     #63 + t2 * t52 * t44) - t122 * t98 * t44 * t49) + t77 * (t108 * t27
     # + t111 * t114 - t2))
      t26 = t9 * t41
      t28 = -t26 - t3
      t31 = t29 * t58 * t54 * (t13 * (t12 * (t76 + t78) * t50 + t15 * (t
     #103 * t28 * t42 - t70) * t48) + t93 * (t65 * (-t100 * t30 * t42 * 
     #t47 - t73 * t134 * t66) + t105 * t2 * t30 - t94 * t15 * t47 * t45)
     # * t21 * t53 + t11 * t23 * t41 * (-t80 * t10 * t9 * t63 + t65 * t1
     #3 * t66 * t47))
      t30 = t31 + t13 * (t83 * t106 * t34 * t21 * t63 * t58 + t9 * (t118
     # * t32 * t92 * t64 - t119 * t132 * t46) * t29 * t59) + t54 * (t85 
     #* t116 * t79 * t53 * t52 * t29 * t34 * t49 * t21 + t103 * t88 * t8
     #1 * t21 * t53 * (t16 * t34 * t49 + (-t47 + t49) * t30 * t2) * t29 
     #+ (t11 * (t81 * t18 * (t26 + t3) * t47 * t66 * t38 * t65 - t52 * t
     #49 * t63 * (t62 * t117 + t116)) + t9 * (-t81 * t3 * t65 * t38 * t6
     #6 * t47 + t118 * t49 * t43 * t63)) * t29 * t59) + t56 * (t102 * t9
     #9 * t12 * (-t79 * t61 - t81 * t61) * t58 + (t9 * (-t129 * t3 * t65
     # * t48 + (t117 + t46) * t44 * t50 * t52) - t91 * t97 * t93 * t53) 
     #* t29 * t59) + t100 * t53 * t68 * t40 * t39 * (t95 + t96)
      t31 = t112 + t86
      t3 = t26 + t3
      t11 = t11 * t33 * t58 * (t65 * (-t81 * t38 * t3 * t27 * t66 + t13 
     #* (t28 * t81 - t41) * t47 * t42) + t44 * t32 * t52 * (t86 * t46 + 
     #t79 * (t9 * t31 * t46 + t31 * t6)))
      t4 = -t9 * t40 - t4
      t3 = t45 * t3
      t26 = t27 * t13
      t3 = t143 * t29 * (t52 * (t112 * t111 * (t4 * t46 - t40 * t6) * t4
     #9 - t40 * t23 * t92 * t63 * (t9 * t46 + t6) * t32) + t27 * t15 * (
     #t26 * t1 * (t79 * (t9 * t45 + t7) + t45) + t103 * (t3 * t66 * t132
     # + t26 * (t41 * t7 + t3) * t42 * t38)) + t70 * t13 * (t81 * t4 - t
     #40) * t49)
      t4 = t121 * t107 * t103 * t29 * t69 * t2 * (t40 * (t70 * t47 * t45
     # + (t76 * t50 * t45 * t23 + t85 * t12 * t50) * t16 * t14) * t13 - 
     #t25 * t48 * t45 * (t103 * t42 + t90) * t54 + t82 * t18 * t45 * (-t
     #49 + t47) - t76 * t16 * t40 * t45 * t49)
      ret = -4 * t37 - 16 * t5 - 64 * t51 - 48 * t33 * (t20 * t23 * (-t8
     #9 * t96 * t45 * t27 + t10 * t44 * t32) + t93 * (t112 * t71 * t34 *
     # t40 * t32 + t8 * t27 * (-t42 * (t103 * t38 + t34) - t82) + t40 * 
     #(t1 * t34 * t45 + t85) * t32 * t12) * t21 * t53) - 24 * t20 * t22 
     #* (t17 * t104 * t40 * (-t67 * t98 + (-t62 - t16) * t44 * t23 * t12
     #) * t49 + t73 * t41 * t15 * (t13 * (t120 * t38 * t42 * t114 + t82 
     #* t113 * t23) + t54 * (t120 * t42 * t27 + t90 * t27 * t23) - t82 *
     # t2) + t74 * (-t115 * t58 * t43 * t63 + t100 * t1) * t46 * t32) - 
     #8 * t19 + 192 * t24 + 32 * t30 + 96 * t11 + 12 * t3 + 128 * t4

      hjetmass_triangle_pmpm_s234_0_mhsq = ret/32q0/(0,1q0)
      return

      end function
