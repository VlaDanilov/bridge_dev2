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
 

      complex*32 function hjetmass_triangle_pmpm_s124_0_mhsq
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


      t1 = za(i1, i2)
      t2 = zb(i2, i1)
      t3 = za(i1, i4)
      t4 = zb(i4, i1)
      t5 = za(i2, i4)
      t6 = zb(i4, i2)
      t7 = t1 * t2
      t8 = t3 * t4
      t9 = t5 * t6
      t10 = t7 + t8 + t9
      t11 = za(i1, i3)
      t12 = zb(i3, i1)
      t13 = za(i2, i3)
      t14 = zb(i3, i2)
      t15 = za(i3, i4)
      t16 = zb(i4, i3)
      t17 = t11 * t12
      t18 = t13 * t14
      t19 = t15 * t16
      t20 = t17 + t18 + t19
      if ( real(t20) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t20 = cg * sqrt(t20 ** 2) + t17 + t18 + t19
      t21 = 0.1q1 / t20
      t22 = -2 * t17 * t10 * t21 + t7 + t8
      t23 = t1 * t14
      t24 = -t16 * t3 + t23
      t25 = 2 * t11 * t10
      t26 = t25 * t16 * t21 - t1 * t6
      t27 = t11 * t22
      t28 = t12 * (-t15 * t26 + t27)
      t25 = -t25 * t14 * t21 + t3 * t6
      t29 = t13 * t25
      t27 = t12 * (t27 + t29)
      t30 = t1 * t12
      t31 = t5 * t16
      t32 = t31 + t30
      t33 = t5 * t14
      t34 = t12 * t3 + t33
      t35 = -2 * t13 * t12 * t10 * t21 + t4 * t5
      t36 = t12 * t22
      t37 = t11 * (t14 * t35 + t36)
      t38 = 2 * t12 * t15 * t10 * t21 - t2 * t5
      t39 = t11 * (t16 * t38 - t36)
      t20 = 0.1q1 / t20
      t40 = 0.1q1 / t4
      t41 = 0.1q1 / t3
      t42 = 0.1q1 / t2
      t43 = 0.1q1 / t26
      t44 = 0.1q1 / t13
      t45 = 0.1q1 / t15
      t46 = 0.1q1 / t1
      t47 = 0.1q1 / t25
      t27 = 0.1q1 / t27
      t28 = 0.1q1 / t28
      t48 = t12 * t28
      t49 = t43 * t45
      t50 = -t49 + t48
      t51 = t9 * t42
      t52 = t51 + t1
      t53 = t47 ** 2
      t54 = t47 * t53
      t55 = t12 ** 2
      t56 = t12 * t55
      t57 = t22 ** 2
      t58 = t22 * t57
      t59 = t27 ** 2
      t60 = t27 * t59
      t61 = t13 ** 2
      t62 = t15 ** 2
      t63 = t40 * t41
      t64 = t63 * t52
      t65 = t64 * t44
      t66 = t48 * t40
      t67 = t42 * t46
      t68 = t28 * t43
      t69 = t68 * t61
      t70 = t20 * t24
      t71 = t70 * t57 * t55
      t72 = t9 * t40
      t73 = t72 + t3
      t72 = -t72 - t3
      t74 = -t9 * t41 - t4
      t75 = t45 ** 2
      t76 = t44 ** 2
      t77 = t43 ** 2
      t78 = t43 * t77
      t79 = t68 * t12
      t80 = t22 * t27
      t81 = t80 * t41
      t82 = t28 * t13
      t83 = t12 * t5
      t84 = t36 * t24 * t5
      t85 = 0.1q1 / t16
      t37 = 0.1q1 / t37
      t39 = 0.1q1 / t39
      t86 = 0.1q1 / t35
      t87 = 0.1q1 / t38
      t88 = 0.1q1 / t14
      t89 = t5 ** 2
      t90 = t5 * t89
      t91 = t35 ** 2
      t92 = t28 ** 2
      t93 = t28 * t92
      t94 = t16 * t22
      t95 = t81 * t47 * t44
      t96 = t83 * t63 * t47 * t44
      t97 = t67 * t73
      t98 = t67 * t89
      t99 = t98 * t6
      t100 = t12 * t26
      t101 = t11 * t55
      t102 = t67 + t63
      t103 = t33 * t41 + t12
      t104 = t5 * t46
      t105 = t12 * t38
      t106 = t41 * t37
      t107 = t3 * t55
      t108 = (-t67 - t63) * t12
      t109 = t42 * t41
      t110 = t22 * t89
      t111 = t42 * t40
      t112 = t51 * t46
      t113 = t14 * t22
      t114 = t12 * t25
      t115 = t24 * t5
      t116 = t67 * t3
      t103 = t12 * (t15 * (t110 * (-t94 * t102 + t108 * t26) * t59 + t10
     #9 * t5 * (t36 * t40 + t104 * (t94 * t40 + t24)) * t27) + t111 * t1
     #1 * (-t104 * t103 * t39 * t35 + t106 * t38 * t5 * (t16 * (t105 * t
     #86 * t88 - t104) - t12) + t107 * t91 * t46 * t87 * t85 * t39)) + t
     #13 * (t104 * t66 * (t103 * t42 * t22 - t115 * t41) + t36 * t89 * (
     #t114 * t102 + t113 * (t63 * (-t112 - 1) - t67)) * t92) + t116 * t5
     # * t57 * t77 * t75 * (-t113 + t114) * t28 * t61
      t117 = t63 * t9 + 1
      t118 = t41 * (-t27 + t28)
      t110 = t46 * (-t107 * t111 * t69 * t57 * t45 - t110 * t68 * t13 * 
     #t24 * t45 + t84 * t69 * t45 * t40 + t118 * t24 * t90) + t62 * (-t1
     #00 * t63 * t89 * t6 * t57 * t76 * t42 * t53 * t27 + t83 * t59 * t4
     #2 * t57 * (-t100 * t117 - t94) * t47 * t44) - t95 * t89 * t15 * t2
     #4
      t33 = t12 * t110 + t12 * t103 + t12 * (t61 * (t57 * (t83 * t67 * t
     #68 * t40 * t45 * (t49 * t9 * t25 - t14) + t49 * t25 * t55 * t5 * (
     #t67 * t72 - t40) * t92) + t58 * (-t99 * t14 * t40 * t77 * t75 * t2
     #8 + t49 * t83 * t14 * (t97 + t40) * t92)) + t62 * (t96 * t57 * (-t
     #30 * t26 + t94 * (-t51 - t1)) * t59 + t95 * (t83 * t24 * t42 + t36
     # * (t42 * (t31 + t30) - t1 * t5 * t26 * t44 * t47) * t40 - t31 * t
     #44 * t47 * t52 * t40 * t57)) + t101 * t42 * t40 * (t30 * t38 ** 2 
     #* t41 * t88 * t86 * t37 + t33 * t91 * t46 * t87 * t85 * t39) + t29
     # * t67 * t63 * t90 * t55 * t6 * t22 * t92 - t36 * t67 * t63 * t90 
     #* t6 * (t100 + t94) * t59 * t15)
      t49 = t113 - t114
      t86 = t18 + t17
      t87 = t17 * t44
      t91 = t87 + t14
      t103 = mt ** 2
      t110 = t17 * t45
      t119 = t17 * t22
      t68 = t68 * t13
      t120 = t41 * t47
      t121 = t14 * t57
      t122 = t82 * t46
      t123 = t39 * t45
      t124 = t15 * t41
      t125 = t111 * t55
      t126 = t11 * t56
      t127 = t16 * t43
      t128 = t69 * t46
      t129 = t128 * t45
      t26 = t111 * t36 * t5 * ((t121 * t62 * t16 * t41 * t44 * t53 * t27
     # + t129 * t22 * (t126 * t25 * t28 + t127 * t49) + t124 * t83 * t46
     # * t16 * (-t29 * t12 * t92 - t27 + t113 * t13 * (t92 + t59))) * t2
     #0 * t10 + t83 * t103 * t46 * t41 * (t27 + t28)) + t125 * t5 * (t10
     #3 * (t11 * (t123 * t46 * (-t12 * t35 * t85 + t5 * t38 * t41) + t10
     #6 * (t104 * t35 - t105 * t88) * t44) + t57 * (t124 * t44 * t47 * t
     #27 - t68 * t45 * t46)) + (t15 * (-t83 * t26 * t41 * t46 * t27 + t1
     #04 * t36 * t41 * (t94 * t11 + t26 * t86) * t59) + t62 * (t121 * t2
     #6 * t41 * t44 * t53 * t27 + t120 * t57 * (t100 * t91 + t94 * t91) 
     #* t59) + t122 * ((t119 * t49 * t28 - t113 + t114) * t41 * t5 + t68
     # * t57 * (t114 * t16 + t113 * (-t110 - t16)))) * t20 * t10)
      t29 = t112 + 1
      t49 = -t93 + t60
      t68 = t16 ** 2
      t91 = t14 ** 2
      t105 = t63 * t51
      t106 = t13 * t55
      t114 = t106 * t5
      t130 = t62 * t58
      t131 = t70 * t12 * t10
      t132 = 0.1q1 / t10
      t133 = t15 * t32
      t134 = t133 * t44 - t34
      t135 = t20 ** 2
      t136 = t10 ** 2
      t137 = t16 * t40
      t138 = t11 ** 2 * t55
      t139 = t61 * t91
      t140 = t62 * t68
      t141 = t55 * t61
      t115 = t106 * t115
      t142 = t13 * t34
      t143 = t142 * t45
      t23 = t10 * (t106 * t99 * t63 * t70 * t15 * (t92 - t59) * t22 + t7
     #0 * (t41 * (t30 * t14 * t62 * t40 * t59 * t53 + t23 * t62 * t40 * 
     #t44 * t27 * t54 + t62 * t40 * t55 * (t18 * t1 + t17 * t52) * t60 *
     # t47) + t69 * (t116 * t16 * t77 * t45 - t137 * t112 * t79 + t126 *
     # t116 * t92 + t19 * (t116 + t40) * t92 * t55)) * t58 + t67 * t71 *
     # t5 * t13 * t15 * t60 * t86) + t136 * (t115 * t67 * t63 * t15 * (t
     #138 * t93 - t139 * t60) * t135 * t57 + t111 * (-t68 * t61 * t78 * 
     #t46 * t28 - t141 * (t140 + t138) * t93 * t46 * t43 - t120 * t55 * 
     #t62 * t60 * (t139 + t138)) * t135 * t24 * t58) + t46 * t21 * t41 *
     # t132 * t89 * t103 * (t134 * t42 + t143 * t40)
      t99 = t16 * t77 * t45
      t138 = t83 * t40 * (t10 * (t13 * (t130 * t109 * t70 * t6 * t55 * t
     #14 * t47 * t60 + t124 * t112 * t71 * (t17 * t49 - t19 * t93)) + t6
     #7 * t70 * t69 * t6 * (t19 * t55 * t92 + t126 * t92 + t99) * t58) +
     # t106 * t124 * t67 * (-t138 * t60 + t140 * t93) * t135 * t24 * t57
     # * t136 - t104 * t103 * t32 * t21 * t41 * t132)
      t23 = t138 + t12 * t23 + t131 * (t22 * (-t116 * t61 * t12 * t57 * 
     #t16 * t77 * t92 + t114 * (t63 * (t59 * (t18 * t80 * t29 - 1) + t92
     # + t119 * t49) + t67 * (t92 * (-t119 * t28 + 1) - t59)) * t15 + t1
     #26 * t61 * t57 * t40 * t43 * t93 + t22 * (t22 * (t105 * t12 * t14 
     #* t59 * t53 + t105 * t14 * t44 * t27 * t54 + t55 * t42 * t86 * t60
     # * t47) - t106 * t31 * t93 * t102) * t62) + t111 * (t91 * (-t13 * 
     #t12 * t62 * t58 * t41 * t59 * t53 - t130 * t41 * t27 * t54) - t101
     # * t14 * t62 * t58 * t41 * t53 * t59 + t13 * t46 * t12 * (t13 * (t
     #15 * t58 * t77 * t92 * t68 + t17 * t58 * t77 * t92 * t16) + t118 *
     # t83 * t15)) * t20 * t10)
      t31 = t6 ** 2
      t68 = t63 * t89
      t91 = t1 ** 2 * t2
      t105 = t91 * t63
      t116 = t68 * t31
      t118 = t89 * t31
      t119 = t44 * t53
      t126 = t12 * t15
      t138 = t104 * t41
      t139 = t3 ** 2 * t4
      t140 = t139 * t67
      t144 = t24 * t22
      t4 = t22 * (t144 * (t107 * t62 * t4 * t22 * t42 * t47 * t60 + t114
     # * t15 * (t7 * t63 * t49 + t67 * (-t116 * t93 + t8 * t49)) + t69 *
     # t22 * (t140 * t75 * t77 + (t7 * t40 + t140) * t92 * t55 - t79 * t
     #45 * t73) + t119 * t36 * t1 * t62 * t59 + t105 * t80 * t62 * t76 *
     # t54) + t111 * t83 * (t120 * t100 * t62 * t44 * t27 - t129 * t113)
     # * t20 * t10) + t111 * t12 * (t41 * (t17 * t34 * t38 * t37 * t88 +
     # t104 * (-t133 * t22 * t28 + t142 * (t123 * t11 * t38 - t80))) + t
     #17 * t39 * t46 * (-t143 + t32) * t85 * t35) * t21 * t103
      t4 = t12 * t4 + t12 * (t43 * (t141 * t67 * t10 * t20 * t5 * t22 * 
     #t25 * t45 * t40 * t28 + t141 * t98 * t31 * t58 * t24 * t40 * t93) 
     #+ t57 * (t106 * t67 * t63 * t90 * t15 * t31 * t24 * t60 + t96 * t1
     #0 * t20 * t62 * t16 * t42 * t27) + t58 * (t119 * t24 * t62 * t12 *
     # (t116 * t42 + t105 + t51) * t59 + t63 * t55 * t62 * t24 * t47 * (
     #t118 * t42 + t91) * t60 + t62 * t68 * t31 * t24 * t76 * t42 * t54 
     #* t27) - t67 * t12 * t61 * t58 * t24 * t45 * (t118 * t40 + t139) *
     # t92 * t77 + t63 * t17 * t42 * (-t126 * t34 * t37 * t44 * t88 + t1
     #04 * (t134 * t37 - t34 * t39)) * t21 * t35 * t103 + t98 * t61 * t3
     #1 * t58 * t24 * t75 * t40 * t78 * t28 + t132 * t89 * (-t138 * (t14
     # * t42 + t137) + t108))
      t7 = t70 * t10
      t8 = t22 * t55 * (t111 * (t124 * (t22 * t34 * t47 - t104 * t32) * 
     #t27 + t122 * (-t5 * t34 * t41 + (t143 - t32) * t43 * t22) - t95 * 
     #t62 * t32) * t21 * t103 + t7 * t5 * (t109 * t62 * t44 * t47 * t27 
     #+ t129 * t40))
      t11 = t7 * t83 * (t42 * (t124 * (-t119 * t121 * t15 + t83 * t46) *
     # t27 + t124 * t36 * (-t104 * t86 + (-t87 - t14) * t47 * t22 * t15)
     # * t59) + t122 * t40 * (t13 * (t99 * t57 - t79 * t57 * (t110 + t16
     #)) + t83 * t41 * (t22 * t28 * (t19 + t17) - 1)))
      t20 = t9 * t27
      t25 = t103 * t42 * t21
      t31 = t55 * t24
      t35 = t67 * t103
      t37 = t61 * t40 * t43
      t38 = t13 * t15
      t6 = t63 * t103 * t90 * t46 * t45 * t132 + t31 * t57 * (-t38 * t93
     # * t102 * t6 * t89 + t37 * t35 * t22 * t21 * t92 + t9 * t22 * (t62
     # * t42 * t47 * t60 + t37 * t93)) + t126 * t42 * t57 * (-t7 * t113 
     #* t15 * t53 * t59 + t63 * (-t133 * t113 * t47 * t59 + t104 * (-t14
     # * t32 * t59 - t16 * t34 * t92) * t13) * t21 * t103) + t112 * t3 *
     # t61 * t58 * t24 * t75 * t78 * t28
      t6 = t12 * t6 + t12 * (t3 * (-t112 * t12 * t61 * t58 * t24 * t45 *
     # t92 * t77 + t141 * t58 * t24 * t93 * t29 * t43) + t115 * (t117 * 
     #t60 - t93 + t67 * (t9 * t60 + t63 * (-t92 + t59) * t21 * t103)) * 
     #t57 * t15 + t131 * t58 * (t137 * t61 * t77 * t92 + t17 * (t97 * t9
     #2 * t45 * t77 * t61 - t65 * t62 * t53 * t59)) + t137 * t67 * t69 *
     # t103 * t58 * t21 * t50 * t34 + t35 * t90 * t41 * t44 * t132 + t13
     #0 * t27 * t47 * (t1 * t55 * t24 * t59 + t63 * (t9 * t1 * t24 * t76
     # * t53 + (-t25 * t14 * t32 + t20 * t30 * t24) * t47 * t44 + t31 * 
     #t27 * (t20 * t1 + t25))))
      t7 = t38 * t144 * t111 * t135 * t56 * t136 * (t138 * (-t18 * t59 +
     # t19 * t92) + t17 * (-t138 * t59 + t113 * t41 * (t15 * t22 * t47 +
     # t104 * t13) * t60 + t94 * t46 * (t22 * t13 * t43 - t124 * t5) * t
     #93 + t138 * t92))
      ret = 96 * t71 * t10 * (t62 * (t65 * t27 * t53 + t12 * (t64 + t42)
     # * t59 * t47) + t69 * (t67 * (t9 * t50 * t40 + t3 * t50) + t66)) +
     # 12 * t84 * (t62 * (-t81 * t76 * t53 * t52 + t36 * t47 * t44 * (-t
     #1 * t41 + t42 * t74) * t59) + t82 * (t48 * t5 * (t63 * (t9 * t46 +
     # t2) + t46) + (t46 * t73 * t75 * t77 + t79 * (-t2 * t40 + t46 * t7
     #2) * t45) * t22 * t13) + t83 * (t67 * t74 - t41) * t59 * t15) - 4 
     #* t33 - 8 * t26 - 64 * t23 + 16 * t4 - 48 * t8 - 192 * t144 * t125
     # * t135 * t136 * (t41 * (t47 * (t36 * t62 * t86 * t59 - t12 * t62 
     #* t27) + t80 * t14 * t62 * t53) + t128 * (t22 * (t101 * t28 + t19 
     #* t48 - t127) - t12)) - 24 * t11 + 32 * t6 + 128 * t7

      hjetmass_triangle_pmpm_s124_0_mhsq = ret/32q0/(0,1q0)
      return

      end function
