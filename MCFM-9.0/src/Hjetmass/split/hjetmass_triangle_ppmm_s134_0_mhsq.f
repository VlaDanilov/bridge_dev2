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
 

      complex*32 function hjetmass_triangle_ppmm_s134_0_mhsq
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

      t1 = zb(i2, i1)
      t2 = za(i1, i4)
      t3 = za(i3, i4)
      t4 = zb(i3, i2)
      t5 = t3 * t4
      t6 = t1 * t2 - t5
      t7 = zb(i3, i1)
      t8 = za(i2, i4)
      t9 = za(i1, i3)
      t10 = zb(i4, i1)
      t11 = zb(i4, i3)
      t12 = t9 * t7
      t13 = t2 * t10
      t14 = t3 * t11
      t15 = t13 + t14 + t12
      t16 = za(i1, i2)
      t17 = za(i2, i3)
      t18 = zb(i4, i2)
      t19 = t16 * t1
      t20 = t17 * t4
      t21 = t8 * t18
      t22 = t21 + t19 + t20
      if ( real(t22) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t22 = cg * sqrt(t22 ** 2) + t19 + t20 + t21
      t23 = 0.1q1 / t22
      t24 = t3 * t7
      t25 = 2 * t8 * t1 * t15 * t23 - t24
      t26 = -2 * t19 * t15 * t23 + t12 + t13
      t27 = 2 * t16
      t28 = t27 * t4 * t15 * t23 + t11 * t2
      t29 = t9 * t4
      t30 = t18 * t2 + t29
      t31 = 2 * t17 * t1 * t15 * t23 + t10 * t3
      t27 = t27 * t18 * t15 * t23 - t11 * t9
      t32 = t8 * t27
      t33 = t1 * (t17 * t28 + t32)
      t34 = t18 * t25
      t35 = t16 * (t31 * t4 + t34)
      t32 = t1 * (-t16 * t26 + t32)
      t36 = t9 * t1
      t37 = t18 * t3 + t36
      t38 = t1 * t26
      t34 = t16 * (-t38 + t34)
      t39 = 0.1q1 / t10
      t32 = 0.1q1 / t32
      t22 = 0.1q1 / t22
      t40 = 0.1q1 / t11
      t41 = 0.1q1 / t7
      t42 = 0.1q1 / t26
      t43 = 0.1q1 / t16
      t44 = 0.1q1 / t2
      t45 = 0.1q1 / t28
      t46 = 0.1q1 / t17
      t33 = 0.1q1 / t33
      t47 = t4 * t37
      t48 = t6 * t18
      t49 = t48 + t47
      t50 = mt ** 2
      t51 = t1 ** 2
      t52 = t1 * t51
      t53 = t45 ** 2
      t54 = t45 * t53
      t55 = t32 ** 2
      t56 = t32 * t55
      t57 = t33 ** 2
      t58 = t33 * t57
      t59 = t26 ** 2
      t60 = t26 * t59
      t61 = t46 ** 2
      t62 = t43 ** 2
      t63 = t43 * t62
      t64 = t8 ** 2
      t65 = t8 * t64
      t66 = t21 * t43
      t67 = t43 * t32
      t68 = t67 * t17
      t69 = t57 * t53
      t70 = t12 * t17
      t71 = t8 * t41
      t72 = t14 * t71
      t73 = t9 * t44
      t74 = t51 * t59
      t75 = t74 * t58
      t76 = t1 * t30
      t77 = t17 * t39
      t78 = t40 * t57
      t79 = t45 * t46
      t80 = t1 * t33
      t81 = t77 * t76
      t82 = t33 * t45
      t83 = t82 * t41
      t84 = t83 * t39
      t85 = t41 * t1
      t86 = t25 * t41
      t87 = t40 * t42
      t88 = t30 * t8
      t35 = 0.1q1 / t35
      t89 = t43 * t9
      t90 = t85 * t15 * t22 + t89
      t91 = t14 * t41
      t92 = t91 + t9
      t93 = t51 * t26
      t94 = t93 * t43
      t95 = t52 * t59
      t96 = t1 * t62
      t97 = t4 * t41
      t98 = t97 * t69
      t99 = t9 * t40
      t100 = t3 * t41
      t101 = t40 * t62
      t102 = t12 * t43
      t103 = t12 * t39 * t44
      t104 = t17 * t55
      t105 = t104 * t62
      t106 = t35 * t40
      t107 = t88 * t17
      t90 = (t17 * (t95 * t43 * t92 * t56 + t96 * t90 * t32 - t94 * t90 
     #* t55) + t60 * (t98 * t51 * t15 * t22 - t52 * t9 * t45 * t58) + t1
     #01 * (t85 * t50 * t44 * t23 + t89 * t10 * t42) + t95 * (t10 * (t99
     # + t100) + t73 * t3) * t58) * t30 * t64
      t108 = t93 * t18
      t47 = t1 * (t88 * (t75 * (-t72 * t26 * t45 + t70 * t40) + t73 * t8
     # * (t69 * t1 * t60 * t46 * (t21 * t15 * t22 - t14) + t68 * t1 * (t
     #14 * t1 * t59 * t55 + ((t66 - t1) * t32 * t26 + t43) * t22 * t15) 
     #- t14 * t51 * t60 * t45 * t58 - t14 * t60 * t54 * t61 * t33) * t39
     #) + t44 * (t64 * (t81 * t32 * t41 * t62 * t26 + t84 * (-t79 * t47 
     #+ t80 * (-t48 - t47)) * t60 + t85 * (t78 * t49 + t77 * (-t76 + t48
     #) * t43 * t55) * t59) + t78 * t77 * t1 * t49 * t59 * t8 - t87 * t8
     #6 * t6) * t23 * t50) + t107 * (t95 * (t103 + 1) * t58 * t3 + t101 
     #* (t50 * t1 * t39 * t44 * t23 + t102 * t42)) - t106 * t86 * t50 * 
     #t16 * t51 * t6 * t23 * t44 + t108 * t15 * t22 * t41 * (t69 * t14 *
     # t59 * t39 * t44 * t46 + t105) * t30 * t65 + t90
      t48 = t10 * t41
      t49 = t44 * (-t100 - t99) - t48 * t40
      t90 = t18 ** 2
      t99 = t4 ** 2
      t101 = t73 * t39
      t109 = t101 * t53
      t110 = (t91 + t9) * t44
      t111 = t18 * t39
      t112 = t77 * t51
      t113 = t112 * t44
      t114 = t77 * t1
      t115 = t17 * t1
      t116 = t26 * t39
      t117 = t116 * t45
      t118 = t77 * t56
      t119 = t118 * t43
      t120 = t77 * t41
      t121 = t44 * t8
      t122 = t76 * t22 * t15
      t123 = 0.1q1 / t15
      t124 = 0.1q1 / t8
      t125 = t48 + t73
      t126 = t22 ** 2
      t127 = t15 ** 2
      t128 = t17 ** 2
      t129 = t3 ** 2
      t130 = t71 * t4
      t131 = t8 * t125
      t132 = t93 * t71
      t133 = t26 * t8
      t134 = t80 * t26
      t135 = t44 * t40
      t72 = t30 * (t121 * (t75 * t99 * t17 * t128 * t39 * t40 + t71 * (t
     #51 * t33 * (t99 * t128 * t59 * t57 + 1) * t40 + t59 * (t51 ** 2 * 
     #t16 * t17 * t56 - t26 * t99 * t54 * t33) * t39)) * t126 * t127 + t
     #133 * (t72 * t4 * t59 * t44 * t46 * t39 * t54 * t33 - t113 * t72 *
     # t43 * t55 - t132 * t17 * (t14 * t66 * t39 * t44 + t1) * t56 + t93
     # * t17 * (-t131 * t4 * t40 + t44 * (t39 * (t8 * (t97 * t11 * t26 *
     # t45 - t18) - t20) - t130) * t3) * t58 + t1 * (t110 * t4 * t8 * t5
     #9 * t39 * t53 - t115 * t40) * t57) * t22 * t15) + t135 * t50 * (-t
     #129 * t7 * t124 * t39 * t123 + t71 * t6 * t23 * (t134 + t43))
      t136 = t78 * t17
      t49 = t88 * t33 * t22 * t15 * t52 * (t135 * t77 * t15 * t22 + t133
     # * t49 * t33 + t136 * (t103 * (-t21 - t20) - t20 - t21) * t59 + t6
     #4 * t18 * t41 * t45 * (t14 * t39 * t44 + 1) * t57 * t60) + t1 * t7
     #2 + t122 * t8 * (t26 * (t64 * (t111 * t75 * t73 * t45 + t108 * (t4
     #9 * t58 + (-t101 - t41) * t43 * t56 * t17)) + t8 * (-t110 * t77 * 
     #t52 * t56 * t26 + t82 * t4 * (t109 * t46 + (t101 + t41) * t57 * t1
     #7 * t51) * t59) - t113 * t57 * (t12 * t40 + t3)) + t121 * t1 * (t8
     # * (t114 * t59 * t40 * t58 * t90 - t111 * t98 * t60) + t85 * t59 *
     # (t58 * (-t117 + t40) + t119) * t90 * t64 + t120 * (t67 * t1 + (-t
     #115 * t58 * t45 - t69) * t60 * t99)) * t22 * t15)
      t34 = 0.1q1 / t34
      t14 = t14 + t12
      t72 = t10 ** 2
      t98 = t9 ** 2
      t99 = t129 * t11 ** 2
      t108 = t39 * t44
      t115 = t108 * t99
      t99 = t99 * t30
      t121 = t51 * t30
      t137 = t7 * t98
      t11 = t129 * t11
      t138 = t11 * t41
      t139 = t27 * t39
      t140 = t79 * t26 * t33
      t141 = t6 * t31 * t42
      t142 = t64 * t26
      t143 = t34 * t42
      t144 = 0.1q1 / t31
      t145 = 0.1q1 / t4
      t146 = t17 * t32
      t147 = t87 * t63
      t94 = -t87 * t64 * t72 * t30 * t63 * t41 + t107 * t43 * (t94 * t71
     # * t55 - t95 * t71 * t56 - t87 * t62) * t10
      t148 = t86 * t1
      t149 = t148 * t145
      t98 = t7 ** 2 * t98
      t150 = t123 * t40 * t3
      t151 = t108 * t51 * (-t133 * t37 * t40 * t33 + (t31 * t40 + t149 *
     # (t25 * t17 * t144 * t124 + 1)) * t35 * t16 * t6) * t23 * t50
      t56 = (t137 * (t69 * t51 * t60 * t39 * t46 - t147 + (t60 * t54 * t
     #61 * t33 - t146 * t63) * t39 * t1 + t95 * (-t56 * t17 * t43 + t26 
     #* t45 * t58) * t39) * t44 + t93 * (t69 * t92 * t46 * t59 + t105 * 
     #t91)) * t30 * t64
      t2 = t2 * t94 + t1 * (t64 * (t59 * (t100 * t111 * t80 * t79 * t15 
     #* t22 * t44 + t121 * (t44 * (-t137 * t40 - t138) - t2 * t72 * t40 
     #* t41) * t58) + t60 * (t99 * t39 * t41 * t44 * t54 * t61 * t33 + t
     #121 * t41 * t45 * (t115 + t13) * t58 + t115 * t69 * t76 * t41 * t4
     #6) + t139 * t100 * t93 * t79 * t15 * t22 * t44 * t33) - t48 * t2 *
     # t17 * t64 * t30 * t63 * t32 - t13 * t107 * t75 * t40 + t73 * t81 
     #* t64 * t26 * t14 * t62 * t55 - t99 * t119 * t74 * t64 * t44 * t41
     # + t135 * t129 * t4 * t123 + t108 * (t40 * (-t8 * t37 * t43 + t141
     #) + t25 * t6 * (t19 * (-t143 * t86 + t106) + t87) * t124 * t17 + t
     #85 * (-t141 * t16 * t25 * t34 + t142 * t37 * (t140 - t67))) * t23 
     #* t50) + t150 * t1 * (-t108 * t37 * t7 - t1) + t88 * t77 * t44 * (
     #-t98 * t147 + t95 * (-t98 * t40 - t11) * t58) + t151 + t56
      t13 = t77 * t55
      t37 = t13 * t43
      t54 = t18 * t26
      t56 = t77 * t7
      t63 = t1 * t44
      t58 = t63 * (t51 * (t64 * (t111 * t4 * t128 * t59 * t30 * t127 * t
     #126 * t40 * t58 + t136 * (t111 + t97) * t126 * t127 * t30 * t26) +
     # t65 * (-t97 * t77 * t18 * t60 * t30 * t127 * t126 * t45 * t58 + t
     #54 * t41 * (t37 + t78) * t126 * t127 * t30 + t97 * t17 * t18 * t59
     # * t30 * t127 * t126 * t40 * t58) + t116 * t88 * t78 * t4 * t128 *
     # t127 * t126) + t52 * (t118 * t18 * t59 * t30 * t127 * t126 * t41 
     #* t65 + t142 * t120 * t30 * t127 * t55 * t126) - t150 * t50 * t6 *
     # t23 * (t56 * t124 + 1))
      t20 = t21 + t20
      t60 = t1 * t28
      t69 = t4 * t26
      t72 = t69 + t60
      t75 = t25 ** 2
      t81 = t20 * t41
      t94 = t1 * t27
      t95 = t79 * t41
      t98 = t95 * t33
      t99 = t8 * t1
      t105 = t63 * t3
      t107 = t54 + t94
      t115 = t26 * t33
      t118 = t80 * t40
      t13 = t105 * (t95 * t111 * t74 * t15 * t22 * t65 * t27 * t57 - t11
     #6 * t80 * t50 * t40 + t116 * (t104 * t85 * t18 * t72 * t43 + t115 
     #* (t97 * t107 * t46 * t53 - t80 * t90 * t40)) * t22 * t15 * t64 + 
     #t118 * (-t139 * t1 + t41 * t72) * t22 * t15 * t8) + t105 * (t50 * 
     #t39 * (-t40 * t43 + t19 * t41 * (t1 * t35 * t144 * t145 - t143) * 
     #t124 * t75 + t99 * (-t67 * t41 * t26 + t98 * t59) + t40 * (t19 * t
     #35 + t42) * t124 * t25) + t99 * (-t111 * t26 * t40 * t33 + t13 * t
     #85 * t26 * t72 + t120 * t67 * t72 + t26 * (t40 * (t1 * (-t139 * t2
     #0 + t81 * t28) + t69 * (-t77 * t18 + t81)) + t117 * t71 * (t133 * 
     #t90 * t46 + t54 * t4 + t94 * t4)) * t57) * t22 * t15)
      t19 = t80 + t79
      t54 = t54 + t94
      t60 = t69 + t60
      t65 = t5 * t108
      t69 = t76 * t8 * t3
      t72 = t30 * t3
      t81 = t44 * t43
      t75 = t16 * (t111 * t52 * t35 * t144 * t145 * t41 * t44 * t75 + t1
     #06 * t51 * t44 * (t111 - t97) * t25) + t134 * t79 * t8 * (t134 * t
     #71 * t107 + t72 * t44) + t74 * t4 * t8 * (t104 * t41 * t43 + t78 *
     # t125) + t109 * t80 * t64 * t59 * t107 * t61 + t38 * (t17 * (t32 *
     # (-t130 * t62 + t65 * t43) + t71 * t51 * t28 * t43 * t55) + t118 *
     # t8 * (t44 * (-t111 + t97) + t80 * t125 * t28)) + t81 * (-t102 * t
     #77 * t4 * t40 - t112 * t88 * t32 - t69 * t32) + t96 * t87 * t28 * 
     #(-t131 - t17)
      t71 = t40 * (t4 * (-t63 * (t31 * t39 + t86) * t42 + t63 * (t77 + t
     #71) * t43 - t8 * (t48 + t73) * t62) - t81 * t76 * t42 * (t56 + t8)
     #) - t133 * t78 * t52 * t27 - t74 * t21 * t78 - t146 * t51 * t8 * t
     #62 * (t101 + t41) * t28
      t67 = t51 * t44 * (t8 * (t38 * (t28 * t41 - t139) * t57 + (-t111 +
     # t97) * t57 * t59) + t67 * t77 * t28) * t129
      t5 = t3 * t71 + t3 * t75 + t40 * (t17 * (-t51 * t39 * t43 - t5 * t
     #62) + t86 * t16 * t52 * t35) + t138 * t108 * t80 * t64 * t59 * t54
     # * t61 * t53 + t108 * t74 * t79 * t33 * t64 * t3 * (t18 * (t115 * 
     #t92 + t41) + t80 * t92 * t27) + t1 * (t148 * (t65 * t16 * t31 * t3
     #4 + t40) + t39 * t40 * (-t70 * t3 * t28 * t62 * t44 + t1 * t31)) *
     # t42 + t99 * (t37 * t11 * t85 * t26 * t60 * t44 + t116 * (t17 * (t
     #32 * (t85 * t4 * t43 - t29 * t62) + t89 * t1 * t60 * t55) - t12 * 
     #t78 * t1 * t54) * t44 * t3 - t85 * t40 * (t134 + t43)) + t39 * t25
     # * t52 * (-t143 * t31 * t41 + t73 * (t149 * t144 + t40) * t35) * t
     #16 - t116 * t8 * t52 * (t73 * t40 * t33 + t68 * t41) + t67 + t98 *
     # t93 * t64 * t44 * (t36 * t116 - t72)
      t7 = t72 * (t43 * (t24 * t114 * t44 * t32 + t87 * (-t10 * t8 - t17
     # * t7 - t12 * (t56 + t8) * t44) * t43) + t59 * (t110 * t80 * t64 *
     # t61 * t53 + t79 * t64 * t51 * (t110 + t48) * t57))
      ret = -32 * t47 - 64 * t49 + 24 * t122 * t83 * t44 * t59 * t64 * t
     #3 * (t79 * t4 + t80 * (t21 * t46 + t4)) + 16 * t2 - 128 * t58 + 8 
     #* t13 + 192 * t93 * t84 * t44 * t126 * t127 * t30 * t64 * (t26 * (
     #t80 * t20 + t4 * t45) + t1) - 96 * t74 * t82 * t22 * t15 * t30 * t
     #64 * (t80 * t41 + t108 * (t19 * t9 + t91 * t19)) + 56 * t113 * t88
     # * t22 * t32 * t15 * t3 * ((t66 + t1) * t32 * t26 + t43) - 4 * t5 
     #+ 80 * t108 * t132 * t50 * t6 * t23 * (t82 * t26 - t68) + 28 * t69
     # * t68 * (t43 * (t103 + 1) + t38 * (-t108 * t14 - 1) * t32) - 12 *
     # t7 + 48 * t63 * (t140 * t100 * t122 * t64 + t77 * (-t134 - t43) *
     # t23 * t40 * t6 * t50)

      hjetmass_triangle_ppmm_s134_0_mhsq = ret/32q0/(0,1q0)
      return

      end function
