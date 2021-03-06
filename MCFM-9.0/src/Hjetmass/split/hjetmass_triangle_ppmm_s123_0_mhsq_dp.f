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
 

      double complex function hjetmass_triangle_ppmm_s123_0_mhsq_dp 
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

      t1 = za(i1, i2)
      t2 = zb(i2, i1)
      t3 = za(i1, i3)
      t4 = zb(i3, i1)
      t5 = za(i2, i3)
      t6 = zb(i3, i2)
      t7 = t1 * t2
      t8 = t3 * t4
      t9 = t5 * t6
      t10 = t7 + t8 + t9
      t11 = za(i1, i4)
      t12 = zb(i4, i1)
      t13 = za(i2, i4)
      t14 = zb(i4, i2)
      t15 = za(i3, i4)
      t16 = zb(i4, i3)
      t17 = t11 * t12
      t18 = t13 * t14
      t19 = t15 * t16
      t20 = t18 + t19 + t17
      if ( dreal(t20) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t19 = cg * cdsqrt(t20 ** 2) + t17 + t18 + t19
      t20 = t13 * t2 + t15 * t4
      t21 = 0.1D1 / t19
      t22 = -2 * t15 * t10 * t12 * t21 - t2 * t5
      t23 = -2 * t13 * t10 * t12 * t21 + t4 * t5
      t24 = t23 * t14
      t25 = t11 * (t16 * t22 + t24)
      t26 = -2 * t17 * t10 * t21 + t7 + t8
      t24 = t11 * (t12 * t26 + t24)
      t27 = t11 * t2
      t28 = -t15 * t6 + t27
      t29 = -2 * t11 * t10
      t30 = t29 * t16 * t21 - t1 * t6
      t29 = t29 * t14 * t21 + t3 * t6
      t31 = t13 * t29
      t32 = t12 * (t15 * t30 + t31)
      t31 = t12 * (t11 * t26 + t31)
      t33 = t11 * t4
      t34 = t13 * t6 + t33
      t35 = 0.1D1 / t3
      t24 = 0.1D1 / t24
      t36 = 0.1D1 / t6
      t37 = 0.1D1 / t26
      t38 = 0.1D1 / t5
      t19 = 0.1D1 / t19
      t39 = 0.1D1 / t12
      t40 = t37 * t39
      t41 = t11 * t24
      t42 = t41 + t40
      t43 = t14 ** 2
      t44 = t14 * t43
      t45 = t22 ** 2
      t46 = t22 * t45
      t47 = t41 * t35
      t48 = t41 * t19
      t49 = t48 * t37
      t50 = 0.1D1 / t1
      t25 = 0.1D1 / t25
      t51 = 0.1D1 / t16
      t52 = t12 * t38
      t53 = t25 * t51
      t54 = t53 * t52
      t55 = t50 * t24
      t56 = t55 + t54
      t57 = t18 * t24
      t58 = t57 + t37
      t59 = t7 * t35
      t60 = t59 + t4
      t61 = t15 ** 2
      t62 = t15 * t61
      t63 = t24 ** 2
      t64 = t24 * t63
      t65 = t37 ** 2
      t66 = t37 * t65
      t67 = t13 ** 2
      t68 = t12 ** 2
      t69 = t25 ** 2
      t70 = t25 * t69
      t71 = t11 ** 2
      t72 = t71 ** 2
      t73 = t11 * t71
      t74 = t64 * t50
      t75 = t74 * t13
      t76 = t15 * t70
      t77 = t24 * t60
      t78 = t35 * t24
      t79 = t37 * t63
      t80 = t79 * t43
      t81 = t12 * t14
      t82 = t8 * t36
      t83 = t82 * t38
      t84 = t63 * t50
      t85 = t84 * t12
      t86 = t85 * t15
      t87 = t35 * t10
      t88 = t87 * t37 * t19
      t89 = t87 * t19
      t90 = t19 * t20
      t91 = t90 * t14
      t92 = 0.1D1 / t14
      t93 = 0.1D1 / t10
      t94 = t14 * t35
      t95 = t94 + t52
      t96 = t52 * t3
      t97 = -t14 - t96
      t98 = t10 ** 2
      t99 = t19 ** 2
      t100 = mt ** 2
      t101 = t2 ** 2
      t102 = t59 * t22
      t103 = t18 * t51
      t104 = t12 * t35
      t105 = t4 * t36
      t106 = t105 * t97
      t107 = t36 * t95
      t108 = t94 * t1
      t109 = t108 * t69 * t51
      t110 = t24 * t22
      t111 = t38 * t36
      t112 = t111 * t35
      t113 = t112 * t45
      t16 = t107 * t74 * t68 * t14 * t45 * t99 * t20 * t11 * t72 + t55 *
     # t17 * t99 * t61 * t20 * t14 * t38 * t36 + t113 * t43 * (t110 * t6
     #6 + t12 * t70 * (t67 * t43 * t51 + t16 * t61)) * t99 * t20 * t73
      t114 = t94 * t28 * t21
      t115 = t100 * t15
      t116 = t115 * t50
      t117 = t116 * t36
      t118 = t5 * t35
      t119 = t35 * t2 * t36
      t120 = t13 * t22
      t121 = t41 * t22
      t122 = t121 * t14
      t123 = t122 * (t10 * (t71 * (-t105 * t91 * t45 * t24 * t65 * t38 +
     # t120 * t90 * (t50 * (-t105 - t118) - t119) * t63 * t43) + t48 * t
     #15 * (t50 * (t105 + t118) + t119) * t20 * t14 - t118 * t85 * t91 *
     # t22 * t73) + t117 * t28 * t21 * t35 + t84 * t98 * t99 * t71 * t67
     # * t20 * t22 * t44 * t35 * t36)
      t16 = t123 + t10 * (t90 * t45 * t43 * (-t104 * t70 * (t103 + t15) 
     #+ (-t102 * t63 * t65 + (-t2 * t64 + (-t59 - t4) * t51 * t70 * t14)
     # * t13 * t12) * t38 * t36) * t73 + t91 * t52 * t2 * t15 * t22 * t3
     #6 * (t63 + t109) * t71 + t91 * t64 * t45 * t12 * (t50 * (-t12 + t1
     #06) - t107 * t2) * t72) + t16 * t98 + t91 * t11 * t10 * (t11 * (-t
     #77 * t14 * t36 * t66 * t38 * t39 * t46 + t86 * (t83 + 1) * t22) + 
     #t71 * (t81 * ((t74 * t10 * t67 * t14 * t19 - t76 * t60 - t8 * t75)
     # * t38 * t36 - t75) * t45 + t80 * t13 * ((t35 * t58 * t19 * t10 - 
     #t77) * t38 * t36 - t78) * t46) + t79 * t14 * t12 * ((t24 * (-t59 -
     # t4) + t88) * t38 * t36 - t78) * t46 * t73 + t89 * t61 * t14 * t36
     # * t56 + t88 * t72 * t68 * t46 * t14 * t38 * t36 * t64) - t117 * (
     #t3 * t101 * t38 * t93 * t92 + t114 * t51)
      t31 = 0.1D1 / t31
      t30 = 0.1D1 / t30
      t77 = 0.1D1 / t11
      t90 = t15 * t23
      t107 = t90 - t120
      t117 = t105 * t38
      t123 = t35 + t117
      t124 = t51 ** 2
      t125 = t51 * t124
      t126 = t29 ** 2
      t127 = t23 * t38
      t128 = t15 * t51
      t129 = t84 * t15
      t130 = t105 * t52
      t131 = t35 * t50
      t132 = t131 * t61
      t133 = t20 * t38
      t134 = t133 * t22
      t135 = t2 * t20
      t136 = t15 * t35
      t137 = t61 * t38
      t32 = 0.1D1 / t32
      t138 = 0.1D1 / t22
      t139 = t105 + t118
      t140 = t15 * t14
      t141 = t104 * t25
      t142 = t45 * t14
      t143 = t14 * t22
      t57 = t143 * (t63 * (-t40 * t117 * t18 * t45 + t90 * t50) + t141 *
     # t124) * t71 + t36 * (t57 * t15 * t22 * t38 * t50 + t128 * t50 * t
     #95 + t15 * (t131 * t29 + (-t104 * t15 * t29 * t32 + t50) * t38 * t
     #26) * t30 - t52 * t4 * (t140 * t26 * t25 + t3 * t50) * t124) * t11
     # + t142 * (t104 * t69 * t51 + t84 * t139) * t73
      t95 = t111 * t22
      t144 = t136 * t29
      t145 = t26 * t35
      t32 = t32 * t30
      t146 = t52 * t29
      t27 = t101 * (t109 * t52 * t73 * t45 * t36 + t95 * (t12 * (-t109 *
     # t15 * t26 - t53) + t143 * t63 * (t40 * t90 * t108 - t13)) * t71) 
     #+ t2 * t57 + t14 * (t2 * (-t111 * t53 * t17 * t15 * t20 + t130 * t
     #45 * t69 * t51 * t73 + t22 * (t12 * (-t128 * t26 * t123 * t69 + t1
     #17 * t124 * t25) + t129 * (-t118 * t26 + t105 * (t127 * t3 - t26))
     #) * t71) - t132 * t51) + t50 * (t12 * (-t27 * t124 - t137 * t51) -
     # t137 * t26 * t30) + t40 * t121 * (t136 * t36 * (t135 + t134) + t1
     #21 * t2 * (-t111 * t102 * t13 + t90 * t123) + t40 * t111 * t22 * t
     #2 * (t107 * t4 + t59 * t107)) * t43 - t132 * t29 * t30 + t119 * t5
     #2 * t13 * t62 * t126 * t77 * t37 * t31 - t135 * t128 * t36 * t50 *
     # t97 * t138 + t146 * t62 * (t32 * t145 + t105 * (t144 * t37 * t77 
     #+ t50) * t31)
      t57 = t83 + 1
      t62 = t36 * t15
      t83 = t11 * t35
      t97 = t13 * t38
      t102 = t135 * t36
      t105 = t50 * t31
      t106 = (t94 * t5 - t106 + t12) * t2
      t27 = t11 * (t110 * t14 * t61 * t50 * t123 + t14 * t2 * (-t141 * t
     #15 * t26 - t139 * t50) * t124 + t53 * (t12 * (t111 * t101 * t26 * 
     #t15 + t137 * t94 * t22) - t101 * t14 * t20 * t36)) + t27 - t40 * t
     #122 * t2 * (t18 * t47 * t45 + t102) + t94 * t73 * t101 * t45 * t36
     # * t63 + t105 * t61 * t12 * ((-t83 + t97) * t36 * t2 + t136) * t29
     # + t106 * t124 * t50 * t26 * t15 * t138 + t143 * t2 * (-t120 * t84
     # * t57 + t62 * (t24 * (t2 * (-t145 + t127) * t24 - t131) - t54 * t
     #35)) * t71
      t54 = t32 * t29 * t35 + t105
      t107 = t4 ** 2
      t108 = t3 ** 2
      t109 = t1 ** 2
      t119 = t39 ** 2
      t122 = t136 * t120
      t132 = t3 * t107
      t137 = t109 * t101
      t139 = t20 * t45
      t1 = t1 * t101 * t71
      t141 = t138 * t125
      t147 = t28 * t29
      t148 = t115 * t111
      t149 = t41 + t40
      t150 = t45 * t63 * t65 * t39
      t151 = t12 * t69
      t152 = t111 * t3
      t153 = t71 * t46
      t154 = t45 * t73
      t155 = t154 * t64
      t156 = t22 * t11
      t6 = (t156 * (t107 * (t152 * (t151 * t124 + t150) * t11 + t3 * t22
     # * t36 * (t52 * t70 * t51 + t74) * t71) - t115 * t112 * t53 * t21 
     #+ t1 * t22 * t35 * t64 * t36 + t33 * t150 + t137 * t113 * t79 * t1
     #1 * t149) + t118 * t11 * (t12 * t25 * t125 + t153 * t64 * t37) * t
     #6 + t131 * (-t141 + t155) * t6 * t5 ** 2) * t20 * t43
      t33 = t15 * t50
      t57 = t33 * (t57 * t15 * t93 * t2 + t36 * (t97 * t3 - t11) * t93 *
     # t101 - t111 * t100 * t28 * t21 * t30 * (t12 * t29 * t92 + t26))
      t74 = t50 * (-t148 * t121 * t21 + t148 * t21 * t51 + t155 * t12 * 
     #(t111 * t108 * t107 + t9)) * t20 * t14
      t1 = t14 * (t11 * (t94 * t22 * t20 * t11 * (t12 * (t9 * t11 * t22 
     #* t51 * t70 + (t7 + t9) * t124 * t69) + t7 * t45 * t39 * t65 * t63
     #) + t111 * (t81 * t109 * t101 * t71 * t20 * t45 * t35 * t70 * t51 
     #+ t143 * t7 * t17 * t4 * t20 * t69 * t124 + t132 * t81 * t20 * t25
     # * t125 + t110 * (t14 * (t37 * (t132 * t71 * t20 * t45 * t63 + (t2
     #3 * t35 * t19 * t39 * t61 - t122 * t19 * t39) * t10 * t2) + t139 *
     # t119 * (t137 * t35 + t132) * t66) + t1 * t12 * t20 * t22 * t63)))
     # + t141 * t20 * (-t52 * t108 * t107 * t36 - t132 * t14 * t36 - t9 
     #* t12) * t50) + t148 * (t12 * (t147 * t35 * t31 * t77 * t61 + t26 
     #* t54 * t28 * t15) + t68 * (t28 * t126 * t35 * t92 * t37 * t31 * t
     #77 * t61 + t147 * t92 * t54 * t15) - t40 * t47 * t43 * t20 * t45) 
     #* t21 + t57 + t74 + t6
      t6 = t7 * t123 + t4
      t9 = t13 * t28
      t23 = t9 * t84
      t54 = t71 * t45
      t57 = t43 * t20
      t3 = t3 * t12
      t5 = t5 * t14
      t66 = t4 * (t3 + t5) * t50
      t74 = t14 * t11
      t101 = t36 * (t136 * t50 * (t43 * t34 * t124 - t147 * t30) + (t80 
     #* t46 * t35 * t38 + t142 * t84 * t52) * t20 * t73) * t21 * t100
      t6 = t36 * (t71 * (t78 * t43 * t37 * t38 * (t40 * t20 - t9 * t24) 
     #* t46 + t142 * (-t23 * t52 + t94 * (-t23 + t52 * (-t15 * t34 - t9)
     # * t51 * t69))) + t84 * t139 * t73 * t43 * t35 - t136 * t17 * t43 
     #* t22 * t34 * t38 * t25 * t124 + t33 * t12 * (t14 * t34 * t38 * t1
     #24 + t144 * t28 * t31)) * t21 * t100 + t57 * t11 * (t12 * (t69 * (
     #-t156 * t18 * t10 * t19 * t123 * t124 + t156 * t128 * t10 * t19 * 
     #t123) + t89 * t15 * t124 * t25 + t54 * t6 * t51 * t70) + t117 * t3
     #9 * t65 * t24 * t46 * (-t18 * t48 * t10 + t7 * t149)) + t101 + t14
     # * t20 * (t81 * t71 * t4 * t22 * t124 * t69 + t17 * t4 * t14 * t12
     #5 * t25 - t66 * t141 + t155 * (t143 * t6 * t37 + t106 + t66) + t74
     # * (-t111 * t59 * t18 * t11 * t46 * t63 * t65 * t39 + t130 * t15 *
     # t25 * t124 - t153 * t35 * t63 * t65) * t19 * t10)
      t9 = t60 * t36
      t23 = t52 * t13
      t34 = t84 * t13
      t48 = t52 * t69
      t59 = t48 * t51
      t23 = t36 * (t98 * (t71 * (t57 * t52 * t15 * (t136 * t69 + t34) * 
     #t99 * t22 + t122 * (t59 + t84) * t99 * t20 * t44) + t72 * (-t133 *
     # t75 * t99 * t68 * t45 * t43 + (-t23 * t46 * t35 * t64 * t37 * t99
     # - t104 * t75 * t45 * t99) * t20 * t44) + t73 * (-t23 * t139 * t76
     # * t44 * t35 * t99 + t129 * t134 * t68 * t14 * t99 + t57 * t86 * t
     #22 * t35 * t99)) - t116 * t2 * t28 * t21 * t93 * (t96 * t92 + 1))
      t60 = t90 - t120
      t63 = t136 * t127
      t64 = t89 * t38
      t66 = t15 * t37
      t68 = t36 * t2
      t13 = t68 * (t10 * (t74 * (-t127 * t55 * t61 + t48 * t54 * t136 + 
     #t154 * t85 * t35 - t33 * t47 * t22) * t19 + t34 * t54 * (-t97 + t8
     #3) * t19 * t43 + t63 * t54 * t79 * t19 * t13 * t44 * t39) + t116 *
     # t38 * (-t29 * t92 * t30 + t121 - t51)) + t68 * (t71 * (-t64 * t79
     # * t67 * t46 * t39 * t44 - t140 * t89 * t52 * t22 * t25 * (t15 * t
     #26 * t25 + t51) + t19 * t22 * t10 * (-t78 * t97 * t45 * t65 * t39 
     #+ t63 * t110 * t65 * t39 + t15 * t13 * (-t59 * t145 + t84 * (-t145
     # + t127))) * t43) + t73 * (t64 * t45 * (t151 * t13 * t51 + t79 * t
     #60) * t43 + t84 * t81 * t19 * t22 * t10 * (-t136 * t26 + t38 * t60
     #)) + t146 * t92 * t61 * t100 * ((t66 * t31 * t77 + t32) * t35 * t2
     #9 + t105) + t140 * (t100 * t22 * t35 * t38 * (t40 * t110 + t53) + 
     #(t97 * t110 * t50 + t136 * t56 * t26) * t19 * t10) * t11)
      t19 = t156 * t25
      t26 = t111 * t91 * t17 * t25 * t10 * t2 * (t19 * (-t103 - t15) + t
     #128)
      t21 = t62 * (t135 * t87 * t49 * t22 * t43 * t39 + t52 * (t121 - t5
     #1) * t21 * t50 * t28 * t100)
      ret = -96 * t49 * t45 * t20 * t15 * t43 * t10 * (t47 + (t7 * t42 *
     # t35 + t4 * t42) * t38 * t36) - 64 * t16 - 4 * t27 - 16 * t1 - 32 
     #* t6 - 12 * t135 * (t45 * (t9 * t41 * t43 * t119 * t65 + t80 * t71
     # * (t9 + t118) * t39) + t51 * (t152 * t17 * t2 * t25 + (t82 * (t96
     # + t14) + t3 + t5) * t51 * t138 * t50)) + 128 * t23 + 8 * t13 - 80
     # * t95 * t115 * t114 * t11 * (t110 * t37 + t53 * t12) - 192 * t66 
     #* t57 * t95 * t47 * t99 * t98 * (t156 * (-t24 * (t18 + t17) - t37)
     # + t15) + 24 * t54 * t102 * t88 * t24 * t43 * (t39 * t58 + t41) + 
     #28 * t135 * t53 * t17 * t14 * (t111 * (t8 * t51 + t156 * (t8 + t7)
     # * t25) + t51 + t19) + 56 * t26 - 48 * t21

      hjetmass_triangle_ppmm_s123_0_mhsq_dp = ret/32d0/(0,1d0)
      return

      end function
