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
 

      double complex function hjetmass_triangle_pmpm_s34_mhsq_s12_dp 
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double complex ret
          double precision mt

          parameter (cg = 1d0)

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
      t11 = t3 * t4
      t12 = t5 * t6
      t13 = t7 * t8
      t14 = t9 * t10
      t15 = t11 + t12 + t13 + t14
      t16 = za(i3, i4)
      t17 = zb(i4, i3)
      t15 = -4 * t1 * t16 * t2 * t17 + t15 ** 2
      t15 = cg * cdsqrt(t15) + t11 + t12 + t13 + t14
      t18 = 2 * t1 * t2
      t19 = t18 + t15
      t20 = (0.1D1 / 0.2D1)
      t18 = t20 * t15 ** 2 - t18 * t16 * t17
      t21 = 2 * t16 * t17 + t15
      t22 = (0.1D1 / 0.4D1)
      t23 = t15 ** 2
      t24 = t1 * t16
      t25 = t24 * t2 * t17
      t26 = t22 * t23
      t27 = t26 - t25
      t27 = 0.1D1 / t27
      t28 = t11 + t13
      t29 = t15 * t1
      t30 = t29 * t2 * t27
      t31 = t26 * t1 * t2 * t27
      t32 = -t20 * t30 * t28 + t31
      t33 = t12 + t14
      t31 = -t20 * t30 * t33 + t31
      t34 = t15 * t22
      t35 = t20 * t7 * t2 * t17 + t34 * t6
      t29 = t29 * t27
      t36 = -t29 * t35
      t23 = t23 * t27
      t25 = t25 * t20 * t15 * t27
      t33 = t22 * t23 * t33 - t25
      t37 = t20 * t9 * t2 * t17 - t34 * t4
      t38 = t15 * t16 * t27
      t39 = t38 * t37
      t40 = t4 * t5 + t8 * t9
      t41 = t23 * t40
      t42 = t10 * t7 + t3 * t6
      t43 = t30 * t42
      t44 = t20 * t5 * t2 * t17 + t34 * t8
      t45 = t38 * t44
      t46 = t29 * (t20 * t3 * t2 * t17 - t34 * t10)
      t28 = t22 * t23 * t28 - t25
      t47 = t4 * t7 + t6 * t9
      t48 = t23 * t47
      t24 = t24 * t20
      t49 = t24 * t4 - t34 * t9
      t50 = t15 * t2 * t27
      t51 = t50 * t49
      t37 = -t29 * t37
      t5 = t24 * t8 + t34 * t5
      t8 = -t50 * t5
      t42 = t23 * t42
      t30 = t30 * t40
      t6 = t24 * t6 + t34 * t7
      t7 = t15 * t17 * t27
      t40 = t7 * t6
      t13 = t13 + t14
      t14 = t22 * t23 * t13 - t25
      t3 = t7 * (t24 * t10 - t34 * t3)
      t10 = -t7 * t49
      t24 = t38 * t17
      t34 = t24 * t47
      t26 = t26 * t16 * t17 * t27
      t13 = -t20 * t24 * t13 + t26
      t11 = t11 + t12
      t12 = t22 * t23 * t11 - t25
      t11 = -t20 * t24 * t11 + t26
      t23 = t29 * t44
      t6 = t50 * t6
      t5 = -t7 * t5
      t7 = 0.1D1 / t32
      t24 = 0.1D1 / t42
      t3 = 0.1D1 / t3
      t25 = 0.1D1 / t43
      t26 = 0.1D1 / t45
      t27 = 0.1D1 / t2
      t29 = 0.1D1 / t8
      t42 = 0.1D1 / t36
      t43 = 0.1D1 / t40
      t15 = 0.1D1 / t15
      t44 = 0.1D1 / t1
      t47 = 0.1D1 / t16
      t49 = 0.1D1 / t17
      t50 = t28 * t34
      t52 = t50 * t30
      t53 = t41 * t51
      t54 = t52 * t43 + t53
      t55 = t9 ** 2
      t56 = t32 ** 2
      t18 = 0.1D1 / t18 ** 2
      t57 = mt ** 2
      t58 = t39 ** 2
      t59 = t4 ** 2
      t60 = t25 ** 2
      t61 = t28 ** 2
      t62 = t51 ** 2
      t63 = t19 ** 2
      t64 = t26 ** 2
      t65 = t26 * t64
      t66 = t29 ** 2
      t67 = t29 * t66
      t68 = t39 * t36
      t69 = t32 * t34
      t1 = t1 * t2
      t2 = t1 * t21 ** 2
      t70 = t1 * t21
      t71 = t70 * t41 * t47 * t49
      t72 = t19 * t30
      t73 = t19 * t12
      t74 = t73 * t41
      t75 = t18 * t58
      t76 = t75 * t32
      t77 = t21 * t11
      t78 = t57 * t39
      t79 = t61 * t48 * t51
      t80 = t39 * t28
      t81 = t19 * t51
      t82 = t81 * t21
      t83 = t21 * t33
      t84 = t82 * t18
      t85 = t53 * t40
      t86 = t49 * t47
      t87 = t16 * t17
      t88 = t87 * t63
      t89 = t86 * t2
      t90 = t86 * t51 * t28
      t91 = t18 * (-t88 * t28 * t48 * t51 * t30 * t27 * t44 + t89 * (t85
     # + t52) * t39) * t29 + t90 * (-t57 * t45 * t48 * t30 * t27 * t15 *
     # t44 - t2 * t54 * t18 * t11) * t66
      t92 = t36 * t51
      t93 = t57 * t51
      t94 = t44 * t27
      t95 = t94 * (t93 * t30 * t10 * t40 * t15 * t47 * t49 * t7 * t3 + (
     #t39 * (-t56 * t31 * t48 * t26 * t18 * t60 + t92 * t30 * t26 * t18 
     #* t25) - t76 * t30 * t12 * t25 * t64) * t63 * t17 * t16)
      t16 = t24 * t91 + t24 * (t29 * (t28 * (-t55 * t51 * t44 * t47 + (t
     #78 * t48 * t30 * t15 * t44 * t47 - t51 * t59) * t49 * t27) + t21 *
     # (-t79 * t10 * t43 + t80 * t48 * t30 - t41 * t62 * t40) * t18 * t1
     #9 + t79 * t16 * t17 * t63 * t37 * t27 * t43 * t44 * t18) + t82 * t
     #28 * t12 * t18 * t54 * t66) + t25 * (t76 * t21 * (t11 * (t72 + t71
     #) - t74) * t64 + t77 * t74 * t75 * t56 * t42 * t65 + t39 * (t2 * t
     #56 * t10 * t34 * t42 * t47 * t49 * t18 + t19 * (-t56 * t37 * t34 *
     # t42 - t68 * t30 + t69 * t30) * t18 * t21 + t69 * t57 * t30 * t27 
     #* t15 * t44 * t47 * t49) * t26) + t83 * t19 * t18 * t48 * t56 * t3
     #9 * t60 * t26 + t84 * t30 * t10 * t40 * t7 * t3 + t95
      t17 = 0.1D1 / t46
      t46 = 0.1D1 / t28
      t54 = 0.1D1 / t51
      t63 = t19 * t37
      t76 = t21 * t10
      t79 = -t63 * t94 * t87 + t76
      t59 = t59 * t27 * t49
      t91 = t55 * t44 * t47 + t59
      t95 = t69 - t68
      t96 = t7 ** 2
      t97 = t30 ** 2
      t98 = t40 ** 2
      t99 = t3 ** 2
      t100 = t3 * t99
      t101 = t94 * t86
      t102 = t101 * t57
      t103 = t102 * t4 * t9
      t104 = t94 * t88
      t105 = t104 * t51 * t37
      t106 = t63 * t21
      t107 = t2 * t10
      t108 = t94 * t57
      t109 = t108 * t34
      t110 = t19 * t14
      t111 = t21 * t13
      t112 = t86 * t28
      t113 = t30 * t28
      t114 = t21 * t41
      t115 = t29 * t24
      t116 = t115 * t51
      t117 = t68 * t19
      t118 = t89 * t39
      t119 = t32 * t39
      t93 = t101 * t93
      t70 = t86 * t70
      t101 = t70 * t10
      t120 = t101 - t63
      t121 = t11 ** 2
      t122 = t12 ** 2
      t123 = t48 * t37
      t124 = t28 * t48
      t125 = t68 * t41
      t126 = t39 * t41
      t127 = t86 * t15
      t128 = t127 * t57
      t129 = t19 * t48
      t130 = t18 * t21
      t131 = t94 * (t128 * (t37 * (-t125 * t17 * t46 - t52 * t7 * t3) + 
     #(t126 * t8 * t34 * t64 - t53 * t34 * t26) * t25 * t32 - t124 * t11
     #5 * t53) + t88 * (t56 * (-t41 * t122 * t25 * t65 * t42 * t58 - t12
     #3 * t25 * t26 * t42 * t39) + t113 * (t124 * t51 * t12 * t24 * t66 
     #* t43 + t37 * t40 * t14 * t7 * t99) + t119 * t53 * t12 * t25 * t64
     #) * t18)
      t45 = t18 * (t117 * (-t114 * t37 + t30 * t79) * t46 * t17 + t116 *
     # t28 * t21 * (-t71 * t34 + (-t113 * t11 * t29 * t43 - t41) * t48 *
     # t19)) + t7 * (t112 * t30 * t40 * (t107 * t13 * t18 + t109 * t23 *
     # t15) * t99 - t111 * t110 * t18 * t97 * t28 * t98 * t54 * t100 + t
     #30 * (-t105 * t40 * t18 - t103 * t40 + t106 * (t39 * t40 - t50) * 
     #t18) * t3) + (t118 * t95 * t18 * t41 - t84 * t95 * t41 + t119 * t9
     #1) * t26 * t25 + t93 * t97 * t98 * t15 * t96 * t3 + t131 + t130 * 
     #(t40 * (t7 * (-t101 * t39 * t30 * t3 - t113 * t19 * (t10 * t14 + t
     #13 * t37) * t99) + t113 * t51 * t120 * t3 * t96) + t56 * (t129 * t
     #10 * t25 * t26 * t42 * t39 - t71 * t121 * t25 * t65 * t42 * t58) -
     # t119 * t81 * t41 * t11 * t25 * t64 - t113 * t81 * t48 * t45 * t66
     # * t24)
      t65 = 0.1D1 / t39
      t71 = t17 ** 2
      t95 = t17 * t71
      t101 = t36 ** 2
      t131 = t41 ** 2
      t132 = t46 ** 2
      t133 = t32 * t46
      t134 = t133 * t101 * t131
      t135 = t41 * t30
      t136 = t135 * t14
      t137 = t133 * t41
      t138 = t119 * t12 * t26
      t139 = t89 * t18
      t27 = t59 * t30 * t40 * t7 * t3 + (t55 * t30 * t40 * t7 * t3 + t49
     # * t15 * t27 * t46 * t17 * t41 * t57 * (t126 * t101 * t46 + (t36 *
     # t5 * t17 - t10) * t48 * t32)) * t47 * t44 + t139 * ((t7 * (t50 * 
     #t13 * t54 * t99 - t39 * t34 * t54 * t3) + t50 * t3 * t96) * t40 * 
     #t97 + (t7 * (-t126 * t46 * t3 + t41 * t13 * t99) + t53 * t3 * t96)
     # * t98 * t30 - t69 * t41 * t10 * t17 * t46) + t104 * (t101 * (t17 
     #* (-t53 * t30 * t7 * t46 + t126 * t30 * t132) + t136 * t46 * t71) 
     #+ t48 * (-t113 * t37 * t7 * t3 + t32 * t41 * t25 * t26 * (t138 * t
     #42 - t51))) * t18 + t130 * (-t126 * t56 * t48 * t11 * t25 * t64 * 
     #t42 - t50 * t97 * t40 * t14 * t7 * t99 * t54 - t135 * t101 * t13 *
     # t71 * t46 - t134 * t14 * t13 * t65 * t95 - t137 * t48 * t10 * t17
     # - t136 * t98 * t7 * t99) * t19
      t44 = t104 * t14 ** 2 + t89 * t13 ** 2
      t47 = t94 * t87
      t49 = t21 * t39
      t55 = t47 * t81 - t49
      t59 = -t47 * t110 + t111
      t87 = t19 * t101 * t131 * t18
      t126 = t70 * t13
      t135 = t129 * t21
      t136 = t28 * t13
      t140 = t18 * t3
      t70 = t70 * t39
      t141 = t70 - t81
      t142 = t30 * t40
      t143 = t130 * t19
      t144 = t130 * t40 * t30
      t145 = t114 * t18
      t146 = t130 * t69
      t147 = t119 * t26
      t148 = t147 * t25
      t149 = t47 * t19
      t150 = t86 * t77 * t1
      t8 = t98 * (-t136 * t89 * t97 * t99 * t96 + t110 * t97 * (-t49 * t
     #54 + t149) * t99 * t7) + t114 * t69 * t39 * t25 * t64 * (t150 * t3
     #2 * t42 - t19 * t8) - t89 * t85 * t10 * t7 * t3
      t4 = t17 * (-t145 * t117 * t32 * t10 * t132 + t94 * t41 * (-t92 * 
     #t88 * t37 * t18 + t86 * (-t69 * t37 * t15 - t4 * t9 * t36) * t57) 
     #* t46) + t18 * t8 + t17 * (t125 * t104 * t18 * t32 * t37 * t132 + 
     #t145 * (-t141 * t10 * t36 - t63 * t69) * t46) + t7 * (t3 * (-t102 
     #* t30 * t15 * (t124 * t10 + t142 * t34) - t139 * t52 * t10 + t143 
     #* (t40 * (-t34 * t97 + t53 * t37) - t113 * t48 * t10)) + t144 * ((
     #t70 * t54 - t19) * t13 * t40 * t30 - t50 * t19 * t23) * t99 - t97 
     #* t28 * t98 * t54 * t18 * t44 * t100) + t96 * (t130 * t110 * t97 *
     # t28 * t98 * t99 + t130 * t97 * t98 * t141 * t3) + t148 * (t94 * (
     #-t88 * t30 * t18 - t128 * t41) * t48 - t146 * t74 * t42 * t26) - t
     #140 * t90 * t2 * t97 * t98 * t7 * t96 + t133 * t114 * t36 * t13 * 
     #t18 * t120 * t71 - t137 * t110 * t36 * t18 * t79 * t71
      t8 = 0.1D1 / t30
      t9 = 0.1D1 / t41
      t13 = t104 * t31 ** 2 + t89 * t33 ** 2
      t23 = t149 * t31
      t49 = -t83 + t23
      t1 = t86 * t83 * t1
      t52 = t19 * t31
      t57 = -t52 + t1
      t64 = t24 ** 2
      t74 = t24 * t64
      t85 = t21 * t57
      t90 = t18 * t8
      t75 = t61 * (t93 * t6 * t10 * t15 * t29 * t64 + t90 * (t85 * t10 +
     # t63 * t49) * t29 * t64 * t62) - t19 * t40 * t18 * t49 * t29 * t64
     # * t62 * t28 - t75 * t56 * t36 * t9 * t25 * t60 * t26 * t13
      t35 = -t32 * t37 * t38 * t35
      t38 = t30 * t39
      t93 = t32 * t33
      t79 = t73 * t79
      t96 = t77 * t120
      t97 = t43 * t18 * t62
      t11 = t61 * (t24 * (-t77 * t73 * t62 * t30 * t43 * t18 * t67 + t97
     # * (-t79 + t96) * t66) + t81 * t48 * t18 * t49 * t29 * t64) + t148
     # * (t94 * (t88 * ((t92 - t138) * t25 * t31 - t38 * t37 * t9) * t18
     # - t127 * t78 * t37) + t130 * (t19 * (t25 * (t147 * (t11 * t31 + t
     #12 * t33) - t35 - t92 * t33) + t38 * t10 * t9) + t70 * (-t93 * t11
     # * t25 * t26 + t10))) + (t51 * (-t118 * t10 * t18 + t106 * t39 * t
     #18 - t103) * t29 + t72 * t62 * t18 * (-t47 * t73 + t77) * t66) * t
     #24 * t28
      t12 = t50 * t43 + t53 * t8
      t38 = t19 * t18
      t12 = t24 * (t51 * (t112 * (t107 * t12 * t18 + t109 * t30 * t15) +
     # t38 * (-t21 * t28 * t37 * t12 + t142 * t55)) * t29 + t97 * t61 * 
     #t30 * (t104 * t122 + t89 * t121) * t67 + t130 * t80 * t51 * t30 * 
     #(t73 - t150) * t66) + t148 * (t105 * t18 - t84 * t10 + t103 - t146
     # * t57 * t25 + t147 * t18 * (t79 - t96) * t42)
      t24 = -t52 + t1
      t1 = t18 * (t61 * (-t82 * t10 * t6 * t64 * t29 + (t73 * (t83 - t23
     #) + t77 * (t52 - t1)) * t66 * t64 * t62) + t9 * t26 * t60 * t58 * 
     #t56 * (t76 * t24 + t63 * (-t83 + t23)) + t80 * t21 * t40 * t24 * t
     #29 * t64 * t51)
      t6 = t130 * t34 * t57 * t29 * t64 * t51 * t61 + (t145 * t40 * t8 *
     # t57 * t29 * t64 + t115 * t94 * (t88 * t37 * t18 - t128 * t10)) * 
     #t62 * t28 + t147 * t60 * (t68 * (t72 * t49 * t9 - t85) * t18 + t35
     # * t102 * t15)
      t23 = t145 * t69
      t2 = t41 * (t36 * (t17 * (t46 * (t38 * (-t47 * t129 * t30 + t114 *
     # (t34 * t51 * t65 - t48)) - t86 * t41 * (t108 * t48 * t15 + t2 * t
     #34 * t18)) - t23 * t19 * t132) + t23 * t65 * (-t110 + t126) * t46 
     #* t71) + t144 * t129 * t7 * t3)
      t15 = t129 * t18 * t46 * t17 * t131 * t36 * (t111 * t32 * t17 * t6
     #5 + t149 * (t65 * (-t32 * t14 * t17 + t51) - t133))
      ret = -16 * t16 - t20 * t2 - t15 * t22 - 8 * t45 - 512 * t75 + 64 
     #* t11 + 32 * t12 + 2048 * t90 * t62 * t61 * t40 * t29 * t74 * t13 
     #+ 256 * t1 - 128 * t6 - 2 * t27 + t17 * (t41 * (-t123 * t104 * t32
     # * t18 + t72 * t130 * t36 * t34 + t36 * t91) * t46 + t87 * t55 * t
     #132 - t119 * t104 * t18 * t131 * t101 * t46 * t132) + t71 * (t18 *
     # t41 * t36 * (-t135 * t32 * t5 + (t21 * (-t110 + t126) - t81 * t59
     # * t65) * t41 * t36) * t46 + t87 * t32 * t59 * t132) - t134 * t65 
     #* t18 * t44 * t95 + t140 * t7 * t40 * t30 * (-t89 * t41 * t34 + t1
     #04 * (t28 * t14 * t3 * t54 - 1) * t30 * t48 + t135 * (-t28 * t7 + 
     #t54 * (-t136 * t3 + t39)) * t30) + 4 * t4 + t143 * (t28 * (192 * t
     #115 * t10 * t62 + 48 * t116 * t30 * t34) + t148 * (t39 * (-1024 * 
     #t93 * t31 * t36 * t60 * t9 - 96 * t37) - 6 * t48 * t41) - 4096 * t
     #62 * t31 * t61 * t33 * t40 * t8 * t29 * t74)

      hjetmass_triangle_pmpm_s34_mhsq_s12_dp = ret/32d0/(0,1d0)
      return

      end function
