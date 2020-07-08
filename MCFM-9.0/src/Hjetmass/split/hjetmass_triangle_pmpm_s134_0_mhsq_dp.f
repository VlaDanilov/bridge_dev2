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
 

      double complex function hjetmass_triangle_pmpm_s134_0_mhsq_dp 
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


      t1 = za(i1, i3)
      t2 = zb(i3, i1)
      t3 = za(i1, i4)
      t4 = zb(i4, i1)
      t5 = za(i3, i4)
      t6 = zb(i4, i3)
      t7 = t1 * t2
      t8 = t3 * t4
      t9 = t5 * t6
      t10 = t7 + t8 + t9
      t11 = za(i1, i2)
      t12 = zb(i2, i1)
      t13 = za(i2, i3)
      t14 = zb(i3, i2)
      t15 = za(i2, i4)
      t16 = zb(i4, i2)
      t17 = t11 * t12
      t18 = t13 * t14
      t19 = t15 * t16
      t20 = t17 + t18 + t19
      if ( dreal(t20) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t20 = cg * cdsqrt(t20 ** 2) + t17 + t18 + t19
      t21 = t13 * t2
      t22 = t15 * t4
      t23 = t22 + t21
      t24 = 0.1D1 / t20
      t25 = 2 * t15
      t26 = t25 * t12 * t10 * t24 - t2 * t5
      t27 = 2 * t13 * t12 * t10 * t24 + t4 * t5
      t28 = t16 * t26
      t29 = t11 * (t14 * t27 + t28)
      t30 = -2 * t17 * t10 * t24 + t7 + t8
      t31 = t11 * (-t12 * t30 + t28)
      t32 = -t25 * t14 * t10 * t24 + t2 * t3
      t25 = -t25 * t16 * t10 * t24 + t8 + t9
      t33 = 2 * t11
      t34 = t33 * t14 * t10 * t24 + t3 * t6
      t33 = t15 * (t33 * t16 * t10 * t24 - t1 * t6)
      t35 = t12 * (-t11 * t30 + t33)
      t36 = t13 * t34
      t33 = t12 * (t36 + t33)
      t37 = t15 * t6
      t38 = t11 * t2
      t39 = t38 + t37
      t40 = 0.1D1 / t6
      t41 = 0.1D1 / t5
      t42 = 0.1D1 / t13
      t43 = 0.1D1 / t34
      t44 = 0.1D1 / t14
      t20 = 0.1D1 / t20
      t45 = 0.1D1 / t30
      t46 = 0.1D1 / t27
      t47 = 0.1D1 / t3
      t31 = 0.1D1 / t31
      t33 = 0.1D1 / t33
      t35 = 0.1D1 / t35
      t29 = 0.1D1 / t29
      t48 = 0.1D1 / t11
      t49 = 0.1D1 / t4
      t50 = 0.1D1 / t12
      t51 = t11 * t26
      t52 = t15 * t30
      t53 = t52 + t51
      t54 = t19 * t31
      t55 = t54 - t45
      t56 = t19 * t50 + t11
      t57 = mt ** 2
      t58 = t29 ** 2
      t59 = t29 * t58
      t60 = t34 ** 2
      t61 = t14 ** 2
      t62 = t15 ** 2
      t63 = t15 * t62
      t64 = t31 ** 2
      t65 = t31 * t64
      t66 = t46 ** 2
      t67 = t46 * t66
      t68 = t26 ** 2
      t69 = t26 * t68
      t70 = t11 ** 2
      t71 = t70 ** 2
      t72 = t11 * t71
      t73 = t11 * t70
      t74 = t12 ** 2
      t75 = t13 * t26
      t76 = t15 * t27
      t77 = t12 * t15
      t78 = t14 * t15
      t79 = t11 * t31
      t80 = t30 * t49
      t81 = t80 * t11
      t82 = t2 * t40
      t83 = t82 * t49
      t84 = t15 * t48
      t80 = t80 * t33
      t85 = t50 * t45
      t86 = t85 * t14
      t87 = t86 * t40
      t88 = t12 * t44
      t89 = t15 * t57
      t90 = t47 * t41
      t91 = t90 * t2
      t92 = t19 * t44
      t93 = t92 + t13
      t94 = t45 ** 2
      t95 = t45 * t94
      t96 = t29 * t26
      t97 = t13 * t68
      t98 = t97 * t44 * t66
      t99 = t82 * t19
      t100 = t49 * t29
      t101 = t100 * t12
      t102 = t91 * t20
      t103 = t102 * t26
      t53 = t103 * t70 * t10 * (t13 * t61 * t68 * t50 * t40 * t94 * t31 
     #+ t78 * t83 * (t27 * (-t19 - t17) + t75 * t16) * t64 + t101 * (t12
     # * (t98 + t96 * (t52 * t13 + t51 * t93) * t46) + t99 * t29 * t53))
     # + t91 * (t11 * (t40 * (t79 * t45 * t68 * t61 * (-t75 * t31 * t56 
     #+ t76 * (t50 * t55 + t79)) + (t75 * t17 * t14 * t53 * t58 - t77 * 
     #t53 * t29 + t78 * (t76 - t75) * t31 + t18 * t70 * t12 * t68 * t64)
     # * t49 * t2) + t44 * (t52 * t74 * t13 * t68 * t29 * t49 * t66 + t8
     #1 * t74 * t62 * t16 * t68 * t58 * t46)) * t20 * t10 + t89 * (t11 *
     # (t68 * (-t88 * t46 * t29 * t49 + t87 * t31) - t83 * (t31 + t29) *
     # t26) + t15 * (t12 * (t84 * t40 * t44 * t45 * t35 * t60 - t83 * t4
     #4 * t35 * t34) + t80 * (-t52 * t42 * t43 + t82))))
      t78 = 0.1D1 / t10
      t104 = t19 + t18
      t105 = t7 * t41
      t106 = t105 + t6
      t107 = t7 * t47
      t108 = t107 + t4
      t109 = t64 - t58
      t110 = -t65 + t59
      t111 = t2 ** 2
      t112 = t2 * t111
      t113 = t20 ** 2
      t114 = t10 ** 2
      t115 = t19 * t65
      t116 = t115 * t61
      t117 = t61 * t41
      t118 = t49 * t47
      t119 = t18 * t59
      t120 = t90 * t77
      t121 = t90 * t82
      t122 = t121 * t20
      t123 = t11 * t47
      t124 = t123 * t100 * t74
      t125 = t62 * t16 ** 2
      t126 = t47 * t50
      t127 = t89 * t111 * t24 * t49 * t40 * t78
      t128 = t41 * t40
      t129 = t128 + t118
      t130 = t13 ** 2
      t131 = t23 * t14
      t132 = t131 * t12 * t2
      t133 = t4 * t40
      t134 = t133 * t117
      t135 = t128 * t4
      t136 = t12 * t61 * t23
      t137 = t70 * (t134 * t20 * t23 * t69 * t50 * t95 * t31 + t77 * t23
     # * t14 * t2 * (t118 * t109 + t128 * t109) * t20 * t26) + t71 * (-t
     #2 * t74 * t14 * t23 * t65 * t129 * t20 * t68 + t136 * t45 * t65 * 
     #(t47 + t135) * t20 * t69) + t73 * (t132 * (-t128 * t115 + t119 * t
     #129) * t20 * t68 + t23 * (t104 * t59 * t46 * t41 * t74 + t116 * t4
     #5 * t47) * t20 * t69)
      t138 = t61 * t40
      t139 = t138 * t95
      t140 = t90 * t23
      t141 = t140 * t113
      t136 = t136 * t90 * t113 * t71 * t69 * t40 * t94 * t64 - t141 * t7
     #2 * t74 * t61 * t69 * t40 * t45 * t65 - t141 * t100 * t11 * t130 *
     # t74 * t69 * t67 + t140 * ((-t130 * t61 - t125) * t49 * t59 * t46 
     #* t74 - t139 * t31) * t113 * t69 * t73 - t131 * t90 * t113 * t70 *
     # t130 * t74 * t69 * t49 * t66 * t58
      t142 = t2 * t14
      t143 = t19 * t10 * t20
      t99 = t17 * t20 * t23 * t10 * (t99 * t70 * t14 * t68 * t41 * t59 +
     # t118 * (t70 * (t119 * t6 * t12 * t46 * t69 + t142 * t59 * (-t117 
     #* t10 * t130 * t40 * t20 + t19) * t68) + t82 * t10 * t14 * t62 * t
     #41 * t20 * (t31 - t29) + t17 * t58 * t66 * t13 * (t41 * (t7 - t143
     #) + t6) * t69 - t128 * t1 * t111 * t73 * t12 * t14 * t68 * t65))
      t71 = t99 + t10 * t137 + t114 * t136 + t10 * (t70 * (t120 * t1 * t
     #111 * t14 * t23 * t40 * t49 * t109 * t20 * t26 + t122 * t1 * t61 *
     # t23 * t69 * t50 * t95 * t31) + t73 * (t23 * (t40 * (t116 * t41 * 
     #t108 * t45 - t117 * t108 * t64 * t94) + t118 * t59 * t46 * t74 * (
     #t105 * t104 + t19 * t6)) * t20 * t69 + t118 * t23 * t14 * t12 * t2
     # * (-t115 + t82 * (t19 * t110 + t119) * t41 * t1) * t20 * t68) + t
     #122 * t71 * t1 * t12 * t61 * t23 * t69 * t45 * t65 + t124 * t13 * 
     #t23 * t44 * t67 * t106 * t20 * t69) + t114 * (t73 * (-t125 * t90 *
     # t83 * t12 * t14 * t23 * t110 * t113 * t68 + t90 * t19 * t61 * t23
     # * t40 * t45 * t64 * (t45 - t54) * t113 * t69) + t90 * t83 * t113 
     #* t72 * t12 * t74 * t14 * t23 * t68 * t65) - t127 * (t126 * t14 * 
     #t23 + t88 * t39 * t41) + t127 * (t23 * t41 + t39 * t47)
      t72 = t1 ** 2 * t111
      t99 = t3 * t4 ** 2
      t104 = t30 ** 2
      t108 = t50 ** 2
      t109 = t101 * t2
      t115 = t85 * t61
      t116 = t115 * t68
      t119 = t80 * t12
      t125 = t40 * t23
      t127 = t125 * t11
      t130 = t5 * t6 ** 2
      t136 = t128 * (-t72 * t47 - t99)
      t137 = t74 * t46
      t144 = t137 * t59
      t145 = t9 * t118
      t135 = t135 * t3
      t146 = t61 * t94
      t147 = t146 * t64
      t148 = t147 * t50
      t149 = t90 * t89
      t150 = t62 * t111
      t151 = t44 ** 2
      t152 = t46 * t29
      t153 = t152 * t11
      t81 = t103 * t81 * t10 * t62 * t44 * t46 * t29 + t152 * t23 * t11 
     #* (t153 * (t6 * (-1 - t145) - t105) * t44 - t118 * (t72 * t41 + t1
     #30) * t66 * t151 - t72 * t90 * t70 * t58 * t49) * t69
      t103 = t13 * t41
      t152 = t149 * t49 * (t137 * t51 * t23 * t32 * t44 * t29 + t52 * t8
     #2 * (t12 * t39 * t35 + t131 * t33)) * t24
      t95 = t99 * t125 * t117 * t79 * t69 * t108 * t95
      t99 = t115 * t122 * t79 * t10 * t15 * t26 * (t76 - t75)
      t8 = t81 * t74 + t70 * (t100 * t102 * t10 * t74 * t15 * t44 * t46 
     #* t68 + t107 * t23 * (-t105 * t74 * t44 * t66 * t58 * t49 + t148) 
     #* t69) + t73 * (t23 * (-t144 * (t130 * t118 + t8 * t41) + (-t9 * t
     #47 + t136) * t65 * t45 * t61) * t69 + t132 * t65 * (t135 + t145) *
     # t68) + t149 * (t62 * (-t131 * t42 * t43 * t33 * t49 * t104 + t119
     # * t23 * t42) - t77 * t39 * t40 * (t88 * t2 * t49 + t84) * t35 * t
     #34 + t127 * (t142 * t26 * t31 * t49 + t109 * t32 - t116 * t31) + t
     #74 * t62 * t60 * t39 * t40 * t44 * t45 * t35 * t48 - t52 * t83 * t
     #12 * t39 * t33) * t24 - t150 * t78 * t129 - t72 * t139 * t140 * t7
     #9 * t69 * t108 + t148 * t69 * t23 * t70 * (t4 - t136) + t132 * t73
     # * (-t135 * t59 + t118 * (-t72 * t128 * t110 - t9 * t59)) * t68 + 
     #t152 - t15 * t112 * t49 * t40 * t78 * (t103 + t123) - t95 + t99
      t9 = t90 * t40 * t49
      t39 = t14 * t26
      t72 = t41 * t44
      t81 = t72 * t137
      t84 = t107 + t4
      t95 = t105 + t6
      t99 = t15 * t26
      t102 = t2 * t15
      t22 = t14 * (t51 * t111 * (-t118 * t76 + t128 * (t107 * t75 * t49 
     #- t76)) * t64 - t118 * t111 * t15 * t40 * (t75 * t41 + t23) * t31)
     # + t61 * (t85 * t40 * t26 * (t97 * t85 * t2 * t84 * t41 + t102 * t
     #23 * t47 + t99 * (t47 * (t22 + t21) - t85 * t2 * t4 * t27) * t41) 
     #* t31 + t38 * t85 * t68 * (t121 * t76 * t1 + t75 * (t128 * (-t107 
     #- t4) - t47)) * t64) + t109 * (t26 * (t23 * (t77 * t41 - t2) + t12
     # * t95 * t29 * t47 * t68 * t70 + t96 * t52 * t17 * t6 * t47) * t46
     # * t44 + t41 * t15 * (t51 * t1 * t111 * t30 * t47 * t40 * t29 + t8
     #2 * t23 - t99 * t47))
      t27 = t128 + t118
      t109 = t64 * t14
      t117 = t2 * t12
      t21 = t70 * (t100 * t91 * t1 * t74 * t151 * t66 * t69 - t39 * t9 *
     # t76 * t1 * t111 * t64 + t109 * (t76 * t87 * t4 * t41 + t21 * t27)
     # * t68) + t73 * (t117 * t58 * t27 * t68 + t81 * t58 * t69) + t90 *
     # t63 * t12 * (-t40 * t34 * t35 + t80)
      t9 = t2 * t21 + t11 * t22 + t111 * (t40 * (-t146 * t90 * t76 * t1 
     #* t68 * t31 * t108 - t86 * t23 * t26 * t31 + t119 * t90 * t62) * t
     #11 + t96 * t77 * (t129 * t29 * t30 - t9) * t70) + t112 * (t9 * t73
     # * t1 * t12 * t68 * t58 + t127 * t49 * (t31 - t29)) + t2 * (t116 *
     # t76 * t70 * t47 * t64 - t90 * t62 * (t77 * t104 * t42 * t43 * t33
     # * t49 + t39 * t40 * t31) * t11) + t124 * t68 * t2 * (t52 * t106 +
     # t51 * t6) * t66 * t151 - t90 * t62 ** 2 * t12 * t6 * t104 * t42 *
     # t49 * t43 * t33 - t36 * t9 * t150 * t12 * t35 + t140 * t12 * t63 
     #* t40 * t45 * t48 * t35 * t60 + t81 * t29 * t68 * t15 * t11 * (t38
     # * t30 * t29 + t118 * (t38 * (t7 * t30 * t29 - 1) - t37))
      t21 = -t65 + t59
      t22 = t118 * t6
      t27 = t31 * t45
      t30 = t27 * t61
      t33 = t12 * t14
      t34 = t51 * t64
      t35 = t15 * t25
      t36 = t27 * t138
      t37 = t51 * t23
      t42 = (t69 * (-t22 * t92 * t10 * t74 * t20 * t66 * t58 + t134 * (-
     #t7 + t143) * t50 * t64 * t94) + t137 * t97 * t90 * t57 * t32 * t24
     # * t58 * t49) * t23 * t70
      t21 = t37 * (t26 * (t1 * (t26 * (t74 * (t100 * t6 * t151 * t47 * t
     #67 + t118 * t11 * t44 * (-t143 * t41 + t6) * t58 * t66 + t70 * (t2
     #2 + t41) * t59 * t46) + t30 * (t70 * t47 * t64 + t128 * (t123 * t8
     #5 * t54 * t10 * t20 + t4 * t94 * t108 + t4 * t70 * t64))) * t2 + t
     #33 * t70 * (t118 * t21 + t128 * t21) * t111) + t51 * t10 * t20 * (
     #-t103 * t74 * t66 * t58 + t147 * t123)) + t90 * t11 * (t99 * t137 
     #* t25 * t58 * t49 + t36 * t26 * (t31 * (t35 - t51) + t85 * t26) + 
     #t33 * t83 * (t13 * t32 * t58 + t34 + t35 * (-t64 + t58))) * t24 * 
     #t57) - t89 * t112 * t49 * t40 * t78 * (t72 + t126) + t42 + (t69 * 
     #(t4 * t61 * t45 * t65 + t144 * t6) + t142 * t12 * t110 * t68) * t2
     #3 * t73 + t140 * t98 * t100 * t57 * t11 * t74 * t32 * t24
      t22 = t38 * t23
      t25 = t85 - t79
      t33 = t19 + t18
      t35 = t14 * t40
      t17 = t22 * t20 * t10 * (t47 * (t35 * (-t11 * t14 * t68 * t94 * t5
     #0 + t102 * t49) * t31 + t109 * t51 * t40 * (t39 * t56 * t45 - t2 *
     # t49 * (t19 + t17))) + t101 * t41 * (-t82 * t15 + t51 * (t93 * t46
     # * t26 * t12 + t82 * t33) * t29 + t97 * t88 * t66))
      t42 = t49 * t95
      t1 = t82 * t1
      t1 = t37 * t2 * (t61 * (-t26 * t40 * t94 * t108 * t84 * t31 + t34 
     #* t85 * (t47 * (t1 + t5) + t133)) + t142 * t11 * (t118 * (-t1 - t5
     #) - t40) * t64 + t12 * t29 * (t38 * t29 * (t128 * (t7 * t49 + t3) 
     #+ t49) + (t42 * t66 * t151 + t153 * (t3 * t41 + t42) * t44) * t26 
     #* t12))
      t2 = t140 * t51 * t89 * t24 * (-t96 * t12 * t49 * t46 + t35 * t27 
     #* t32)
      t3 = t100 * t74
      t3 = t141 * t51 * t15 * t114 * (t46 * (t51 * t74 * t33 * t49 * t58
     # - t3 * t15) + t36 * (t12 * t70 * t26 * t31 + t51 * t55 - t15) + t
     #3 * t75 * t66)
      ret = -32 * t21 - 8 * t53 + 64 * t71 + 16 * t8 + 4 * t9 - 48 * t22
     # * t15 * (-t81 * t96 * t10 * t20 * t49 + (-t115 * t10 * t26 * t20 
     #* t31 + (t12 * t32 * t31 + t96 * t14) * t49 * t24 * t41 * t57) * t
     #47 * t40) + 96 * t20 * t68 * t23 * t15 * t11 * t10 * (t74 * (-t100
     # * t44 * t47 * t95 * t66 + t11 * (t118 * (-t105 - t6) - t41) * t58
     # * t46) + t30 * (-t79 * t47 + t128 * (t107 * t25 + t25 * t4))) - 2
     #4 * t17 - 128 * t131 * t120 * t113 * t26 * t70 * t114 * (t28 * t40
     # * t65 * (-t117 * t49 + t39 * t45) * t70 + t49 * (t82 * t12 * t64 
     #+ (t12 * t68 * t46 * t59 + t39 * t82 * t59) * t16 * t13) * t11 + t
     #83 * t19 * (t64 - t58) - t18 * t83 * t58) + 12 * t1 + 80 * t2 + 19
     #2 * t3

      hjetmass_triangle_pmpm_s134_0_mhsq_dp = ret/32d0/(0,1d0)
      return

      end function
