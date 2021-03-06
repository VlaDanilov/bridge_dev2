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
 

      double complex function hjetmass_qqbgg_triangle_pmmp_s34_mhsq_s12_dp
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double complex ret
          double precision mt

      t1 = za(i2, i3)
      t2 = za(i1, i3)
      t3 = zb(i3, i1)
      t4 = zb(i3, i2)
      t5 = za(i1, i4)
      t6 = zb(i4, i1)
      t7 = za(i2, i4)
      t8 = zb(i4, i2)
      t9 = t2 * t3
      t10 = t1 * t4
      t11 = t5 * t6
      t12 = t7 * t8
      t13 = t11 + t12 + t9 + t10
      t14 = za(i1, i2)
      t15 = za(i3, i4)
      t16 = zb(i2, i1)
      t17 = zb(i4, i3)
      t13 = -4 * t14 * t15 * t16 * t17 + t13 ** 2
      t13 = t11 + t12 + cdsqrt(t13) + t9 + t10
      t18 = (0.1D1 / 0.4D1)
      t19 = t13 ** 2
      t20 = t14 * t15
      t21 = t20 * t16 * t17
      t22 = t19 * t18 - t21
      t22 = 0.1D1 / t22
      t23 = (0.1D1 / 0.2D1)
      t24 = t13 * t18
      t25 = t23 * t2 * t16 * t17 - t24 * t8
      t26 = t13 * t14 * t22
      t27 = t26 * t25
      t28 = t23 * t1 * t16 * t17 + t24 * t6
      t29 = t13 * t15 * t22
      t30 = t29 * t28
      t31 = t1 * t3 + t6 * t7
      t19 = t19 * t22
      t32 = t19 * t31
      t33 = t2 * t4 + t5 * t8
      t34 = t26 * t16
      t35 = t11 + t9
      t36 = t19 * t18
      t37 = t36 * t14 * t16
      t38 = -t23 * t34 * t35 + t37
      t39 = t21 * t23 * t13 * t22
      t35 = t18 * t19 * t35 - t39
      t40 = t38 * t35
      t41 = -t32 * t34 * t33 / 8
      t42 = -t30 * t26 * (t23 * t5 * t16 * t17 + t24 * t4) + t40 + t41
      t41 = t27 * t29 * (t23 * t7 * t16 * t17 - t24 * t3) + t40 + t41
      t43 = 2 * t14 * t16 + t13
      t21 = t23 * t13 ** 2 - 2 * t21
      t44 = t12 + t10
      t37 = -t23 * t34 * t44 + t37
      t9 = t9 + t10
      t10 = t18 * t19 * t9 - t39
      t45 = 2 * t15 * t17 + t13
      t31 = t34 * t31
      t20 = t20 * t23
      t34 = t24 * t1 + t20 * t6
      t46 = t13 * t16 * t22
      t47 = -t46 * t34
      t48 = -t24 * t2 + t20 * t8
      t22 = t13 * t17 * t22
      t49 = t22 * t48
      t2 = t1 * t8 + t2 * t6
      t8 = t29 * t17
      t50 = t8 * t2
      t33 = -t31 * t19 * t33 / 8
      t4 = t47 * t22 * (t20 * t4 + t24 * t5) + t40 + t33
      t5 = t26 * t28
      t3 = t49 * t46 * (t20 * t3 - t24 * t7) + t40 + t33
      t7 = t18 * t19 * t44 - t39
      t11 = t11 + t12
      t12 = t36 * t15 * t17
      t20 = -t23 * t8 * t11 + t12
      t2 = t19 * t2
      t22 = -t22 * t34
      t8 = -t23 * t8 * t9 + t12
      t9 = t18 * t19 * t11 - t39
      t11 = -t29 * t25
      t4 = 0.1D1 / t4
      t12 = 0.1D1 / t41
      t18 = 0.1D1 / t21
      t3 = 0.1D1 / t3
      t19 = 0.1D1 / t17
      t21 = 0.1D1 / t42
      t24 = 0.1D1 / t15
      t25 = 0.1D1 / t14
      t26 = 0.1D1 / t16
      t28 = t21 ** 2
      t29 = t21 * t28
      t33 = t12 ** 2
      t34 = t12 * t33
      t36 = t28 - t33
      t39 = t7 ** 2
      t41 = t20 ** 2 + t39
      t42 = t37 + t9
      t44 = -t7 - t8
      t51 = t4 - t3
      t52 = t35 + t7
      t53 = t38 + t37
      t54 = t35 ** 2
      t55 = t35 * t54
      t56 = t27 ** 2
      t57 = t30 ** 2
      t58 = t3 ** 2
      t59 = t3 * t58
      t60 = t1 ** 2
      t61 = t8 ** 2
      t62 = t6 ** 2
      t63 = t4 ** 2
      t64 = t4 * t63
      t65 = t18 ** 2
      t66 = t45 ** 2
      t67 = t20 * t2
      t68 = t50 * t53
      t69 = t49 * t47
      t70 = t30 * t27
      t71 = t70 * t32
      t72 = t71 * t38
      t73 = t72 * t5 * t11
      t74 = t73 * t33
      t75 = t69 * t35
      t76 = t75 * t31
      t77 = t76 * t50 * t58
      t78 = t56 * t57 * t32
      t79 = t24 * t19
      t80 = t79 * t16 * t14
      t81 = t52 * t33
      t82 = t51 * t25 * t24
      t13 = 0.1D1 / t13
      t83 = t35 + t7 + t8
      t84 = t38 + t37 + t10
      t85 = -t35 - t7 - t20
      t86 = t38 ** 2
      t87 = t10 ** 2
      t88 = mt ** 2
      t89 = t30 * t38
      t90 = t35 * t28
      t91 = t38 * (t29 - t34)
      t92 = t85 * t33
      t93 = t83 * t28
      t94 = t18 * t45
      t95 = t18 * t47
      t96 = t95 * t43
      t97 = t80 * t30
      t15 = t15 * t17
      t17 = t15 * t43 ** 2
      t98 = t17 * t65
      t14 = t14 * t16
      t16 = t79 * t35
      t99 = t89 * t21
      t100 = t35 * t47
      t101 = t100 * t4
      t102 = t38 + t37 + t9
      t103 = t37 + t10
      t104 = t9 ** 2
      t105 = t20 * t58
      t106 = t44 * t63
      t107 = t94 * t47
      t108 = t107 * t54 * t63
      t109 = t30 * t4
      t110 = t109 * t94
      t111 = t30 * t2
      t112 = t50 * t9
      t113 = t5 * t28
      t114 = t42 * t33
      t115 = t28 * t103
      t116 = t5 * t12
      t117 = t98 * t26 * t25
      t118 = t117 * t30 * (-t116 * t38 * t2 + (t102 * t33 * t47 + t89 * 
     #(-t104 * t34 + t29 * t86)) * t32 * t56)
      t11 = t45 * (t38 * (t111 * t22 * t12 - t77) + t78 * (-t114 + t115)
     # + t72 * (t113 * t11 + t112 * t33)) * t65 * t43
      t119 = t94 * t80 * (-t94 * t78 * t38 * t54 * t34 + (t100 * (t94 * 
     #(t105 + t106) + t63) - t108 + t110) * t50 * t49 * t31)
      t120 = t70 * t60 * t25 * t24
      t121 = t37 ** 2
      t112 = (-t69 * t62 * t3 + t70 * (t12 * (t112 * t88 * t32 * t38 * t
     #13 * t24 * t25 * t12 - t62) + t21 * t62)) * t26 * t19
      t122 = t117 * (-t100 * t5 * t2 * t3 + t91 * t78 * t121)
      t11 = t45 * (t80 * (t78 * t36 - t77) * t18 + t43 * (t69 * (t58 * (
     #-t42 * t50 + t67) + t68 * t63) * t31 * t35 - t74) * t65) + t69 * (
     #t62 * t26 * t19 * t4 + t82 * t60) + t80 * (t69 * t31 * t54 * t50 *
     # t58 + t78 * (-t38 * t41 * t34 + t28 * t44 + t81 + (t54 + t39 + t6
     #1) * t38 * t29)) * t65 * t66 + t98 * t2 * t5 * t25 * t26 * (t99 + 
     #t101) + (-t97 * t66 * t50 * t65 * t3 + t16 * t47 * (t67 * t88 * t1
     #3 * t25 * t26 + t14 * t66 * t50 * t7 * t65) * t58) * t49 * t31 + t
     #30 * (t96 * (-t28 + t33 + t94 * (t92 + t93)) + t97 * (t91 + (t20 *
     # t33 - t90) * t65 * t66) + t98 * (-t84 * t28 * t47 + t89 * (t29 * 
     #t87 - t34 * t86)) * t26 * t25) * t32 * t56 + t120 * t21 + t11 + t1
     #18 + t119 - t99 * t43 * t65 * t45 * t2 * t22 + t122 + t112 - t120 
     #* t12
      t44 = t49 ** 2
      t56 = t15 * t43
      t60 = t56 * t37 * t25 * t26
      t62 = t94 * t35
      t67 = t80 * t45
      t77 = t18 * t43
      t91 = t43 * t63
      t46 = -t22 * t46 * t48
      t48 = t45 * t65
      t97 = t78 * t18
      t112 = t97 * t38
      t118 = -t59 + t64
      t119 = -t7 - t20
      t120 = t28 * t38
      t122 = t78 * t38
      t123 = t79 * t88
      t124 = t123 * t13
      t125 = t124 * t25 * t26
      t126 = t79 * t18
      t127 = t48 * t43
      t39 = -t69 * t50 * t3 * (t125 + t127) + t95 * (t14 * (t79 * t30 * 
     #t63 * t45 + t126 * (-t30 * t8 * t63 - t55 * t47 * t59 + t100 * (t3
     #9 * t64 - t41 * t59)) * t66) + t17 * t95 * t26 * t25 * (t38 * t58 
     #+ (t87 + t86 + t121) * t64 * t35)) * t44
      t41 = t88 * t13
      t87 = t125 * t2
      t24 = t78 * (t38 * (t19 * (-t41 * t24 * t25 * t26 * t28 + t94 * t1
     #4 * (t35 + t8) * t24 * t29) + t77 * t33 * (t12 * t42 + t94)) + t77
     # * (t34 * (t56 * t9 * t18 * t25 * t26 + 1) - t29) * t86)
      t128 = t100 * (-t48 * t51 * t5 * t50 * t43 - t87 * t22 * t3)
      t24 = t31 * t39 + t31 * (t48 * t69 * t43 * t50 * t4 + t47 * (t80 *
     # (-t63 * t52 * t65 * t66 * t30 + t100 * t118) + t65 * t47 * t43 * 
     #(t45 * (t119 * t58 - t106) - t91 * t15 * t84 * t26 * t25)) * t44) 
     #+ t43 * (t48 * (t100 * t2 * t22 * t51 + t78 * (t86 * (t119 * t34 +
     # t29 * t83) - t120)) - t112 * t29 * t103) - t79 * t73 * t88 * t13 
     #* t25 * t26 * t28 + t48 * (t31 * (t47 * t30 * (-t67 * t85 * t58 + 
     #t91 * t53) * t44 + t91 * t75 * (-t10 * t50 + t46)) + t67 * t22 * t
     #50 * (-t100 * t3 + t99)) + t112 * (t67 * (t29 * t7 + t34 * (t7 * (
     #-1 + t62) - t20 - t35)) + t77 * (t34 * (t9 * (-t20 * t45 + t60) - 
     #t40 * t45) + t45 * (t10 * t83 + t37 * t83) * t29)) + t125 * (t101 
     #* t2 * t22 + t122 * t33) + t24 + t128
      t39 = t80 * t66
      t40 = t39 * t65
      t73 = -t37 - t10
      t75 = t47 ** 2
      t78 = t39 * t35
      t106 = t12 - t21
      t112 = t17 * t5 * t25 * t26
      t128 = t77 * t31
      t129 = t38 * t22
      t130 = t31 * t54
      t131 = t38 * t57
      t132 = t131 * t27
      t133 = t63 - t58
      t134 = t35 * t22
      t135 = t132 * t22
      t136 = t69 * t30
      t41 = t41 * t35
      t137 = t35 * t5
      t81 = t79 * (-t41 * t31 * t75 * t44 * t58 * t25 * t26 + t14 * (t94
     # * (-t130 * t118 * t44 * t75 - t134 * t133 * t49 * t75 - t135 * t2
     #8) + (t75 * (t130 * (t119 * t59 + t64 * t7) * t44 - t134 * t58 * t
     #52 * t49) + t136 * t22 * t3 + t135 * (t90 - t81)) * t65 * t66))
      t62 = t77 * (t27 * (t95 * (t56 * t106 * t26 * t25 * t5 - t45 * t22
     # * t12) * t30 + t131 * (t33 * (t129 * t94 - t5) + t113)) + t62 * t
     #49 * t75 * (t133 * t49 * t31 + t137 * t58))
      t90 = -t38 - t37 - t9
      t118 = t7 + t8
      t119 = t4 * t35
      t131 = t31 * t35
      t134 = (t37 + t9) * t33
      t136 = t75 * (t43 * (t56 * t18 * t4 * t5 * (t119 * (t38 + t37 + t1
     #0) - 1) * t26 * t25 + t94 * (t4 - t3) * t22 + t137 * (t63 - t58)) 
     #* t49 + t131 * (t43 * (t38 * t64 + t59 * t90) + (t80 * t59 - t77 *
     # t64 * (t38 + t37)) * t7 * t45) * t44) - t136 * t94 * t43 * t5 * t
     #3 + t18 * t30 * (-t45 * t43 * t38 * t5 * t50 * t21 + t70 * (t43 * 
     #(t45 * (t21 * (-t86 * t22 * t21 + t5) - t116) - t134 * t56 * t38 *
     # t5 * t26 * t25) + t129 * t39 * t28 * t118))
      t137 = t69 * t4
      t138 = t38 * t50
      t139 = t49 * t31
      t140 = t49 * t75
      t141 = t20 * t59
      t142 = t8 * t64
      t143 = t38 * t37
      t105 = t94 * (t77 * (t30 * (t5 * (t138 * t12 + t137) + t27 * t47 *
     # t22 * t21) + t140 * (t139 * t9 * t59 - t5 * t63) * t54 + t140 * (
     #t5 * (-t63 * t7 + t105) + t58 * t42 * t22) * t35) + t80 * (t75 * (
     #t94 * t54 * t22 * t63 * t49 + t131 * (t141 + t94 * (-t141 + t142) 
     #* t7) * t44) - t110 * t69 * t22 + t135 * t33))
      t36 = (-t124 * t99 * t50 * t5 + t98 * (t57 * t27 * t86 * t5 * t36 
     #+ t131 * (t59 * (-t53 * t9 - t143) + t64 * (t10 * t53 + t143)) * t
     #44 * t75)) * t26 * t25
      t14 = t140 * t18 * t35 * (t112 * t58 * t18 * t90 + t14 * (t66 * (t
     #126 * t22 * t118 * t63 + t142 * t139 * t16 * t18) - t139 * t79 * t
     #118 * t64 * t45))
      t1 = t14 + t18 * t136 + t75 * (t128 * t35 * (t103 * t64 + t94 * (t
     #59 * (t102 * t20 + t102 * t7) + t64 * (-t10 * t7 - t8 * t84))) * t
     #44 + t65 * t3 * (-t78 * t20 * t22 * t3 + t112) * t49) + t123 * t69
     # * t26 * t25 * (t46 * t35 * t31 * t63 * t13 - t51 * t6 * t1) + t70
     # * (t123 * t1 * t6 * t25 * t26 * t106 + (t113 * t17 * t38 * t25 * 
     #t26 * t103 + t39 * (t12 * (-t38 * t20 * t12 + 1) - t21) * t22) * t
     #65 * t30) + t125 * (t76 * (-t46 * t58 + t69 * t63) + t116 * t89 * 
     #t50) + t127 * (t75 * (t35 * (t58 * (t5 * t7 + t129) + t63 * (-t22 
     #* t84 - t5 * t8)) * t49 + t130 * (t53 * t59 - t64 * t84) * t44) + 
     #t132 * (t22 * (t114 - t115) + t5 * (-t92 - t93))) + t81 + t62 + t3
     #6 + t105
      t5 = t40 * t50
      t6 = t111 * t106
      t5 = t71 * (t45 * (t138 * t80 * t28 * t18 + t43 * t38 * (-t8 * t2 
     #* t28 - t68 * t33) * t65) + t5 * (t21 * (-t38 * t21 * t83 + 1) - t
     #12) - t120 * t87 * t8) + t2 * t31 * (t41 * t82 * t50 * t26 * t19 +
     # t137 * t117 * (-t119 * t84 + 1) + (-t100 * t63 * t18 + t48 * (t10
     #0 * t7 * t133 - t47 * t54 * t58 - t109)) * t49 * t43) + t128 * t2 
     #* t49 * (t15 * t96 * t3 * (t3 * t35 * t102 - 1) * t26 * t25 + t108
     # + t94 * t30 * t3 + t100 * (t94 * t8 * t63 + t58)) + (t45 * (-t89 
     #* t80 * t33 * t50 * t18 + t43 * (t50 * (t106 * t47 + t30 * (t120 *
     # t103 + t28 * t86)) - t6) * t65) - t6 * t125 - t92 * t5 * t89) * t
     #32 * t27
      t6 = t33 - t28
      t14 = -t12 + t21
      t3 = t2 * (t43 * (t72 * t6 * t18 + t48 * (t72 * (t28 * (t35 + t7) 
     #+ t33 * (-t35 - t7 - t20)) + t131 * t50 * (-t4 + t3))) + t117 * t3
     #2 * t27 * (t14 * t47 + t30 * (t38 * (t28 * t73 + t134) + t6 * t86)
     #) + t138 * t125 * t32 * t14)
      t6 = 8
      ret = t23 * t138 * t127 * t32 * t2 * t106 + t6 * (t18 * (t31 * (t9
     #1 * t107 * t30 * t10 * t44 + (t43 * t58 + t63 * (t78 * t61 * t4 * 
     #t18 - t43) + t17 * t59 * t18 * t35 * (-t104 - t86 - t121) * t26 * 
     #t25) * t44 * t75) + t97 * t17 * t26 * t25 * t86 * (t29 * t73 + t34
     # * t37)) + t31 * (t95 * (t67 * (t107 * t55 * t64 - t30 * t58) + t7
     #7 * t58 * (t56 * t42 * t26 * t25 * t47 - t45 * t30 * t102)) * t44 
     #+ t69 * (-t46 * t127 * t35 * t58 + (-t16 * t10 * t13 * t25 * t26 *
     # t63 + t79 * t26 * t25 * t13 * t4) * t50 * t88)) + t24 + t125 * t7
     #4 + t122 * t65 * (t43 * (-t60 * t10 * t29 + (t37 * t85 - t52 * t9)
     # * t34 * t45) + t39 * (t29 * (-t35 * t7 - t52 * t8) + t52 * t34 * 
     #t20)) + t40 * t22 * t50 * (-t89 * t12 + t101)) + 16 * t1 - 4 * t11
     # + 2 * t5 - t3

      hjetmass_qqbgg_triangle_pmmp_s34_mhsq_s12_dp = ret/32d0/(0,1d0)
      return

      end function
