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
 

      double complex function hjetmass_qqbgg_bubble_pmpm_s124_dp 
     &     (i1,i2,i3,i4,za,zb,mt,p,flip)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double complex ret
          double precision mt
          double precision p(mxpart,4)
          double complex alpha
          logical flip

      alpha = (za(i1,i2)*zb(i2,i1) + za(i1,i4)*zb(i4,i1) +
     & za(i2,i4)*zb(i4,i2))/(za(i1,i2)*zb(i2,i1) + za(i1,i4)*zb(i4,i1))

      p(5,:) = dreal(alpha)*p(i1,:)
      p(6,:) = (1d0-dreal(alpha))*p(i1,:) + p(i2,:) + p(i4,:)
      if (flip .eqv. .true.) then
          call spinoru_dp(6,p,zb,za)
      else
          call spinoru_dp(6,p,za,zb)
      end if

      t1 = za(6, i1)
      t2 = za(i2, i4)
      t3 = zb(i2, 5)
      t4 = zb(i3, i1)
      t5 = za(6, i3)
      t6 = za(i1, i2)
      t7 = za(i1, i3)
      t8 = zb(i2, i1)
      t9 = zb(i3, 5)
      t10 = za(i1, i4)
      t11 = zb(i4, i1)
      t12 = t6 * t8
      t13 = t10 * t11
      t14 = t2 * zb(i4, i2)
      t15 = t12 + t13 + t14
      t16 = zb(i4, 5)
      t17 = za(6, i2)
      t18 = za(5, i2)
      t19 = zb(i2, 6)
      t20 = za(5, i3)
      t21 = zb(i3, 6)
      t22 = za(5, i4)
      t23 = zb(i4, 6)
      t24 = t18 * t19
      t25 = t20 * t21
      t26 = t22 * t23
      t27 = t26 + t24 + t25
      t28 = zb(i1, 6)
      t29 = za(6, i4)
      t30 = t17 * t3
      t31 = t5 * t9
      t32 = t29 * t16
      t33 = t31 + t32 + t30
      t34 = t1 * t28
      t35 = t18 * t3
      t36 = t17 * t19
      t37 = t20 * t9
      t38 = t5 * t21
      t39 = t22 * t16
      t40 = t29 * t23
      t41 = -t39 + t40 - t37 + t38 + t34 - t35 + t36
      t41 = t41 ** 2
      t42 = 4 * t33 * t27 + t41
      t42 = cdsqrt(t42)
      t43 = t42 + t39 - t40 + t37 - t38 - t34 + t35 - t36
      t42 = t42 - t39 + t40 - t37 + t38 + t34 - t35 + t36
      t44 = 0.1D1 / t27
      t45 = (0.1D1 / 0.2D1)
      t46 = t45 * t22
      t47 = -t46 * t43 * t44 - t29
      t48 = t45 * t21
      t49 = -t48 * t43 * t44 + t9
      t41 = 4 * t33 * t27 + t41
      t41 = cdsqrt(t41)
      t50 = t41 + t39 - t40 + t37 - t38 - t34 + t35 - t36
      t50 = 0.1D1 / t50
      t51 = 2
      t52 = t45 * t43 * t44
      t53 = t51 * t33
      t54 = -t53 * t50 - t52
      t55 = t44 * (t43 + t42)
      t46 = t46 * t42 * t44 - t29
      t48 = t48 * t42 * t44 + t9
      t41 = t41 - t39 + t40 - t37 + t38 + t34 - t35 + t36
      t41 = 0.1D1 / t41
      t56 = t45 * t42 * t44
      t53 = t53 * t41 + t56
      t57 = 0.1D1 / t20
      t58 = t18 * t5
      t59 = -t58 * t57 + t17
      t60 = 0.1D1 / t21
      t61 = t9 * t60
      t62 = t5 * t57
      t63 = t61 + t62
      t64 = 0.1D1 / t23
      t65 = t16 * t64
      t66 = t65 + t62
      t67 = t65 * t18 + t17
      t68 = 0.1D1 / t22
      t69 = t29 * t68
      t70 = t65 + t69
      t71 = -t65 + t52
      t72 = -t65 - t56
      t73 = zb(i4, i3)
      t74 = t61 - t65
      t75 = -t69 * t18 + t17
      t76 = t69 + t52
      t77 = t69 - t56
      t78 = zb(i3, i2)
      t79 = -t52 * t18 - t17
      t80 = t56 * t18 - t17
      t81 = t5 * t59
      t82 = t20 * t17
      t58 = -t58 + t82
      t83 = t64 ** 2
      t84 = t16 ** 2
      t85 = t22 * t5
      t86 = -t85 * t57 + t29
      t87 = t62 - t69
      t88 = -t65 * t21 + t9
      t89 = -t52 * t19 + t3
      t90 = t56 * t19 + t3
      t91 = 0.1D1 / t27
      t92 = t7 * t4
      t93 = za(i3, i4) * t73
      t94 = za(i2, i3) * t78
      t95 = t26 + t24 + t25
      t96 = t91 ** 2
      t97 = t91 * t96
      t44 = (t40 + t92 + t93 + t94 + t12 + t13 + t14 + t34 + t36 + t38) 
     #* t44
      t98 = t95 * t96
      t99 = t98 * t43 + t44
      t100 = t39 - t92 - t93 - t94 - t12 - t13 - t14 + t35 + t37
      t101 = t43 ** 2
      t102 = t96 * t101
      t103 = t97 * t43 * t101
      t104 = (t39 + t35 + t37) * t91
      t105 = t104 * t43 + t30 + t31 + t32
      t106 = (0.1D1 / 0.4D1)
      t107 = (0.1D1 / 0.8D1)
      t40 = t40 + t92 + t93 + t94 + t12 + t13 + t14 + t34 + t36 + t38
      t108 = t40 * t17
      t109 = t18 ** 2 * t3
      t110 = (-t12 - t13 - t14 - t92 - t93 - t94 + t39 + t37) * t91
      t111 = t42 ** 2
      t112 = t98 * t111
      t113 = t95 * t91
      t114 = t113 * t42 + t35 + t37 + t39
      t115 = (0.3D1 / 0.4D1)
      t116 = t109 * t42 * t91 + t110 * t42 * t18 + t112 * t115 * t18 - t
     #51 * t17 * t114 + t108
      t117 = -t12 - t13 - t14 - t92 - t93 - t94 + t37 + t35
      t40 = t40 * t29
      t118 = t117 * t91
      t119 = t22 ** 2 * t16
      t26 = -t26 - t24 - t25
      t44 = t26 * t96 * t42 + t44
      t120 = t96 * t111
      t97 = t97 * t42 * t111
      t32 = t104 * t42 - t30 - t31 - t32
      t104 = t106 * t120 * t18 * t100 + t107 * t97 * t18 * t95 + t45 * t
     #42 * t17 * t44 - t17 * t32
      t98 = t98 * t101
      t26 = -t109 * t43 * t91 - t110 * t43 * t18 + t115 * t98 * t18 - t5
     #1 * t17 * (t26 * t91 * t43 + t35 + t37 + t39) + t108
      t101 = -t106 * t102 * t18 * (-t39 + t92 + t93 + t94 + t12 + t13 + 
     #t14 - t35 - t37) - t107 * t103 * t18 * t95 - t45 * t43 * t17 * t99
     # + t17 * t105
      t32 = t106 * t120 * t22 * t100 + t107 * t97 * t22 * t95 + t45 * t4
     #2 * t29 * t44 - t29 * t32
      t15 = 0.1D1 / t15
      t6 = 0.1D1 / t6
      t44 = 0.1D1 / t55
      t55 = 0.1D1 / t42
      t97 = 0.1D1 / t43
      t108 = 0.1D1 / t8
      t109 = 0.1D1 / t11
      t54 = 0.1D1 / t54
      t53 = 0.1D1 / t53
      t7 = 0.1D1 / t7
      t110 = t47 * t108
      t111 = t79 * t109
      t120 = t10 * t6
      t121 = t46 * t108
      t122 = t80 * t109
      t123 = t53 * t41
      t124 = t123 * t55
      t125 = t124 * t48
      t126 = t50 * t54
      t127 = t7 * t44
      t128 = -0.1D1 / t71
      t129 = -0.1D1 / t72
      t130 = t120 * t108
      t131 = t130 + t109
      t132 = t2 ** 2
      t133 = t78 * t108
      t134 = t123 * t48
      t135 = t29 * t108
      t136 = t17 * t109
      t137 = t97 * t55
      t138 = 0.1D1 / t63
      t72 = 0.1D1 / t72
      t71 = 0.1D1 / t71
      t139 = 0.1D1 / t5
      t140 = 0.1D1 / t66
      t141 = 0.1D1 / t70
      t66 = 0.1D1 / t66
      t142 = t46 * t73
      t143 = t78 * t80
      t144 = t126 * (-t47 * t73 + t78 * t79)
      t145 = t123 * (t143 - t142) + t144
      t146 = t135 + t136
      t147 = t108 * t5 * t86 + t109 * t81
      t148 = t138 ** 2
      t149 = t57 ** 2
      t150 = t148 * t140
      t69 = t69 * t141
      t151 = t91 * t44 * t33
      t152 = t34 * t64 * t6
      t153 = t137 * t27
      t154 = t109 * t108
      t155 = t154 * t28
      t12 = t7 * (t155 * (t152 * (t69 * (t88 * (-t51 * t65 * t17 * (t65 
     #* (t25 + t24) - t35 - t37) + t17 * (t31 + t30) + t117 * t18 * t83 
     #* t84 - t65 * (t12 + t13 + t14 + t92 + t93 + t94 + t38 + t34 + t36
     #) * t17 - t18 * (t25 + t24) * t64 * t83 * t16 * t84) * t72 ** 2 * 
     #t71 ** 2 * t96 - (t18 * t20 * t84 * t83 + t82 * t65 * t51 + t17 * 
     #t5) * t149 * t66 ** 2) + (-t140 ** 2 * t149 * t60 * t138 - t150 * 
     #t149 * t60) * t81 * t9) + t151 * t4 * (t123 * (t90 * t64 * t129 + 
     #t120) * t80 + t126 * t79 * (t89 * t64 * t128 + t120))) + (t151 * t
     #28 * (t109 * t145 + t130 * t145) + t34 * (-t6 * t147 * t60 * t139 
     #** 2 * t138 - t6 * t57 * t147 * t60 * t139 * t148) * t9 + t153 * t
     #4 * (t120 * t146 * t16 + t146 * t3)) * t15 * t2)
      t13 = 0.1D1 / t77
      t14 = 0.1D1 / t87
      t18 = 0.1D1 / t9
      t24 = 0.1D1 / t76
      t25 = t28 ** 2
      t30 = t60 ** 2
      t31 = t141 ** 2
      t36 = 0.1D1 / t70 ** 2
      t38 = 0.1D1 / t63 ** 2
      t63 = t68 ** 2
      t18 = t4 * t18
      t70 = t133 * t28
      t81 = t16 * t109
      t82 = t1 * t2
      t84 = t59 * t109
      t92 = t84 * t28 * t57
      t93 = t1 ** 2 * t2 * t4 * t6
      t94 = t57 * t66
      t117 = t82 * t6
      t3 = t64 * (t94 * t120 * t29 * t28 * t67 ** 2 * t63 * t31 * t109 -
     # t94 * t69 * t28 * t67 * t109) + t83 * (-t81 * t120 * t75 ** 2 * t
     #28 * t57 * t68 * t14 * t36 + t81 * t135 * t1 * t4 * t25 * t36 * t1
     #3 * t6 * t91 * t24 * (t10 * t75 - t82) * t63) + t9 * (t92 * t138 *
     # t60 * t140 * t64 + t82 * t30 * t6 * (t155 * t4 * t57 + (t70 * t57
     # + ((t62 * t23 + t16) * t109 - (t62 * t19 + t3) * t108) * t139 * t
     #4) * t15 * t2) * t148) + t69 * t155 * t64 * t91 * (t4 * t67 * (-t6
     #5 * t19 + t3) + t65 * (-t67 * t78 + t73 * (t65 * t22 + t29) + t120
     # * t1 * t4 * t67 * t141 * t68 - t93 * t141 * t68) * t28) * t71 * t
     #72 + t117 * (-t155 * t5 * t4 * t149 * t38 * t60 + ((-t70 * t60 + t
     #18 * ((-t61 * t23 + t16) * t109 - (-t61 * t19 + t3) * t108)) * t14
     #9 * t38 * t5 + t18 * t139 * (t108 * t3 - t81)) * t15 * t2)
      t13 = 0.1D1 / t74
      t18 = 0.1D1 / t74
      t19 = 0.1D1 / t87
      t24 = t15 * t2
      t3 = t7 * t3 + t7 * t109 * t28 * (t6 * (t57 * (t64 * (-t82 * t29 *
     # t67 * t31 * t66 * t63 + t61 * t59 * t138 * t140 * t19 * (t10 * t5
     #9 - t82) * t68) + t82 * t16 * (t135 * t28 * t73 * t13 * t141 * t66
     # * t60 + t75 * t36 * t14) * t68 * t83 + t132 * t1 * t9 * t73 * t14
     #8 * t15 * t30) + t82 * t60 * t73 * t5 * (-t24 * t38 + t61 * (-t18 
     #* t38 + t150) * t108 * t64 * t28) * t149) + t69 * t34 * t4 * t16 *
     # t108 * t83 * t91 * t71 * t72)
      t5 = -0.1D1 / t77
      t13 = -0.1D1 / t76
      t14 = t9 ** 2
      t18 = t123 * t42
      t19 = t126 * t43 * t128
      t36 = t18 * t129
      t38 = t88 * t72 * t71 * t91
      t59 = t6 * t1
      t62 = t155 * t64
      t66 = t62 * t96
      t70 = t7 * t28
      t5 = t70 * (t66 * (t36 * (-t143 + t142) + t93 * t68 * (-t18 * t5 *
     # t129 + t19 * t13) + (-t36 * (t120 * t80 * t5 * t68 + 1) + t19 * (
     #t120 * t79 * t13 * t68 + 1)) * t4 * t1 + t144 * t128 * t43) * t44 
     #* t33 + t59 * (t24 * ((t108 * t86 + t84) * t30 * t148 * t14 + t136
     # + t135) * t139 + t108 * t64 * (-t69 * t2 * t88 * t91 * t71 * t72 
     #- t28 * t67 * t109 * (t38 + t94) * t63 * t31 * t29 ** 2 + t92 * t1
     #50 * t14 * t30)))
      t13 = t34 * t91
      t14 = t13 * t108
      t18 = (-t110 - t111) * t97
      t30 = (t121 + t122) * t55
      t31 = t13 * t109
      t63 = t134 * t129
      t68 = t126 * t49
      t74 = t68 * t128
      t34 = t34 * t6
      t4 = t7 * (t34 * t9 * (t62 * t58 * t149 * t138 * t60 * t140 + t24 
     #* (-t153 * t146 + (t58 * t109 - (-t20 * t29 + t85) * t108) * t60 *
     # t57 * t139 * t138)) + t2 * (t152 * t91 * t108 * (t74 + t63) + (t1
     #26 * (t18 * t89 + t31) + t120 * (t123 * (t30 * (t56 * t23 + t16) +
     # t14) + t126 * (t18 * (-t52 * t23 + t16) + t14)) + t123 * (t30 * t
     #90 + t31)) * t15 * t4) * t44 * t33)
      t14 = t38 * t69 * t67
      t10 = t70 * t64 * (t14 * t109 + (t14 * t10 + t82 * (-t61 * t138 * 
     #t140 * t57 + t69 * (-t65 * t28 * t73 * t72 * t71 * t91 * t109 + t9
     #4))) * t108 * t6)
      t14 = t104 * t109 + t108 * t32
      t16 = (t100 * t102 * t106 * t22 - t103 * t107 * t22 * t95 - t45 * 
     #t43 * t29 * t99 + t29 * t105) * t108 + t101 * t109
      t20 = t44 ** 2
      t23 = t97 ** 2
      t31 = t55 ** 2
      t38 = t54 ** 2
      t22 = t34 * t24 * t7 * (t9 * t27 ** 2 * t23 * t31 * (t108 * t29 * 
     #t33 + t109 * t17 * t33) + (t124 * (t14 * t21 + t48 * (t108 * (t115
     # * t112 * t22 - t51 * t29 * t114 + t118 * t42 * t22 + t119 * t42 *
     # t91 - t32 * t53 + t40) + t109 * (-t104 * t53 + t116))) + (t54 * (
     #t16 * t21 + t49 * (t26 * t109 + (t115 * t98 * t22 - t118 * t43 * t
     #22 - t119 * t43 * t91 + t51 * t29 * (t113 * t43 - t35 - t37 - t39)
     # + t40) * t108)) - t49 * t16 * t38) * t97 * t50) * t91 * t20 * t33
     #)
      t32 = t41 ** 2
      t35 = t50 ** 2
      t37 = t53 ** 2
      t39 = t48 * t104
      t40 = t24 * t33
      t18 = t34 * t127 * t33 * (t40 * (t18 * t38 * t35 * t49 - t30 * t48
     # * t32 * t37) + t66 * (t41 * (t53 * (t129 * (t104 * t21 + t116 * t
     #48) - t39 * t129 ** 2) - t39 * t129 * t37) + t126 * t128 * (t101 *
     # (t49 * (t128 + t54) - t21) - t26 * t49)) * t44)
      t21 = t34 * t7 * t20 * t33 * (t66 * t44 * (t74 * t101 + t63 * t104
     #) + t24 * (t134 * t14 * t31 - t68 * t23 * t16))
      t23 = t59 * t154 * t127 * t25 * t33 ** 2 * t64 * t91 * (-t79 * t49
     # * t35 * t128 * t38 + t48 * t80 * t32 * t129 * t37)
      ret = 96 * t127 * t15 * t33 * t2 * (t126 * (-t111 * t8 - t47 + t12
     #0 * (-t110 * t11 - t79)) * t97 * t49 + t125 * (t122 * t8 + t46 + t
     #120 * (t121 * t11 + t80))) + 48 * t7 * ((t64 * (t126 * t131 * t128
     # * t79 * t49 + t134 * t80 * t129 * t131) + t6 * (-t133 * (t123 + t
     #126) + (-t123 - t126) * t109 * t73) * t15 * t132 * t1) * t91 * t44
     # * t33 * t28 + t137 * t15 * t27 * t9 * t2 * (t136 * t8 + t29 + t12
     #0 * (t135 * t11 + t17))) - 16 * t12 + 8 * t5 - 32 * t4 - 12 * t10 
     #+ 24 * t117 * t154 * t127 * t25 * t33 * t73 * t64 * t96 * (-t36 + 
     #t19) - 256 * t22 - 128 * t18 + 512 * t21 + 1024 * t40 * t13 * t7 *
     # t44 * t20 * t6 * (-t68 * t97 * t16 + t125 * t14) + 64 * t23 + 4 *
     # t3

      hjetmass_qqbgg_bubble_pmpm_s124_dp = -ret/16d0*(0,1d0)
      return

      end function
