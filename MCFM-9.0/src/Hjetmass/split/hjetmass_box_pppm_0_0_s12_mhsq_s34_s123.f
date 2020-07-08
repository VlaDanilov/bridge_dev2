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
 

      complex*32 function hjetmass_box_pppm_0_0_s12_mhsq_s34_s123
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          
          complex*32 ret
          double precision mt

      t1 = za(i2, i4)
      t2 = zb(i2, i1)
      t3 = zb(i3, i2)
      t4 = za(i1, i2)
      t5 = za(i1, i3)
      t6 = zb(i4, i1)
      t7 = za(i3, i4)
      t8 = zb(i3, i1)
      t9 = za(i1, i4)
      t10 = za(i2, i3)
      t11 = zb(i4, i3)
      t12 = t1 * t3
      t13 = t9 * t8
      t14 = t13 + t12
      t15 = t7 * t11
      t16 = t15 * t14
      t17 = t4 * t2
      t18 = t5 * t8
      t19 = t10 * t3
      t20 = t15 * (t19 + t17 + t18)
      t21 = zb(i4, i2)
      t22 = t5 * t6
      t23 = t10 * t21
      t24 = t23 + t22
      t25 = mt ** 2
      t26 = 16 * t25 * t24 * t16 + 4 * t20 ** 2
      t26 = sqrt(t26)
      t16 = 0.1q1 / t16
      t27 = t13 + t12
      t28 = t15 * t16
      t27 = t28 * (t17 * t27 + t18 * t27 + t19 * t27)
      t29 = 2
      t17 = t29 * (t19 + t17 + t18)
      t30 = t26 * t16 * t14 / 2
      t31 = t17 + t30 - t27
      t20 = t29 * t20
      t32 = -t20 + t26
      t24 = 0.1q1 / t24
      t22 = t22 / 2
      t33 = t13 / 4
      t34 = -t33 * t32 * t16 + t22 * t31 * t24
      t35 = t21 * t31 * t24 / 2
      t36 = t12 / 4
      t37 = t10 * (-t35 + t3) + t36 * t32 * t16
      t17 = t17 - t30 - t27
      t20 = -t20 - t26
      t22 = -t33 * t20 * t16 + t22 * t17 * t24
      t26 = t21 * t17
      t27 = t26 * t24 / 2
      t30 = t10 * (-t27 + t3) + t36 * t20 * t16
      t33 = t20 * t16 / 4
      t27 = -t33 * t9 * t3 + t27 * t5
      t33 = t17 * t24 * t10 * t6 / 2 - t33 * t1 * t8
      t36 = t32 * t16 / 4
      t35 = -t36 * t9 * t3 + t35 * t5
      t36 = t31 * t24 * t10 * t6 / 2 - t36 * t1 * t8
      t38 = t8 * t21
      t39 = t3 * t6
      t40 = -t38 + t39
      t41 = t20 * t17
      t42 = t32 * t31
      t5 = 0.1q1 / t5
      t37 = 0.1q1 / t37
      t43 = 0.1q1 / t6
      t30 = 0.1q1 / t30
      t4 = 0.1q1 / t4
      t44 = t41 + t42
      t45 = t31 ** 2
      t46 = t17 ** 2
      t47 = t46 + t45
      t48 = t24 ** 2
      t49 = t8 ** 2
      t50 = t31 * t45 * t37
      t51 = t17 * t46 * t30
      t52 = t35 * t37
      t53 = t27 * t30
      t54 = t16 * t49
      t55 = t54 * t5 * t43
      t56 = t15 * t4 * t24
      t57 = 0.1q1 / t9
      t58 = t31 * t37
      t59 = t17 * t30
      t60 = t59 + t58
      t61 = t45 * t37
      t62 = t30 * t46
      t63 = t62 + t61
      t64 = t22 * t17
      t65 = t31 * t34
      t66 = t65 + t64
      t67 = t36 * t43
      t68 = t67 + t1
      t69 = t34 ** 2
      t70 = t2 ** 2
      t71 = t4 ** 2
      t72 = t35 * t5
      t73 = t3 * t7
      t74 = t4 * t69
      t75 = t22 * t5
      t76 = t2 * t11
      t77 = t2 * t30
      t78 = t1 * t21
      t79 = t22 * t4
      t80 = t45 * t34
      t81 = t80 * t78
      t82 = t46 * t22
      t83 = t25 * t43
      t84 = t83 * t8
      t85 = t1 * t27
      t86 = t85 * t5
      t87 = t76 * t4
      t88 = t24 * t7
      t89 = t27 * t5
      t90 = t89 - t3
      t91 = t72 - t3
      t92 = t31 + t17
      t93 = t22 ** 2
      t94 = t7 ** 2
      t95 = t12 * t57
      t96 = t52 * t31
      t97 = t53 * t17
      t98 = t17 * t27
      t99 = t31 * t35
      t100 = t1 * t4
      t101 = t7 * t5
      t102 = t53 * t2 * t17
      t103 = t96 * t34
      t104 = t18 * t4
      t105 = t7 * t8
      t106 = t11 * (t4 * t94 * (t52 * t45 * t24 * t57 + t96 * t5 + t97 *
     # (t17 * t24 * t57 + t5)) * t8 + t7 * (t24 * (-t25 * t4 * t63 + (t6
     #1 * t91 + t62 * t90) * t57 * t7) + t95 * t4 * t92 * t43) * t2 + t1
     #00 * t83 * t60 * t70 - t101 * t43 * t71 * (t98 * t33 + t99 * t36))
     # - t105 * t79 * t53 * t17 * t5 + t105 * (t5 * (t83 * t4 * t92 * t2
     #1 - t103 * t4 + t102) - t104 * t60 * t3) + t73 * t59 * t93 * t4 * 
     #t5
      t107 = t2 * t3
      t108 = -t107 * (t59 + t58) + (t96 + t97) * t4 * t8
      t109 = t2 * t37
      t110 = t59 * t22
      t111 = t17 * t33
      t26 = t24 * (t108 * t8 + t3 * t71 * (t31 * t36 + t111) * t43 * t11
     # + t43 * (t5 * (t21 * (t3 * (t31 * (t34 * (-t109 - t4) + t74 * t37
     #) + t110 * (t79 - t2)) + t8 * (t98 * (t77 + t4) + t96 * t2)) + t87
     # * (t97 * t22 + t103)) - t96 * t87 * t8 + t104 * t3 * (t58 * (t76 
     #- t38) - t26 * t8 * t30)) * t57 * t1) + t87 * t81 * t37 * t43 * t5
     #7 * t48
      t96 = t2 * t34
      t69 = t7 * (t5 * (-t110 * t107 + t58 * (t3 * (-t96 + t74) + t2 * t
     #8 * t35) + t4 * (t4 * (t17 * t93 + t31 * t69) + t102 * t33) * t43 
     #* t11) + t43 * (t38 * (t31 * (-t107 * t37 + (t37 * (-t34 * t5 + t8
     #) + t5) * t4 * t35) + t59 * (t27 * t4 * (-t75 + t8) - t107)) + t76
     # * (t4 * (t5 * (-t99 - t98) + t59 * t18 * t3) - t102 * t5)) * t57 
     #* t1) + t76 * t83 * t1 * t71 * t92 + t3 * t11 * t5 * t43 * t57 * (
     #t59 * t93 * t4 - t96 * t58) * t94
      t18 = t24 * t69 + t7 * t26 + t24 * t106 + t88 * (t57 * (t43 * (t31
     # * (-t38 * t12 * t4 + (t74 * t73 * t5 + (-t72 + t3) * t70 * t1) * 
     #t37 * t11) + t79 * t78 * t11 * t24 * (t77 + t4) * t46 + t1 * (-t76
     # * t53 * t8 * t4 - t3 * t21 * t4 * (t75 + t8) + t70 * t3 * t11 * t
     #30) * t17 + t81 * t11 * t24 * t71) - t56 * t18 * t3 * t63 - t56 * 
     #t5 * (t80 * t52 + t82 * t53)) + t19 * t11 * t5 * t43 * t71 * t66 -
     # t39 * t25 * t5 * t4 * t60 * t9 + t87 * (t58 * (t72 * t68 + t84) +
     # t59 * (t86 + t84))) + t84 * t1 * t71 * t5 * (t27 + t35) - t75 * t
     #76 * t59 * t24 * t94 * t3 * t57 * t43
      t26 = (t37 + t30) * t8
      t27 = t12 * t43
      t39 = t3 * t5
      t59 = t34 * t43
      t63 = t87 * t9
      t69 = t24 * t11
      t70 = t4 * t24
      t74 = t9 ** 2
      t75 = t2 * (-t83 + t9)
      t49 = t7 * t49
      t77 = t49 * t21 * t43
      t79 = t1 * t2
      t15 = t15 * t24
      t67 = t70 * (t17 * (t105 * t90 - t77 + t5 * (t7 * (t85 * t4 - t73 
     #+ t75) + t2 * t74 * t10 * t4) * t11) + t31 * (t11 * (t67 * t2 * t9
     # * t4 + t75 * t101 - t39 * t94) - t77 + t105 * t91) - t51 * t48 * 
     #t94 * t22 * t21 * t11 * t57 + t15 * (t21 * (-t105 * t37 - t100) + 
     #t34 * (t109 + t4)) * t45 + t15 * (t4 * (-t78 + t22) + t21 * (t105 
     #* t22 * t43 * t57 - t79) * t30) * t46)
      t63 = t4 * (-t73 * t55 * t25 * t9 * (t20 * t30 + t32 * t37) + t69 
     #* (t4 * (t2 * (t74 * t10 * t31 * t5 - t43 * t66 * t1 + t111 * t9 *
     # t43) + t101 * t9 * t66) + t105 * t2 * t43 * t92)) + t67 + t70 * (
     #t7 * (-t39 * t66 + t11 * (t1 * (t17 * t3 + t31 * (t72 + t3)) - t19
     # * t5 * t92 * t9 - t23 * t5 * (t82 + t80) * t24 * t43) * t4 + t76 
     #* t24 * (-t61 * t21 * t68 + t62 * (-t21 * t33 * t43 + t22))) - t63
     # * t1 * t92 - t63 * t83 * t10 * t5 * t92 + t69 * (-t50 * t21 * t34
     # * t24 * t57 + t62 * (t76 * t22 * t43 * t57 - t89 * t6 - t38) + t6
     #1 * (-t72 * t6 + t59 * t57 * (t38 + t76))) * t94)
      t19 = t19 * t9
      t31 = -t31 * t42 * t16 * t24 * t7 * t40
      t40 = t17 * t41 * t16 * t24 * t7 * t40
      t11 = t70 * t43 * t11
      t67 = t52 + t53
      t68 = t105 + t79
      t69 = t37 + t30
      t70 = t34 + t22
      t74 = t8 * t4
      t75 = t25 * t5
      t76 = t7 * t69
      t3 = t56 * (t64 * (-t86 * t43 * t4 + (-t88 * t17 + t43 * t68) * t3
     #0 * t3) + t59 * t58 * t3 * t68 - t80 * t88 * t3 * t37) * t57 + t75
     # * (t43 * (t4 * (t3 * (t2 * (-t30 * t33 - t36 * t37) + t4 * (t33 +
     # t36)) - t4 * t70 * t8) * t9 + t74 * t73 + t79 * (t107 * t69 - t74
     # * t67)) + t88 * t108 - t9 * t71 * t14) + t75 * (-t49 * t43 * t4 *
     # t67 - t12 * t2 * t9 * t4 * t69 + t3 * (t9 * (-t10 * t43 * t71 - t
     #76 * t4) + t76 * t2 * t43) * t8) + t15 * t71 * t43 * (t65 * (-t72 
     #* t1 * t57 + t8) + t64 * t8)
      t14 = t40 - t31
      ret = t29 * t18 - 8 * t4 * (t25 * (t8 * (t13 * t4 * t43 + t73 * (-
     #t24 * t60 + t26 * t43)) + t27 * (t26 + t5) * t2) + t95 * t56 * t43
     # * t66) + t11 * (t8 * (t39 * t44 + t38 * (t62 * t20 + t61 * t32) *
     # t24) * t16 * t94 + t100 * (-t40 + t31) + t101 * t4 * (t42 * (t8 *
     # (t1 * t35 + t19) - t12 * t34) + t41 * (t8 * (t85 + t19) - t12 * t
     #22)) * t16) / 4 - (0.3q1 / 0.4q1) * t27 * t28 * t24 * t8 * t71 * (
     #t41 + t42) + t11 * t5 * (t14 * t4 * t10 * t9 + t14 * t7) / 8 + 6 *
     # t83 * t12 * t71 * t5 * t70 - t56 * (t7 * (t21 * (-t5 * t47 * t24 
     #- t6 * (t51 + t50) * t48) + t55 * (t53 * t41 + t52 * t42)) + (t24 
     #* (t47 * t6 + t23 * (-t46 - t45) * t5) + t54 * t44 * t43) * t4 * t
     #9) / 2 + t63 - 4 * t3 - 12 * t84 * t12 * t71

      hjetmass_box_pppm_0_0_s12_mhsq_s34_s123 = ret/32q0/(0,1q0)
      return

      end function
