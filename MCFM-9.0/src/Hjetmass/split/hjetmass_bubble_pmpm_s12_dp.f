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
 

      double complex function hjetmass_bubble_pmpm_s12_dp 
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double complex ret
          double precision mt

      t1 = za(i1, i4)
      t2 = za(i2, i3)
      t3 = zb(i3, i1)
      t4 = za(i1, i3)
      t5 = za(i3, i4)
      t6 = zb(i3, i2)
      t7 = 0.1D1 / t6
      t8 = 0.1D1 / t4
      t9 = t2 * t8
      t10 = t3 * t7
      t11 = t9 + t10
      t12 = zb(i4, i1)
      t13 = zb(i4, i2)
      t14 = 0.1D1 / t13
      t15 = t12 * t14
      t16 = t9 + t15
      t17 = zb(i4, i3)
      t18 = za(i2, i4)
      t19 = t1 * t2
      t20 = t19 * t8
      t21 = t18 - t20
      t22 = t2 ** 2
      t23 = t8 ** 2
      t24 = t8 * t23
      t25 = t18 * t3
      t26 = t1 * t22
      t27 = 0.2D1
      t28 = 0.1D1 / t1
      t29 = t18 * t28 + t15
      t30 = t12 ** 2
      t31 = t14 ** 2
      t32 = t1 * t6
      t33 = t1 * t3
      t34 = -t15 * t6 + t3
      t35 = t15 * t1 + t18
      t36 = t2 * t6
      t37 = t36 * t8 + t3
      t38 = t2 * t3
      t39 = t18 * t12
      t40 = t38 + t39
      t41 = t4 * t6
      t42 = t1 * t13 + t41
      t43 = t4 * t3
      t44 = t1 * t12
      t45 = t18 * t13
      t46 = t43 - t36 + t44 - t45
      t46 = t46 ** 2
      t47 = 0.4D1 * t40 * t42 + t46
      t47 = cdsqrt(t47)
      t48 = t43 - t36 + t44 - t45 + t47
      t49 = 0.1D1 / t42
      t50 = 0.1D1 / 0.2D1
      t51 = t50 * t1
      t52 = -t51 * t48 * t49 - t18
      t53 = 0.1D1 / t42
      t54 = t5 * t17
      t55 = t3 * t13
      t56 = t6 * t12
      t57 = t56 + t55
      t58 = t48 ** 2
      t59 = t53 ** 2
      t60 = t59 ** 2
      t61 = t53 * t59
      t62 = t18 * t57
      t63 = t62 * t59
      t64 = t43 * t6 * t61
      t65 = t25 * t49 * (t54 + t43)
      t66 = t44 + t43
      t67 = t45 + t36
      t68 = t1 ** 2
      t69 = t6 * t18
      t70 = t69 * t67
      t71 = t58 * t59
      t72 = t4 * t18
      t73 = t54 * t1
      t74 = t6 * (t72 - t19) + t73
      t75 = t61 * t48 * t58
      t76 = (-t39 - t38) * t18
      t77 = 0.3D1 / 0.4D1
      t78 = 0.1D1 / 0.4D1
      t79 = -0.1D1 / 0.8D1
      t80 = 0.1D1 / 0.16D2
      t81 = 0.3D1 / 0.2D1
      t82 = t39 * t27 * t53 + t38 * t81 * t49
      t57 = t50 * t48 * (t1 * (-t63 * t48 - t64 * t58) + t65) + t77 * t3
     #3 * t58 * t59 * t66 - t78 * t71 * (t57 * t53 * t48 * t68 + t70 + t
     #54 * (-t33 + t69)) + t79 * t75 * t6 * t74 + t80 * t32 * t58 ** 2 *
     # t60 * t42 - t3 * (t71 * t19 * t6 + t76) + t33 * t48 * t82
      t46 = 0.4D1 * t40 * t42 + t46
      t46 = cdsqrt(t46)
      t83 = t43 - t36 + t44 - t45 + t46
      t84 = t50 * t48 * t49
      t85 = -t9 - t84
      t83 = 0.1D1 / t83
      t86 = t27 * t40
      t87 = -t86 * t83 - t84
      t47 = t43 - t36 + t44 - t45 - t47
      t88 = t49 * (t47 - t48)
      t89 = -t50 * t6 * t47 * t49 + t3
      t90 = t67 * t53
      t91 = t47 ** 2
      t92 = t47 * t91
      t93 = t1 * t91 * t59
      t94 = -t43 + t54
      t95 = (t54 + t36 + t45) * t18
      t96 = t68 * t12
      t97 = t94 * t53
      t46 = t43 - t36 + t44 - t45 - t46
      t98 = t50 * t47 * t49
      t99 = t15 - t98
      t46 = 0.1D1 / t46
      t86 = -t86 * t46 - t98
      t45 = t49 * t18 * (t54 + t36 + t45)
      t67 = t67 * t59
      t100 = t18 * t40
      t101 = t40 * t53
      t98 = t9 + t98
      t102 = -t15 + t84
      t84 = -t84 * t6 + t3
      t49 = -t51 * t47 * t49 - t18
      t51 = t56 + t55
      t103 = t1 * t18
      t104 = t18 ** 2
      t105 = t104 * t13
      t106 = t54 * t33
      t107 = t2 * t6 ** 2 * t18
      t108 = (-t54 * t18 - t105) * t6
      t109 = t25 * (t54 + t43)
      t110 = t91 * t59
      t111 = t110 * t6
      t55 = -t56 - t55
      t56 = t50 * t47 * (t1 * (-t63 * t47 - t64 * t91) + t65) + t77 * t1
     #10 * t33 * t66 + t78 * t110 * (t55 * t53 * t47 * t68 - t70 + t54 *
     # (t33 - t69)) + t79 * t92 * t61 * t6 * t74 + t80 * t32 * t91 ** 2 
     #* t60 * t42 - t3 * (t111 * t19 + t76) + t33 * t47 * t82
      t60 = t71 * t1
      t63 = t101 * t48 * t1 - t50 * t48 * (t67 * t48 * t1 + t45) + t78 *
     # t60 * (-t54 + t43 + t44) + t79 * t75 * t1 * t42 + t100
      t64 = 0.1D1 / t98
      t11 = 0.1D1 / t11
      t17 = 0.1D1 / t17
      t5 = 0.1D1 / t5
      t65 = 0.1D1 / t16
      t16 = 0.1D1 / t16
      t66 = -0.1D1 / t85
      t70 = 0.1D1 / t102
      t29 = 0.1D1 / t29
      t76 = -0.1D1 / t99
      t80 = t33 * t7 + t18
      t82 = t34 ** 2
      t112 = t21 * t53 * t64 * t66
      t65 = t65 * t14
      t113 = t21 * t3
      t114 = t18 * t37
      t115 = t39 * t31 * t29
      t5 = t5 * t17
      t17 = t38 * t23
      t4 = t5 * (t17 * (t21 * (-t27 * t26 * t23 * (t20 * t55 + t62) + t5
     #4 * t32 * t2 * t22 * t24 + (-t105 * t6 + t54 * (t33 - t69)) * t23 
     #* t22 + t104 * t3 * t12 - t9 * t54 * t25 + t68 * t22 ** 2 * t6 * t
     #13 * t23 ** 2 + t20 * t3 * t12 * (-0.4D1 * t18 + 0.3D1 * t20)) * t
     #59 * t64 ** 2 * t66 ** 2 + t65 * (t20 * t27 * t3 + t26 * t6 * t23 
     #- t25)) * t11 * t7 + t115 * t28 * (t82 * (t27 * t15 * t19 * t34 - 
     #t32 * t4 * t12 * t30 * t14 * t31 - t94 * t1 * t31 * t30 + t25 * t2
     # - t15 * (t36 + t54) * t18) * t59 * t70 ** 2 * t76 ** 2 + (-t15 * 
     #t33 * t27 + t32 * t30 * t31 - t25) * t8 * t16))
      t9 = 0.1D1 / t88
      t12 = 0.1D1 / t85
      t13 = -0.1D1 / t102
      t20 = 0.1D1 / t99
      t25 = -0.1D1 / t98
      t26 = 0.1D1 / t86
      t30 = 0.1D1 / t87
      t31 = t9 ** 2
      t55 = t49 * t56
      t62 = t55 * t8 * t25
      t45 = (-t50 * t47 * (t67 * t47 * t1 + t45) - t78 * t93 * (t54 - t4
     #3 - t44) + t79 * t1 * t92 * t61 * t42 + t101 * t47 * t1 + t100) * 
     #t20
      t54 = t52 * t57 * t12 * t8
      t67 = t30 * t83
      t78 = t5 * t31
      t79 = t33 + t69
      t85 = t5 * t8
      t86 = t52 ** 2
      t87 = t47 * t49
      t88 = t85 * t9 * t59 * t40 * (t48 * (t67 * t18 * t84 * t12 * t52 +
     # t67 * t3 * t12 * t86) - t87 * t46 * t25 * t26 * (t18 * t89 + t3 *
     # t49))
      t94 = t89 ** 2
      t98 = t84 ** 2
      t99 = t45 * t94 * t14
      t100 = t47 * t46
      t31 = t5 * t9 * t31 * t61 * t40 * (-t100 * t26 * (t62 + t99) + t67
     # * (t98 * t63 * t13 * t14 + t54) * t48)
      t101 = t26 ** 2
      t102 = t12 * t8
      t19 = t67 * t48 * (t98 * (-t63 * t14 * t13 ** 2 + (-t63 * t30 + t9
     #7 * t48 * t1 - t96 * t48 * t53 + t27 * t1 * (t90 * t48 - t38 - t39
     #) + t77 * t60 * t42 + t95) * t14 * t13) + t102 * (t52 * (t27 * t10
     #3 * t48 * t53 * t51 - t50 * t75 * t32 * t42 - t77 * t71 * t6 * (t6
     # * (-t72 + t19) - t73) + t81 * t71 * t68 * t51 - 0.3D1 * t33 * (-t
     #41 * t71 + t38 + (t44 + t43) * t53 * t48) - t109 + (-t106 + t107 -
     # t108) * t53 * t48 + 0.4D1 * t33 * (t36 * t48 * t53 - t39)) + t57 
     #* (-t52 * (t12 + t30) + t1)))
      t1 = t78 * t61 * t40 * (t19 + t100 * (t101 * (-t62 - t99) + t26 * 
     #(t8 * (t25 * (t1 * t56 + t49 * (t27 * t103 * t47 * t53 * t51 - t50
     # * t32 * t92 * t61 * t42 + t77 * t111 * t74 + t81 * t110 * t68 * t
     #51 + 0.3D1 * t33 * (t41 * t91 * t59 - t38 + (-t44 - t43) * t53 * t
     #47) - t109 - (t106 - t107 + t108) * t53 * t47 + 0.4D1 * t33 * (t36
     # * t47 * t53 - t39))) - t55 * t25 ** 2) + t94 * t20 * t14 * (-t45 
     #+ t27 * t1 * (t90 * t47 - t38 - t39) + t77 * t93 * t42 + t97 * t47
     # * t1 - t96 * t47 * t53 + t95))))
      t19 = t30 ** 2
      t27 = t83 ** 2
      t19 = t5 * t9 * t59 * t40 ** 2 * (t48 * (t98 * t27 * t13 * t14 * t
     #19 * t52 + t102 * t84 * t27 * t19 * t86) - t87 * t89 * t46 ** 2 * 
     #t101 * (t89 * t14 * t20 + t49 * t8 * t25))
      ret = -0.8D1 * t5 * (t2 * (t3 ** 2 * t21 * t37 * t7 ** 2 * t23 * (
     #t65 + t112) * t11 ** 2 + t65 * t10 * t23 * (t114 - t113) * t11) + 
     #t65 * t3 * t24 * t80 * t11 * t22 + t115 * (t18 * t35 * t82 * t28 *
     #* 2 * t29 * t53 * t70 * t76 + (t28 * (-t3 * t35 + (t35 * t28 * t29
     # + 0.1D1) * t34 * t18) + t15 * (-t69 * t28 - t3)) * t8 * t16)) + 0
     #.16D2 * t4 - 0.128D3 * t78 * t59 * t40 * ((t45 * t89 * t14 * (t6 *
     # t47 * t53 - t89) - t62) * t26 * t46 + t67 * (t84 * t63 * t13 * t1
     #4 * (t6 * t48 * t53 - t84) - t54)) - 0.12D2 * t5 * t17 * t112 * t1
     #1 * t7 * (t114 - t113) - 0.20D2 * t85 * t53 * (t113 * t22 * t23 * 
     #t11 * t66 * t64 * t80 + (-t91 * t49 * t25 * t46 * t26 * t79 + t67 
     #* t79 * t12 * t52 * t58) * t9 * t59 * t40) + 0.24D2 * t88 + 0.256D
     #3 * t31 - 0.64D2 * t1 + 0.32D2 * t19

      hjetmass_bubble_pmpm_s12_dp = ret/16d0*(0,1d0)
      return

      end function
