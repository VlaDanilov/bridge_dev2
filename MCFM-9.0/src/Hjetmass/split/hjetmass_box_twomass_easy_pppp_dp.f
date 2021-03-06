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
 

      double complex function hjetmass_box_twomass_easy_pppp_dp 
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
      t3 = zb(i3, i2)
      t4 = zb(i4, i1)
      t5 = za(i1, i2)
      t6 = za(i2, i4)
      t7 = zb(i2, i1)
      t8 = zb(i4, i2)
      t9 = za(i1, i3)
      t10 = zb(i3, i1)
      t11 = za(i3, i4)
      t12 = zb(i4, i3)
      t13 = t5 * t8
      t14 = t9 * t12
      t15 = t14 + t13
      t16 = t4 ** 2
      t17 = t1 ** 2
      t18 = t17 * t16
      t19 = t18 * t15
      t20 = t6 * t8
      t21 = t11 * t12
      t22 = t1 * t4
      t23 = t22 * t2 * t3
      t24 = t9 * t10
      t25 = t22 * ((-t21 - t20) * t7 * t5 + t23 - t24 * (t21 + t20))
      t11 = t11 * t10 + t7 * t6
      t26 = mt ** 2
      t27 = 16 * t22 * t26 * t11 * t19 + 4 * t25 ** 2
      t27 = cdsqrt(t27)
      t28 = t21 + t20
      t29 = -t28 * t7 * t5 - t24 * t28 + t23
      t19 = 0.1D1 / t19
      t30 = t14 + t13
      t31 = -t5 * t30 * t7 - t24 * t30
      t15 = t22 * t27 * t19 * t15 / 4
      t20 = t18 * t19 * (t20 * t31 + t21 * t31 + t23 * t30) / 2
      t21 = -t29 + t15 + t20
      t11 = t22 * t11
      t15 = -t29 - t15 + t20
      t20 = -2 * t25
      t23 = t20 + t27
      t11 = 0.1D1 / t11
      t25 = 0.1D1 / t4
      t29 = 0.1D1 / t1
      t28 = t28 * t29
      t25 = t28 * t25 * t5
      t30 = t15 * t11
      t31 = t30 * t6
      t32 = t31 - t25
      t33 = t23 * t19 * t5 / 4
      t34 = -t10 * t32 - t33 * t12
      t32 = -t7 * t32 - t33 * t8
      t20 = t20 - t27
      t27 = t21 * t11
      t6 = t27 * t6
      t25 = t6 - t25
      t33 = t20 * t19 * t5 / 4
      t35 = -t10 * t25 - t33 * t12
      t25 = -t7 * t25 - t33 * t8
      t28 = t28 * t5
      t31 = t31 * t4 - t28
      t6 = t6 * t4 - t28
      t28 = 0.1D1 / t9
      t33 = 0.1D1 / t5
      t36 = 0.1D1 / t25
      t2 = 0.1D1 / t2
      t37 = 0.1D1 / t32
      t38 = t21 ** 2
      t39 = t15 ** 2
      t40 = t39 * t37
      t41 = t38 * t36
      t42 = t40 + t41
      t43 = -t21 - t15
      t44 = t36 + t37
      t45 = t5 ** 2
      t46 = t10 ** 2
      t47 = t11 ** 2
      t48 = t7 ** 2
      t49 = t26 * t33
      t50 = t26 * t7
      t51 = t24 * t33
      t52 = t3 * t33
      t53 = t10 * t2
      t54 = t3 * t28
      t55 = t28 * t7
      t56 = t50 * t3
      t57 = t21 * t36
      t58 = t15 * t37
      t59 = t58 + t57
      t60 = t26 * t4
      t61 = t4 * t2
      t62 = t3 * t26
      t63 = t10 * t7
      t64 = t5 * t48
      t65 = t3 * t4
      t66 = t8 * t10
      t67 = -t65 - t66
      t68 = t33 ** 2
      t69 = t66 * t59
      t70 = t59 * t12
      t71 = t24 * t12
      t72 = t60 * t12
      t73 = t55 * t33
      t38 = t10 * (t38 + t39)
      t74 = t7 * t31
      t75 = t3 * t36
      t76 = t3 * t37
      t77 = t3 * t7
      t78 = t8 * t3
      t79 = t11 * t10
      t80 = t11 * t1
      t81 = t26 * t28
      t82 = t31 + t6
      t83 = t58 * t34
      t84 = t34 * t29
      t85 = t7 * t36
      t86 = t14 * t7 * t37 * (-t84 * t3 + t30 * t46) + t85 * (t3 * t12 *
     # t6 + t66 * t27 * t35) + t80 * (t53 * t43 + t54 * t43) * t16 + t65
     # * (t26 * t12 * t36 + (t57 * t35 + t83) * t11 * t7)
      t87 = t4 * t28
      t88 = t7 * t6
      t89 = t26 * t10
      t90 = t4 * t34
      t26 = t78 * (t36 * (t28 * (t60 + t88) + t89 * t29) + t37 * (t26 * 
     #(t10 * t29 + t87) + t7 * (t28 * t31 - t84)))
      t36 = t54 * t16
      t84 = t22 * t33
      t91 = t33 * t10
      t84 = t7 * (-t84 * t40 * t46 * t11 + t57 * (t1 * (-t27 * t46 * t33
     # * t4 + t36) - t81 * t66) + t58 * (t36 * t1 + t66 * (t33 * t34 - t
     #81))) + t28 * (-t79 * t40 * t22 + t57 * (t10 * (-t22 * t27 + t13) 
     #+ t65 * t5) - t58 * t5 * t67) * t48 + t84 * (t91 * (t15 * t31 + t2
     #1 * t6) + t87 * t21 * t25) * t2 + t91 * t4 * (t43 * t28 * t8 * t1 
     #- t62 * t59)
      t92 = t14 * t52 * (-t85 * t35 + t89 * t44) * t29
      t16 = t11 * t84 + t33 * t86 + t4 * (t55 * (-t30 * t7 * t2 - t78 * 
     #t44) * t5 + t79 * (t43 * t68 * t12 * t1 + t51 * t43 * t2 + t77 * t
     #59)) + t11 * (t10 * (t7 * (t70 * t7 + t69) + t52 * t43 * t4 + t71 
     #* t43 * t68) + t73 * (t58 * (t31 * t67 - t72) + t57 * (t6 * t67 - 
     #t72)) * t1) + t35 * (-t75 * t7 * t8 * t29 - t14 * t61 * t68) + t22
     # * t68 * t10 * (t55 * (t40 * t31 + t41 * t6) * t1 + t38) * t47 + t
     #76 * t12 * (t33 * (t60 + t74) - t64 * t29) + t77 * (-t55 * t44 * t
     #8 * t45 - t66 * t44 * t5 - t71 * t44) * t29 + t80 * (t2 * (t68 * (
     #t15 * t34 + t21 * t35) + t55 * t43) + t38 * t28 * t68 * t11 * t1) 
     #* t16 - t81 * t70 * t11 * t48 + t26 - t64 * (t75 * t12 * t29 + t61
     # * t27 * t28) + t14 * t2 * (-t10 * t82 - t90) * t68 + t92
      t22 = t41 * t4
      t25 = t32 + t25
      t6 = t2 * (t33 * (t12 * (-t25 * t4 - t7 * t82) + t8 * (-t10 * t31 
     #- t4 * (t34 + t35))) - t87 * t8 * t25) + t8 * (t11 * t43 * t33 * t
     #46 - t55 * t2 * t82 - t53 * t6 * t33) + t16 + t51 * t11 * t12 * t7
     # * (t83 * t33 + t57 * (t33 * t35 + t10)) + t80 * (t69 * t55 * t4 +
     # t63 * (-t57 * t12 * t6 - t58 * (t12 * t31 + t90 * t30) - t22 * t3
     #5 * t11) * t68 + t4 * (t70 * t63 + (t15 * (t32 * t4 + t74) + t88 *
     # t21) * t2 * t28) * t33) - t73 * t18 * t10 * t47 * t42 - t36 * t80
     # * t49 * t59 - t77 * t4 * t12 * t44
      t16 = t28 * (t66 + t65) + t91 * t12
      t18 = -t75 + t2
      t25 = -t76 + t2
      t26 = t61 * (t91 + t55)
      t14 = t53 * (t14 * t33 + t8)
      t13 = t13 * t55
      t24 = t29 * t2 * (t12 * (t24 * (t49 - t7) + t50) + t8 * ((t81 - t1
     #0) * t7 * t5 + t89))
      t30 = 8
      ret = t30 * (t12 * (-t64 * t2 * t29 + t49 * (t63 * t11 * t59 + t61
     #)) + t8 * (t2 * (-t9 * t46 * t29 + t60 * t28) + t55 * (-t62 * t29 
     #* t44 - t61) * t5) + t29 * (-t8 * (t45 * t48 * t28 * t2 + t10 * t3
     #) - t56 * t44 * t12 - t51 * t3 * t12 - t9 ** 2 * t46 * t12 * t33 *
     # t2) + t4 * (t12 * (t2 * (-t7 - t51) - t52) - t8 * (t53 + t54) + t
     #2 * (t10 * (t49 * t43 - t43 * t7) + t50 * t43 * t28) * t11 - t55 *
     # t49 * t1 * t10 * t42 * t47)) + 12 * t87 * t56 * t11 * t59 - 4 * t
     #6 - 2 * t19 * t11 * t4 * (-t73 * t23 * t15 * t39 * t47 * t17 * t10
     # * t4 * t37 + t20 * t21 * (-t22 * t73 * t47 * t17 * t10 + t12 * t1
     #8 * t7 + t13 * t18 + t14 + t27 * (t85 * t16 - t26) * t1) + (t12 * 
     #t25 * t7 + t13 * t25 + t14) * t23 * t15 + t80 * (t7 * t16 * t37 - 
     #t26) * t23 * t39) + 16 * t24

      hjetmass_box_twomass_easy_pppp_dp = ret/32d0/(0,1d0)
      return

      end function
