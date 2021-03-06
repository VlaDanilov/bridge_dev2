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
 

      double complex function hjetmass_triangle_pppm_s123_0_mhsq_dp 
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
      t2 = za(i1, i3)
      t3 = zb(i2, i1)
      t4 = zb(i3, i1)
      t5 = za(i2, i4)
      t6 = za(i3, i4)
      t7 = t5 * t3
      t8 = t6 * t4
      t9 = t8 + t7
      t10 = zb(i4, i2)
      t11 = za(i1, i4)
      t12 = za(i2, i3)
      t13 = zb(i4, i1)
      t14 = zb(i3, i2)
      t15 = t1 * t3
      t16 = t2 * t4
      t17 = t12 * t14
      t18 = t17 + t15 + t16
      t19 = zb(i4, i3)
      t20 = t11 * t13
      t21 = t5 * t10
      t22 = t6 * t19
      t23 = t20 + t21 + t22
      if ( dreal(t23) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t22 = cg * cdsqrt(t23 ** 2) + t20 + t21 + t22
      t23 = 0.1D1 / t22
      t24 = -2 * t20 * t18 * t23 + t15 + t16
      t25 = t11 * (t10 * (-2 * t5 * t18 * t13 * t23 + t12 * t4) + t13 * 
     #t24)
      t26 = t11 * t3
      t27 = -t14 * t6 + t26
      t28 = 0.1D1 / t13
      t29 = 0.1D1 / t5
      t24 = 0.1D1 / t24
      t25 = 0.1D1 / t25
      t30 = t16 + t15
      t31 = t30 * t28
      t32 = t31 * t24 + (t17 + t15 + t16) * t25 * t11
      t33 = t9 ** 2
      t34 = t9 * t33
      t35 = t19 * t29
      t36 = t10 * t28
      t37 = t36 * t35
      t38 = t37 * t25
      t39 = 0.1D1 / t12
      t40 = 0.1D1 / t19
      t41 = 0.1D1 / t1
      t42 = 0.1D1 / t6
      t43 = t13 * t14
      t44 = t3 * t19
      t45 = t44 + t43
      t46 = t36 * t4
      t47 = -t46 + t14
      t46 = t46 - t14
      t48 = t4 ** 2
      t49 = t24 ** 2
      t50 = t24 * t49
      t51 = t19 ** 2
      t52 = t3 ** 2
      t53 = t28 ** 2
      t54 = t25 ** 2
      t55 = t25 * t54
      t56 = t14 ** 2
      t57 = t16 * t3
      t58 = t9 * t24
      t59 = t5 * t41
      t60 = t10 * t48
      t61 = t1 * t52
      t62 = t61 * t51 * t28
      t63 = t60 * t2
      t64 = t6 * t29
      t65 = t28 * t24
      t66 = t11 * t29
      t67 = t11 ** 2
      t68 = t11 * t67
      t69 = t1 * t10
      t70 = t3 * t28
      t71 = t70 * t51
      t72 = t4 * t10
      t73 = t72 * t11
      t74 = t17 * t11
      t75 = t14 * t40
      t76 = t25 * t24
      t77 = t76 * t29
      t45 = t9 * (t4 * (-t67 * t4 * t41 * t29 * t25 + t24 * t39 * (t75 -
     # t70)) + t77 * (t65 * ((t57 * t19 + t61 * t19 - t63) * t42 * t11 +
     # t62) + t74 * t25 * (-t73 * t42 + t71) + t26 * t4 * (-t69 * t11 * 
     #t42 + t2 * t51 * t28) * t25) * t33 + t7 * t39 * t40 * t24 * t42 * 
     #t47 + t58 * (t25 * (t66 * t36 * t48 + t56) + (t15 * t40 * t46 - t1
     #6 * (t75 + t70)) * t42 * t24 * t39)) + t9 * (t25 * (t26 * t58 * t4
     #2 * t46 - t59 * t56 + t29 * (t57 * t51 * t53 + (t16 * t14 + t15 * 
     #t47) * t42 * t11) * t49 * t33) + t65 * (-t48 * (t10 * t40 + t64) +
     # t58 * (t63 * t40 - t61) * t42) * t39 + t66 * t24 * t33 * ((t2 * (
     #t4 * t45 - t60) + t17 * t45 + t15 * t45) * t42 * t11 + t62) * t54)
      t22 = 0.1D1 / t22
      t46 = 0.1D1 / t11
      t47 = t10 * t42
      t60 = t35 + t47
      t62 = t20 * t29
      t63 = t62 + t10
      t70 = mt ** 2
      t75 = t47 * t35
      t78 = t42 * t54
      t79 = t42 * t55
      t80 = t79 * t19 * t67
      t81 = t58 * t10
      t82 = t81 * t29
      t83 = t70 * t23
      t84 = t17 * t3
      t85 = t78 * t19 * t67
      t86 = t47 * t19
      t87 = t75 * t33
      t79 = t83 * (t11 * (-t38 * t33 * t49 * t42 + t86 * t27 * (-t59 + t
     #58) * t54) - t87 * t67 * t24 * t54 + t78 * t35 * t9 * t68 * t13 * 
     #t41 + t9 * t28 * t39 * t49 * t60) + (t10 * (t7 * t42 + t4) + t44) 
     #* t49 * t39 * t22 * t28 * t18 + t84 * t79 * t35 * t68 * t9 * t13
      t7 = t33 * t79 + t33 * (t19 * (t11 * (-t36 * t3 * t9 * t49 * t42 *
     # t25 + t10 * (t47 * t31 * t49 * t33 - t7 * t58 * t36 * t42 + t59 *
     # t14) * t54) + t4 * t68 * t13 * t29 * t41 * t54 + t8 * t28 * t39 *
     # t29 * t49 + (t41 * (t72 + t43) + t47 * t17 * t33 * t29 * t49) * t
     #54 * t67) * t22 * t18 + t85 * (-t77 * t15 * t17 * t33 * t10 - t70 
     #* t27 * t13 * t23 * t41 + (t84 * t25 + t83 * t41) * t10 * t9) + t1
     #6 * t9 * (t3 * (t1 * (-t78 * t37 * t9 * t11 * t49 - t47 * t58 * t3
     #5 * t67 * t55 + t53 * (-t75 * t9 * t25 + t60 * t46 * t39) * t50) +
     # t80 * t63) + t80 * t17 * (t41 * t63 - t82)))
      t23 = t41 * (t17 + t16) + t3 + t81 * t18 * t22
      t27 = -t17 - t15 - t16
      t37 = t13 ** 2
      t59 = t10 ** 2
      t63 = t10 * t9
      t70 = t22 * t18
      t75 = t42 * t25
      t77 = t70 * t11
      t16 = t70 * t24 * t34 * (t70 * t80 * t5 * t9 * t10 * t59 + t35 * t
     #49 * t39 * (-t77 + t31) + t47 * t49 * (t70 * t35 * t9 * t67 * t25 
     #+ t31 * t39 - t77 * t39) + t70 * t85 * t58 * t59) + t75 * t22 * t1
     #8 * t19 * t11 * t34 * (t11 * (-t63 * t29 * t49 * t30 * t25 + t59 *
     # (t5 * (t41 * (-t21 * t18 * t22 + t16 + t17) + t3) + t58 * t27) * 
     #t54) + t67 * (t63 * t13 * t18 * t22 * t29 * t49 * t25 + t82 * t13 
     #* t27 * t54) + t29 * t37 * t23 * t54 * t68 - t63 * t31 * t29 * t50
     # - t70 * t67 ** 2 * t13 * t37 * t41 * t29 * t54)
      t31 = t72 - t43
      t37 = t35 * t67
      t58 = t42 * t11
      t59 = t4 * t14 * t41
      t63 = t9 * t11
      t78 = t2 ** 2 * t48
      t79 = t12 ** 2 * t56
      t80 = t79 + t78
      t46 = t46 * t39
      t52 = t1 ** 2 * t52
      t81 = t54 * t19 * t11
      t82 = t52 + t78
      t1 = t33 * (t2 * (t48 * (t37 * t41 * t54 + t46 * (t35 * t6 + t10) 
     #* t49 * t53) + t81 * t59) + t75 * t19 * (-t13 * t29 * t41 * t80 * 
     #t54 * t68 - t57 * t10 * t53 * t49 - t61 * t10 * t67 * t54 + t76 * 
     #t26 * t36 * t27) * t9 + t87 * t76 * (t65 * t82 * t25 * t11 + t49 *
     # t53 * t82 + (t79 + t78 + t52) * t54 * t67) + t44 * (t46 * t15 * t
     #53 * t49 + t11 * t14 * t54)) + t33 * (t28 * (t87 * t74 * t54 * t30
     # * t49 + t76 * t14 * t19 * (t8 * t29 + t3)) + t53 * (t3 * (-t86 * 
     #t15 * t9 * t25 + t46 * (t4 * (t19 * (t64 * t1 + t2) + t69) + t47 *
     # t30 * t5)) * t49 + t46 * t9 * (-t78 * t60 + t52 * (-t35 - t47)) *
     # t50) + t81 * (t12 * t56 * t41 + t66 * (t17 * t41 + t3) * t4 + t63
     # * t75 * (-t10 * t41 * t80 - t62 * t61)))
      ret = 40 * t76 * t35 * t14 * t34 * t32 + 24 * t38 * t24 * t34 * t4
     # * t32 + 8 * t45 + 64 * t7 - 48 * t76 * t73 * t22 * t18 * t19 * t3
     #4 * (t25 * (t66 + t36) + t65 * t29) - 128 * t16 + 16 * t63 * (-t59
     # * t25 + t70 * t42 * t39 * (t40 * (-t72 + t43) + t3) * t49 * t9 + 
     #t70 * (t54 * (t3 * (t51 * (-t66 - t36) - t13 * t37 * t42) + t58 * 
     #(t10 * t31 + t62 * t31)) * t24 + t25 * t29 * (t58 * (t72 - t43 - t
     #44) - t71) * t49) * t33) - 32 * t1 - 256 * t70 * t86 * t55 * t13 *
     # t68 * t34 * t23 + 384 * t86 * t18 ** 2 * t22 ** 2 * t68 * t34 * t
     #13 * t41 * t55 * (t20 + t21) - 80 * t77 * t76 * t19 * t34 * (t14 *
     # (t20 * t25 + t24) * t29 + t25 * (t26 * t42 + t14) * t10)

      hjetmass_triangle_pppm_s123_0_mhsq_dp = ret/32d0/(0,1d0)
      return

      end function
