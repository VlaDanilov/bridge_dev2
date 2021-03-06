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
 


      double complex function hjetmass_triangle_pppm_s34_s134_0_dp 
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

      t1 = za(i3, i4)
      t2 = zb(i4, i3)
      t3 = za(i1, i3)
      t4 = zb(i3, i1)
      t5 = za(i1, i4)
      t6 = zb(i4, i1)
      t7 = t3 * t4
      t8 = t5 * t6
      t9 = t8 + t7
      if ( dreal(t9) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t9 = cg * cdsqrt(t9 ** 2) + t7 + t8
      t10 = za(i1, i2)
      t11 = zb(i2, i1)
      t12 = za(i2, i3)
      t13 = za(i2, i4)
      t14 = t12 * t4 + t13 * t6
      t15 = zb(i3, i2)
      t16 = zb(i4, i2)
      t17 = t12 * t15
      t18 = t13 * t16
      t19 = 0.1D1 / t9
      t20 = t10 * t11
      t21 = t1 * t2
      t22 = -2 * t21 * t20 * t19 + t17 + t18
      t23 = -2 * t10 * t1
      t24 = t2 * (t23 * t4 * t19 + t13)
      t23 = t2 * (t23 * t6 * t19 - t12)
      t25 = t10 * t22
      t26 = t3 * t24
      t27 = (t23 * t5 + t25 + t26) * t11
      t28 = t1 * (-2 * t5 * t11 * t2 * t19 - t15)
      t29 = t11 * t22
      t30 = (t28 * t6 + t4 * t1 * (-2 * t3 * t11 * t2 * t19 + t16) + t29
     #) * t10
      t31 = 0.1D1 / t10
      t23 = 0.1D1 / t23
      t32 = 0.1D1 / t22
      t33 = 0.1D1 / t13
      t9 = 0.1D1 / t9
      t27 = 0.1D1 / t27
      t34 = t15 * t14
      t35 = t11 * t24
      t36 = t35 + t34
      t37 = t4 * t22
      t38 = t34 - t37
      t39 = -t8 - t7
      t40 = t8 + t7
      t41 = t4 ** 2
      t42 = t24 ** 2
      t43 = t11 ** 2
      t44 = t43 ** 2
      t45 = t11 * t43
      t46 = t27 ** 2
      t47 = t27 * t46
      t48 = t32 ** 2
      t49 = t4 * t33
      t50 = t49 * t5
      t51 = t14 * t5
      t52 = t23 * t6
      t53 = t52 * t22
      t54 = t53 * t24
      t55 = t24 * t14
      t56 = t24 * t22
      t57 = t27 * t42
      t58 = t11 * t42
      t59 = t58 * t31
      t60 = t31 * t32
      t61 = t41 * t22
      t62 = t9 * t23
      t63 = t62 * t5
      t64 = t14 ** 2
      t65 = t52 * t42
      t66 = t4 * t24
      t67 = t52 * t24
      t68 = t49 * t1
      t69 = t24 * t43
      t26 = t63 * (t31 * (-t36 * t4 + t65 * (-t68 + t11)) + t69 * (-t20 
     #* t5 * t64 * t15 * t33 + t35 * t51 * (-t8 - t20 - t7) * t33 + t37 
     #* (t22 * (t8 + t20 + t7) - t26 * t11)) * t46 + t11 * (t22 * (t38 *
     # t4 + t52 * (-t34 + t37) * t24) + t51 * (t11 * (t66 - t65) - t61 +
     # t34 * (-t67 + t4)) * t33) * t27) + t63 * (t1 * (t60 * t52 * t42 *
     # t33 * t36 + t59 * t33 * t36 * t48) + t31 * (t52 * t38 * t24 + t61
     #) + t43 * (t24 * (-t54 + t34) * t27 + t55 * (t22 * (t15 * t39 + t5
     #0 * t40) - t51 * t15 * t33 * t40) * t46) + t45 * (t57 + t56 * (-t8
     # * t24 + (t50 - t15) * t14 * t10) * t46) + t54 * t49 * t51 * t11 *
     # t27 - t25 * t44 * t42 * t46)
      t36 = 0.1D1 / t2
      t54 = 0.1D1 / t6
      t65 = 0.1D1 / t1
      t30 = 0.1D1 / t30
      t70 = t13 * t22 + t51
      t71 = t5 * t43
      t72 = t11 * t15
      t73 = t11 * t65
      t54 = t54 * t36
      t74 = t54 * t28
      t75 = t16 * t36
      t76 = t75 * t65 + t33
      t77 = t5 ** 2
      t78 = t17 * t36
      t55 = t55 * t46
      t18 = t11 * (t5 * (t60 * t36 * (t73 - t49) * t23 * t42 + t27 * t36
     # * (-t61 * t33 + t73 * t38) * t23 * t24) + t77 * (t43 * t14 * t46 
     #* t76 * t23 * t42 + t55 * t33 * t11 * (t37 * (-t78 * t65 - 1) + t3
     #4) * t23) - t74 * t10 * t4 * t41 * t33 * t30) + t71 * t24 * t42 * 
     #t31 * t33 * t36 * t48 * t23 + t58 * (t72 * t36 * t27 + t60 * t15 *
     # t36 * (t51 * t32 * t33 - 1) + t71 * ((t18 * t22 + t17 * (t51 * t3
     #3 + t22)) * t65 * t36 + t22) * t46) * t23 + t74 * t41 * (-t10 * t4
     #3 * t65 * t30 + t49 * t32 + t73 * t32) + (t71 * (t22 * (t34 - t37)
     # + (t16 * (t34 * t70 - t37 * t70) + t17 * (t5 * t64 * t15 * t33 - 
     #t4 * t22 ** 2 + t34 * t22)) * t65 * t36) * t46 - t37 * t11 * t15 *
     # t36 * t27 + t36 * t4 * t31 * (t50 + t15)) * t23 * t24
      t37 = t17 * t33
      t38 = t37 + t16
      t48 = mt ** 2
      t58 = t19 ** 2
      t61 = t23 ** 2
      t64 = t6 ** 2
      t70 = t21 * t9
      t73 = t70 * t5
      t74 = t73 * t42 * t33
      t79 = t42 * t22
      t80 = t1 * t33
      t81 = t9 ** 2
      t82 = t56 * t46 * t43
      t83 = t33 * t23 * t14 * t77
      t2 = t83 * (t67 * t48 * t4 * t58 * t31 + t21 * (t81 * (-t31 * t41 
     #+ t57 * t43 * (t11 * (t77 * t64 * t22 * t46 - t8 * t27) - t52 + t1
     #0 ** 2 * t45 * t22 * t46)) + t82 * (-t69 * t10 * t27 + t35 * t39 *
     # t27 + t4) * t9)) + t51 * t23 * (t9 * (t44 * (-t25 * t5 * t42 * t3
     #8 * t47 - t74 * t10 * t46) + t45 * (t79 * t5 * (t70 * t3 ** 2 * t4
     #1 * t33 - t8 * t38 + t7 * (-t37 - t16)) * t47 + t74 * (t53 * t10 -
     # t7) * t46) - t73 * t64 * t42 * t31 * t33 * t61 + t71 * t56 * (t38
     # * t4 + t21 * (t63 * t24 * t33 * t64 + t67 * t7 * t33 * t9)) * t46
     # + t80 * (t66 * t15 + (t64 * t42 * t61 + t41) * t9 * t22 * t5 * t2
     #) * t27 * t11) + t35 * t50 * (t27 * ((-t5 * t11 * t27 - t23) * t22
     # * t6 + t11) - t60) * t58 * t48)
      t3 = (t21 + t17) * t33
      t10 = t46 * t9
      t35 = t80 * t72 * t6
      t8 = t83 * (t21 * (t66 * (t43 * (-t22 * t40 * t46 + t27) - t25 * t
     #45 * t46 + t52 * t31 - t53 * t11 * t27) * t81 + t79 * t47 * t45 * 
     #(t8 * t20 + t7 * (t8 + t20)) * t81) + t48 * t41 * t58 * (t29 * t27
     # - t31))
      ret = -4 * t18 - 64 * t2 - 8 * t26 - 32 * t14 * (t77 * (t42 * (t11
     # * (t43 * (t46 * (t9 * (t3 + t16) + t48 * t22 * t36 * t33 * t65 * 
     #t19) + t22 * (t17 * t76 + t16) * t47) - t35 * t10 - t48 * t36 * t3
     #1 ** 2 * t32 * t33 * t65 * t19) * t23 + t10 * t22 * t43 * t6 * (-t
     #3 - t16) * t61) + t82 * t49 * t48 * t21 * (-2 * t7 * t19 + 1) * t3
     #6 * t23 * t65 * t19) + (-t35 * t9 * t27 * t61 - t80 * t62 * t43 * 
     #t15 * t46 * (t7 + t20)) * t42 * t5 + t54 * t48 * t41 * t28 ** 2 * 
     #t19 * t33 * t65 * t32 * (-t20 * t30 + t32)) - 16 * t23 * t24 * t5 
     #* (t68 * t9 * t31 * (-t34 * t32 + t4) + t55 * t15 * (t33 * (t78 + 
     #t1) + t75) * t43 + t56 * t51 * (t21 * t33 + (t12 ** 2 * t15 ** 2 *
     # t33 + t13 * t16 ** 2) * t65 * t36) * t47 * t45) - 128 * t8 + 24 *
     # t59 * t62 * t50 * t1 * t32

      hjetmass_triangle_pppm_s34_s134_0_dp = ret/32d0/(0,1d0)
      return

      end function
