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
 

      complex*32 function hjetmass_triangle_pppp_s234_0_mhsq_rat
     &     (i1,i2,i3,i4,za,zb)
      implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          
          complex*32 ret
          real*16 cg

      t1 = zb(i2, i1)
      t2 = za(i2, i4)
      t3 = zb(i4, i3)
      t4 = za(i1, i2)
      t5 = zb(i3, i1)
      t6 = za(i2, i3)
      t7 = zb(i3, i2)
      t8 = zb(i4, i2)
      t9 = za(i3, i4)
      t10 = t6 * t7
      t11 = t2 * t8
      t12 = t9 * t3
      t13 = t11 + t12 + t10
      t14 = za(i1, i3)
      t15 = za(i1, i4)
      t16 = zb(i4, i1)
      t17 = t4 * t1
      t18 = t14 * t5
      t19 = t15 * t16
      t20 = t18 + t19 + t17
      if ( real(t20) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t20 = cg * sqrt(t20 ** 2) + t17 + t18 + t19
      t20 = 0.1q1 / t20
      t21 = -2 * t4
      t22 = t21 * t5 * t13 * t20 + t2 * t3
      t23 = t21 * t16 * t13 * t20 - t3 * t6
      t24 = t14 * t22
      t25 = t15 * t23
      t26 = t1 * (t25 + t24)
      t27 = t16 * t2 + t5 * t6
      t21 = t21 * t1 * t13 * t20 + t10 + t11
      t28 = -2 * t14
      t29 = t28 * t1 * t13 * t20 + t8 * t9
      t30 = -2 * t15
      t31 = t30 * t1 * t13 * t20 - t7 * t9
      t32 = t1 * (t21 * t4 + t24)
      t33 = t28 * t5 * t13 * t20 + t10 + t12
      t11 = t30 * t16 * t13 * t20 + t11 + t12
      t34 = t5 * t29
      t35 = t21 * t1
      t36 = t4 * (t34 + t35)
      t37 = t16 * t31
      t38 = t4 * (t34 + t37)
      t30 = t30 * t5 * t13 * t20 + t2 * t7
      t28 = t28 * t16 * t13 * t20 + t6 * t8
      t39 = 0.1q1 / t15
      t40 = 0.1q1 / t2
      t41 = 0.1q1 / t23
      t36 = 0.1q1 / t36
      t6 = 0.1q1 / t6
      t42 = 0.1q1 / t21
      t9 = 0.1q1 / t9
      t26 = 0.1q1 / t26
      t13 = 0.1q1 / t13
      t38 = 0.1q1 / t38
      t43 = 0.1q1 / t4
      t32 = 0.1q1 / t32
      t44 = -t26 - t38
      t45 = t15 ** 2
      t46 = t6 ** 2
      t47 = t31 ** 2
      t48 = t1 ** 2
      t49 = t43 ** 2
      t50 = t9 ** 2
      t51 = t31 * t42
      t52 = t1 * t11
      t53 = t15 * t42
      t54 = t1 * t45
      t55 = t51 * t5 * t40
      t56 = t6 * t1
      t57 = t15 * t28
      t58 = t4 * t16
      t59 = t58 * t39
      t60 = t4 * t29
      t61 = t60 * t44
      t62 = t14 * t21
      t63 = t7 * t39
      t64 = t4 * t31 * t39
      t65 = t3 * t43
      t66 = t5 * t8
      t67 = t7 * t16
      t68 = t9 * t43
      t69 = 0.1q1 / t31
      t70 = t39 ** 2
      t71 = t4 ** 2
      t72 = t14 * t30
      t73 = t11 * t32
      t74 = t5 * t40
      t75 = t1 * t9
      t76 = t6 * t8
      t77 = t76 * t9
      t34 = (t74 * t54 * t28 * t50 * t43 * t32 - t24 * t53 * t1 * t40 * 
     #t50 * t49 + t40 * (t39 * (t1 * (t4 * (-t18 * t21 * t26 + 1) - t34 
     #* t26 * t71) - t18) + t72 * t58 * t70 * t41 + t35 * (1 - t5 * t38 
     #* (t62 + t60)) * t69) * t46 + t9 * (t75 * t58 * t51 * t29 * t36 + 
     #t74 * (t14 * (t1 * (t21 * t44 + t73) + t59 * t47 * t36 * t42) + t6
     #1 * t1 + t57 * t1 * t32)) * t6 + t16 * t47 * t40 * t42 ** 2 * t50 
     #* (t60 * t5 * t36 - 1)) * t20 * t27
      t58 = t32 ** 2
      t78 = t22 * t40
      t79 = t21 * t9
      t80 = t28 * t6
      t81 = t78 * t23
      t82 = t40 * t6
      t83 = t51 * t39
      t84 = t7 * t40
      t85 = t13 * t9
      t86 = t17 * t39
      t87 = t19 * t43
      t88 = t16 * t6
      t89 = t20 * t27
      t90 = t30 * t16 * t40
      t8 = t89 * (t13 * (t88 * (t6 * (-t12 * t40 - t8) - t84) + t75 * (t
     #76 + t84)) + t51 * t13 * (t82 * (t18 + t17) * t39 * t3 + t1 * t50 
     #* (t10 * t40 + t8)) + t56 * (t35 * t4 * (t16 * (t11 * t6 * t9 + (t
     #4 * t39 * t6 + t9) * t40 * t30) + t9 * (t33 * t40 + t80) * t5) * t
     #26 ** 2 + t6 * (-t75 * t23 - t78 * (t21 * t16 * t41 + t1) * t39 * 
     #t4 + t90 * t71 * t21 * t70 * t41) * t26 - t25 * t48 * t21 * t50 * 
     #t58 + t90 * t9 * t32) * t14)
      t8 = t8 + t89 * (t14 * (t48 * (t15 * (t73 * t78 * t42 * t50 * t43 
     #- t80 * t5 * t9 * (t79 + t78) * t58) + t45 * (-t74 * t22 * t28 * t
     #50 * t43 * t58 - t81 * t42 * t50 * t49 * t32) + t82 * t26 * (t71 *
     # t5 * t21 * t33 * t39 * t6 * t26 - t22 * t9)) + t16 * t50 * t32 * 
     #(t15 * t30 * t40 * t43 + t31 * t6) * t1 - t81 * t15 * t9 * t58 * (
     #t68 * t15 + t6) * t1 * t48 + t85 * t5 * (-t84 * t43 + t76 * (t83 -
     # t43))) + t51 * t40 * t13 * (-t88 + t75) * t3 + t85 * (-t84 * t19 
     #* t43 + t76 * (-t87 + t51 * (t86 - t16))))
      t25 = t51 * t9
      t4 = t1 * (t40 * (t26 * (t86 * t22 * t6 - t7 * t27 * t9) + t5 * (t
     #21 ** 2 * t69 * t6 * t38 + t79 * t38) + t6 * (t16 * t22 + t27 * t3
     #) * t32) + t79 * t88 * (t38 + t32) + t77 * t27 * (-t26 + t32)) + t
     #13 * (t74 * (t65 + t63) + t88 * t65) + t64 * t16 * (t74 * (t6 + t2
     #5) + t75 * t6) * t36 + t75 * (t87 * t78 * t32 + t63 * t13) + t40 *
     # (-t5 * t69 * t20 * (t62 * t43 + t29) * t46 + t20 * t42 * (t5 * (t
     #14 * t39 * t42 * t47 + t51 * t29) + t1 * t43 * (t15 * t33 + t72)) 
     #* t50 - t17 * t63 * t6 * t26 - t39 * t46 * t20 * (t33 * t4 + t24) 
     #* t41 * t16 + t75 * t65 * t15 * t32) * t27 + t26 * (t23 * t6 + t78
     #) * t9 * t48
      t22 = t89 * t40 * t13 * (t67 * t25 - t56 * t3)
      t24 = 16
      ret = t24 * (t13 * (t77 * (t1 * t39 + t16 * t43) * t2 + t74 * (t12
     # * t39 * t6 + t10 * t68)) + t13 * ((t1 * t3 + t66) * t6 * t39 + t6
     #8 * (t66 + t67)) + (t50 * (t40 * (t43 * (t18 * (t52 * t15 * t32 + 
     #t51) - t53 * (t52 + t37)) + t54 * t23 * t42 * t49) + t59 * t14 * t
     #47 * t42 * (t56 + t55) * t36 + t48 * t6 * t32 * (t11 * t14 + t57))
     # + t9 * (t16 * (t60 * t55 * t36 * t6 + t1 * (t62 * t44 + t61) * t4
     #6) + t63 * t51 * t40 * t13 * (t18 + t17)) - t65 * t40 * t6 * t13 *
     # (t18 + t19) + t59 * t40 * t41 * t46 * (t21 - t64)) * t20 * t27 + 
     #t34) + 32 * t8 + 4 * t40 * (t27 * (-t39 * t3 * t41 * t6 + t68 * t7
     # * t42) + t5 * (t6 * (t21 * t69 * t43 + t39) + t9 * (t83 + t43))) 
     #+ 8 * t4 - 48 * t22

      hjetmass_triangle_pppp_s234_0_mhsq_rat = ret/32q0/(0,1q0)
      return

      end function
