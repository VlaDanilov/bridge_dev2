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
 

      complex*32 function hjetmass_triangle_pmpm_s34_s134_0
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          
          complex*32 ret
          double precision mt
          real*16 cg

      t1 = za(i3, i4)
      t2 = zb(i4, i3)
      t3 = za(i1, i3)
      t4 = zb(i3, i1)
      t5 = za(i1, i4)
      t6 = zb(i4, i1)
      t7 = t3 * t4
      t8 = t5 * t6
      t9 = t8 + t7
      if ( real(t9) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t9 = cg * sqrt(t9 ** 2) + t7 + t8
      t10 = za(i1, i2)
      t11 = zb(i2, i1)
      t12 = za(i2, i3)
      t13 = za(i2, i4)
      t14 = t12 * t4
      t15 = t13 * t6
      t16 = t15 + t14
      t17 = 0.1q1 / t9
      t18 = t1 * t10
      t19 = -2 * t18
      t20 = t2 * (t19 * t4 * t17 + t13)
      t21 = t2 * (t19 * t6 * t17 - t12)
      t22 = zb(i3, i2)
      t23 = zb(i4, i2)
      t24 = t12 * t22
      t25 = t13 * t23
      t19 = t19 * t11 * t2 * t17 + t24 + t25
      t26 = t3 * t20
      t27 = (t10 * t19 + t21 * t5 + t26) * t11
      t28 = t1 * t2
      t21 = 0.1q1 / t21
      t9 = 0.1q1 / t9
      t19 = 0.1q1 / t19
      t29 = 0.1q1 / t11
      t27 = 0.1q1 / t27
      t30 = -t28 - t24 - t25
      t31 = t20 ** 2
      t32 = t6 ** 2
      t33 = mt ** 2
      t34 = t27 ** 2
      t35 = t27 * t34
      t36 = t17 ** 2
      t37 = t10 ** 2
      t38 = t9 ** 2
      t39 = t19 ** 2
      t40 = t11 ** 2
      t41 = t11 * t40
      t42 = t4 ** 2
      t43 = t5 * t32
      t44 = t4 * t9
      t45 = t33 * t4
      t46 = t45 * t6 * t36
      t47 = t32 * t31 * t21 ** 2
      t48 = t28 * t38
      t49 = t20 * t11
      t50 = t9 * t31
      t51 = t42 * t29
      t52 = t19 * t4
      t53 = t21 * t16 ** 2 * t5
      t54 = t25 + t24
      t36 = t33 * t42 * t36
      t55 = 0.1q1 / t1
      t56 = 0.1q1 / t10
      t57 = 0.1q1 / t2
      t58 = t5 * t12
      t59 = -t18 + t58
      t60 = t15 + t14
      t61 = t8 + t7
      t62 = t13 ** 2
      t63 = t2 * t16
      t64 = t63 * t5
      t65 = t64 * t13
      t66 = t18 * t4
      t67 = t10 * t35
      t68 = t19 * t56 * t29
      t69 = t16 * t21
      t12 = t69 * (t27 * (t40 * (t64 * t18 * t31 * t34 + t9 * t20 * t10 
     #* (t20 * t59 * t4 + t5 * (t20 * t6 + t63) * t13) * t27) + t9 * (t6
     # * (t59 * t4 * t21 * t31 + t65 * t21 * t20) - t65 * t4 + t43 * t13
     # * t31 * t21) + t49 * t9 * (t20 * (t5 * (t7 * t60 + t8 * t60) - t6
     #6 * t61) + t65 * t61) * t27) + t20 * (t5 * (-t49 * t16 * t19 * t39
     # + t67 * (t12 ** 2 * t22 ** 2 + t62 * t23 ** 2) * t20 * t16 * t40)
     # + t45 * t13 * (-t68 + t27)) * t57 * t55)
      t18 = t4 * t11
      t22 = t53 * t48 * t20 * (t10 * (t8 * t7 * t40 * t20 * t35 - t4 * t
     #6 * t21 * t27 - t18 * t61 * t34) + t37 * (t41 * t20 * t61 * t35 - 
     #t4 * t40 * t34) + t52 * (t6 * t29 * t21 + t19))
      t23 = t9 * t21 * t6
      t45 = t16 * t5
      t17 = t69 * t20 * (t45 * (t25 * t24 * t67 * t40 * t20 + (t27 * (t1
     #0 * (t20 * t27 * t40 + t18 * t28 * (-2 * t7 * t17 + 1) * t27) - t4
     #) + t68 * t4 - t20 * t39 * t56) * t17 * t33) * t57 * t55 + t44 * t
     #27 * (t5 * (-t15 - t14) + t66) + t5 * (t10 * (t23 * t30 * t34 * t1
     #1 + t54 * t35 * t40) - t23 * t39) * t20 * t16)
      t1 = t69 * t49 * t34 * (t20 * (t5 * ((t14 * t54 + t15 * t54) * t57
     # * t55 + t14 + t15) + t4 * (-t54 * t57 - t1) * t10) + t45 * t13 * 
     #(t54 * t55 + t2))
      t2 = t20 * (t27 * (t4 * (-t13 * t4 + t21 * (t4 * t10 * t16 + t15 *
     # t20)) * t57 + t69 * t13 * (t11 * t13 + (t13 * t3 - t58) * t56 * t
     #4) * t55) - t69 * (t62 * t55 * t56 + t51 * t57) * t19)
      t4 = 64
      ret = t4 * (t53 * (t10 * (t50 * t40 * (t7 * t30 - t8 * t54) * t35 
     #+ t44 * t49 * (t26 * t28 * t6 * t21 * t9 + t24 + t25) * t34 + t36 
     #* t27) - t36 * t29 * t19) + t53 * (t10 * (t27 * (-t46 * t20 * t21 
     #+ t48 * (t47 + t42)) + t49 * (-t46 * t5 + t28 * (t43 * t20 * t21 *
     # t38 + t44)) * t34 + t50 * t28 * t40 * (t3 ** 2 * t42 * t9 + t5 **
     # 2 * t32 * t9 - t8) * t35) + t19 * (t20 * (t9 * (-t49 * t39 + t52)
     # + t46 * t29 * t21) + t48 * (-t6 * t31 * t21 * t19 - t11 * t31 * t
     #39 - t47 * t29 - t51)) + t37 * (t48 * t40 * t6 * t31 * t21 * t34 +
     # t50 * t41 * t30 * t35) + t48 * t10 * t37 * t40 ** 2 * t31 * t35))
     # + 16 * t12 + 128 * t22 + 32 * t17 - 8 * t1 - 4 * t2

      hjetmass_triangle_pmpm_s34_s134_0 = ret/32q0/(0,1q0)
      return

      end function
