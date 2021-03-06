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
 

      double complex function hjetmass_triangle_pmpm_s124_0_mhsq_rat_dp 
     &     (i1,i2,i3,i4,za,zb)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double complex ret
          double precision cg

      t1 = za(i3, i4)
      t2 = zb(i3, i1)
      t3 = za(i1, i4)
      t4 = za(i2, i4)
      t5 = zb(i3, i2)
      t6 = t2 * t3 + t4 * t5
      t7 = za(i1, i2)
      t8 = zb(i2, i1)
      t9 = zb(i4, i1)
      t10 = za(i1, i3)
      t11 = zb(i4, i2)
      t12 = t7 * t8
      t13 = t3 * t9
      t14 = t11 * t4 + t12 + t13
      t15 = za(i2, i3)
      t16 = zb(i4, i3)
      t17 = t10 * t2
      t18 = t15 * t5
      t19 = t1 * t16
      t20 = t17 + t18 + t19
      if ( dreal(t20) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t18 = cg * cdsqrt(t20 ** 2) + t17 + t18 + t19
      t18 = 0.1D1 / t18
      t12 = -2 * t17 * t14 * t18 + t12 + t13
      t13 = 2 * t10
      t17 = -t13 * t5 * t14 * t18 + t11 * t3
      t19 = t10 * t12
      t20 = t2 * (t15 * t17 + t19)
      t21 = -t16 * t3 + t5 * t7
      t22 = t16 * t4 + t2 * t7
      t23 = 2 * t1 * t2 * t14 * t18 - t4 * t8
      t24 = t2 * t12
      t25 = t10 * (t16 * t23 - t24)
      t26 = -2 * t15 * t2 * t14 * t18 + t4 * t9
      t11 = t13 * t14 * t16 * t18 - t11 * t7
      t13 = t2 * (t1 * t11 - t19)
      t19 = t10 * (t26 * t5 + t24)
      t7 = 0.1D1 / t7
      t24 = 0.1D1 / t15
      t8 = 0.1D1 / t8
      t27 = 0.1D1 / t1
      t20 = 0.1D1 / t20
      t3 = 0.1D1 / t3
      t13 = 0.1D1 / t13
      t11 = 0.1D1 / t11
      t17 = 0.1D1 / t17
      t9 = 0.1D1 / t9
      t28 = t1 ** 2
      t29 = t2 ** 2
      t30 = t6 * t15
      t31 = t30 * t27
      t32 = t4 * t22 * t7
      t33 = t13 * t7 * t15
      t34 = t24 * t17
      t35 = t9 * t18 * t8
      t14 = 0.1D1 / t14
      t36 = t13 ** 2
      t37 = t20 ** 2
      t38 = t4 ** 2
      t39 = t12 ** 2
      t40 = t6 * t16
      t41 = t2 * t21
      t42 = t5 * t22
      t43 = t1 * t22
      t5 = 0.1D1 / t5
      t19 = 0.1D1 / t19
      t25 = 0.1D1 / t25
      t16 = 0.1D1 / t16
      t16 = t2 * t26 * t16
      t44 = t4 * t26 * t7
      t45 = t2 * t23 * t5
      t23 = t7 * (t4 * t23 * t3 - t16) * t25
      t46 = t24 * t19 * t3
      t47 = t4 * t7 * t3
      t5 = t35 * t29 * (t10 * (t16 * t22 * t7 * t25 + t23 * t31 + (t45 *
     # t19 + t44 * (-t19 - t25)) * t3 * t6 + t46 * t26 * (-t2 * t6 * t5 
     #+ t32) * t1) + t47 * (t43 * t13 - t30 * t20) * t12)
      ret = -48 * t35 * t12 * t29 * (-t34 * t28 * t12 * t22 * t3 * t20 +
     # t3 * (t6 * t12 * t17 - t32) * t20 * t1 + t33 * (t4 * t6 * t3 + (-
     #t31 + t22) * t11 * t12)) + 32 * t2 * (t4 * t38 * t7 * t3 * t14 * (
     #t24 * t8 + t27 * t9) + t35 * t3 * t7 * t15 * t4 * t2 * t1 * (-t40 
     #* t36 - t42 * t37 + t41 * (t37 - t36)) * t39 + t35 * ((t40 * t27 *
     # t11 ** 2 * t13 + t2 * (t41 + t40) * t11 * t36) * t7 * t15 ** 2 + 
     #t17 * t3 * t20 * t28 * (t29 * t21 * t20 + t42 * (-t2 * t20 - t34))
     #) * t12 * t39) - 64 * t3 * t18 * t14 * t7 * t38 * t2 * (t8 * (t43 
     #* t24 - t6) + t9 * (t31 - t22)) - 8 * t9 * t8 * t4 * t29 * (t10 * 
     #(t23 * t27 + t46 * (-t45 + t44)) + t39 * (t34 * t1 * t20 * t3 + t3
     #3 * t27 * t11) + t47 * (t20 - t13) * t12) + 16 * t5

      hjetmass_triangle_pmpm_s124_0_mhsq_rat_dp = ret/32d0/(0,1d0)
      return

      end function
