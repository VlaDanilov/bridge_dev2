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
 

      complex*32 function hjetmass_triangle_pppm_s34_s234_0_rat
     &     (i1,i2,i3,i4,za,zb)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          
          complex*32 ret
          real*16 cg

      t1 = zb(i3, i1)
      t2 = za(i1, i2)
      t3 = zb(i2, i1)
      t4 = za(i1, i3)
      t5 = t2 * t3
      t6 = t4 * t1
      t7 = t6 + t5
      if ( real(t7) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t7 = cg * sqrt(t7 ** 2) + t5 + t6
      t8 = za(i3, i4)
      t9 = zb(i3, i2)
      t10 = za(i1, i4)
      t11 = za(i2, i3)
      t12 = (0.1q1 / 0.2q1)
      t13 = t12 * t7
      t14 = t11 * t3
      t15 = -t9 * (t14 * t10 + t13 * t8)
      t16 = za(i2, i4)
      t17 = zb(i4, i2)
      t18 = zb(i4, i1)
      t19 = zb(i4, i3)
      t20 = t2 * t9
      t21 = -t11 * (t13 * t19 + t20 * t18)
      t22 = t11 * t9
      t23 = t22 * (t13 - t5)
      t24 = t20 * t14 * (t10 * t18 + t5 + t6)
      t19 = t12 * t14 * t7 * (-t10 * t19 + t20) - t24
      t12 = -t12 * t20 * t7 * (t18 * t8 - t14) - t24
      t20 = t10 * t11
      t24 = t9 * (-t20 * t1 + t13 * t16)
      t8 = 0.1q1 / t8
      t25 = 0.1q1 / t2
      t17 = 0.1q1 / t17
      t26 = 0.1q1 / t23
      t16 = 0.1q1 / t16
      t27 = 0.1q1 / t18
      t12 = 0.1q1 / t12
      t21 = 0.1q1 / t21
      t7 = 0.1q1 / t7
      t19 = 0.1q1 / t19
      t28 = t18 * t24
      t6 = t1 * t22 * (t13 - t6) + t28
      t13 = t1 * t16
      t29 = t3 * t8
      t30 = t29 + t13
      t31 = t10 ** 2
      t32 = t15 ** 2
      t33 = t3 ** 2
      t34 = t19 ** 2
      t35 = t15 * t8
      t36 = t22 * t1
      t14 = t14 * t1
      t37 = t1 * t15
      t38 = t19 * t16
      t1 = t7 * t17 * t1
      t7 = t24 * t16
      t17 = t3 * t24
      t24 = t10 * t3 * t19 - t21
      t39 = -32
      ret = t39 * (t1 * (t21 * (t17 * t16 + (-t7 - t35) * t21 * t23 * t1
     #8) + t3 * t15 * t30 * t27 * t26 - t38 * t3 * (t17 + t37) * t10 + t
     #31 * t3 * t33 * t23 ** 2 * t8 * t34) + t1 * (t14 * t32 * t8 * t16 
     #* (-t2 * t15 * t12 + t27) * t26 ** 2 + t3 * (t13 * t11 * t8 * t25 
     #* (-t15 * t19 + t22 * t21) - t36 * t34 * (t16 * t2 + t4 * t8) * t2
     #3 * t33 + t34 * (t16 * t6 + t35 * t18) * t23 * t3) * t31 + t20 * t
     #8 * t25 * t16 * t21 * (-t28 * t23 * t21 + t37) + t38 * t8 * t11 * 
     #t33 * (-t36 * t25 + (-t14 * t9 + t25 * t6) * t19 * t23) * t10 * t3
     #1 + t37 * t16 * t21 - t5 * t32 * t12 * t30 * t26)) + 64 * t29 * t1
     # * (t20 * t7 * t24 * t25 + t15 * t24)

      hjetmass_triangle_pppm_s34_s234_0_rat = ret/32q0/(0,1q0)
      return

      end function
