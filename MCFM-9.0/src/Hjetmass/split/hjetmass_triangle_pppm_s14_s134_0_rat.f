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
 

      complex*32 function hjetmass_triangle_pppm_s14_s134_0_rat
     &     (i1,i2,i3,i4,za,zb)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          
          complex*32 ret
          real*16 cg

      t1 = za(i2, i4)
      t2 = zb(i3, i2)
      t3 = za(i2, i3)
      t4 = zb(i4, i2)
      t5 = t3 * t2
      t6 = t1 * t4
      t7 = t6 + t5
      if ( real(t7) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t6 = cg * sqrt(t7 ** 2) + t5 + t6
      t7 = za(i1, i2)
      t8 = za(i1, i4)
      t9 = zb(i4, i3)
      t10 = zb(i2, i1)
      t11 = za(i1, i3)
      t12 = t11 * t2 + t4 * t8
      t13 = zb(i3, i1)
      t14 = zb(i4, i1)
      t15 = za(i3, i4)
      t16 = t11 * t13 + t14 * t8
      t17 = 0.1q1 / t6
      t18 = 2 * t7 * t15
      t19 = -t18 * t10 * t9 * t17 + t16
      t20 = (0.1q1 / 0.2q1)
      t21 = t20 * t6
      t22 = t1 * t10
      t23 = t22 * t9
      t24 = t15 * (-t21 * t13 + t23)
      t6 = -t7 * t15 * t10 * t9 + t20 * t6 * t16
      t16 = t3 * t10 * t9
      t20 = t15 * (t21 * t14 + t16)
      t21 = t9 * (t18 * t2 * t17 + t8)
      t13 = t15 * (2 * t23 * t17 - t13)
      t14 = t15 * (2 * t16 * t17 + t14)
      t16 = t10 * t19
      t23 = (-t13 * t4 - t14 * t2 + t16) * t7
      t11 = (t1 * t9 * (t18 * t4 * t17 - t11) - t19 * t7 + t21 * t3) * t
     #10
      t18 = 0.1q1 / t15
      t11 = 0.1q1 / t11
      t6 = 0.1q1 / t6
      t20 = 0.1q1 / t20
      t25 = 0.1q1 / t7
      t14 = 0.1q1 / t14
      t23 = 0.1q1 / t23
      t8 = 0.1q1 / t8
      t26 = 0.1q1 / t9
      t27 = 0.1q1 / t19
      t3 = 0.1q1 / t3
      t28 = t16 * t11 + t25
      t29 = t1 ** 2
      t30 = t11 ** 2
      t19 = t10 ** 2 * t19
      ret = -32 * t26 * t17 * t8 * ((t2 * t4 * t28 * t3 ** 2 + t18 * t12
     # * t10 * (t21 * t25 ** 2 * t27 - t19 * t21 * t30 + t16 * t2 * t15 
     #* t9 * (-2 * t5 * t17 + 1) * t30) * t3) * t1 * t29 + t2 ** 2 * (t7
     # * t10 * t13 ** 2 * t14 * t23 - t24 ** 2 * t6 * t20)) - 64 * t3 * 
     #t17 * t8 * t2 * t29 * (-t19 * t29 * t4 * t12 * t17 * t30 + t2 * t2
     #6 * t28 - t22 * (t10 * t11 + t25 * t27) * t17 * t12)

      hjetmass_triangle_pppm_s14_s134_0_rat = ret/32q0/(0,1q0)
      return

      end function
