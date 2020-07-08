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
 

      complex*32 function hjetmass_triangle_pppm_s34_s134_0_rat
     &     (i1,i2,i3,i4,za,zb)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          
          complex*32 ret
          real*16 cg

      t1 = za(i1, i4)
      t2 = za(i2, i4)
      t3 = zb(i2, i1)
      t4 = zb(i4, i1)
      t5 = za(i1, i2)
      t6 = zb(i4, i2)
      t7 = t5 * t3
      t8 = t2 * t6
      t9 = t8 + t7
      if ( real(t9) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t9 = cg * sqrt(t9 ** 2) + t7 + t8
      t10 = za(i3, i4)
      t11 = (0.1q1 / 0.2q1)
      t12 = t11 * t9
      t13 = t1 * t4
      t14 = t13 * (t12 - t7)
      t15 = za(i2, i3)
      t16 = t1 * t3
      t17 = t4 * (t12 * t10 + t16 * t15)
      t18 = zb(i3, i2)
      t19 = t5 * t4
      t20 = t15 * t18
      t7 = t7 * t13 * (t20 + t7 + t8)
      t21 = t11 * t19 * t9 * (-t10 * t18 + t16) - t7
      t22 = zb(i4, i3)
      t23 = t1 * (t12 * t22 + t19 * t18)
      t24 = zb(i3, i1)
      t12 = t1 * (-t2 * t18 * t4 + t12 * t24)
      t13 = t11 * t9 * (t10 * t22 + t24 * za(i1, i3)) - t20 * t13
      t7 = t11 * t16 * t9 * (-t15 * t22 + t19) - t7
      t7 = 0.1q1 / t7
      t11 = 0.1q1 / t14
      t16 = 0.1q1 / t23
      t19 = 0.1q1 / t4
      t10 = 0.1q1 / t10
      t17 = 0.1q1 / t17
      t20 = 0.1q1 / t5
      t21 = 0.1q1 / t21
      t9 = 0.1q1 / t9
      t15 = 0.1q1 / t15
      t22 = t3 ** 2
      t24 = t7 ** 2
      t7 = t3 * t7
      t25 = t16 * t15
      t9 = t9 * t10 * t2 ** 2
      ret = 32 * t9 * (t19 * (-t8 * t3 * (t7 * t23 + t15) * t20 ** 2 + t
     #14 * (t2 * (t6 * (-t18 * t13 * t16 ** 2 * t15 ** 2 + t22 * t18 * t
     #13 * t24 - t3 * t22 * t23 * t24) + t12 * t22 * t24 * t6 ** 2) + t3
     # * t18 * (t25 + t7)) * t20) + t1 ** 2 * t22 ** 2 * t4 * t11 * (t5 
     #* t18 * t21 + t17)) - 64 * t9 * t3 * t12 * t6 * t20 * t19 * (t25 +
     # t7)

      hjetmass_triangle_pppm_s34_s134_0_rat = ret/32q0/(0,1q0)
      return

      end function
