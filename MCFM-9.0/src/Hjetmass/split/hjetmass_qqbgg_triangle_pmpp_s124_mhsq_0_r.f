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
 

      double complex function hjetmass_qqbgg_triangle_pmpp_s124_mhsq_0_rat_dp
     &     (i1,i2,i3,i4,za,zb)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double complex ret
          double precision cg

      t1 = za(i2, i3)
      t2 = zb(i3, i1)
      t3 = za(i1, i2)
      t4 = za(i2, i4)
      t5 = zb(i4, i3)
      t6 = t2 * t3 + t4 * t5
      t7 = za(i1, i3)
      t8 = zb(i3, i2)
      t9 = za(i3, i4)
      t10 = t7 * t2
      t11 = t1 * t8
      t12 = t9 * t5
      t13 = t11 + t12 + t10
      if ( dreal(t13) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t11 = cg * cdsqrt(t13 ** 2) + t10 + t11 + t12
      t12 = zb(i2, i1)
      t13 = za(i1, i4)
      t14 = zb(i4, i1)
      t15 = zb(i4, i2)
      t16 = t3 * t12
      t17 = t13 * t14
      t18 = t15 * t4 + t16 + t17
      t11 = 0.1D1 / t11
      t19 = 0.2D1 * t1
      t20 = -t19 * t2 * t18 * t11 + t14 * t4
      t21 = 0.2D1 * t9 * t2 * t18 * t11 - t12 * t4
      t10 = -0.2D1 * t10 * t18 * t11 + t16 + t17
      t16 = t7 * (t10 * t2 + t20 * t8)
      t17 = 0.2D1 * t7
      t22 = t2 * (t1 * (-t17 * t8 * t18 * t11 + t13 * t15) + t10 * t7)
      t15 = t17 * t18 * t5 * t11 - t15 * t3
      t17 = -t13 * t5 + t3 * t8
      t19 = t19 * t18 * t5 * t11 + t14 * t3
      t18 = 0.1D1 / t18
      t23 = 0.1D1 / t3
      t12 = 0.1D1 / t12
      t24 = 0.1D1 / t4
      t25 = 0.1D1 / t13
      t14 = t14 * t12
      t9 = 0.1D1 / t9
      t16 = 0.1D1 / t16
      t26 = 0.1D1 / t15
      t21 = 0.1D1 / t21
      t22 = 0.1D1 / t22
      t27 = t2 * t22
      ret = 0.128D3 * t11 * t18 * t6 * (t5 * (t14 * t23 + t25) + (t25 * 
     #t3 + t14) * t24 * t2) + 0.32D2 * t24 * t12 * t11 * t23 * (t1 * (t2
     # * (-t27 * t17 * t19 + t6 * t9) - t6 * (t26 * t9 + t27) * t10 * t5
     #) - t2 * t6 * (t5 * t7 * t16 + t21) * t20 + t10 * (-t5 * (t13 * t2
     # + t4 * t8) * t9 ** 2 * t26 + t2 ** 2 * t17 * t22 ** 2 * (t15 * t2
     # + t19 * t8)) * t1 ** 2)

      hjetmass_qqbgg_triangle_pmpp_s124_mhsq_0_rat_dp = ret/32d0*(0,1d0)
      return

      end function
