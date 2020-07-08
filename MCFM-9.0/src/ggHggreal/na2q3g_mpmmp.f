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
 
      function na2q3g_mpmmp(j2,j3,j4,j5,j1,zb,za)
c--- Calculation of the amplitudes using results of S. Badger
c--- as adapted from routines written by F. Caola      
C  Arguments in call are strange because 
C this has been created from A0Hqbqggg_mp_mpp
      
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      complex(dp)::na2q3g_mpmmp
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer::j1,j2,j3,j4,j5
      real(dp)::qsq
      real(dp)::sH1, sH4, sH5, s125, s123, s234, s345
      complex(dp)::zab_1_pH_4, zab_5_pH_1, zab_5_pH_4, zab_3_pH_4, 
     & zab_1_pH5_4, zab_3_pH4_5,zab_1_23_4, zab_3_45_2,  zab_5_34_2
      complex(dp)::zbb_2_pH34_pH_4, zbb_4_pH_15_2
      complex(dp)::zab2, zab3
      real(dp):: s3, s4

c--- begin statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      zab3(j1,j2,j3,j4,j5)=
     & za(j1,j2)*zb(j2,j5)+za(j1,j3)*zb(j3,j5)+za(j1,j4)*zb(j4,j5)

      s3(j1,j2,j3)=s(j1,j2)+s(j1,j3)+s(j2,j3)
      s4(j1,j2,j3,j4)=s(j1,j2)+s(j1,j3)+s(j1,j4)
     &               +s(j2,j3)+s(j2,j4)+s(j3,j4)
c--- end statement functions

      qsq=
     & +s(j1,j2)+s(j1,j3)+s(j1,j4)+s(j1,j5)
     &          +s(j2,j3)+s(j2,j4)+s(j2,j5)
     &                   +s(j3,j4)+s(j3,j5)
     &                            +s(j4,j5)


      sH1 = s4(j2,j3,j4,j5)
      sH4 = s4(j1,j2,j3,j5)
      sH5 = s4(j1,j2,j3,j4)
      s125 = s3(j1,j2,j5)
      s123 = s3(j1,j2,j3)
      s234 = s3(j2,j3,j4)
      s345 = s3(j3,j4,j5)

      zab_1_pH_4 = -zab3(j1,j2,j3,j5,j4)
      zab_5_pH_1 = -zab3(j5,j2,j3,j4,j1)
      zab_5_pH_4 = -zab3(j5,j1,j2,j3,j4)
      zab_3_pH_4 = -zab3(j3,j1,j2,j5,j4)
      zab_1_pH5_4 = -zab2(j1,j2,j3,j4)
      zab_3_pH4_5 = -zab2(j3,j1,j2,j5)
      zab_1_23_4 = zab2(j1,j2,j3,j4)
      zab_3_45_2 = zab2(j3,j4,j5,j2)
      zab_5_34_2 = zab2(j5,j3,j4,j2)
      
      zbb_4_pH_15_2 = -zb(j4,j2) * s125 - zb(j4,j3)*zab2(j3,j1,j5,j2)
      zbb_2_pH34_pH_4 = zbb_4_pH_15_2

      na2q3g_mpmmp=
     & -za(j1,j3)**3/(za(j1,j2)*za(j1,j5)*za(j3,j4)*za(j4,j5))
     & -zab_3_45_2**3/(s345*za(j3,j4)*za(j4,j5)*zab_5_34_2*zb(j2,j1))
     & -sH1**2*zb(j4,j2)**3
     & /(s234*zab_5_pH_1*zab_5_34_2*zb(j3,j2)*zb(j4,j3))
     & -qsq**2*zb(j4,j1)*zb(j4,j2)**3
     & /(sH5*zab_5_pH_1*zab_5_pH_4*zb(j2,j1)*zb(j3,j2)*zb(j4,j3))
     & +za(j1,j3)**3*zb(j5,j4)**3
     & /(s123*za(j1,j2)*zab_1_pH5_4*zab_3_pH4_5)
     & -zab_3_pH_4**3*zb(j5,j2)**3
     & /(sH4*s125*zab_3_pH4_5*zb(j2,j1)*zbb_2_pH34_pH_4)
     & -zab_1_pH_4**3*zb(j4,j2)**3
     & /(za(j1,j5)*zab_1_23_4*zab_5_pH_4*zb(j3,j2)
     & *zb(j4,j3)*zbb_4_pH_15_2)

      return
      end
