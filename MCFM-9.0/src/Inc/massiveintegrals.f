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
 
      complex(dp):: 
     & I41x2x3x4(2),I4m1x2x3x4(2),I41x2x4x3,
     & I3m12x3x4(2),I312x3x4(2),
     & I3m13x2x4(2),I313x2x4(2),
     & F41x2x3x4(2),F4m1x2x3x4(2),
     & I461x2x3x4(2),I46m1x2x3x4(2),F212(2),
     & I3m1x23x4,I3m2x3x41,I32x3x41,I31x23x4,I2m,I2h23,F2m23,Im2m,I23
      common/massiveintegrals/
     & I41x2x3x4,I4m1x2x3x4,I41x2x4x3,
     & I3m12x3x4,I312x3x4,
     & I3m13x2x4,I313x2x4,
     & F41x2x3x4,F4m1x2x3x4,
     & I461x2x3x4,I46m1x2x3x4,F212,
     & I3m1x23x4,I3m2x3x41,I32x3x41,I31x23x4,
     & I2m,I2h23,F2m23,Im2m,I23
!$omp threadprivate(/massiveintegrals/)
