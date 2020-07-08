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
 
      DATA (NHEL(IHEL,  1),IHEL=1,5) / +1, +1, +1, +1, +1/
      DATA (NHEL(IHEL,  2),IHEL=1,5) / -1, +1, +1, +1, +1/
      DATA (NHEL(IHEL,  3),IHEL=1,5) / +1, -1, +1, +1, +1/
      DATA (NHEL(IHEL,  4),IHEL=1,5) / +1, +1, -1, +1, +1/
      DATA (NHEL(IHEL,  5),IHEL=1,5) / +1, +1, +1, -1, +1/
      DATA (NHEL(IHEL,  6),IHEL=1,5) / +1, +1, +1, +1, -1/
      DATA (NHEL(IHEL,  7),IHEL=1,5) / -1, -1, +1, +1, +1/
      DATA (NHEL(IHEL,  8),IHEL=1,5) / -1, +1, -1, +1, +1/
      DATA (NHEL(IHEL,  9),IHEL=1,5) / -1, +1, +1, -1, +1/
      DATA (NHEL(IHEL, 10),IHEL=1,5) / -1, +1, +1, +1, -1/
      DATA (NHEL(IHEL, 11),IHEL=1,5) / +1, -1, -1, +1, +1/
      DATA (NHEL(IHEL, 12),IHEL=1,5) / +1, -1, +1, -1, +1/
      DATA (NHEL(IHEL, 13),IHEL=1,5) / +1, -1, +1, +1, -1/
      DATA (NHEL(IHEL, 14),IHEL=1,5) / +1, +1, -1, -1, +1/
      DATA (NHEL(IHEL, 15),IHEL=1,5) / +1, +1, -1, +1, -1/
      DATA (NHEL(IHEL, 16),IHEL=1,5) / +1, +1, +1, -1, -1/
      save NHEL
!$omp threadprivate(NHEL)
