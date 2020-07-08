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
 
      function DOTKS(I,J)
      implicit none
      include 'types.f'
      integer:: I,J
      real(dp):: DOTKS
      real(dp):: PLAB(4,10)
      COMMON/MOM/PLAB
!$omp threadprivate(/MOM/)  

         DOTKS=PLAB(4,I)*PLAB(4,J)-PLAB(3,I)*PLAB(3,J)
     1        -PLAB(2,I)*PLAB(2,J)-PLAB(1,I)*PLAB(1,J)
      RETURN
      END

