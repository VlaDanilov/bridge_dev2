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
 
      integer:: colourchoice
!--- 'colourchoice' allows calculation by colour structure
!--- For Gflag=.true. [QQGG, QQGGG processes]
!--- 1) Only leading colour ( NCF . N )
!--- 2) Only sub-leading ( NCF . 1/N )
!--- 3) Only sub-sub-leading ( NCF . 1/N**3 )  [QQGGG only]
!--- 0) The total
!--- For Qflag=.true. [QQBQQB process]
!--- 1) Only leading colour ( NCF . 1 )
!--- 2) Only sub-leading ( NCF . 1/N ) [Identical quarks only]
!--- 0) The total
      common/ColC/colourchoice
!$omp threadprivate(/ColC/)
