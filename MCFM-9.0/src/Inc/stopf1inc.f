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
 
      complex(dp):: lc,ls,lVc,lVs,lp
     &,LsA,LsB1,LsB2,tr1Xfc,tr1Xfs,tr2fu,tr3Xc
     &,tr3c00fs,tr3c001fs,tr3c002fs,tr3Xs,tr3s00ft,tr3s001ft
     &,tr3s002ft,tr4Xc,tr4Xs,tr5Xc,tr5Xs,B0csf,B0cgsf
     &,lRc1,lRc2,lRs1,lRs2,lRcs,BfunX
      common/stopf1inc/ lc,ls,lVc,lVs,lp
     &,LsA,LsB1,LsB2,tr1Xfc,tr1Xfs,tr2fu,tr3Xc
     &,tr3c00fs,tr3c001fs,tr3c002fs,tr3Xs,tr3s00ft,tr3s001ft
     &,tr3s002ft,tr4Xc,tr4Xs,tr5Xc,tr5Xs,B0csf,B0cgsf
     &,lRc1,lRc2,lRs1,lRs2,lRcs,BfunX
!$omp threadprivate(/stopf1inc/)
      include 'sck.f'

