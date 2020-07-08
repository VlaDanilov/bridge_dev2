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
 
      subroutine sumamp(coeff,scints,amp,str)
      implicit none
      include 'types.f'
c--- routine to multiply scalar integrals (in scints) by the computed
c--- coefficients (in coeff) and return the result in amp
c--- the 6-character string "str" is only used as output when checking primitives
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'Wbbmlabels.f'
      integer:: iep,j,k
      character*6 str
      complex(dp):: amp(-2:0)
      logical:: numcheck
      common/numcheck/numcheck
!$omp threadprivate(/numcheck/)

c--- multiply scalar integrals by coefficients
c---  NB: only need sum over finite pieces, poles handled separately
c      do iep=-2,0
      do iep=0,0
      amp(iep)=czip
      do j=1,4
      do k=1,20
      amp(iep)=amp(iep)+coeff(j,k)*scints(j,k,iep)
c      write(6,*) iep,j,k,coeff(j,k),scints(j,k,iep)
      enddo
      enddo      
      enddo     
c--- add purely rational term
      amp(0)=amp(0)+coeff(0,irat)

      if (numcheck) then
        write(6,89) str//'(-2) =',amp(-2) 
        write(6,89) str//'(-1) =',amp(-1) 
        write(6,89) str//'( 0) =',amp( 0)
        write(6,*)
      endif

      return
      
   89 format(SP,a13,2e20.11,' i')   
      
      end
      
