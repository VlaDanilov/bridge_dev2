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
 
      function Rgen(pjet,i,p,j)
      implicit none
      include 'types.f'
      real(dp):: Rgen
c--- Calculate the angular separation between pjet(i) and p(j)
c--- This routine is a generalization of R.f: Rgen(p,i,p,j) == R(p,i,j) 
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: pjet(mxpart,4),p(mxpart,4),
     & r1,r2,dely,delphi,ei,ej
      integer:: i,j

      ei=sqrt(pjet(i,1)**2+pjet(i,2)**2+pjet(i,3)**2)
      ej=sqrt(p(j,1)**2+p(j,2)**2+p(j,3)**2)

      r1= (ei+pjet(i,3))*(ej-p(j,3))/
     &     ((ej+p(j,3))*(ei-pjet(i,3)))
      dely=0.5_dp*log(r1)

      r2= (pjet(i,1)*p(j,1)+pjet(i,2)*p(j,2))
     &     /sqrt((pjet(i,1)**2+pjet(i,2)**2)*(p(j,1)**2+p(j,2)**2))
      if (r2 > +0.9999999_dp) r2=+1._dp
      if (r2 < -0.9999999_dp) r2=-1._dp
      delphi=acos(r2)

      Rgen=sqrt(dely**2+delphi**2)
      
      return
      end
      
