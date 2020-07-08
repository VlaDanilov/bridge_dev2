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
 
      subroutine hwwdecay(p,j1,j2,j3,j4,msq)
      implicit none
      include 'types.f'
************************************************************************
*     Author: R.K. Ellis, September 2012                               *
*                                                                      *
*     matrix element squared for the process of                        *
*     Higgs decay  H --> WW --> ne(j1)+e+(j2)+e-(j3)+nu~(j4)           *
************************************************************************
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ewcouple.f'
      integer:: j1,j2,j3,j4
      real(dp):: p(mxpart,4),msq,dot,hdecay,s12,s13,s24,s34
      
      s12=2._dp*dot(p,j1,j2)
      s13=2._dp*dot(p,j1,j3)
      s24=2._dp*dot(p,j2,j4)
      s34=2._dp*dot(p,j3,j4)

      hdecay=gwsq**3*wmass**2*s13*s24
      hdecay=hdecay
     &  /(((s12-wmass**2)**2+(wmass*wwidth)**2)
     &   *((s34-wmass**2)**2+(wmass*wwidth)**2))
      msq=hdecay
      
      return
      end
      
      
