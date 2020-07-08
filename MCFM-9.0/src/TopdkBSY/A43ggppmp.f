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
 
      function A43ggppmp(e1,p2,p3,e4,za,zb,zab,zba)
      implicit none
      include 'types.f'
      complex(dp):: A43ggppmp
      
C-----Authors: John Campbell and Keith Ellis, November 2011
C---- taken from arXiv:1101.5947 [hep-ph], Eq. (5,53)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'zprods_decl.f'
      include 'zabprods_decl.f'
      include 'sprods_com.f'
      include 'epinv.f'
      include 'scale.f'
      real(dp):: mt2,beta,xlog
      complex(dp):: BSYALggppmphp,BSYARggppmphp,BSYALslggpppm
      complex(dp):: BSYA0ggppmp,lnrat,Vslc
      integer:: e1,p2,p3,e4,e1p,e4p
       
      mt2=mt**2
      beta=sqrt(1d0-4d0*mt2/s(p2,p3))
      xlog=log((1d0-beta)/(1d0+beta))/beta
C---- arXiv:1101.5947 [hep-ph], Eq. (52)
c----   with an extra minus sign to account for A4;3 -> -Vslc, c.f. Eq. (53))
      Vslc=-cplx1(epinv)*
     & (+s(1,p2)/s(p2,p3)*(lnrat(musq,-s(1,p2))+lnrat(mt2,-s(1,p2)))
     &  +s(1,p3)/s(p2,p3)*(lnrat(musq,-s(1,p3))+lnrat(mt2,-s(1,p3)))
     & +lnrat(musq,-s(p2,p3))+(1d0-2d0*mt2/s(p2,p3))*xlog)

c--- perform the swaps 5 <-> 7 and 6 <-> 8     
      e1p=12-e1
      e4p=14-e4
      
      A43ggppmp=(
     &  Vslc*(BSYA0ggppmp(e1,p2,p3,e4,za,zb,zab,zba)
     &       -BSYA0ggppmp(e1p,p3,p2,e4p,zb,za,zba,zab))
     & +BSYALggppmphp(e1,p2,p3,e4,za,zb,zab,zba)
     & -BSYALggppmphp(e1p,p3,p2,e4p,zb,za,zba,zab)
     & +BSYARggppmphp(e1,p2,p3,e4,za,zb,zab,zba)
     & -BSYARggppmphp(e1p,p3,p2,e4p,zb,za,zba,zab)
     & +BSYALslggpppm(e1,p2,e4,p3,za,zb,zab,zba)
     & -BSYALslggpppm(e1p,p3,e4p,p2,zb,za,zba,zab))

      return
      end

