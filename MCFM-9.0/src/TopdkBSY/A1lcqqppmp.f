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
 
      function A1lcqqppmp(e1,p2,p3,e4,za,zb,zab,zba)
      implicit none
      include 'types.f'
      complex(dp):: A1lcqqppmp
      
C-----Authors: John Campbell and Keith Ellis, November 2011
C---- taken from arXiv:1101.5947 [hep-ph], Eq. (55)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'zabprods_decl.f'
      include 'epinv.f'
      include 'scale.f'
      real(dp):: mt2
      complex(dp):: BSYA0qqppmp,BSYA1lcqqppmp,Vlc,lnrat
      integer:: e1,p2,p3,e4
      mt2=mt**2 
c--- NB: added a minus sign in front of the double pole compared to Eq. (55)
      Vlc=-cplx1(epinv**2)+cplx1(epinv)*(8d0/3d0
     & -(lnrat(musq,-s(1,p2))+lnrat(mt2,-s(1,p2))))
      A1lcqqppmp=Vlc*BSYA0qqppmp(e1,p2,p3,e4,za,zb,zab,zba)
     &              +BSYA1lcqqppmp(e1,p2,p3,e4,za,zb,zab,zba)  
      
      return
      end
