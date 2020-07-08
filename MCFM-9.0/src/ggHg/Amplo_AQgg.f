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
 
      subroutine Amplo_AQgg(p1,p2,p3,p4,ab,ba)
      implicit none
      include 'types.f'
c--- this routine is a wrapper to the new versions of the leading
c--- order amplitudes that are calculated in the same way as the
c--- new virtual amplitudes
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_com.f'
      integer:: p1,p2,p3,p4
      complex(dp):: ab(2,2,2),ba(2,2,2),
     & A0HAQggmppp,A0HAQggmpmm,A0HAQggmpmp,A0HAQggmppm
      
      ab(1,2,2)=A0HAQggmppp(p1,p2,p3,p4,za,zb)
      ab(1,1,1)=A0HAQggmpmm(p1,p2,p3,p4,za,zb)
      ab(1,1,2)=A0HAQggmpmp(p1,p2,p3,p4,za,zb)
      ab(1,2,1)=A0HAQggmppm(p1,p2,p3,p4,za,zb)

      ba(1,2,2)=A0HAQggmppp(p1,p2,p4,p3,za,zb)
      ba(1,1,1)=A0HAQggmpmm(p1,p2,p4,p3,za,zb)
      ba(1,1,2)=A0HAQggmppm(p1,p2,p4,p3,za,zb)
      ba(1,2,1)=A0HAQggmpmp(p1,p2,p4,p3,za,zb)

c--- Obtain opposite helicities through parity relation, c.f. Eq.(4.6) of DS 0906.0008
      ab(2,1,1)=-A0HAQggmppp(p1,p2,p3,p4,zb,za)
      ab(2,2,2)=-A0HAQggmpmm(p1,p2,p3,p4,zb,za)
      ab(2,2,1)=-A0HAQggmpmp(p1,p2,p3,p4,zb,za)
      ab(2,1,2)=-A0HAQggmppm(p1,p2,p3,p4,zb,za)

      ba(2,1,1)=-A0HAQggmppp(p1,p2,p4,p3,zb,za)
      ba(2,2,2)=-A0HAQggmpmm(p1,p2,p4,p3,zb,za)
      ba(2,2,1)=-A0HAQggmppm(p1,p2,p4,p3,zb,za)
      ba(2,1,2)=-A0HAQggmpmp(p1,p2,p4,p3,zb,za)

      return
      end
      
