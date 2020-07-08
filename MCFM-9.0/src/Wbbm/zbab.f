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
 
      function zbab(j1,p1,p2,j2)
      implicit none
      include 'types.f'
      complex(dp):: zbab
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'momwbbm.f'
      integer:: j1,j2

      real(dp):: p1(4),p2(4),E,px,py,pz
      complex(dp):: Up(4),Vp(4),zp1(4),zp2(4),spstrng2,rtpp

C---Calculate spinor 1
      E=mom(j1,4)      
      px=+mom(j1,3)      
      py=-mom(j1,2)      
      pz=+mom(j1,1)      
      rtpp=sqrt(cplx1(E+pz))
      Up(1)=czip
      Up(2)=czip
      Up(3)=rtpp
      Up(4)=cplx2(px,-py)/rtpp

C---Calculate spinor 2
      E=mom(j2,4)      
      px=+mom(j2,3)      
      py=-mom(j2,2)      
      pz=+mom(j2,1)      
      rtpp=sqrt(cplx1(E+pz))
      Vp(1)=czip
      Vp(2)=czip
      Vp(3)=cplx2(px,-py)/rtpp
      Vp(4)=-rtpp

      zp1(1)=cplx1(p1(1))
      zp1(2)=cplx1(p1(2))
      zp1(3)=cplx1(p1(3))
      zp1(4)=cplx1(p1(4))

      zp2(1)=cplx1(p2(1))
      zp2(2)=cplx1(p2(2))
      zp2(3)=cplx1(p2(3))
      zp2(4)=cplx1(p2(4))

      zbab=spstrng2(Up,zp1,zp2,Vp)
      return
      end
