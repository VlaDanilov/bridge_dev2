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
 
      function Finite0x2qqgagas(s,t)
!     Results taken from hep-ph/0201274
! \bibitem{Anastasiou:2002zn} 
!  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
!  %``Two loop QED and QCD corrections to massless fermion boson scattering,''
!  Nucl.\ Phys.\ B {\bf 629}, 255 (2002)
!  [hep-ph/0201274].
      implicit none
      include 'types.f'
      real(dp)::Finite0x2qqgagas
      include 'constants.f'
      include 'nf.f'
      include 'scale.f'
      real(dp)::s,t,u,x,y,z,zinv,Lx,Ly,Ls,
     & Li4x,Li4y,Li4z,Li4zinv,
     & Li3x,Li3y,
     & Li2x,Li2y,
     & Li4,Li3,Li2,
     & AGTYAs,xAs,AGTYBs,xBs,AGTYD2s,xD2s,AGTYE3s,xE3s
      real(dp),parameter::sumq=11._dp/9._dp
      u=-s-t
      x=-t/s
      y=-u/s
      z=-y/x
      zinv=-x/y

      Lx=log(x)
      Ly=log(y)
      Ls=log(s/musq)

      Li4x=Li4(x)
      Li4y=Li4(y)
      Li4z=Li4(z)
      Li4zinv=Li4(zinv)

      Li3x=Li3(x)
      Li3y=Li3(y)

      Li2x=Li2(x)
      Li2y=Li2(y)
      xAs=AGTYAs(s,t,u,Lx,Ly,Li2x,Li2y,Li3x,Li3y,Li4x,Li4y,
     & Li4z,Li4zinv)
      xBs=AGTYBs(s,t,u,Lx,Ly,Ls,Li2x,Li3x,Li4x,Li2y,Li3y,Li4y,
     & Li4z,Li4zinv)
      xD2s=AGTYD2s(s,t,u,Lx,Ly,Ls,Li2x,Li3x,Li4x,Li2y,Li3y,Li4y,
     & Li4z,Li4zinv)
      xE3s=AGTYE3s(t,u,Lx,Ly,Ls)

      Finite0x2qqgagas=two*xn
     & *(sumq*TR*CF*xAs+CF**2*xBs+CF*CA*xD2s+nf*CF*xE3s)
      return
      end





