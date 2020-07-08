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
 
      function A53(j1,j2,j3,j4,j5,za,zb,mtsq)
!  Amplitudes taken from Appendix IV of
!  Z.~Bern, L.~J.~Dixon and D.~A.~Kosower,
!  %``One loop amplitudes for e+ e- to four partons,''
!  Nucl.\ Phys.\ B {\bf 513}, 3 (1998)
!  [hep-ph/9708239].
!  Eq.(IV.11)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      complex(dp):: A53,L1,kinstructure,F1anom
      integer:: j1,j2,j3,j4,j5
      real(dp):: s45,s12,mtsq,musq
      logical,parameter::largemassexp=.false.
      logical,parameter::tandbtogether=.true.
      s12=s(j1,j2)
      s45=s(j1,j2)+s(j2,j3)+s(j3,j1)
      kinstructure=two*zb(j5,j3)*zb(j3,j1)*za(j2,j4)/s45
      A53=czip
      if (tandbtogether) A53=-kinstructure*L1(-s12,-s45)/two/s45
      if (mtsq == 0._dp) then
         A53=+kinstructure*L1(-s12,-s45)/two/s45
      else
         if (largemassexp) then
            A53=A53+kinstructure/(24._dp*mtsq)
     &      *(1._dp+(2._dp*s12+s45)/15._dp/mtsq*zip
     &      +(2._dp*s12*s45+3._dp*s12**2+s45**2)/140._dp/mtsq**2*zip)
         else
!     musq is irrelevant, set to some value
            musq=abs(s45)
            A53=A53+kinstructure*F1anom(s12,s45,mtsq,musq)
         endif
      endif
      return  
      end
