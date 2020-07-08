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
 
      function assemble(order,taucut,beama0,beamb0,beama1,beamb1,
     & beama2,beamb2,soft1,soft2,hard)
      implicit none
      include 'types.f'
c---- Given beam, soft and hard functions, calculates
c---- O(alphas^order) correction after all convolutions
      include 'constants.f'
      include 'nf.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'scet_const.f'
      include 'kpart.f'
      integer:: order
      real(dp):: assemble,soft1(-1:1),soft2(-1:3),hard(2),
     & beama0,beamb0,
     & beama1(-1:1),beamb1(-1:1),
     & beama2(-1:3),beamb2(-1:3),
     & Lt1cut,full(-1:3),taucut

      Lt1cut = log(taucut/scale)

      if (coeffonly) then
        assemble=zip
      else
        assemble=beama0*beamb0
      endif
      if ( (order ==1) .or.
     &    ((order ==2) .and. (coeffonly .eqv. .false.)) ) then
      full(3)=zip
      full(2)=zip
      full(1)=beama1(1)*beamb0 + beamb1(1)*beama0 + beama0*beamb0*
     & soft1(1)

      full(0)=beama1(0)*beamb0 + beamb1(0)*beama0 + beama0*beamb0*
     & soft1(0)

      full(-1)=beama1(-1)*beamb0 + beamb1(-1)*beama0 + beama0*beamb0*
     & soft1(-1) + beama0*beamb0*hard(1)
      assemble=assemble+ason2pi*(
     & +full(-1)
     & +full(0)*Lt1cut
     & +full(1)*Lt1cut**2/two)
      endif
      if (order > 1) then
      full(3)=beama2(3)*beamb0 + beamb2(3)*beama0 + beama1(1)*beamb1(1)
     &  + beama1(1)*beamb0*soft1(1) + beamb1(1)*beama0*soft1(1) + 
     & beama0*beamb0*soft2(3)

      full(2)=beama2(2)*beamb0 + beamb2(2)*beama0 + 3._dp/2._dp*beama1(0)
     & *beamb1(1) + 3._dp/2._dp*beama1(0)*beamb0*soft1(1) + 3._dp/2._dp*
     & beama1(1)*beamb1(0) + 3._dp/2._dp*beama1(1)*beamb0*soft1(0) + 3._dp/
     & 2._dp*beamb1(0)*beama0*soft1(1) + 3._dp/2._dp*beamb1(1)*beama0*
     & soft1(0) + beama0*beamb0*soft2(2)

      full(1)=beama2(1)*beamb0 + beamb2(1)*beama0 + beama1(-1)*beamb1(1
     & ) + beama1(-1)*beamb0*soft1(1) + 2._dp*beama1(0)*beamb1(0) + 2._dp
     & *beama1(0)*beamb0*soft1(0) + beama1(1)*beamb1(-1) - 2._dp*beama1(
     & 1)*beamb1(1)*zeta2 + beama1(1)*beamb0*soft1(-1) - 2._dp*beama1(1)
     & *beamb0*soft1(1)*zeta2 + beama1(1)*beamb0*hard(1) + beamb1(-1)*
     & beama0*soft1(1) + 2._dp*beamb1(0)*beama0*soft1(0) + beamb1(1)*
     & beama0*soft1(-1) - 2._dp*beamb1(1)*beama0*soft1(1)*zeta2 + 
     & beamb1(1)*beama0*hard(1) + beama0*beamb0*soft1(1)*hard(1) + 
     & beama0*beamb0*soft2(1)

      full(0)=beama2(0)*beamb0 + beamb2(0)*beama0 + beama1(-1)*beamb1(0
     & ) + beama1(-1)*beamb0*soft1(0) + beama1(0)*beamb1(-1) - beama1(0
     & )*beamb1(1)*zeta2 + beama1(0)*beamb0*soft1(-1) - beama1(0)*
     & beamb0*soft1(1)*zeta2 + beama1(0)*beamb0*hard(1) - beama1(1)*
     & beamb1(0)*zeta2 + 2._dp*beama1(1)*beamb1(1)*zeta3 - beama1(1)*
     & beamb0*soft1(0)*zeta2 + 2._dp*beama1(1)*beamb0*soft1(1)*zeta3 + 
     & beamb1(-1)*beama0*soft1(0) + beamb1(0)*beama0*soft1(-1) - 
     & beamb1(0)*beama0*soft1(1)*zeta2 + beamb1(0)*beama0*hard(1) - 
     & beamb1(1)*beama0*soft1(0)*zeta2 + 2._dp*beamb1(1)*beama0*soft1(1)
     & *zeta3 + beama0*beamb0*soft1(0)*hard(1) + beama0*beamb0*soft2(0)

      full(-1)=beama2(-1)*beamb0 + beamb2(-1)*beama0 + beama1(-1)*
     & beamb1(-1) + beama1(-1)*beamb0*soft1(-1) + beama1(-1)*beamb0*
     & hard(1) - beama1(0)*beamb1(0)*zeta2 + beama1(0)*beamb1(1)*zeta3
     &  - beama1(0)*beamb0*soft1(0)*zeta2 + beama1(0)*beamb0*soft1(1)*
     & zeta3 + beama1(1)*beamb1(0)*zeta3 - 1._dp/10._dp*beama1(1)*beamb1(
     & 1)*zeta2**2 + beama1(1)*beamb0*soft1(0)*zeta3 - 1._dp/10._dp*
     & beama1(1)*beamb0*soft1(1)*zeta2**2 + beamb1(-1)*beama0*soft1(-1)
     &  + beamb1(-1)*beama0*hard(1) - beamb1(0)*beama0*soft1(0)*zeta2
     &  + beamb1(0)*beama0*soft1(1)*zeta3 + beamb1(1)*beama0*soft1(0)*
     & zeta3 - 1._dp/10._dp*beamb1(1)*beama0*soft1(1)*zeta2**2 + beama0*
     & beamb0*soft1(-1)*hard(1) + beama0*beamb0*soft2(-1) + beama0*
     & beamb0*hard(2)
      assemble=assemble+ason2pi**2*(
     & +full(-1)
     & +full(0)*Lt1cut
     & +full(1)*Lt1cut**2/two
     & +full(2)*Lt1cut**3/three
     & +full(3)*Lt1cut**4/four)
      endif

      return
      end
