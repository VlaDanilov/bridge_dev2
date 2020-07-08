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
 
      subroutine gen4_intf(r,p,wt,*)
      implicit none
      include 'types.f'
c--- Generates phase space for gg -> VV | gg-> H -> VV interference processes.
c--- If the mass of the Higgs boson is above (or close) to the minimum
c--- invariant mass required by the cuts (m34min and m56min),
c--- the phase space is generated using a B.W. for the Higgs boson (gen4h.f);
c--- if the Higgs mass is below this threshold then the phase space is
c--- generated using the standard VV routine gen4.f
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'mxdim.f'
      include 'limits.f'
      include 'first.f'
      logical, save::useHiggsBW
      real(dp):: r(mxdim),p(mxpart,4),wt,threshold

c--- initialization on first call      
      if (first) then
        threshold=wsqmin+bbsqmin        
c--- generate using B.W. if Higgs mass above or close to threshold      
        if (hmass > threshold-5._dp*hwidth) then
c--- otherwise generate without B.W.
          useHiggsBW=.true.
        else
          useHiggsBW=.false.
        endif
        write(6,*)
        write(6,*) 'gen4_intf: useHiggsBW = ',useHiggsBW
        first=.false.
      endif

c--- switch between gen4h and gen4
      if (useHiggsBW) then
        call gen4h(r,p,wt,*999)
      else
        call gen4(r,p,wt,*999)
      endif
      
      return

 999  wt=0._dp
      return 1

      end
      
