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
 
      function vetow_2gam(p)
       implicit none
      include 'types.f'
      logical:: vetow_2gam
c--- returns TRUE if the momenta passed into the function
c--- should be vetoed according to the current value of "ipsgen"
      include 'mxpart.f'
      include 'ipsgen.f'
      include 'masses.f'
      include 'kpart.f'
      real(dp):: p(mxpart,4),dot,s345,s346,s3456,xwid
c--- note: parameter "xwid" controls the number of widths away
c--- from the peak to generate according to a BW
      parameter(xwid=5._dp)
      
      vetow_2gam=.false.
      
c--- veto for fragmentation contribution
      if (kpart==kfrag) then
        if (ipsgen == 1) then
          s345=2._dp*(dot(p,3,4)+dot(p,3,5)+dot(p,4,5))
          if (abs(sqrt(s345)-wmass) < xwid*wwidth) then
            vetow_2gam=.true.
          endif
        endif
        return
      endif

c--- veto for everything else     
      if (ipsgen == 1) then
        s345=2._dp*(dot(p,3,4)+dot(p,3,5)+dot(p,4,5))
        if (abs(sqrt(s345)-wmass) < xwid*wwidth) then
          vetow_2gam=.true.
          return
        endif
      endif
      if ((ipsgen == 1) .or. (ipsgen == 3)) then
        s346=2._dp*(dot(p,3,4)+dot(p,3,6)+dot(p,4,6))
        if (abs(sqrt(s346)-wmass) < xwid*wwidth) then
          vetow_2gam=.true.
          return
        endif
      endif
      if ((ipsgen == 1) .or. (ipsgen == 3) .or. (ipsgen == 4))then
        s3456=2._dp*(dot(p,3,4)+dot(p,3,5)+dot(p,3,6)
     &            +dot(p,4,5)+dot(p,4,6)+dot(p,5,6))
        if (abs(sqrt(s3456)-wmass) < xwid*wwidth) then
          vetow_2gam=.true.
          return
        endif
      endif
      
      return
      end
      
