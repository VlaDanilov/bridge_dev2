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
 
      subroutine setnotag()
      implicit none
      include 'types.f'
      
      include 'removebr.f'
      include 'nproc.f'
      include 'notag.f'
c--- this routine sets the value of "notag", the number of jets
c--- that may be safely ignored without affecting finiteness of result;
c--- the minimum number of jets allowed by the code is equal to
c---    nqcdjets - notag
c--- where nqcdjets is itself process-dependent

c--- Modifying this routine to allowe *larger* values of notag
c--- than the defaults below (or >0 for processes not listed here)
c--- should be done with care

c--- Modifying this routine to allow *smaller* values of notag,
c--- i.e. stricter constraints on the number of jets observed,
c--- should not cause problems

      if     ((nproc == 62) .or. (nproc == 63)
     &    .or.(nproc == 64) .or. (nproc == 65)) then
c---- WW production, hadronic W decays
        notag=2

      elseif(((nproc == 161) .or. (nproc == 162)
     &    .or.(nproc == 163) .or. (nproc == 166)
     &    .or.(nproc == 167) .or. (nproc == 168))
     &   .and.(removebr)) then
c---- 5-flavor t-channel single top, top BR removed
c---    the calculation is inclusive of all additional jets; 
c---    can set notag=0 to explicitly require an additional jet at LO
        notag=1

      elseif ((nproc == 231) .or. (nproc == 232)
     &    .or.(nproc == 233) .or. (nproc == 234)
     &    .or.(nproc == 235) .or. (nproc == 236)
     &    .or.(nproc == 237) .or. (nproc == 238)
     &    .or.(nproc == 239) .or. (nproc == 240)) then
c---- 4-flavor t-channel single top
c---    the calculation requires the presence of two light jets that
c---    are present; to compute an inclusive cross section, one
c---    can set notag=1, or use jet cuts that have no effect
c---    (for the studies in arXiv:1204.1513, FERMILAB-PUB-12-078-T
c---      we have set notag=1)    
        notag=0

      elseif (nproc == 280) then
c---- direct photon production, presence of jet not required
        notag=1

      elseif ((nproc == 503) .or. (nproc == 506)
     &   .or. (nproc == 513) .or. (nproc == 516)) then
c---- ttW production, hadronic W in top decay
        notag=2

      elseif ((nproc == 532) .or. (nproc == 533)) then
c---- ttZ production, Z-> bb~ and hadronic W in top decay
        notag=4

      elseif ((nproc == 564) .or. (nproc == 567)) then
c---- tZ production with top decay
        notag=2
      endif
      
      return
      end
      
