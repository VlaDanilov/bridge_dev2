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
 
      subroutine setupscet(nprocabove)
      implicit none
c---- Routine to setup variables for performing SCET calculation;
c---- inspects value of nproc (passed via common block) and
c---- returns nprocabove, the corresponding process with one
c---- additional jet;  stops with error message if appropriate.
      include 'types.f'
      include 'nproc.f'
      include 'taucut.f'
      integer:: nprocabove

      if     ((nproc == 1) .or. (nproc == 6)
     &   .or. (nproc == 31) .or. (nproc == 32)) then
        nprocabove=nproc+10
        ntau=0        
      elseif ((nproc == 11) .or. (nproc == 16)) then
        nprocabove=nproc+11
        ntau=1       
      elseif (nproc == 41) then
        nprocabove=44
        ntau=1 
      elseif (nproc == 42) then
        nprocabove=46
        ntau=1 
      elseif ((nproc == 91) .or. (nproc == 96)) then
        nprocabove=nproc+519
        ntau=0
      elseif ((nproc == 92) .or. (nproc == 97)) then
        nprocabove=nproc+519
        ntau=0
      elseif ((nproc == 93) .or. (nproc == 98)) then
        nprocabove=nproc+519
        ntau=0
      elseif ((nproc == 94) .or. (nproc == 99)) then
        nprocabove=nproc+519
        ntau=0
      elseif (nproc == 101) then
        nprocabove=621
        ntau=0
      elseif (nproc == 104) then
        nprocabove=622
        ntau=0
      elseif (nproc == 106) then
        nprocabove=623
        ntau=0
      elseif (nproc == 110) then
        nprocabove=620
        ntau=0
      elseif (nproc == 111) then
        nprocabove=203
        ntau=0
      elseif (nproc == 112) then
        nprocabove=204
        ntau=0
      elseif (nproc == 119) then
        nprocabove=210
        ntau=0
      elseif (nproc == 120) then
        nprocabove=205
        ntau=0
      elseif (nproc == 203) then
        nprocabove=272
        ntau=1
      elseif (nproc == 204) then
        nprocabove=271
        ntau=1
      elseif (nproc == 210) then
        nprocabove=270
        ntau=1
      elseif ((nproc == 285)) then
        nprocabove=nproc+1
        ntau=0
      elseif (nproc == 300) then
        nprocabove=302
        ntau=0 
      elseif (nproc == 305) then
        nprocabove=307
        ntau=0
      elseif (nproc ==290 ) then
        nprocabove=292
        ntau=0
      elseif (nproc ==295 ) then
        nprocabove=297
        ntau=0
      elseif (nproc == 3000) then
        nprocabove=3002
        ntau=0
      else
        write(6,*) 'This process cannot be computed at NNLO'
        write(6,*) 'or at NLO using SCET.'
        stop
      endif
      
      return
      end
      
