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
 
      subroutine read_jetcuts(read_ptmin,read_etamin,read_etamax)
      implicit none
      include 'types.f'
      
      include 'clustering.f'
      include 'jetcuts.f'
      integer:: nargs
      character*72 jetcutsfile
      real(dp):: read_ptmin,read_etamin,read_etamax
      real(dp):: ptmin,etamax,ptmin_tev,etamax_tev,
     & ptmin_lhc,etamax_lhc
      logical:: useTevcuts,useLHCcuts
      
      read_ptmin=ptjetmin
      read_etamin=etajetmin
      read_etamax=etajetmax

      return
      end
      
      
      
