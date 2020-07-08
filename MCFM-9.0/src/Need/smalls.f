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
 
      subroutine smalls(s,npart,*)
      implicit none
      include 'types.f'
c    cut if radiated parton too close
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'cutoff.f'
      include 'kprocess.f'
      integer:: npart,j,k
      real(dp):: s(mxpart,mxpart)
      
      do j=3,npart+2
      if ((-s(1,j) < cutoff_s) .or. (-s(2,j) < cutoff_s)) return 1
        do k=j+1,npart+2
        if (s(j,k) < cutoff_s) return 1
        enddo
      enddo
      
      return
      
      
      
      if (npart == 2) then
      if ( 
     &      (-s(1,4) < cutoff_s)
     & .or. (-s(2,4) < cutoff_s)
     & .or. (+s(3,4) < cutoff_s)
     & .or. (-s(1,3) < cutoff_s)
     & .or. (-s(2,3) < cutoff_s)
     & ) return 1
      
      elseif (npart == 3) then
      if ( 
     &      (-s(1,5) < cutoff_s)
     & .or. (-s(2,5) < cutoff_s)
     & .or. (-s(1,4) < cutoff_s)
     & .or. (-s(2,4) < cutoff_s)
     & .or. (+s(4,5) < cutoff_s)
     & .or. (+s(3,4) < cutoff_s)
     & .or. (+s(3,5) < cutoff_s)
     & .or. (-s(1,3) < cutoff_s)
     & .or. (-s(2,3) < cutoff_s)
     & ) return 1

      elseif (npart == 4) then
      if ( 
     &      (-s(1,6) < cutoff_s)
     & .or. (-s(2,6) < cutoff_s)
     & .or. (-s(1,5) < cutoff_s)
     & .or. (-s(2,5) < cutoff_s)
     & .or. (-s(1,4) < cutoff_s)
     & .or. (-s(2,4) < cutoff_s)
     & .or. (-s(1,3) < cutoff_s)
     & .or. (-s(2,3) < cutoff_s)
     & .or. (+s(3,4) < cutoff_s)
     & .or. (+s(3,5) < cutoff_s)
     & .or. (+s(3,6) < cutoff_s)
     & .or. (+s(4,5) < cutoff_s)
     & .or. (+s(4,6) < cutoff_s)
     & .or. (+s(5,6) < cutoff_s)
     & ) return 1
c      if (  (kcase==kqq_tbg) .or. (kcase==kqqtbgg)
c     & .or. (kcase==kepem3j) .or. (kcase==kW_tndk)
c     & .or. (kcase==kZ_tjet)) then
c      if ( 
c     &      (+s(3,4) < cutoff_s)
c     & .or. (+s(3,5) < cutoff_s)
c     & .or. (+s(3,6) < cutoff_s)
c     & .or. (+s(4,5) < cutoff_s)
c     & .or. (+s(4,6) < cutoff_s)
c     & ) return 1
c      endif
     
      elseif (npart == 5) then
      if ( 
     &      (-s(1,5) < cutoff_s)
     & .or. (-s(2,5) < cutoff_s)
     & .or. (-s(1,6) < cutoff_s)
     & .or. (-s(2,6) < cutoff_s)
     & .or. (-s(1,7) < cutoff_s)
     & .or. (-s(2,7) < cutoff_s)
     & .or. (-s(1,4) < cutoff_s)
     & .or. (-s(2,4) < cutoff_s)
     & .or. (-s(1,3) < cutoff_s)
     & .or. (-s(2,3) < cutoff_s)
     & .or. (+s(3,4) < cutoff_s)
     & .or. (+s(3,5) < cutoff_s)
     & .or. (+s(3,6) < cutoff_s)
     & .or. (+s(3,7) < cutoff_s)
     & .or. (+s(4,5) < cutoff_s)
     & .or. (+s(4,6) < cutoff_s)
     & .or. (+s(4,7) < cutoff_s)
     & .or. (+s(5,6) < cutoff_s)
     & .or. (+s(5,7) < cutoff_s)
     & .or. (+s(6,7) < cutoff_s)
     & ) return 1

      elseif (npart == 6) then
      if ( 
     &      (-s(1,5) < cutoff_s)
     & .or. (-s(2,5) < cutoff_s)
     & .or. (-s(1,6) < cutoff_s)
     & .or. (-s(2,6) < cutoff_s)
     & .or. (-s(1,7) < cutoff_s)
     & .or. (-s(2,7) < cutoff_s)
     & .or. (-s(1,8) < cutoff_s)
     & .or. (-s(2,8) < cutoff_s)
     & .or. (+s(5,6) < cutoff_s)
     & .or. (+s(5,7) < cutoff_s)
     & .or. (+s(5,8) < cutoff_s)
     & .or. (+s(6,7) < cutoff_s)
     & .or. (+s(6,8) < cutoff_s)
     & .or. (+s(7,8) < cutoff_s)
     & ) return 1
      
      elseif (npart == 7) then
      if ( 
     &      (-s(1,5) < cutoff_s)
     & .or. (-s(2,5) < cutoff_s)
     & .or. (-s(1,6) < cutoff_s)
     & .or. (-s(2,6) < cutoff_s)
     & .or. (-s(1,7) < cutoff_s)
     & .or. (-s(2,7) < cutoff_s)
     & .or. (-s(1,8) < cutoff_s)
     & .or. (-s(2,8) < cutoff_s)
     & .or. (-s(1,9) < cutoff_s)
     & .or. (-s(2,9) < cutoff_s)
     & .or. (+s(5,6) < cutoff_s)
     & .or. (+s(5,7) < cutoff_s)
     & .or. (+s(5,8) < cutoff_s)
     & .or. (+s(5,9) < cutoff_s)
     & .or. (+s(6,7) < cutoff_s)
     & .or. (+s(6,8) < cutoff_s)
     & .or. (+s(6,9) < cutoff_s)
     & .or. (+s(7,8) < cutoff_s)
     & .or. (+s(7,9) < cutoff_s)
     & .or. (+s(8,9) < cutoff_s)
     & ) return 1
      
      endif      
      
      return
      end
