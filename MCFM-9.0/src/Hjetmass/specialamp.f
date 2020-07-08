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
 
      subroutine specialamp(amp)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'qqbggnames.f'
      real(dp)::p(mxpart,4)
      complex(dp)::amp(2,2)
      complex(dp)::bubints(Nbubints),triints(Ntriints),boxints(Nboxints)
      integer::ord
      include 'specialmom.f'
      write(6,*) 'p1',p(1,4),p(1,1),p(1,2),p(1,3)
      write(6,*) 'p2',p(2,4),p(2,1),p(2,2),p(2,3)
      write(6,*) 'p3',p(3,4),p(3,1),p(3,2),p(3,3)
      write(6,*) 'p4',p(4,4),p(4,1),p(4,2),p(4,3)
      call spinoru(4,p,za,zb)
      ord=0
      call hjetmass_qqbgg_bubints_init_dp(1,2,3,4,bubints,ord,s)
      call hjetmass_qqbgg_triints_init_dp(1,2,3,4,triints,ord,s)
      call hjetmass_qqbgg_boxints_init_dp(1,2,3,4,boxints,ord,s)
      call specialpoint(boxints,triints,bubints,amp)
!      write(6,*) 'amp(1,1)',amp(1,1)/2._dp,abs(amp(1,1))/2._dp
!      write(6,*) 'amp(1,2)',amp(1,2)/2._dp,abs(amp(1,2))/2._dp
!      write(6,*) 'amp(2,1)',amp(2,1)/2._dp,abs(amp(2,1))/2._dp
!      write(6,*) 'amp(2,2)',amp(2,2)/2._dp,abs(amp(2,2))/2._dp
      return
      end
