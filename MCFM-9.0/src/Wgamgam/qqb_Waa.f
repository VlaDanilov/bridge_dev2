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
 
      subroutine qqb_Waa(p,msq)
      implicit none
      include 'types.f'
      
c---  matrix element squared for the process
c---  q(p1)+qbar(p2)+W(l(p3)+lbar(p4))+gamma(p5)+gamma(p6)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_com.f'
      include 'nwz.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'sprods_com.f' 
      real(dp):: udb,dbu,dub,ubd,msqWaa,fac,stat
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf),msqWaa_frm
      complex(dp):: temp
      parameter(stat=0.5_dp)
      msq(:,:)=0._dp 
      fac=stat*aveqq*xn*4._dp*esq**2*gwsq**2

      call spinoruorz(6,p,za,zb)
c      call spinoru(6,p,za,zb)
      
      if (nwz == 1) then
        udb=fac*msqWaa(1,2,3,4,5,6,za,zb)
        dbu=fac*msqWaa(2,1,3,4,5,6,za,zb)
        msq(2,-1)=udb
        msq(4,-3)=udb
        msq(-1,2)=dbu
        msq(-3,4)=dbu
      elseif (nwz == -1) then
        dub=fac*msqWaa(2,1,4,3,5,6,zb,za)
        ubd=fac*msqWaa(1,2,4,3,5,6,zb,za)
        msq(1,-2)=dub
        msq(3,-4)=dub
        msq(-2,1)=ubd
        msq(-4,3)=ubd
      endif
      
      return
      end

