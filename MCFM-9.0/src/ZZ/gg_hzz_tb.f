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
 
      subroutine gg_hZZ_tb(p,msq)
        use mod_qcdloop_c
      implicit none
      include 'types.f'
      
c--- Author: J. M. Campbell, September 2013
c--- Matrix element squared for gg -> H -> ZZ signal process
c--- The exact result for massive bottom and top quark loops is included
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'interference.f'
      integer:: h1,h2,h34,h56
      real(dp):: p(mxpart,4),msq(fn:nf,fn:nf),msqgg,fac,
     & pswap(mxpart,4),oprat
      complex(dp):: ggH_bquark(2,2,2,2),ggH_tquark(2,2,2,2),Ahiggs,
     & ggH_bquark_swap(2,2,2,2),ggH_tquark_swap(2,2,2,2),Ahiggs_swap

      msq(:,:)=0._dp
      
      call getggHZZamps(p,ggH_bquark,ggH_tquark)
      
      if (interference) then
c--- for interference, compute amplitudes after 4<->6 swap
       pswap(1,:)=p(1,:)
       pswap(2,:)=p(2,:)
       pswap(3,:)=p(3,:)
       pswap(4,:)=p(6,:)
       pswap(5,:)=p(5,:)
       pswap(6,:)=p(4,:)
       call getggHZZamps(pswap,ggH_bquark_swap,ggH_tquark_swap)
      endif
      
      msqgg=0._dp
      do h1=1,2
      h2=h1
      do h34=1,2
      do h56=1,2
      
c--- compute total Higgs amplitude
      AHiggs=
     &  +ggH_bquark(h1,h2,h34,h56)
     &  +ggH_tquark(h1,h2,h34,h56)
     
      if (interference .eqv. .false.) then
c--- normal case
        msqgg=msqgg+abs(AHiggs)**2
      else
c--- with interference
        AHiggs_swap=
     &  +ggH_bquark_swap(h1,h2,h34,h56)
     &  +ggH_tquark_swap(h1,h2,h34,h56)
        if (h34 == h56) then
          oprat=1._dp-two*real(conjg(AHiggs)*AHiggs_swap)
     &                 /(abs(AHiggs)**2+abs(AHiggs_swap)**2)
        else
          oprat=1._dp
        endif
        if (bw34_56) then
          msqgg=msqgg+two*abs(AHiggs)**2*oprat
        else
          msqgg=msqgg+two*abs(AHiggs_swap)**2*oprat
        endif
      endif
      enddo
      enddo
      enddo

c--- overall factor extracted (c.f. getggHZZamps.f)
      fac=avegg*V*(four*esq*gsq/(16._dp*pisq)*esq)**2
      
      msq(0,0)=msqgg*fac*vsymfact

      return
      end
      
      
