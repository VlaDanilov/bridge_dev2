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
 
!========= new style of routine in which integrated frag and dipoles are combined together 
!==== C Williams April 2013
      subroutine qqb_dirgam_frag_combo(p,p_phys,msq) 
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'ewcharge.f'
      include 'frag.f'
      include 'msqbits.f'
      include 'sprods_com.f'
      real(dp):: p(mxpart,4),p_phys(mxpart,4)
      real(dp):: msq(-nf:nf,-nf:nf)
      real(dp):: msq34_1(-nf:nf,-nf:nf),msqbits34_1(12)
      real(dp):: msq34_1_swap(-nf:nf,-nf:nf),msqbits34_1_swap(12)
      integer:: j,k,i
      real(dp):: virt_dip,xl,dot,fsq 
      real(dp):: aewo2pi,fi_gaq
      real(dp):: D(0:5),d_sum_q
      real(dp):: smallc,ss,tt,uu,smalld
      real(dp):: ga_ga,ag_ga,gq_gq,qg_gq,gg_gg,qa_gg,aq_gg
      real(dp):: fac
      common/D/D
!$omp threadprivate(/D/)
   
      aewo2pi=esq/(fourpi*twopi)
      fsq=frag_scale**2
      msq(:,:)=0._dp 

!======= integrated dipoles 
      xl=log(-two*dot(p_phys,1,3)/fsq)
      virt_dip=+aewo2pi*(fi_gaq(z_frag,p_phys,xl,3,1,2))
!======== fragmentation functions 
c---- Generate array D(j) corresponding to MCFM notation 0=gluon 1=down 2=up ....
      do i=0,5
         D(i)=0._dp
         if     (fragset == 'BFGset_I') then
            call get_frag(z_frag,fsq,1,i,D(i))   
         elseif (fragset == 'BFGsetII') then  
            call get_frag(z_frag,fsq,2,i,D(i))   
         elseif (fragset == 'GdRG__LO') then 
            call GGdR_frag(z_frag,i,D(i),0)
         else
            write(6,*) 'Unrecognized fragmentation set name: ',fragset
            stop        
         endif
      enddo

!============== matrix element functions 
      
      call qqb_2jet(p,msq34_1) 
      msqbits34_1(:)=msqbits(:) 
      call qqb_2jet_swap(p,msq34_1_swap)
!      call fill_dirgam_swap(7,msq34_1_swap)
      msqbits34_1_swap(:)=msqbits(:) 
!===== higher order pieces 
      call dotem(4,p,s)
      fac=gsq**2
     
      ss=s(1,2)
      tt=s(1,3)
      uu=s(2,3)
      qa_gg=+fac*aveqq*smallc(ss,tt,uu)*half
      aq_gg=+fac*aveqq*smallc(ss,tt,uu)*half

      gq_gq=-fac*aveqg*smallc(tt,ss,uu)
      ga_ga=-fac*aveqg*smallc(tt,uu,ss)
      qg_gq=-fac*aveqg*smallc(uu,ss,tt)
      ag_ga=-fac*aveqg*smallc(uu,tt,ss)

      gg_gg=fac*avegg*smalld(ss,tt,uu)*half



!========= pieces with dipoles and fragmentation 
!========== factors of two arise from pieces with both 7 and 8 (34 and 35) singular regions, pieces with only 0
      do j=-nf,nf 
      do k=-nf,nf

!=========
!========= quark -> photon fragmentation contributions
!=========
! quark-gluon
        if ((k == 0).and.(j.ne.0)) then
           msq(j,k)=(Q(j)**2*virt_dip+D(abs(j)))*msq34_1(j,k)
! gluon-quark
        elseif((j==0).and.(k.ne.0)) then 
           msq(j,k)=(Q(k)**2*virt_dip+D(abs(k)))*msq34_1(j,k)
! gluon-gluon
        elseif((j==0).and.(k==0)) then 
           msq(j,k)=((3._dp*Q(1)**2+2._dp*Q(2)**2)*virt_dip+d_sum_q(0))*msq34_1(j,k)*two
! quark-quark
        elseif((j>0).and.(k>0)) then 
           msq(j,k)=(Q(j)**2*virt_dip+D(j))*msq34_1(j,k)
     &             +(Q(k)**2*virt_dip+D(k))*msq34_1_swap(j,k)
! antiquark-antiquark
        elseif((j<0).and.(k<0)) then 
           msq(j,k)=(Q(j)**2*virt_dip+D(-j))*msq34_1(j,k)
     &             +(Q(k)**2*virt_dip+D(-k))*msq34_1_swap(j,k)
! quark-antiquark
        elseif((j>0).and.(k<0)) then 
          if(j==-k) then 
            if  ((abs(j) == 2) .or. (abs(j) == 4)) then
              msq(j,k)=
     &          ((Q(j)**2*virt_dip+D(j))*msqbits34_1(uub_uub)
     &         + (Q(2)**2*virt_dip)*msqbits34_1(uub_ccb)
     &         + (3._dp*Q(1)**2*virt_dip)*msqbits34_1(uub_ddb)
     &           +d_sum_q(j)*msqbits34_1(uub_ddb))
     &          +((Q(k)**2*virt_dip+D(-k))*msqbits34_1_swap(uub_uub)
     &         + (Q(2)**2*virt_dip)*msqbits34_1_swap(uub_ccb)
     &         + (3._dp*Q(1)**2*virt_dip)*msqbits34_1_swap(uub_ddb)
     &           +d_sum_q(-k)*msqbits34_1_swap(uub_ddb))
            else
              msq(j,k)=
     &            ((Q(j)**2*virt_dip+D(j))*msqbits34_1(ddb_ddb)
     &           + (2._dp*Q(2)**2*virt_dip)*msqbits34_1(ddb_uub)
     &           + (2._dp*Q(1)**2*virt_dip)*msqbits34_1(ddb_ssb)
     &           +d_sum_q(j)*msqbits34_1(ddb_ssb))
     &          +((Q(k)**2*virt_dip+D(-k))*msqbits34_1_swap(ddb_ddb)
     &           + (2._dp*Q(2)**2*virt_dip)*msqbits34_1_swap(ddb_uub)
     &           + (2._dp*Q(1)**2*virt_dip)*msqbits34_1_swap(ddb_ssb)
     &           +d_sum_q(-k)*msqbits34_1_swap(ddb_ssb))
            endif
          else
              msq(j,k)=(Q(j)**2*virt_dip+D(j))*msq34_1(j,k)
     &                +(Q(k)**2*virt_dip+D(-k))*msq34_1_swap(j,k)
  
          endif
! antiquark-quark
        elseif((j<0).and.(k>0)) then 
          if(j==-k) then 
            if  ((abs(j) == 2) .or. (abs(j) == 4)) then
              msq(j,k)=
     &          ((Q(j)**2*virt_dip+D(-j))*msqbits34_1(ubu_uub)
     &         + (Q(2)**2*virt_dip)*msqbits34_1(ubu_ccb)
     &         + (3._dp*Q(1)**2*virt_dip)*msqbits34_1(ubu_ddb)
     &           +d_sum_q(-j)*msqbits34_1(ubu_ddb))
     &          +((Q(k)**2*virt_dip+D(k))*msqbits34_1_swap(ubu_uub)
     &         + (Q(2)**2*virt_dip)*msqbits34_1_swap(ubu_ccb)
     &         + (3._dp*Q(1)**2*virt_dip)*msqbits34_1_swap(ubu_ddb)
     &           +d_sum_q(k)*msqbits34_1_swap(ubu_ddb))
            else
              msq(j,k)=
     &            ((Q(j)**2*virt_dip+D(-j))*msqbits34_1(dbd_ddb)
     &           + (2._dp*Q(2)**2*virt_dip)*msqbits34_1(dbd_uub)
     &           + (2._dp*Q(1)**2*virt_dip)*msqbits34_1(dbd_ssb)
     &           +d_sum_q(-j)*msqbits34_1(ddb_ssb))
     &         + ((Q(k)**2*virt_dip+D(k))*msqbits34_1_swap(dbd_ddb)
     &           + (2._dp*Q(2)**2*virt_dip)*msqbits34_1_swap(dbd_uub)
     &           + (2._dp*Q(1)**2*virt_dip)*msqbits34_1_swap(dbd_ssb)
     &           +d_sum_q(k)*msqbits34_1_swap(dbd_ssb))
            endif
          else
            msq(j,k)=(Q(j)**2*virt_dip+D(-j))*msq34_1(j,k)
     &              +(Q(k)**2*virt_dip+D(k))*msq34_1_swap(j,k)
  
          endif
        endif

!=========
!========= gluon -> photon fragmentation contributions (formally, higher order)
!=========
        if((j==0).and.(k<0)) then 
! gluon-antiquark
          msq(j,k)=msq(j,k)+D(0)*ga_ga
        elseif((j==0).and.(k>0)) then 
! gluon-quark
          msq(j,k)=msq(j,k)+D(0)*gq_gq
        elseif((j<0).and.(k==0)) then 
! antiquark-gluon
          msq(j,k)=msq(j,k)+D(0)*ag_ga
        elseif((j>0).and.(k==0)) then 
! quark-gluon
          msq(j,k)=msq(j,k)+D(0)*qg_gq
        elseif((j>0).and.(k<0).and.(j==-k)) then 
! quark-antiquark
          msq(j,k)=msq(j,k)+D(0)*qa_gg*two
        elseif((j<0).and.(k>0).and.(j==-k)) then 
! antiquark-quark
          msq(j,k)=msq(j,k)+D(0)*aq_gg*two
        elseif((j==0).and.(k==0)) then 
! gluon-gluon
          msq(j,k)=msq(j,k)+D(0)*gg_gg*two
        endif
        
      enddo
      enddo
      
      
      return 
      end
      

      function d_sum_q(j)
      implicit none
      include 'types.f'
      real(dp):: d_sum_q
!-----sum over i=1,5 i!=j   
       
      real(dp):: D(0:5) 
      common/D/D
      integer:: i,j
!$omp threadprivate(/D/)

      d_sum_q=0._dp
      do i=1,5 
         if(i.ne.j) then 
            d_sum_q=d_sum_q+D(i) 
         endif
      enddo 
      return 
      end
