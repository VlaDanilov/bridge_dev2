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
 
!--------------------------------------------------------------- 
!   This subroutine checks the number of external dipoles----
!---absorbing the correct number into the fragmenation functions 
!   it then returns the finite (msq_qcd*dip) ---------------------
!--------------------------------------------------------------- 
!--- Author C. Williams Nov 2012
!-----------------------------------------------------------------
! modified under frag re-write March 2013

      subroutine qqb_dm_monophot_fragdips(p,p_phys,qcd_tree,msq_out) 
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'frag.f'
      real(dp):: p(mxpart,4),p_phys(mxpart,4)
      real(dp):: msq_qcd(-nf:nf,-nf:nf),msq_out(-nf:nf,-nf:nf)
      integer:: j,k
      real(dp):: virt_dips,xl,dot,fsq 
      real(dp):: aewo2pi,fi_gaq
      external qcd_tree
   
      aewo2pi=esq/(fourpi*twopi)      
      
      fsq=frag_scale**2

      xl=log(-two*dot(p_phys,2,5)/fsq)
      
      virt_dips=+aewo2pi*(fi_gaq(z_frag,p_phys,xl,5,2,2))
      
      do j=-nf,nf
         do k=-nf,nf
            msq_qcd(j,k)=0d0
            msq_out(j,k)=0d0
         enddo
      enddo
      
      call qcd_tree(p,msq_qcd) 

      do j=-nf,nf
         do k=-nf,nf
            
           if((j==0).and.(k.ne.0)) then 
              msq_out(j,k)=Q(k)**2*msq_qcd(j,k)*virt_dips
           elseif((j.ne.0).and.(k==0)) then 
              msq_out(j,k)=Q(j)**2*msq_qcd(j,k)*virt_dips
           endif
            
         enddo
      enddo
     
     
      return 
      end subroutine


        
