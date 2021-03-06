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
 
      subroutine qqb_fourgam_v(p,msq)
      implicit none
      include 'types.f'
************************************************************************
*    Author: C Williams  and T. Dennen                                     *
*    Sept, 2014.                                                      *
*    LO matrix element squared, averaged over initial colors         *
*    and spins (plus all crossings)                                    *
c     q(-p1)+qbar(-p2) --> gam(p3) + gam(p4) + gam(p5) + gam(p6)         *
************************************************************************
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcharge.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zprods_decl.f'
      include 'epinv.f' 
      include 'epinv2.f'
      include 'scheme.f' 
      integer:: j
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf),
     & qqb,qbq,fac,cfac
      real(dp):: symfac 
      parameter (symfac=one/24d0)
      complex(dp):: qqb_aaaa_nlo_e2(2,2,2,2,2)
     &     ,qqb_aaaa_lo(2,2,2,2,2)
      complex(dp):: qqb_aaaa_nlo_e1(2,2,2,2,2)
      complex(dp):: qqb_aaaa_nlo_e0(2,2,2,2,2)
      integer:: i1,i2,i3,i4,i5
      complex(dp):: qqb_aaaa_nlo(2,2,2,2,2)

      scheme='dred'

      call spinoru(6,p,za,zb)

      qqb_aaaa_nlo_e2(:,:,:,:,:)=czip
      qqb_aaaa_nlo_e1(:,:,:,:,:)=czip
      qqb_aaaa_nlo_e0(:,:,:,:,:)=czip

      qqb_aaaa_lo(:,:,:,:,:)=czip
 
      if(epinv.ne.0d0) then 
!--------epinv^2 pieces
      call aaaa_fill(1,2,3,4,5,6,za,zb,qqb_aaaa_nlo_e2,qqb_aaaa_lo,-2)
!-------- epinv^1 pieces
      call aaaa_fill(1,2,3,4,5,6,za,zb,qqb_aaaa_nlo_e1,qqb_aaaa_lo,-1)
       endif
!-------- epinv^0 pieces
      call aaaa_fill(1,2,3,4,5,6,za,zb,qqb_aaaa_nlo_e0,qqb_aaaa_lo,0)

      qqb_aaaa_nlo(:,:,:,:,:)=qqb_aaaa_nlo_e2(:,:,:,:,:)*epinv*epinv2
     &     +qqb_aaaa_nlo_e1(:,:,:,:,:)*epinv+qqb_aaaa_nlo_e0(:,:,:,:,:)

      qqb=0d0
      
c--- overall coupling, LO fac
      fac=4d0*esq**4*xn*symfac*aveqq

C------ NLO fac
      fac=-2d0*fac*Cf*ason2pi
      
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i4=1,2
      do i5=1,2
         qqb=qqb+Dble(im*qqb_aaaa_nlo(i1,i2,i3,i4,i5)*
     &        conjg(qqb_aaaa_lo(i1,i2,i3,i4,i5)))
      enddo
      enddo
      enddo
      enddo
      enddo
 
      qbq=qqb

      msq(:,:)=0d0
      do j=1,nf
        cfac=Q(j)**8
        msq(j,-j)=cfac*qqb*fac
        msq(-j,j)=cfac*qbq*fac
      enddo
      
      return
      end
      
      
