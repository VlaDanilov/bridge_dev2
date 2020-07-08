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
 
      subroutine qqb_dirgam_gvec(p,n,in,msq)
      implicit none
      include 'types.f'
C*********************************************************************** 
c     Author: R.K. Ellis                                               *
c     October, 2002.                                                   *
c     Matrix element for gamma production                              *
c     averaged over initial colours and spins                          *
c     contracted with the vector n(mu) (orthogonal to p4)              *
c     q(-p1)+qbar(-p2)--> gamma(p3)+ g(p4)                             *
C*********************************************************************** 
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: j,k,in
C--in is the label of the parton dotted with n
      real(dp):: msq(-nf:nf,-nf:nf),msqa(-nf:nf,-nf:nf),
     & p(mxpart,4),n(4),nDn

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo
 
      call qqb_dirgam(p,msqa)
      nDn=n(4)**2-n(1)**2-n(2)**2-n(3)**2

      call checkndotp(p,n,in)

      do j=-nf,nf

      if (in == 1) then
        msq(0,j)=-0.5_dp*nDn*msqa(0,j)
      elseif (in == 2) then
        msq(j,0)=-0.5_dp*nDn*msqa(j,0)
      elseif (in == 4) then      
        msq(j,-j)=-0.5_dp*nDn*msqa(j,-j)
      endif

      enddo

      return
      end
