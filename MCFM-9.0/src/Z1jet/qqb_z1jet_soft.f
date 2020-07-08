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
 
      subroutine qqb_z1jet_soft(P,msq)
      implicit none
      include 'types.f'
c---Soft matrix element squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) --> Z +g(p5) +g(p6)
c                          |
c                          --> l(p3)+a(p4)
c                            
c   with p6 soft
c--all momenta incoming
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      integer:: j,k
      real(dp):: P(mxpart,4),msq(-nf:nf,-nf:nf),
     & msq0(-nf:nf,-nf:nf),s(mxpart,mxpart),
     & facqqb,facqbq,facgqb,facqbg,facgq,facqg
      real(dp):: eik12,eik15,eik25
      
      call dotem(6,p,s)

      eik12=s(1,2)/(s(1,6)*s(2,6))
      eik15=s(1,5)/(s(1,6)*s(6,5))
      eik25=s(2,5)/(s(2,6)*s(6,5))

      facqqb=two*gsq*(xn*(eik15+eik25)-eik12/xn)
      facqbq=facqqb
      facqg=two*gsq*(xn*(eik12+eik25)-eik15/xn)
      facgq=two*gsq*(xn*(eik12+eik15)-eik25/xn)
      facgqb=two*gsq*(xn*(eik15+eik12)-eik25/xn)
      facqbg=two*gsq*(xn*(eik25+eik12)-eik15/xn)

c--- SUB-LEADING ONLY
c      facqqb=two*gsq*(-eik12/xn)
c      facqbq=facqbq
c      facqg=two*gsq*(-eik15/xn)
c      facgq=two*gsq*(-eik25/xn)
c      facgqb=two*gsq*(-eik25/xn)
c      facqbg=two*gsq*(-eik15/xn)

      call qqb_z1jet(p,msq0)
      
      do j=-nf,nf
      do k=-nf,nf
            if ((j > 0) .and. (k<0)) then
               msq(j,k)=facqqb*msq0(j,k)
            elseif ((j < 0) .and. (k>0)) then
               msq(j,k)=facqbq*msq0(j,k)
            elseif ((j > 0) .and. (k==0)) then
               msq(j,k)=facqg*msq0(j,k)
            elseif ((j < 0) .and. (k==0)) then
               msq(j,k)=facqbg*msq0(j,k)
            elseif ((j == 0) .and. (k>0)) then
               msq(j,k)=facgq*msq0(j,k)
            elseif ((j == 0) .and. (k<0)) then
               msq(j,k)=facgqb*msq0(j,k)
            endif
      enddo
      enddo
      return
      end

 
