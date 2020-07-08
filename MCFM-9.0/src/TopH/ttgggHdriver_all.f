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
 
      subroutine ttgggHdriver(q,ampsq)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: q(mxpart,4),ampsq
      integer:: h1,h2,h3,h4,h6
      complex(dp):: m(6),mqed,
     & m126(2,2,2,2,2),m261(2,2,2,2,2),m612(2,2,2,2,2),
     & m621(2,2,2,2,2),m216(2,2,2,2,2),m162(2,2,2,2,2),
     & x126(2,2,2,2,2),x261(2,2,2,2,2),x612(2,2,2,2,2),
     & x621(2,2,2,2,2),x216(2,2,2,2,2),x162(2,2,2,2,2),
     & dum(2,2,2,2,2)
      call ttgggHamp(q,3,4,1,2,6,1,1,x126,x261,x612,x621,x216,x162)
      call ttgggHamp(q,3,4,1,2,6,1,1,m126,dum,dum,dum,dum,dum)
      call ttgggHamp(q,3,4,2,6,1,1,1,m261,dum,dum,dum,dum,dum)
      call ttgggHamp(q,3,4,6,1,2,1,1,m612,dum,dum,dum,dum,dum)
      call ttgggHamp(q,3,4,6,2,1,1,1,m621,dum,dum,dum,dum,dum)
      call ttgggHamp(q,3,4,2,1,6,1,1,m216,dum,dum,dum,dum,dum)
      call ttgggHamp(q,3,4,1,6,2,1,1,m162,dum,dum,dum,dum,dum)
      ampsq=0d0
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
      do h6=1,2

      m(1)=m126(h1,h2,h3,h4,h6)
      m(2)=m261(h1,h2,h4,h6,h3)
      m(3)=m612(h1,h2,h6,h3,h4)
      m(4)=m621(h1,h2,h6,h4,h3)
      m(5)=m216(h1,h2,h4,h3,h6)
      m(6)=m162(h1,h2,h3,h6,h4)
      write(6,*) 'm(1)',m(1)/x126(h1,h2,h3,h4,h6)
      write(6,*) 'm(2)',m(2)/x261(h1,h2,h3,h4,h6)
      write(6,*) 'm(3)',m(3)/x612(h1,h2,h3,h4,h6)
      write(6,*) 'm(4)',m(4)/x621(h1,h2,h3,h4,h6)
      write(6,*) 'm(5)',m(5)/x216(h1,h2,h3,h4,h6)
      write(6,*) 'm(6)',m(6)/x162(h1,h2,h3,h4,h6)
      pause
      mqed=m(1)+m(2)+m(3)+m(4)+m(5)+m(6)

c      ampsq=ampsq+(xnsq-2d0)*(
c     & +real(m(1)*conjg(m(1)))+real(m(2)*conjg(m(2)))
c     & +real(m(3)*conjg(m(3)))+real(m(4)*conjg(m(4)))
c     & +real(m(5)*conjg(m(5)))+real(m(6)*conjg(m(6))))
c       ampsq=ampsq
c     & +real(m(1)*conjg(m(4)-m(6)-m(5)))
c     & +real(m(2)*conjg(m(6)-m(5)-m(4)))
c     & +real(m(3)*conjg(m(5)-m(4)-m(6)))
c     & +real(m(4)*conjg(m(1)-m(2)-m(3)))
c     & +real(m(5)*conjg(m(3)-m(1)-m(2)))
c     & +real(m(6)*conjg(m(2)-m(3)-m(1)))
c      ampsq=ampsq+real(mqed*conjg(mqed))/xn**2

      ampsq=ampsq+xnsq*(
     & +abs(m(1))**2+abs(m(2))**2+abs(m(3))**2
     & +abs(m(4))**2+abs(m(5))**2+abs(m(6))**2)
     
      ampsq=ampsq-(
     & +abs(m(1)+m(6)+m(3))**2
     & +abs(m(5)+m(2)+m(4))**2
     & +abs(m(1)+m(5)+m(6))**2
     & +abs(m(4)+m(2)+m(3))**2
     & +abs(m(1)+m(5)+m(2))**2
     & +abs(m(6)+m(3)+m(4))**2)

      ampsq=ampsq+(1d0+1d0/xnsq)*abs(mqed)**2

      enddo
      enddo
      enddo
      enddo
      enddo

      return
      end
