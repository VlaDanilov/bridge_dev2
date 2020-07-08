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
 
      subroutine dittdrein(p,l1,l2,costhdd)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: j,l1,l2
      real(dp):: costhdd,p(mxpart,4),pt,ptem,ptep,xl
      real(dp):: pem(4),pep(4),psum(4)
      real(dp):: p_cm(4)

      do j=1,4
      pem(j)=p(l1,j)
      pep(j)=p(l2,j)
      psum(j)=pem(j)+pep(j)
      enddo

      ptem=pt(l1,p)
      ptep=pt(l2,p)


      if (ptem > ptep) then
      call boosta(psum,pem,p_cm)
      else
      call boosta(psum,pep,p_cm)
      endif

      xl=sqrt((psum(1)**2+psum(2)**2+psum(3)**2)
     &       *(p_cm(1)**2+p_cm(2)**2+p_cm(3)**2))
      costhdd=(p_cm(1)*psum(1)+p_cm(2)*psum(2)+p_cm(3)*psum(3))/xl

      return
      end


 


