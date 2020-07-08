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
 
      subroutine VKlSt(q,m,aux,i,f)
      implicit none
      include 'types.f'
      include 'constants.f'
c     subroutine for massive V spinor alla Kleiss and Stirling
C     q complex vector q(4),q(1),q(2),q(3)=E_q,qx,qy,qz
C     'aux' massless auxiliary complex vector 
C     i polarization 1 or -1
C     f output spinor
      
      integer:: i,j
      real(dp):: m
      complex(dp):: q(4),qt(4),aux(4),f(4),fqt(4),eta(4),den,
     & qDaux,qsq
      
      if (m == 0._dp) then
        call uspinor0(q,-i,f)
        return
      else 
        qDaux=q(4)*aux(4)-q(1)*aux(1)-q(2)*aux(2)-q(3)*aux(3)
        qsq=q(4)**2-q(1)**2-q(2)**2-q(3)**2
        qt(:)=q(:)-0.5_dp*qsq/qDaux*aux(:)
        call uspinor0(aux,i,eta)
        call ubarspinor0(qt,-i,fqt)
        den=eta(1)*fqt(1)+eta(2)*fqt(2)+eta(3)*fqt(3)+eta(4)*fqt(4)
      endif

      call kslashU(q,eta,f)
      do j=1,4
      f(j)=(f(j)-eta(j)*m)/den
      enddo
      return
      end


