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
 
      subroutine qqWZggampf(i1,i2,i3,i4,i5,i6,i7,i8,b7,b8,za,zb,msq)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zprods_decl.f'
      integer:: i1,i2,i3,i4,i5,i6,i7,i8,b7,b8,h5,h7,h8
      real(dp):: msq
      complex(dp):: a78(2,2,2),a87(2,2,2),aq(2,2,2)

      call qqWZggamp(i1,i2,i3,i4,i5,i6,i7,i8,b7,b8,za,zb,a78,a87)
      aq(:,:,:)=a78(:,:,:)+a87(:,:,:)

      msq=zip
      do h5=1,2
      do h7=1,2
      do h8=1,2
      msq=msq+gsq**2*esq**4*V*xn*(
     & +real(a78(h5,h7,h8)*conjg(a78(h5,h7,h8)))
     & +real(a87(h5,h7,h8)*conjg(a87(h5,h7,h8)))
     & -real(aq(h5,h7,h8)*conjg(aq(h5,h7,h8)))
     &  /xn**2)
      enddo
      enddo
      enddo
      return
      end
