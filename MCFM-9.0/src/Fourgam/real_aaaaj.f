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
 
! T. Dennnen
! returns the squared amplitude for qQaaaag
      function real_aaaaj(j1,j2,j3,j4,j5,j6,j7,za,zb)
      implicit none
      include 'types.f'
      real(dp)::real_aaaaj
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer::j1,j2,j3,j4,j5,j6,j7
      integer::h5,h6
      complex(dp)::real_aaajj1(2,2,2,2,2,2),
     & real_aaajj2(2,2,2,2,2,2),real_aaaaj_hel(2,2,2,2,2,2)
      integer::h1,h2,h3,h4
      
      call real_aaajj_fill(j1,j2,j3,j4,j5,j6,j7,za,zb,real_aaajj1)
      call real_aaajj_fill(j1,j2,j3,j4,j5,j7,j6,za,zb,real_aaajj2)

      
      do h5=1,2
      do h6=1,2
      real_aaaaj_hel(:,:,:,:,h5,h6) = real_aaajj1(:,:,:,:,h5,h6)
     & +real_aaajj2(:,:,:,:,h6,h5) 
      enddo
      enddo
      
      real_aaaaj=zip

      do h1=1,2
         do h2=1,2
            do h3=1,2
               do h4=1,2
                  do h5=1,2
                     do h6=1,2
                        real_aaaaj=real_aaaaj
     &                    +abs(real_aaaaj_hel(h1,h2,h3,h4,h5,h6))**2
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      return
      end

