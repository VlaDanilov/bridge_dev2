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
 
!===== T. Dennen, May 2014
!===== Initialise array of bubble integral values.
!===== For use with qqb->4gamma and qqb->2j2gamma
!===== coefficients by same author.
      subroutine Bubint_init(i1,i2,i3,i4,i5,i6,Bubint,ord)
        use mod_qcdloop_c
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'scale.f'
      integer:: i1,i2,i3,i4,i5,i6
      integer:: ord
      complex(dp):: Bubint(25)
      real(dp):: t

      t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

      Bubint(:) = czip

      Bubint(1) = qlI2(t(i1,i2,i3),0d0,0d0,musq,ord)
      Bubint(2) = qlI2(t(i1,i2,i4),0d0,0d0,musq,ord)
      Bubint(3) = qlI2(t(i1,i2,i5),0d0,0d0,musq,ord)
      Bubint(4) = qlI2(t(i1,i2,i6),0d0,0d0,musq,ord)
      Bubint(5) = qlI2(t(i1,i3,i4),0d0,0d0,musq,ord)
      Bubint(6) = qlI2(t(i1,i3,i5),0d0,0d0,musq,ord)
      Bubint(7) = qlI2(t(i1,i3,i6),0d0,0d0,musq,ord)
      Bubint(8) = qlI2(t(i1,i4,i5),0d0,0d0,musq,ord)
      Bubint(9) = qlI2(t(i1,i4,i6),0d0,0d0,musq,ord)
      Bubint(10) = qlI2(t(i1,i5,i6),0d0,0d0,musq,ord)
      Bubint(11) = qlI2(s(i1,i2),0d0,0d0,musq,ord)
      Bubint(12) = qlI2(s(i5,i6),0d0,0d0,musq,ord)
      Bubint(13) = qlI2(s(i4,i6),0d0,0d0,musq,ord)
      Bubint(14) = qlI2(s(i4,i5),0d0,0d0,musq,ord)
      Bubint(15) = qlI2(s(i3,i6),0d0,0d0,musq,ord)
      Bubint(16) = qlI2(s(i3,i5),0d0,0d0,musq,ord)
      Bubint(17) = qlI2(s(i3,i4),0d0,0d0,musq,ord)
      Bubint(18) = qlI2(s(i1,i3),0d0,0d0,musq,ord)
      Bubint(19) = qlI2(s(i2,i6),0d0,0d0,musq,ord)
      Bubint(20) = qlI2(s(i2,i5),0d0,0d0,musq,ord)
      Bubint(21) = qlI2(s(i2,i4),0d0,0d0,musq,ord)
      Bubint(22) = qlI2(s(i1,i4),0d0,0d0,musq,ord)
      Bubint(23) = qlI2(s(i2,i3),0d0,0d0,musq,ord)
      Bubint(24) = qlI2(s(i1,i5),0d0,0d0,musq,ord)
      Bubint(25) = qlI2(s(i1,i6),0d0,0d0,musq,ord)

      return
      end
