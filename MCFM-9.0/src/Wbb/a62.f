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
 
      function a62(st,j1,j2,j3,j4,j5,j6,za,zb) 
      implicit none
      include 'types.f'
      complex(dp):: a62
      
************************************************************************
*     Author: R.K. Ellis                                               *
*     July, 1998.                                                      *
*     Updated 8/16/01: UV subtraction included                         *
************************************************************************
*     implementation of Eqs. (2.7) and (2.8) of BDKW hep-ph/9610370
*     with ns=0
*     character string st can take the value pp or pm
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6
      character*2 st
      complex(dp):: a6,aa6sf,aa6tp,aa6uv

c----Includes ultraviolet subtraction aa6uv (with +ve sign since
c--- it appears in interference terms,e.g. 1234 x 3214)
      call a6routine(st,j1,j2,j3,j4,j5,j6,za,zb,aa6sf,aa6tp,aa6uv) 
      a62=a6(st,j1,j2,j3,j4,j5,j6,za,zb)/xnsq
     & +(aa6sf*real(nf,dp)-aa6tp)/xn+aa6uv

      if (st == 'pp') then
      a62=a62+a6('pm',j1,j3,j2,j4,j5,j6,za,zb)*(one+one/xnsq)
     &       -a6('sl',j2,j3,j1,j4,j5,j6,za,zb)/xnsq
      elseif (st == 'pm') then
      a62=a62+a6('pp',j1,j3,j2,j4,j5,j6,za,zb)*(one+one/xnsq)
     &       +a6('sl',j3,j2,j1,j4,j5,j6,za,zb)/xnsq
      else
      write(6,*) 'Unimplemented st in a62',st 
      stop
      endif
      return
      end

