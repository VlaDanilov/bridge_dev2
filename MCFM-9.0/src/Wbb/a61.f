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
 
      function a61(st,j1,j2,j3,j4,j5,j6,za,zb) 
      implicit none
      include 'types.f'
      complex(dp):: a61
************************************************************************
*     Author: R.K. Ellis                                               *
*     July, 1998.                                                      *
************************************************************************
*     implementation of Eqs. (2.7) and (2.8) of BDKW hep-ph/9610370
*     with ns=0
*     character string st can take the value pp or pm
*     a61(pp,j1,j2,j3,j4,j5,j6,za,zb) corresponds to 
*     q(j1,+)+Q(j3,-)+e(j6,+)+q~(j4)+Q~(j2)+e~(j5)
*     a61(pm,j1,j2,j3,j4,j5,j6,za,zb) corresponds to 
*     q(j1,+)+Q(j3,+)+e(j6,+)+q~(j4)+Q~(j2)+e~(j5)
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6
      character*2 st
      complex(dp):: a6,aa6sf,aa6tp,aa6uv

c----Includes ultraviolet subtraction aa6uv
      call a6routine(st,j1,j2,j3,j4,j5,j6,za,zb,aa6sf,aa6tp,aa6uv) 
      a61=(one-two/xnsq)*a6(st,j1,j2,j3,j4,j5,j6,za,zb)
     & -(real(nf,dp)*aa6sf-aa6tp)/xn-aa6uv

      if (st == 'pp') then
      a61=a61+(-two*a6('pm',j1,j3,j2,j4,j5,j6,za,zb)
     &             +a6('sl',j2,j3,j1,j4,j5,j6,za,zb))/xnsq
      elseif (st == 'pm') then
      a61=a61+(-two*a6('pp',j1,j3,j2,j4,j5,j6,za,zb)
     &             -a6('sl',j3,j2,j1,j4,j5,j6,za,zb))/xnsq
      else
      write(6,*) 'Unimplemented st in a61',st 
      stop
      endif
      return
      end

