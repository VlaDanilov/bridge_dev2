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
 
      function virt5_VH(ip,za,zb)
      implicit none
      include 'types.f'
      real(dp):: virt5_VH
************************************************************************ 
*     Author: R.K. Ellis                                               *
*     July, 1999.                                                      *
*   Given za and zb calculate the                                      *
*   the interference of the amplitude for the process                  *
*   0--> qb_R(1)+q_L(2)+l_L(3)+a_R(4)+g_L/R(5)                         *
*   at one loop with the corresponding lowest order amplitude          *
*   summed over the polarizations of the emitted gluon                 *
*   Virtual terms are in units of 
*   (as/4/pi) (4 pi)^ep Gamma(1+ep)*Gamma(1-ep)^2/Gamma(1-2*ep)
************************************************************************
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'zprods_decl.f'
      integer:: ip(5)
      complex(dp):: A5LOm,A5NLOm,A5LOp,A5NLOp

c   0--> qb_R(1)+q_L(2)+l_L(3)+a_R(4)+g_L(5)
      call A5NLO_VH(ip(1),ip(2),ip(3),ip(4),ip(5),za,zb,A5LOm,A5NLOm)
c   0--> qb_R(1)+q_L(2)+l_L(3)+a_R(4)+g_R(5)
      call A5NLO_VH(ip(2),ip(1),ip(4),ip(3),ip(5),zb,za,A5LOp,A5NLOp)

      virt5_VH=
     & +ason2pi*(real(conjg(A5LOp)*A5NLOp,dp)
     &          +real(conjg(A5LOm)*A5NLOm,dp))

      return
      end

