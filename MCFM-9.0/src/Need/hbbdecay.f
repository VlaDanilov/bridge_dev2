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
 
      subroutine hbbdecay(p,ib,ibb,msq)
      implicit none
      include 'types.f'
************************************************************************
*     Author: J.M. Campbell, June 2012                                 *
*                                                                      *
*     matrix element squared for the process of                        *
*     Higgs decay  H --> b(ib)+b~(ibb)                                 *
*     with bottom mass included                                        *
************************************************************************
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      integer:: ib,ibb
      real(dp):: p(mxpart,4),s56,msq,msqhbb
      
      s56=2._dp*(p(ib,4)*p(ibb,4)-p(ib,1)*p(ibb,1)
     &        -p(ib,2)*p(ibb,2)-p(ib,3)*p(ibb,3))+2._dp*mb**2
      
      msq=msqhbb(s56)
      
      return
      end
      
      
      
      function msqhbb(s)
      implicit none
      include 'types.f'
      real(dp):: msqhbb
      real(dp):: s
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'couple.f'
      include 'kpart.f'
      include 'msbarmasses.f'
      include 'ewcouple.f'
      include 'first.f'
      include 'hbbparams.f'

      msqhbb=xn*gwsq*mb_eff**2/(4._dp*wmass**2)*2._dp*(s-4._dp*mb**2)

      return
      end
      
      
