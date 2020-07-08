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
 
      subroutine htautaudecay(p,jm,jp,msq)
      implicit none
      include 'types.f'
************************************************************************
*     Author: J.M. Campbell, August 2012                               *
*                                                                      *
*     matrix element squared for the process of                        *
*     Higgs decay  H --> tau^-(jm)+tau^+(jp)                           *
*     with tau mass included                                           *
************************************************************************
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      integer:: jm,jp
      real(dp):: p(mxpart,4),s56,msq,msqhtautau
      
      s56=two*(p(jm,4)*p(jp,4)-p(jm,1)*p(jp,1)
     &        -p(jm,2)*p(jp,2)-p(jm,3)*p(jp,3))+two*mtau**2
      
      msq=msqhtautau(s56)
      
      return
      end
      
      
      
      function msqhtautau(s)
      implicit none
      include 'types.f'
      real(dp):: msqhtautau
      
      real(dp):: s
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ewcouple.f'
 
      msqhtautau=gwsq*mtausq/(four*wmass**2)*two*(s-four*mtau**2)

      return
      end
      
      
