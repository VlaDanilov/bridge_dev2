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
 
      subroutine extract(target_struc,target_struc_v)
      implicit none
      include 'types.f'
************************************************************************
*     Author: J.M. Campbell                                            *
*     February, 2005.                                                  *
*                                                                      *
*     Extracts the matrix element colour structure from the common     *
*     include file 'msq_struc' and returns it in the arrays            *
*     target_struc (lowest order) and target_struc_v (lowest order     *
*     contracted with a vector)                                        *
*                                                                      *
*     For efficiency, this is only done for indices (-1,0,+1)          *
*                                                                      *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'msq_struc.f'
      integer::istruc,j,k      
      real(dp)::target_struc(8,-1:1,-1:1),
     &               target_struc_v(8,-1:1,-1:1)
      
      do j=-1,1
      do k=-1,1
      do istruc=1,8
        target_struc(istruc,j,k)=msq_struc(istruc,j,k)
c--- dummy for now: needs to be filled once the appropriate
c--- modifications to the gvec routine have been made       
        target_struc_v(istruc,j,k)=msq_strucv(istruc,j,k)
      enddo
      enddo
      enddo

      return
      end
      
