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
 
!--- 
      subroutine xI1gg(z,I1gg)
      implicit none
      include 'types.f'       
      include 'constants.f'       
       
      real(dp),intent(in) :: z
      real(dp) :: I1gg(0:2),omz
       
      real(dp) :: fpgg
       
      omz=one-z
      fpgg = 2*(omz+z**2)**2/z/omz
       
      I1gg(0) = -pisq/6.0_dp
      I1gg(1) = 2.0_dp*(omz+z**2)**2/z
      I1gg(2) = - fpgg*log(z)
      if(z.eq.1.0_dp) I1gg(2) = 2.0_dp
       
      I1gg(:) = CA*I1gg(:)
       
      return
      end
       
!--- 
      function I1gqi(z)
      implicit none
      include 'types.f'       
      include 'constants.f'       
       
      real(dp), intent(in) :: z
      real(dp) :: I1gqi
       
      real(dp) :: fpgq
       
      fpgq = (1.0_dp+(1.0_dp-z)**2)/z
       
      I1gqi = fpgq*log((1.0_dp-z)/z) + z
       
      I1gqi = I1gqi * CF
       
      return
      end
