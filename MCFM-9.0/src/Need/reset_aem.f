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
 

      subroutine reset_aem(ae_in)
      implicit none
      include 'types.f'
      include 'ewcouple.f'
      include 'first.f'
      include 'constants.f'
      include 'mpicommon.f'
      real(kind=dp), intent(in) :: ae_in
      real(kind=dp) :: esq_old,ae_old

      esq_old=esq
      ae_old=esq/(fourpi)

      esq=ae_in*fourpi
      
      if(first) then
         first=.false.
         if (rank == 0) then
         write(6,*) '*********** Changed alpha_EM **********'
         write(6,88) '*  old value of esq :',esq_old,'     *'
         write(6,88) '*  new value of esq :',esq,'     *'
         write(6,88) '*  new 1/alpha_EM   :',fourpi/esq,'     *'
         write(6,*) '***************************************'
         endif
      endif

 88   format(1x,a21,f12.8,a6)
      return
      end
        
