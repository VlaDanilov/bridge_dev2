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
 
      subroutine higgsw(br)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'masses.f'
      real(dp):: wff,mfsq,br
c-----approximate form for the width of the standard model higgs
c-----valid for low masses
      wff(mfsq)=sqrt(2._dp)/8._dp/pi*gf*hmass*mfsq
     & *(1._dp-4._dp*mfsq/hmass**2)**1.5_dp

      hwidth=3._dp*(wff(mbsq)+wff(mcsq))+wff(mtausq)
      write(6,*) 'hmass,hwidth',hmass,hwidth
      write(6,*) 'mtausq,mcsq,mbsq',mtausq,mcsq,mbsq
      write(6,*) 
      br=3._dp*wff(mbsq)/hwidth
      return
      end
