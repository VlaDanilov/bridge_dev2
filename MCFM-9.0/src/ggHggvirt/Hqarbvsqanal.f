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
 
      function Hqarbvsqanal(i1,i2,i3,i4)
      implicit none
      include 'types.f'
      real(dp):: Hqarbvsqanal
c--- This routine is simply a wrapper to the (non-identical) four
c--- quark virtual routines. By changing "imode" below, one can
c--- switch between a squared ME calculation and one using amplitudes
      
      integer:: i1,i2,i3,i4,imode
      real(dp):: Hqarbvsq,Ampvirtsq_AQaq_nonid,res_sq,res_amp
      
      imode=2
c--- imode=1   ! compute square using amplitudes from H4pCode      
c--- imode=2   ! compute square directly using results from EGZ
c--- imode=3   ! compare the two calculations 
      
      if ((imode == 1) .or. (imode == 3)) then
        res_amp=Ampvirtsq_AQaq_nonid(i1,i2,i3,i4)
      endif
      
      if ((imode == 2) .or. (imode == 3)) then
        res_sq=Hqarbvsq(i1,i2,i3,i4)
      endif

c--- Checked that Hqarbvsq == Ampvirtsq_AQaq_nonid on 23/10/09
      if (imode == 3) then
       write(6,99) 'i1,i2,i3,i4 Hqarbvsq',i1,i2,i3,i4,1._dp-res_sq/res_amp
      endif
      
      if (imode == 2) then
        Hqarbvsqanal=res_sq
      else
        Hqarbvsqanal=res_amp
      endif
      
      return
      
   99 format(a20,4i3,e21.12) 
      
      end
      
