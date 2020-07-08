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
 
      subroutine breitw(x1,mminsq,mmaxsq,rmass,rwidth,msq,wt)       
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
c---- Given a number 0<x<1 generate a mass-squared msq and a weight wt 
c---- such that mminsq<msq<mmaxsq
c---- points are generated around resonance position rmass, but 
c---- breit-wigner should still be included in the matrix element
c     wt is the jacobian between integration in msq and integration in x1
      real(dp):: x1,mminsq,mmaxsq,rmass,rwidth,msq,wt
      real(dp):: almin,almax,al,tanal
      include 'zerowidth.f'

c--- in case the maximum msq is very small, just generate linearly for safety
      if ((mmaxsq < rmass*1.e-3_dp) .and. (zerowidth .eqv. .false.)) then
        msq=mminsq+x1*(mmaxsq-mminsq)
        wt=mmaxsq-mminsq
        return
      endif

      if (zerowidth) then
          tanal=0._dp
          almax=+pi/two
          almin=-pi/two
      else
c--- if generating off-resonance then just do logarithmic sampling in msq
          if (mminsq > rmass**2) then
            wt=one
            call pick(2,msq,mminsq,mmaxsq,x1,wt)
            return
          endif          
          almin=atan((mminsq-rmass**2)/rmass/rwidth)
          almax=atan((mmaxsq-rmass**2)/rmass/rwidth)
          al=(almax-almin)*x1+almin
          tanal=tan(al)
      endif

      msq=rmass**2+rmass*rwidth*tanal
c---- bw=(1._dp+tanal**2)*rmass**2*rwidth**2
      wt=(almax-almin)*rmass*rwidth*(1._dp+tanal**2)

      if (msq < 0._dp) then
        msq=mminsq
        wt=0._dp
      endif
      
      return
      end

