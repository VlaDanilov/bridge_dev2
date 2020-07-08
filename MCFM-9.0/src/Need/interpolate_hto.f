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
 
      subroutine interpolate_hto(mh,hwidth)
      implicit none
      include 'types.f'
c--- This routine reads in the file hto_output.dat, prepared
c--- by the HTO program of G. Passarino, and uses linear interpolation
c--- to obtain the Higgs width parameter CPGH for any input Higgs
c--- mass (mh) between 50 GeV and 1500 GeV
      
      logical:: first
      real(dp):: mh,hwidth,m1,g1,m2,g2,eps
      data first/.true./
      save first
      parameter(eps=1.e-8_dp)
      
      if ((mh < 50._dp) .or. (mh > 1500._dp)) then
        write(6,*) 'Higgs mass outside HTO interpolation range;'
        write(6,*) 'requires  50 GeV < mH < 1500 GeV'  
        stop
      endif

      if (first) then
        first=.false.
        write(6,*)
        write(6,*) 'Higgs width parameter obtained by interpolating:'
        write(6,*) '   '
        write(6,*) '*'
        write(6,*) ' *'
        write(6,*) '  * --- cpHTO v 1.1 (May 2012) by Giampiero '
        write(6,*) ' *'
        write(6,*) '*'      
      endif
      
      open(unit=12,file='hto_output.dat',status='unknown')
      read(12,*)
      read(12,*)
   77 continue
      read(12,99) m1,g1
      if (mh-m1 > 0.5_dp+eps) goto 77
      read(12,99) m2,g2
      close(12)
      
      hwidth=(g2-g1)*(mh-m1)/(m2-m1)+g1

      return

   99 format(f10.1,f21.12)
      end
      
