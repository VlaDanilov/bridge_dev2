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
 
      subroutine ggZZparsecheck(helstring,h1,h2,resarr,string)
      implicit none
      include 'types.f'
      include 'cplx.h'
      
      logical:: found
      integer:: h1,h2,h3,h4,j3,j4
      character*(*) string
      character*160 line
      character*2 helstring
      character*21 filename
      real(dp):: xr,xc,xabs,ratio
      complex(dp):: resarr(2,2)

      do h3=-1,1,2
      do h4=-1,1,2     
      
      found=.false. 

      write(filename,91) h1,h2,h3,h4,helstring
   91 format(SP,'out.KC_',i2,'_',i2,'_',i2,'_',i2,'_',a2)
   
      open(unit=21,file=filename,status='unknown')
   
   77 continue   
      read(21,'(a160)',end=99) line
      if (index(line,trim(string)) == 0) goto 77
      close(unit=21)
      found=.true.
      
      if (index(line,'xe') == 0) then
        read(line,71) xr,xc,xabs 
   71   format(60x,3f24.15)   
      else
        read(line,72) xr,xc,xabs 
   72   format(18x,3f24.15)   
      endif
      
      j3=1+(h3+1)/2
      j4=1+(h4+1)/2
      ratio=real(resarr(j3,j4)/cplx2(xr,xc))
   99 continue  
      if (found .eqv. .false.) then
        if (abs(resarr(j3,j4)) < 1d-12) then
          ratio=1._dp
        else
          ratio=-77._dp ! signals zero in KC but non-zero here
        endif 
      endif
      write(6,*) string,': ',helstring,h1,h2,h3,h4,ratio
      
      enddo
      enddo
      
      return
      end
               
