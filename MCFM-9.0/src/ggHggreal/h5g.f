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
 
      SUBROUTINE h5g(ANS)
C  
C RETURNS AMPLITUDE SQUARED SUMMED OVER COLORS
C AND HELICITIES FOR THE POINT IN PHASE SPACE P(mxpart,4)
C FOR PROCESS : g g -> g g g h 

C overall coupling factors have been removed
C  
      
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
C  
C CONSTANTS
C  
      integer,parameter:: NEXTERNAL=5, NCOMB= 16

C  
C ARGUMENTS 
C  
      real(dp)::ANS
C  
C LOCAL VARIABLES 
C  
      integer::NHEL(NEXTERNAL,NCOMB)
      real(dp)::T,GG_GGG
      integer:: IHEL
      include 'hels.f'

C ----------
C BEGIN CODE
C ----------


      ANS=zip

c--  sum over helicities
      DO IHEL=1,16
              T=GG_GGG(NHEL(1,IHEL))            
              ANS=ANS+T
      ENDDO

c--  Multiply by two to account for the other 16 helicity configs
      ANS=ANS*two

      RETURN
      END
       
       

