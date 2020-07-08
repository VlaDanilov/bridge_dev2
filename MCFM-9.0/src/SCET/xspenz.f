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
 
      FUNCTION XSPENZ(Z)
C
C CALCULATES THE EULER DILOG FOR ARBITRARY COMPLEX ARGUMENT Z
C CALLS THE SERIES CALCULATION XCDIL 
C T. RIEMANN, ABOUT 1984, FOLLOWING:  
C t Hooft, G. and Veltman, M., SCALAR ONE LOOP INTEGRALS, 
C Nucl. Phys.B153, 365, App. A
C 
      implicit none
      include 'types.f'
      include 'constants.f'
      complex(dp)::xspenz,xcdil,z,cspenz
      real(dp)::REZ,AIZ,AAZ,AAS,AIS,RE1,RES
      integer::n,jprint
c      logical::first
c      IMPLICIT COMPLEX*16(X,Y)
c      IMPLICIT REAL*8(A-H,O-W,Z)
c      COMPLEX*16 Z,XCDIL,CSPENZ,LOG
c      COMMON/CDZPIF/PI,F1
      EXTERNAL XCDIL
      DATA N/0/

C
      JPRINT=0
C

      REZ=REAL(Z)
      AIZ=AIMAG(Z)
      AAZ=ABS(Z)
      IF(AAZ) 11,9,11
 9    CSPENZ=DCMPLX(0.D0,0.D0)
      GOTO 99
 11   IF(AAZ-one) 6,4,1
 1    RE1=REAL(one/Z)
      IF(RE1-half) 3,3,2
2     CONTINUE
      CSPENZ=XCDIL(one-one/Z)-2.*pisqo6-LOG(Z)*LOG(one-one/Z)
     U      -half*(LOG(-Z))**2
      GOTO 99
3     CONTINUE
      CSPENZ=-XCDIL(one/Z)-pisqo6-half*LOG(-Z)**2
      GOTO 99
 4    IF(REZ-one) 6,5,1
 5    CSPENZ=DCMPLX(pisqo6,0.D0)
      GOTO 99
 6    IF(REZ-half) 7,7,8
7     CONTINUE
      CSPENZ=XCDIL(Z)
      GOTO 99
8     CONTINUE
      CSPENZ=-XCDIL(one-Z)+pisqo6-LOG(Z)*LOG(one-Z)
 99   CONTINUE
      AAS= ABS(CSPENZ)
      RES=REAL(CSPENZ)
      AIS=AIMAG(CSPENZ)
      IF(JPRINT) 97,97,98
 98   CONTINUE
 97   CONTINUE
      XSPENZ=CSPENZ
      END

C================================================================ 
