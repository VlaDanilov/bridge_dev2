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
 
*********************************************************
*AUTHOR: FABIO MALTONI                                  *
*DATE  : 10/10/2001                                     *
*NOTES : PROGRAM GENERATED BY BHIGGS.M                  *
*********************************************************

      function  BBBBH(I1,I2,I3,I4)                                    
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_com.f'

      real(dp)::  BBBBH
* ---------------------------------------------------------------------
* returns the matrix element squared for the process                   
*                                                                      
*          0 -> bbar1 b2 bbar3 b4 h                                    
*                                                                      
* All momenta outgoing.                                                
* No averaging is performed for initial spins or colors.               
* No symmetry factors are included.                                    
* No STRONG coupling is included.                                      
* No Bottom-Higgs coupling is included.                                
* ---------------------------------------------------------------------

*     ARGUMENTS 
      integer:: I1,I2,I3,I4

*     LOCAL 
      real(dp)::  T123,T124,T134,T234
      complex(dp):: AMP(2,4)
      integer:: L


      t123=s(i1,i2)+s(i2,i3)+s(i1,i3)
      t124=s(i1,i2)+s(i2,i4)+s(i1,i4)
      t134=s(i1,i3)+s(i1,i4)+s(i3,i4)
      t234=s(i2,i3)+s(i2,i4)+s(i3,i4)

      amp(1,1)=
     &   (-two*(za(i4,i2)*zb(i2,i1) + za(i4,i3)*zb(i3,i1))*zb(i3,i2))/
     -   (T234*s(i3,i4)) + (two*zb(i3,i1)*
     -     (za(i4,i1)*zb(i2,i1) - za(i4,i3)*zb(i3,i2)))/
     -   (T134*s(i3,i4))

      amp(2,1)=
     &   (-two*zb(i1,i2)*(za(i4,i1)*zb(i1,i3) + za(i4,i2)*zb(i2,i3)))/
     -   (T124*s(i1,i4)) + (two*zb(i1,i3)*
     -     (-(za(i4,i1)*zb(i1,i2)) + za(i4,i3)*zb(i2,i3)))/
     -   (T134*s(i1,i4))

      amp(1,2)=
     &   (-two*(za(i3,i2)*zb(i2,i1) + za(i3,i4)*zb(i4,i1))*zb(i4,i2))/
     -   (T234*s(i3,i4)) + (two*zb(i1,i4)*
     -     (za(i1,i3)*zb(i2,i1) - za(i4,i3)*zb(i4,i2)))/
     -   (T134*s(i3,i4))

      amp(2,2)=
     &        (two*zb(i1,i2)*(-(za(i2,i3)*zb(i2,i4)) + 
     -       za(i1,i3)*zb(i4,i1)))/(T123*s(i2,i3)) - 
     -  (two*zb(i2,i4)*(za(i3,i2)*zb(i2,i1) + za(i3,i4)*zb(i4,i1)))/
     -   (T234*s(i2,i3))

      amp(1,3)=
     &    (two*zb(i1,i4)*(za(i2,i1)*zb(i3,i1) - za(i2,i4)*zb(i4,i3)))/
     -   (T124*s(i1,i2)) + (two*
     -     (za(i2,i1)*zb(i1,i3)*zb(i4,i1) + 
     -       za(i3,i2)*zb(i3,i1)*zb(i4,i3)))/(T123*s(i1,i2))

      amp(2,3)=
     &   (two*zb(i3,i4)*(za(i2,i3)*zb(i1,i3) - za(i2,i4)*zb(i4,i1)))/
     -   (T234*s(i2,i3)) + (two*
     -     (za(i1,i2)*zb(i1,i3)*zb(i4,i1) + 
     -       za(i2,i3)*zb(i3,i1)*zb(i4,i3)))/(T123*s(i2,i3))

      amp(1,4)=
     &   (two*zb(i3,i2)*(za(i2,i1)*zb(i4,i2) + za(i3,i1)*zb(i4,i3)))/
     -   (T123*s(i1,i2)) + (two*zb(i4,i2)*
     -     (za(i2,i1)*zb(i3,i2) - za(i4,i1)*zb(i4,i3)))/
     -   (T124*s(i1,i2))

      amp(2,4)=
     &   (two*(za(i3,i1)*zb(i2,i3) + za(i4,i1)*zb(i2,i4))*zb(i3,i4))/
     -   (T134*s(i1,i4)) + (two*zb(i2,i4)*
     -     (-(za(i2,i1)*zb(i2,i3)) + za(i4,i1)*zb(i3,i4)))/
     -   (T124*s(i1,i4))


       

      BBBBH=0._dp                                        
      do l=1,4                                         
      BBBBH=BBBBH+abs(amp(1,l))**2+abs(amp(2,l))**2
     &     +two/xn*real(amp(1,l)*conjg(amp(2,l)))     
      enddo                                            


      BBBBH=two*BBBBH*V/four                            

      END
