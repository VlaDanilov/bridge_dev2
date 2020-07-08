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
 
************************************************************************
*    Author: J. Campbell, January 2014                                 *
*                                                                      *
*    Virtual amplitude for the tri-photon process                      *
*    q(i1)^- + qb(i2)^+ + gamma(i3)^+ + gamma(i4)^+ + gamma(i5)^-      *
*                                                                      *
*    Adapted from the original routine of C. Williams, March 2013      *
*                                                                      *
************************************************************************
      function virt_trigam(i1,i2,i3,i4,i5,za,zb) 
      implicit none
      include 'types.f'
      complex(dp):: virt_trigam
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      include 'sprods_com.f' 
      include 'zprods_decl.f' 
      include 'epinv.f' 
      include 'epinv2.f' 
      include 'scale.f'
      integer:: i1,i2,i3,i4,i5
      complex(dp):: Vcc,Fcc,A5LO,zab,zab2,trigam
      complex(dp):: Vpole,Alo,Lsm1,L0,L1,lnrat,l12,l25

c--- statement functions
      zab(i1,i2,i3)=+za(i1,i2)*zb(i2,i3)
      zab2(i1,i2,i3,i4)=zab(i1,i2,i4)+zab(i1,i3,i4)

      l12=lnrat(musq,-s(i1,i2))
      l25=lnrat(musq,-s(i2,i5))
      Alo=trigam(i1,i2,i3,i4,i5,za,zb)

      Vpole=-(epinv**2+epinv*l12+0.5_dp*l12**2)
     &     -(3._dp/2._dp*(epinv+l25+2._dp))
    
      virt_trigam=
     & +Vpole*Alo
    
     & +(za(i1,i3)**3*za(i2,i4)*za(i4,i5)**2
     &  -za(i1,i4)**3*za(i2,i3)*za(i3,i5)**2)
     &  /(za(i1,i3)*za(i1,i4)*za(i2,i3)*za(i2,i4)*za(i3,i4)**3)
     &  *Lsm1(-s(i3,i5),-s(i1,i2),-s(i4,i5),-s(i1,i2))
     & -(za(i1,i2)*za(i4,i5))**2/(za(i1,i3)*za(i2,i4)**3*za(i3,i4))
     &  *Lsm1(-s(i4,i5),-s(i1,i3),-s(i2,i5),-s(i1,i3))
     & +(za(i1,i2)*za(i3,i5))**2/(za(i1,i4)*za(i2,i3)**3*za(i3,i4))
     &  *Lsm1(-s(i3,i5),-s(i1,i4),-s(i2,i5),-s(i1,i4))
     & +za(i1,i5)**2/(za(i1,i4)*za(i2,i3)*za(i3,i4))
     &  *Lsm1(-s(i4,i5),-s(i2,i3),-s(i1,i5),-s(i2,i3))
     & -za(i1,i5)**2/(za(i1,i3)*za(i2,i4)*za(i3,i4))
     &  *Lsm1(-s(i3,i5),-s(i2,i4),-s(i1,i5),-s(i2,i4))
     & -(za(i1,i2)*za(i3,i5))**2
     &  /(za(i1,i3)*za(i2,i3)**2*za(i2,i4)*za(i3,i4))
     &  *Lsm1(-s(i1,i3),-s(i4,i5),-s(i1,i2),-s(i4,i5))
     & +(za(i1,i2)*za(i4,i5))**2
     &  /(za(i1,i4)*za(i2,i4)**2*za(i2,i3)*za(i3,i4))
     &  *Lsm1(-s(i1,i4),-s(i3,i5),-s(i1,i2),-s(i3,i5))
     & -za(i1,i5)**2/(za(i2,i4)*za(i1,i3)*za(i3,i4))
     &  *Lsm1(-s(i1,i2),-s(i3,i5),-s(i2,i4),-s(i3,i5))
     & +za(i1,i5)**2/(za(i2,i3)*za(i1,i4)*za(i3,i4))
     &  *Lsm1(-s(i1,i2),-s(i4,i5),-s(i2,i3),-s(i4,i5))
    
     & -za(i1,i2)*za(i2,i5)**2*zb(i3,i2)/(za(i2,i3)*za(i2,i4)**2)
     &  *L0(-s(i1,i3),-s(i4,i5))/s(i4,i5)
     & -za(i1,i2)*za(i2,i5)**2*zb(i4,i2)/(za(i2,i4)*za(i2,i3)**2)
     &  *L0(-s(i1,i4),-s(i3,i5))/s(i3,i5)
     & +za(i1,i2)*za(i4,i5)**2*zb(i4,i3)/(za(i2,i4)**2*za(i3,i4))
     &  *L0(-s(i1,i3),-s(i2,i5))/s(i2,i5)
     & +za(i1,i3)*za(i4,i5)**2*zb(i4,i3)**2/(2._dp*za(i2,i4)*za(i3,i4))
     &  *L1(-s(i1,i3),-s(i2,i5))/s(i2,i5)**2
     & +za(i1,i2)*za(i3,i5)**2*zb(i4,i3)/(za(i2,i3)**2*za(i3,i4))
     &  *L0(-s(i1,i4),-s(i2,i5))/s(i2,i5)
     & -za(i1,i4)*za(i3,i5)**2*zb(i4,i3)**2/(2._dp*za(i2,i3)*za(i3,i4))
     &  *L1(-s(i1,i4),-s(i2,i5))/s(i2,i5)**2
     & -za(i1,i2)*za(i1,i5)*za(i2,i5)/(za(i1,i3)*za(i2,i3)*za(i2,i4)**2)
     &  *lnrat(-s(i4,i5),-s(i2,i5))
     & -za(i1,i2)*za(i1,i5)*za(i2,i5)/(za(i1,i4)*za(i2,i4)*za(i2,i3)**2)
     &  *lnrat(-s(i3,i5),-s(i2,i5))
    
     & +0.5_dp*zb(i3,i4)/zb(i2,i5)*(za(i1,i5)/(za(i2,i5)*za(i3,i4))
     &   *(za(i3,i5)/za(i2,i3)+za(i4,i5)/za(i2,i4))
     &  +(zb(i2,i3)/za(i2,i4)-zb(i2,i4)/za(i2,i3))/zb(i1,i5))
            
c      write(6,*) 'virt_trigam,virt_trigam/Alo         ',
c     & virt_trigam,virt_trigam/Alo
   
      return 
      end 

