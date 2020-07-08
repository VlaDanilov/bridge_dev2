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
 


      subroutine qqb_dm_monojet_Vamps(p,i1,i2,i3,i4,i5,amp) 
      implicit none
      include 'types.f'
       
      include 'dm_params.f' 
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f' 
      real(dp):: p(mxpart,4) 
!----- fills amplitude for q g qb chi,chib 
      complex(dp):: amp(2,2,2,2) 
      real(dp):: q(mxpart,4)
      integer:: i1,i2,i3,i4,i5
!------ z1jet amplitude for testing (note order is q qb l1 l2 g) 
!------ and gluon helicity is summed over 
!------ returns amp squared. 

      amp(:,:,:,:)=czip
      if(xmass>1d-8) then 
!--------- generate massless phase space 
      call gen_masslessvecs(p,q,i4,i5)
!--------- generate spinors 
      call spinoru(5,q,za,zb)
      else
!-------- massless dm can use usual spinoru
         call spinoru(5,p,za,zb) 
         
      endif

      
!------- helicity conserving amplitudes (can be checked against Z) 
      amp(1,1,1,2)=(za(i2,i4)*zb(i5,i3))/zb(i2,i1) + 
     -  (za(i1,i4)*zb(i3,i1)*zb(i5,i3))/(zb(i2,i1)*zb(i3,i2))
      amp(1,2,1,2)=(za(i1,i4)*zb(i5,i2))/za(i2,i3) + 
     -  (za(i1,i3)*za(i1,i4)*zb(i5,i3))/(za(i1,i2)*za(i2,i3))
      amp(1,1,2,1)=(za(i2,i5)*zb(i4,i3))/zb(i2,i1) + 
     -  (za(i1,i5)*zb(i3,i1)*zb(i4,i3))/(zb(i2,i1)*zb(i3,i2))
      amp(1,2,2,1)=(za(i1,i5)*zb(i4,i2))/za(i2,i3) + 
     -  (za(i1,i3)*za(i1,i5)*zb(i4,i3))/(za(i1,i2)*za(i2,i3))
      amp(2,1,1,2)= -((za(i2,i4)*zb(i5,i1))/zb(i3,i2)) - 
     -  (za(i3,i4)*zb(i3,i1)*zb(i5,i1))/(zb(i2,i1)*zb(i3,i2))
      amp(2,2,1,2)=-((za(i1,i3)*za(i3,i4)*zb(i5,i1))
     &/(za(i1,i2)*za(i2,i3))) - (za(i3,i4)*zb(i5,i2))/za(i1,i2)
      amp(2,1,2,1)=-((za(i2,i5)*zb(i4,i1))/zb(i3,i2)) - 
     -  (za(i3,i5)*zb(i3,i1)*zb(i4,i1))/(zb(i2,i1)*zb(i3,i2))
      amp(2,2,2,1)=-((za(i1,i3)*za(i3,i5)*zb(i4,i1))
     &/(za(i1,i2)*za(i2,i3))) - (za(i3,i5)*zb(i4,i2))/za(i1,i2)

!------ helicity violating amplitudes 
      amp(1,1,1,1)=-((xmass*za(i2,i4)*zb(i4,i3))/(zb(i2,i1)*zb(i5,i4)))- 
     -  (xmass*za(i1,i4)*zb(i3,i1)*zb(i4,i3))/
     -   (zb(i2,i1)*zb(i3,i2)*zb(i5,i4)) + 
     -  (xmass*za(i2,i5)*zb(i5,i3))/(zb(i2,i1)*zb(i5,i4)) + 
     -  (xmass*za(i1,i5)*zb(i3,i1)*zb(i5,i3))/
     -   (zb(i2,i1)*zb(i3,i2)*zb(i5,i4))
      amp(1,2,1,1)=-((xmass*za(i1,i4)*zb(i4,i2))/(za(i2,i3)*zb(i5,i4)))- 
     -  (xmass*za(i1,i3)*za(i1,i4)*zb(i4,i3))/
     -   (za(i1,i2)*za(i2,i3)*zb(i5,i4)) + 
     -  (xmass*za(i1,i5)*zb(i5,i2))/(za(i2,i3)*zb(i5,i4)) + 
     -  (xmass*za(i1,i3)*za(i1,i5)*zb(i5,i3))/
     -   (za(i1,i2)*za(i2,i3)*zb(i5,i4))
      amp(2,1,1,1)=(xmass*za(i2,i4)*zb(i4,i1))/(zb(i3,i2)*zb(i5,i4)) + 
     -  (xmass*za(i3,i4)*zb(i3,i1)*zb(i4,i1))/
     -   (zb(i2,i1)*zb(i3,i2)*zb(i5,i4)) - 
     -  (xmass*za(i2,i5)*zb(i5,i1))/(zb(i3,i2)*zb(i5,i4)) - 
     -  (xmass*za(i3,i5)*zb(i3,i1)*zb(i5,i1))/
     -   (zb(i2,i1)*zb(i3,i2)*zb(i5,i4))
      amp(2,2,1,1)=(xmass*za(i1,i3)*za(i3,i4)*zb(i4,i1))/
     -   (za(i1,i2)*za(i2,i3)*zb(i5,i4)) + 
     -  (xmass*za(i3,i4)*zb(i4,i2))/(za(i1,i2)*zb(i5,i4)) - 
     -  (xmass*za(i1,i3)*za(i3,i5)*zb(i5,i1))/
     -   (za(i1,i2)*za(i2,i3)*zb(i5,i4)) - 
     -  (xmass*za(i3,i5)*zb(i5,i2))/(za(i1,i2)*zb(i5,i4))
      amp(1,1,2,2)=(xmass*za(i2,i4)*zb(i4,i3))/(za(i4,i5)*zb(i2,i1)) + 
     -  (xmass*za(i1,i4)*zb(i3,i1)*zb(i4,i3))/
     -   (za(i4,i5)*zb(i2,i1)*zb(i3,i2)) - 
     -  (xmass*za(i2,i5)*zb(i5,i3))/(za(i4,i5)*zb(i2,i1)) - 
     -  (xmass*za(i1,i5)*zb(i3,i1)*zb(i5,i3))/
     -   (za(i4,i5)*zb(i2,i1)*zb(i3,i2))
      amp(1,2,2,2)=(xmass*za(i1,i4)*zb(i4,i2))/(za(i2,i3)*za(i4,i5)) + 
     -  (xmass*za(i1,i3)*za(i1,i4)*zb(i4,i3))/
     -   (za(i1,i2)*za(i2,i3)*za(i4,i5)) - 
     -  (xmass*za(i1,i5)*zb(i5,i2))/(za(i2,i3)*za(i4,i5)) - 
     -  (xmass*za(i1,i3)*za(i1,i5)*zb(i5,i3))/
     -   (za(i1,i2)*za(i2,i3)*za(i4,i5))
      amp(2,1,2,2)=-((xmass*za(i2,i4)*zb(i4,i1))/(za(i4,i5)*zb(i3,i2)))- 
     -  (xmass*za(i3,i4)*zb(i3,i1)*zb(i4,i1))/
     -   (za(i4,i5)*zb(i2,i1)*zb(i3,i2)) + 
     -  (xmass*za(i2,i5)*zb(i5,i1))/(za(i4,i5)*zb(i3,i2)) + 
     -  (xmass*za(i3,i5)*zb(i3,i1)*zb(i5,i1))/
     -   (za(i4,i5)*zb(i2,i1)*zb(i3,i2))
      amp(2,2,2,2)= -((xmass*za(i1,i3)*za(i3,i4)*zb(i4,i1))/
     -     (za(i1,i2)*za(i2,i3)*za(i4,i5))) - 
     -  (xmass*za(i3,i4)*zb(i4,i2))/(za(i1,i2)*za(i4,i5)) + 
     -  (xmass*za(i1,i3)*za(i3,i5)*zb(i5,i1))/
     -   (za(i1,i2)*za(i2,i3)*za(i4,i5)) + 
     -  (xmass*za(i3,i5)*zb(i5,i2))/(za(i1,i2)*za(i4,i5))

     
!------ define square 
!      amp2=0d0 
!      do h1=1,2 
!         do h2=1,2 
!            do h3=1,2 
!               do h4=1,2 
!                  amp2=amp2+abs(dm(L)*amp(1,h2,h3,h4))**2
!          
!               enddo
!            enddo
!         enddo
!      enddo

!      write(6,*) amp2 
!     
!---- check against z 
!      do h1=1,2
!         do h2=1,2
!            amp_tes(h1,h2)=abs(amp(h1,1,h2,3-h2))**2
!     &           +abs(amp(h1,2,h2,3-h2))**2
!         enddo
!      enddo

!      write(6,*) amp_tes(1,1)**2+amp_tes(1,2)**2
!     &+amp_tes(2,1)**2+amp_tes(2,2)**2

!      write(6,*) i1,i2,i3,i4,i5,amp_tes
!      write(6,*) abs(amp(1,1,1,2))**2+abs(amp(1,2,1,2))**2

!      call zgamps2(i1,i3,i5,i4,i2,za,zb,z1jet_tes) 
!      write(6,*) '1, 1',z1jet_tes(1,1) 
!      write(6,*) '2, 2',z1jet_tes(2,2)
!      write(6,*) '2, 1',z1jet_tes(2,1)
!      write(6,*) '1, 2',z1jet_tes(1,2)
!      write(6,*) h1,h2
!      write(6,*) ' 1, 1 mine  ',amp_tes(1,1)
!      write(6,*) ' 1, 2 mine  ',amp_tes(1,2)     
!      write(6,*) ' 2, 1 mine  ',amp_tes(2,1)
!      write(6,*) ' 2, 2 mine  ',amp_tes(2,2)     
!      pause
      return 
      end
