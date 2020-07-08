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
 
      subroutine qqb_dm_monojet_Samps(p,i1,i2,i3,i4,i5,amp) 
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
      real(dp):: bp,beta,s34
      complex(dp):: amp_prod(2,2),amp_dec(2,2) 
      integer:: h1,h2,h3,h4
      complex(dp):: s123 
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
!====== dm decay 

      s34=Dble(za(i4,i5)*zb(i5,i4))
      beta=sqrt(1d0-4d0*xmass**2/s34) 
      bp=0.5d0*(one+beta)
      
!      write(6,*) 'In LO ',i4,i5,bp,s34,xmass
      call dm_scal_decay(i4,i5,za,zb,bp,amp_dec) 

      s123=za(i1,i2)*zb(i2,i1)+za(i2,i3)*zb(i3,i2)+za(i1,i3)*zb(i3,i1) 

      amp_prod(1,1)=-s123/(zb(i1,i2)*zb(i2,i3))
      amp_prod(2,2)=-s123/(za(i1,i2)*za(i2,i3))
      amp_prod(1,2)=-za(i1,i3)**2/(za(i1,i2)*za(i2,i3))
      amp_prod(2,1)=-zb(i1,i3)**2/(zb(i1,i2)*zb(i2,i3))

      do h1=1,2
         do h2=1,2
            do h3=1,2
               do h4=1,2 
                  amp(h1,h2,h3,h4)=amp_prod(h1,h2)*amp_dec(h3,h4) 
               enddo
            enddo
         enddo
      enddo

      
      
      
 
      return 
      end 


      subroutine dm_scal_decay(i2,i1,za,zb,bp,amp_dec) 
      implicit none
      include 'types.f'
       
      include 'dm_params.f' 
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f' 
      integer:: i2,i1 
      real(dp):: bp 
      complex(dp):: amp_dec(2,2) 

      amp_dec(1,2)=czip
      amp_dec(2,1)=czip
      amp_dec(1,1)=za(i2,i1)*bp-xmass**2/(zb(i1,i2)*bp)
      amp_dec(2,2)=zb(i2,i1)*bp-xmass**2/(za(i1,i2)*bp)

      return 
      end 
      
