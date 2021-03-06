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
 


      subroutine qqb_dm_gg_Samps(p,i1,i2,i3,i4,i5,i6,msq1,msq2,msqsl) 
      implicit none
      include 'types.f'
       
      include 'dm_params.f' 
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
!----- vector amplitude for 
!----- q(i1)+g(i2)+g(i3)+qb(i4)+x(i5)+x(i6) 
      real(dp):: p(mxpart,4) 
!-----fills amplitude for q g qb chi,chib  
      real(dp):: q(mxpart,4)
      integer:: i1,i2,i3,i4,i5,i6
!------and gluon helicity is summed over 
!------returns amp squared. 
!      complex(dp):: z2jet_1LL(-1:1,-1:1)
!      real(dp):: amp_tes(2,2) 
      integer:: h1,h2,h3,h4,h5
!      complex(dp):: cor_fac 
      complex(dp):: amp_hconsx1(2,2,2,2,2),amp_hconsx2(2,2,2,2,2)
      complex(dp):: amp_hconsQED(2,2,2,2,2)
      complex(dp):: amp_prod_1(2,2,2),amp_prod_2(2,2,2)
      complex(dp):: amp_dec(2,2)
      real(dp):: s34,bp,beta
!------ index refers to helicit of qqb (sum over everyting else) 
      real(dp):: msq1(2),msq2(2),msqsl(2)
      real(dp):: fac
!---- fac is relative difference w.r.t. vector operator 
      
      fac=1d0 
 
      amp_prod_1(:,:,:)=czip
      amp_prod_2(:,:,:)=czip
      
      if(xmass>1d-8) then 
!--------- generate massless phase space 
      call gen_masslessvecs(p,q,i5,i6)
!--------- generate spinors 
      call spinoru(6,q,za,zb)
      else
!-------- massless dm can use usual spinoru
         call spinoru(6,p,za,zb)       
      endif

!====== helicity amplitudes
!      call dm_gg_helamps(i1,i2,i3,i4,i5,i6,za,zb,amp_hconsx1)
!======= other colour ordering, swap gluons 
!      call dm_gg_helamps(i1,i3,i2,i4,i5,i6,za,zb,amp_hconsx2) 


      call dm_gg_Shelamps(i1,i2,i3,i4,za,zb,amp_prod_1)
      call dm_gg_Shelamps(i1,i3,i2,i4,za,zb,amp_prod_2)
!----- decay amplitudes 
      s34=Dble(za(i5,i6)*zb(i6,i5))
      beta=sqrt(1d0-4d0*xmass**2/s34) 
      bp=0.5d0*(one+beta)
!      write(6,*) 'In gg ',i5,i6,bp,s34,xmass
      call dm_scal_decay(i5,i6,za,zb,bp,amp_dec) 
      do h1=1,2
         do h2=1,2
            do h3=1,2 
               do h4=1,2 
                  do h5=1,2
        amp_hconsx1(h1,h2,h3,h4,h5)=amp_prod_1(h1,h2,h3)*amp_dec(h4,h5)
        amp_hconsx2(h1,h2,h3,h4,h5)=amp_prod_2(h1,h2,h3)*amp_dec(h4,h5)
      enddo
      enddo
      enddo
      enddo
      enddo

      
      do h1=1,2
         do h2=1,2
            do h3=1,2
               do h4=1,2 
                  do h5=1,2
               amp_hconsQED(h1,h2,h3,h4,h5)=amp_hconsx1(h1,h2,h3,h4,h5) 
     &                +amp_hconsx2(h1,h3,h2,h4,h5) 
             enddo
          enddo
       enddo
      enddo
      enddo

      msq1(:)=0d0 
      msq2(:)=0d0 
      msqsl(:)=0d0 
!------ msq 
      do h1=1,2
         do h2=1,2
            do h3=1,2 
               do h4=1,2
                  do h5=1,2 
!      do h1=1,1
!         do h2=2,2
!           do h3=2,2 
!               do h4=1,2
!                  do h5=1,2 
            msq1(h1)=msq1(h1)+fac*abs(amp_hconsx1(h1,h2,h3,h4,h5))**2
            msq2(h1)=msq2(h1)+fac*abs(amp_hconsx2(h1,h2,h3,h4,h5))**2
          msqsl(h1)=msqsl(h1)+fac*abs(amp_hconsQED(h1,h2,h3,h4,h5))**2
                  enddo
               enddo
            enddo
         enddo
      enddo
!----- include -1/N_c^2 for subleading piece 
      do h1=1,2
         msqsl(h1)=-one/9d0*msqsl(h1)
      enddo

      return 
      end

      
      subroutine dm_gg_Shelamps(i1,i2,i3,i4,za,zb,amp)
      implicit none
      include 'types.f'
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      include 'zprods_decl.f' 
      integer:: i1,i2,i3,i4
      complex(dp):: amp(2,2,2) 
      complex(dp):: s1234 
      complex(dp):: s(mxpart,mxpart)
!-------q(i4)+qb(i1)+g(i2)+g(i3) 
      integer:: h1,h2 
       
      do h1=1,6
         do h2=1,6
            s(h1,h2)=za(h1,h2)*zb(h2,h1) 
         enddo
      enddo


      s1234=za(i1,i2)*zb(i2,i1)+za(i1,i3)*zb(i3,i1)+za(i1,i4)*zb(i4,i1)
     &     +za(i2,i3)*zb(i3,i2)+za(i2,i4)*zb(i4,i2)+za(i3,i4)*zb(i4,i3) 

      amp(1,1,1)=s1234/(zb(i2,i1)*zb(i3,i2)*zb(i4,i3))     
      amp(1,1,2)= -((za(i1,i2)*za(i2,i4)**2)/
     -     ((s(i2,i3) + s(i2,i4) + 
     -         s(i3,i4))*za(i2,i3)*
     -       za(i3,i4))) - 
     -  za(i2,i4)**2/
     -   (za(i2,i3)*za(i3,i4)*
     -     zb(i2,i1)) + 
     -  (za(i1,i2)*za(i2,i4)*
     -     zb(i3,i1))/
     -   ((s(i1,i2) + s(i1,i3) + 
     -       s(i2,i3))*za(i2,i3)*
     -     zb(i2,i1)) - 
     -  (za(i1,i4)*za(i2,i4)*
     -     zb(i3,i1))/
     -   (za(i2,i3)*za(i3,i4)*
     -     zb(i2,i1)*zb(i3,i2)) + 
     -  (za(i1,i2)*za(i1,i4)*
     -     zb(i3,i1)**2)/
     -   ((s(i1,i2) + s(i1,i3) + 
     -       s(i2,i3))*za(i2,i3)*
     -     zb(i2,i1)*zb(i3,i2)) + 
     -  (za(i1,i4)*za(i2,i4)**2*
     -     zb(i4,i3))/
     -   ((s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))*za(i2,i3)*
     -     za(i3,i4)*zb(i3,i2))
      amp(1,2,1)= -((za(i1,i3)**2*za(i3,i4))/
     -     ((s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -       za(i1,i2)*za(i2,i3))) + 
     -  (za(i1,i3)**2*za(i1,i4)*zb(i2,i1))/
     -   ((s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     za(i1,i2)*za(i2,i3)*zb(i3,i2)) - 
     -  za(i1,i3)**2/
     -   (za(i1,i2)*za(i2,i3)*zb(i4,i3)) + 
     -  (za(i1,i3)*za(i3,i4)*zb(i4,i2))/
     -   ((s(i2,i3) + s(i2,i4) + s(i3,i4))*
     -     za(i2,i3)*zb(i4,i3)) - 
     -  (za(i1,i3)*za(i1,i4)*zb(i4,i2))/
     -   (za(i1,i2)*za(i2,i3)*zb(i3,i2)*
     -     zb(i4,i3)) + 
     -  (za(i1,i4)*za(i3,i4)*zb(i4,i2)**2)/
     -   ((s(i2,i3) + s(i2,i4) + s(i3,i4))*
     -     za(i2,i3)*zb(i3,i2)*zb(i4,i3))
      amp(1,2,2)=-(za(i1,i4)**2/(za(i1,i2)*za(i2,i3)*za(i3,i4)))

      amp(2,1,1)=zb(i4,i1)**2/
     -  (zb(i2,i1)*zb(i3,i2)*
     -    zb(i4,i3))
      amp(2,1,2)= zb(i3,i1)**2/
     -   (za(i3,i4)*zb(i2,i1)*
     -     zb(i3,i2)) + 
     -  (za(i2,i4)*zb(i3,i1)*
     -     zb(i4,i1))/
     -   (za(i2,i3)*za(i3,i4)*
     -     zb(i2,i1)*zb(i3,i2)) - 
     -  (za(i1,i2)*zb(i3,i1)**2*
     -     zb(i4,i1))/
     -   ((s(i1,i2) + s(i1,i3) + 
     -       s(i2,i3))*za(i2,i3)*
     -     zb(i2,i1)*zb(i3,i2)) - 
     -  (za(i2,i4)*zb(i3,i1)*
     -     zb(i4,i3))/
     -   ((s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))*za(i3,i4)*
     -     zb(i3,i2)) + 
     -  (zb(i3,i1)**2*zb(i4,i3))/
     -   ((s(i1,i2) + s(i1,i3) + 
     -       s(i2,i3))*zb(i2,i1)*
     -     zb(i3,i2)) - 
     -  (za(i2,i4)**2*zb(i4,i1)*
     -     zb(i4,i3))/
     -   ((s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))*za(i2,i3)*
     -     za(i3,i4)*zb(i3,i2))
      amp(2,2,1)=-((za(i1,i3)**2*zb(i2,i1)*
     -       zb(i4,i1))/
     -     ((s(i1,i2) + s(i1,i3) + 
     -         s(i2,i3))*za(i1,i2)*
     -       za(i2,i3)*zb(i3,i2))) - 
     -  (za(i1,i3)*zb(i2,i1)*
     -     zb(i4,i2))/
     -   ((s(i1,i2) + s(i1,i3) + 
     -       s(i2,i3))*za(i1,i2)*
     -     zb(i3,i2)) + 
     -  (za(i1,i3)*zb(i4,i1)*
     -     zb(i4,i2))/
     -   (za(i1,i2)*za(i2,i3)*
     -     zb(i3,i2)*zb(i4,i3)) + 
     -  zb(i4,i2)**2/
     -   (za(i1,i2)*zb(i3,i2)*
     -     zb(i4,i3)) + 
     -  (zb(i2,i1)*zb(i4,i2)**2)/
     -   ((s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))*zb(i3,i2)*
     -     zb(i4,i3)) - 
     -  (za(i3,i4)*zb(i4,i1)*
     -     zb(i4,i2)**2)/
     -   ((s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))*za(i2,i3)*
     -     zb(i3,i2)*zb(i4,i3))
      amp(2,2,2)=-(s1234/(za(i1,i2)*za(i2,i3)*za(i3,i4)))
      return 
      end 

