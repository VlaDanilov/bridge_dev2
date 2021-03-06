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
 

      subroutine qqb_dm_qqb(p,i1,i2,i5,i6,amp_a,amp_b) 
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      include 'dm_params.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      include 'first.f'
!------ amplitudes=> ME's for q(i1)+Qb(i5)+Q(i6)+qb(i2)+x(3)+x~(4) 
      real(dp):: p(mxpart,4),q(mxpart,4)
      integer:: i1,i2,i5,i6 
!      real(dp):: msqA,msqB 
      complex(dp):: amp_a(2,2,2,2),amp_b(2,2,2,2)
      integer:: h1,h2,h3,h4
!      complex(dp):: aqqb_zbb_new,cor
      complex(dp):: qqb_dm_qqb_VLR,qqb_dm_qqb_VRL,qqb_dm_qqb_VLL
     & ,qqb_dm_qqb_VRR
      complex(dp):: qqb_dm_qqb_Ax_VLR,qqb_dm_qqb_Ax_VRL
     & ,qqb_dm_qqb_Ax_VLL,qqb_dm_qqb_Ax_VRR
      complex(dp):: amp_pa(2,2),amp_pb(2,2),amp_dec(2,2)
      real(dp):: bp,beta,s34

!      complex(dp):: tc1(2,2,2,2),Rc2(2,2,2,2)
      
      do h1=1,2
         do h2=1,2
            do h3=1,2
               do h4=1,2
                  amp_a(h1,h2,h3,h4)=czip
                  amp_b(h1,h2,h3,h4)=czip
               enddo
            enddo
         enddo
      enddo

!      ampA(:,:,:,:)=czip 
   
      if(xmass>1d-8) then 
!---------generate massless phase space 
         call gen_masslessvecs(p,q,3,4)
!---------generate spinors 
         call spinoru(6,q,za,zb)
      else
!--------massless dm can use usual spinoru
         call spinoru(6,p,za,zb)       
      endif

      do h1=1,6
         do h2=1,6 
            s(h1,h2)=Dble(za(h1,h2)*zb(h2,h1))
         enddo
      enddo



      if(dm_mediator=='vector') then 
!----- Z => qqb L L (LR,RL,LL,RR)
         amp_a(1,1,1,2)=qqb_dm_qqb_VLR(i1,i6,i5,i2,3,4)
         amp_a(1,1,2,1)=qqb_dm_qqb_VRL(i1,i6,i5,i2,3,4)
         amp_a(1,1,1,1)=qqb_dm_qqb_VLL(i1,i6,i5,i2,3,4)
         amp_a(1,1,2,2)=qqb_dm_qqb_VRR(i1,i6,i5,i2,3,4) 
!-----Z => bb L L (LR,RL,LL,RR)
         amp_b(1,1,1,2)=qqb_dm_qqb_VLR(i6,i1,i2,i5,3,4)
         amp_b(1,1,2,1)=qqb_dm_qqb_VRL(i6,i1,i2,i5,3,4)
         amp_b(1,1,1,1)=qqb_dm_qqb_VLL(i6,i1,i2,i5,3,4)
         amp_b(1,1,2,2)=qqb_dm_qqb_VRR(i6,i1,i2,i5,3,4) 
!-----Z => qqb R R (LR,RL,LL,RR)
         amp_a(2,2,1,2)=-qqb_dm_qqb_VLR(i2,i5,i6,i1,3,4)
         amp_a(2,2,2,1)=-qqb_dm_qqb_VRL(i2,i5,i6,i1,3,4)
         amp_a(2,2,1,1)=-qqb_dm_qqb_VLL(i2,i5,i6,i1,3,4)
         amp_a(2,2,2,2)=-qqb_dm_qqb_VRR(i2,i5,i6,i1,3,4) 
!-----Z => bb R R (LR,RL,LL,RR)
         amp_b(2,2,1,2)=-qqb_dm_qqb_VLR(i5,i2,i1,i6,3,4)
         amp_b(2,2,2,1)=-qqb_dm_qqb_VRL(i5,i2,i1,i6,3,4)
         amp_b(2,2,1,1)=-qqb_dm_qqb_VLL(i5,i2,i1,i6,3,4)
         amp_b(2,2,2,2)=-qqb_dm_qqb_VRR(i5,i2,i1,i6,3,4) 
!-----Z => qqb L R (LR,RL,LL,RR)
         amp_a(1,2,1,2)=qqb_dm_qqb_VLR(i1,i5,i6,i2,3,4)
         amp_a(1,2,2,1)=qqb_dm_qqb_VRL(i1,i5,i6,i2,3,4)
         amp_a(1,2,1,1)=qqb_dm_qqb_VLL(i1,i5,i6,i2,3,4)
         amp_a(1,2,2,2)=qqb_dm_qqb_VRR(i1,i5,i6,i2,3,4) 
!-----Z => bb L R (LR,RL,LL,RR)
         amp_b(1,2,1,2)=-qqb_dm_qqb_VLR(i5,i1,i2,i6,3,4)
         amp_b(1,2,2,1)=-qqb_dm_qqb_VRL(i5,i1,i2,i6,3,4)
         amp_b(1,2,1,1)=-qqb_dm_qqb_VLL(i5,i1,i2,i6,3,4)
         amp_b(1,2,2,2)=-qqb_dm_qqb_VRR(i5,i1,i2,i6,3,4) 
!-----Z => qqb R L (LR,RL,LL,RR)
         amp_a(2,1,1,2)=-qqb_dm_qqb_VLR(i2,i6,i5,i1,3,4)
         amp_a(2,1,2,1)=-qqb_dm_qqb_VRL(i2,i6,i5,i1,3,4)
         amp_a(2,1,1,1)=-qqb_dm_qqb_VLL(i2,i6,i5,i1,3,4)
         amp_a(2,1,2,2)=-qqb_dm_qqb_VRR(i2,i6,i5,i1,3,4) 
!-----Z => bb R L (LR,RL,LL,RR)
         amp_b(2,1,1,2)=qqb_dm_qqb_VLR(i6,i2,i1,i5,3,4)
         amp_b(2,1,2,1)=qqb_dm_qqb_VRL(i6,i2,i1,i5,3,4)
         amp_b(2,1,1,1)=qqb_dm_qqb_VLL(i6,i2,i1,i5,3,4)
         amp_b(2,1,2,2)=qqb_dm_qqb_VRR(i6,i2,i1,i5,3,4) 
      elseif(dm_mediator=='axvect') then 
         !----- Z => qqb L L (LR,RL,LL,RR)
         amp_a(1,1,1,2)=qqb_dm_qqb_Ax_VLR(i1,i6,i5,i2,3,4)
         amp_a(1,1,2,1)=qqb_dm_qqb_Ax_VRL(i1,i6,i5,i2,3,4)
         amp_a(1,1,1,1)=qqb_dm_qqb_Ax_VLL(i1,i6,i5,i2,3,4)
         amp_a(1,1,2,2)=qqb_dm_qqb_Ax_VRR(i1,i6,i5,i2,3,4) 
!-----Z => bb L L (LR,RL,LL,RR)
         amp_b(1,1,1,2)=qqb_dm_qqb_Ax_VLR(i6,i1,i2,i5,3,4)
         amp_b(1,1,2,1)=qqb_dm_qqb_Ax_VRL(i6,i1,i2,i5,3,4)
         amp_b(1,1,1,1)=qqb_dm_qqb_Ax_VLL(i6,i1,i2,i5,3,4)
         amp_b(1,1,2,2)=qqb_dm_qqb_Ax_VRR(i6,i1,i2,i5,3,4) 
!-----Z => qqb R R (LR,RL,LL,RR)
         amp_a(2,2,1,2)=-qqb_dm_qqb_Ax_VLR(i2,i5,i6,i1,3,4)
         amp_a(2,2,2,1)=-qqb_dm_qqb_Ax_VRL(i2,i5,i6,i1,3,4)
         amp_a(2,2,1,1)=-qqb_dm_qqb_Ax_VLL(i2,i5,i6,i1,3,4)
         amp_a(2,2,2,2)=-qqb_dm_qqb_Ax_VRR(i2,i5,i6,i1,3,4) 
!-----Z => bb R R (LR,RL,LL,RR)
         amp_b(2,2,1,2)=-qqb_dm_qqb_Ax_VLR(i5,i2,i1,i6,3,4)
         amp_b(2,2,2,1)=-qqb_dm_qqb_Ax_VRL(i5,i2,i1,i6,3,4)
         amp_b(2,2,1,1)=-qqb_dm_qqb_Ax_VLL(i5,i2,i1,i6,3,4)
         amp_b(2,2,2,2)=-qqb_dm_qqb_Ax_VRR(i5,i2,i1,i6,3,4) 
!-----Z => qqb L R (LR,RL,LL,RR)
         amp_a(1,2,1,2)=qqb_dm_qqb_Ax_VLR(i1,i5,i6,i2,3,4)
         amp_a(1,2,2,1)=qqb_dm_qqb_Ax_VRL(i1,i5,i6,i2,3,4)
         amp_a(1,2,1,1)=qqb_dm_qqb_Ax_VLL(i1,i5,i6,i2,3,4)
         amp_a(1,2,2,2)=qqb_dm_qqb_Ax_VRR(i1,i5,i6,i2,3,4) 
!-----Z => bb L R (LR,RL,LL,RR)
         amp_b(1,2,1,2)=-qqb_dm_qqb_Ax_VLR(i5,i1,i2,i6,3,4)
         amp_b(1,2,2,1)=-qqb_dm_qqb_Ax_VRL(i5,i1,i2,i6,3,4)
         amp_b(1,2,1,1)=-qqb_dm_qqb_Ax_VLL(i5,i1,i2,i6,3,4)
         amp_b(1,2,2,2)=-qqb_dm_qqb_Ax_VRR(i5,i1,i2,i6,3,4) 
!-----Z => qqb R L (LR,RL,LL,RR)
         amp_a(2,1,1,2)=-qqb_dm_qqb_Ax_VLR(i2,i6,i5,i1,3,4)
         amp_a(2,1,2,1)=-qqb_dm_qqb_Ax_VRL(i2,i6,i5,i1,3,4)
         amp_a(2,1,1,1)=-qqb_dm_qqb_Ax_VLL(i2,i6,i5,i1,3,4)
         amp_a(2,1,2,2)=-qqb_dm_qqb_Ax_VRR(i2,i6,i5,i1,3,4) 
!-----Z => bb R L (LR,RL,LL,RR)
         amp_b(2,1,1,2)=qqb_dm_qqb_Ax_VLR(i6,i2,i1,i5,3,4)
         amp_b(2,1,2,1)=qqb_dm_qqb_Ax_VRL(i6,i2,i1,i5,3,4)
         amp_b(2,1,1,1)=qqb_dm_qqb_Ax_VLL(i6,i2,i1,i5,3,4)
         amp_b(2,1,2,2)=qqb_dm_qqb_Ax_VRR(i6,i2,i1,i5,3,4)
      elseif(dm_mediator=='scalar') then 
         call qqb_dm_qqb_Samps(p,i1,i5,i6,i2,za,zb,amp_pa)
         call qqb_dm_qqb_Samps(p,i5,i1,i2,i6,za,zb,amp_pb)
         s34=s(3,4)
         beta=sqrt(1d0-4d0*xmass**2/s34) 
         bp=0.5d0*(one+beta)
         call dm_scal_decay(3,4,za,zb,bp,amp_dec)          
         do h1=1,2
            do h2=1,2 
               do h3=1,2
                  do h4=1,2 
            amp_a(h1,h2,h3,h4)=amp_pa(h1,h2)*amp_dec(h3,h4)
            amp_b(h1,h2,h3,h4)=amp_pb(h1,h2)*amp_dec(h3,h4)
         enddo
      enddo
      enddo
      enddo
      elseif(dm_mediator=='pseudo') then 
         if(first) then 
            first = .false. 
            call check_dmAxC
         endif         

         call qqb_dm_qqb_Samps(p,i1,i5,i6,i2,za,zb,amp_pa)
         call qqb_dm_qqb_Samps(p,i5,i1,i2,i6,za,zb,amp_pb)
         s34=s(3,4)
         beta=sqrt(1d0-4d0*xmass**2/s34) 
         bp=0.5d0*(one+beta)
         
!         write(6,*) 'In QQ ',3,4,bp,s34,xmass
         call dm_Pscal_decay(3,4,za,zb,bp,amp_dec)        
         do h1=1,2
            do h2=1,2 
               do h3=1,2
                  do h4=1,2 
            amp_a(h1,h2,h3,h4)=amp_pa(h1,h2)*amp_dec(h3,h4)
            amp_b(h1,h2,h3,h4)=amp_pb(h1,h2)*amp_dec(h3,h4)
         enddo
      enddo
      enddo
      enddo
      endif
      return 
      end 

 
!--------- function for helicity conserving amplitudes 
      function qqb_dm_qqb_VLR(i1,i2,i3,i4,i5,i6) 
      implicit none
      include 'types.f'
      complex(dp):: qqb_dm_qqb_VLR
        
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      include 'zprods_com.f'
      include 'sprods_com.f'
      integer:: i1,i2,i3,i4,i5,i6
!------ copy exisitng MCFM structure     
c--- This corresponds to A++(1,2,3,4) of eq. (12.3) in BDK
c    The notation of BDK calculates the following amplitude
c
c     q3(L)----<----------q2            q3(L)------<--------q2          
c                 0                             0
c                 0                             0
c                 0                             0
c     q1(R)------<--------q4            q1(R)------<--------q4
c             )                                         )
c            (                                         (
c             )                                         )
c     l5(L)-------<-------l6            l5(L)-------<-------l6
c
c     Note that this function has the property
c     Conjg(aqqb_zbb_new(i1,i2,i3,i4,i5,i6))=
C          -aqqb_zbb_new(i4,i3,i2,i1,i6,i5)
!------ note that the default expresion (i.e. 123456) corresponds 
!------ to amp_A(1,2,1,2) and amp_B(2,1,1,2)

      qqb_dm_qqb_VLR=-((-((s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i1,i3)*
     &     za(i4,i5)*
     -         zb(i2,i1)*zb(i6,i1)) + 
     -      (s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i3,i4)*za(i3,i5)*
     -       zb(i3,i2)*zb(i6,i1) + 
     -      (s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i3,i4)*za(i4,i5)*
     -       zb(i4,i2)*zb(i6,i1) - 
     -      (s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i2,i3)*za(i4,i5)*
     -       zb(i2,i1)*zb(i6,i2))/
     -    (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -      (s(i2,i3) + s(i2,i4) + s(i3,i4))))
      
      return 
      end 

      function qqb_dm_qqb_VRL(i1,i2,i3,i4,i5,i6) 
      implicit none
      include 'types.f'
      complex(dp):: qqb_dm_qqb_VRL
        
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      include 'zprods_com.f'
      include 'sprods_com.f'
      integer:: i1,i2,i3,i4,i5,i6
!------ copy exisitng MCFM structure     
c--- This corresponds to A++(1,2,3,4) of eq. (12.3) in BDK
c    The notation of BDK calculates the following amplitude
c
c     q3(L)----<----------q2            q3(L)------<--------q2          
c                 0                             0
c                 0                             0
c                 0                             0
c     q1(R)------<--------q4            q1(R)------<--------q4
c             )                                         )
c            (                                         (
c             )                                         )
c     l5(R)-------<-------l6            l5(R)-------<-------l6
c
c     Note that this function has the property
c     Conjg(aqqb_zbb_new(i1,i2,i3,i4,i5,i6))=
C          -aqqb_zbb_new(i4,i3,i2,i1,i6,i5)
!------ note that the default expresion (i.e. 123456) corresponds 
!------ to amp_A(1,2,2,1) and amp_B(2,1,2,1)

      qqb_dm_qqb_VRL=-((-((s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i1,i3)
     &     *za(i4,i6)*
     -         zb(i2,i1)*zb(i5,i1)) + 
     -      (s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i3,i4)*za(i3,i6)*
     -       zb(i3,i2)*zb(i5,i1) + 
     -      (s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i3,i4)*za(i4,i6)*
     -       zb(i4,i2)*zb(i5,i1) - 
     -      (s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i2,i3)*za(i4,i6)*
     -       zb(i2,i1)*zb(i5,i2))/
     -    (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -      (s(i2,i3) + s(i2,i4) + s(i3,i4))))
      return 
      end 

      function qqb_dm_qqb_VLL(i1,i2,i3,i4,i5,i6) 
      implicit none
      include 'types.f'
      complex(dp):: qqb_dm_qqb_VLL
        
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      include 'dm_params.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      integer:: i1,i2,i3,i4,i5,i6
!------ copy exisitng MCFM structure     
c--- This corresponds to A++(1,2,3,4) of eq. (12.3) in BDK
c    The notation of BDK calculates the following amplitude
c
c     q3(L)----<----------q2            q3(L)------<--------q2          
c                 0                             0
c                 0                             0
c                 0                             0
c     q1(R)------<--------q4            q1(R)------<--------q4
c             )                                         )
c            (                                         (
c             )                                         )
c     l5(L)-----X-<-------l6(L)            l5(R)-------<-------l6
c
c     Note that this function has the property
c     Conjg(aqqb_zbb_new(i1,i2,i3,i4,i5,i6))=
C          -aqqb_zbb_new(i4,i3,i2,i1,i6,i5)
!------ note that the default expresion (i.e. 123456) corresponds 
!------ to amp_A(1,2,1,1) and amp_B(2,1,1,1)

      qqb_dm_qqb_VLL=(xmass*(-((s(i2,i3) + s(i2,i4) + s(i3,i4))*
     &     za(i1,i3)*za(i4,i5)*
     -         zb(i2,i1)*zb(i5,i1)) + 
     -      (s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i3,i4)*za(i3,i5)*
     -       zb(i3,i2)*zb(i5,i1) + 
     -      (s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i3,i4)*za(i4,i5)*
     -       zb(i4,i2)*zb(i5,i1) - 
     -      (s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i2,i3)*za(i4,i5)*
     -       zb(i2,i1)*zb(i5,i2) + 
     -      (s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i1,i3)*za(i4,i6)*
     -       zb(i2,i1)*zb(i6,i1) - 
     -      (s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i3,i4)*za(i3,i6)*
     -       zb(i3,i2)*zb(i6,i1) - 
     -      (s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i3,i4)*za(i4,i6)*
     -       zb(i4,i2)*zb(i6,i1) + 
     -      (s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i2,i3)*za(i4,i6)*
     -       zb(i2,i1)*zb(i6,i2)))/
     -  (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -    (s(i2,i3) + s(i2,i4) + s(i3,i4))*zb(i6,i5))
      return 
      end 

            function qqb_dm_qqb_VRR(i1,i2,i3,i4,i5,i6) 
            implicit none
      include 'types.f'
      complex(dp):: qqb_dm_qqb_VRR
        
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      include 'dm_params.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      integer:: i1,i2,i3,i4,i5,i6
!------ copy exisitng MCFM structure     
c--- This corresponds to A++(1,2,3,4) of eq. (12.3) in BDK
c    The notation of BDK calculates the following amplitude
c
c     q3(L)----<----------q2            q3(L)------<--------q2          
c                 0                             0
c                 0                             0
c                 0                             0
c     q1(R)------<--------q4            q1(R)------<--------q4
c             )                                         )
c            (                                         (
c             )                                         )
c     l5(R)-----X-<-------l6(R)            l5(R)-------<-------l6
c
c     Note that this function has the property
c     Conjg(aqqb_zbb_new(i1,i2,i3,i4,i5,i6))=
C          -aqqb_zbb_new(i4,i3,i2,i1,i6,i5)
!------ note that the default expresion (i.e. 123456) corresponds 
!------ to amp_A(1,2,2,2) and amp_B(2,1,2,2)

      qqb_dm_qqb_VRR=-((xmass*(-((s(i2,i3) + s(i2,i4) + s(i3,i4))*
     ^za(i1,i3)*
     -           za(i4,i5)*zb(i2,i1)*zb(i5,i1)) + 
     -        (s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i3,i4)*za(i3,i5)*
     -         zb(i3,i2)*zb(i5,i1) + 
     -        (s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i3,i4)*za(i4,i5)*
     -         zb(i4,i2)*zb(i5,i1) - 
     -        (s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i2,i3)*za(i4,i5)*
     -         zb(i2,i1)*zb(i5,i2) + 
     -        (s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i1,i3)*za(i4,i6)*
     -         zb(i2,i1)*zb(i6,i1) - 
     -        (s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i3,i4)*za(i3,i6)*
     -         zb(i3,i2)*zb(i6,i1) - 
     -        (s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i3,i4)*za(i4,i6)*
     -         zb(i4,i2)*zb(i6,i1) + 
     -        (s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i2,i3)*za(i4,i6)*
     -         zb(i2,i1)*zb(i6,i2)))/
     -    (s(i2,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -      (s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i5,i6)))
      return 
      end 
