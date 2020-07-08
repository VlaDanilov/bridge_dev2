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
 


      subroutine qqb_dm_monojet_lc_PSamps(p,i1,i2,i3,i4,i5,amp) 
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      include 'dm_params.f' 
      include 'zprods_decl.f' 
      include 'scale.f'
      include 'epinv.f'
!----- LEADING COLOR 
!----- fills amplitude for q g qb chi,chib 
      complex(dp):: amp(2,2,2,2) 
      complex(dp):: amp_tree(2,2,2,2) 
      real(dp):: p(mxpart,4),q(mxpart,4) 
      integer:: i1,i2,i3,i4,i5 
      integer:: h1,h2,h3,h4 
      complex(dp):: Lsm1
      real(dp):: s23,s123,s34
      complex(dp):: lnrat
      real(dp):: s(mxpart,mxpart)
      complex(dp):: vfac,ffac,rat(2,2),amp_dec(2,2) 
      integer:: j,k
      real(dp):: uvsub
      real(dp):: s13,beta,bp
!---- ps point for test 

!      p(i1,4)=2.18799500337960849d0
!      p(i1,1)=-1.75847163017272490d0
!      p(i1,2)=-0.01742318231181945d0
!      p(i1,3)=-1.30184334441972760d0

!      p(i3,4)=1.55171692881328013d0
!      p(i3,1)= 1.18390845497861889d0
!      p(i3,2)=0.98045174296920807d0
!      p(i3,3)=-0.21189756276205669d0
      
!      p(i2,4)=1.052445768129343192d0
!      p(i2,1)=0.75255756964896313d0
!      p(i2,2)=-0.46085308444383691d0
!      p(i2,3)=0.573509924740636861d0

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
      s(:,:)=0d0
      do j=1,5 
         do k=1,5
            s(j,k)=Dble(za(j,k)*zb(k,j))
         enddo
      enddo
!------ basis integrals 
!      l12=lnrat(musq,-s(i1,i2))
!      l23=lnrat(musq,-s(i2,i3))

      s123=s(i1,i2)+s(i2,i3)+s(i1,i3)
      s13=s(i1,i2) ! s13 = s_qg 
      s23=s(i2,i3) ! s23 = s_gqb

!------ univerisal v and f pieces
      vfac=-2d0*epinv**2+epinv*(-lnrat(musq,-s13) 
     &     - lnrat(musq,-s23))
     & -0.5d0*(lnrat(musq,-s13)**2+lnrat(musq,-s23)**2)

      ffac=-Lsm1(-s23,-s123,-s13,-s123)
!----- additional uv-counterterm from vertex renorm 
      uvsub=-1.5d0*epinv
 !------helicity dependent rational pieces 
      
!      atreet=s123/(zb(i2,i1)*zb(i3,i2))
      rat(1,1)=0.5d0*(s13+s23)/(zb(i2,i1)*zb(i3,i2))
      rat(2,2)=0.5d0*(s13+s23)/(za(i2,i1)*za(i3,i2))
      rat(1,2)=czip
      rat(2,1)=czip
!      write(6,*) atreet,abs(rt2*atreet) 
!      write(6,*) vfac+ffac,abs(two*(vfac+ffac))
!      write(6,*) rat(1,1)*rt2*two,abs(rat(1,1)*rt2*two)
!      write(6,*) abs(rt2*two*((vfac+ffac)*atreet+rat(1,1)))
!      pause 


     
!====== scalar decay 
      s34=Dble(za(i4,i5)*zb(i5,i4))
      beta=sqrt(1d0-4d0*xmass**2/s34) 
      bp=0.5d0*(one+beta)
      
      call dm_Pscal_decay(i4,i5,za,zb,bp,amp_dec) 
      
!======= tree (note swap on gluon quark in these amplitudes) 
      call qqb_dm_monojet_PSamps(p,i1,i2,i3,i4,i5,amp_tree) 


      do h1=1,2 
         do h2=1,2 
            do h3=1,2 
               do h4=1,2 
                amp(h1,h2,h3,h4)=(vfac+ffac+uvsub)*amp_tree(h1,h2,h3,h4)
     &                 +rat(h1,h2)*amp_dec(h3,h4) 
               enddo
            enddo
         enddo
      enddo
      
      return 
      end 
