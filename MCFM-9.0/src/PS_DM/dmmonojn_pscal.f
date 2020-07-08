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
 


      subroutine dmmonojn_Pscal(j1,j2,j5,p,n1,qqbg) 
      implicit none
      include 'types.f'
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      include 'dm_params.f' 
      include 'zprods_decl.f'
!----- ME squared for q(j1)+qb(j2)+x(j3)+x~(j4)+g(j5) X n_nu 
      integer:: j1,j2,j5,j3,j4 
      real(dp):: n1(4),p(mxpart,4),q(mxpart,4) 
      complex(dp)::amp(2,2,2),zab(mxpart,mxpart),zba(mxpart,mxpart)
      real(dp):: qqbg(2) 
!      real(dp):: z1jetn
!----- return as funciton of quark line helicity 
      integer:: h1,h2,h3
      complex(dp):: amp_p(2),amp_dec(2,2) 
      real(dp):: bp,beta,s34
      real(dp):: fac
      complex(dp):: s(mxpart,mxpart) 

      fac=0.25d0

      j3=3
      j4=4 

      call checkndotp(p,n1,j5)

      if(xmass>1d-8) then 
!---------generate massless phase space 
         call gen_masslessvecs(p,q,j3,j4)
!---------generate spinors 
         call spinoru(5,q,za,zb)
         call spinork(5,q,zab,zba,n1)          
      else
!--------massless dm can use usual spinoru
         call spinoru(5,p,za,zb)      
         call spinork(5,p,zab,zba,n1) 

      endif

      
      do h1=1,5 
         do h2=1,5
            s(h1,h2)=za(h1,h2)*zb(h2,h1)
         enddo
      enddo

!bp = 1/2(1+beta) 
! beta = sqrt(1-4xmass**2/s34) 
      s34=Dble(za(j3,j4)*zb(j4,j3))
      beta=sqrt(1d0-4d0*xmass**2/s34)
      bp=0.5d0*(one+beta)
      
!      write(6,*) 'In gvec ',j3,j4,bp,s34,xmass
      call dm_Pscal_decay(j3,j4,za,zb,bp,amp_dec) 
!------ bulid amplitudes     
!------- helicity conserving 
      
      amp_p(1)=(za(j1,j2)*zab(j1,j1))/s(j1,j5) - 
     -  (za(j2,j5)*zab(j1,j5))/s(j1,j5) - 
     -  (za(j1,j2)*zab(j2,j2))/s(j2,j5)
     &     - (za(j1,j5)*zab(j2,j5))/s(j2,j5)
      amp_p(2)=-((zb(j2,j1)*zba(j1,j1))/s(j1,j5)) + 
     -  (zb(j5,j2)*zba(j1,j5))/s(j1,j5) + 
     -  (zb(j2,j1)*zba(j2,j2))/s(j2,j5) 
     &     + (zb(j5,j1)*zba(j2,j5))/s(j2,j5)


    
      do h1=1,2 
         do h2=1,2 
            do h3=1,2 
               amp(h1,h2,h3)=amp_p(h1)*amp_dec(h2,h3) 
            enddo
         enddo
      enddo

      
       qqbg(:)=0d0 
       do h1=1,2 
          do h2=1,2 
             do h3=1,2 
                 qqbg(h1)=qqbg(h1)+fac*abs(amp(h1,h2,h3))**2
              enddo
           enddo
        enddo

!        write(6,*) j1,j2,j5
!        write(6,*) one/16d0*abs(amp(2,2,1))**2/s(j3,j4)**2 
!     &/z1jetn(j1,j2,j5,p,n1) 
!        write(6,*) one/16d0*abs(amp(2,1,2))**2/s(j3,j4)**2 
!     &/z1jetn(j2,j1,j5,p,n1) 
!        pause 

      return 
      end 