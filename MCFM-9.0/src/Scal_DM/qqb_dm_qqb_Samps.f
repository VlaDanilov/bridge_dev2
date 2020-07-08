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
 
      
      subroutine qqb_dm_qqb_Samps(p,i1,i2,i3,i4,za,zb,amp)
      implicit none
      include 'types.f'
       
!---------- amplitudes for q(i1)+Q(i2)+Qb(i3)+q(i4) 
!---------- with i1 i4 coupling to DM 
      include 'dm_params.f' 
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      include 'zprods_decl.f'
      real(dp):: p(mxpart,4),q(mxpart,4) 
      integer:: i1,i2,i3,i4
      complex(dp):: amp(2,2) 
      real(dp):: s(mxpart,mxpart) 
      integer:: h1,h2 
      real(dp):: fac 
    

!====== fac w.r.t. vector !--- note in gg we defined fac as 
!------- at ME level, here we are at amp so sqrt(fac) 
      fac=sqrt(1d0) 

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

      

      amp(1,1)= (za(i1,i2)*(za(i1,i4)*zb(i3,i1) + za(i2,i4)*zb(i3,i2)))/
     -   ((s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i2,i3)*zb(i3,i2)) - 
     -  (za(i2,i4)*(-(za(i1,i2)*zb(i3,i2)) + za(i1,i4)*zb(i4,i3)))/
     -   ((s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i2,i3)*zb(i3,i2))
      amp(1,2)=(za(i1,i3)*(za(i1,i4)*zb(i2,i1) - za(i3,i4)*zb(i3,i2)))/
     -   ((s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i2,i3)*zb(i3,i2)) - 
     -  (za(i3,i4)*(za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2)))/
     -   ((s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i2,i3)*zb(i3,i2))
      amp(2,1)= ((za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))*zb(i4,i3))/
     -   ((s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i2,i3)*zb(i3,i2)) - 
     -  (zb(i3,i1)*(za(i1,i2)*zb(i4,i1) - za(i2,i3)*zb(i4,i3)))/
     -   ((s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i2,i3)*zb(i3,i2))
      amp(2,2)=((-(za(i2,i3)*zb(i2,i1)) + za(i3,i4)*zb(i4,i1))
     &     *zb(i4,i2))/
     -   ((s(i2,i3) + s(i2,i4) + s(i3,i4))*za(i2,i3)*zb(i3,i2)) - 
     -  (zb(i2,i1)*(za(i1,i3)*zb(i4,i1) + za(i2,i3)*zb(i4,i2)))/
     -   ((s(i1,i2) + s(i1,i3) + s(i2,i3))*za(i2,i3)*zb(i3,i2))

      do h1=1,2
         do h2=1,2 
            amp(h1,h2)=fac*amp(h1,h2) 
         enddo
      enddo


      return 
      end 
