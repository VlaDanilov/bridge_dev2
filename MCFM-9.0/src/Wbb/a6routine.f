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
 
      subroutine a6routine(st,j1,j2,j3,j4,j5,j6,za,zb,a6sf,a6tp,a6uv) 
      implicit none
      include 'types.f'
************************************************************************
*     Author: R.K. Ellis                                               *
*     July, 1998.                                                      *
************************************************************************
*     a6sf is the sum of a6s and a6f 
*     it is only the sum which is needed.
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'masses.f'
      include 'epinv.f'
      include 'toploops.f'
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp)::atree,virtsf,virtuv,virttp,Lnrat,a6sf,a6tp,a6uv,
     & tree,A6texact
      character*2 st 
      logical:: msbar
      common/msbar/msbar
      
      if (st == 'sl') then
      write(6,*) 'error in a6routine',st
      stop
      endif
      tree=atree(st,j1,j2,j3,j4,j5,j6,za,zb)
      
      virttp=czip
      
      if (toplight) then
        if     (toploops == tapprox) then
          if (mt .ne. 0._dp) then
            virttp=-2._dp/15._dp*s(j2,j3)/mt**2
          else
            stop 'mt=0 in a6routine'
          endif
        elseif (toploops == texact) then
          virttp=A6texact(s(j2,j3),mt**2)
        elseif (toploops == tnone) then
          virttp=czip
        endif
      else
        virttp=czip
      endif
      
! Check that asymptotic limit is okay
!      do k=1,1000
!      s(j2,j3)=float(k)*1000._dp
!      write(6,*) sqrt(s(j2,j3)),-2._dp/15._dp*s(j2,j3)/mt**2/A6texact(s(j2,j3),mt**2)
!      enddo
!      stop
      
      virtsf=two/three*epinv
     & +two/three*Lnrat(musq,-s(j2,j3))+10._dp/9._dp
      virtuv=(epinv*(11._dp-two/xn*real(nf,dp))-one)/three
      if (msbar) virtuv=virtuv+two*CF/xn
c---virtuv is the infinite and finite renormalization 
c---to get us to MS bar system
c--the term commented is associated with the number of legs
c  and in this program is taken care of in factorization procedure.
      a6sf=tree*virtsf 
      a6tp=tree*virttp 
      a6uv=tree*virtuv 
      
      return
      end

 
