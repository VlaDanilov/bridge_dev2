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
 

      complex*32 function hjetmass_bubble_pppm_mhsq
     &     (i1,i2,i3,i4,za,zb,mt,p,flip)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          
          double precision mt
          real*16 p(mxpart,4)
          complex*32 hjetmass_bubble_pppm_s234
          complex*32 hjetmass_bubble_pppm_s134
          complex*32 hjetmass_bubble_pppm_s124
          complex*32 hjetmass_bubble_pppm_s34
          complex*32 hjetmass_bubble_pppm_s14
          logical flip

      hjetmass_bubble_pppm_mhsq = 
     & -hjetmass_bubble_pppm_s234(i1,i2,i3,i4,za,zb,mt,p,flip)
     & -hjetmass_bubble_pppm_s134(i1,i2,i3,i4,za,zb,mt,p,flip)
     & -hjetmass_bubble_pppm_s124(i1,i2,i3,i4,za,zb,mt,p,flip)
     & -hjetmass_bubble_pppm_s34(i1,i2,i3,i4,za,zb,mt)
     & -hjetmass_bubble_pppm_s14(i1,i2,i3,i4,za,zb,mt)
      return

      end function
