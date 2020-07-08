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
 
      subroutine getvbfpoint(p)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: p(mxpart,4)
      
c--- Stand-alone Madgraph PS point SQRTS=1000d0
      p(           1 ,4)=  -500.00000000000000d0   
      p(           1 ,1)=  -0.0000000000000000d0   
      p(           1 ,2)=  -0.0000000000000000d0   
      p(           1 ,3)=  -500.00000000000000d0   
      p(           2 ,4)=  -500.00000000000000d0   
      p(           2 ,1)=  -0.0000000000000000d0   
      p(           2 ,2)=  -0.0000000000000000d0   
      p(           2 ,3)=   500.00000000000000d0   
      p(           3 ,4)=   85.531224838488669d0   
      p(           3 ,1)=  -8.2219322397786812d0   
      p(           3 ,2)=   36.163783768203281d0   
      p(           3 ,3)=  -77.072504800241362d0   
      p(           4 ,4)=   181.42881161004260d0   
      p(           4 ,1)=  -57.859982948193725d0   
      p(           4 ,2)=  -171.86373408663516d0   
      p(           4 ,3)=  -5.6118589848131073d0   
      p(           5 ,4)=   82.849301077435555d0   
      p(           5 ,1)=  -65.909547623589134d0   
      p(           5 ,2)=  -49.895215719628702d0   
      p(           5 ,3)=   5.5141336005866366d0   
      p(           6 ,4)=   381.47038530081528d0   
      p(           6 ,1)=   190.18527704151884d0   
      p(           6 ,2)=   292.04294098458701d0   
      p(           6 ,3)=  -155.11330013659799d0   
      p(           7 ,4)=   54.231407011799924d0   
      p(           7 ,1)=  -31.133016208179807d0   
      p(           7 ,2)=  -7.9279665679113993d0   
      p(           7 ,3)=   43.691282361116336d0   
      p(           8 ,4)=   214.48887016141768d0   
      p(           8 ,1)=  -27.060798021777508d0   
      p(           8 ,2)=  -98.519808378614982d0   
      p(           8 ,3)=   188.59224795994947d0   
      
      return
      end
      
