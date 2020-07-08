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
 
      integer:: Nbub,NTri,NBox
      parameter(Nbub=13,NTri=3,NBox=28)
      integer:: b15,b25,b16,b26,b125,b126,b12,b34,b56,b156,b256,
     & b234,b134,b345,b346,b12a,b12b
      parameter(b15=1,b25=2,b16=3,b26=4,b125=5,b126=6,
     & b12=7,b34=8,b56=9,b156=10,b256=11,b12a=12,b12b=13) 
      parameter(b234=b156,b134=b256,b345=b126,b346=b125)
      integer:: 
     & c34_15_26,
     & c16_25_34,
     & c34_56_12
      parameter(      
     & c34_15_26=1,
     & c16_25_34=2,
     & c34_56_12=3
     & )
      integer:: 
     & d125x15_12,
     & d125x25_12,
     & d125x12_25,
     & d126x12_26,
     & d126x16_12,
     & d126x26_12,
     & d256x25_56,
     & d256x26_56,
     & d156x15_56,
     & d156x16_56,
     & d15_34x125_156,
     & d25_34x125_256,
     & d16_34x126_156,
     & d26_34x126_256,
     & d12_34x125_126,
     & d12_34x126_125,
     & d25_16x126_156,
     & d26_15x125_156,
     & d15_26x126_256,
     & d16_34x125_156,
     & h15_34x125_26,
     & h25_34x125_16,
     & h16_34x126_25,
     & h26_34x126_15,
     & h12_34x125_56,
     & h12_34x126_56,
     & h56_34x12_156,
     & h56_34x12_256

      parameter(
     & d125x15_12=1,
     & d125x25_12=2,
     & d125x12_25=3,
     & d126x12_26=4,
     & d126x16_12=5,
     & d126x26_12=6,
     & d256x25_56=7,
     & d256x26_56=8,
     & d156x15_56=9,
     & d156x16_56=10,
     & d15_34x125_156=11,
     & d25_34x125_256=12,
     & d16_34x126_156=13,
     & d26_34x126_256=14,
     & d12_34x125_126=15,
     & d12_34x126_125=16,
     & d25_16x126_156=17,
     & d26_15x125_156=18,
     & d15_26x126_256=19,
     & d16_34x125_156=20,
     & h15_34x125_26=21,
     & h25_34x125_16=22,
     & h16_34x126_25=23,
     & h26_34x126_15=24,
     & h12_34x125_56=25,
     & h12_34x126_56=26,
     & h56_34x12_156=27,
     & h56_34x12_256=28
     & )
