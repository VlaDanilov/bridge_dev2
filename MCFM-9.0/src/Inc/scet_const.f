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
 
      real(dp),parameter::
     & zeta3=1.202056903159594285399738161511449990764986292_dp,
     & zeta2=pi**2/6._dp,
     & be0 = 11/three*CA-4/three*TR*NF,
     & be1= 34/three*CA**2-(20/three*CA+4*CF)*TR*NF,
     & Ga0 = four,
     & Ga1 = 4/three*((four-pisq)*CA+5*be0),
     & Gaq0 = Ga0*CF,
     & Gaq1 = Ga1*CF,
     & Gag0 = Ga0*CA,
     & Gag1 = Ga1*CA,
     & gBq0 = 6*CF,
     & gBq1 = CF*(CA*(146/nine-80*zeta3)+CF*(three-4*pisq+48*zeta3)
     &           +be0*(121/nine+2*pisq/three)),
     & gBg0 = 2*be0,
     & gBg1 = CA*(CA*(182/nine-32*zeta3)
     & +be0*(94/nine-2/three*pisq))+2*be1,
     & gams1=CA*(-64/nine+28*zeta3)+be0*(-56/nine+2*zeta2)
