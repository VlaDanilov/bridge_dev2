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
 
      t1 = 0.1D1 / sman
      t2 = mH ** 2
      t3 = 0.1D1 / uman
      t4 = cdlogwrap(tman * t1)
      t5 = cdlogwrap(tman * t3)
      t6 = cdlogwrap(uman * t1)
      t7 = tman ** 2
      t8 = cA * t7
      t9 = LogMums + LogMumt + LogMumu
      t10 = sman * uman
      c3mt0 = t7 * ((0.2D1 / 0.3D1) * LogMuMtop + 7 * cA - (0.20D2 / 0.3
     #D1) * tr) + (0.4D1 / 0.3D1) * t7 * ((dilogc(-t1 * t2 + 1) + dilogc
     #(-t2 / tman + 1) + dilogc(-t2 * t3 + 1)) * cA - LogMuMtop * tr) + 
     #t8 * (-pi ** 2 + t4 ** 2 + t5 ** 2 + t6 ** 2) / 3 + (0.11D2 / 0.9D
     #1) * t8 * t9 - (0.20D2 / 0.9D1) * tr * (t7 * t9 + t10) - 2 * t7 * 
     #(-Kg + cF) + (0.2D1 / 0.9D1) * t10 * cA

