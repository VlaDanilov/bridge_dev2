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
 
      !complex(dp) :: w(29)

      w(1)=N**(-1)
      w(2)=OneMz**(-1)
      w(3)=y**(-1)
      w(4)=log(1 - xv - xu)
      w(5)=log(1 + xu*xv**(-1))
      w(6)=log(1 - xu)
      w(7)=log(xu)
      w(8)=log(xv)
      w(9)=li2(1 - xv**(-1))
      w(10)=li2(1 - xv**(-1) + xu*xv**(-1))
      w(11)=li2( - xu*xv**(-1))
      w(12)=yPz**(-1)
      w(13)=w(7) - w(8)
      w(14)=w(8) + w(5)
      w(13)=w(13)*w(14)
      w(15)=w(14) + 3._dp/2._dp
      w(16)=pi*im
      w(15)=w(15)*w(16)
      w(13)=w(13) + w(15) + 3._dp/2._dp*w(4)
      w(15)=w(4) - w(8)
      w(17)=w(15) + w(16)
      w(18)=1._dp/2._dp*w(2)
      w(19)=w(18) + 1
      w(19)=w(19)*w(17)
      w(19)=w(19) + 1._dp/2._dp
      w(19)=w(19)*w(2)
      w(20)=w(17)*w(6)
      w(21)=w(10) + w(11)
      w(19)=w(19) - w(20) - w(21) + 7._dp/2._dp
      w(22)=w(19) - w(13)
      w(20)=w(20) - w(9) + w(10)
      w(23)=1._dp/2._dp*w(4)
      w(24)=w(23) + w(20) + 1._dp/2._dp*w(16)
      w(25)= - z*w(24)
      w(26)=w(2) - 1
      w(26)=w(26)*w(17)
      w(27)=1._dp/2._dp*w(3)
      w(27)=w(20)*w(27)
      w(28)= - z**2*w(27)
      w(25)=w(28) + 1._dp/2._dp*w(26) + w(25)
      w(25)=w(3)*w(25)
      w(22)=w(25) + 1._dp/2._dp*w(22) + w(9)
      w(22)=w(1)*w(22)
      w(15)= - w(5) + w(15) + w(7)
      w(15)=w(15)*im
      w(25)=pi*im**2
      w(15)=w(15) + w(25)
      w(15)=w(15)*pi
      w(25)=w(7) - 3._dp/4._dp
      w(25)=w(25)*w(4)
      w(14)=w(14) - 3._dp/4._dp
      w(14)=w(14)*w(7)
      w(28)=w(8)*w(5)
      w(29)=Pi**2
      w(14)=w(15) + w(28) + w(25) - w(14) + 1._dp/6._dp*w(29)
      w(15)= - w(19) - w(14)
      w(19)=1._dp/2._dp*N
      w(15)=w(15)*w(19)
      w(25)=Pi*im
      w(23)=w(23) - w(25) + 1._dp/2._dp*w(7) + w(16) - w(8)
      w(25)=1._dp/2._dp*beta0
      w(23)=w(23)*w(25)
      w(15)=w(15) - w(23) + w(22)
      w(18)=w(18) + w(6)
      w(18)=w(18)*w(17)
      w(18)=w(18) + w(21) - 3
      w(14)=w(18) - w(14)
      w(14)=w(14)*w(19)
      w(21)=w(12)*w(8)
      w(22)=w(21) - 1
      w(22)=w(22)*z*w(12)
      w(13)=w(22) - w(21) + w(18) + w(13)
      w(18)= - z + 1
      w(18)=w(18)*w(27)
      w(18)=w(18) - w(24)
      w(18)=z*w(18)
      w(18)=1._dp/2._dp*w(8) + w(18) + w(20)
      w(18)=w(3)*w(18)
      w(13)=w(18) + w(9) - 1._dp/2._dp*w(13)
      w(13)=w(1)*w(13)
      w(13)=w(14) - w(23) + w(13)
      w(14)= - w(2)*w(17)
      w(17)=w(3)*z*w(20)
      w(14)=w(17) + w(14) + w(4) + w(16)
      w(14)=w(3)*w(14)
      w(16)=w(26) + 1
      w(16)=w(16)*w(2)
      w(16)=w(16) - 1
      w(14)=w(14) + w(22) - 1._dp/2._dp*w(16) - w(21)
      w(14)=w(1)*w(14)
      w(16)=w(16)*w(19)
      w(14)=w(14) + w(16)
      w(14)=1._dp/2._dp*w(14)


       c_alphai(1) = w(15)

       c_betai(1) = w(13)

       c_gammai(1) = w(14)

       c_alphai(0) =  1._dp

       c_betai(0) =  1._dp

       c_gammai(0) =  0
