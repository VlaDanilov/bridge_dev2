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
 
      function fc002(msq,taugs,qsq)
      implicit none
      include 'types.f'
      include 'constants.f'
      complex(dp):: fc002
      
      real(dp):: qsq,msq,qsqhat,taugs,qgs,cgs2
      complex(dp):: lnrat,I3me

      qsqhat=qsq-msq
      qgs=qsqhat-taugs
      cgs2=(qsqhat-taugs)**2-four*taugs*msq

      fc002 =
     &  - 1.D0/12.D0*qgs*cgs2**(-1)*msq
     &  + 1.D0/6.D0*cgs2**(-1)*msq*taugs
     &  + 1.D0/12.D0*cgs2**(-1)*msq**2*taugs*qsq**(-1)
     &  + 1.D0/12.D0*cgs2**(-1)*msq**2
     &  - 1.D0/12.D0*cgs2**(-1)*msq**3*qsq**(-1)
     &  - I3me(msq,taugs,qsq)*cgs2**(-2)*msq**2*taugs**3
     &  + 1.D0/2.D0*lnrat( - taugs,msq)*qgs*cgs2**(-2)*msq*taugs**2
     &  - 1.D0/12.D0*lnrat( - taugs,msq)*qgs*cgs2**(-1)*taugs
     &  - 1.D0/2.D0*lnrat( - qsqhat,msq)*qgs*cgs2**(-2)*msq*taugs**2*
     & qsq**(-1)*qsqhat
     &  + 1.D0/12.D0*lnrat( - qsqhat,msq)*qgs*cgs2**(-1)*msq*qsq**(-1)*
     & qsqhat
     &  + 1.D0/12.D0*lnrat( - qsqhat,msq)*qgs*cgs2**(-1)*taugs*
     & qsq**(-1)*qsqhat
     &
      fc002 = fc002 - lnrat( - qsqhat,msq)*cgs2**(-2)*msq**2*taugs**2*
     & qsq**(-1)*qsqhat
     &  - 1.D0/12.D0*lnrat( - qsqhat,msq)*cgs2**(-1)*msq**2*taugs*
     & qsq**(-2)*qsqhat
     &  - 1.D0/12.D0*lnrat( - qsqhat,msq)*cgs2**(-1)*msq**2*qsq**(-1)*
     & qsqhat
     &  + 1.D0/12.D0*lnrat( - qsqhat,msq)*cgs2**(-1)*msq**3*qsq**(-2)*
     & qsqhat
     &

      return
      end

      function fc00(msq,taugs,qsq)
      implicit none
      include 'types.f'
      include 'constants.f'
      complex(dp):: fc00
      
      real(dp):: qsq,msq,qsqhat,taugs,qgs,cgs2
      complex(dp):: lnrat,I3me

      qsqhat=qsq-msq
      qgs=qsqhat-taugs
      cgs2=(qsqhat-taugs)**2-four*taugs*msq

      fc00 =
     &  - 1.D0/2.D0*I3me(msq,taugs,qsq)*cgs2**(-1)*msq*taugs**2
     &  + 1.D0/4.D0*lnrat( - taugs,msq)*qgs*cgs2**(-1)*taugs
     &  - 1.D0/4.D0*lnrat( - qsqhat,msq)*qgs*cgs2**(-1)*taugs*qsq**(-1)
     & *qsqhat
     &  - 1.D0/2.D0*lnrat( - qsqhat,msq)*cgs2**(-1)*msq*taugs*qsq**(-1)
     & *qsqhat
     &

      return
      end
