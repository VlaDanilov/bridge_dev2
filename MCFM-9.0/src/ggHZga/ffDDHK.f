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
 
      function ffDDHK(s12,s34,msq)
      implicit none
      include 'types.f'
      include 'cplx.h'
      complex(dp):: ffDDHK
      
C-----Author: R.K.Ellis July 2012
C-----The top loop function taken from
C-----Djouadi,Driesen,Hollik,Kraft arXiv:hep-ph/9701342v1, Eq(4)
      real(dp):: s12,s34,msq
      complex(dp):: C0DDHK,C2DDHK
      ffDDHK=C0DDHK(s12,s34,msq)+4._dp*C2DDHK(s12,s34,msq)
      return
      end

      function fWDDHK(s12,s34,msq)
      implicit none
      include 'types.f'
      include 'cplx.h'
      complex(dp):: fWDDHK
      
C-----Author: R.K.Ellis July 2012
C-----The W-loop function taken from
C-----Djouadi,Driesen,Hollik,Kraft arXiv:hep-ph/9701342v1, Eq(4)
c----- msq here should only be equal to wmass**2
      real(dp):: s12,s34,msq
      complex(dp):: C0DDHK,C2DDHK
      fWDDHK=2._dp*(s12/msq*(1._dp-2._dp*msq/s34)
     & +2._dp*(1._dp-6._dp*msq/s34))*C2DDHK(s12,s34,msq)
     & +4._dp*(1._dp-4._dp*msq/s34)*C0DDHK(s12,s34,msq)
      return
      end

      function f0DDHK(s12,s34)
      implicit none
      include 'types.f'
      complex(dp):: f0DDHK
C-----The F0-function taken from
C-----Djouadi,Driesen,Hollik,Kraft arXiv:hep-ph/9701342v1, Eq(2)
C-----suitably generalized to allow off-shell Z-line
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'zcouple.f'
      include 'couple.f'
      include 'kpart.f'
      include 'first.f'
      include 'msbarmasses.f'
      integer:: top
      complex(dp):: ffDDHK,fWDDHK
      real(dp):: s12,s34,cotw,mtsq,mwsq,mt_eff,massfrun
      parameter(top=2)
      save mt_eff
!$omp threadprivate(mt_eff)
                       
      if (first) then           
c--- run mt to appropriate scale
        if (kpart==klord) then
          mt_eff=massfrun(mt_msbar,hmass,amz,1)
        else
          mt_eff=massfrun(mt_msbar,hmass,amz,2)
        endif  
        first=.false.
      endif
      
      cotw=sqrt((1._dp-xw)/xw)
      mtsq=mt**2 
      mwsq=wmass**2 
      f0DDHK=s34*(cotw*fWDDHK(s12,s34,mwsq)
     &+2._dp*Q(top)*xn*mt_eff**2/s34*(L(top)+R(top))*ffDDHK(s12,s34,mtsq))
      return
      end
      

      function qlC2DDHK(s12,s34,msq)
        use mod_qcdloop_c
      implicit none
      include 'types.f'
      include 'cplx.h'
      complex(dp):: qlC2DDHK
      
      include 'scale.f'
      real(dp):: s12,s34,msq 
      qlC2DDHK=-cplx1(0.5_dp/(s12-s34))
     & +0.5_dp*s34/(s12-s34)**2
     & *(qlI2(s34,msq,msq,musq,0)-qlI2(s12,msq,msq,musq,0))
     & -msq/(s12-s34)*qlI3(s12,s34,0._dp,msq,msq,msq,musq,0)
      return
      end
c       - 1/2*[s12-s34]^-1
c          + 1/2*B0f(p1,mt,mt)*s34*[s12-s34]^-2
c          - 1/2*B0f(p12,mt,mt)*s34*[s12-s34]^-2
c          - C0DDHK(p1,p2,mt,mt,mt)*mt^2*[s12-s34]^-1


      function qlC0DDHK(s12,s34,msq)
        use mod_qcdloop_c
      implicit none
      include 'types.f'
      complex(dp):: qlC0DDHK
      
      include 'scale.f'
      real(dp):: s12,s34,msq 
      
      qlC0DDHK=qlI3(s12,s34,0._dp,msq,msq,msq,musq,0)
      
      return
      end
