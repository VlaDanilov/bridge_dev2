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
 
      subroutine PSZHHamps(pa,pb,pc,pd,amp)
        use mod_qcdloop_c
      implicit none
      include 'types.f'
C--    Formula taken from Plehn, Spira and Zerwas                        
C--    Nucl. Phys. B479 (1996) 46
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'scale.f'
      real(dp):: pa(4),pb(4),pc(4),pd(4),pcsq,pdsq,
     & ss,tt,uu,S,T,U,T1,T2,U1,U2,rhoc,rhod,tauQ,mQsq,lambdaHHH
      complex(dp):: Ftriangle,Fbox,Gbox,
     & Dabc,Dbac,Dacb,Cab,Cbc,Cac,Ccd,Cbd,Cad,
     & Cdelta,Cbox,amp(2)

      mQsq=mt**2
      ss=(pa(4)+pb(4))**2
     &  -(pa(1)+pb(1))**2-(pa(2)+pb(2))**2-(pa(3)+pb(3))**2
      tt=(pa(4)+pc(4))**2
     &  -(pa(1)+pc(1))**2-(pa(2)+pc(2))**2-(pa(3)+pc(3))**2
      uu=(pb(4)+pc(4))**2
     &  -(pb(1)+pc(1))**2-(pb(2)+pc(2))**2-(pb(3)+pc(3))**2
      S=ss/mQsq
      T=tt/mQsq
      U=uu/mQsq
      lambdaHHH=3._dp*(hmass/zmass)**2
      rhoc=hmass**2/mQsq
      rhod=hmass**2/mQsq
      T1=T-rhoc
      U1=U-rhoc
      T2=T-rhod
      U2=U-rhod
      tauQ=4._dp/S
C      ss=(pa+pb)**2
C      tt=(pa+pc)**2
C      uu=(pb+pc)**2
      pcsq=pc(4)**2-pc(1)**2-pc(2)**2-pc(3)**2
      pdsq=pd(4)**2-pd(1)**2-pd(2)**2-pd(3)**2
      Cab=qlI3(0._dp,0._dp,ss,mQsq,mQsq,mQsq,musq,0)
      Cbc=qlI3(0._dp,pcsq,uu,mQsq,mQsq,mQsq,musq,0)
      Cac=qlI3(0._dp,pcsq,tt,mQsq,mQsq,mQsq,musq,0)
      Ccd=qlI3(pcsq,pdsq,ss,mQsq,mQsq,mQsq,musq,0)
      Cbd=qlI3(0._dp,pdsq,tt,mQsq,mQsq,mQsq,musq,0)
      Cad=qlI3(0._dp,pdsq,uu,mQsq,mQsq,mQsq,musq,0)
      Dabc=qlI4(0._dp,0._dp,pcsq,pdsq,ss,uu,mQsq,mQsq,mQsq,mQsq,musq,0)
      Dbac=qlI4(0._dp,0._dp,pcsq,pdsq,ss,tt,mQsq,mQsq,mQsq,mQsq,musq,0)
      Dacb=qlI4(0._dp,pcsq,0._dp,pdsq,tt,uu,mQsq,mQsq,mQsq,mQsq,musq,0)


      Ftriangle=2._dp/S*(2._dp+(4._dp-S)*mQsq*Cab)

      Fbox=(4._dp*S+8._dp*S*mQsq*Cab-2._dp*S*(S+rhoc+rhod-8._dp)
     & *mQsq**2*(Dabc+Dbac+Dacb)
     & +(rhoc+rhod-8._dp)*mQsq*(T1*Cac+U1*Cbc+U2*Cad+T2*Cbd
     & -(T*U-rhoc*rhod)*mQsq*Dacb))/S**2
      Gbox=((T**2+rhoc*rhod-8._dp*T)*mQsq
     & *(S*Cab+T1*Cac+T2*Cbd-S*T*mQsq*Dbac)
     & +(U**2+rhoc*rhod-8._dp*U)*mQsq
     & *(S*Cab+U1*Cbc+U2*Cad-S*U*mQsq*Dabc)
     & -(T**2+U**2-2._dp*rhoc*rhod)*(T+U-8._dp)*mQsq*Ccd
     & -2._dp*(T+U-8._dp)*(T*U-rhoc*rhod)*mQsq**2
     & *(Dabc+Dbac+Dacb))/(S*(T*U-rhoc*rhod))

      Cdelta=lambdaHHH*zmass**2/cplx2(ss-hmass**2,hmass*hwidth) 
      Cbox=cone

      amp(1)=Cdelta*Ftriangle+Cbox*Fbox
      amp(2)=Cbox*Gbox
      return
      end
