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
 
      function I3m(s1,s2,s3)
      implicit none
      include 'types.f'
      complex(dp):: I3m
C     This is the function I3m, a massless triangle with all three external 
C     lines offshell defined in BDK
C     %\cite{Bern:1997sc}
C     \bibitem{Bern:1997sc}
C     Z.~Bern, L.~J.~Dixon and D.~A.~Kosower,
C     %``One-loop amplitudes for e+ e- to four partons,''
C     Nucl.\ Phys.\ B {\bf 513}, 3 (1998)
C     [arXiv:hep-ph/9708239].
C     %%CITATION = HEP-PH 9708239;%%
C     defined in their equation II.9
C     \int da_1 da_2 da_3 /(-a_1*a_2*s1-a_2*a_3*s2-a_3*a_1*s3)
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: s1,s2,s3,smax,smid,smin,del3,rtdel3
      real(dp):: i3m1a,flag
      complex(dp):: i3m1b

      smax=max(s1,s2,s3)
      smin=min(s1,s2,s3)
      smid=s1+s2+s3-smax-smin
      del3=s1**2+s2**2+s3**2-two*(s1*s2+s2*s3+s3*s1)      

      if (del3 > 0) then
      rtdel3=sqrt(del3)
         if (smax < 0) then
c---case all negative
             flag=0._dp
             i3m=i3m1b(smax,smid,smin,rtdel3,flag)
         elseif (smin > 0) then
c---case all positive
             flag=0._dp
             i3m=-i3m1b(-smin,-smid,-smax,rtdel3,flag)
         elseif ((smid < 0) .and. (smin < 0)) then
c---case two negative and one positive
             flag=+1._dp
             i3m=i3m1b(smin,smid,smax,rtdel3,flag)
         elseif ((smax > 0).and.(smid > 0)) then
c---case two positive and one negative
             flag=-1._dp
             i3m=-i3m1b(-smax,-smid,-smin,rtdel3,flag)
         endif
      elseif (del3 < 0) then 
      rtdel3=sqrt(-del3)
         if (smax < 0) then
c---case all negative
             i3m=+cplx1(i3m1a(+s1,+s2,+s3,rtdel3))
         elseif (smin > 0) then
c---case all positive
             i3m=-cplx1(i3m1a(-s1,-s2,-s3,rtdel3))  
          endif
      endif

      return
      end     


      function I3m1a(s1,s2,s3,rtmdel)
      implicit none
      include 'types.f'
      real(dp):: I3m1a
      
C     symmetric form of Lu and Perez
C     %\cite{Lu:1992ny}
c     \bibitem{Lu:1992ny}
c     H.~J.~Lu and C.~A.~Perez,
c     %``Massless one loop scalar three point integral and associated Clausen,
c     %Glaisher and L functions,''
c     SLAC-PUB-5809
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: s1,s2,s3,d1,d2,d3,rtmdel,arg1,arg2,arg3,dclaus

      d1=s1-s2-s3
      d2=s2-s3-s1
      d3=s3-s1-s2
      
      arg1=two*atan(rtmdel/d1)
      arg2=two*atan(rtmdel/d2)
      arg3=two*atan(rtmdel/d3)
      i3m1a=two/rtmdel*(Dclaus(arg1)+Dclaus(arg2)+Dclaus(arg3))

      end

      
      function I3m1b(s1,s2,s3,rtdel,flag)
      implicit none
      include 'types.f'
      complex(dp):: I3m1b
      
C     form of Ussyukina and Davydychev
C %\cite{Usyukina:1994iw}
C \bibitem{Usyukina:1994iw}
C   N.~I.~Usyukina and A.~I.~Davydychev,
C   %``New results for two loop off-shell three point diagrams,''
C  Phys.\ Lett.\ B {\bf 332}, 159 (1994)
C  [arXiv:hep-ph/9402223].
C  %%CITATION = HEP-PH 9402223;%%

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: s1,s2,s3,d3,temp,ddilog,xlog,ylog,rat
      real(dp):: x,y,rho,rtdel,argx,argy,argdlx,argdly,flag
      d3=s3-s1-s2
      x=s1/s3
      y=s2/s3
      rat=0.5_dp*(d3+rtdel)/s3
      if (abs(rat) < 1.e-3_dp) rat=2._dp*s1*s2/(s3*(d3-rtdel))
      rho=1._dp/rat
      argx=rho*x
      argy=rho*y
      argdlx=-argx
      argdly=-argy

      if ((argdlx > 1._dp) .or. (argdly > 1._dp)) then
      write(6,*) 'problems with call of I3m1b'
      write(6,*) 'argdlx',argdlx
      write(6,*) 'argdly',argdly
      stop
      endif

      xlog=log(abs(argx))
      ylog=log(abs(argy))
      temp=xlog*ylog+pisq/3._dp+(ylog-xlog)*log((one+argy)/(one+argx))
     & +two*(ddilog(argdlx)+ddilog(argdly))
      I3m1b=cplx1(temp-abs(flag)*pisq)
     & +impi*cplx1(flag*(xlog+ylog))
      I3m1b=-I3m1b/cplx1(rtdel)
      end





