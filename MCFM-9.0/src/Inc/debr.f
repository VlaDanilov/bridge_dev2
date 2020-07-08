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
 
      complex(dp):: HQ1H,HQ2H,HQ3H,HQPH,HP1H,HP2H,HP3H,HPQH,H1QH,H2QH,
     & H3QH,H1PH,H2PH,H3PH,H12H,H13H,H23H,H21H,H31H,H32H,H1LH,H2LH,H3LH,
     & HL1H,HL2H,HL3H,H1L_H,H2L_H,H3L_H,HL_1H,HL_2H,HL_3H,HL_QH,HL_PH,
     & HQL_H,HPL_H,
     & TQ1T,TQ2T,TQ3T,TQPT,TP1T,TP2T,TP3T,TPQT,T1QT,T2QT,T3QT,T1PT,T2PT,
     & T3PT,T12T,T13T,T23T,T21T,T31T,T32T,T1L_T,T2L_T,T3L_T,T1LT,T2LT,
     & T3LT,TL_1T,TL_2T,TL_3T,TL_QT,TL_PT,TLQT,TLPT,TQLT,TPLT
      common/debr1/HQ1H,HQ2H,HQ3H,HQPH,HP1H,HP2H,HP3H,HPQH,H1QH,H2QH,
     & H3QH,H1PH,H2PH,H3PH,H12H,H13H,H23H,H21H,H31H,H32H,H1LH,H2LH,H3LH,
     & HL1H,HL2H,HL3H,H1L_H,H2L_H,H3L_H,HL_1H,HL_2H,HL_3H,HL_QH,HL_PH,
     & HQL_H,HPL_H,
     & TQ1T,TQ2T,TQ3T,TQPT,TP1T,TP2T,TP3T,TPQT,T1QT,T2QT,T3QT,T1PT,T2PT,
     & T3PT,T12T,T13T,T23T,T21T,T31T,T32T,T1L_T,T2L_T,T3L_T,T1LT,T2LT,
     & T3LT,TL_1T,TL_2T,TL_3T,TL_QT,TL_PT,TLQT,TLPT,TQLT,TPLT
      real(dp):: S12,S13,S23,SQ1,SQ2,SQ3,SP1,SP2,SP3
      common/debr2/S12,S13,S23,SQ1,SQ2,SQ3,SP1,SP2,SP3
      integer:: q,k1,k2,k3,p,l,l_