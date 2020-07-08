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
 
      logical:: PDFerrors
      integer:: maxPDFsets,currentPDF
c--- 40 is my choice for the maximum number of dipoles
c--- 50 is the old choice for the maximum number of PDF error sets,
c--- that has been increased to 1000 in order to include NNPDF sets
      real(dp):: PDFxsec(0:1000),PDFxsec_nd(0:1000,0:40),
     & PDFwgt(0:1000)
      common/PDFerrors/PDFerrors,maxPDFsets,PDFxsec
      common/PDFweight/PDFwgt  
!$omp threadprivate(/PDFweight/)
