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
 
      function realint(vector,wgt)
          use PDFerrors
#ifdef HAVE_LHAPDF
        use LHAPDF, only: getalphas
#endif
          use Scalevar
          use ieee_arithmetic
          use VVconfig_m
          use singletop2_m, only : singletop2_gs, singletop2_tree, singletop2_gs_light, singletop2_gs_heavy
          use singletop2_realamps_m, only : singletop2_real, singletop2_real_light, singletop2_real_heavy
          use singletop2_scale_m
          use bbfrac_m, only : bbfrac
          use MCFMStorage
          use SCET
          use m_gencuts, only : enable_reweight_user, reweight_user
          use Superhisto, only : shtmpreset, shtmpcommit
      implicit none
      real(dp):: realint
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'debug.f'
      include 'noglue.f'
      include 'nflav.f'
      include 'vegas_common.f'
      include 'ptilde.f'
      include 'phasemin.f'
      include 'new_pspace.f'
      include 'npart.f'
      include 'scale.f'
      include 'facscale.f'
      include 'dynamicscale.f'
      include 'realwt.f'
      include 'maxwt.f'
      include 'kprocess.f'
      include 'masses.f'
      include 'dipolescale.f'
      include 'stopscales.f'
      include 'flags.f'
      include 'frag.f'
      include 'ipsgen.f'
      include 'decay1q2a.f'
      include 'outputoptions.f'
      include 'breit.f'
      include 'dm_params.f' 
      include 'runstring.f'
      include 'energy.f'
      include 'incldip.f'
      include 'nproc.f'
      include 'initialscales.f'
      include 'badpoint.f'
      include 'nqcdjets.f'
      include 'taucut.f'
      include 'qcdcouple.f'
      include 'scalevar.f'
      include 'couple.f'
      include 'nlooprun.f'
      include 'cutoff.f'
      include 'ewcorr.f'

      integer:: ih1,ih2,j,k,nd,nmax,nmin,nvec,ii,t
      integer itrial
      real(dp):: alphas,msqtrial,
     & fx1up(-nf:nf),fx2up(-nf:nf),fx1dn(-nf:nf),fx2dn(-nf:nf),
     & dipfx1up(0:maxd,-nf:nf),dipfx2up(0:maxd,-nf:nf),
     & dipfx1dn(0:maxd,-nf:nf),dipfx2dn(0:maxd,-nf:nf),
     & xmsqvar(2,0:maxd),asorig,scaleup,scaledn
      real(dp):: vector(mxdim),W,ptmp, xmsqsum
      real(dp) :: centralval(0:maxd)
      real(dp):: fx1(-nf:nf),fx2(-nf:nf),
     & dipfx1(0:maxd,-nf:nf),dipfx2(0:maxd,-nf:nf),
     & fx1_H(-nf:nf),fx2_H(-nf:nf),fx1_L(-nf:nf),fx2_L(-nf:nf)
      real(dp):: p(mxpart,4),pjet(mxpart,4),p1ext(4),p2ext(4)
      real(dp):: pswt
      real(dp):: s(mxpart,mxpart),wgt,msq(-nf:nf,-nf:nf)
      ! to distinguish for singletop2 between light and heavy line emissions
      real(dp) :: msq_light(-nf:nf,-nf:nf), msq_heavy(-nf:nf,-nf:nf)
      real(dp) :: msqc_light(maxd,-nf:nf,-nf:nf), msqc_heavy(maxd,-nf:nf,-nf:nf)

      real(dp):: msqc(maxd,-nf:nf,-nf:nf),xmsq(0:maxd)
      ! ,msqc_new(maxd,-nf:nf,-nf:nf),bit1,bit2
      real(dp):: flux,BrnRat
      real(dp):: xx1,xx2,q(mxpart,4)
      real(dp):: m3,m4,m5,R,Rbbmin,dot
      real(dp):: xmsqjk, plo(mxpart,4),pswtdip
      integer:: sgnj,sgnk
      common/Rbbmin/Rbbmin
      logical:: bin,failed
      logical:: includedipole,includereal
      real(dp):: QandGint
      external qqb_w2jet_g,qqb_w2jet_gs,qqb_z2jet_g,qqb_z2jet_gs,
     & qqb_w2jet,qqb_w1jet_gs,qqb_z2jet,qqb_z1jet_gs,qqb_Hg_g,qqb_Hg_gs,
     & qqb_hww_g,qqb_hww_gs,qqb_zbb_g,qqb_zbb_gs,
     & qqb_wh_ww,qqb_wh_ww_gs,
     & qqb_wh_zz,qqb_wh_zz_gs,
     & qqb_wbb_g,qqb_wbb_gs,
     & qqb_dirgam_g,qqb_dirgam_gs,qqb_hflgam_g,qqb_hflgam_gs,
     & qqb_trigam_g,qqb_trigam_gs,qqb_gmgmjt_g,qqb_gmgmjt_gs,
     & qqb_w_g,qqb_w_gs,qqb_z1jet,qqb_z_gs,qqb_ww_g,qqb_ww_gs,
     & qqb_wz_g,qqb_wz_gs,qqb_zz_g,qqb_zz_gs,qqb_wgam_g,qqb_wgam_gs,
     & qqb_QQb_g,qqb_QQb_gs,
     & VV_Hqq_g,VV_Hqq_gs,VV_HWW_g,VV_HWW_gs,
     & gg_Hg,gg_H_gs,
     & gg_HWWgg,gg_HWWg_gs,gg_HZZgg,gg_HZZg_gs,
     & gg_Hgg,gg_Hg_gs,
     & gQ_zQ_g,gQ_zQ_gs,qqb_tbb_g,qqb_tbb_gs,
     & qqb_w_tndk_g,qqb_w_tndk_gs,
     & qqb_w_twdk_g,qqb_w_twdk_gs,qqb_w_twdk_gdk,qqb_w_twdk_gsdk,
     & qqb_zbjet_g,qqb_zbjet_gs,qqb_w_cjet_g,qqb_w_cjet_gs,
     & qqb_wbfromc_g,qqb_wbfromc_gs,
     & gg_hggg,gg_hgg_gs,qq_tchan_htqg,qq_tchan_htq_gs,
     & qg_tbq_g,qg_tbq_gs,qq_tbg_g,qq_tbg_gs,epem3j_g,epem3j_gs,
     & qq_tchan_ztqg_dk,qq_tchan_ztq_dk_gs,
     & qq_tchan_ztqg,qq_tchan_ztq_gs,
     & qqb_QQbdk_g,qqb_QQbdk_gs,qqb_gamgam_g,qqb_gamgam_gs,
     & qqb_zaj_gs,qqb_zaj_g,
     & qqb_gmgmjt_g_mad,
     & qqb_zaa_g,qqb_zaa_gs,
     & qqb_waa_g,qqb_waa_gs,
     & qqb_WH1jet_g,qqb_WH1jet_gs,
     & qqb_ZH1jet_g,qqb_ZH1jet_gs,
     & qqb_waa_g_mad,qqb_tottth_g,qqb_tottth_g_mad
      common/density/ih1,ih2
      common/bin/bin
      common/Pext/p1ext,p2ext
!$omp threadprivate(/Pext/)
      common/nmax/nmax
      common/BrnRat/BrnRat
      common/nmin/nmin

      external :: qqb_zgam_new, qqb_zgam

      QandGflag = .false.
      p(1:npart+2,:) = 0._dp
      pjet(1:npart+2,:) = 0._dp
      pswt = 0._dp
      realint = 0._dp      
c--- ensure isolation code does not think this is fragmentation piece
      z_frag=0._dp

      W=sqrts**2

      currentPDF = 0

c-- note: new_pspace now signifies multi-channel integration
      if (new_pspace) then
        call gen_lops(vector,plo,pswt,*999)
!        call writeout(plo)
        npart=npart+1
        call multichan(vector(ndim-2),vector(ndim-1),vector(ndim),
     &                 vector(ndim+2),plo,p,pswtdip,*999)

        pswt=pswt*pswtdip
      else
        call gen_realps(vector,p,pswt,*999)
      endif

      if (all(.not. ieee_is_nan(p(1:npart+2,4))) .eqv. .false.) then
          if (debug) then
              write(6,*) 'Discarding NaN or infinite phase space point'
          endif
          goto 999
      endif

      nvec=npart+2
      call dotem(nvec,p,s)
      
c----calculate the x's for the incoming partons from generated momenta
      xx1=two*(p(1,4)*p2ext(4)-p(1,3)*p2ext(3))/W
      xx2=two*(p(2,4)*p1ext(4)-p(2,3)*p1ext(3))/W

c      if (debug) write(*,*) 'Reconstructed x1,x2 ',xx1,xx2
      
      if ((xx1 >  1._dp) .or. (xx2 >  1._dp)
     &.or.(xx1 < xmin) .or. (xx2 < xmin)) then
         goto 999
      endif

c--- (moved to includedipole) impose cuts on final state
c      call masscuts(p,*999)

! small safety cuts
      call smallnew(p,npart,*999)
     
c--- extra cut to divide WQj/ZQj regions
      if ( (nproc == 312) .or. (nproc == 317)
     & .or.(nproc == 322) .or. (nproc == 327)
     & .or.(nproc == 342) .or. (nproc == 352)) then
        if (R(p,5,6) < Rbbmin) goto 999
      endif

c--- see whether this point will pass cuts - if it will not, do not
c--- bother calculating the matrix elements for it, instead set to zero
      includereal=includedipole(0,p)
      incldip(0)=includereal

c     it is always a good idea to initialize all elements as zero
      msq(:,:) = 0
      msqc(:,:,:) = 0
 
      if ((kcase==kqg_tbq) .or. (kcase==k4ftwdk)) then
          msqLH(:,:) = 0._dp
          msqHL(:,:) = 0._dp
      endif
      
c--- test to see whether we need Gflag and Qflag together
      if ( ((kcase==kW_2jet) .or. (kcase==kZ_2jet))
     &.and. (Qflag) .and. (Gflag) ) then
        QandGflag=.true.
        QandGint=0._dp
c--- first pass: Gflag
        Gflag=.true.
        Qflag=.false.
      endif
      
c--- restart from here when calculating with Qflag and Gflag
c--- (W+2 jet and Z+2 jet processes only)
   44 continue   

  771 continue
      if (dynamicscale) then
          call scaleset(initscale,initfacscale,p)
          dipscale(0)=facscale
      else
          call usescales(initscale,initfacscale)
          dipscale(0)=facscale
      endif

      if (doPDFAlphas) then
          call updateAlphas(scale)
      endif

      if (kcase==ktopanom .and. use_DDIS) then
        ! we store an array for dipole modified factorization scales, dipole pdfs and pdfs here, which need to be reset
        call singletop2_scale_reset()
      endif
      
      if ((doScalevar .and. currentPDF == 0 .and. bin) .and. (foundpow .eqv. .false.)) then
        itrial=1
      endif
   66 continue

c--- Calculate the required matrix elements
      if     (kcase==kW_only) then
c        call singcheck(qqb_w_g,qqb_w_gs,p)         ! Checked 11/30/01
        if (includereal) call qqb_w_g(p,msq)
        call qqb_w_gs(p,msqc)
      elseif (kcase==kW_1jet) then
c        call singcheck(qqb_w2jet,qqb_w1jet_gs,p)   ! Checked 11/16/01
        if (includereal) call qqb_w2jet(p,msq)
        call qqb_w1jet_gs(p,msqc)  
      elseif (kcase==kWgamma) then
!         if(includereal) call singcheck(qqb_wgam_g,qqb_wgam_gs,p) ! Checked 08/27/02
        if (includereal) call qqb_wgam_g(p,msq)
        call qqb_wgam_gs(p,msqc)  
      elseif (kcase==kWbfrmc) then
        if (includereal) call qqb_wbfromc_g(p,msq)
        call qqb_wbfromc_gs(p,msqc)  
      elseif (kcase==kW_cjet) then
c        call singcheck(qqb_w_cjet_g,qqb_w_cjet_gs,p) ! Checked 15/05/07
        if (includereal) call qqb_w_cjet_g(p,msq)
        call qqb_w_cjet_gs(p,msqc)
      elseif (kcase==kZgamma) then
        if (decayChannel() == decayQuarks) then
          !if(includereal) call singcheck(qqb_zaj_vdecay,qqb_zgam_gs,p)
          if (includereal) call qqb_zaj_vdecay(p,msq)
          call qqb_zgam_gs(p,msqc)
        else
          !if(includereal) call singcheck(qqb_zaj,qqb_zgam_gs,p)
          call set_anomcoup(p)
          if (includereal) call qqb_zaj(p,msq)
          call qqb_zgam_gs(p,msqc)
        endif
      elseif (kcase==kZ_2gam) then
c         if(includereal) call singcheck(qqb_zaa_g,qqb_zaa_gs,p)
        if (includereal) call qqb_zaa_g(p,msq)
        call qqb_zaa_gs(p,msqc)
      elseif (kcase==kW_2gam) then
c        if(includereal) call singcheck(qqb_waa_g,qqb_waa_gs,p)
c        call compare_madgraph(p,qqb_waa_g,qqb_waa_g_mad)
         stop
c        if (includereal) call qqb_waa_g(p,msq)      
c        call qqb_waa_gs(p,msqc)
      elseif (kcase==kZgajet) then
          if (decayChannel() == decayQuarks) then
c             if(includereal) call singcheck(qqb_zaj_g_vdecay,qqb_zaj_gs_vdecay,p)
              if (includereal) call qqb_zaj_g_vdecay(p,msq)
              call qqb_zaj_gs_vdecay(p,msqc)
          else
              !if(includereal) call singcheck(qqb_zaj_g,qqb_zaj_gs,p)
              call set_anomcoup(p)
              if (includereal) call qqb_zaj_g(p,msq)
              call qqb_zaj_gs(p,msqc)
          endif
      elseif (kcase==kWbbmas) then
c        call singcheck(qqb_wbbm_g,qqb_wbbm_gs,p)     ! Checked 10/21/10
      if (includereal) call qqb_wbbm_g(p,msq)      
        call qqb_wbbm_gs(p,msqc)      
      elseif (kcase==kWttmas) then
c        call singcheck(qqb_wbbm_g,qqb_wbbm_gs,p)     ! Checked 10/21/10
      if (includereal) call qqb_wbbm_g(p,msq)      
        call qqb_wbbm_gs(p,msqc)      
      elseif (kcase==kWbbbar) then
c        call singcheck(qqb_wbb_g,qqb_wbb_gs,p)     ! Checked 11/30/01
        if (includereal) call qqb_wbb_g(p,msq)      
        call qqb_wbb_gs(p,msqc)      
      elseif (kcase==kW_2jet) then  
c        if (includereal) call singcheck(qqb_w2jet_g,qqb_w2jet_gs_new,p) ! Re-checked June 09
        if (includereal)  call qqb_w2jet_g(p,msq)
        call qqb_w2jet_gs_new(p,msqc)
      elseif (kcase==kZ_only) then
c        call singcheck(qqb_z1jet,qqb_z_gs,p)         ! Checked 11/30/01
        if (includereal) call qqb_z1jet(p,msq)      
        call qqb_z_gs(p,msqc)     
      elseif (kcase==kZ_1jet) then
c        call singcheck(qqb_z2jet,qqb_z1jet_gs,p)   ! Checked 11/16/01
        if (includereal) call qqb_z2jet(p,msq)      
        call qqb_z1jet_gs(p,msqc)  
      elseif (kcase==kZ_2jet) then
c        if (includereal) call singcheck(qqb_z2jet_g,qqb_z2jet_gs_new,p) ! Re-checked June 09
        if (includereal) call qqb_z2jet_g(p,msq)
        call qqb_z2jet_gs_new(p,msqc)
c        call qqb_z2jet_gs(p,msqc) 
c        write(6,*) 'Gflag,Qflag',Gflag,Qflag
c        do j=-nf,nf
c        do k=-nf,nf
c        bit1=0._dp
c        bit2=0._dp
c        bit1=sum(msqc(1:ndmax,j,k))
c        bit2=sum(msqc_new(1:ndmax,j,k))
c        if ((abs(bit1) .gt. 1.e-100_dp) .or. (abs(bit2) .gt. 1.e-100_dp)) then
c          if (abs(bit1/bit2-1._dp) .gt. 1.e-12_dp) then
c            write(6,*) 'j,k,bit1,bit1/bit2',j,k,bit1,bit2,bit1/bit2-1._dp
c          endif
c        endif
c        enddo
c        enddo
c        pause
      elseif (kcase==kZbbbar) then
c        call singcheck(qqb_zbb_g,qqb_zbb_gs,p)     ! Checked 11/30/01
        if (includereal) call qqb_zbb_g(p,msq)
        call qqb_zbb_gs(p,msqc) 
      elseif (kcase==kWWqqbr) then
c        call singcheck(qqb_ww_g,qqb_ww_gs,p)       ! Checked 11/30/01
        if (includereal) call qqb_ww_g(p,msq)
        call qqb_ww_gs(p,msqc)      
      elseif (kcase==kWWqqdk) then
        if (includereal) call dkqqb_ww_g(p,msq)
        call dkqqb_ww_gs(p,msqc)      
      elseif (kcase==kWZbbar) then
c        call singcheck(qqb_wz_g,qqb_wz_gs,p)       ! Checked 12/05/01
        if (includereal) call qqb_wz_g(p,msq)
        call qqb_wz_gs(p,msqc)      
      elseif (kcase==kZZlept) then
c        call singcheck(qqb_zz_g,qqb_zz_gs,p)       ! Checked 12/05/01
        if (includereal) call qqb_zz_g(p,msq)
        call qqb_zz_gs(p,msqc)      
      elseif (kcase==kWHbbar) then
c        call singcheck(qqb_wh_g,qqb_wh_gs,p)
        if (includereal) call qqb_wh_g(p,msq)
        call qqb_wh_gs(p,msqc)
      elseif (kcase==kWHbbdk) then
        if (includereal) call dkqqb_wh_g(p,msq)
        call dkqqb_wh_gs(p,msqc)
      elseif (kcase==kWH1jet) then
c     call singcheck(qqb_WH1jet_g,qqb_WH1jet_gs,p)
         if (toponly) then
            msq=zip
            msqc=zip
         else
            if (includereal) call qqb_WH1jet_g(p,msq)
            call qqb_WH1jet_gs(p,msqc)
         endif
      elseif (kcase==kWH__WW) then
c        call singcheck(qqb_wh_ww_g,qqb_wh_ww_gs,p)
        if (includereal) call qqb_wh_ww_g(p,msq)
        call qqb_wh_ww_gs(p,msqc)
      elseif (kcase==kWH__ZZ) then
c        call singcheck(qqb_wh_zz_g,qqb_wh_zz_gs,p)
        if (includereal) call qqb_wh_zz_g(p,msq)
        call qqb_wh_zz_gs(p,msqc)
      elseif (kcase==kWHgaga) then
c        call singcheck(qqb_wh_gaga_g,qqb_wh_gaga_gs,p)
        if (includereal) call qqb_wh_gaga_g(p,msq)
        call qqb_wh_gaga_gs(p,msqc) 
      elseif (kcase==kZHbbar) then
c        call singcheck(qqb_zh_g,qqb_zh_gs,p)
        if (includereal) call qqb_zh_g(p,msq)
        call qqb_zh_gs(p,msqc)     
      elseif (kcase==kZHbbdk) then
         if (includereal) call dkqqb_zh_g(p,msq)
         call dkqqb_zh_gs(p,msqc)
         
      elseif (kcase==kZHgaga) then
c        call singcheck(qqb_zh_gaga_g,qqb_zh_gaga_gs,p)
        if (includereal) call qqb_zh_gaga_g(p,msq)
        call qqb_zh_gaga_gs(p,msqc)     
      elseif (kcase==kZH__WW) then
c        call singcheck(qqb_zh_ww_g,qqb_zh_ww_gs,p)
        if (includereal) call qqb_zh_ww_g(p,msq)
        call qqb_zh_ww_gs(p,msqc)     
      elseif (kcase==kZH__ZZ) then
c        call singcheck(qqb_zh_zz_g,qqb_zh_zz_gs,p)
        if (includereal) call qqb_zh_zz_g(p,msq)
        call qqb_zh_zz_gs(p,msqc)     
      elseif (kcase==kZH1jet) then
c        call singcheck(qqb_ZH1jet_g,qqb_ZH1jet_gs,p)
c        if (includereal) call qqb_ZHgg_mad(p,msqa)
c         do j=-nf,nf
c         do k=-nf,nf
c         write(6,*) 'j,k,msqa',j,k,msq(j,k),msqa(j,k)
c         enddo
c         enddo
c         pause
        if (toponly) then
          msq=zip
          msqc=zip
        else
          if (includereal) call qqb_ZH1jet_g(p,msq)
          call qqb_ZH1jet_gs(p,msqc)
        endif
      elseif (kcase==ktwo_ew) then
c        if (includereal) call singcheck(qqb_twojet_mix_g,qqb_twojet_mix_gs,p) 
        if (  (-s(1,3) < cutoff) .or. (-s(2,3) < cutoff)
     &   .or. (-s(1,4) < cutoff) .or. (-s(2,4) < cutoff)
     &   .or. (-s(1,5) < cutoff) .or. (-s(2,5) < cutoff)
     &   .or. ( s(3,4) < cutoff) .or. ( s(3,5) < cutoff)
     &   .or. ( s(4,5) < cutoff)) goto 999
        if (includereal) call qqb_twojet_mix_g(p,msq)
        call qqb_twojet_mix_gs(p,msqc)
        wt_noew=zip ! non-EW calculation has no real contribution
      elseif (kcase==kdirgam) then
!        if (includereal) call singcheck(qqb_dirgam_g,qqb_dirgam_gs,p) 
        if (includereal) call qqb_dirgam_g(p,msq)
        call qqb_dirgam_gs(p,msqc)
      elseif (kcase==khflgam) then
c        if (includereal) call singcheck(qqb_hflgam_g,qqb_hflgam_gs,p)
        if (includereal) call qqb_hflgam_g(p,msq)
        call qqb_hflgam_gs(p,msqc)
      elseif (kcase==kgamgam) then
c        if (includereal) call singcheck(qqb_gamgam_g,qqb_gamgam_gs,p)
        if (includereal) call qqb_gamgam_g(p,msq)
        call qqb_gamgam_gs(p,msqc)
      elseif (kcase==kgg2gam) then
c        if (includereal) call singcheck(gg_2gam_g,gg_2gam_gs,p)
        if (includereal) call gg_2gam_g(p,msq)
        call gg_2gam_gs(p,msqc)
      elseif (kcase==kgmgmjt) then
c        if (includereal) call singcheck(qqb_gmgmjt_g,qqb_gmgmjt_gs,p)
        if (includereal) call qqb_gmgmjt_g(p,msq)
        call qqb_gmgmjt_gs(p,msqc)
      elseif (kcase==ktrigam) then
c        if (includereal) call singcheck(qqb_trigam_g,qqb_trigam_gs,p)
        if (includereal) call qqb_trigam_g(p,msq)
        call qqb_trigam_gs(p,msqc)
      elseif (kcase==kfourga) then 
c        if (includereal) call singcheck(qqb_fourgam_g,qqb_fourgam_gs,p)       
        if (includereal) call qqb_fourgam_g(p,msq)
        call qqb_fourgam_gs(p,msqc)
       elseif (kcase==kggfus0) then
c         call singcheck(gg_hg,gg_h_gs,p)       ! Checked 28/02/03
         if (includereal) call gg_hg(p,msq)
         call gg_h_gs(p,msqc)
       elseif (kcase==kHigaga) then
c         call singcheck(gg_hgamgamg,gg_hgamgam_gs,p)
         if (includereal) call gg_hgamgamg(p,msq)
         call gg_hgamgam_gs(p,msqc)
       elseif (kcase==kHi_Zga) then
c         call singcheck(gg_hzgamg,gg_hzgam_gs,p)
         if (includereal) call gg_hzgamg(p,msq)
         call gg_hzgam_gs(p,msqc)
      elseif ((kcase==kHWW_4l) .or. (kcase==kHWW2lq)) then
c        call singcheck(qqb_hww_g,qqb_hww_gs,p)
        if (includereal) call qqb_hww_g(p,msq)      
        call qqb_hww_gs(p,msqc)      
      elseif (kcase==kHWWdkW) then
        if (includereal) call dkqqb_hww_g(p,msq)      
        call dkqqb_hww_gs(p,msqc)      
      elseif (kcase==kHWWdkW) then
        if (includereal) call dkqqb_hww_g(p,msq)      
        call dkqqb_hww_gs(p,msqc)      
      elseif (kcase==kHZZ_4l) then
c        call singcheck(qqb_hzz_g,qqb_hzz_gs,p)
        if (includereal) call qqb_hzz_g(p,msq)      
        call qqb_hzz_gs(p,msqc)  
      elseif (kcase==kH_1jet) then
c        call singcheck(qqb_Hg_g,qqb_Hg_gs,p)       ! Checked 19/02/02
        if (includereal) call qqb_Hg_g(p,msq)  
        call qqb_Hg_gs(p,msqc) 
      elseif ((kcase==ktt_bbl) .or. (kcase==ktt_bbh)) then
c         call singcheck(qqb_QQbdk_g,qqb_QQbdk_gs,p) ! Checked 15/8/08
        if (includereal) call qqb_QQbdk_g(p,msq)  
        call qqb_QQbdk_gs(p,msqc) 
      elseif ((kcase==ktt_ldk) .or. (kcase==ktt_hdk)) then
        if     (decay1q2a == 1) then
c------ radiation in the decay of the top quark
        if (includereal) call dk1qqb_QQb_g(p,msq)
          call dk1qqb_QQb_gs(p,msqc)
      elseif (decay1q2a == 2) then 
c------ radiation in the decay of the anti-top quark
        if (includereal) call dk2qqb_QQb_g(p,msq)
          call dk2qqb_QQb_gs(p,msqc)
      else
        write(6,*) 'Problem generating PS, decay1q2a=',decay1q2a
        stop
      endif
      elseif (kcase==ktt_udk) then
        if     (decay1q2a == 1) then
c------ radiation in the decay of the top quark
        if (includereal) call dk1uqqb_QQb_g(p,msq)
          call dk1uqqb_QQb_gs(p,msqc)
      elseif (decay1q2a == 2) then 
c------ radiation in the decay of the anti-top quark
        if (includereal) call dk2uqqb_QQb_g(p,msq)
          call dk2uqqb_QQb_gs(p,msqc)
      else
        write(6,*) 'Problem generating PS, decay1q2a=',decay1q2a
        stop
      endif
      elseif (kcase==ktthWdk) then
        if     (decay1q2a == 1) then
c------ radiation in the decay of the top quark
        if (includereal) call dkW1qqb_QQb_g(p,msq)
          call dkW1qqb_QQb_gs(p,msqc)
      elseif (decay1q2a == 2) then 
c------ radiation in the decay of the anti-top quark
        if (includereal) call dkW2qqb_QQb_g(p,msq)
          call dkW2qqb_QQb_gs(p,msqc)
      else
        write(6,*) 'Problem generating PS, decay1q2a=',decay1q2a
        stop
      endif
      elseif (kcase==ktt_bbu) then
c         call singcheck(qqb_QQbdku_g,qqb_QQbdku_gs,p) !
c         pause
        if (includereal) call qqb_QQbdku_g(p,msq)  
        call qqb_QQbdku_gs(p,msqc) 
      elseif (kcase==ktt_mix) then
         if ((s(3,5) < cutoff) .or. (s(4,5) < cutoff)) goto 999
         if (abs(two*dot(p,3,4)+two*mt**2-zmass**2) .lt. 0.1_dp) then
           msq(:,:)=zip
           msqc(:,:,:)=zip
         else
c           call singcheck(qqb_QQb_mix_g,qqb_QQb_mix_gs,p)
           if (includereal) call qqb_QQb_mix_g(p,msq)
           call qqb_QQb_mix_gs(p,msqc)
           wt_noew=zip ! non-EW calculation has no real contribution
         endif
      elseif ((kcase==ktt_tot) .or. (kcase==kcc_tot)
     &   .or. (kcase==kbb_tot)) then
c        call singcheck(qqb_QQb_g,qqb_QQb_gs,p)
        if (includereal) call qqb_QQb_g(p,msq)
        call qqb_QQb_gs(p,msqc)
      elseif (kcase==ktopanom) then
        bbfrac = 0._dp

        msq_light = 0._dp
        msq_heavy = 0._dp
        msqc_light = 0._dp
        msqc_heavy = 0._dp
        ! this is required even when use_DDIS is .false.
        call singletop2_scale_setup(p)
        if (includereal) then
          ! setup once for real emission matrix elements
          call singletop2_real(p,msq_light,.true.,.false.)
          call singletop2_real(p,msq_heavy,.false.,.true.)
          msq = msq_light + msq_heavy
        endif
        ! scale setup for use_DDIS case is called within _gs and dipolesub
        call singletop2_gs(p,msqc_light,.true.,.false.)
        call singletop2_gs(p,msqc_heavy,.false.,.true.)
        msqc = msqc_light + msqc_heavy

        ! reset scale setup before singcheck, which calls the emission routine first
c       call singletop2_scale_setup(p)
c       call singcheck(singletop2_real_heavy, singletop2_gs_heavy, p)
      elseif (kcase==kbq_tpq) then
        bbfrac = 0._dp
c        call singcheck(qqb_tbb_g,qqb_tbb_gs,p)
        if (includereal) call qqb_tbb_g(p,msq)
        call qqb_tbb_gs(p,msqc)

      elseif (kcase==kt_bbar) then
c        call singcheck(qqb_tbb_g,qqb_tbb_gs,p)
        if (includereal) call qqb_tbbdk_g(p,msq)
        call qqb_tbbdk_gs(p,msqc)
      elseif (kcase==kttdkay) then
c        call singcheck(qqb_tbb_g_new,qqb_tbb_gs,p)
        if (includereal) call bq_tpq_gdk(p,msq)
        call bq_tpq_gsdk(p,msqc)
      elseif (kcase==ktdecay) then
c        call singcheck(qqb_tbb_g_new,qqb_tbb_gs,p)
      if (includereal) call dkqqb_tbbdk_g(p,msq)
        call dkqqb_tbbdk_gs(p,msqc)       
       elseif (kcase==kW_tndk) then
c        call singcheck(qqb_w_tndk_g,qqb_w_tndk_gs,p)      ! Checked 12/3/04
        if (includereal) call qqb_w_tndk_g(p,msq)
        call qqb_w_tndk_gs(p,msqc)
      elseif (kcase==kW_twdk) then
c        call singcheck(qqb_w_twdk_g,qqb_w_twdk_gs,p)            ! Checked 2/4/05
        if (includereal) call qqb_w_twdk_g(p,msq)
        call qqb_w_twdk_gs(p,msqc)
      elseif (kcase==kWtdkay) then
        if (includereal) call dkqqb_w_twdk_g(p,msq)
        call dkqqb_w_twdk_gs(p,msqc)
      elseif ( (kcase==kqq_ttw)) then 
        if (includereal) call qqb_ttw_g(p,msq)
        call qqb_ttw_gs(p,msqc)
      elseif ( (kcase==kttwldk)) then
        if     (decay1q2a == 1) then
c------ radiation in the decay of the top quark
        if (includereal) call dk1qqb_ttw_g(p,msq)
          call dk1qqb_ttw_gs(p,msqc)
      elseif (decay1q2a == 2) then 
c------ radiation in the decay of the anti-top quark
        if (includereal) call dk2qqb_ttw_g(p,msq)
          call dk2qqb_ttw_gs(p,msqc)
      else
        write(6,*) 'Problem generating PS, decay1q2a=',decay1q2a
        stop
      endif
      elseif (kcase==kggfus1) then
!        if (includereal) call singcheck(gg_hgg,gg_hg_gs,p)
        if (includereal) call gg_hgg(p,msq)
        call gg_hg_gs(p,msqc)
      elseif (kcase==kHi_Zaj) then
!        if (includereal) call singcheck(gg_hgg_zgam,gg_hg_zgam_gs,p)
        if (includereal) call gg_hgg_zgam(p,msq)
        call gg_hg_zgam_gs(p,msqc)
      elseif (kcase==khjetma) then
        badpoint = .false.

        if (includereal) call hjetmass_r(p,msq)
        if (badpoint .eqv. .false.) then
          call hjetmass_gs(p,msqc)
        else
            msq = 0._dp
            msqc = 0._dp
        endif
      elseif (kcase==kHgagaj) then
c        call singcheck(gg_hgagagg,gg_hgagag_gs,p)
        if (includereal) call gg_hgagagg(p,msq)
        call gg_hgagag_gs(p,msqc)
      elseif (kcase==kHWWjet) then
c        call singcheck(gg_hWWgg,gg_hWWg_gs,p)
        if (includereal) call gg_hWWgg(p,msq)
        call gg_hWWg_gs(p,msqc)
      elseif (kcase==kHWW2jt) then
c        call singcheck(gg_hWWgg,gg_hWWg_gs,p)
        if (includereal) call gg_hWWggg(p,msq)
        call gg_hWWgg_gs(p,msqc)
      elseif (kcase==kHZZjet) then
c        call singcheck(gg_hZZgg,gg_hZZg_gs,p)
        if (includereal) call gg_hZZgg(p,msq)
        call gg_hZZg_gs(p,msqc)
      elseif (kcase==kHZZ2jt) then
c        call singcheck(gg_hZZgg,gg_hZZg_gs,p)
        if (includereal) call gg_hZZggg(p,msq)
        call gg_hZZgg_gs(p,msqc)
c      write(6,*) msq
c      pause
      elseif (kcase==kqq_Hqq) then
c        call singcheck(VV_Hqq_g,VV_Hqq_gs,p)
        if (includereal) call VV_Hqq_g(p,msq)
        call VV_Hqq_gs(p,msqc)
      elseif (kcase==kqq_Hgg) then
c        call singcheck(VV_Hgaga_g,VV_Hgaga_gs,p)
        if (includereal) call VV_Hgaga_g(p,msq)
        call VV_Hgaga_gs(p,msqc)
      elseif (kcase==kqq_HWW) then
c        call singcheck(VV_HWW_g,VV_HWW_gs,p)
       if (includereal) call VV_HWW_g(p,msq)
        call VV_HWW_gs(p,msqc)
      elseif (kcase==kqq_HZZ) then
c        call singcheck(VV_HZZ_g,VV_HZZ_gs,p)
        if (includereal) call VV_HZZ_g(p,msq)
        call VV_HZZ_gs(p,msqc)
      elseif (kcase==kggfus2) then
c        call singcheck(gg_hggg,gg_hgg_gs,p)            ! Checked 10/29/09
        if (includereal) call gg_hggg(p,msq)
        call gg_hgg_gs(p,msqc)
      elseif (kcase==kgagajj) then
c        call singcheck(gg_hggg,gg_hgg_gs,p)
        if (includereal) call gg_hggg(p,msq)
        call gg_hgg_gs(p,msqc)
      elseif (kcase==kqg_tbq) then
c        call singcheck(qg_tbq_g,qg_tbq_gs,p)       ! Checked 2/4/08
       if (includereal) call qg_tbq_g(p,msq)
       call qg_tbq_gs(p,msqc)
      elseif (kcase==k4ftwdk) then
c        call singcheck(qg_tbq_g,qg_tbq_gs,p)
       if (includereal) call qg_tbqdk_g(p,msq)
       call qg_tbqdk_gs(p,msqc)
      elseif (kcase==kdk_4ft) then
c       if (includereal) call dkqg_tbqdk_g_old(p,msq)  
       if (includereal) call dkqg_tbqdk_g(p,msq)  
       call dkqg_tbqdk_gs(p,msqc)
      elseif (kcase==kqq_tbg) then
c        call singcheck(qq_tbg_g,qq_tbg_gs,p)       ! Checked 8/9/08
       if (includereal) call qq_tbg_g(p,msq)
       call qq_tbg_gs(p,msqc)
      elseif (kcase==kepem3j) then
c        call singcheck(epem3j_g,epem3j_gs,p)       ! Checked 17/11/08
       if (includereal) call epem3j_g(p,msq)
       call epem3j_gs(p,msqc)
      elseif (kcase==kgQ__ZQ) then
c        call singcheck(gQ_zQ_g,gQ_zQ_gs,p)
        if (includereal) call gQ_zQ_g(p,msq)
        call gQ_zQ_gs(p,msqc)
      elseif (kcase==kZ_bjet) then
c        call singcheck(qqb_zbjet_g,qqb_zbjet_gs,p)      ! Checked 07/18/05
        if (includereal) call qqb_zbjet_g(p,msq)
        call qqb_zbjet_gs(p,msqc)
      elseif (kcase==kW_bjet) then
c        call singcheck(qqb_wbjet_g,qqb_wbjet_gs,p) ! Rechecked 14/3/08
        if (includereal) call qqb_wbjet_g(p,msq)
        call qqb_wbjet_gs(p,msqc)
      elseif (kcase==kWcsbar) then
c        call singcheck(qqb_w_g,qqb_w_gs,p)         ! Checked 11/30/01
        if (includereal) call qqb_w_g(p,msq)      
        call qqb_w_gs(p,msqc)     
      elseif (kcase==kWcs_ms) then
        if (includereal) call qqb_w_cjet(p,msq)
        ndmax=0
      elseif (kcase==kZ_tjet) then
c        if (includereal) call singcheck(qq_tchan_ztqg,qq_tchan_ztq_gs,p)
        if (includereal) call qq_tchan_ztqg(p,msq)
        call qq_tchan_ztq_gs(p,msqc)
      elseif (kcase==kZ_tdkj) then
c        if (includereal) call singcheck(qq_tchan_ztqg_dk,
c     &                                  qq_tchan_ztq_dk_gs,p)
        if (includereal) call qq_tchan_ztqg_dk(p,msq)
        call qq_tchan_ztq_dk_gs(p,msqc)
      elseif (kcase==kH_tjet) then
c        if (includereal) call singcheck(qq_tchan_htqg,qq_tchan_htq_gs,p)
        if (includereal) call qq_tchan_htqg(p,msq)
        call qq_tchan_htq_gs(p,msqc)
      elseif (kcase==kH_tdkj) then
        if (includereal) call qq_tchan_htqg_dk(p,msq)
        call qq_tchan_htq_dk_gs(p,msqc)
      elseif (kcase==ktottth) then
c        if (includereal) call qqb_tottth_g(p,msq)
c        call qqb_tottth_gs(p,msqsc)
c        call singcheck(qqb_tottth_g,qqb_tottth_gs,p)
c        call compare_madgraph(p,qqb_tottth_g,qqb_tottth_g_mad)
      elseif (kcase==kdm_jet) then 
         if(includereal) call qqb_dm_monojet_g(p,msq)       
         call qqb_dm_monojet_gs(p,msqc)
      elseif (kcase==kdm_gam) then
!         if(includereal) call singcheck(qqb_dm_monophot_g
!     &        ,qqb_dm_monophot_gs,p)
         if(includereal) call qqb_dm_monophot_g(p,msq) 
         call qqb_dm_monophot_gs(p,msqc)

      endif

! code to find power of alpha-s for scale variation
      if ((doScalevar .and. currentPDF == 0 .and. bin) .and. (foundpow .eqv. .false.)) then
        if (itrial == 1) then
          msqtrial=maxval(msq)
          if (msqtrial == 0) goto 999
          if (dynamicscale) call scaleset(initscale,initfacscale,p)
          as=as*two
          ason2pi=ason2pi*two
          ason4pi=ason4pi*two
          gsq=gsq*two
          itrial=itrial+1
          goto 66
        endif
        msqtrial=maxval(msq)/msqtrial
        alphaspow=-1
        if (abs(msqtrial-one) < 1.e-8) alphaspow=0
        if (abs(msqtrial-two) < 1.e-8) alphaspow=1
        if (abs(msqtrial-four) < 1.e-8) alphaspow=2
        if (abs(msqtrial-eight) < 1.e-8) alphaspow=3
        if (abs(msqtrial-16._dp) < 1.e-8) alphaspow=4
        if (abs(msqtrial-32._dp) < 1.e-8) alphaspow=5
        if (alphaspow == -1) then
          write(6,*) 'Unable to determine power of alpha-s for scale variation'
          write(6,*) ' msqtrial = ',msqtrial
          stop
        endif
        as=as/two
        ason2pi=ason2pi/two
        ason4pi=ason4pi/two
        gsq=gsq/two
        foundpow=.true.
!        write(6,*) 'Found alpha-s power: ',alphaspow
        goto 66
      endif

      xmsq = 0._dp
      bbfrac = 0._dp

      flux=fbGeV2/(two*xx1*xx2*W)

c--- initialize a PDF set here, if calculating errors
  777 continue

      if (maxPDFsets > 0) then
          pdfreweight(:) = 0._dp
          if (currentPDF == 0) then
              ! saves PDF central value xmsq
              centralval(:) = 0._dp
          endif
      endif

      xmsq = 0._dp
      bbfrac = 0._dp

c--- calculate PDF's

      if ((maxPDFsets > 0) .and. bin) then
! PDF errors and dynamic scale
          if (dynamicscale) then
            do nd=ndmax,0,-1  ! so that fx1,fx2 correct for real kinematics
              !--- in case dipole is not used, set up dummy value of scale for safety
              !--- and set all PDF entries to zero
              if (dipscale(nd) < 1.e-8_dp) then        
                dipscale(nd)=dipscale(0)
                fx1(:)=0._dp
                fx2(:)=0._dp
              else
                call fdist(ih1,xx1,dipscale(nd),fx1)
                call fdist(ih2,xx2,dipscale(nd),fx2)
                dipfx1(nd,:)=fx1(:)
                dipfx2(nd,:)=fx2(:)
              endif
            enddo
            if ((kcase==kqg_tbq) .or. (kcase==k4ftwdk)) then
              fx1_H(:)=fx1(:)
              fx1_L(:)=fx1(:)
              fx2_H(:)=fx2(:)
              fx2_L(:)=fx2(:)
            endif
! PDF errors and fixed scale
          else ! no dynamic scale
            !--- for single top + b, make sure to use two different scales
            if ((kcase==kqg_tbq) .or. (kcase==k4ftwdk)) then
              call fdist(ih1,xx1,facscale_H,fx1_H)
              call fdist(ih2,xx2,facscale_H,fx2_H)
              call fdist(ih1,xx1,facscale_L,fx1_L)
              call fdist(ih2,xx2,facscale_L,fx2_L)
              do j=-nf,nf
                if (j == 0) then  ! heavy quark line has gluon init. state
                  fx1(j)=fx1_H(j)
                  fx2(j)=fx2_H(j)
                else
                  fx1(j)=fx1_L(j)
                  fx2(j)=fx2_L(j)
                endif
              enddo
            else
              call fdist(ih1,xx1,facscale,fx1)
              call fdist(ih2,xx2,facscale,fx2)
            endif
          endif
      else ! case of no PDFerrors
! NO PDF errors and dynamic scale
        if (dynamicscale) then
            do nd=ndmax,0,-1  ! so that fx1,fx2 correct for real kinematics
              !--- in case dipole is not used, set up dummy value of scale for safety
              !--- and set all PDF entries to zero
              if (dipscale(nd) < 1.e-8_dp) then        
                dipscale(nd)=dipscale(0)
                fx1(:)=0._dp
                fx2(:)=0._dp
              else
                if (doscalevar .and. currentPDF == 0 .and. bin) then
                  call fdist(ih1,xx1,dipscale(nd)*two,fx1)
                  call fdist(ih2,xx2,dipscale(nd)*two,fx2)
                  dipfx1up(nd,:)=fx1(:)
                  dipfx2up(nd,:)=fx2(:)
                  call fdist(ih1,xx1,dipscale(nd)/two,fx1)
                  call fdist(ih2,xx2,dipscale(nd)/two,fx2)
                  dipfx1dn(nd,:)=fx1(:)
                  dipfx2dn(nd,:)=fx2(:)
                  xmsqvar(:,nd)=zip
                endif
                if (kcase == ktopanom .and. use_DDIS) then
                  ! first part: pdf's for dipole contributions
                  if (nd > 0) then
                      if (singletop2_dipscale(nd,1) /= 0._dp .and. singletop2_dipscale(nd,2) /= 0._dp) then
                          call fdist(ih1,xx1,singletop2_dipscale(nd,1),singletop2_dipole_pdfs(nd,1,:))
                          call fdist(ih2,xx2,singletop2_dipscale(nd,2),singletop2_dipole_pdfs(nd,2,:))
                      else
                          if (sum(abs(msqc(nd,:,:))) > 0._dp) then
                              write (*,*) "CRITICAL WARNING: dipole contribution ", sum(abs(msqc(nd,:,:))), " but scale zero"
                          endif
                      endif
                  endif
                else ! usual case
                  call fdist(ih1,xx1,dipscale(nd),fx1)
                  call fdist(ih2,xx2,dipscale(nd),fx2)
                  dipfx1(nd,:)=fx1(:)
                  dipfx2(nd,:)=fx2(:)
                endif
              endif
            enddo
            if (kcase == ktopanom .and. use_DDIS) then
              ! scales were previously set by some dipole routine, reset with unmodified momenta
              call singletop2_scale_setup(p)
              ! in addition we need all the PDF sets for the real emission matrix elements
              call fdist(ih1,xx1,facscale_beam1_islight_onlight, singletop2_pdfs(i_beam1_islight_onlight,:))
              call fdist(ih2,xx2,facscale_beam2_isheavy_onlight, singletop2_pdfs(i_beam2_isheavy_onlight,:))

              call fdist(ih1,xx1,facscale_beam1_isheavy_onlight, singletop2_pdfs(i_beam1_isheavy_onlight,:))
              call fdist(ih2,xx2,facscale_beam2_islight_onlight, singletop2_pdfs(i_beam2_islight_onlight,:))

              call fdist(ih1,xx1,facscale_beam1_islight_onheavy, singletop2_pdfs(i_beam1_islight_onheavy,:))
              call fdist(ih2,xx2,facscale_beam2_isheavy_onheavy, singletop2_pdfs(i_beam2_isheavy_onheavy,:))

              call fdist(ih1,xx1,facscale_beam1_isheavy_onheavy, singletop2_pdfs(i_beam1_isheavy_onheavy,:))
              call fdist(ih2,xx2,facscale_beam2_islight_onheavy, singletop2_pdfs(i_beam2_islight_onheavy,:))
            endif
            if ((kcase==kqg_tbq) .or. (kcase==k4ftwdk)) then
              fx1_H(:)=fx1(:)
              fx1_L(:)=fx1(:)
              fx2_H(:)=fx2(:)
              fx2_L(:)=fx2(:)
            endif
! NO PDF errors and NO dynamic scale
        else ! no dynamic scale
              !--- for single top + b, make sure to use two different scales
              if ((kcase==kqg_tbq) .or. (kcase==k4ftwdk)) then
                  call fdist(ih1,xx1,facscale_H,fx1_H)
                  call fdist(ih2,xx2,facscale_H,fx2_H)
                  call fdist(ih1,xx1,facscale_L,fx1_L)
                  call fdist(ih2,xx2,facscale_L,fx2_L)
                  do j=-nf,nf
                    if (j == 0) then  ! heavy quark line has gluon init. state
                      fx1(j)=fx1_H(j)
                      fx2(j)=fx2_H(j)
                    else
                      fx1(j)=fx1_L(j)
                      fx2(j)=fx2_L(j)
                    endif
                  enddo
              else
                  call fdist(ih1,xx1,facscale,fx1)
                  call fdist(ih2,xx2,facscale,fx2)
                  if (doScalevar .and. currentPDF == 0 .and. bin) then
                    call fdist(ih1,xx1,facscale*two,fx1up)
                    call fdist(ih2,xx2,facscale*two,fx2up)
                    call fdist(ih1,xx1,facscale/two,fx1dn)
                    call fdist(ih2,xx2,facscale/two,fx2dn)
                    xmsqvar(:,:)=zip
                  endif
              endif
        endif

      endif

      do j=-nflav,nflav
      do k=-nflav,nflav

c     ! only gluon quark, this selects only the bbfrac piece
c     if (.not. ((j==0 .and. abs(k) /= 5) .or. (abs(j) /= 5 .and. k==0))) cycle

      if (ggonly) then
      if ((j.ne.0) .or. (k.ne.0)) cycle
      endif

      if (gqonly) then
      if (((j==0).and.(k==0)) .or. ((j.ne.0).and.(k.ne.0))) cycle
      endif      
      
      if (noglue) then 
      if ((j==0) .or. (k==0)) cycle
      endif

      if (omitgg) then 
      if ((j==0) .and. (k==0)) cycle
      endif

      if ((kcase==kWcsbar).and.(j .ne. 4).and.(k .ne. 4)) cycle

      if     (j > 0) then
        sgnj=+1
      elseif (j < 0) then
        sgnj=-1
      else
        sgnj=0
      endif
      if     (k > 0) then
        sgnk=+1
      elseif (k < 0) then
        sgnk=-1
      else
        sgnk=0
      endif

c--- for single top + b, make sure to use two different scales
      if ((kcase==kqg_tbq) .or. (kcase==k4ftwdk)) then
        xmsqjk=fx1_L(j)*fx2_H(k)*msqLH(j,k)
     &        +fx1_H(j)*fx2_L(k)*msqHL(j,k)
      else
        xmsqjk = 0._dp

        if (kcase == ktopanom .and. use_DDIS) then
          if (abs(k) == 5) then ! qb, gb
            xmsqjk = xmsqjk + singletop2_pdfs(i_beam1_islight_onlight,j) *
     &              singletop2_pdfs(i_beam2_isheavy_onlight,k) * msq_light(j,k)
            xmsqjk = xmsqjk + singletop2_pdfs(i_beam1_islight_onheavy,j) *
     &              singletop2_pdfs(i_beam2_isheavy_onheavy,k) * msq_heavy(j,k)
          elseif (abs(j) == 5) then ! bq, bg
            xmsqjk = xmsqjk + singletop2_pdfs(i_beam1_isheavy_onlight,j) *
     &              singletop2_pdfs(i_beam2_islight_onlight,k) * msq_light(j,k)
            xmsqjk = xmsqjk + singletop2_pdfs(i_beam1_isheavy_onheavy,j) *
     &              singletop2_pdfs(i_beam2_islight_onheavy,k) * msq_heavy(j,k)
          elseif (abs(j) /= 5 .and. k == 0) then ! this is qg
            xmsqjk = xmsqjk + singletop2_pdfs(i_beam1_islight_onheavy,j) *
     &              singletop2_pdfs(i_beam2_isheavy_onheavy,k) * msq_heavy(j,k)
          elseif (j == 0 .and. abs(k) /= 5) then ! this is crossed qg
            xmsqjk = xmsqjk + singletop2_pdfs(i_beam1_isheavy_onheavy,j) *
     &              singletop2_pdfs(i_beam2_islight_onheavy,k) * msq_heavy(j,k)
          else ! no other cases for our process
            if (msq(j,k) /= 0._dp) then
              write (*,*) "CRITICAL WARNING: unexpected channel contribution", j,k, msq(j,k)
            endif
            xmsqjk = 0._dp
          endif
        else ! usual case
          xmsqjk=fx1(j)*fx2(k)*msq(j,k)
        endif
      
      endif
      
      xmsq(0)=xmsq(0)+xmsqjk

      if (doscalevar .and. currentPDF == 0 .and. bin) then
        if (dynamicscale) then
          xmsqvar(1,0)=xmsqvar(1,0)+dipfx1up(0,j)*dipfx2up(0,k)*msq(j,k)
          xmsqvar(2,0)=xmsqvar(2,0)+dipfx1dn(0,j)*dipfx2dn(0,k)*msq(j,k)
        else
          xmsqvar(1,0)=xmsqvar(1,0)+fx1up(j)*fx2up(k)*msq(j,k)
          xmsqvar(2,0)=xmsqvar(2,0)+fx1dn(j)*fx2dn(k)*msq(j,k)
        endif
      endif

      ! ktopanom: two b's, real emission contribution
      if (nproc==169 .or. nproc==164 .or. nproc==161 .or. nproc==166) then
        if ((j == 0 .and. abs(k) /= 5) .or. (abs(j) /= 5 .and. k ==0)) then
            bbfrac(0) = bbfrac(0) + xmsqjk
        endif
      endif

      ! add dipole contributions
      do nd=1,ndmax

        if (dynamicscale) then
            if (kcase == ktopanom .and. use_DDIS) then
                xmsqjk = singletop2_dipole_pdfs(nd,1,j) * singletop2_dipole_pdfs(nd,2,k) * (-msqc(nd,j,k))
c               write (*,*) j,k,singletop2_dipole_pdfs(nd,1,j),singletop2_dipole_pdfs(nd,2,k),(-msqc(nd,j,k))
            else
              xmsqjk=dipfx1(nd,j)*dipfx2(nd,k)*(-msqc(nd,j,k))
              if (doscalevar .and. currentPDF == 0 .and. bin) then
                xmsqvar(1,nd)=xmsqvar(1,nd)+dipfx1up(nd,j)*dipfx2up(nd,k)*(-msqc(nd,j,k))
                xmsqvar(2,nd)=xmsqvar(2,nd)+dipfx1dn(nd,j)*dipfx2dn(nd,k)*(-msqc(nd,j,k))
              endif
            endif
        else ! no dynamicscale
            xmsqjk=fx1(j)*fx2(k)*(-msqc(nd,j,k))
            if (doscalevar .and. currentPDF == 0 .and. bin) then
              xmsqvar(1,nd)=xmsqvar(1,nd)+fx1up(j)*fx2up(k)*(-msqc(nd,j,k))
              xmsqvar(2,nd)=xmsqvar(2,nd)+fx1dn(j)*fx2dn(k)*(-msqc(nd,j,k))
            endif
        endif ! dynamicscale

      ! ktopanom: two b's, dipole contribution
        if (nproc==169 .or. nproc==164 .or. nproc==161 .or. nproc==166) then
c         if (nd == 13 .or. nd == 14) the
            if ((j == 0 .and. abs(k) /= 5) .or. (abs(j) /= 5 .and. k ==0)) then
                bbfrac(nd) = bbfrac(nd) + xmsqjk
            endif
c         endif
        endif

        xmsq(nd)=xmsq(nd)+xmsqjk

      enddo

      enddo
      enddo

      ! for PDF errors we want to save the central result as normalization
      ! for pdfreweight
      if (currentPDF == 0) then
          realint = 0._dp
      endif

      ! summation and histogramming for all dipole contributions
      do nd=0,ndmax
          if (kcase == ktopanom .or. kcase == kbq_tpq) then
            if (xmsq(nd) /= 0._dp) then
                bbfrac(nd) = bbfrac(nd)/xmsq(nd)
            endif
          endif

          xmsq(nd)=xmsq(nd)*flux*pswt/BrnRat

          if (xmsq(nd) /= 0._dp) then
              q(:,:) = ptilde(nd,:,:) 
              ! check cuts
              if (incldip(nd)) incldip(nd)=includedipole(nd,q)
              if (incldip(nd) .eqv. .false.) then
                  call dotem(nvec,p,s)
                  xmsq(nd) = 0._dp
              endif
          endif

          if (ieee_is_nan(xmsq(nd))) then
            if (debug) write(6,*) 'discarding point with weight NaN: pswt=',pswt
            if (bin) call threadStorageOp(shtmpreset)
            goto 999
          endif 

          if (.not. ieee_is_finite(xmsq(nd))) then
            if (debug) write(6,*) 'discarding point with infinite weight: pswt=',pswt
            if (bin) call threadStorageOp(shtmpreset)
            goto 999
          endif

          ! histogramming
          if (xmsq(nd) /= 0._dp) then
              if (currentPDF == 0) then
                  if (includeTaucutgrid(nd)) then
                      realint = realint + xmsq(nd)
                  endif
              endif

              if (bin) then
                  if (maxPDFsets > 0) then
                      if (currentPDF == 0) then
                          centralval(:) = xmsq(:)
                      endif

                      if (currentPDF > 0) then
                          pdfreweight(currentPDF) = (centralval(nd) - xmsq(nd))*wgt
                      endif

                  endif

                  if (doscalevar .and. currentPDF == 0 .and. bin) then
                      scalereweight(:) = 0._dp
                      if (dynamicscale) then
                        asorig=alphas(dipscale(nd)*initscale/initfacscale,amz,nlooprun)
                        scaleup=dipscale(nd)*initscale/initfacscale*two
                        scaledn=dipscale(nd)*initscale/initfacscale/two
                      else
                        asorig=as
                        scaleup=scale*two
                        scaledn=scale/two
                      endif
                      scalereweight(1)=(alphas(scaleup,amz,nlooprun)/asorig)**alphaspow
                      scalereweight(2)=(alphas(scaledn,amz,nlooprun)/asorig)**alphaspow
                      scalereweight(1)=scalereweight(1)*xmsqvar(1,nd)*flux*pswt/BrnRat/xmsq(nd)
                      scalereweight(2)=scalereweight(2)*xmsqvar(2,nd)*flux*pswt/BrnRat/xmsq(nd)
                      if (maxscalevar == 6) then
                        scalereweight(3)=scalereweight(1)*xmsq(nd)/xmsqvar(1,nd)/flux/pswt*BrnRat
                        scalereweight(4)=scalereweight(2)*xmsq(nd)/xmsqvar(2,nd)/flux/pswt*BrnRat
                        scalereweight(5)=xmsqvar(1,nd)*flux*pswt/BrnRat/xmsq(nd)
                        scalereweight(6)=xmsqvar(2,nd)*flux*pswt/BrnRat/xmsq(nd)
                      endif
                  endif

                  call getptildejet(nd,pjet)
                  call dotem(nvec,pjet,s)
                  call nplotter(pjet,xmsq(nd)*wgt,(xmsq(nd)*wgt)**2,nd)
              endif
          endif
      enddo ! loop over dipole contributions

      if (bin) then
         call threadStorageOp(shtmpcommit)

#ifdef HAVE_LHAPDF
         if (maxPDFsets > 0) then
            currentPDF=currentPDF+1
            if (currentPDF <= maxPDFsets) then
                if (doPDFAlphas) then
                    ! only if as(mz) has changed recompute matrix element
                    if ( abs(getalphas(zmass) - amz) > 1d-5 ) then
                        goto 771
                    endif
                endif
                goto 777
            endif
         endif
#endif
      endif

      if (abs(realint*wgt) > wtmax) then
        wtmax=abs(realint*wgt)
      endif

c--- handle special case of Qflag and Gflag
      if (QandGflag) then
        QandGint=QandGint+realint
        if ((Gflag) .and. (.not.(Qflag))) then
c--- go back for second pass (Qflag), calling includedipole first to reset
c--- the values of "jets" and "jetlabel"
          includereal=includedipole(0,p)
          Qflag=.true.
          Gflag=.false.
          goto 44
        else
c--- return both to .true. and assign value to realint (to return to VEGAS)
          Qflag=.true.
          Gflag=.true.
          realint=QandGint
        endif
      endif

      if (enable_reweight_user .and. includereal) then
          realint = realint * reweight_user(p)
      endif

      return

      ! any errors where the whole point is set to zero
 999  continue

      realint = 0._dp

      if (QandGflag) then
        Qflag=.true.
        Gflag=.true.
      endif

      return
      end
