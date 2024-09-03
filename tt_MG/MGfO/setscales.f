      subroutine set_alphaS(xp)
      implicit none
      include "genps.inc"
      include "nexternal.inc"
      include "run.inc"
      include "coupl.inc"
      include "timing_variables.inc"
      
      double precision xp(0:3,nexternal)
      double precision dummy,dummyQES,dummies(2)
      integer i,j

      character*80 muR_id_str,muF1_id_str,muF2_id_str,QES_id_str
      common/cscales_id_string/muR_id_str,muF1_id_str,
     #                         muF2_id_str,QES_id_str
      double precision PP(0:3,max_particles)
      COMMON /MOMENTA_PP/PP

      logical firsttime
      data firsttime/.true./

      call cpu_time(tBefore)

      if (firsttime) then
        firsttime=.false.
        call set_ren_scale(xp,dummy)
        if(dummy.lt.0.2d0)then
          write(*,*)'Error in set_alphaS: muR too soft',dummy
          stop
        endif
        call set_fac_scale(xp,dummies)
        if(dummies(1).lt.0.2d0.or.dummies(2).lt.0.2d0)then
          write(*,*)'Error in set_alphaS: muF too soft',
     #              dummies(1),dummies(2)
          stop
        endif
        call set_QES_scale(xp,dummyQES)
        if(scale.lt.0.2d0)then
          write(*,*)'Error in set_alphaS: QES too soft',dummyQES
          stop
        endif
        write(*,*)'Scale values (may change event by event):'
        write(*,200)'muR,  muR_reference: ',dummy,
     #              dummy/muR_over_ref,muR_over_ref
        write(*,200)'muF1, muF1_reference:',dummies(1),
     #              dummies(1)/muF1_over_ref,muF1_over_ref
        write(*,200)'muF2, muF2_reference:',dummies(2),
     #              dummies(2)/muF2_over_ref,muF2_over_ref
        write(*,200)'QES,  QES_reference: ',dummyQES,
     #              dummyQES/QES_over_ref,QES_over_ref
        write(*,*)' '
        write(*,*)'muR_reference [functional form]:'
        write(*,*)'   ',muR_id_str(1:len_trim(muR_id_str))
        write(*,*)'muF1_reference [functional form]:'
        write(*,*)'   ',muF1_id_str(1:len_trim(muF1_id_str))
        write(*,*)'muF2_reference [functional form]:'
        write(*,*)'   ',muF2_id_str(1:len_trim(muF2_id_str))
        write(*,*)'QES_reference [functional form]: '
        write(*,*)'   ',QES_id_str(1:len_trim(QES_id_str))
        write(*,*)' '
        write(*,*) 'alpha_s=',g**2/(16d0*atan(1d0))
        do i=0,3
          do j=1,max_particles
            pp(i,j) = 0d0
          enddo
        enddo
      endif
      call set_QES_scale(xp,dummyQES)
      call set_fac_scale(xp,dummies)
      call set_ren_scale(xp,dummy)
      if ( .not.fixed_ren_scale.or.
     &         .not.fixed_couplings.or.
     &             .not.fixed_QES_scale) then
        if (.not.fixed_couplings)then
          do i=0,3
            do j=1,nexternal
              PP(i,j)=xp(i,j)
            enddo
          enddo
        endif
      endif

      call cpu_time(tAfter)
      t_coupl=t_coupl+(tAfter-tBefore)
      
 200  format(1x,a,2(1x,d12.6),2x,f4.2)

      return
      end

      subroutine set_ren_scale(pp,muR)
      implicit none
      include 'genps.inc'
      include 'nexternal.inc'
      include 'run.inc'
      include 'coupl.inc'
      double precision pp(0:3,nexternal),muR
      double precision mur_temp,mur_ref_dynamic,alphas
      double precision pi
      parameter (pi=3.14159265358979323846d0)
      character*80 muR_id_str,muF1_id_str,muF2_id_str,QES_id_str
      common/cscales_id_string/muR_id_str,muF1_id_str,
     #                         muF2_id_str,QES_id_str
      character*80 temp_scale_id
      common/ctemp_scale_id/temp_scale_id
      double precision minscaleR
      parameter (minscaleR=2d0)
      logical calculatedBorn
      common/ccalculatedBorn/calculatedBorn
      temp_scale_id='  '
      if(fixed_ren_scale)then
        mur_temp=muR_ref_fixed
        temp_scale_id='fixed'
      else
        mur_temp=max(minscaleR,muR_ref_dynamic(pp))
      endif
      muR=muR_over_ref*mur_temp
      muR2_current=muR**2
      muR_id_str=temp_scale_id
      mu_r = muR
      scale=muR
      g=sqrt(4d0*pi*alphas(scale))
      call update_as_param()
      calculatedBorn=.false.
      return
      end


      function muR_ref_dynamic(pp)
      use extra_weights
      implicit none
      include 'genps.inc'
      include 'nexternal.inc'
      include 'run.inc'
      include 'cuts.inc'
      include 'orders.inc'
      double precision muR_ref_dynamic,pp(0:3,nexternal)
      double precision tmp,scale_global_reference,pt,et,dot,sumdot
      external pt,et,dot,sumdot
      character*80 temp_scale_id
      common/ctemp_scale_id/temp_scale_id
      integer i,imurtype
      parameter (imurtype=1)
      integer j
      LOGICAL  IS_A_J(NEXTERNAL),IS_A_LP(NEXTERNAL),IS_A_LM(NEXTERNAL)
      LOGICAL  IS_A_PH(NEXTERNAL)
      COMMON /TO_SPECISA/IS_A_J,IS_A_LP,IS_A_LM,IS_A_PH
      integer NN,NJET,JET(nexternal),iqcd
      double precision pQCD(0:3,nexternal),PJET(0:3,nexternal)
      double precision rfj,sycut,palg,amcatnlo_fastjetdmergemax
     &     ,tmp1,tmp2,xm2
      integer nFxFx_ren_scales
      double precision FxFx_ren_scales(0:nexternal),FxFx_fac_scale(2)
      common/c_FxFx_scales/FxFx_ren_scales,nFxFx_ren_scales
     $     ,FxFx_fac_scale
      integer bpower
      tmp=0
      if (nincoming.eq.1) then
         tmp=pp(0,1) ! mass of the decaying particle
         temp_scale_id='Mass of decaying particle'
      elseif(ickkw.eq.3)then
         bpower=born_orders(qcd_pos)/2
         if (bpower.gt.nFxFx_ren_scales) then
            tmp=FxFx_ren_scales(0)**
     &           (bpower-(nFxFx_ren_scales))
         elseif(bpower.eq.0) then
            tmp=FxFx_ren_scales(0)
         else
            tmp=1d0
         endif
         do i=1,nFxFx_ren_scales
            tmp=tmp*FxFx_ren_scales(i)
         enddo
         tmp=tmp**(1d0/max(dble(bpower),1d0))
         temp_scale_id='FxFx merging scale'
      elseif(imurtype.eq.1)then
        tmp=scale_global_reference(pp)
      elseif(imurtype.eq.2)then
        do i=nincoming+1,nexternal
          tmp=tmp+pt(pp(0,i))
        enddo
        temp_scale_id='sum_i pT(i), i=final state'
      elseif(imurtype.eq.3)then

         write (*,*) "imurtype=3 not possible in setscales.f: "/
     $        /"need to check number of Born orders."
         stop 1
         tmp1=0d0
         tmp2=1d0
         iqcd=0
         if (nint(wgtbpower).eq.0) then
            do i=nincoming+1,nexternal
               xm2=dot(pp(0,i),pp(0,i))
               if(xm2.le.0.d0)xm2=0.d0
               tmp1=tmp1+sqrt(pt(pp(0,i))**2+xm2)
            enddo
            tmp=tmp1
         else
            nn=0
            do i=1,nexternal
               if (is_a_j(i)) then
                  nn=nn+1
                  do j=0,3
                     pQCD(j,nn)=pp(j,i)
                  enddo
               endif
            enddo
            palg=1.d0
            rfj=0.4d0
            sycut=0d0
            call amcatnlo_fastjetppgenkt_timed(pQCD,NN,rfj,sycut,palg,
     $           pjet,njet,jet)
            if (nn-1.gt.nint(wgtbpower)) then
               write (*,*) 'More Born QCD partons than Born QCD '/
     &              /'couplings: cannot used this scale choice',imurtype
               stop
            elseif (nn-1.eq.nint(wgtbpower)) then
               do i=1,nint(wgtbpower)
                  tmp2=tmp2*sqrt(amcatnlo_fastjetdmergemax(nn-i-1))
               enddo
               tmp=tmp2**(1d0/wgtbpower)
            elseif (nn-1.lt.nint(wgtbpower)) then
               do i=nincoming+1,nexternal-1
                  if (.not.is_a_j(i)) then
                     xm2=dot(pp(0,i),pp(0,i))
                     if(xm2.le.0.d0)xm2=0.d0
                     tmp1=tmp1+sqrt(pt(pp(0,i))**2+xm2)
                  else
                     iqcd=iqcd+1
                     tmp2=tmp2*sqrt(amcatnlo_fastjetdmergemax(nn-iqcd-1))
                  endif
               enddo
               tmp=tmp2*tmp1**(nint(wgtbpower)-iqcd)
               tmp=tmp**(1d0/wgtbpower)
            endif
         endif
         temp_scale_id='geometric mean #3'
      else
        write(*,*)'Unknown option in muR_ref_dynamic',imurtype
        stop
      endif
      muR_ref_dynamic=tmp
c
      return
      end

      subroutine set_fac_scale(pp,muF)
      implicit none
      include 'genps.inc'
      include 'nexternal.inc'
      include 'run.inc'
      include 'coupl.inc'
      double precision pp(0:3,nexternal),muF(2)
      double precision muf_temp(2),muF_ref_dynamic
      character*80 muR_id_str,muF1_id_str,muF2_id_str,QES_id_str
      common/cscales_id_string/muR_id_str,muF1_id_str,
     #                         muF2_id_str,QES_id_str
      character*80 temp_scale_id,temp_scale_id2
      common/ctemp_scale_id/temp_scale_id
      double precision minscaleF
      parameter (minscaleF=2d0)
      temp_scale_id='  '
      temp_scale_id2='  '
      if(fixed_fac_scale)then
        muf_temp(1)=muF1_ref_fixed
        muf_temp(2)=muF2_ref_fixed
        temp_scale_id='fixed'
        temp_scale_id2='fixed'
      else
        muf_temp(1)=max(minscaleF,muF_ref_dynamic(pp))
        muf_temp(2)=muf_temp(1)
        temp_scale_id2=temp_scale_id
      endif
      muF(1)=muF1_over_ref*muf_temp(1)
      muF(2)=muF2_over_ref*muf_temp(2)
      muF12_current=muF(1)**2
      muF22_current=muF(2)**2
      muF1_id_str=temp_scale_id
      muF2_id_str=temp_scale_id2
      if(muF(1).le.0.d0.or.muF(2).le.0.d0)then
        write(*,*)'Error in set_fac_scale: muF(*)=',muF(1),muF(2)
        stop
      endif
      q2fact(1)=muF12_current
      q2fact(2)=muF22_current
      return
      end


      function muF_ref_dynamic(pp)
      implicit none
      include 'genps.inc'
      include 'nexternal.inc'
      include 'run.inc'
      double precision muF_ref_dynamic,pp(0:3,nexternal)
      double precision tmp,scale_global_reference,pt,et,dot,xm2,sumdot
      external pt,et,dot,sumdot
      character*80 temp_scale_id
      common/ctemp_scale_id/temp_scale_id
      integer i,imuftype
      parameter (imuftype=1)
      integer nFxFx_ren_scales
      double precision FxFx_ren_scales(0:nexternal),FxFx_fac_scale(2)
      common/c_FxFx_scales/FxFx_ren_scales,nFxFx_ren_scales
     $     ,FxFx_fac_scale
      tmp=0
      if(ickkw.eq.3)then
        tmp=(FxFx_fac_scale(1)+FxFx_fac_scale(2))/2d0
        temp_scale_id='FxFx merging scale'
      elseif(imuftype.eq.1)then
        tmp=scale_global_reference(pp)
      elseif(imuftype.eq.2)then
        do i=nincoming+1,nexternal
          tmp=tmp+pt(pp(0,i))**2
        enddo
        tmp=sqrt(tmp)
        temp_scale_id='Sqrt[sum_i pT(i)**2], i=final state'
      else
        write(*,*)'Unknown option in muF_ref_dynamic',imuftype
        stop
      endif
      muF_ref_dynamic=tmp
c
      return
      end


      subroutine set_QES_scale(pp,QES)
      implicit none
      include 'genps.inc'
      include 'nexternal.inc'
      include 'run.inc'
      include 'coupl.inc'
      include 'q_es.inc'
      double precision pp(0:3,nexternal),QES
      double precision QES_temp,QES_ref_dynamic
      double precision pi
      parameter (pi=3.14159265358979323846d0)
      character*80 muR_id_str,muF1_id_str,muF2_id_str,QES_id_str
      common/cscales_id_string/muR_id_str,muF1_id_str,
     #                         muF2_id_str,QES_id_str
      character*80 temp_scale_id
      common/ctemp_scale_id/temp_scale_id
      double precision minscaleES
      parameter (minscaleES=2d0)
      temp_scale_id='  '
      if(fixed_QES_scale)then
        QES_temp=QES_ref_fixed
        temp_scale_id='fixed'
      else
        QES_temp=max(minscaleES,QES_ref_dynamic(pp))
      endif
      QES=QES_over_ref*QES_temp
      QES2_current=QES**2
      QES_id_str=temp_scale_id
      if(QES.le.0.d0)then
        write(*,*)'Error in set_QES_scale: QES=',QES
        stop
      endif
      QES2=QES2_current
      return
      end


      function QES_ref_dynamic(pp)
      implicit none
      include 'genps.inc'
      include 'nexternal.inc'
      double precision QES_ref_dynamic,pp(0:3,nexternal)
      double precision tmp,scale_global_reference,pt
      external pt
      character*80 temp_scale_id
      common/ctemp_scale_id/temp_scale_id
      integer i,iQEStype
      parameter (iQEStype=1)
      tmp=0
      if (nincoming.eq.1) then
         tmp=pp(0,1) ! mass of the decaying particle
         temp_scale_id='Mass of decaying particle'
      elseif(iQEStype.eq.1)then
        tmp=scale_global_reference(pp)
      elseif(iQEStype.eq.2)then
        do i=nincoming+1,nexternal
          tmp=tmp+pt(pp(0,i))
        enddo
        temp_scale_id='sum_i pT(i), i=final state'
      else
        write(*,*)'Unknown option in QES_ref_dynamic',iQEStype
        stop
      endif
      QES_ref_dynamic=tmp
      return
      end


      function scale_global_reference(pp)
      implicit none
      include 'genps.inc'
      include 'nexternal.inc'
      include 'run.inc'
      include 'cuts.inc'
      double precision scale_global_reference,pp(0:3,nexternal)
      double precision tmp,pt,et,dot,xm2,sumdot,xmt2,ptmp(0:3),mf,mi,rdn
      external pt,et,dot,sumdot
      integer i,j,nbin
      character*80 temp_scale_id
      common/ctemp_scale_id/temp_scale_id
      tmp=0
      if(ickkw.eq.-1)then
         tmp=ptj
         temp_scale_id='NLO+NNLL veto scale: ptj_max'
      elseif(dynamical_scale_choice.eq.0) then
          tmp=muR_ref_fixed
          temp_scale_id='fixed scale'
      elseif(dynamical_scale_choice.eq.10) then
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc      USER-DEFINED SCALE: ENTER YOUR CODE HERE                                 cc
cc      to use this code you must set                                            cc
cc                 dynamical_scale_choice = 10                                   cc
cc      in the run_card (run_card.dat)                                           cc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         tmp=0d0
c$$$         mf=5000d0
c$$$         mi=1000d0
c$$$         nbin=40
         tmp=(pp(0,3)+pp(0,4))**2-(pp(1,3)+pp(1,4))**2
     $        -(pp(2,3)+pp(2,4))**2-(pp(3,3)+pp(3,4))**2
         tmp=dsqrt(tmp)/2d0
c$$$         rdn=tmp
c$$$         tmp=mi+(mf-mi)/nbin*(floor(nbin*(tmp-mi)/(mf-mi))+0.5d0)
c         temp_scale_id='scale = local mean M_bin '
c         call random_number(rdn)
c$$$         if (tmp.gt.mf) then
c$$$            write(*,*)'m min = ', mi,'mmax = ', mf,' mu = ',tmp
c$$$            write(*,*)' tmp = ',rdn,' xmul = ',floor(nbin*(rdn-mi)/(mf-mi))
c$$$         endif
c$$$  tmp=0d0
c$$$         tmp=(pp(0,3)+pp(0,4))**2-(pp(1,3)+pp(1,4))**2
c$$$     $        -(pp(2,3)+pp(2,4))**2-(pp(3,3)+pp(3,4))**2
c$$$         tmp=dsqrt(tmp)
         temp_scale_id='scale = M_inv/2 '
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc      USER-DEFINED SCALE: END OF USER CODE                                     cc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      else
        write(*,*)'Unknown option in scale_global_reference',dynamical_scale_choice
        stop
      endif
      scale_global_reference=tmp
      return
      end

