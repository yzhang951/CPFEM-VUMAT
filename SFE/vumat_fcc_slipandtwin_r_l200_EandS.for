c===================================================================
c                                                                  c
c VUMAT: Rate-dependent crystal plasticity with back stresses      c
c Basic crystal plasticity model for beginners                     c
c Developed by: Yin Zhang, Peking University                       c
c Contact: yinzhang@pku.edu.cn                                     c 
c Homepage: https://yzhang951.github.io/                           c
c Runs on Abaqus/Explicit (working on 6.13 and 2022)               c
c Latest version: 09/24/2024                                       c
c zzrewrite 11/21/2024                                                                 c
c-------------------------------------------------------------------
c This VUMAT models cyclic loading in FCC crystal plasticity
c ISVs: gamma_dot, back stress 1, back stress 2, threshold stress
c===================================================================

C Important: Double precision recommended
C num_el, path of 'aeuler', both need to be modified for each run

c-------------------------------------------------------------------
c VUSDFLD user redfine field variables
c Subroutine to define Euler angles and other field variables
c-------------------------------------------------------------------

       subroutine vusdfld(
c Read only variables -
     1 nblock, nstatev, nfieldv, nprops, ndir, nshr,
     2 jElem, kIntPt, kLayer, kSecPt,
     3 stepTime, totalTime, dt, cmname,
     4 coordMp, direct, dir_cos, charLength, props,
     5 stateOld,
c Wirte only variables
     6 stateNew, field)
c
       include 'vaba_param.inc'
       dimension jElem(nblock), coordMp(nblock,*),
     1           direct(nblock,3,3), dir_cos(nblock,3,3),
     2           charLength(nblock), props(nprops),
     3           stateOld(nblock,nstatev),
     4           stateNew(nblock,nstatev),
     5           field(nblock,nfieldv)
       character*80 cmname, FILE1, FILE2
c      
c     Local arrays from vgetvrm are dimensioned to
c     maximum block size (maxblk)
       parameter(nrData=6, num_el=14758)
       character*3 cData(maxblk*nrData)
       dimension rData(maxblk*nrData), jData(maxblk*nrData)
       real*8    psi(num_el,3), F_p(3,3), Weibull(num_el), wb
       PI = 4.D0*ATAN(1.D0)

       if(totalTime.le.dt.or.dt.eq.1.D0) then
c     Read in Euler angles from external file and store in STATENEW           
         write(*,*) "Read in Euler Angles"
         FILE1 = 'D:\dya\work1\abaqus\single\SFE\e14758
     &\aeuler'
         FILE2 = 'D:\dya\work1\Normal_Distribution
     &\normal14758l200sig5e_3.txt'
c         FILE2 = 'D:\twip_dis_channel\work_modify_v2\disturb
c     &\normal3D_N14758_ell200.txt'
         open(unit=65,file=FILE1,status='UNKNOWN')
         read(65,*)
         read(65,*)
         read(65,*) ICRYS_FLAG
         read(65,*)
         read(65,*)
         read(65,*)
         read(65,*)
         do i = 1,num_el
           read(65,*) ph, th, om 
           psi(i,1) = ph*PI/180.0
           psi(i,2) = th*PI/180.0
           psi(i,3) = om*PI/180.0
         end do
         close(unit=65)
         do k = 1,nblock
           i = jElem(k)
           stateNew(k,1) = psi(i,1)
           stateNew(k,2) = psi(i,2)
           stateNew(k,3) = psi(i,3)
c           field(k,1) = 1.0
         end do
         open(unit=60,file=FILE2,status='UNKNOWN')
         do i = 1,num_el
           read(60,*) wb
           Weibull(i) = wb
c           write(*,*) wb
         end do
         close(unit=60)
         do k = 1,nblock
             i = jElem(k)
           stateNew(k,121) = Weibull(i)
           stateOld(k,121) = Weibull(i)
         end do
       end if

       return
       end

c===================================================================
c     User defined materials, custom consitutive models here!      c
c           Detailed formulas in Equation.pdf                      c
c===================================================================
c-------------------------------------------------------------------
       subroutine vumat(
C Read only (unmodifiable) variables -
     1 nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2 stepTime, totalTime, dt, cmname, coordMp, charLength,
     3 props, density, strainInc, relSpinInc,
     4 tempOld, stretchOld, defgradOld, fieldOld,
     5 stressOld, stateOld, enerInternOld, enerInelasOld,
     6 tempNew, stretchNew, defgradNew, fieldNew,
C Write only (modifiable) variables -
     7 stressNew, stateNew, enerInternNew, enerInelasNew)
C
       include 'vaba_param.inc'
       dimension props(nprops), density(nblock), coordMp(nblock,*),
     1 charLength(nblock), strainInc(nblock,ndir+nshr),
     2 relSpinInc(nblock,nshr), tempOld(nblock),
     3 stretchOld(nblock,ndir+nshr), 
     4 defgradOld(nblock,ndir+nshr+nshr),
     5 fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6 stateOld(nblock,nstatev), enerInternOld(nblock),
     7 enerInelasOld(nblock), tempNew(nblock),
     8 stretchNew(nblock,ndir+nshr),
     9 defgradNew(nblock,ndir+nshr+nshr),
     1 fieldNew(nblock,nfieldv),
     2 stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     3 enerInternNew(nblock), enerInelasNew(nblock)
C       

       parameter(num_slip_sys = 12,   ! Total number of slip systems
     &           num_twin_sys = 12,   ! Total number of twin systems
     &           num_el = 14758,       ! Total number of elements, need to modify for different simulations!
     &           B_k = 1.38064852E-23)! Boltzmann constant in SI unit

       character*80 cmname
c-------------------------------------------------------------------
c  Dimension other arrays used in this UMAT sub
c-------------------------------------------------------------------
       real*8
     & array1(3,3),
     & array2(3,3),
     & array3(num_slip_sys),
     & del(3,3),            ! Kronecker delta tensor
     & dir_cos_0(3,3),      ! Direction cosine of original Euler angles
     & dir_cos(3,3),        ! Direction cosine of current Euler angles
     & rel_spin(3,3),       ! relative spin at start of inc
     & rel_spin_t(3,3),     ! relative spin at end of inc
     & U_stretch(3,3), U_stretch_inv(3,3), ! Stretch tensor and its inverse
     & C0(3,3,3,3),         ! Forth order elastic tensor at crystal basis
     & C(3,3,3,3),          ! Forth order elastic tensor at labortory basis
     & F_0(3,3), F_t(3,3),  ! Deformation gradient at start and end of inc
     & F_p_0(3,3),          ! Plastic deformation gradient at start of inc
     & F_p_t(3,3),          ! Plastic deformation gradient at end of inc
     & E_el(3,3),           ! Elastic Green strain tensor (reference config)
     & F_el_0(3,3),         ! Elastic deformation gradient at start of inc
     & F_el_inv_0(3,3),     ! Inverse of last one
     & F_el_t(3,3),         ! Elastic deformation gradient at end of inc
     & F_el_inv_t(3,3),     ! Inverse of last one
     & xL_p(3,3),           ! Plastic velocity gradient tensor
     & tau(num_slip_sys),   ! Resovled shear stress on each slip system
     & rho_0(num_slip_sys), ! Dislocation density at each slip system (optional)
     & rho_t(num_slip_sys), ! Last one at end of inc
     & back_stress_one_0(num_slip_sys), ! Back stress 1, start of inc
     & back_stress_one_t(num_slip_sys), ! Back stress 1, end of inc
     & back_stress_two_0(num_slip_sys), ! Back stress 2, start of inc
     & back_stress_two_t(num_slip_sys), ! Back stress 2, end of inc
     & threshold_0(num_slip_sys),       ! Threshold resistance, start of inc
     & threshold_t(num_slip_sys),       ! Threshold resistance, end of inc
     & gamma_dot_0(num_slip_sys),       ! Shear rate, start of inc
     & gamma_dot_t(num_slip_sys),       ! Shear rate, end of inc
     & psi_t(3),                        ! Euler angles
     & spk2(3,3),                       ! Second Piola-Kirchoff stress (sym and ref config)
     & sig(3,3), sig_star(3,3),         ! Cauchy stress in lab basis and corotational basis
     & xs0(3,num_slip_sys),             ! Slip direction in lab basis before elastic deformation
     & xs0_xtal(3,num_slip_sys),        ! Slip direction in crystal basis
     & xs(3,num_slip_sys),              ! Slip direction in lab basis after elastic deformation
     & xm0(3,num_slip_sys),             ! Slip plane in lab basis before elastic deformation
     & xm0_xtal(3,num_slip_sys),        ! Slip plane in crystal basis
     & xm(3,num_slip_sys),              ! Slip plane in lab basis after elastic deformation
     & Aii(num_slip_sys,num_slip_sys)   ! Dislocation forest hardening coeff (also called latent hardening)

       real*8 C11, C12, C44, shear_mod_0, ! Elastic moduli
     & psi_ang, theta_ang, phi_ang,
     & g0, am, s0, h0, s_sat, ah, A00, A01,
     & dtime, PI, one, plastic_work_inc, dPEEQ,
     & T, mu, burgers
       
C --- Variable definition for forest dislocation hardening -----
       real*8 h_dis, omega1, omega2, tau_dis(num_slip_sys)
       real*8 k1, k2, rho
C --- Variable definition for the Hall-Petch effect ----------------
       real*8 k_hp, d_grain, tau_hp
C --- Variable definition for irradiation defects hardening --------
       integer num_def_sys
       real*8 h_def, omega3, omega4, tau_def(num_slip_sys)
       real*8 rho_def_0, rho_def_t, d_def, rho_def, P_anni
C --- Variable definition for lattice resistance -------------------
       real*8 tau_0, tau_c(num_slip_sys)
C --- Variable definition for twin systems -------------------------
       real*8 xs0_xtal_tw(3,num_twin_sys), xm0_xtal_tw(3,num_twin_sys)
       real*8 xs0_tw(3,num_twin_sys), xm0_tw(3,num_twin_sys)
C --- Var def for interactive relationship between slip and twin ---
       real*8 affiliation_sl_tw(num_slip_sys,num_twin_sys)
       real*8 affiliation_tw_tw(num_twin_sys,num_twin_sys)
       real*8 xi_alpbet_cop,xi_alpbet_noncop
       real*8 xi_alpbet(num_slip_sys,num_twin_sys)
       real*8 xi_betbet_cop,xi_betbet_noncop
       real*8 xi_betbet(num_twin_sys,num_twin_sys)
C --- Variable definition for deformation twinning -----------------
       real*8 gamma_dot_tw_0(num_twin_sys), f_beta_0(num_twin_sys)
       real*8 gamma_dot_tw_t(num_twin_sys), f_beta_t(num_twin_sys)
       real*8 tau_tw(num_twin_sys), mfp_twin(num_twin_sys)
       real*8 volume_beta(num_twin_sys), nu, x_c
       real*8 gamma_sf, tau_r_tw, tau_f_tw, L_tw
       real*8 nucrhorate_tw(num_twin_sys), Vcs, m_tw
       real*8 A_tw, f_beta_max, f_beta_dot(num_twin_sys)
       real*8 gamma_tw
C --- Variable definition for slip resistance caused by twins ------
       real*8 N_alpha_0(num_slip_sys), N_alpha_t(num_slip_sys)
       real*8 t_tw, N_max, tau_TBs(num_slip_sys), l_avedistance
       real*8 mfp_slip(num_slip_sys)
C --- Var def for the energy of dislocation cross-slip through def -    
       real*8 A_cs, omega_cs1, omega_cs2, omega_cs3, p_cs, q_cs
C --- Var def for the energy of dislocation cutting through def ----
       real*8 A_cut, alpha_cut, k_cut, m_cut
C --- Variable definition for the probability of cutting def -------
       real*8 P_cut(num_slip_sys)
       integer affiliation_sl_cs(num_slip_sys)
       real*8 tauce(num_slip_sys)
       real*8 tauge(num_slip_sys), taucs(num_slip_sys)
       real*8 E_cut(num_slip_sys)
       real*8 E_cs_positive(num_slip_sys), E_cs_negative(num_slip_sys)
C --- Variable definition for introducing randomness  --------------       
       real*8 Weibull(num_el), wb
       real*8 d_def_rand
       
       PI = 4.D0*DATAN(1.D0)
       one = 1.D0
       dtime = dt                       ! Timestep
c===================================================================
c      TOTAL ISVs:
c===================================================================
c      3  - Euler angles
c      9  - relative Spin Matrix
c      9  - F_p(i,j)
c      12 - dummy variables for 12 slip systems
c      12 - back_stress(12) for 12 slip systems
c      12 - threshold(12) for 12 slip systems
c      12 - gamma_dot(12) for 12 slip systems
c===================================================================

C       write(*,*) "Material Properties"
    
c===================================================================
c       Assign props() array to logical variable names
c===================================================================
c-------------------- Macroscopic Properties  ----------------------
       C11          = props(1)  ! Elastic constants C11
       C12          = props(2)  ! Elastic constants C12
       C44          = props(3)  ! Elastic constants C44, shear modulus
       psi_ang      = props(4)  ! Euler angles
       theta_ang    = props(5)
       phi_ang      = props(6)
c-------------------- Rate-dependent flow rule ---------------------       
       g0           = props(7)  ! Reference shear rate
       am           = props(8)  ! Strain rate sensitivity
       s0           = props(9)  ! Initial slip resistance
c-------------------- Hardening ------------------------------------
       h0           = props(10) ! Initial hardening rate
       s_sat        = props(11) ! Saturated slip resistance
       ah           = props(12) ! Hardening exponent 
c-------------------- Threshold Stress Evolution -------------------       
       A00          = props(13) ! Dislocation self interaction coeff
       A01          = props(14) ! Dislocation interaction coeff
c-------------------- Backstress Evolution -------------------------
       C_back_one   = props(15) ! Backstress 1 generation parameter
       C_dyn_one    = props(16) ! Backstress 1 dynamic recovery para
       C_back_two   = props(17) ! Backstress 2 generation parameter
       C_dyn_two    = props(18) ! Backstress 2 dynamic recovery para
C-------------------- Dislocation density evolution ----------------
       k1           = props(19) ! Dislocation mutiplication 
       k2           = props(20) ! Dislocation annihilation
C-------------------- Heat activates -------------------------------
       T            = props(21) ! Temperature
C-------------------- Additional Macroscopic Properties ------------
       burgers      = props(22) ! Magnitude of the Burgers vector
       rho          = props(23) ! Initial dislocation density
C-------------------- Forest dislocation hardening -----------------
       h_dis        = props(24) ! Dislocation hardening coefficient
       omega1       = props(25) ! Interaction coefficient
       omega2       = props(26) ! Interaction coefficient
C-------------------- The Hall-Petch effect ------------------------
       k_hp         = props(27) ! HP coefficient
       d_grain      = props(28) ! Average grain size
C-------------------- Irradiation defects hardening ----------------
       h_def        = props(29) ! HP coefficient
       omega3       = props(30) ! Interaction coefficient
       omega4       = props(31) ! Interaction coefficient
       rho_def      = props(32) ! Initial irr def density
       d_def        = props(33) ! Average irr def size
C-------------------- Lattice resistance ---------------------------
       tau_0        = props(34) 
C-------------------- Mean free paths ------------------------------
       xi_alpbet_cop    = props(35) ! Interaction coefficient between slip and twin systems (coplanar)
       xi_alpbet_noncop = props(36) ! Interaction coefficient between slip and twin systems (cross-slip)
       xi_betbet_cop    = props(37) ! Interaction coefficient between twin systems (coplanar) 
       xi_betbet_noncop = props(38) ! Interaction coefficient between twin systems (non-coplanar) 
       t_tw             = props(39) ! Average thickness of twin lamellas
C-------------------- TBs hardening for slip -----------------------
       N_max            = props(40) ! The saturated number of piled-up dislocations at boundaries
       l_avedistance    = props(41) ! Average distance between slip bands
C-------------------- Twin nucleation ------------------------------
       nu           = props(42) ! Poisson ratio
       x_c          = props(43) ! Critical distance to form the twin nucleus
       gamma_sf     = props(44) ! Stacking fault energy
       L_tw         = props(45) ! Width of twin embryo
       Vcs          = props(46) ! Activation volume for cross-slip
       A_tw         = props(47) ! Parameter in the probability of forming twin nuclei
       f_beta_max   = props(48) ! The maximum of twin volume fraction in each twin system
       m_tw         = props(49) ! Twinning rate sensitivity
C-------------------- irr def evolution------------------------------       
       A_cs         = props(50) ! The fitting parameters of the energy barrier 
       omega_cs1    = props(51) ! for cross-slip when encountering irr defects
       omega_cs2    = props(52)
       omega_cs3    = props(53)
       p_cs         = props(54)
       q_cs         = props(55)
       A_cut        = props(56) ! The fitting parameters of the energy barrier 
       alpha_cut    = props(57) ! for cutting though when encountering irr defects
       k_cut        = props(58)
       m_cut        = props(59)
       P_anni       = props(60) ! The probability of annihilation of the defect.
       
       
       
       mu = dsqrt((C11-C12)*C44/2)
       num_def_sys = 4                  
       gamma_tw = sqrt(2.d0)/2.d0       ! The characteristic twinning shear strain in FCC
c------------ Euler angles for single xtal simulation --------------
       psi_t(1) = psi_ang*PI/180.
       psi_t(2) = theta_ang*PI/180.
       psi_t(3) = phi_ang*PI/180.

c--------------------- Disl interaction coeff matrix ---------------
       do i = 1,num_slip_sys
         do j = 1,num_slip_sys
           Aii(i,j) = A01
         end do
         Aii(i,i) = A00
       end do

C===================================================================
C===================================================================
C
C
C     Abaqus Explicit sends in data in blocks of NBLOCK=128
C     elements at a time. 
C
C===================================================================       
C=================================================================== 
C===================================================================       
C      Start loop over Nblock       
C===================================================================
       do km = 1,nblock
c===================================================================       
c      Initialize internal state variables for the first time step
c===================================================================

         if (totalTime.le.dt.or.dt.eq.one) then
c           write(*,*) "Initializing ISVs, first time step"
c--------------- Read in Euler angles from external file -----------
c--- Important!: comment this section for poly xtal simulations ----
c--- the Euler angles will be read from file 'aeuler' through VUSDFLD
c           psi_t(1) = psi_ang*PI/180.0
c           psi_t(2) = theta_ang*PI/180.0
c           psi_t(3) = phi_ang*PI/180.0
           psi_t(1) = stateNew(km,1)
           psi_t(2) = stateNew(km,2)
           psi_t(3) = stateNew(km,3)
c--------------------- Initialize F_p and F_p_inv. -----------------
           do i = 1,3
             do j = 1,3
               F_p_0(i,j) = 0.D0
               F_el_0(i,j) = 0.D0
               rel_spin(i,j) = 0.D0
             end do
             F_p_0(i,i) = 1.D0
             F_el_0(i,i) = 1.D0
             rel_spin(i,i) = 1.D0
           end do
c------------ Initialize shear rate and back stress ----------------
           do i = 1,num_slip_sys
             gamma_dot_0(i) = 0.D0
             back_stress_one_0(i) = 0.D0
             back_stress_two_0(i) = 0.D0
           end do
c-------------------- Initialize threshold stress -------------------
           do i = 1,num_slip_sys
             threshold_0(i) = s0
           end do
c-------------------- Initialize dislocation density ----------------           
           do i = 1,num_slip_sys
C             stateOld(km,i) = rho
              rho_0(i) = rho
           end do   
c-------------------- Initialize def tw shear rate ------------------         
           do i = 1,num_twin_sys
               gamma_dot_tw_0(i) = 0.d0
           enddo
c-------------------- Initialize twin volume fraction ---------------           
           do i = 1,num_twin_sys
               f_beta_0(i) = 0.d0
           enddo
c-------------------- Initialize piled-up dislocation ---------------           
           do i = 1,num_slip_sys
               N_alpha_0(i) = 0.d0
           enddo
c-------------------- Initialize def den ----------------------------  
           rho_def_0 = statenew(km,121)*rho_def           
c           rho_def_0 = 1.d0*rho_def       
c-------------------- Initialize def size ---------------------------  
c           d_def_rand = statenew(km,121)*d_def      
           d_def_rand = 1.d0*d_def    
C           write(*,*) "Initialization finished!"
         end if
c====================================================================         
c       End of initialization. 
c       Read in internal state variables for nonzero time step
c====================================================================
         if (totalTime.gt.dt.and.dt.ne.one) then
c           write(*,*) "Read in ISVs, nonzero time step"
           n = 0
c-------------------- Read in Euler Angles 1-3 ----------------------
           do i = 1,3
             n = n + 1
             psi_t(i) = stateOld(km,n) 
           end do 
c-------------------- Read in Relative Spin Matrix 4-12--------------
           do i = 1,3
             do j = 1,3
               n = n + 1
               rel_spin(i,j) = stateOld(km,n)
             end do
           end do 
c-------------------- Read in F_p 13-21 -----------------------------
           do i = 1,3
             do j = 1,3
               n = n + 1
               F_p_0(i,j) = stateOld(km,n)
             end do
           end do
c-------------------- Read in dislocation shear rate 22-33 ----------
           do i = 1, num_slip_sys
             n = n + 1
             gamma_dot_0(i) = stateOld(km,n)
           end  do
c-------------------- Read in accumlated shear 34-45 ----------------
           do i = 1, num_slip_sys
             n = n + 1
             rho_0(i) = stateOld(km,n)
           end do
c-------------------- Read in threshold stress 46-57 ----------------
           do i = 1, num_slip_sys
             n = n + 1
c             threshold_0(i) = stateOld(km,n)
           end do
c-------------------- Read in back stress one 58-69 -----------------
           do i = 1, num_slip_sys
             n = n + 1
c             back_stress_one_0(i) = stateOld(km,n)
           end do
c-------------------- Read in back stress two 70-81 -----------------
           do i = 1, num_slip_sys
             n = n + 1
c             back_stress_two_0(i) = stateOld(km,n)
           end do
c-------------------- Read in PEEQ 82 -------------------------------
           n = n + 1
           n = n + 1 
c-------------------- Read in def density 84 ------------------------ 
           n = n + 1
           rho_def_0 = stateold(km,n)
c-------------------- Read in def tw shear rate 85-96 ---------------
           do i = 1,num_twin_sys
               n = n + 1
               gamma_dot_tw_0(i) = stateold(km,n)
           enddo
c-------------------- Read in twin volume fraction 97-108 -----------
           do i = 1,num_twin_sys
               n = n + 1
               f_beta_0(i) = stateold(km,n)
           enddo      
c-------------------- Read in piled-up dislocation 109-120 ----------
           do i = 1,num_twin_sys
               n = n + 1
               N_alpha_0(i) = stateold(km,n)
           enddo  
c-------------------- Read in random def size -----------------------
c           d_def_rand = stateold(km,121)*d_def
           d_def_rand = 1.d0*d_def    
C           write(*,*) "ISVs read in finished!"
         end if

c====================================================================         
c        End of ISVs read in. 
c====================================================================

c====================================================================       
c       Define Dir_cos and Slip System
c====================================================================
C         write(*,*) "Define Dir_cos and Slip system!"
         call calc_Dir_Cos(psi_t, dir_cos_0)
         call aa_dot_bb(3, dir_cos_0, rel_spin, dir_cos)
C		 Define slip systems for FCC materials
C        Tricks: xs and xm can be swapped for BCC materials
         call def_Slip_Sys(xm0_xtal, xs0_xtal, num_slip_sys)
         call def_Twin_Sys(xm0_xtal_tw, xs0_xtal_tw, num_twin_sys)
         call calc_affiliation_matrix(xm0_xtal, xs0_xtal, 
     &   xm0_xtal_tw, xs0_xtal_tw, num_twin_sys, num_slip_sys,
     &   affiliation_sl_tw)
         call calc_affiliation_matrix(xm0_xtal_tw, xs0_xtal_tw, 
     &   xm0_xtal_tw, xs0_xtal_tw, num_twin_sys, num_twin_sys,
     &   affiliation_tw_tw)
         call calc_xi(num_slip_sys, num_twin_sys, xi_alpbet_cop,
     &   xi_alpbet_noncop, affiliation_sl_tw, xi_alpbet)
         call calc_xi(num_twin_sys, num_twin_sys, xi_betbet_cop,
     &   xi_betbet_noncop, affiliation_tw_tw, xi_betbet)
         
C        Rotate the slip direction and plane to the lab basis
         do n = 1, num_slip_sys
           do i = 1,3
             xs0(i,n) = 0.D0
             xm0(i,n) = 0.D0
             do j = 1,3
               xs0(i,n) = xs0(i,n) + dir_cos_0(i,j)*xs0_xtal(j,n)
               xm0(i,n) = xm0(i,n) + dir_cos_0(i,j)*xm0_xtal(j,n)
             end do
           end do
         end do
      
         do n = 1, num_twin_sys
           do i = 1,3
             xs0_tw(i,n) = 0.D0
             xm0_tw(i,n) = 0.D0
             do j = 1,3
               xs0_tw(i,n) = 
     &         xs0_tw(i,n) + dir_cos_0(i,j)*xs0_xtal_tw(j,n)
               xm0_tw(i,n) = 
     &         xm0_tw(i,n) + dir_cos_0(i,j)*xm0_xtal_tw(j,n)
             end do
           end do
         end do      

c====================================================================       
c       Define Kron and Elastic Tensor
c====================================================================
C         write(*,*) "Define Kron and Elastic Tensor!"
         call def_Kron_Del(del)
         call calc_4th_C(C11, C12, C44, del, C0)
c====================================================================       
c       Read in deformation gradient tensor
c====================================================================
C         write(*,*) "Read in deformation gradient tensor!"
C         write(*,*) "Read in deformation gradient tensor!"
         if (nshr.eq.3) then
           F_0(1,1) = defgradOld(km,1)
           F_0(2,2) = defgradOld(km,2)
           F_0(3,3) = defgradOld(km,3)
           F_0(1,2) = defgradOld(km,4)
           F_0(2,3) = defgradOld(km,5)
           F_0(3,1) = defgradOld(km,6)
           F_0(2,1) = defgradOld(km,7)
           F_0(3,2) = defgradOld(km,8)
           F_0(1,3) = defgradOld(km,9)


           F_t(1,1) = defgradNew(km,1)
           F_t(2,2) = defgradNew(km,2)
           F_t(3,3) = defgradNew(km,3)
           F_t(1,2) = defgradNew(km,4)
           F_t(2,3) = defgradNew(km,5)
           F_t(3,1) = defgradNew(km,6)
           F_t(2,1) = defgradNew(km,7)
           F_t(3,2) = defgradNew(km,8)
           F_t(1,3) = defgradNew(km,9)
         elseif (nshr.eq.1) then
c        Plane Strain problem         
           F_0(1,1) = defgradOld(km,1)
           F_0(2,2) = defgradOld(km,2)
           F_0(3,3) = defgradOld(km,3)
           F_0(1,2) = defgradOld(km,4)
           F_0(2,3) = 0.d0
           F_0(3,1) = 0.d0
           F_0(2,1) = defgradOld(km,5)
           F_0(3,2) = 0.d0
           F_0(1,3) = 0.d0


           F_t(1,1) = defgradNew(km,1)
           F_t(2,2) = defgradNew(km,2)
           F_t(3,3) = defgradNew(km,3)
           F_t(1,2) = defgradNew(km,4)
           F_t(2,3) = 0.d0
           F_t(3,1) = 0.d0
           F_t(2,1) = defgradNew(km,5)
           F_t(3,2) = 0.d0
           F_t(1,3) = 0.d0
         endif

c====================================================================       
c       Rotate elastic tensor and slip direction and normal direction
c====================================================================
C         write(*,*) "Rotate Elastic tensor and slip system!"
         call rotate_4th(dir_cos_0, C0, C)
         call calc_F_el(F_p_0, F_0, F_el_inv_0, F_el_0)
         call calc_Rot_Slip(xs0, xs, F_el_0, num_slip_sys)
         call calc_Rot_Norm(xm0, xm, F_el_inv_0, num_slip_sys)
c====================================================================       
c        Calculate the plastic Velocity Gradient, Green Strain, 
C        SPK2 stress and Cauchy Stress
c        See 'Equation.pdf'
c====================================================================
C         write(*,*) "Calculate the L_p, E_el, Spk2, Sig"
c         call calc_L_p(num_slip_sys, gamma_dot_0, xs0, xm0, xL_p)
         call calc_L_p_new(num_slip_sys, gamma_dot_0, xs0, xm0, xL_p,
     &   num_twin_sys, gamma_dot_tw_0, xs0_tw, xm0_tw, f_beta_0)
         call calc_F_p_subend(xL_p, dtime, F_p_0, F_p_t)
         call calc_F_el(F_p_t, F_t, F_el_inv_t, F_el_t)
         call calc_E_el(F_el_t, E_el)
         call calc_Spk2(C, E_el, spk2)
         call calc_Sig(F_el_t, spk2, sig)
         call calc_Tau(num_slip_sys, tau, spk2, xs0, xm0)
         call calc_Tau(num_twin_sys, tau_tw, spk2,
     &   xs0_tw, xm0_tw)

c==================================================================== 
c        Update the rotation matrix
c====================================================================
c         write(*,*) "Updating the rotation matrix", F_p_0, F_p_t
         if (nshr.eq.3) then
           U_stretch(1,1) = stretchNew(km,1)
           U_stretch(2,2) = stretchNew(km,2)
           U_stretch(3,3) = stretchNew(km,3)
           U_stretch(1,2) = stretchNew(km,4)
           U_stretch(2,3) = stretchNew(km,5)
           U_stretch(3,1) = stretchNew(km,6)
           U_stretch(3,2) = U_stretch(2,3)
           U_stretch(1,3) = U_stretch(3,1)
           U_stretch(2,1) = U_stretch(1,2)
         elseif (nshr.eq.1) then
c        Plane Strain problem   
           U_stretch(1,1) = stretchNew(km,1)
           U_stretch(2,2) = stretchNew(km,2)
           U_stretch(3,3) = stretchNew(km,3)
           U_stretch(1,2) = stretchNew(km,4)
           U_stretch(2,3) = 0.d0
           U_stretch(3,1) = 0.d0
           U_stretch(3,2) = U_stretch(2,3)
           U_stretch(1,3) = U_stretch(3,1)
           U_stretch(2,1) = U_stretch(1,2)  
         endif
         call inverse_3x3(U_stretch, U_stretch_inv)
         call aa_dot_bb(3, F_t, U_stretch_inv, rel_spin_t)

c====================================================================       
c        Calculate the evolution of internal state variables
c====================================================================

C         write(*,*) "Evolution of internal state variables"
c         call calc_back(num_slip_sys, C_back_one, C_dyn_one, 
c     &   gamma_dot_0, back_stress_one_0, dtime, back_stress_one_t)
c         call calc_back(num_slip_sys, C_back_two, C_dyn_two, 
c     &   gamma_dot_0, back_stress_two_0, dtime, back_stress_two_t)
C         write(*,*) "Evolution of back stress: ", back_stress_t
c         call calc_threshold(num_slip_sys, Aii, h0, threshold_0, 
c     &   s_sat, ah, gamma_dot_0, dtime, threshold_t)
c         write(*,*) "Evolution of threshold stress: ", threshold_t
c------- Dislocation slipping ---------------------------------------
         call calc_tau_dis(num_slip_sys, mu, burgers, h_dis,
     &   omega1, omega2, rho_0, tau_dis)
         call calc_tau_hp(k_hp, d_grain, tau_hp)
         call calc_tau_def(num_def_sys, num_slip_sys, mu, burgers,
     &   h_def, omega3, omega4, rho_def_0, d_def_rand, tau_def)
         call calc_tau_TBs(num_slip_sys, num_twin_sys, mu, burgers,
     &   d_grain, xi_alpbet, f_beta_0, t_tw, N_alpha_0, N_max, tau_TBs,
     &   mfp_slip)
         call calc_crss(num_slip_sys, tau_0, tau_hp, tau_dis, 
     &   tau_def, tau_TBs, tau_c)
c         call f(num_slip_sys, g0, tau, back_stress_one_t, 
c     &   back_stress_two_t, threshold_t, am, gamma_dot_t)
         call fnew(num_slip_sys, g0, tau, tau_TBs, tau_c, am,
     &   gamma_dot_t)
         call calc_P_cut(num_slip_sys, xs0_xtal, xs0, xm0, spk2,
     &   tau_def, gamma_sf, mu, burgers, A_cs, omega_cs1, omega_cs2,
     &   omega_cs3, p_cs, q_cs, rho_def_0, d_def_rand, A_cut, alpha_cut,
     &   k_cut, m_cut, B_k, T, P_cut, affiliation_sl_cs,
     &   tauce, tauge, taucs, E_cs_positive, E_cs_negative, E_cut)
         call calc_rho(num_slip_sys, gamma_dot_t, rho_0, 
     &   dtime, rho_t, k1, k2, N_alpha_0, N_max, mfp_slip, burgers,
     &   rho_def_0, d_def_rand, P_cut, affiliation_sl_cs)
         call calc_rho_def(num_slip_sys, gamma_dot_t, 
     &   burgers, d_def_rand, rho_def_0, dtime, rho_def_t, 
     &   P_cut, tau, P_anni)
         call calc_rho_piledup(num_slip_sys, N_alpha_0, l_avedistance,
     &   burgers, N_max, gamma_dot_t, dtime, N_alpha_t)
C         write(*,*) "Evolution of gamma_dot: ", gamma_dot_t
c------- Dubug ------------------------------------------------------
         do i = 1,num_slip_sys
             if ((abs(tau(i))-tau_TBs(i))/tau_c(i).gt.4.d0
     &        .and. totaltime.gt.0.d0)then
                 write(*,*) 'RSS is too large'
                 write(*,*) 'eular:' , psi_t
                 call abort
             endif
         enddo
      
c------- Deformation twinning ---------------------------------------         
         call calc_newtw_volume(num_twin_sys, d_grain, 
     &   xi_betbet, f_beta_0, t_tw, mfp_twin, volume_beta)
         call calc_tau_r_tw(mu, burgers, nu, x_c,
     &   gamma_sf, tau_r_tw)
         call calc_tau_formation_tw(mu, burgers, gamma_sf,
     &   d_grain, tau_f_tw)
         call calc_nucrhorate_tw(num_slip_sys, num_twin_sys, L_tw,
     &   affiliation_sl_tw, gamma_dot_0, rho_0, tau, tau_TBs, tau_c,
     &   nucrhorate_tw, m_tw, totaltime)
         call calc_rate_ftw(num_twin_sys, volume_beta, tau_r_tw,
     &   tau_f_tw, nucrhorate_tw, f_beta_0, Vcs, B_k, T, A_tw,
     &   f_beta_max, tau_tw, f_beta_dot)
         call fnew_tw(num_twin_sys, gamma_tw, f_beta_dot, tau_tw,
     &   gamma_dot_tw_t)
         call calc_ftw(num_twin_sys, f_beta_0, f_beta_dot,
     &   dtime, f_beta_t)
c==================================================================== 
c        Update the stress tensor in corotational basis system
c        sig_star = R^T * sig * R, Eq. (8)
c====================================================================
         call aa_dot_bb(3, sig, rel_spin_t, array1)
         call transpose(3, rel_spin_t, array2)
         call aa_dot_bb(3, array2, array1, sig_star)
         do i = 1,ndir
           stressNew(km,i) = sig_star(i,i)
         end do
         stressNew(km,4) = sig_star(1,2)
         stressNew(km,5) = sig_star(2,3)
         stressNew(km,6) = sig_star(3,1)

c==================================================================== 
c        Update the internal energy
c        dE = sig * de / density
c====================================================================
C         write(*,*) "Update the internal energy"
         enerInternNew(km) = enerInternOld(km)
         do i = 1, ndir+nshr
           enerInternNew(km) = enerInternNew(km) + 
     &     0.5*(stressOld(km,i)+stressNew(km,i))
     &     *strainInc(km,i)/density(km)
           if(i.gt.ndir) then
             enerInternNew(km) = enerInternNew(km) + 
     &       0.5*(stressOld(km,i)+stressNew(km,i))
     &       *strainInc(km,i)/density(km)
           end if
         end do

         plastic_work_inc = 0.D0
         do i = 1, num_slip_sys
           plastic_work_inc = plastic_work_inc +
     &       DABS(tau(i)*gamma_dot_t(i)*dtime) +
     &       DABS(tau_tw(i)*gamma_dot_tw_t(i)*dtime)
         end do
         enerInelasNew(km) = enerInelasOld(km) + 
     &                       plastic_work_inc/density(km)

         call calc_dPEEQ(sig, plastic_work_inc, dPEEQ)


c====================================================================    
c        Store in internal state variables for nonzero time step
c====================================================================
!         write(*,*) "Store in state variables "

         n = 0
c-------------------- Store in Euler Angles 1-3 ---------------------
         do i = 1,3
           n = n + 1
           stateNew(km,n) = psi_t(i)
         end do 
c-------------------- Store in Relative Spin Matrix 4-12-------------
         do i = 1,3
           do j = 1,3
             n = n + 1
             stateNew(km,n) = rel_spin_t(i,j)
           end do
         end do 
c-------------------- Store in F_p 13-21 ----------------------------
         do i = 1,3
           do j = 1,3
             n = n + 1
             stateNew(km,n) = F_p_t(i,j) 
           end do
         end do
c-------------------- Store in dislocation shear rate 22-33 ---------
         do i = 1, num_slip_sys
           n = n + 1
           stateNew(km,n) = gamma_dot_t(i)
         end  do
c-------------------- Store in dislocation density 34-45 ------------
         do i = 1, num_slip_sys
           n = n + 1
           stateNew(km,n) = rho_t(i)
         end do
c-------------------- Store in threshold stress 46-57 ---------------
         do i = 1, num_slip_sys
           n = n + 1
c           stateNew(km,n) = threshold_t(i)
c           stateNew(km,n) = one
           stateNew(km,n) = tau_def(i)
c           stateNew(km,n) = mfp_slip(i)
c           stateNew(km,n) = abs(tau_tw(i))
c           stateNew(km,n) = stateOld(km,n)
c     &       +dabs(gamma_dot_t(i))* dtime*
c     &      (1.d0-N_alpha_t(i)/N_max)/(burgers*mfp_slip(i))
         end do
c-------------------- Store in Back stress 58-69 --------------------
         do i = 1, num_slip_sys
           n = n + 1
c           stateNew(km,n) = back_stress_one_t(i)
c           stateNew(km,n) = one
c           stateNew(km,n) = 
c     &            1.d0-exp(-Vcs/(B_k*T)*(tau_r_tw-abs(tau_tw(i))))
           stateNew(km,n) = tau_dis(i)
         end do
c-------------------- Store in Back stress 70-81 --------------------
         do i = 1, num_slip_sys
           n = n + 1
c           stateNew(km,n) = back_stress_two_t(i)
c           stateNew(km,n) = one
c           stateNew(km,n) = exp(-(tau_f_tw/abs(tau_tw(i)))**A_tw)
c           stateNew(km,n) = tau_f_tw
c           stateNew(km,n) = mfp_slip(i)
           stateNew(km,n) = 1.d0 - P_cut(i)
         end do
c-------------------- Store in PEEQ 82 ------------------------------
         n = n + 1
         stateNew(km,n) = stateOld(km,n) + dPEEQ
c-------------------- Store in plastic disspation 83-----------------
         n = n + 1
         stateNew(km,n) = enerInelasNew(km)
c-------------------- Store in irr def density 84--------------------
         n = n + 1
         stateNew(km,n) = rho_def_t  
c-------------------- Store in def tw shear rate 85-96 --------------
         do i = 1,num_twin_sys
             n = n + 1
             stateNew(km,n) = gamma_dot_tw_t(i)
         enddo
c-------------------- Store in twin volume fraction 97-108 ----------
         do i = 1,num_twin_sys
             n = n + 1
             stateNew(km,n) = f_beta_t(i)
         enddo      
c-------------------- Store in piled-up dislocation 109-120 ---------
         do i = 1,num_twin_sys
             n = n + 1
             stateNew(km,n) = N_alpha_t(i)
         enddo     
c-------------------- Store in Weibull number 121 -------------------
         n = n + 1
         stateNew(km,n) = stateold(km,n)
c         write(*,*) "Next Element!"
c-------------------- sdv used to debug -----------------------------
c-------------------- 122 -------------------------------------------         
         n = n + 1
         stateNew(km,n) = sum(f_beta_t)
c         stateNew(km,n) = d_def_rand
c-------------------- 123-134 ---------------------------------------         
         do i = 1, num_slip_sys
           n = n + 1
c           stateNew(km,n) = back_stress_two_t(i)
c           stateNew(km,n) = one
c           stateNew(km,n) = exp(-(tau_f_tw/abs(tau_tw(i)))**A_tw)
c           stateNew(km,n) = tau_f_tw
           stateNew(km,n) = tau_c(i)
         end do
c-------------------- 135-146 ---------------------------------------         
         do i = 1, num_slip_sys
           n = n + 1
c           stateNew(km,n) = nucrhorate_tw(i)
           stateNew(km,n) = abs(tau(i))
         end do  
c-------------------- 147-158 ---------------------------------------               
         do i = 1, num_slip_sys
           n = n + 1
           stateNew(km,n) = (abs(tau(i))-tau_TBs(i))/tau_c(i)
         end do   
c-------------------- 159-170 ---------------------------------------       
         do i = 1, num_slip_sys
           n = n + 1
c           stateNew(km,n) = mfp_slip(i)SIGN(1.D0, tau(i)-B1(i)-B2(i))
           stateNew(km,n) = SIGN(1.d0,tau_c(i) - stateOld(km,122+i))
         end do
c-------------------- 171 ---------------------------------------      
         n = n + 1
         stateNew(km,n) = statenew(km,121)*rho_def - rho_def_t
              
       end do
     
       return
       end 

c=====================================================================
c=====================================================================
c      SUBROUTINE:  
c      Calculate the function f = [shear rate expression]
c      Eq. (11)
c=====================================================================
c=====================================================================	
       subroutine f(num_slip_sys, g0, tau, B1, B2, S, am, gamma_dot)
       implicit double precision (a-h,o-z,k)
       real*8 tau(num_slip_sys), B1(num_slip_sys), B2(num_slip_sys),
     & S(num_slip_sys), gamma_dot(num_slip_sys)


       do i = 1,num_slip_sys
         tau_eff = DABS(tau(i)-B1(i)-B2(i))
         sgn = SIGN(1.D0, tau(i)-B1(i)-B2(i))
         gamma_dot(i) = sgn*g0*((tau_eff/S(i))**(1.D0/am))
         if (isnan(gamma_dot(i))) then
            write(*,*) tau(i), tau_eff, S(i), sgn, am
            write(*,*) 'gamma_dot is a NaN '
            call abort
         end if

       end do
       return

       end
      
c====================================================================
c====================================================================
c      SUBROUTINE: Calculate back stress, B
c      Eq. (9)
c====================================================================
c====================================================================	
       subroutine calc_back(num_slip_sys, k_back, k_dyn, gamma_dot, 
     & B_0, dtime, B_1)
       implicit double precision (a-h,o-z,k)
       real*8    B_0(num_slip_sys),
     &           B_1(num_slip_sys),
     &           dB_dt(num_slip_sys),
     &           gamma_dot(num_slip_sys)

       do i = 1,num_slip_sys
         dB_dt(i) = k_back*gamma_dot(i) 
     &            - k_dyn*B_0(i)*DABS(gamma_dot(i))
         B_1(i) = B_0(i) + dB_dt(i)*dtime
         if (isnan(B_1(i))) then
           write(*,*) 'B_1 is a NaN '
           call abort
         end if
       end do
       
       return

       end
c====================================================================
c====================================================================
c      SUBROUTINE: Calculate threshold stress, S
c      Eq. (14)
c====================================================================
c====================================================================	
       subroutine calc_threshold(num_slip_sys, Aii, h0, threshold_0, 
     & s_sat, ah, gamma_dot, dtime, threshold_t)
       implicit double precision (a-h,o-z,k)
       real*8    threshold_0(num_slip_sys),
     &           Aii(num_slip_sys, num_slip_sys),
     &           gamma_dot(num_slip_sys),
     &           threshold_t(num_slip_sys)
       real*8    sfrac, dS

       do i = 1,num_slip_sys
         dS = 0.D0
         do j = 1,num_slip_sys
            sgn = SIGN(1.D0, 1.D0-threshold_0(j)/s_sat)
            sfrac = DABS(1.D0-threshold_0(j)/s_sat)
            dS = dS + DABS(gamma_dot(j))*Aii(i,j)*h0*sgn*(sfrac**ah)
         end do

         threshold_t(i) = threshold_0(i) + dS*dtime
         if (threshold_t(i).lt.0.D0) then
           write(*,*) dS, h0, s_sat, ah, gamma_dot
           write(*,*) 'S is negative '
           call abort
         end if
       end do

       return

       end



c====================================================================
c====================================================================
c      SUBROUTINE: Calculate the SPK2 = C:E_el
c      Eq. (5)
c====================================================================
c====================================================================	

       subroutine calc_Spk2(C, E_el, spk2)
       implicit double precision (a-h,o-z)
       real*8 C(3,3,3,3), E_el(3,3), spk2(3,3)

       do i = 1,3
         do j = 1,3
           spk2(i,j) = 0.D0
           do k = 1,3
             do l = 1,3
               spk2(i,j) = spk2(i,j) + C(i,j,k,l)*E_el(k,l)
               if (isnan(spk2(i,j))) then
                 write(*,*) 'Spk2 is a NaN '
                 call abort
               end if
               if (abs(spk2(i,j))-1.d12 .gt. 0.d0) then
                 write(*,*) 'Spk2 is a too large '
                 write(*,*) E_el(k,l)
                 call abort
               end if
             end do
           end do 
         end do
       end do
       return

       end

c====================================================================
c====================================================================
c      SUBROUTINE: Calculate the L_p
c      Eq. (2)
c====================================================================
c====================================================================	

       subroutine calc_L_p(num_slip_sys, gamma_dot, xs, xm, xL_p)
       implicit double precision (a-h,o-z)

       real*8 xL_p(3,3), xs(3,num_slip_sys), 
     & xm(3,num_slip_sys), gamma_dot(num_slip_sys)

       do i = 1,3
         do j = 1,3
           xL_p(i,j) = 0.D0
           do m = 1,num_slip_sys
             xL_p(i,j) = xL_p(i,j) + gamma_dot(m)*xs(i,m)*xm(j,m)
             if (isnan(xL_p(i,j))) then
               write(*,*) 'xL_p is a NaN '
               call abort
             end if
           end do
         end do
       end do

       return
      end
      

c====================================================================
c====================================================================
c      SUBROUTINE: Calculate the F_p
c      Eq. (1)
c====================================================================
c====================================================================	

       subroutine calc_F_p_subend(xL_p,dtime,F_p0,F_p1)
       
       implicit double precision (a-h,o-z)

       real*8 xL_p(3,3),F_p0(3,3),xL_p_dot_L_p(3,3),
     & exp_L_p_dtime(3,3),F_p1(3,3),c_identity(3,3)
       real*8 temp_var, det

c--------------------------------------------------------------------
c Determine omega for calculation of: exp(xL_p*dtime)
c-------------------------------------------------------------------

       call aa_dot_dot_bb(3,xL_p,xL_p,temp_var) ! Double dot product
       
       omega = dsqrt(0.5*temp_var)*dtime   ! Used in next calculation
      
       if(omega.eq.0.D0) then
         do i = 1,3
           do j = 1,3
             F_p1(i,j) = F_p0(i,j)
           end do
         end do
         return
       end if
           
c--------------------------------------------------------------------
       
       call def_Identity(c_identity) ! Define identity matrix
       call aa_dot_bb(3,xL_p,xL_p,xL_p_dot_L_p) ! Dot product
       
       do i = 1,3
         do j = 1,3
           exp_L_p_dtime(i,j) = c_identity(i,j) + xL_p(i,j)*dtime
         end do
       end do
       call calc_Determinant(exp_L_p_dtime,temp_var)
       temp_var = 1.D0/(temp_var**(1.D0/3.D0))
       do i = 1,3
         do j = 1,3
           exp_L_p_dtime(i,j) = exp_L_p_dtime(i,j)*temp_var
         end do
       end do
       
c--------------------------------------------------------------------
c Determine F_p_n+1 = exp(xL_p*dtime) * F_p_n
c-------------------------------------------------------------------        
       
       call aa_dot_bb(3,exp_L_p_dtime,F_p0,F_p1) ! Dot product
       
       if (minval(F_p1) .lt. -1.d1 .or. maxval(F_p1) .gt. 1.d1) then
          write(*,*) 'F_p1 is too large'
          write(*,*) F_p1
          write(*,*) F_p0
          write(*,*) xL_p
          call abort
      endif
       
       return
       end
      
c====================================================================
c====================================================================
c      SUBROUTINE: Calculate the Green Strain Tensor
c      Eq. (4)
c====================================================================
c====================================================================

       subroutine calc_E_el(F_el,E_el)
       
       implicit double precision (a-h,o-z)

       real*8 F_el(3,3), E_el(3,3), F_el_t(3,3)
      
       call transpose(3,F_el,F_el_t)      ! Transpose
       call aa_dot_bb(3,F_el_t,F_el,E_el) ! Dot tensors
     
       do i = 1,3
         E_el(i,i) = E_el(i,i) - 1.D0 ! Finish calculation
           do j = 1,3
             E_el(i,j) = 0.5D0*E_el(i,j)
             if (E_el(i,j) .gt. 1.d2)then
                 write(*,*) 'E_el is too large'
                 write(*,*) E_el(i,j)
                 write(*,*) F_el
                 call abort
             endif
             
           end do
       end do
            
       return
       end
      
c====================================================================
c====================================================================
c      SUBROUTINE: Calculate the elastic deformation gradient
c      Eq. (0), F=Fe*Fp
c====================================================================
c====================================================================	

       subroutine calc_F_el(F_p,F,F_el_inv,F_el)
       
       implicit double precision (a-h,o-z)
       
       real*8 F_p(3,3),F(3,3),F_el(3,3), 
     & F_el_inv(3,3),F_p_inv(3,3)
       
      call inverse_3x3(F_p,F_p_inv)       ! Perform calculations
      call aa_dot_bb(3,F,F_p_inv,F_el)
      call inverse_3x3(F_el,F_el_inv)
      
      if (minval(F_el) .lt. -5.d1 .or. maxval(F_el) .gt. 5.d1) then
          write(*,*) 'F_el is too large'
          write(*,*) F_el
          write(*,*) F_p
          call abort
      endif
      

      return 
      end
 
       
c====================================================================
c====================================================================
c      SUBROUTINE: Calculate Cauchy Stress
c      Eq. (7)
c====================================================================
c====================================================================

       subroutine calc_Sig(F_el,spk2,sig)
       
       implicit double precision (a-h,o-z)
       
       real*8 F_el(3,3),spk2(3,3),sig(3,3), 
     & F_el_t(3,3), a1(3,3)
       
       call aa_dot_bb(3,F_el,spk2,a1)
       call transpose(3,F_el,F_el_t)
       call aa_dot_bb(3,a1,F_el_t,sig)
       
       call calc_Determinant(F_el,det)


       do i = 1,3
         do j = 1,3
           sig(i,j) = sig(i,j)/det
         end do
       end do
       
       return 
       end
c====================================================================
c====================================================================
c      SUBROUTINE: Calculate the resolved shear stress
c      Eq. (3)
c====================================================================
c====================================================================

       subroutine calc_Tau(num_slip_sys, tau, sig, xs, xm)
       
       implicit double precision (a-h,o-z)
       
       real*8 tau(num_slip_sys),sig(3,3),
     & xs(3,num_slip_sys), xm(3,num_slip_sys), dv1(3), dv2(3), dv3(3)
     
       do i = 1,num_slip_sys
         
         dv1(1) = xs(1,i)   ! Dummy vectors
         dv1(2) = xs(2,i)
         dv1(3) = xs(3,i)
         
         dv2(1) = xm(1,i)   ! Dummy vectors
         dv2(2) = xm(2,i)
         dv2(3) = xm(3,i)
         
         call a_dot_bb(dv1,sig, dv3)      ! Dot product
         call a_dot_b(dv3,dv2,tau(i))     ! Dot product
         
c         if (abs(tau(i))-1.d12 .gt. 0.d0) then
c             write(*,*) 'tau is too large'
c             write(*,*) sig
c             call abort
c         endif
            
       end do
       
       return
       end
         

c====================================================================
c====================================================================
c      SUBROUTINE: Rotate the Norm direction
c====================================================================
c====================================================================	

       subroutine calc_Rot_Norm(xm0,xm,F_el_inv, num_slip_sys)
       
       implicit double precision (a-h,o-z)
       
       real*8 
     & xm0(3,num_slip_sys),xm(3,num_slip_sys),
     & F_el_inv(3,3), d_vec_1(3),d_vec_2(3), F_el_inv_t(3,3)
     
       call transpose(3,F_el_inv,F_el_inv_t)
       
       do i = 1, num_slip_sys
       
         d_vec_1(1) = xm0(1,i)
         d_vec_1(2) = xm0(2,i)
         d_vec_1(3) = xm0(3,i)
         
         call calc_Rot_Vec(d_vec_1,d_vec_2,F_el_inv_t)
         R2   = d_vec_2(1)*d_vec_2(1) + 
     &          d_vec_2(2)*d_vec_2(2) + 
     &          d_vec_2(3)*d_vec_2(3) 

        
         xm(1,i) = d_vec_2(1)/DSQRT(R2)
         xm(2,i) = d_vec_2(2)/DSQRT(R2)
         xm(3,i) = d_vec_2(3)/DSQRT(R2)
         
       end do
       
       return
       end
         

c====================================================================
c====================================================================
c      SUBROUTINE: Rotate the Slip direction
c====================================================================
c====================================================================	

       subroutine calc_Rot_Slip(xs0,xs,F_el, num_slip_sys)
       
       implicit double precision (a-h,o-z)
       
       real*8 
     & xs0(3,num_slip_sys),xs(3,num_slip_sys),
     & F_el(3,3), d_vec_1(3),d_vec_2(3)
       
       do i = 1, num_slip_sys
       
         d_vec_1(1) = xs0(1,i)
         d_vec_1(2) = xs0(2,i)
         d_vec_1(3) = xs0(3,i)
         
         call calc_Rot_Vec(d_vec_1,d_vec_2,F_el) ! Calculation
         R2   = d_vec_2(1)*d_vec_2(1) + 
     &          d_vec_2(2)*d_vec_2(2) + 
     &          d_vec_2(3)*d_vec_2(3) 

         xs(1,i) = d_vec_2(1)/DSQRT(R2)
         xs(2,i) = d_vec_2(2)/DSQRT(R2)
         xs(3,i) = d_vec_2(3)/DSQRT(R2)
         
       end do
       
       return
       end
         
c====================================================================
c====================================================================
c      SUBROUTINE: Subroutine Rotate a Vector
c====================================================================
c====================================================================

       subroutine calc_Rot_Vec(vector1,vector2,rotation)
       
       implicit double precision (a-h,o-z)

       real*8 vector1(3),vector2(3),rotation(3,3)
       
       do i = 1,3
         sum = 0.0
         do j = 1,3
           sum = sum + rotation(i,j)*vector1(j)
         end do
         vector2(i) = sum
       end do
       
       return
       end
 

c====================================================================
c====================================================================
c      SUBROUTINE: Definition of identity matrix
c====================================================================
c====================================================================

       subroutine def_Identity(c_Identity)

       implicit double precision (a-h,o-z)
       
       real*8 c_Identity(3,3)
       
       do i = 1,3
         do j = 1,3
          c_Identity(i,j) = 0.0
         end do
          c_Identity(i,i) = 1.0
       end do  
       
       return
       end    

 
c====================================================================
c====================================================================
c       SUBROUTINE: Define the slip systems
c====================================================================
c====================================================================	 

        subroutine def_Slip_Sys(slip_direction,
     &  perp_direction,num_slip_sys)
        implicit double precision (a-h,o-z) 
c-------------------------------------------------------------------
c  Dimension arrays 
c-------------------------------------------------------------------

        real*8
     &  slip_direction(3,num_slip_sys),   ! Miller indices of slip plane normals
     &  perp_direction(3,num_slip_sys),   ! Miller indices of slip plane directions
     &  dumb_v1(3),                       ! Dumb vector
     &  dumb_v2(3)                        ! Dumb vector
     
c-------------------------------------------------------------------
c  Assign slip system normals and slip directions for an FCC.
c-------------------------------------------------------------------

c     !   plane      dir
c     /  1, 1, 1,  0, 1,-1 /  
c     /  1, 1, 1, -1, 0, 1 /  
c     /  1, 1, 1,  1,-1, 0 /  
c     /  1,-1,-1,  0,-1, 1 /  
c     /  1,-1,-1, -1, 0,-1 /  
c     /  1,-1,-1,  1, 1, 0 /  
c     / -1,-1, 1,  0,-1,-1 /
c     / -1,-1, 1,  1, 0, 1 /
c     / -1,-1, 1, -1, 1, 0 /
c     / -1, 1,-1,  0, 1, 1 /  
c     / -1, 1,-1,  1, 0,-1 /  
c     / -1, 1,-1, -1,-1, 0 / 


        slip_direction(1,1) =  1.0
        slip_direction(2,1) =  1.0
        slip_direction(3,1) =  1.0
     
        slip_direction(1,2) =  1.0
        slip_direction(2,2) =  1.0
        slip_direction(3,2) =  1.0
        
        slip_direction(1,3) =  1.0
        slip_direction(2,3) =  1.0
        slip_direction(3,3) =  1.0
        
        slip_direction(1,4) =  1.0
        slip_direction(2,4) = -1.0
        slip_direction(3,4) = -1.0
       
        slip_direction(1,5) =  1.0
        slip_direction(2,5) = -1.0
        slip_direction(3,5) = -1.0
        
        slip_direction(1,6) =  1.0
        slip_direction(2,6) = -1.0
        slip_direction(3,6) = -1.0
        
        slip_direction(1,7) = -1.0
        slip_direction(2,7) = -1.0
        slip_direction(3,7) =  1.0

        slip_direction(1,8) = -1.0
        slip_direction(2,8) = -1.0
        slip_direction(3,8) =  1.0
        
        slip_direction(1,9) = -1.0
        slip_direction(2,9) = -1.0
        slip_direction(3,9) =  1.0
        
        slip_direction(1,10)= -1.0
        slip_direction(2,10)=  1.0
        slip_direction(3,10)= -1.0
        
        slip_direction(1,11)= -1.0
        slip_direction(2,11)=  1.0
        slip_direction(3,11)= -1.0
        
        slip_direction(1,12)= -1.0
        slip_direction(2,12)=  1.0
        slip_direction(3,12)= -1.0
      

        perp_direction(1,1) =  0.0               
        perp_direction(2,1) =  1.0
        perp_direction(3,1) = -1.0
                     
        perp_direction(1,2) = -1.0               
        perp_direction(2,2) =  0.0
        perp_direction(3,2) =  1.0              
        
        perp_direction(1,3) =  1.0
        perp_direction(2,3) = -1.0               
        perp_direction(3,3) =  0.0               
        
        perp_direction(1,4) =  0.0               
        perp_direction(2,4) = -1.0               
        perp_direction(3,4) =  1.0
        
        perp_direction(1,5) = -1.0               
        perp_direction(2,5) =  0.0               
        perp_direction(3,5) = -1.0               
        
        perp_direction(1,6) =  1.0
        perp_direction(2,6) =  1.0
        perp_direction(3,6) =  0.0               
        
        perp_direction(1,7) =  0.0               
        perp_direction(2,7) = -1.0
        perp_direction(3,7) = -1.0

        perp_direction(1,8) =  1.0
        perp_direction(2,8) =  0.0               
        perp_direction(3,8) =  1.0               
        
        perp_direction(1,9) = -1.0               
        perp_direction(2,9) =  1.0               
        perp_direction(3,9) =  0.0               
        
        perp_direction(1,10)=  0.0               
        perp_direction(2,10)=  1.0               
        perp_direction(3,10)=  1.0               
        
        perp_direction(1,11)=  1.0
        perp_direction(2,11)=  0.0               
        perp_direction(3,11)= -1.0
        
        perp_direction(1,12)= -1.0               
        perp_direction(2,12)= -1.0 
        perp_direction(3,12)=  0.0 
        
c-------------------------------------------------------------------
c  Normalize Miller Indices
c-------------------------------------------------------------------
        
        do i = 1,num_slip_sys
          dumb_v1(1) = perp_direction(1,i)
          dumb_v1(2) = perp_direction(2,i)
          dumb_v1(3) = perp_direction(3,i)
          
          dumb_v2(1) = slip_direction(1,i)
          dumb_v2(2) = slip_direction(2,i)
          dumb_v2(3) = slip_direction(3,i)
            
          call normalize_vector(dumb_v1)
          call normalize_vector(dumb_v2)
          
          perp_direction(1,i) = dumb_v1(1)
          perp_direction(2,i) = dumb_v1(2)
          perp_direction(3,i) = dumb_v1(3)
                                
          slip_direction(1,i) = dumb_v2(1)
          slip_direction(2,i) = dumb_v2(2)
          slip_direction(3,i) = dumb_v2(3)
          
        end do
        
        
c-------------------------------------------------------------------
c  Check for normality
c-------------------------------------------------------------------
        
        do i = 1,num_slip_sys
            dumb_v1(1) = perp_direction(1,i)
            dumb_v1(2) = perp_direction(2,i)
            dumb_v1(3) = perp_direction(3,i)
            
            dumb_v2(1) = slip_direction(1,i)
            dumb_v2(2) = slip_direction(2,i)
            dumb_v2(3) = slip_direction(3,i)
            
            prod = 0.0
            
            call a_dot_b(dumb_v1,dumb_v2,prod)
              
            if (prod .gt.  1E-5) then
              write(7,*) 'ERROR: slip sys not normal'
            else
C              write(*,*) 'Normal'
            end if
                          
          
        end do
                   
        return
        end


      
c====================================================================
c====================================================================
c       SUBROUTINE: Define the twin systems
c====================================================================
c====================================================================	 

        subroutine def_Twin_Sys(slip_direction,
     &  perp_direction,num_twin_sys)
        implicit double precision (a-h,o-z) 
c-------------------------------------------------------------------
c  Dimension arrays 
c-------------------------------------------------------------------

        real*8
     &  slip_direction(3,num_twin_sys),   ! Miller indices of twin plane normals
     &  perp_direction(3,num_twin_sys),   ! Miller indices of twin plane directions
     &  dumb_v1(3),                       ! Dumb vector
     &  dumb_v2(3)                        ! Dumb vector
     
c-------------------------------------------------------------------
c  Assign twin system normals and twin directions for an FCC.
c-------------------------------------------------------------------

c     !   plane      dir
c     /  1, 1, 1, -2, 1, 1 /  
c     /  1, 1, 1,  1,-2, 1 /  
c     /  1, 1, 1,  1, 1,-2 /  
c     / -1,-1, 1,  2,-1, 1 /  
c     / -1,-1, 1, -1, 2, 1 /  
c     / -1,-1, 1, -1,-1,-2 /  
c     / -1, 1, 1,  2, 1, 1 /
c     / -1, 1, 1, -1,-2, 1 /
c     / -1, 1, 1, -1, 1,-2 /
c     /  1,-1, 1, -2,-1, 1 /  
c     /  1,-1, 1,  1, 2, 1 /  
c     /  1,-1, 1,  1,-1,-2 / 


        slip_direction(1,1) =  1.0
        slip_direction(2,1) =  1.0
        slip_direction(3,1) =  1.0
     
        slip_direction(1,2) =  1.0
        slip_direction(2,2) =  1.0
        slip_direction(3,2) =  1.0
        
        slip_direction(1,3) =  1.0
        slip_direction(2,3) =  1.0
        slip_direction(3,3) =  1.0
        
        slip_direction(1,4) = -1.0
        slip_direction(2,4) = -1.0
        slip_direction(3,4) =  1.0
       
        slip_direction(1,5) = -1.0
        slip_direction(2,5) = -1.0
        slip_direction(3,5) =  1.0
        
        slip_direction(1,6) = -1.0
        slip_direction(2,6) = -1.0
        slip_direction(3,6) =  1.0
        
        slip_direction(1,7) = -1.0
        slip_direction(2,7) =  1.0
        slip_direction(3,7) =  1.0

        slip_direction(1,8) = -1.0
        slip_direction(2,8) =  1.0
        slip_direction(3,8) =  1.0
        
        slip_direction(1,9) = -1.0
        slip_direction(2,9) =  1.0
        slip_direction(3,9) =  1.0
        
        slip_direction(1,10)=  1.0
        slip_direction(2,10)= -1.0
        slip_direction(3,10)=  1.0
        
        slip_direction(1,11)=  1.0
        slip_direction(2,11)= -1.0
        slip_direction(3,11)=  1.0
        
        slip_direction(1,12)=  1.0
        slip_direction(2,12)= -1.0
        slip_direction(3,12)=  1.0
      

        perp_direction(1,1) = -2.0               
        perp_direction(2,1) =  1.0
        perp_direction(3,1) =  1.0
                     
        perp_direction(1,2) =  1.0               
        perp_direction(2,2) = -2.0
        perp_direction(3,2) =  1.0              
        
        perp_direction(1,3) =  1.0
        perp_direction(2,3) =  1.0               
        perp_direction(3,3) = -2.0               
        
        perp_direction(1,4) =  2.0               
        perp_direction(2,4) = -1.0               
        perp_direction(3,4) =  1.0
        
        perp_direction(1,5) = -1.0               
        perp_direction(2,5) =  2.0               
        perp_direction(3,5) =  1.0               
        
        perp_direction(1,6) = -1.0
        perp_direction(2,6) = -1.0
        perp_direction(3,6) = -2.0               
        
        perp_direction(1,7) =  2.0               
        perp_direction(2,7) =  1.0
        perp_direction(3,7) =  1.0

        perp_direction(1,8) = -1.0
        perp_direction(2,8) = -2.0               
        perp_direction(3,8) =  1.0               
        
        perp_direction(1,9) = -1.0               
        perp_direction(2,9) =  1.0               
        perp_direction(3,9) = -2.0               
        
        perp_direction(1,10)= -2.0               
        perp_direction(2,10)= -1.0               
        perp_direction(3,10)=  1.0               
        
        perp_direction(1,11)=  1.0
        perp_direction(2,11)=  2.0               
        perp_direction(3,11)=  1.0
        
        perp_direction(1,12)=  1.0               
        perp_direction(2,12)= -1.0 
        perp_direction(3,12)= -2.0 
        
c-------------------------------------------------------------------
c  Normalize Miller Indices
c-------------------------------------------------------------------
        
        do i = 1,num_twin_sys
          dumb_v1(1) = perp_direction(1,i)
          dumb_v1(2) = perp_direction(2,i)
          dumb_v1(3) = perp_direction(3,i)
          
          dumb_v2(1) = slip_direction(1,i)
          dumb_v2(2) = slip_direction(2,i)
          dumb_v2(3) = slip_direction(3,i)
            
          call normalize_vector(dumb_v1)
          call normalize_vector(dumb_v2)
          
          perp_direction(1,i) = dumb_v1(1)
          perp_direction(2,i) = dumb_v1(2)
          perp_direction(3,i) = dumb_v1(3)
                                
          slip_direction(1,i) = dumb_v2(1)
          slip_direction(2,i) = dumb_v2(2)
          slip_direction(3,i) = dumb_v2(3)
          
        end do
        
        
c-------------------------------------------------------------------
c  Check for normality
c-------------------------------------------------------------------
        
        do i = 1,num_twin_sys
            dumb_v1(1) = perp_direction(1,i)
            dumb_v1(2) = perp_direction(2,i)
            dumb_v1(3) = perp_direction(3,i)
            
            dumb_v2(1) = slip_direction(1,i)
            dumb_v2(2) = slip_direction(2,i)
            dumb_v2(3) = slip_direction(3,i)
            
            prod = 0.0
            
            call a_dot_b(dumb_v1,dumb_v2,prod)
              
            if (prod .gt.  1E-5) then
              write(7,*) 'ERROR: twin sys not normal'
            else
C              write(*,*) 'Normal'
            end if
                          
          
        end do
                   
        return
        end
      

c====================================================================
c====================================================================
c      Calculate the spatial relationship 
c      between slip sys and twin sys
c      0:non-coplanar 
c      1:Coplanar and acute angles 
c      -1:Coplanar and obtuse
c====================================================================
c====================================================================	 

       subroutine calc_affiliation_matrix(xm_xtal, xs_xtal, 
     & xm_xtal_tw, xs_xtal_tw, num_twin_sys, num_slip_sys,
     & affiliation_sl_tw)
       implicit double precision (a-h,o-z)
       real*8 xm_xtal(3,num_slip_sys), xs_xtal(3,num_slip_sys)
       real*8 xm_xtal_tw(3,num_twin_sys), xs_xtal_tw(3,num_twin_sys)
       real*8 affiliation_sl_tw(num_slip_sys,num_twin_sys)
       real*8 temp_m, temp_s, norm_slip, norm_twin
       
       
       do i = 1,num_slip_sys
           do j = 1,num_twin_sys
               affiliation_sl_tw(i,j) = 0.d0
               temp_m = xm_xtal(1,i)*xm_xtal_tw(1,j)
     &                + xm_xtal(2,i)*xm_xtal_tw(2,j)
     &                + xm_xtal(3,i)*xm_xtal_tw(3,j)
               norm_slip = sqrt(xm_xtal(1,i)**2+xm_xtal(2,i)**2
     &                   + xm_xtal(3,i)**2)
               norm_twin = sqrt(xm_xtal_tw(1,i)**2+xm_xtal_tw(2,i)**2
     &                   + xm_xtal_tw(3,i)**2)
               if (abs(1.d0-abs(temp_m)/(norm_slip*norm_twin))
     &             .lt.1.0d-2)then
                   temp_s = xs_xtal(1,i)*xs_xtal_tw(1,j)
     &                    + xs_xtal(2,i)*xs_xtal_tw(2,j)
     &                    + xs_xtal(3,i)*xs_xtal_tw(3,j)
                   if (abs(temp_s) .gt. 1e-2) then
                       affiliation_sl_tw(i,j) = 1.d0
                   else
                       affiliation_sl_tw(i,j) = -1.d0
                   endif
               endif
           enddo
       enddo
                          
       return
       end

      
c====================================================================
c====================================================================
c      Calculate the Interaction coefficient between slip and twin 
c      systems
c====================================================================
c====================================================================	 

       subroutine calc_xi(num_slip_sys, num_twin_sys, xi_alpbet_cop,
     & xi_alpbet_noncop, affiliation_sl_tw, xi_alpbet)
       implicit double precision (a-h,o-z)
       real*8 xi_alpbet_cop, xi_alpbet_noncop
       real*8 xi_alpbet(num_slip_sys,num_twin_sys)
       real*8 affiliation_sl_tw(num_slip_sys,num_twin_sys)

       
       do i = 1,num_slip_sys
           do j = 1,num_twin_sys
               xi_alpbet(i,j) = xi_alpbet_noncop
               if (abs(affiliation_sl_tw(i,j)) .eq. 1.d0) then
                   xi_alpbet(i,j) = xi_alpbet_cop
               endif
           enddo
       enddo
                          
       return
       end

      
            
c====================================================================
c      Subroutine: Calculate the direction cosines 
c====================================================================
c====================================================================	 

       subroutine calc_Dir_Cos(psi,dir_cos)
       
       implicit double precision (a-h,o-z)
       
       real*8 psi(3),dir_cos(3,3)
       
       s1 = dsin(psi(1))     ! Perform calculation
       c1 = dcos(psi(1))
       s2 = dsin(psi(2))
       c2 = dcos(psi(2))
       s3 = dsin(psi(3))
       c3 = dcos(psi(3))

       dir_cos(1,1) = c1*c3 - s1*s3*c2
       dir_cos(1,2) = s1*c3 + c1*s3*c2
       dir_cos(1,3) = s3*s2

       dir_cos(2,1) = -c1*s3 - s1*c3*c2
       dir_cos(2,2) = -s1*s3 + c1*c3*c2
       dir_cos(2,3) = c3*s2

       dir_cos(3,1) = s1*s2
       dir_cos(3,2) = -c1*s2
       dir_cos(3,3) = c2
       
       return
       end     
c====================================================================
c====================================================================
c      Define the elastic tensor
c====================================================================
c====================================================================	 

       subroutine calc_4th_C(C_11, C_12, C_44,del,C)
       
       implicit double precision (a-h,o-z)
       
       real*8 del(3,3), C(3,3,3,3)
       
       do i = 1,3
        do j = 1,3
         do k = 1,3
          do l = 1,3
           C(i,j,k,l) = C_12 * del(i,j) * del(k,l) +
     &      C_44 * (del(i,k)*del(j,l)+del(i,l)*del(k,j))
          end do
         end do
        end do
       end do
       C(1,1,1,1) = C_11
       C(2,2,2,2) = C_11
       C(3,3,3,3) = C_11
       
       return 
       end
 
c====================================================================
c====================================================================
c      Calculate the Equivalent plastic strain inc
c====================================================================
c====================================================================	 

       subroutine calc_dPEEQ(sig, plastic_work_inc ,dPEEQ)
       implicit double precision (a-h,o-z)
       real*8    sig(3,3), sig_dev(3,3)
       real*8    von_Mise, hydro
       hydro = (sig(1,1)+sig(2,2)+sig(3,3))/3.D0
       sig_dev = sig
       sig_dev(1,1) = sig(1,1) - hydro
       sig_dev(2,2) = sig(2,2) - hydro
       sig_dev(3,3) = sig(3,3) - hydro
       call aa_dot_dot_bb(3, sig_dev, sig_dev, von_Mise) 
       von_Mise = DSQRT(von_Mise*1.5D0)
       dPEEQ = plastic_work_inc/von_Mise
       if (plastic_work_inc .eq. 0.D0) then
         dPEEQ = 0.D0
       end if
       return
       end

c====================================================================
c====================================================================
c      SUBROUTINE: Calculate dislocation density evolution 
c====================================================================
c====================================================================	

       subroutine calc_rho(num_slip_sys, gamma_dot, rho_0, 
     & dtime, rho_1, k1, k2, N_alpha, N_max, mfp_slip, burgers,
     & rho_def_0, d_def, P_cut, affiliation_sl_cs)
       implicit double precision (a-h,o-z,k)
       real*8    rho_0(num_slip_sys),
     &           rho_1(num_slip_sys),
     &           gamma_dot(num_slip_sys),
     &           N_alpha(num_slip_sys), N_max,
     &           mfp_slip(num_slip_sys)
       real*8    burgers, rho_def_0, d_def
       real*8    P_cut(num_slip_sys)
       integer   affiliation_sl_cs(num_slip_sys)
       

       do i = 1,num_slip_sys
           j = affiliation_sl_cs(i)
           
           rho_1(i)=rho_0(i)
     &     + dabs(gamma_dot(i))* dtime*
     &     (k1*dsqrt(rho_0(i)) 
     &     - k2*rho_0(i)
     &     + (1.d0-N_alpha(i)/N_max)/(burgers*mfp_slip(i)))
c     &     - 2.d0*sqrt(6.d0)/3.d0*d_def*rho_def_0*
c     &     (dabs(gamma_dot(i))*(1.d0-P_cut(i)) 
c     &     - dabs(gamma_dot(j))*(1.d0-P_cut(j)))*dtime 
         
           if (gamma_dot(i) .gt. 1d20) then
               write(*,*) "k:",k1,k2
               write(*,*) "dtime",dtime,dabs(gamma_dot(i))
               write(*,*) rho_1(i),rho_0(i)
               write(*,*) dabs(gamma_dot(i))*k1*dsqrt(rho_0(i))*dtime
               write(*,*) "------"
           endif
           if (isnan(rho_1(i))) then
               write(*,*) gamma_dot(i),rho_0(i),P_cut(i)
               stop 'rho is a NaN '
           endif
           if (rho_1(i) .lt. 0.d0) then
               rho_1(i) = 0.d0
           endif
           
       end do
      
       return
       end

c====================================================================
c====================================================================
c     SUBROUTINE:
c     Calculate the dislocation forest hardening 
c     tau_dis^alpha = mu*b*sqrt(h_dis*sum(K^alpha_beta*rho^alpha))
c
c--------------------------------------------------------------------

       subroutine calc_tau_dis(num_slip_sys, mu, burgers, h_dis,
     &     omega1, omega2, rho_0, tau_dis)
       
       implicit double precision (a-h,o-z)
       
       real*8 rho_0(num_slip_sys), tau_dis(num_slip_sys)
       real*8 mu, burgers, h_dis, omega1, omega2
       real*8 term_var(num_slip_sys), del(num_slip_sys,num_slip_sys)
 
       call def_Kron_Del_higher(del,num_slip_sys)
      
       do i = 1,num_slip_sys
           term_var(i) = 0
           do j = 1,num_slip_sys
               term_var(i) = term_var(i) + 
     &         (omega1 + (1 - omega2)*del(i, j))*rho_0(j)
           end do
           tau_dis(i) = mu*burgers*dsqrt(h_dis*term_var(i))
       end do
       

       return
       end

c====================================================================
c====================================================================
c     SUBROUTINE:
c     Calculate the Hall¨CPetch hardening 
c     tau_hp = k_hp*1/sqrt(d)
c
c--------------------------------------------------------------------

       subroutine calc_tau_hp(k_hp, d_grain, tau_hp)
       
       implicit double precision (a-h,o-z)
       
       real*8 tau_hp
       real*8 k_hp, d_grain
       
       
       tau_hp = k_hp/dsqrt(d_grain)
       
      
       return
       end
      
c====================================================================
c====================================================================
c     SUBROUTINE:
c     Calculate the irradiation defects hardening 
c     tau_def^alpha = mu*b*sqrt(h_def*sum(L^alpha_beta*N_def*d_def))
c
c--------------------------------------------------------------------

       subroutine calc_tau_def(num_def_sys, num_slip_sys, mu, burgers,
     & h_def, omega3, omega4, rho_def_0, d_def, tau_def)
       
       implicit double precision (a-h,o-z)
       
       integer index
       real*8 tau_def(num_slip_sys)
       real*8 mu, burgers, h_def, omega3, omega4
       real*8 rho_def_0, d_def
       real*8 term_var(num_slip_sys), del(num_def_sys,num_def_sys)
       
       call def_Kron_Del_higher(del,num_def_sys)
       
       do i = 1,num_slip_sys
           term_var(i) = 0
c           index = ceiling(i/3)
           index = int((i + 2) / 3) 
           do j = 1,num_def_sys
               term_var(i) = term_var(i) + 
     &         (omega3 + (1 - omega4)*del(index, j))*
     &         rho_def_0/4.d0*d_def 
           end do
           tau_def(i) = mu*burgers*dsqrt(h_def*term_var(i))
       end do

       
       return
       end
      
c====================================================================
c====================================================================
c     SUBROUTINE:
c     Calculate the resistance caused by the dislocation pile-ups
c     at the GBs and TBs 
c     tau_TBs^alpha = mu*b/Lambda*N^alpha(1-N^alpha/N^max)
c     It should be noted that if it is a single crystal, 
c     then the grain size should be set to infinity
c--------------------------------------------------------------------

       subroutine calc_tau_TBs(num_slip_sys, num_twin_sys, mu, burgers,
     & d_grain, xi_alpbet, f_beta, t_tw, N_alpha, N_max, tau_TBs, 
     & mfp_slip)
       
       implicit double precision (a-h,o-z)
       real*8 tau_TBs(num_slip_sys)
       real*8 mu, burgers, d_grain
       real*8 xi_alpbet(num_slip_sys,num_twin_sys), f_beta(num_twin_sys)
       real*8 t_tw, N_alpha(num_slip_sys), N_max
       real*8 mfp_slip(num_slip_sys)
       
       call calc_meanfreepath(num_slip_sys,num_twin_sys, d_grain, 
     & xi_alpbet, f_beta, t_tw, mfp_slip)
       
       do i = 1,num_slip_sys
           tau_TBs(i) = mu*burgers/mfp_slip(i)
     &       *N_alpha(i)
       end do
       
       return
       end
      
       
c====================================================================
c====================================================================
c     SUBROUTINE:
c     Calculate the irradiation defects hardening 
c     tau_c^alpha = tau_0+tau_hp+tau_dis^alpha+tau_def^alpha
c
c--------------------------------------------------------------------

       subroutine calc_crss(num_slip_sys, tau_0, tau_hp, tau_dis, 
     & tau_def, tau_TBs, tau_c)
       
       implicit double precision (a-h,o-z)
       
       real*8 tau_0, tau_hp
       real*8 tau_dis(num_slip_sys)
       real*8 tau_def(num_slip_sys)
       real*8 tau_TBs(num_slip_sys)
       real*8 tau_c(num_slip_sys)
       
       do i = 1,num_slip_sys
           tau_c(i) = tau_0 
     &     + tau_hp
     &     + tau_dis(i) 
     &     + tau_def(i)
c     &     + tau_TBs(i)
       end do
       
       return
       end
      
      
      
c====================================================================
c====================================================================
c      SUBROUTINE:  
c      Calculate the function f = [shear rate expression]
c      a new function
c====================================================================
c====================================================================	
       subroutine fnew(num_slip_sys, g0, tau, tau_TBs, tau_c, am,
     & gamma_dot)
       
       implicit double precision (a-h,o-z,k)
       
       real*8 tau(num_slip_sys), tau_c(num_slip_sys),
     & gamma_dot(num_slip_sys), tau_TBs(num_slip_sys)
       real*8 temp_ratio, tau_eff
       
       do i = 1,num_slip_sys
           tau_eff = dabs(tau(i)) - tau_TBs(i)
           ! if back stress is large than RSS
           if (tau_eff .lt. 0.d0) then
               tau_eff = 0.d0
           endif           
C		 tau_crss = 
         if (dabs(tau_eff/tau_c(i)) .lt. 1.2d0) then
               temp_ratio = dabs(tau_eff/tau_c(i))
c               write(*,*) "temp_ratio", temp_ratio
           else
             temp_ratio = 1.2d0
c			 write(80,*) "Time step too large!",temp_ratio
c              call abort
         endif
         sgn = sign(1.D0, tau(i))
         gamma_dot(i) = sgn*g0*(temp_ratio**(1.D0/am))
         
c         if(dabs(tau_eff) .lt. tau_c(i)) then
c             gamma_dot(i) = 0.d0
c         endif
         
         if (isnan(gamma_dot(i))) then
           write(*,*) tau(i), sgn, am
           write(*,*) temp_ratio, tau_c(i), tau_eff       
           call exit(26)
c            stop 'gamma_dot is a NaN '
         endif

       end do
       return

       end      
      
c====================================================================
c====================================================================
c      SUBROUTINE: Calculate defect density evolution 
c====================================================================
c====================================================================	

       subroutine calc_rho_def(num_slip_sys, gamma_dot, burgers,
     & d_def, rho_def_0, dtime, rho_def_1, P_cut, tau, P_anni)
       implicit double precision (a-h,o-z,k)
       real*8    rho_def_0,
     &           rho_def_1,
     &           gamma_dot(num_slip_sys),
     &           temp_var,  
     &           P_cut(num_slip_sys),  
     &           tau(num_slip_sys),
     &           P_anni 
       temp_var = 0.d0
       
       do i = 1,num_slip_sys
           temp_var = temp_var + dabs(gamma_dot(i))
     &     *(0.5d0*P_cut(i)+0.5d0)*dtime 
       enddo
       
       rho_def_1 = rho_def_0 - 2.d0*sqrt(6.d0)/(3.d0*burgers)
     & *rho_def_0*d_def*P_anni*temp_var  
       
       if (rho_def_1 .lt. 0.d0) then
           rho_def_1 = 0.d0
       endif

         if (isnan(rho_def_1)) then
            write(*,*) gamma_dot
            stop 'rho_def is a NaN '
         endif
       
       return

       end
      
c====================================================================
c====================================================================
c      SUBROUTINE: Calculate dislocation piled-up before the GBs and
c      TBs evolution 
c====================================================================
c====================================================================	

       subroutine calc_rho_piledup(num_slip_sys, N_alpha_0, 
     & l_avedistance, burgers, N_max, gamma_dot, dtime, N_alpha_t)
       
       implicit double precision (a-h,o-z,k)
       
       real*8 N_alpha_0(num_slip_sys), l_avedistance
       real*8 burgers, N_max, gamma_dot(num_slip_sys)
       real*8 dtime, N_alpha_t(num_slip_sys)
       
       do i = 1,num_slip_sys
           N_alpha_t(i) = N_alpha_0(i) + l_avedistance/burgers
     &                  *(1.d0-N_alpha_0(i)/N_max)*abs(gamma_dot(i))
     &                  *dtime
           if (N_alpha_t(i) .ge. N_max) then
               N_alpha_t(i) = N_max
           endif
       enddo
       return

       end
      
            
c====================================================================
c====================================================================
c      SUBROUTINE: Calculate the L_p
c      slip+twin
c====================================================================
c====================================================================	

       subroutine calc_L_p_new(num_slip_sys, gamma_dot, xs, xm, xL_p,
     & num_twin_sys, gamma_dot_tw, xs_tw, xm_tw, f_beta )
       implicit double precision (a-h,o-z)

       real*8 xL_p(3,3), xs(3,num_slip_sys), 
     & xm(3,num_slip_sys), gamma_dot(num_slip_sys),
     & xs_tw(3,num_twin_sys), xm_tw(3,num_twin_sys),
     & gamma_dot_tw(num_twin_sys), f_beta(num_twin_sys),
     & ft 
       
       ft = 0.d0
       
       do i = 1,num_twin_sys
           ft = ft + f_beta(i)
       enddo
       
       do i = 1,3
         do j = 1,3
           xL_p(i,j) = 0.D0
           do m = 1,num_slip_sys
             xL_p(i,j) = xL_p(i,j) 
     &       + (1.d0-ft)*gamma_dot(m)*xs(i,m)*xm(j,m)
           end do
           do m = 1,num_twin_sys
             xL_p(i,j) = xL_p(i,j) 
     &       + gamma_dot_tw(m)*xs_tw(i,m)*xm_tw(j,m)
             if (isnan(xL_p(i,j))) then
               write(*,*) 'xL_p is a NaN '
               write(*,*) m, gamma_dot_tw(m)
               call abort
             endif
             if (abs(xL_p(i,j)) .gt. 1.d6 ) then
               write(*,*) 'L_p is too large'
               write(*,*) xL_p(i,j)
               write(*,*) gamma_dot_tw(m)
               write(*,*) xs_tw(i,m), xm_tw(j,m)
               call abort
             endif
           end do
         end do
       end do

       return
       end
      
      
c====================================================================
c====================================================================
c     SUBROUTINE:
c     Calculate the mean free path
c     It should be noted that if it is a single crystal, 
c     then the grain size should be set to infinity
c--------------------------------------------------------------------

       subroutine calc_meanfreepath(num_slip_sys,num_twin_sys, d_grain, 
     & xi_alpbet, f_beta, t_tw, mfp)
       
       implicit double precision (a-h,o-z)
       
       real*8 d_grain
       real*8 xi_alpbet(num_slip_sys,num_twin_sys), f_beta(num_twin_sys)
       real*8 t_tw
       real*8 mfp(num_slip_sys), temp_var(num_slip_sys),ft
       real*8 temp_var2(num_slip_sys)
       
       ft = 0.d0
       
       do i = 1,num_twin_sys
           ft = ft + f_beta(i)
       enddo
       do i = 1,num_slip_sys
           temp_var(i) = 0.d0
           do j = 1,num_twin_sys
               temp_var(i) = temp_var(i)
     &         + xi_alpbet(i,j)*f_beta(j)/(t_tw*(1.d0-ft))
           enddo
       end do
       do i = 1,num_slip_sys
           temp_var2(i) = 1.d0/temp_var(i)
       enddo
       do i = 1,num_slip_sys
           mfp(i) = 1.d0/(1.d0/d_grain+1.d0/temp_var2(i))
       enddo
       
       return
       end
      

c====================================================================
c====================================================================
c     SUBROUTINE:
c     Calculate the the volume of a newly generated twin
c     V^beta = pi/4*mfp_twin^beta^2*t
c--------------------------------------------------------------------

       subroutine calc_newtw_volume(num_twin_sys, d_grain, 
     & xi_betbet, f_beta, t_tw, mfp_twin, volume_beta)
       
       implicit double precision (a-h,o-z)
       
       real*8 d_grain
       real*8 xi_betbet(num_twin_sys,num_twin_sys), f_beta(num_twin_sys)
       real*8 t_tw
       real*8 mfp_twin(num_twin_sys), volume_beta(num_twin_sys), pi

       PI = 4.D0*DATAN(1.D0)
       
       call calc_meanfreepath(num_twin_sys,num_twin_sys, d_grain, 
     & xi_betbet, f_beta, t_tw, mfp_twin)

       do i = 1,num_twin_sys
           volume_beta(i) = pi/4.d0*mfp_twin(i)**2*t_tw
       enddo
       
       return
       end
   
      
c====================================================================
c====================================================================
c     SUBROUTINE:
c     Calculate the stress needed to bring two partial dislocations 
c     from the equilibrium distance x0 to a critical distance xc 
c     in which the twin nucleation can be facilitated
c     tau_r = mu*b/(2pi)*(1/(x0+xc)+1/2x0)
c--------------------------------------------------------------------

       subroutine calc_tau_r_tw(mu, burgers, nu, x_c,
     & gamma_sf, tau_r_tw)
       
       implicit double precision (a-h,o-z)
       
       real*8 mu, burgers, nu, x_c
       real*8 gamma_sf, tau_r_tw
       real*8 pi, x_0
       
       PI = 4.D0*DATAN(1.D0)
c       x_0 = mu/gamma_sf*(burgers**2)/(8.d0*pi)*(2.d0+nu)/(1.d0-nu)
c       tau_r_tw = mu*burgers/(2.d0*pi)*(1.d0/(x_0+x_c)+1.d0/(2.d0*x_0))
       x_0 = mu/gamma_sf*(burgers**2)/(24.d0*pi)*(2.d0+nu)/(1.d0-nu)
       tau_r_tw = mu*burgers/(2.d0*pi)*(1.d0/(x_0+x_c)+1.d0/(2.d0*x_0))
c       write(*,*) gamma_sf,x_0,tau_r_tw*1d-6
       return
       end
      
  
c====================================================================
c====================================================================
c     SUBROUTINE:
c     Calculate the critical stress for twin formation
c     tau_f_tw = gammasf/bp+Gbp/dgrain
c--------------------------------------------------------------------

       subroutine calc_tau_formation_tw(mu, burgers, gamma_sf,
     & d_grain, tau_f_tw)
       
       implicit double precision (a-h,o-z)
       
       real*8 mu, burgers, d_grain
       real*8 gamma_sf, tau_f_tw
 
       tau_f_tw = gamma_sf/(burgers/sqrt(3.d0))
     & + mu*(burgers/sqrt(3.d0))/d_grain 
       
       return
       end
      
c====================================================================
c====================================================================
c     SUBROUTINE:
c     Calculate the number density of potential twin nuclei per 
c     unit time
c--------------------------------------------------------------------

       subroutine calc_nucrhorate_tw(num_slip_sys, num_twin_sys, L_tw,
     & affiliation_sl_tw, gamma_dot, rho, tau, tau_TBs, tau_c,
     & nucrhorate_tw, m_tw, totaltime)
       
       implicit double precision (a-h,o-z)
       
       real*8 L_tw, affiliation_sl_tw(num_slip_sys,num_twin_sys)
       real*8 gamma_dot(num_slip_sys), rho(num_slip_sys)
       real*8 tau(num_slip_sys), tau_TBs(num_slip_sys)
       real*8 tau_c(num_slip_sys)
       real*8 nucrhorate_tw(num_twin_sys)
       real*8 gamma_dot_sum, rho_sum, m_tw, temp_flag(num_twin_sys)
       real*8 tau_eff2, tau_eff3
       integer temp_flag2, temp_flag3
       
       do i = 1,num_twin_sys
           temp_flag(i) = 0
           temp_flag2 = 0
           temp_flag3 = 0
           do j = 1,num_slip_sys
               if (affiliation_sl_tw(j,i) .eq. 1.d0) then
                   if (temp_flag2 .gt. 0)then
                       temp_flag3 = j
                   else
                       temp_flag2 = j
                   endif
               endif
           enddo
c           if (abs(tau(temp_flag2))-tau_c(temp_flag2).ge.0.d0
c     &     .and.abs(tau(temp_flag3))-tau_c(temp_flag3).ge.0.d0 )then
c               temp_flag(i) = 1.d0
c           endif
           tau_eff2 = dabs(tau(temp_flag2)) - tau_TBs(temp_flag2)
           tau_eff3 = dabs(tau(temp_flag3)) - tau_TBs(temp_flag3)
           
           if (tau_eff2 .lt. 0.d0) tau_eff2 = 0.d0
           if (tau_eff3 .lt. 0.d0) tau_eff3 = 0.d0
           
           temp_flag(i) = ( (tau_eff2 * tau_eff3) / 
     &                    (tau_c(temp_flag2) * tau_c(temp_flag3)) )
     &                    **(1.d0/m_tw)
           if (isnan(temp_flag(i)) .or.
     &       abs(temp_flag(i)).gt.huge(1.d0)) then
               write(*,*) 'temp_flag is a NaN or Inf'
               write(*,*) i, temp_flag(i)
               write(*,*) temp_flag2, temp_flag3
               write(*,*) abs(tau(temp_flag2)*tau(temp_flag3))
     &     /(tau_c(temp_flag2)*tau_c(temp_flag3))
               write(*,*) tau(temp_flag2), tau(temp_flag3)
               write(*,*) tau_c(temp_flag2), tau_c(temp_flag3)
               call abort
            end if
            if (abs(temp_flag(i)) .gt. 1.d5 .and. 
     &           totaltime .gt. 0.d0) then
               write(*,*) 'temp_flag is too large'
               write(*,*) i, temp_flag(i), m_tw
               write(*,*) tau(temp_flag2), tau(temp_flag3)
               write(*,*) tau_c(temp_flag2), tau_c(temp_flag3)
               write(*,*) abs(tau(temp_flag2)/tau_c(temp_flag2))
               write(*,*) abs(tau(temp_flag3)/tau_c(temp_flag3))
               write(*,*) temp_flag2, temp_flag3
               call abort
            end if
            if (abs(temp_flag(i)) .gt. 1.2d0) then
               temp_flag(i) = 1.2d0
            end if
       enddo
             
       gamma_dot_sum = abs(gamma_dot(temp_flag2))
     &               + abs(gamma_dot(temp_flag3))
       rho_sum = rho(temp_flag2) + rho(temp_flag3)
       
       
       do i = 1,num_twin_sys
           nucrhorate_tw(i) = temp_flag(i)*gamma_dot_sum
     &                      *rho_sum/(3.d0*L_tw)
           if (nucrhorate_tw(i) .lt. 0.d0)then
               write(*,*) 'rate of twin nuc less than 0'
               write(*,*) nucrhorate_tw(i), temp_flag(i)
               write(*,*) gamma_dot_sum, rho_sum
               write(*,*) L_tw
               call abort
           endif
           if (isnan(nucrhorate_tw(i)) .or.
     &       abs(nucrhorate_tw(i)).gt.huge(1.d0)) then
               write(*,*) 'nucrhorate_tw is a NaN or Inf'
               write(*,*) i, nucrhorate_tw(i), temp_flag(i)
               write(*,*) gamma_dot_sum, rho_sum, L_tw
               call abort
           end if
           if (abs(nucrhorate_tw(i)) .gt. 1.d30) then
               write(*,*) 'nucrhorate_tw is too large'
               write(*,*) i, nucrhorate_tw(i), temp_flag(i)
               write(*,*) gamma_dot_sum, rho_sum, L_tw
               call abort
           end if
             
       enddo
       
       return
       end

c====================================================================
c====================================================================
c     SUBROUTINE:
c     Calculate the change rate of twin fraction
c     f_beta_dot = (1-ft)*V_beta*N_t_beta_dot
c--------------------------------------------------------------------

       subroutine calc_rate_ftw(num_twin_sys, volume_beta, tau_r_tw,
     & tau_f_tw, nucrhorate_tw, f_beta, Vcs, B_k, T, A_tw, f_beta_max,
     & tau_tw, f_beta_dot)
       
       implicit double precision (a-h,o-z)
       
       real*8 volume_beta(num_twin_sys), tau_r_tw
       real*8 tau_f_tw, nucrhorate_tw(num_twin_sys)
       real*8 f_beta(num_twin_sys), Vcs, B_k, T, A_tw, f_beta_max
       real*8 tau_tw(num_twin_sys), f_beta_dot(num_twin_sys)
       real*8 pncs, ptw, ft
       
       ft = 0.d0
       
       do i = 1,num_twin_sys
           ft = ft + f_beta(i)
       enddo
       
       do i = 1,num_twin_sys
           pncs = 1.d0-exp(-Vcs/(B_k*T)*(tau_r_tw-abs(tau_tw(i))))
           if (pncs .lt. 0.d0)then
               pncs = 0.d0
           endif
           if(tau_tw(i).gt.0.d0) then
               ptw = exp(-(tau_f_tw/abs(tau_tw(i)))**A_tw)
           else
               ptw = 0.d0
           endif
           f_beta_dot(i) = (1.d0-ft)*volume_beta(i)*nucrhorate_tw(i)
     &                    *pncs*ptw
           
           if(f_beta_dot(i) .lt. 0.d0) then
               write(*,*) 'the change rate of twin fraction less than 0'
               write(*,*) ft, volume_beta(i), nucrhorate_tw(i),pncs,ptw
               write(*,*) tau_r_tw, tau_tw(i)
               call abort
           endif
           if (isnan(f_beta_dot(i))) then
               write(*,*) 'f_beta_dot is a NaN '
               write(*,*) i, volume_beta(i), nucrhorate_tw(i)
               write(*,*) ft, pncs, ptw
               call abort
           end if
           if (abs(f_beta_dot(i)) .gt. 1.d6 ) then
               write(*,*) 'f_beta_dot is too large'
               write(*,*) f_beta_dot(i)
               write(*,*) ft
               write(*,*) volume_beta(i), pncs, ptw
               write(*,*) nucrhorate_tw(i)
               call abort
           endif
           if(f_beta(i) .gt. f_beta_max) then
               f_beta_dot(i) = 0.d0
           endif
           
      
       enddo
       
       
       return
       end

      
c====================================================================
c====================================================================
c     SUBROUTINE:
c     Calculate the shear rate contributed by a twin system
c     gamma_beta_dot = gamma_tw*f_beta_dot
c--------------------------------------------------------------------

       subroutine fnew_tw(num_twin_sys, gamma_tw, f_beta_dot, tau_tw,
     & gamma_dot_tw)
       
       implicit double precision (a-h,o-z)
       
       real*8 gamma_tw, f_beta_dot(num_twin_sys)
       real*8 tau_tw(num_twin_sys) 
       real*8 gamma_dot_tw(num_twin_sys)
  
       do i = 1,num_twin_sys
           gamma_dot_tw(i) = gamma_tw*f_beta_dot(i)
           if (isnan(gamma_dot_tw(i))) then
               write(*,*) 'gamma_dot_tw is a NaN '
               write(*,*) i, f_beta_dot(i)
               call abort
           end if
           if (abs(gamma_dot_tw(i)) .gt. 1.d6 ) then
               write(*,*) 'gamma_dot_tw is too large'
               write(*,*) gamma_dot_tw(i)
               write(*,*) gamma_tw
               write(*,*) f_beta_dot(i)
               call abort
             endif
       enddo
       
       return
       end
       
c====================================================================
c====================================================================
c     SUBROUTINE:
c     Calculate the twin volume fraction evolution 
c--------------------------------------------------------------------

       subroutine calc_ftw(num_twin_sys, f_beta_0, f_beta_dot,
     & dtime, f_beta_t)
       
       implicit double precision (a-h,o-z)
       
       real*8 f_beta_0(num_twin_sys)
       real*8 f_beta_dot(num_twin_sys), dtime
       real*8 f_beta_t(num_twin_sys)
  
       do i = 1,num_twin_sys
           f_beta_t(i) = f_beta_0(i)+f_beta_dot(i)*dtime
       enddo
       
       return
       end

      
c====================================================================
c====================================================================
c     SUBROUTINE:
c     Calculate the probability of dislocation cutting through 
c     irradiation defects. 
c--------------------------------------------------------------------

       subroutine calc_P_cut(num_slip_sys, xs_xtal, xs0, xm0, spk2,
     & tau, gamma_sf, mu, burgers, A_cs, omega_cs1, omega_cs2,
     & omega_cs3, p_cs, q_cs, rho_def, d_def, A_cut, alpha_cut, k_cut, 
     & m_cut, B_k, T, P_cut, affiliation_sl_cs, tauce, tauge, taucs,
     & E_cs_positive, E_cs_negative, E_cut)
       
       implicit double precision (a-h,o-z)
       
       real*8 P_cut(num_slip_sys), E_cut(num_slip_sys)
       real*8 E_cs_positive(num_slip_sys), E_cs_negative(num_slip_sys)
       
       real*8 xs_xtal(3,num_slip_sys)
       integer affiliation_sl_cs(num_slip_sys)
       
       real*8 xs0(3,num_slip_sys), xm0(3,num_slip_sys), spk2(3,3)
       real*8 tau(num_slip_sys), tauce(num_slip_sys)
       real*8 tauge(num_slip_sys), taucs(num_slip_sys)
       
       real*8 gamma_sf, mu, burgers, A_cs, omega_cs1, omega_cs2
       real*8 omega_cs3, p_cs, q_cs
       
       real*8 rho_def, d_def, A_cut, alpha_cut, k_cut, m_cut
       
       real*8 B_k, T

       
       call calc_rela_gc(num_slip_sys, xs_xtal, affiliation_sl_cs)
       call calc_four_tau_rss(num_slip_sys, xs0, xm0, spk2,
     & affiliation_sl_cs, tauce, tauge, taucs)
       call calc_E_cs(num_slip_sys, abs(tau), abs(tauce), abs(tauge), 
     & abs(taucs), gamma_sf, mu, burgers, A_cs, omega_cs1, omega_cs2, 
     & omega_cs3, p_cs, q_cs, E_cs_positive)
       call calc_E_cs(num_slip_sys, abs(tau), -abs(tauce), -abs(tauge), 
     & abs(taucs), gamma_sf, mu, burgers, A_cs, omega_cs1, omega_cs2, 
     & omega_cs3, p_cs, q_cs, E_cs_negative)
       call calc_E_cut(num_slip_sys, tau, rho_def, d_def,
     & mu, burgers, A_cut, alpha_cut, k_cut, m_cut, E_cut, totaltime)
       
       do i = 1,num_slip_sys
c           P_cut(i) = 1.d0
           P_cut(i) = 5.d-1*(
     &     1.d0/(1.d0+exp((E_cut(i)-E_cs_positive(i))/(B_k*T)))
     &     + 1.d0/(1.d0+exp((E_cut(i)-E_cs_negative(i))/(B_k*T))))
           if (P_cut(i).gt.1.d0)then
               write(*,*) "The possiblity greater than 1"
               call abort
           endif
           if (isnan(P_cut(i))) then
               write(*,*) 'P_cut is a NaN '
               write(*,*) E_cut(i),E_cs_positive(i),E_cs_negative(i)
               write(*,*)  abs(tau(i)), -abs(tauce(i)), -abs(tauge(i)),
     &         abs(taucs(i)), mu, burgers
               write(*,*) 'spk2:'
               write(*,*) spk2
               call abort
           endif
       enddo
       
       return
       end      
      
      
c====================================================================
c====================================================================
c     SUBROUTINE:
c     Calculate the relationship between the original glide plane 
c     and the cross-slip plane.
c     affiliation_sl_cs(alpha)=eta means that original glide plane
c     is alpha,  the cross-slip plane is eta
c--------------------------------------------------------------------

       subroutine calc_rela_gc(num_slip_sys, xs_xtal, affiliation_sl_cs)
       
       implicit double precision (a-h,o-z)
       
       real*8 xs_xtal(3,num_slip_sys)
       integer affiliation_sl_cs(num_slip_sys)
       real*8 dot_prod, magnitude1, magnitude2
       real*8 tolerance
       
       tolerance = 1.d-6
       
       do i = 1, num_slip_sys
           affiliation_sl_cs(i) = 0
       enddo
       
       do i = 1, num_slip_sys
           do j = i + 1, num_slip_sys
               ! Compute the dot product of the i-th and j-th vectors
               dot_prod = 0.d0
               do k = 1, 3
                   dot_prod = dot_prod + xs_xtal(k, i) * xs_xtal(k, j)
               end do
        
               ! Compute the magnitudes of the i-th and j-th vectors
               magnitude1 = 0.d0
               magnitude2 = 0.d0
               do k = 1, 3
                   magnitude1 = magnitude1 + xs_xtal(k, i) ** 2
                   magnitude2 = magnitude2 + xs_xtal(k, j) ** 2
               end do
               magnitude1 = sqrt(magnitude1)
               magnitude2 = sqrt(magnitude2)
        
               ! Check if vectors are parallel (dot product is proportional to magnitudes)
               if (abs(dot_prod) .gt. 
     &         (magnitude1 * magnitude2 - tolerance)) then
                   ! If they are parallel, set the relationship in A
                   affiliation_sl_cs(i) = j
                   affiliation_sl_cs(j) = i
               endif
           enddo
       end do

       
       return
       end         

           
c====================================================================
c====================================================================
c      SUBROUTINE: Calculate the resolved shear stress
c      on original glide plane and cross-slip plane
c      Perpendicular to slip dir and Parallel to slip dir
c====================================================================
c====================================================================

       subroutine calc_four_tau_rss(num_slip_sys, xs, xm, sig,
     & affiliation_sl_cs, tauce, tauge, taucs)
       
       implicit double precision (a-h,o-z)
       
       integer affiliation_sl_cs(num_slip_sys)
       real*8 sig(3,3), xs(3,num_slip_sys), 
     & xm(3,num_slip_sys), dv1(3), dv2(3), dv3(3)
       real*8 tauce(num_slip_sys), tauge(num_slip_sys)
       real*8 taucs(num_slip_sys)
       real*8 dv4(3), k
            
       do i = 1,num_slip_sys
         k = affiliation_sl_cs(i)
         dv1(1) = xs(1,k)   ! Dummy vectors
         dv1(2) = xs(2,k)
         dv1(3) = xs(3,k)
         
         dv2(1) = xm(1,k)   ! Dummy vectors
         dv2(2) = xm(2,k)
         dv2(3) = xm(3,k)
         
         call cross_product(dv1(1),dv1(2),dv1(3),
     &   dv2(1),dv2(2),dv2(3),dv4(1),dv4(2),dv4(3))
         call a_dot_bb(dv4,sig, dv3)      ! Dot product
         call a_dot_b(dv3,dv2,tauce(i))     ! Dot product
         
       end do
       do i = 1,num_slip_sys
         
         dv1(1) = xs(1,i)   ! Dummy vectors
         dv1(2) = xs(2,i)
         dv1(3) = xs(3,i)
         
         dv2(1) = xm(1,i)   ! Dummy vectors
         dv2(2) = xm(2,i)
         dv2(3) = xm(3,i)
         call cross_product(dv1(1),dv1(2),dv1(3),
     &   dv2(1),dv2(2),dv2(3),dv4(1),dv4(2),dv4(3))
         call a_dot_bb(dv4,sig, dv3)      ! Dot product
         call a_dot_b(dv3,dv2,tauge(i))     ! Dot product
         
       end do       
       do i = 1,num_slip_sys
         k = affiliation_sl_cs(i)
         dv1(1) = xs(1,k)   ! Dummy vectors
         dv1(2) = xs(2,k)
         dv1(3) = xs(3,k)
         
         dv2(1) = xm(1,k)   ! Dummy vectors
         dv2(2) = xm(2,k)
         dv2(3) = xm(3,k)
         
         call a_dot_bb(dv1,sig, dv3)      ! Dot product
         call a_dot_b(dv3,dv2,taucs(i))     ! Dot product
         
       end do       
       
       return
       end
         
      
c====================================================================
c====================================================================
c     SUBROUTINE:
c     Calculate the energy of dislocation cross-slip through irr def 
c--------------------------------------------------------------------

       subroutine calc_E_cs(num_slip_sys, sig_gs, sig_ce, sig_ge, 
     & sig_cs, gamma_sf, mu, burgers, A_cs, omega_cs1, omega_cs2, 
     & omega_cs3, p_cs, q_cs, E_cs)

       
       implicit double precision (a-h,o-z)
       
       real*8 sig_gs(num_slip_sys), sig_ce(num_slip_sys)
       real*8 sig_ge(num_slip_sys), sig_cs(num_slip_sys)
       real*8 Gs(num_slip_sys), Ce(num_slip_sys)
       real*8 Ge(num_slip_sys), Cs(num_slip_sys)
       real*8 gamma_sf, mu, burgers
       real*8 A_cs, omega_cs1, omega_cs2, omega_cs3, p_cs, q_cs
       real*8 E_cs(num_slip_sys)
       real*8 temp
       
       do i = 1,num_slip_sys
           Gs(i) = sig_gs(i)/mu
           Cs(i) = sig_cs(i)/mu
           Ce(i) = (3.d0*sqrt(2.d0)*gamma_sf + burgers*sig_ce(i))
     &     /(3.d0*sqrt(2.d0)*mu*burgers) 
           Ge(i) = (3.d0*sqrt(2.d0)*gamma_sf + burgers*sig_ge(i))
     &     /(3.d0*sqrt(2.d0)*mu*burgers)
       enddo
       
       do i = 1,num_slip_sys
           temp = Cs(i)+omega_cs1*Ce(i)+omega_cs2*Gs(i)
     &     +omega_cs3*Ge(i)
           if (temp .le. 0.d0)then
               temp = 0.d0
           else
               temp = temp**p_cs
           endif
           
           if (1.d0-temp .lt. 0.d0) then
               write(*,*) 'RSS is too large'
               E_cs(i) = 0.d0
c               call abort
           else
               E_cs(i) = A_cs*mu*burgers**3.d0*(1.d0 - temp)**q_cs
           endif
       enddo
       
       return
       end      
     
       
c====================================================================
c====================================================================
c     SUBROUTINE:
c     Calculate the energy of dislocation cutting through irr def 
c--------------------------------------------------------------------

       subroutine calc_E_cut(num_slip_sys, sig_gs, rho_def, d_def,
     & mu, burgers, A_cut, alpha_cut, k_cut, m_cut, E_cut)

       
       implicit double precision (a-h,o-z)
       
       real*8 sig_gs(num_slip_sys), tau_cut(num_slip_sys)
       real*8 rho_def, d_def, L_def
       real*8 mu, burgers
       real*8 A_cut, alpha_cut, k_cut, m_cut
       real*8 E_cut(num_slip_sys)
       
       L_def = rho_def**(-1.d0/3.d0)*0.3
c      high-density limit: spacing-controlled cutting model no longer valid
c      use constant E_cut = A_cut * mu * b^3
       if (L_def - d_def .le. 0.d0) then
           E_cut_const = mu * burgers**3.d0 * A_cut
           do i = 1, num_slip_sys
               E_cut(i) = E_cut_const
           enddo
           return
       endif
       
       do i = 1,num_slip_sys
           tau_cut(i) = k_cut*mu*burgers/(L_def-d_def)
     &     +m_cut*mu
           E_cut(i) = mu*burgers**3.d0*A_cut
     &     *(1.d0-abs(sig_gs(i))/tau_cut(i))**alpha_cut
           
           if(abs(sig_gs(i)) .ge. tau_cut(i))then
               E_cut(i) = 0.d0
           endif
       enddo
              
       return
       end             
       
       
c====================================================================
c====================================================================
c====================== Utility  Subroutines ========================
c====================================================================
c====================================================================

      
c====================================================================
c====================================================================
c      SUBROUTINE: Define Kronecker Delta tensor
c====================================================================
 

       subroutine def_Kron_Del(del)

       implicit double precision (a-h,o-z)
       
       real*8 del(3,3)
       
       do i = 1,3
         do j = 1,3
          del(i,j) = 0.0
         end do
          del(i,i) = 1.0
       end do  
       
       return
       end      

c====================================================================
c====================================================================
c      SUBROUTINE: Define Kronecker Delta tensor(higher order)
c====================================================================
 

       subroutine def_Kron_Del_higher(del, order)

       implicit double precision (a-h,o-z)
       
       integer order
       real*8 del(order, order)
       
       do i = 1,order
         do j = 1,order
          del(i,j) = 0.d0
         end do
          del(i,i) = 1.d0
       end do  
       
       return
       end      
      
      
c====================================================================
c====================================================================
c      SUBROUTINE: dot product between a and b
c====================================================================
c====================================================================	 

       subroutine a_dot_b(v1,v2,prod)
       
       implicit double precision (a-h,o-z)
       
       real*8 v1(3),v2(3)
       
       prod = 0.0    ! Initialize 
       
       do i = 1,3
         prod = prod + v1(i)*v2(i) ! Perform calculation
       end do
      
       return
       end

c--------------------------------------------------------------------
c  SUBROUTINE: Calculate a vector dot product.
c
c  c = a dot bb
c
c--------------------------------------------------------------------


       subroutine a_dot_bb(a,b,c)
       
       implicit double precision (a-h,o-z)
       
       real*8 a(3),b(3,3),c(3)
       
       do i = 1,3
         c(i) = 0.0  ! Initialize
         do j = 1,3
           c(i) = c(i)+a(j)*b(j,i) ! Perform calculation
         end do
       end do
       
       return
       end
        
c--------------------------------------------------------------------
c  SUBROUTINE: Calculate a vector cross product.
c
c  c = a X b
c
c--------------------------------------------------------------------

      subroutine cross_product(a1,a2,a3,b1,b2,b3,c1,c2,c3)

      implicit double precision (a-h,o-z)

      c1 = a2*b3 - a3*b2
      c2 = a3*b1 - a1*b3
      c3 = a1*b2 - a2*b1

      return
      end

c====================================================================
c====================================================================
c
c  SUBROUTINE: Calculate sin(theta) between two vectors.
c
c====================================================================
c====================================================================      

      subroutine calc_sin(x1,x2,x3,y1,y2,y3,z)

      implicit double precision (a-h,o-z)

      z = sqrt(1 - (x1*y1 + x2*y2 + x3*y3)**2)

      return
      end

c====================================================================
c====================================================================
c
c  SUBROUTINE: Normalize the length of a vector to one.
c
c--------------------------------------------------------------------

      subroutine normalize_vector(v1)

      implicit double precision (a-h,o-z)
      
      Dimension v1(3)

      xlength = sqrt(v1(1)*v1(1)+v1(2)*v1(2)+v1(3)*v1(3))
      v1(1) = v1(1) / xlength
      v1(2) = v1(2) / xlength
      v1(3) = v1(3) / xlength

      return
      end

c====================================================================
c====================================================================
c
c  SUBROUTINE: Transpose a (n x n) tensor.
c
c--------------------------------------------------------------------

      subroutine transpose(n,a,b)

      implicit double precision (a-h,o-z)

      real*8 a(n,n), b(n,n)

      do i = 1,n
        do j = 1,n
          b(i,j) = a(j,i)
        end do
      end do

      return
      end


c====================================================================
c====================================================================
c  SUBROUTINE: 
c  Calculate the dot product of two 2nd rank tensors.
c  Result is stored in cc(i,j)
c
c--------------------------------------------------------------------

      subroutine aa_dot_bb(n,a,b,c)

      implicit double precision (a-h,o-z)

      real*8 a(n,n), b(n,n), c(n,n)

      do i = 1,n
        do j = 1,n
          c(i,j) = 0
          do k = 1,n
            c(i,j) = c(i,j) + a(i,k)*b(k,j)
          end do
        end do
      end do

      return
      end

c====================================================================
c====================================================================
c
c  SUBROUTINE: Calculate the double dot product of two 2nd rank tensors.
c
c--------------------------------------------------------------------

      subroutine aa_dot_dot_bb(n,a,b,sum)

      implicit double precision (a-h,o-z)

      real*8 a(n,n), b(n,n)

      sum = 0.0
      do i = 1,n
        do j = 1,n
          sum = sum + a(i,j)*b(i,j)
        end do
      end do

      return
      end

c====================================================================
c====================================================================
c
c  SUBROUTINE: Calculate the double dot product of two 4th rank tensors.
c  Result is stored in c(i,j,k,l)
c
c--------------------------------------------------------------------

      subroutine aaaa_dot_dot_bbbb(n,a,b,c)

      implicit double precision (a-h,o-z)

      real*8 a(n,n,n,n), b(n,n,n,n), c(n,n,n,n)

      do i = 1,n
       do j = 1,n
        do k = 1,n
         do l = 1,n
          c(i,j,k,l) = 0
          do m1 = 1,n
           do m2 = 1,n
            c(i,j,k,l) = c(i,j,k,l) + a(i,j,m1,m2)*b(m1,m2,k,l)
           end do !m2
          end do !m1
         end do !l
        end do !k
       end do !j
      end do !i

      return
      end

c====================================================================
c====================================================================
c
c  SUBROUTINE: Calculate the double dot product of a 4th rank tensor and
c  a 2nd rank tensor.  Result is stored in c(i,j).
c
c--------------------------------------------------------------------

      subroutine aaaa_dot_dot_bb(n,a,b,c)

      implicit double precision (a-h,o-z)

      real*8 a(n,n,n,n), b(n,n), c(n,n)

      do i = 1,n
       do j = 1,n
        c(i,j) = 0
        do k = 1,n
         do l = 1,n
          c(i,j) = c(i,j) + a(i,j,k,l)*b(k,l)
         end do !l
        end do !k
       end do !j
      end do !i

      return
      end

c====================================================================
c====================================================================
c
c  SUBROUTINE: Calculate the double dot product of a 2nd rank tensor and
c  a 4th rank tensor.  Result is stored in c(i,j).
c
c--------------------------------------------------------------------

      subroutine aa_dot_dot_bbbb(n,a,b,c)

      implicit double precision (a-h,o-z)

      real*8 a(n,n), b(n,n,n,n), c(n,n)

      do i = 1,n
       do j = 1,n
        c(i,j) = 0
        do k = 1,n
         do l = 1,n
          c(i,j) = c(i,j) + a(k,l) * b(k,l,i,j)
         end do !l
        end do !k
       end do !j
      end do !i

      return
      end

c====================================================================
c====================================================================
c
c  SUBROUTINE: Rotates any 3x3x3x3 tensor by a rotation matrix.
c
c  c(i,j,k,l) = a(i,m) * a(j,n) * a(k,p) * a(l,q) * b(m,n,p,q)
c
c--------------------------------------------------------------------

      subroutine rotate_4th(a,b,c)

      implicit double precision (a-h,r-z)

      real*8 a(3,3), b(3,3,3,3), c(3,3,3,3)
      do i = 1,3
        do j = 1,3
          do k = 1,3
            do l = 1,3
              c(i,j,k,l) = 0.D0
            end do
          end do
        end do
      end do

      do i = 1,3
        do j = 1,3
          do k = 1,3
            do l = 1,3
              do m = 1,3
                do n = 1,3
                  do p = 1,3
                    do q = 1,3
                      c(i,j,k,l) = c(i,j,k,l) + 
     &                a(i,m)*a(j,n)*a(k,p)*a(l,q)*b(m,n,p,q)

                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
      return
      end

c====================================================================
c====================================================================
c
c  SUBROUTINE: Calculate the inverse of a 3 x 3 matrix.
c
c--------------------------------------------------------------------

      subroutine inverse_3x3(a,b)

      implicit double precision (a-h,o-z)

      real*8 a(3,3), b(3,3)

      b(1,1) = a(2,2) * a(3,3) - a(3,2) * a(2,3)
      b(1,2) = a(3,2) * a(1,3) - a(1,2) * a(3,3)
      b(1,3) = a(1,2) * a(2,3) - a(2,2) * a(1,3)
      b(2,1) = a(3,1) * a(2,3) - a(2,1) * a(3,3)
      b(2,2) = a(1,1) * a(3,3) - a(3,1) * a(1,3)
      b(2,3) = a(2,1) * a(1,3) - a(1,1) * a(2,3)
      b(3,1) = a(2,1) * a(3,2) - a(3,1) * a(2,2)
      b(3,2) = a(3,1) * a(1,2) - a(1,1) * a(3,2)
      b(3,3) = a(1,1) * a(2,2) - a(2,1) * a(1,2)

      det = a(1,1) * b(1,1) + a(1,2) * b(2,1) + a(1,3) * b(3,1)

      do i = 1,3
         do j = 1,3
            b(i,j) = b(i,j) / det
         end do
      end do

      return
      end


c====================================================================
c====================================================================
c  SUBROUTINE: 
c  Solve simultaneous equations using LU decomposition (Crout's method)
c  Result is stored in b(i)
c
c  Reference: "Numerical Recipes" Section 2.3  p. 31
c
c--------------------------------------------------------------------

      subroutine simeq(n,a,b)

      implicit double precision (a-h,o-z)

      real*8 a(n,n), b(n), index(n)

      call LU_Decomp(n,a,index)
      call LU_BackSub(n,a,index,b)

      return
      end

c====================================================================
c====================================================================
c  SUBROUTINE: 
c  Calculate the inverse of a matrix using
c  LU decomposition (Crout's method)
c
c  Reference: "Numerical Recipes" Section 2.3  p. 31
c
c--------------------------------------------------------------------

      subroutine inverse(n,a,b)

      implicit double precision (a-h,o-z)

      real*8 a(n,n), b(n,n), c(n,n), index(n)

      do i = 1,n
         do j = 1,n
            c(i,j) = a(i,j)
         end do
      end do

      do i = 1,n
         do j = 1,n
            b(i,j) = 0.0
         end do
         b(i,i) = 1.0
      end do

      call LU_Decomp(n,c,index)
      do j = 1,n
         call LU_BackSub(n,c,index,b(1,j))
      end do

      return
      end

c====================================================================
c====================================================================
c  SUBROUTINE: 
c  This sub performs an LU Decomposition (Crout's method) on the
c  matrix "a". It uses partial pivoting for stability. The index()
c  vector is used for the partial pivoting.  The v() vector is
c  a dummy work area.
c
c  Reference: "Numerical Recipes" Section 2.3  p. 31
c
c--------------------------------------------------------------------

      subroutine LU_Decomp(n,a,index)

      implicit double precision (a-h,o-z)

      real*8 a(n,n), index(n), v(n)

      tiny = 1.0e-20

c--------------------------------------------------------------------
c  Loop over the rows to get the implicit scaling info.
c--------------------------------------------------------------------

      do i = 1,n
         a_max = 0.0
         do j = 1,n
            a_max = max(a_max,abs(a(i,j)))
         end do !j
         v(i) = 1.0 / a_max
      end do !i

c--------------------------------------------------------------------
c  Begin big loop over all the columns.
c--------------------------------------------------------------------

      do j = 1,n

         do i = 1,j-1
            sum = a(i,j)
            do k = 1,i-1
               sum = sum - a(i,k) * a(k,j)
            end do
            a(i,j) = sum
         end do

         a_max = 0.0
         do i = j,n
            sum = a(i,j)
            do k = 1,j-1
               sum = sum - a(i,k) * a(k,j)
            end do
            a(i,j) = sum
            dummy = v(i) * abs(sum)
            if ( dummy .gt. a_max ) then
               imax = i
               a_max = dummy
            end if
         end do

c--------------------------------------------------------------------
c  Pivot rows if necessary.
c--------------------------------------------------------------------

         if ( j .ne. imax ) then
            do k = 1,n
               dummy = a(imax,k)
               a(imax,k) = a(j,k)
               a(j,k) = dummy
            end do
            v(imax) = v(j)
         end if
         index(j) = imax

c--------------------------------------------------------------------
c  Divide by the pivot element.
c--------------------------------------------------------------------

         if ( a(j,j) .eq. 0.0 ) a(j,j) = tiny
         if ( j .ne. n ) then
            dummy = 1.0 / a(j,j)
            do i = j+1,n
               a(i,j) = a(i,j) * dummy
            end do
         end if

      end do !j

      return
      end

c====================================================================
c====================================================================
c  SUBROUTINE: 
c  Solves a set of simultaneous equations by doing back substitution.
c  The answer in returned in the b() vector.  The a(,) matrix
c  must have already been "LU Decomposed" by the above subroutine.
c
c  Reference: "Numerical Recipes" Section 2.3  p. 31
c
c--------------------------------------------------------------------

      subroutine LU_BackSub(n,a,index,b)

      implicit double precision (a-h,o-z)

      real*8 a(n,n), index(n), b(n)

      ii = 0

c--------------------------------------------------------------------
c  Do the forward substitution.
c--------------------------------------------------------------------

      do i = 1,n
         m = index(i)
         sum = b(m)
         b(m) = b(i)
         if ( ii .ne. 0 ) then
            do j = ii,i-1
               sum = sum - a(i,j) * b(j)
            end do
         else if ( sum .ne. 0.0 ) then
            ii = i
         end if
         b(i) = sum
      end do

c--------------------------------------------------------------------
c  Do the back substitution.
c--------------------------------------------------------------------

      do i = n,1,-1
         sum = b(i)
         if ( i .lt. n ) then
            do j = i+1,n
               sum = sum - a(i,j) * b(j)
            end do
         end if
         b(i) = sum / a(i,i)
      end do

      return
      end

c====================================================================
c====================================================================
c  SUBROUTINE: 
c  Restore a symmetric 4th rank tensor stored in Voigt notation
c  back to its 4th rank form.
c
c--------------------------------------------------------------------

      subroutine Voigt_to_forth(b,a)

      implicit double precision (a-h,o-z)

      real*8 a(3,3,3,3), b(6,6)

      do i = 1,3
       do j = 1,3
        ia = i
        if (i.ne.j) ia=9-i-j
        do k = 1,3
         do l = 1,3
          ib = k
          if (k.ne.l) ib=9-k-l
          a(i,j,k,l) = b(ia,ib)
          if (ia.gt.3) a(i,j,k,l) = a(i,j,k,l) / 2
          if (ib.gt.3) a(i,j,k,l) = a(i,j,k,l) / 2
         end do
        end do
       end do
      end do

      return
      end


c====================================================================
c====================================================================
c  SUBROUTINE: 
c  Store a SYMMETRIC 4th rank tensor in Voigt notation.
c
c--------------------------------------------------------------------

      subroutine forth_to_Voigt(a,b)

      implicit double precision (a-h,o-z)

      real*8 a(3,3,3,3), b(6,6)

      do i = 1,3
       do j = i,3   ! not 1 !!!
        ia = i
        if (i.ne.j) ia=9-i-j
        do k = 1,3
         do l = k,3 ! not 1 !!!
          ib = k
          if (k.ne.l) ib=9-k-l
          b(ia,ib) = a(i,j,k,l)
         end do
        end do
       end do
      end do

      return
      end



c====================================================================
c====================================================================
c  FUNCTION: 
c  Perform x**y but while retaining the sign of x.
c
c--------------------------------------------------------------------

      function power(x,y)

      implicit double precision (a-h,o-z)

      if (x.eq.0.0) then
        if (y.gt.0.0) then
          power = 0.0
        else if (y .lt. 0.0) then
          power = 1.0d+300
        else
          power = 1.0
        end if
      else
         power = y * log10(abs(x))
         if (power .gt. 300.) then
           power = 1.d+300
         else
           power = 10.d0 ** power
         end if
         if (x .lt. 0.0) power = -power
      end if

      return
      end

c===================================================================
c===================================================================
c  SUBROUTINE:
c  Print out Euler angles in Kocks notation.
c
c-------------------------------------------------------------------

      subroutine kocks_angles(angle ,array1)

      implicit double precision (a-h,o-z)

      real*8 angle(3), array1(3,3)

      pi = 4 * DATAN(1.D0)

      if (DABS(array1(3,3)) .gt. 0.99999) then
        psi   = DATAN2(array1(2,1),array1(1,1))
        theta = 0.0
        phi   = 0.0
      else
        psi   = DATAN2(array1(2,3),array1(1,3))
        theta = DACOS(array1(3,3))
        phi   = DATAN2(array1(3,2),-array1(3,1))
      end if
      angle(1) = psi
      angle(2) = theta
      angle(3) = phi
      return
      end

c====================================================================
c====================================================================
c  SUBROUTINE:
c  Calculate the determinant of a 3 x 3 matrix.
c
c--------------------------------------------------------------------

       subroutine calc_Determinant(a,det)
       
       implicit double precision (a-h,o-z)
       
       real*8 a(3,3)
       
       b1 = a(2,2) * a(3,3) - a(3,2) * a(2,3)    ! Perform calculation
       b2 = a(3,1) * a(2,3) - a(2,1) * a(3,3)
       b3 = a(2,1) * a(3,2) - a(3,1) * a(2,2)

       det = a(1,1) * b1 + a(1,2) * b2 + a(1,3) * b3

       return
       end

