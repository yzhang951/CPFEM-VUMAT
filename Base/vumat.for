c===================================================================
c                                                                  c
c VUMAT: Rate-dependent crystal plasticity with back stresses      c
c Basic crystal plasticity model for beginners                     c
c Developed by: Yin Zhang, Peking University                       c
c Contact: yinzhang@pku.edu.cn                                     c 
c Homepage: https://yzhang951.github.io/                           c
c Runs on Abaqus/Explicit (working on 6.13 and 2022)               c
c Latest version: 09/24/2024                                       c
c                                                                  c
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
       character*80 cmname, FILE1
c      
c     Local arrays from vgetvrm are dimensioned to
c     maximum block size (maxblk)
       parameter(nrData=6, num_el=125)
       character*3 cData(maxblk*nrData)
       dimension rData(maxblk*nrData), jData(maxblk*nrData)
       real*8    psi(num_el,3), F_p(3,3)
       PI = 4.D0*ATAN(1.D0)

       if(totalTime.le.dt.or.dt.eq.1.D0) then
c     Read in Euler angles from external file and store in STATENEW           
         write(*,*) "Read in Euler Angles"
         FILE1 = '/data/yinzhang/GitHub/CPFEM-VUMAT/Base/aeuler'
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

       parameter(num_slip_sys = 12, ! Total number of slip systems
     &           num_el = 125)      ! Total number of elements, need to modify for different simulations!

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
     & dtime, PI, one, plastic_work_inc, dPEEQ

       PI = 4.D0*DATAN(1.D0)
       one = 1.D0
       dtime = dt                       ! Timestep
c====================================================================
c      TOTAL ISVs:
c====================================================================
c      3  - Euler angles
c      9  - relative Spin Matrix
c      9  - F_p(i,j)
c      12 - dummy variables for 12 slip systems
c      12 - back_stress(12) for 12 slip systems
c      12 - threshold(12) for 12 slip systems
c      12 - gamma_dot(12) for 12 slip systems
c====================================================================

C       write(*,*) "Material Properties"
    
c====================================================================
c       Assign props() array to logical variable names
c====================================================================
c-------------------- Macroscopic Properties  -----------------------
       C11          = props(1)  ! Elastic constants C11
       C12          = props(2)  ! Elastic constants C12
       C44          = props(3)  ! Elastic constants C44, shear modulus
       psi_ang      = props(4)  ! Euler angles
       theta_ang    = props(5)
       phi_ang      = props(6)
c-------------------- Rate-dependent flow rule ----------------------       
       g0           = props(7)  ! Reference shear rate
       am           = props(8)  ! Strain rate sensitivity
       s0           = props(9) ! Initial slip resistance
c-------------------- Hardening -------------------------------------
       h0           = props(10) ! Initial hardening rate
       s_sat        = props(11) ! Saturated slip resistance
       ah           = props(12) ! Hardening exponent 
c-------------------- Threshold Stress Evolution --------------------       
       A00          = props(13) ! Dislocation self interaction coeff
       A01          = props(14) ! Dislocation interaction coeff
c-------------------- Backstress Evolution --------------------------
       C_back_one   = props(15) ! Backstress 1 generation parameter
       C_dyn_one    = props(16) ! Backstress 1 dynamic recovery para
       C_back_two   = props(17) ! Backstress 2 generation parameter
       C_dyn_two    = props(18) ! Backstress 2 dynamic recovery para


c------------ Euler angles for single xtal simulation ---------------
       psi_t(1) = psi_ang*PI/180.
       psi_t(2) = theta_ang*PI/180.
       psi_t(3) = phi_ang*PI/180.

c--------------------- Disl interaction coeff matrix ----------------
       do i = 1,num_slip_sys
         do j = 1,num_slip_sys
           Aii(i,j) = A01
         end do
         Aii(i,i) = A00
       end do

C====================================================================
C====================================================================
C
C
C     Abaqus Explicit sends in data in blocks of NBLOCK=128
C     elements at a time. 
C
C====================================================================       
C==================================================================== 
C====================================================================       
C      Start loop over Nblock       
C====================================================================
       do km = 1,nblock
c====================================================================       
c      Initialize internal state variables for the first time step
c====================================================================

         if (totalTime.le.dt.or.dt.eq.one) then
c           write(*,*) "Initializing ISVs, first time step"
c--------------- Read in Euler angles from external file ------------
c--- Important!: comment this section for poly xtal simulations -----
c--- the Euler angles will be read from file 'aeuler' through VUSDFLD
           psi_t(1) = psi_ang*PI/180.0
           psi_t(2) = theta_ang*PI/180.0
           psi_t(3) = phi_ang*PI/180.0
c--------------------- Initialize F_p and F_p_inv. ------------------
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
c------------ Initialize shear rate and back stress -----------------
           do i = 1,num_slip_sys
             gamma_dot_0(i) = 0.D0
             back_stress_one_0(i) = 0.D0
             back_stress_two_0(i) = 0.D0
           end do
c-------------------- Initialize threshold stress --------------------
           do i = 1,num_slip_sys
             threshold_0(i) = s0
           end do
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
c             rho_0(i) = stateOld(km,n)
           end do
c-------------------- Read in threshold stress 46-57 ----------------
           do i = 1, num_slip_sys
             n = n + 1
             threshold_0(i) = stateOld(km,n)
           end do
c-------------------- Read in back stress one 58-69 -----------------
           do i = 1, num_slip_sys
             n = n + 1
             back_stress_one_0(i) = stateOld(km,n)
           end do
c-------------------- Read in back stress two 70-81 -----------------
           do i = 1, num_slip_sys
             n = n + 1
             back_stress_two_0(i) = stateOld(km,n)
           end do
c-------------------- Read in PEEQ 82 -------------------------------
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
         if (nshr.eq.3) then
C        3D problem
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
           F_0(2,3) = 0.D0
           F_0(3,1) = 0.D0
           F_0(2,1) = defgradOld(km,5)
           F_0(3,2) = 0.D0
           F_0(1,3) = 0.D0


           F_t(1,1) = defgradNew(km,1)
           F_t(2,2) = defgradNew(km,2)
           F_t(3,3) = defgradNew(km,3)
           F_t(1,2) = defgradNew(km,4)
           F_t(2,3) = 0.D0
           F_t(3,1) = 0.D0
           F_t(2,1) = defgradNew(km,5)
           F_t(3,2) = 0.D0
           F_t(1,3) = 0.D0
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
         call calc_L_p(num_slip_sys, gamma_dot_0, xs0, xm0, xL_p)
         call calc_F_p_subend(xL_p, dtime, F_p_0, F_p_t)
         call calc_F_el(F_p_t, F_t, F_el_inv_t, F_el_t)
         call calc_E_el(F_el_t, E_el)
         call calc_Spk2(C, E_el, spk2)
         call calc_Sig(F_el_t, spk2, sig)
         call calc_Tau(num_slip_sys, tau, spk2, xs0, xm0)

c==================================================================== 
c        Update the rotation matrix
c====================================================================
c         write(*,*) "Updating the rotation matrix", F_p_0, F_p_t
         if (nshr.eq.3) then
c        3D problem
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
           U_stretch(2,3) = 0.D0
           U_stretch(3,1) = 0.D0
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
         call calc_back(num_slip_sys, C_back_one, C_dyn_one, 
     &   gamma_dot_0, back_stress_one_0, dtime, back_stress_one_t)
         call calc_back(num_slip_sys, C_back_two, C_dyn_two, 
     &   gamma_dot_0, back_stress_two_0, dtime, back_stress_two_t)
C         write(*,*) "Evolution of back stress: ", back_stress_t
         call calc_threshold(num_slip_sys, Aii, h0, threshold_0, 
     &   s_sat, ah, gamma_dot_0, dtime, threshold_t)
c         write(*,*) "Evolution of threshold stress: ", threshold_t
         call f(num_slip_sys, g0, tau, back_stress_one_t, 
     &   back_stress_two_t, threshold_t, am, gamma_dot_t)
C         write(*,*) "Evolution of gamma_dot: ", gamma_dot_t

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
     &       DABS(tau(i)*gamma_dot_t(i)*dtime)
         end do
         enerInelasNew(km) = enerInelasOld(km) + 
     &                       plastic_work_inc/density(km)

         call calc_dPEEQ(sig, plastic_work_inc, dPEEQ)


c=====================================================================    
c        Store in internal state variables for nonzero time step
c=====================================================================
!         write(*,*) "Store in state variables "

         n = 0
c-------------------- Store in Euler Angles 1-3 ----------------------
         do i = 1,3
           n = n + 1
           stateNew(km,n) = psi_t(i)
         end do 
c-------------------- Store in Relative Spin Matrix 4-12---------------
         do i = 1,3
           do j = 1,3
             n = n + 1
             stateNew(km,n) = rel_spin_t(i,j)
           end do
         end do 
c-------------------- Store in F_p 13-21 -----------------------------
         do i = 1,3
           do j = 1,3
             n = n + 1
             stateNew(km,n) = F_p_t(i,j) 
           end do
         end do
c-------------------- Store in dislocation shear rate 22-33 ----------
         do i = 1, num_slip_sys
           n = n + 1
           stateNew(km,n) = gamma_dot_t(i)
         end  do
c-------------------- Store in accumlated shear34-45 -----------------
         do i = 1, num_slip_sys
           n = n + 1
           stateNew(km,n) = stateOld(km,n)+DABS(gamma_dot_t(i))*dtime
         end do
c-------------------- Store in threshold stress 46-57 ----------------
         do i = 1, num_slip_sys
           n = n + 1
           stateNew(km,n) = threshold_t(i)
         end do
c-------------------- Store in back stress one 58-69 -----------------
         do i = 1, num_slip_sys
           n = n + 1
           stateNew(km,n) = back_stress_one_t(i)
         end do
c-------------------- Store in back stress two 70-81 -----------------
         do i = 1, num_slip_sys
           n = n + 1
           stateNew(km,n) = back_stress_two_t(i)
         end do
c-------------------- Store in PEEQ 82 -------------------------------
         n = n + 1
         stateNew(km,n) = stateOld(km,n) + dPEEQ
c-------------------- Store in plastic disspation 83------------------
         n = n + 1
         stateNew(km,n) = enerInelasNew(km)
c         write(*,*) "Next Element!"
       
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

