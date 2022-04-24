C***********************************************************************
C
C                          ABAQUS  VUMAT  FOR
C
C                                   
C                  RATE-DEPENDENT POWER-LAW CRYSTAL PLASTICITY.
C      
C
C **
C
C  TO BE USED WITH ABAQUS VERSION 6.13
C
C  LAST MODIFIED June 7 2020
C
C   This VUMAT uses an explicit Euler-forward method for integrating the
C   flow rule and the evolution equations for the slip resistances
C
C                
C
C
C****************************************************************************
C       The model contains 12 material constants, and calculation parameters
C       which are:
C
C      C11       = props(1) -- Elastic modulus, Pa
C      C12       = props(2) -- Elastic modulus, Pa
C      C44       = props(3) -- Elastic modulus, Pa
C      GDOT0     = props(4) -- Reference shearing rate sec^{-1}
C      AM        = props(5) -- Power law
C      S0        = props(6) -- Initial val. of the  def. resis., Pa
C      AH0       = props(7) -- Value of the initial hardening rate, Pa
C      SSAT      = props(8) -- Saturation value for the def. resis., Pa
C      RHARD     = props(9) -- Exponent in hardening equation
C      QL        = props(10) -- Latent hardening parameter
C      TEXCALTIME= props(21) -- Set equal to the time at which you wish to
C                               first calculate the texture. 
C      TEXFREQ   = props(22) -- Set equal to the time increment 
C                               to reset TEXCALTIME = TEXCALTIME + TEXFREQ
C
C                                 If texture calculation is to be performed
C                                 only at the end of the step, then set
C                                 TEXCALTIME very close to the final time
C                                 and set TEXFREQ to a large number.
C
C            These properties are specified by
C
C      *USER MATERIAL,CONSTANTS=12
C      170.E9,124.E9,75.E9,0.001,0.012,16.E6,250.E6,190.E6,
C      2.25,1.4,0.499, 1.
C
C      in the main input deck.
C
C
C      The model also use 3 solution dependent variables 
C          SDV1  = EBARP_T    ---  Equiv. tensile plastic strain
C          SDV2  = EBARDOTP_T ---  Equiv. tensile plastic strain rate 
C          SDV3  = EQSTRESS_T      ---  Equiv. tensile stress 
C       SDV15-20  = Backstress  --- 6 Backstress component 
C      which are specified by
C       
C       *DEPVAR
C      18 
C      in the main input deck.
C 
C      All other state variables are stored in large matrices.
C
C**************************************************************
C THE MODEL NEEDS THE FOLLOWING EXTERNAL FILES
C 1.      vumatcomm:
C
C         This file contains the required information for the size of 
C         of the finite element problem: 
C
C         noel = # of finite elements.
C         ncrys= # of different crystal orientations at each integ. point.
C         nslip= # of slip systems in each crystal
C
C          PARAMETER(NOEL=480,NCRYS=1,NSLIP=12)
C
C 2.      aeuler:
C
C           This file contains information concerning the euler angles
C           for the orientations of crystals making up a materail point.
c            
C            For a full finite element calculation of a single crystal, set
C
C            icrys_flag=1, 
C
C            and give the euler angles 
C
C            theta, phi, omega, b11, b22, b33, phase
C
C            for the orientation of the crystal.
C
C            For a full finite element calculation in which each finite element potentially
C            has a different orientation, set 
C
C            icrys_flag=2, 
C
C            and give the euler angles
C
C            theta, phi, omega, b11, b22, b33, phase
C
C            for each finite element. The number of
C            entries for theta,phi,omega, need to be as many as total
C            the number of elements, NOEL
C
C            For a taylor model calculation in which each integration point has the same initial
C            set of crystal orientations, set 
C
C            icrys_flag=1,
C            theta,phi,omega 
C            
C            the number of entries for theta,phi,omega need  to be as many as ncrys.
C
C
C            Example: Entries for a finite element calculation of a 110 oriented sxal
C
C            ICRYS_FLAG : 1 --- TAYLOR MODEL, SAME SET OF CRYSTALS AT  INTEGRATION POINTS
C                         2 --- FEM MODEL,    DIFFERENT   CRYSTALS AT INTEGRATION POINTS
C            1
C            EULER ANGLES FOR EACH CRYSTAL AT AN INTEGRATION POINT
C            THETA PHI  OMEGA
C            -45. 0. .0
C
C
C  3. slipsys:
C
C      This file contains the maximum number of slipsystems and the
C      Miller indices {hkl}<uvw> for the slip systems
C
C      EXAMPLE: for fcc materials:
C
C      MAXIMUM NUMBER OF SLIP SYSTSEMS
C      12
C      SLIP SYSTEMS
C      1  1  1  1 -1  0
C      1  1  1  1  0 -1
C      1  1  1  0  1 -1
C      -1  1  1  1  0  1
C      -1  1  1  1  1  0
C      -1  1  1  0  1 -1
C       1 -1  1  1  0 -1
C       1 -1  1  0  1  1
C       1 -1  1  1  1  0
C       1  1 -1  1 -1  0
C       1  1 -1  1  0  1
C       1  1 -1  0  1  1
C
C*******************************************************************************
      SUBROUTINE VUMAT (
C      Read only (unmodifiable) variables :-
     +                    NBLOCK, NDIR, NSHR, NSTATEV, NFIELDV,
     +                    NPROPS, LANNEAL, STEP_TIME, TOTAL_TIME,
     +                    DT, CMNAME, COORD_MP, CHAR_LENGTH, PROPS,
     +                    DENSITY, STRAIN_INC, REL_SPIN_INC,
     +                    TEMP_OLD, STRETCH_OLD, DEFGRAD_OLD,
     +                    FIELD_OLD, STRESS_OLD, STATE_OLD, 
     +                    ENER_INTERN_OLD, ENER_INELAS_OLD, TEMP_NEW,
     +                    STRETCH_NEW, DEFGRAD_NEW, FIELD_NEW, 
C      Read and write (modifiable) variables :-
     +                    STRESS_NEW, STATE_NEW, ENER_INTERN_NEW,
     +                    ENER_INELAS_NEW)
     
      INCLUDE 'VABA_PARAM.INC'
      
      DIMENSION COORD_MP(NBLOCK,*), CHAR_LENGTH(NBLOCK), 
     +            PROPS(NPROPS), DENSITY(NBLOCK), 
     +            STRAIN_INC(NBLOCK,NDIR+NSHR), 
     +              REL_SPIN_INC(NBLOCK,NSHR), TEMP_OLD(NBLOCK), 
     +              STRETCH_OLD(NBLOCK,NDIR+NSHR),
     +              DEFGRAD_OLD(NBLOCK,NDIR+NSHR+NSHR), 
     +              FIELD_OLD(NBLOCK,NFIELDV), 
     +            STRESS_OLD(NBLOCK,NDIR+NSHR),
     +              STATE_OLD(NBLOCK,NSTATEV), 
     +            ENER_INTERN_OLD(NBLOCK),
     +              ENER_INELAS_OLD(NBLOCK), TEMP_NEW(NBLOCK), 
     +              STRETCH_NEW(NBLOCK,NDIR+NSHR), 
     +              DEFGRAD_NEW(NBLOCK,NDIR+NSHR+NSHR), 
     +              FIELD_NEW(NBLOCK,NFIELDV), 
     +            STRESS_NEW(NBLOCK,NDIR+NSHR),
     +              STATE_NEW(NBLOCK,NSTATEV), 
     +            ENER_INTERN_NEW(NBLOCK),
     +              ENER_INELAS_NEW(NBLOCK)

      CHARACTER*8 CMNAME
            
      CHARACTER*80 WRKDIR, FILE1, FILE2, FILE3, FILE4, FILE5, FILE6
      CHARACTER*3 MESSAGES
C****************************************************************
C         Specify the number of elements, the number of crystals
C         at each integration point, and the number of slip systems
C         in each crystal.
C****************************************************************
C     BH : backstress hardening
      PARAMETER(NOEL=990,NCRYS=1,NSLIP=12,BH=0.0)
      
      PARAMETER(ZERO=0.D0, ONE=1.D0, ONE_HALF=0.5D0, TWO=2.D0,
     +              ONE_THIRD=1.D0/3.D0, TWO_THIRD=2.D0/3.D0,
     +              THREE_HALF=1.5D0, THREE=3.D0)

      INTEGER I, J, K, L, II, JJ, KK, LL, KM, 
     +            INOEL,ICRYS,ISLIP,JSLIP, NOEL_PT,
     +        ICRYS_FLAG,IEULERERR,IPHASE
      INTEGER PHASE(NOEL)

      REAL*8 C11,C12,C44,GDOT0,AM,S0,AH0,SSAT,RHARD
      REAL*8 QL,TEXCALTIME,TEXFREQ
      
      REAL*8 DTIME,EBARP_T,EBARDOTP_T,EQSTRESS_T,
     +       PI,TH, PH, OM, THK, PHIK, PSIK, VOLUME,
     +            AMV_C_MAG, ANV_C_MAG, EIG11,EIG22,EIG33,
     +       DET_FP_TAU_INV,DET_FSTAR_TAU, DET_FSTAR_TAU_INV,
     +       DET_U_TAU_INV,VAL,
     +       PBAR, EQSTRESS_TAU, PLASTIC_WORK_INC,DEQPS,
     +       TAU_STAR,RATIO,FRAC,DELTAG,
     +       STRESS_POWER,SFRAC,TAU_ALPHA,S_ALPHA,Interval  
      REAL*8 C11_FCC,C12_FCC,C44_FCC,
     +       GDOT0_FCC,AM_FCC,S0_FCC,AH0_FCC,SSAT_FCC,
     +       RHARD_FCC,QL_FCC,
     +       C11_BCC,C12_BCC,C44_BCC,
     +       GDOT0_BCC,AM_BCC,S0_BCC,AH0_BCC,SSAT_BCC,
     +       RHARD_BCC,QL_BCC

      REAL*8 AIDENT(3,3),VOL(NCRYS), 
     +       Q(3,3),Q_CRYS(NCRYS,3,3), Q_ALL(NOEL,NCRYS,3,3),
     +       VOL_FRAC_ALL(NOEL,NCRYS), EIGEN(NOEL,3),
     +       AMV_C(NSLIP,3), ANV_C(NSLIP,3),
     +       SMAT_C(NSLIP,3,3),SMAT_FCC(NSLIP,3,3),SMAT_BCC(NSLIP,3,3),
     +       SMAT_0_ALL(NOEL,NCRYS,NSLIP,3,3),     
     +       FP_T_INV_ALL(NOEL,NCRYS,3,3),
     +       FP_T_INV(3,3), FP_TAU_INV(3,3),
     +       S_T_ALL(NOEL,NCRYS,NSLIP), 
     +       DGAMMA_T_ALL(NOEL,NCRYS,NSLIP),  
     +       DGAMMA_T(NSLIP), DGAMMA_TAU(NSLIP),
     +       ACCGAM_T_ALL(NOEL,NCRYS,NSLIP),
     +       ACCGAM_T(NSLIP), ACCGAM_TAU(NSLIP),
     +       QLAT(NSLIP,NSLIP),
     +       F_T(3,3), F_TAU(3,3), 
     +       R_TAU(3,3), U_TAU(3,3),U_TAU_INV(3,3),
     +       ELMAT_C(3,3,3,3), ELMAT_G_ALL(NOEL,NCRYS,3,3,3,3),
     +       FSTAR_TAU(3,3), RSTAR_TAU(3,3), USTAR_TAU(3,3),
     +       ESTAR_TAU(3,3), TSTAR_TAU(3,3), Bij(3,3), 
     +       T_TAU(3,3), TBAR_TAU(3,3),TBARDEV_TAU(3,3),STRESS(3,3), 
     +       S_T(NSLIP), AH_T(NSLIP), AHMAT_T(NSLIP,NSLIP),      
     +       S_TAU(NSLIP),RSS_TAU(NSLIP),RSS_TAU_EFF(NSLIP),
     +       QTEX(3,3), AUX1(3,3), AUX2(3,3) 
C**********************************************************************
C                Set up common blocks to store the relevant
C                internal variables of the model
C**********************************************************************
        COMMON/VUMATPRP/C11_FCC,C12_FCC,C44_FCC,
     +                GDOT0_FCC,AM_FCC,S0_FCC,AH0_FCC,SSAT_FCC,
     +                RHARD_FCC,QL_FCC,
     +                C11_BCC,C12_BCC,C44_BCC,
     +                GDOT0_BCC,AM_BCC,S0_BCC,AH0_BCC,SSAT_BCC,
     +                RHARD_BCC,QL_BCC,
     +                TEXCALTIME,TEXFREQ
        COMMON/VUMATCST1/Q_ALL, SMAT_0_ALL,VOL_FRAC_ALL
        COMMON/VUMATCST2/AIDENT, PI, QLAT, VOLUME
        COMMON/VUMATMSGS/MESSAGES     
        COMMON/VUMATVAR/ELMAT_G_ALL,FP_T_INV_ALL,S_T_ALL,DGAMMA_T_ALL,
     +                  ACCGAM_T_ALL, EIGEN
        COMMON/VUMATFLG/NOEL_PT, PHASE
C**********************************************************************
C             BEGIN COMPUTATION
C**********************************************************************
        DTIME = DT
C**********************************************************************
C         Abaqus Explicit sends in data in blocks of NBLOCK=128 elements
C         (or since these are reduced integration elements the same number 
C          of integration or material points) at a time. Set up an element
C         pointer NOEL_PT which is intialized in a common block. This pointer 
C         is incremented by NBLOCK after every loop over NBLOCK elements.
C**********************************************************************
      IF (NOEL_PT .EQ. NOEL) NOEL_PT = 0
        
C**********************************************************************
C       write(*,*) 'NBLOCK',NBLOCK
C       write(*,*) 'density',density(NBlock)
       DO 1000 KM = 1,NBLOCK      
C
C     During the data-check phase ABAQUS
C     calls the VUMAT with a fictitious set of deformation gradients
C     and totalTime, stepTime, and DT all equal to 1.0.
C     During this intial dummy step  initialize some variables and 
C     open some files for writing messages 
C
C           if (NOEL_PT+KM.eq.1) then
C           write(*,*) 'read', DGAMMA_T_ALL(1,1,4),DGAMMA_T_ALL(1,1,5)
C           endif       
      
      IF (TOTAL_TIME .EQ. 0.0D0 .OR. (TOTAL_TIME-DT) .EQ. 0.0D0) THEN

      FILE1 = '~/AlHEA/response/aeuler_BCC'
      FILE2 = '~/AlHEA/3/slipsys'
      FILE3 = '~/AlHEA/response/VU_MESSAGES'
      FILE4 = '~/AlHEA/response/VU_TEXTURES'
      FILE5 = '~/AlHEA/response/VU_SLIPS'
      
      
      DTIME = DT
C       OPEN(UNIT=61,FILE='temp',STATUS='UNKNOWN')
C       write(61,*) NCRYS
C       close(61)
C       stop          
       PI = 4.D0*DATAN(1.D0)
          
         AIDENT(1,1) = ONE
         AIDENT(1,2) = ZERO
         AIDENT(1,3) = ZERO
         AIDENT(2,1) = ZERO
         AIDENT(2,2) = ONE
         AIDENT(2,3) = ZERO
         AIDENT(3,1) = ZERO
         AIDENT(3,2) = ZERO
         AIDENT(3,3) = ONE
C
C              Store the material parameters from the props vector
C
       C11_FCC   = props(1)   
       C12_FCC   = props(2)   
       C44_FCC   = props(3)   
       GDOT0_FCC = props(4)   
       AM_FCC    = props(5)     
       S0_FCC    = props(6)  
       AH0_FCC   = props(7)  
       SSAT_FCC  = props(8)  
       RHARD_FCC = props(9)  
       QL_FCC    = props(10)  

       C11_BCC   = props(11)   
       C12_BCC   = props(12)   
       C44_BCC   = props(13)   
       GDOT0_BCC = props(14)   
       AM_BCC    = props(15)     
       S0_BCC    = props(16)  
       AH0_BCC   = props(17)  
       SSAT_BCC  = props(18)  
       RHARD_BCC = props(19)  
       QL_BCC    = props(20)  


       TEXCALTIME= props(21)
       TEXFREQ   = props(22)
C       write(*,*) props(6),props(7),props(8),props(16),props(17),props(18)
         EBARP_T    = ZERO
         EBARDOTP_T = ZERO
         EQSTRESS_T = ZERO
          
         STATE_OLD(KM,1)  = EBARP_T    
         STATE_OLD(KM,2)  = EBARDOTP_T 
         STATE_OLD(KM,3)  = EQSTRESS_T
C      Back Stress      
       STATE_OLD(KM,10) = ZERO
       STATE_OLD(KM,11) = ZERO
       STATE_OLD(KM,12) = ZERO
       STATE_OLD(KM,13) = ZERO
       STATE_OLD(KM,14) = ZERO
       STATE_OLD(KM,15) = ZERO
        
C
C              Read in the Euler angles for all the crystals
C
C         FILE1 = '/usr/home/cricket/tzhu/SXAL_VUMAT/AEULER'

C           ICRYS_FLAG = 1
           OPEN(UNIT=65,FILE=FILE1,STATUS='UNKNOWN')
           READ(65,*)
           READ(65,*)
           READ(65,*) ICRYS_FLAG
           READ(65,*) 
           READ(65,*) 
           READ(65,*) 
           READ(65,*) 
      
          IF (ICRYS_FLAG .EQ. 1) THEN
C
C          Then, this is a Taylor model calculation,
C          and the same set of crystals are used at all the
C          integration points in all the elements. 
C
c          If ICRYS=2,
C          then it is a full finite-element calculation
C          of a polycrystal; that is, each finite element
C          is a separate crystal.
C
C           The euler angles are in two possible formats:
C           The first format is the Kalidindi format in
C           which all orientations are weighted equally, with unit weight.
C           The second is the Kocks format, where the different
C           orientations are given different weights. The latter
C           method can sometimes save computation time in complex
C           calculations. The Kalidindi orientations are listed as
C           (TH,PH OM), whereas the Kocks orientations are denoted as
C           (THK,PHIK,PSIK). The relationship between the two is
C
C            TH = THK
C              PH = 90. + PHIK
C               OM = 90. - PSIK
C
C            ALL CALCULATIONS ARE DONE IN THE (TH, PHI, OM) SYSTEM.
C
C            The quantity VOLUME used below is the ``volume fraction''
C            of each orientation in  Kocks' representation. The same 
C            quantity in the Kalidindi format is the simply a counter 
C            for the number of grains used to represent a material point 
C            in the Taylor model.
C
C            The quantity NCRYS is the total number of crystal orientations
C            used to represent the texture.
C
C        The e11 e22 e33 represent the eigenstrain for each element.
C        Modified by Yin Zhang, 27th Nov 2017
C
         VOLUME = ZERO
         DO 10 ICRYS = 1,NCRYS
C             pause
              READ(65,*) PH, TH, OM, EIG11, EIG22, EIG33, IPHASE

                VOL(ICRYS) = ONE
C            
C               CONVERT THE ANGLES INTO RADIANS
C
               TH = TH*PI/180.
              PH = PH*PI/180.
              OM = OM*PI/180.
            EIGEN(ICRYS,1) = 0
            EIGEN(ICRYS,2) = 0
            EIGEN(ICRYS,3) = 0
            
                VOLUME = VOLUME + VOL(ICRYS)
            
C
C                  Calculate the rotation matrix [Q] for each crystal
C
              CALL ROTMAT(TH,PH,OM,Q)
C
C                  Store the rotation matrix for all the crystals
C                  at an integration point in the matrix Q_CRYS
C                  which is indexed by the counter ICRYS = 1 TO NCRYS
C              
              Q_CRYS(ICRYS,1,1) = Q(1,1)
              Q_CRYS(ICRYS,1,2) = Q(1,2)
              Q_CRYS(ICRYS,1,3) = Q(1,3)
              Q_CRYS(ICRYS,2,1) = Q(2,1)
              Q_CRYS(ICRYS,2,2) = Q(2,2)
              Q_CRYS(ICRYS,2,3) = Q(2,3)
              Q_CRYS(ICRYS,3,1) = Q(3,1)
              Q_CRYS(ICRYS,3,2) = Q(3,2)
              Q_CRYS(ICRYS,3,3) = Q(3,3)
10            CONTINUE              

C
C                  Store the rotation matrix for all the elements
C                  in the matrix Q_ALL which is indexed by the counters
C                  INOEL = 1 TO NOEL and ICRYS=1 TO NCRYS
C
C                  CALCULATE THE VOLUME FRACTION OF EACH CRYSTAL ORIENTATION
C                  IN EACH ELEMENT AND STORE IT IN THE MATRIX
C                  VOL_FRAC_ALL
C
            DO 11 INOEL = 1,NOEL
            DO 11 ICRYS = 1,NCRYS
              Q_ALL(INOEL,ICRYS,1,1) = Q_CRYS(ICRYS,1,1)
              Q_ALL(INOEL,ICRYS,1,2) = Q_CRYS(ICRYS,1,2)
              Q_ALL(INOEL,ICRY0,1,3) = Q_CRYS(ICRYS,1,3)
              Q_ALL(INOEL,ICRYS,2,1) = Q_CRYS(ICRYS,2,1)
              Q_ALL(INOEL,ICRYS,2,2) = Q_CRYS(ICRYS,2,2)
              Q_ALL(INOEL,ICRYS,2,3) = Q_CRYS(ICRYS,2,3)
              Q_ALL(INOEL,ICRYS,3,1) = Q_CRYS(ICRYS,3,1)
              Q_ALL(INOEL,ICRYS,3,2) = Q_CRYS(ICRYS,3,2)
              Q_ALL(INOEL,ICRYS,3,3) = Q_CRYS(ICRYS,3,3)
              VOL_FRAC_ALL(INOEL,ICRYS) = VOL(ICRYS)/VOLUME
11            CONTINUE

          ELSE
C
C              Then there are different crystals in
C              different elements. Proceed accordingly.
C          
            DO 20 INOEL = 1,NOEL
              VOLUME = ZERO
            DO 20 ICRYS = 1,NCRYS

              READ(65,*) PH, TH, OM, EIG11, EIG22, EIG33, IPHASE
                 
                VOL(ICRYS) = ONE

              TH = TH*PI/180.
              PH = PH*PI/180.
              OM = OM*PI/180.
            EIGEN(INOEL,1) = EIG11
            EIGEN(INOEL,2) = EIG22
            EIGEN(INOEL,3) = EIG33
            PHASE(INOEL) = IPHASE
     
              VOLUME = VOLUME + VOL(ICRYS)

              CALL ROTMAT(TH,PH,OM,Q)

              Q_ALL(INOEL,ICRYS,1,1) = Q(1,1)
              Q_ALL(INOEL,ICRYS,1,2) = Q(1,2)
              Q_ALL(INOEL,ICRYS,1,3) = Q(1,3)
              Q_ALL(INOEL,ICRYS,2,1) = Q(2,1)
              Q_ALL(INOEL,ICRYS,2,2) = Q(2,2)
              Q_ALL(INOEL,ICRYS,2,3) = Q(2,3)
              Q_ALL(INOEL,ICRYS,3,1) = Q(3,1)
              Q_ALL(INOEL,ICRYS,3,2) = Q(3,2)
              Q_ALL(INOEL,ICRYS,3,3) = Q(3,3)
            VOL_FRAC_ALL(INOEL,ICRYS) = VOL(ICRYS)/VOLUME            
      
20            CONTINUE

          ENDIF
          CLOSE(UNIT=65)
C
C              Read the slip systems for the crystalline material
C              Store the slip plane normal in ANV_C, and the slip
C              direction in AMV_C
C
          OPEN(UNIT=70,FILE=FILE2,STATUS='UNKNOWN')
          READ(70,*)
          READ(70,*)
          READ(70,*)
          DO 30 ISLIP = 1,NSLIP
            READ(70,*) (ANV_C(ISLIP,J),J=1,3),(AMV_C(ISLIP,J),J=1,3)
30          CONTINUE
          CLOSE(UNIT=70)
C
C               Normalize the slip system vectors to unit magnitude
C
          DO 40 I = 1,NSLIP
          AMV_C_MAG = DSQRT(AMV_C(I,1)**TWO + AMV_C(I,2)**TWO +
     +                            AMV_C(I,3)**TWO)
            ANV_C_MAG = DSQRT(ANV_C(I,1)**TWO + ANV_C(I,2)**TWO +
     +                          ANV_C(I,3)**TWO)

            AMV_C(I,1) = AMV_C(I,1)/AMV_C_MAG
            AMV_C(I,2) = AMV_C(I,2)/AMV_C_MAG
            AMV_C(I,3) = AMV_C(I,3)/AMV_C_MAG

            ANV_C(I,1) = ANV_C(I,1)/ANV_C_MAG
            ANV_C(I,2) = ANV_C(I,2)/ANV_C_MAG
            ANV_C(I,3) = ANV_C(I,3)/ANV_C_MAG
40          CONTINUE
C
C            Calculate the components of the Schmid tensors in the crystal
C            basis and store them in SMAT_C = [S^{\alpha}_0] in equation
C            A1 of Kalidindi
C
          DO 45 ISLIP = 1,NSLIP
            SMAT_FCC(ISLIP,1,1) = AMV_C(ISLIP,1)*ANV_C(ISLIP,1)
            SMAT_FCC(ISLIP,1,2) = AMV_C(ISLIP,1)*ANV_C(ISLIP,2)
            SMAT_FCC(ISLIP,1,3) = AMV_C(ISLIP,1)*ANV_C(ISLIP,3)
            SMAT_FCC(ISLIP,2,1) = AMV_C(ISLIP,2)*ANV_C(ISLIP,1)
            SMAT_FCC(ISLIP,2,2) = AMV_C(ISLIP,2)*ANV_C(ISLIP,2)
            SMAT_FCC(ISLIP,2,3) = AMV_C(ISLIP,2)*ANV_C(ISLIP,3)
            SMAT_FCC(ISLIP,3,1) = AMV_C(ISLIP,3)*ANV_C(ISLIP,1)
            SMAT_FCC(ISLIP,3,2) = AMV_C(ISLIP,3)*ANV_C(ISLIP,2)
            SMAT_FCC(ISLIP,3,3) = AMV_C(ISLIP,3)*ANV_C(ISLIP,3)
            SMAT_BCC(ISLIP,1,1) = ANV_C(ISLIP,1)*AMV_C(ISLIP,1)
            SMAT_BCC(ISLIP,1,2) = ANV_C(ISLIP,1)*AMV_C(ISLIP,2)
            SMAT_BCC(ISLIP,1,3) = ANV_C(ISLIP,1)*AMV_C(ISLIP,3)
            SMAT_BCC(ISLIP,2,1) = ANV_C(ISLIP,2)*AMV_C(ISLIP,1)
            SMAT_BCC(ISLIP,2,2) = ANV_C(ISLIP,2)*AMV_C(ISLIP,2)
            SMAT_BCC(ISLIP,2,3) = ANV_C(ISLIP,2)*AMV_C(ISLIP,3)
            SMAT_BCC(ISLIP,3,1) = ANV_C(ISLIP,3)*AMV_C(ISLIP,1)
            SMAT_BCC(ISLIP,3,2) = ANV_C(ISLIP,3)*AMV_C(ISLIP,2)
            SMAT_BCC(ISLIP,3,3) = ANV_C(ISLIP,3)*AMV_C(ISLIP,3)


45          CONTINUE
C
C            Calculate the components of the Schmid tensors in the crystal
C            basis and store them in SMAT_0_ALL = [Q][S^{\alpha}_0][Q]^T in 
C            equation A1 of Kalidindi
C                   SMAT_FCC: FCC slip system
C                   SMAT_BCC: BCC slip system
      DO 46 INOEL = 1,NOEL
              IF (PHASE(INOEL).EQ.1) THEN
                DO 47 ISLIP = 1,NSLIP
                  DO 47 I = 1,3
                  DO 47 J = 1,3
                     SMAT_C(ISLIP,I,J)=SMAT_FCC(ISLIP,I,J)      
47                  CONTINUE
              ELSE
                DO 48 ISLIP = 1,NSLIP
                  DO 48 I = 1,3
                  DO 48 J = 1,3
                     SMAT_C(ISLIP,I,J)=SMAT_BCC(ISLIP,I,J)      
48                  CONTINUE
              END IF
      DO 46 ICRYS = 1,NCRYS
      DO 46 ISLIP = 1,NSLIP
          
       SMAT_0_ALL(INOEL,ICRYS,ISLIP,1,1) = 
     + Q_ALL(INOEL,ICRYS,1,1)*SMAT_C(ISLIP,1,1)*Q_ALL(INOEL,ICRYS,1,1) +
     + Q_ALL(INOEL,ICRYS,1,1)*SMAT_C(ISLIP,1,2)*Q_ALL(INOEL,ICRYS,1,2) +
     + Q_ALL(INOEL,ICRYS,1,1)*SMAT_C(ISLIP,1,3)*Q_ALL(INOEL,ICRYS,1,3) +
     + Q_ALL(INOEL,ICRYS,1,2)*SMAT_C(ISLIP,2,1)*Q_ALL(INOEL,ICRYS,1,1) +
     + Q_ALL(INOEL,ICRYS,1,2)*SMAT_C(ISLIP,2,2)*Q_ALL(INOEL,ICRYS,1,2) +
     + Q_ALL(INOEL,ICRYS,1,2)*SMAT_C(ISLIP,2,3)*Q_ALL(INOEL,ICRYS,1,3) +
     + Q_ALL(INOEL,ICRYS,1,3)*SMAT_C(ISLIP,3,1)*Q_ALL(INOEL,ICRYS,1,1) +
     + Q_ALL(INOEL,ICRYS,1,3)*SMAT_C(ISLIP,3,2)*Q_ALL(INOEL,ICRYS,1,2) +
     + Q_ALL(INOEL,ICRYS,1,3)*SMAT_C(ISLIP,3,3)*Q_ALL(INOEL,ICRYS,1,3) 
      
       SMAT_0_ALL(INOEL,ICRYS,ISLIP,1,2) = 
     + Q_ALL(INOEL,ICRYS,1,1)*SMAT_C(ISLIP,1,1)*Q_ALL(INOEL,ICRYS,2,1) +
     + Q_ALL(INOEL,ICRYS,1,1)*SMAT_C(ISLIP,1,2)*Q_ALL(INOEL,ICRYS,2,2) +
     + Q_ALL(INOEL,ICRYS,1,1)*SMAT_C(ISLIP,1,3)*Q_ALL(INOEL,ICRYS,2,3) +
     + Q_ALL(INOEL,ICRYS,1,2)*SMAT_C(ISLIP,2,1)*Q_ALL(INOEL,ICRYS,2,1) +
     + Q_ALL(INOEL,ICRYS,1,2)*SMAT_C(ISLIP,2,2)*Q_ALL(INOEL,ICRYS,2,2) +
     + Q_ALL(INOEL,ICRYS,1,2)*SMAT_C(ISLIP,2,3)*Q_ALL(INOEL,ICRYS,2,3) +
     + Q_ALL(INOEL,ICRYS,1,3)*SMAT_C(ISLIP,3,1)*Q_ALL(INOEL,ICRYS,2,1) +
     + Q_ALL(INOEL,ICRYS,1,3)*SMAT_C(ISLIP,3,2)*Q_ALL(INOEL,ICRYS,2,2) +
     + Q_ALL(INOEL,ICRYS,1,3)*SMAT_C(ISLIP,3,3)*Q_ALL(INOEL,ICRYS,2,3) 
    
       SMAT_0_ALL(INOEL,ICRYS,ISLIP,1,3) = 
     + Q_ALL(INOEL,ICRYS,1,1)*SMAT_C(ISLIP,1,1)*Q_ALL(INOEL,ICRYS,3,1) +
     + Q_ALL(INOEL,ICRYS,1,1)*SMAT_C(ISLIP,1,2)*Q_ALL(INOEL,ICRYS,3,2) +
     + Q_ALL(INOEL,ICRYS,1,1)*SMAT_C(ISLIP,1,3)*Q_ALL(INOEL,ICRYS,3,3) +
     + Q_ALL(INOEL,ICRYS,1,2)*SMAT_C(ISLIP,2,1)*Q_ALL(INOEL,ICRYS,3,1) +
     + Q_ALL(INOEL,ICRYS,1,2)*SMAT_C(ISLIP,2,2)*Q_ALL(INOEL,ICRYS,3,2) +
     + Q_ALL(INOEL,ICRYS,1,2)*SMAT_C(ISLIP,2,3)*Q_ALL(INOEL,ICRYS,3,3) +
     + Q_ALL(INOEL,ICRYS,1,3)*SMAT_C(ISLIP,3,1)*Q_ALL(INOEL,ICRYS,3,1) +
     + Q_ALL(INOEL,ICRYS,1,3)*SMAT_C(ISLIP,3,2)*Q_ALL(INOEL,ICRYS,3,2) +
     + Q_ALL(INOEL,ICRYS,1,3)*SMAT_C(ISLIP,3,3)*Q_ALL(INOEL,ICRYS,3,3)
          
       SMAT_0_ALL(INOEL,ICRYS,ISLIP,2,1) = 
     + Q_ALL(INOEL,ICRYS,2,1)*SMAT_C(ISLIP,1,1)*Q_ALL(INOEL,ICRYS,1,1) +
     + Q_ALL(INOEL,ICRYS,2,1)*SMAT_C(ISLIP,1,2)*Q_ALL(INOEL,ICRYS,1,2) +
     + Q_ALL(INOEL,ICRYS,2,1)*SMAT_C(ISLIP,1,3)*Q_ALL(INOEL,ICRYS,1,3) +
     + Q_ALL(INOEL,ICRYS,2,2)*SMAT_C(ISLIP,2,1)*Q_ALL(INOEL,ICRYS,1,1) +
     + Q_ALL(INOEL,ICRYS,2,2)*SMAT_C(ISLIP,2,2)*Q_ALL(INOEL,ICRYS,1,2) +
     + Q_ALL(INOEL,ICRYS,2,2)*SMAT_C(ISLIP,2,3)*Q_ALL(INOEL,ICRYS,1,3) +
     + Q_ALL(INOEL,ICRYS,2,3)*SMAT_C(ISLIP,3,1)*Q_ALL(INOEL,ICRYS,1,1) +
     + Q_ALL(INOEL,ICRYS,2,3)*SMAT_C(ISLIP,3,2)*Q_ALL(INOEL,ICRYS,1,2) +
     + Q_ALL(INOEL,ICRYS,2,3)*SMAT_C(ISLIP,3,3)*Q_ALL(INOEL,ICRYS,1,3)
    
       SMAT_0_ALL(INOEL,ICRYS,ISLIP,2,2) = 
     + Q_ALL(INOEL,ICRYS,2,1)*SMAT_C(ISLIP,1,1)*Q_ALL(INOEL,ICRYS,2,1) +
     + Q_ALL(INOEL,ICRYS,2,1)*SMAT_C(ISLIP,1,2)*Q_ALL(INOEL,ICRYS,2,2) +
     + Q_ALL(INOEL,ICRYS,2,1)*SMAT_C(ISLIP,1,3)*Q_ALL(INOEL,ICRYS,2,3) +
     + Q_ALL(INOEL,ICRYS,2,2)*SMAT_C(ISLIP,2,1)*Q_ALL(INOEL,ICRYS,2,1) +
     + Q_ALL(INOEL,ICRYS,2,2)*SMAT_C(ISLIP,2,2)*Q_ALL(INOEL,ICRYS,2,2) +
     + Q_ALL(INOEL,ICRYS,2,2)*SMAT_C(ISLIP,2,3)*Q_ALL(INOEL,ICRYS,2,3) +
     + Q_ALL(INOEL,ICRYS,2,3)*SMAT_C(ISLIP,3,1)*Q_ALL(INOEL,ICRYS,2,1) +
     + Q_ALL(INOEL,ICRYS,2,3)*SMAT_C(ISLIP,3,2)*Q_ALL(INOEL,ICRYS,2,2) +
     + Q_ALL(INOEL,ICRYS,2,3)*SMAT_C(ISLIP,3,3)*Q_ALL(INOEL,ICRYS,2,3)
    
       SMAT_0_ALL(INOEL,ICRYS,ISLIP,2,3) = 
     + Q_ALL(INOEL,ICRYS,2,1)*SMAT_C(ISLIP,1,1)*Q_ALL(INOEL,ICRYS,3,1) +
     + Q_ALL(INOEL,ICRYS,2,1)*SMAT_C(ISLIP,1,2)*Q_ALL(INOEL,ICRYS,3,2) +
     + Q_ALL(INOEL,ICRYS,2,1)*SMAT_C(ISLIP,1,3)*Q_ALL(INOEL,ICRYS,3,3) +
     + Q_ALL(INOEL,ICRYS,2,2)*SMAT_C(ISLIP,2,1)*Q_ALL(INOEL,ICRYS,3,1) +
     + Q_ALL(INOEL,ICRYS,2,2)*SMAT_C(ISLIP,2,2)*Q_ALL(INOEL,ICRYS,3,2) +
     + Q_ALL(INOEL,ICRYS,2,2)*SMAT_C(ISLIP,2,3)*Q_ALL(INOEL,ICRYS,3,3) +
     + Q_ALL(INOEL,ICRYS,2,3)*SMAT_C(ISLIP,3,1)*Q_ALL(INOEL,ICRYS,3,1) +
     + Q_ALL(INOEL,ICRYS,2,3)*SMAT_C(ISLIP,3,2)*Q_ALL(INOEL,ICRYS,3,2) +
     + Q_ALL(INOEL,ICRYS,2,3)*SMAT_C(ISLIP,3,3)*Q_ALL(INOEL,ICRYS,3,3)
    
       SMAT_0_ALL(INOEL,ICRYS,ISLIP,3,1) = 
     + Q_ALL(INOEL,ICRYS,3,1)*SMAT_C(ISLIP,1,1)*Q_ALL(INOEL,ICRYS,1,1) +
     + Q_ALL(INOEL,ICRYS,3,1)*SMAT_C(ISLIP,1,2)*Q_ALL(INOEL,ICRYS,1,2) +
     + Q_ALL(INOEL,ICRYS,3,1)*SMAT_C(ISLIP,1,3)*Q_ALL(INOEL,ICRYS,1,3) +
     + Q_ALL(INOEL,ICRYS,3,2)*SMAT_C(ISLIP,2,1)*Q_ALL(INOEL,ICRYS,1,1) +
     + Q_ALL(INOEL,ICRYS,3,2)*SMAT_C(ISLIP,2,2)*Q_ALL(INOEL,ICRYS,1,2) +
     + Q_ALL(INOEL,ICRYS,3,2)*SMAT_C(ISLIP,2,3)*Q_ALL(INOEL,ICRYS,1,3) +
     + Q_ALL(INOEL,ICRYS,3,3)*SMAT_C(ISLIP,3,1)*Q_ALL(INOEL,ICRYS,1,1) +
     + Q_ALL(INOEL,ICRYS,3,3)*SMAT_C(ISLIP,3,2)*Q_ALL(INOEL,ICRYS,1,2) +
     + Q_ALL(INOEL,ICRYS,3,3)*SMAT_C(ISLIP,3,3)*Q_ALL(INOEL,ICRYS,1,3)
    
       SMAT_0_ALL(INOEL,ICRYS,ISLIP,3,2) = 
     + Q_ALL(INOEL,ICRYS,3,1)*SMAT_C(ISLIP,1,1)*Q_ALL(INOEL,ICRYS,2,1) +
     + Q_ALL(INOEL,ICRYS,3,1)*SMAT_C(ISLIP,1,2)*Q_ALL(INOEL,ICRYS,2,2) +
     + Q_ALL(INOEL,ICRYS,3,1)*SMAT_C(ISLIP,1,3)*Q_ALL(INOEL,ICRYS,2,3) +
     + Q_ALL(INOEL,ICRYS,3,2)*SMAT_C(ISLIP,2,1)*Q_ALL(INOEL,ICRYS,2,1) +
     + Q_ALL(INOEL,ICRYS,3,2)*SMAT_C(ISLIP,2,2)*Q_ALL(INOEL,ICRYS,2,2) +
     + Q_ALL(INOEL,ICRYS,3,2)*SMAT_C(ISLIP,2,3)*Q_ALL(INOEL,ICRYS,2,3) +
     + Q_ALL(INOEL,ICRYS,3,3)*SMAT_C(ISLIP,3,1)*Q_ALL(INOEL,ICRYS,2,1) +
     + Q_ALL(INOEL,ICRYS,3,3)*SMAT_C(ISLIP,3,2)*Q_ALL(INOEL,ICRYS,2,2) +
     + Q_ALL(INOEL,ICRYS,3,3)*SMAT_C(ISLIP,3,3)*Q_ALL(INOEL,ICRYS,2,3)
    
       SMAT_0_ALL(INOEL,ICRYS,ISLIP,3,3) = 
     + Q_ALL(INOEL,ICRYS,3,1)*SMAT_C(ISLIP,1,1)*Q_ALL(INOEL,ICRYS,3,1) +
     + Q_ALL(INOEL,ICRYS,3,1)*SMAT_C(ISLIP,1,2)*Q_ALL(INOEL,ICRYS,3,2) +
     + Q_ALL(INOEL,ICRYS,3,1)*SMAT_C(ISLIP,1,3)*Q_ALL(INOEL,ICRYS,3,3) +
     + Q_ALL(INOEL,ICRYS,3,2)*SMAT_C(ISLIP,2,1)*Q_ALL(INOEL,ICRYS,3,1) +
     + Q_ALL(INOEL,ICRYS,3,2)*SMAT_C(ISLIP,2,2)*Q_ALL(INOEL,ICRYS,3,2) +
     + Q_ALL(INOEL,ICRYS,3,2)*SMAT_C(ISLIP,2,3)*Q_ALL(INOEL,ICRYS,3,3) +
     + Q_ALL(INOEL,ICRYS,3,3)*SMAT_C(ISLIP,3,1)*Q_ALL(INOEL,ICRYS,3,1) +
     + Q_ALL(INOEL,ICRYS,3,3)*SMAT_C(ISLIP,3,2)*Q_ALL(INOEL,ICRYS,3,2) +
     + Q_ALL(INOEL,ICRYS,3,3)*SMAT_C(ISLIP,3,3)*Q_ALL(INOEL,ICRYS,3,3)

46      CONTINUE 

C                  Material property based on phase ID
        IF (PHASE(INOEL).EQ.1) THEN   
          S0 = S0_FCC
          QL = QL_FCC
        ELSE
          S0 = S0_BCC
          QL = QL_BCC
        END IF

C
C           The calculation procedure needs the inverse of the plastic
C           deformation gradient for all elements and all crystals. This 
C           is stored and updated in a  matrix of internal variables.
C           Initialize these state variables at the beginning of a 
C           calculation.
C        
      DO 60 INOEL = 1,NOEL
      DO 60 ICRYS = 1,NCRYS
            FP_T_INV_ALL(INOEL,ICRYS,1,1) = ONE
            FP_T_INV_ALL(INOEL,ICRYS,1,2) = ZERO
            FP_T_INV_ALL(INOEL,ICRYS,1,3) = ZERO
            FP_T_INV_ALL(INOEL,ICRYS,2,1) = ZERO
            FP_T_INV_ALL(INOEL,ICRYS,2,2) = ONE
            FP_T_INV_ALL(INOEL,ICRYS,2,3) = ZERO
            FP_T_INV_ALL(INOEL,ICRYS,3,1) = ZERO
            FP_T_INV_ALL(INOEL,ICRYS,3,2) = ZERO
            FP_T_INV_ALL(INOEL,ICRYS,3,3) = ONE
C
C             Initailize the slip system deformation resistances, 
C             the shear increments, and the accumulated shears
C             for all slip systems in every crystal.
C              
           DO 55 ISLIP = 1,NSLIP
           S_T_ALL(INOEL,ICRYS,ISLIP)       = S0
             DGAMMA_T_ALL(INOEL,ICRYS,ISLIP)  = ZERO
           ACCGAM_T_ALL(INOEL,ICRYS,ISLIP)  = ZERO
55           CONTINUE
60      CONTINUE
C
C            Initialize the latent hardening matrix
C
      DO 64 ISLIP = 1,NSLIP
      DO 64 JSLIP = 1,NSLIP
        IF (ISLIP .NE. JSLIP) THEN
              QLAT(ISLIP,JSLIP) = QL
          ELSE
              QLAT(ISLIP,JSLIP) = ONE
        ENDIF
        IF (ISLIP .NE. JSLIP) QLAT(JSLIP,ISLIP) = QLAT(ISLIP,JSLIP)
64    CONTINUE

C         FILE3 = '/usr/home/cricket/tzhu/SXAL_VUMAT/VUMAT_MESSAGES'
C         FILE4 = '/usr/home/cricket/tzhu/SXAL_VUMAT/VUMAT_TEXTURES'
C         FILE5 = '/usr/home/cricket/tzhu/SXAL_VUMAT/VUMAT_SLIPS'

          OPEN(UNIT=80,FILE=FILE3,STATUS='UNKNOWN')
           OPEN(UNIT=90,FILE=FILE4,STATUS='UNKNOWN')
           OPEN(UNIT=100,FILE=FILE5,STATUS='UNKNOWN')
          
      ENDIF
C**************************************************************************
C              End of initialization. Proceed with computation.
C**************************************************************************

C                  Material property based on phase ID
      IF (PHASE(NOEL_PT+KM).EQ.1) THEN   
         C11   = C11_FCC
         C12   = C12_FCC
         C44   = C44_FCC
         GDOT0 = GDOT0_FCC
         AM    = AM_FCC
         S0    = S0_FCC
         AH0   = AH0_FCC
         SSAT  = SSAT_FCC
         RHARD = RHARD_FCC
         QL    = QL_FCC
      ELSE
         C11   = C11_BCC
         C12   = C12_BCC
         C44   = C44_BCC
         GDOT0 = GDOT0_BCC
         AM    = AM_BCC
         S0    = S0_BCC
         AH0   = AH0_BCC
         SSAT  = SSAT_BCC
         RHARD = RHARD_BCC
         QL    = QL_BCC
      END IF

C
C      Copy the old and new Deformation gradients into F_T and F_TAU 
C      respectively
C         write(*,*) 'RHARD',am,S0,ah0,SSAT
      
      F_T(1,1) = DEFGRAD_OLD(KM,1)
      F_T(2,2) = DEFGRAD_OLD(KM,2)
      F_T(3,3) = DEFGRAD_OLD(KM,3)
      F_T(1,2) = DEFGRAD_OLD(KM,4)

      F_TAU(1,1) = DEFGRAD_NEW(KM,1)
      F_TAU(2,2) = DEFGRAD_NEW(KM,2)
      F_TAU(3,3) = DEFGRAD_NEW(KM,3)
      F_TAU(1,2) = DEFGRAD_NEW(KM,4)

      U_TAU(1,1) = STRETCH_NEW(KM,1)
      U_TAU(2,2) = STRETCH_NEW(KM,2)
      U_TAU(3,3) = STRETCH_NEW(KM,3)
      U_TAU(1,2) = STRETCH_NEW(KM,4)
        
      IF (NSHR .EQ. 1) THEN
        F_T(2,3) = ZERO
        F_T(3,1) = ZERO
        F_T(2,1) = DEFGRAD_OLD(KM,5)
        F_T(3,2) = ZERO
        F_T(1,3) = ZERO

        F_TAU(2,3) = ZERO
        F_TAU(3,1) = ZERO
        F_TAU(2,1) = DEFGRAD_NEW(KM,5)
        F_TAU(3,2) = ZERO
        F_TAU(1,3) = ZERO

        U_TAU(2,3) = ZERO
        U_TAU(3,1) = ZERO
        U_TAU(2,1) = U_TAU(1,2)
        U_TAU(3,2) = ZERO
        U_TAU(1,3) = ZERO
      ELSE
        F_T(2,3) = DEFGRAD_OLD(KM,5)
        F_T(3,1) = DEFGRAD_OLD(KM,6)
        F_T(2,1) = DEFGRAD_OLD(KM,7)
        F_T(3,2) = DEFGRAD_OLD(KM,8)
        F_T(1,3) = DEFGRAD_OLD(KM,9)

        F_TAU(2,3) = DEFGRAD_NEW(KM,5)
        F_TAU(3,1) = DEFGRAD_NEW(KM,6)
        F_TAU(2,1) = DEFGRAD_NEW(KM,7)
        F_TAU(3,2) = DEFGRAD_NEW(KM,8)
        F_TAU(1,3) = DEFGRAD_NEW(KM,9)

        U_TAU(2,3) = STRETCH_NEW(KM,5)
        U_TAU(3,1) = STRETCH_NEW(KM,6)
        U_TAU(2,1) = U_TAU(1,2)
        U_TAU(3,2) = U_TAU(2,3)
        U_TAU(1,3) = U_TAU(3,1)
       ENDIF

C
C            CALCULATE THE ELASTIC STIFFNESS MATRIX
C         FOR EACH CRYSTAL
C
          
      DO 66 I = 1,3
      DO 66 J = 1,3
      DO 66 K = 1,3
      DO 66 L = 1,3
            ELMAT_C(I,J,K,L) = ZERO
66      CONTINUE

          ELMAT_C(1,1,1,1) = C11
          ELMAT_C(1,1,2,2) = C12
          ELMAT_C(1,1,3,3) = C12
          ELMAT_C(2,2,1,1) = C12
          ELMAT_C(2,2,2,2) = C11
          ELMAT_C(2,2,3,3) = C12
          ELMAT_C(3,3,1,1) = C12
          ELMAT_C(3,3,2,2) = C12
          ELMAT_C(3,3,3,3) = C11
          ELMAT_C(1,2,1,2) = C44
          ELMAT_C(1,2,2,1) = C44
          ELMAT_C(2,1,1,2) = C44
          ELMAT_C(2,1,2,1) = C44
          ELMAT_C(1,3,1,3) = C44
          ELMAT_C(1,3,3,1) = C44
          ELMAT_C(3,1,1,3) = C44
          ELMAT_C(3,1,3,1) = C44
          ELMAT_C(2,3,2,3) = C44
          ELMAT_C(2,3,3,2) = C44
          ELMAT_C(3,2,2,3) = C44
          ELMAT_C(3,2,3,2) = C44
C
C          Initialize the volume averaged Cauchy stress
C          {\bar T}_{\tau} at the end of the step
C
            
        TBAR_TAU(1,1) = ZERO
        TBAR_TAU(1,2) = ZERO
        TBAR_TAU(1,3) = ZERO
        TBAR_TAU(2,1) = ZERO
        TBAR_TAU(2,2) = ZERO
        TBAR_TAU(2,3) = ZERO
        TBAR_TAU(3,1) = ZERO
        TBAR_TAU(3,2) = ZERO
        TBAR_TAU(3,3) = ZERO
C
C          Initialize the plastic work increment 
c
        PLASTIC_WORK_INC = ZERO
C
C        Start the calculation for the stress {\bar T}_{\tau}
C        and the update of other variables
C
C******************************************************************
C          Begin loop over all crystals at an integration point
C******************************************************************
        
  
      
      DO 500 ICRYS = 1,NCRYS
C
C              First read in the stored value of {\fpt}^{-1}
C              for the particular crystal at an integration point
C              in an element
C
            
          FP_T_INV(1,1) = FP_T_INV_ALL(NOEL_PT+KM,ICRYS,1,1)
          FP_T_INV(1,2) = FP_T_INV_ALL(NOEL_PT+KM,ICRYS,1,2)
          FP_T_INV(1,3) = FP_T_INV_ALL(NOEL_PT+KM,ICRYS,1,3)
          FP_T_INV(2,1) = FP_T_INV_ALL(NOEL_PT+KM,ICRYS,2,1)
          FP_T_INV(2,2) = FP_T_INV_ALL(NOEL_PT+KM,ICRYS,2,2)
          FP_T_INV(2,3) = FP_T_INV_ALL(NOEL_PT+KM,ICRYS,2,3)
          FP_T_INV(3,1) = FP_T_INV_ALL(NOEL_PT+KM,ICRYS,3,1)
          FP_T_INV(3,2) = FP_T_INV_ALL(NOEL_PT+KM,ICRYS,3,2)
          FP_T_INV(3,3) = FP_T_INV_ALL(NOEL_PT+KM,ICRYS,3,3)
           
           
C
C              Read in the stored value of the orientation matrix [Q] 
C              for the particular crystal at an integration point
C              in an element
C
          Q(1,1) = Q_ALL(NOEL_PT+KM,ICRYS,1,1)
          Q(1,2) = Q_ALL(NOEL_PT+KM,ICRYS,1,2)
          Q(1,3) = Q_ALL(NOEL_PT+KM,ICRYS,1,3)
          Q(2,1) = Q_ALL(NOEL_PT+KM,ICRYS,2,1)
          Q(2,2) = Q_ALL(NOEL_PT+KM,ICRYS,2,2)
          Q(2,3) = Q_ALL(NOEL_PT+KM,ICRYS,2,3)
          Q(3,1) = Q_ALL(NOEL_PT+KM,ICRYS,3,1)
          Q(3,2) = Q_ALL(NOEL_PT+KM,ICRYS,3,2)
          Q(3,3) = Q_ALL(NOEL_PT+KM,ICRYS,3,3)

C
C              Calculate the elastic stiffness matrix
C              in the global basis
C          
          DO 68 I = 1,3
          DO 68 J = 1,3
          DO 68 K = 1,3
          DO 68 L = 1,3
            ELMAT_G_ALL(NOEL_PT+KM,ICRYS,I,J,K,L) = ZERO
            DO 68 II = 1,3
            DO 68 JJ = 1,3
            DO 68 KK = 1,3
            DO 68 LL = 1,3
              ELMAT_G_ALL(NOEL_PT+KM,ICRYS,I,J,K,L) = 
     +                 ELMAT_G_ALL(NOEL_PT+KM,ICRYS,I,J,K,L) +
     +                 Q(I,II)*Q(K,KK)*
     +               Q(L,LL)*Q(J,JJ)*
     +               ELMAT_C(II,JJ,KK,LL)
68          CONTINUE         
C          write(*,*) ELMAT_G_ALL(NOEL_PT+KM,ICRYS,1,1,1,1)
C          write(*,*) ELMAT_G_ALL(NOEL_PT+KM,ICRYS,1,2,1,2)
C
C          Read in the stored value of the slip system resistances
C          and the shear strain increments at the beginning of the
C          step for the slip systems of  a particular crystal 
C          at an integration point in an element
C
   
      
          DO 69 ISLIP = 1,NSLIP
            S_T(ISLIP)      = S_T_ALL(NOEL_PT+KM,ICRYS,ISLIP)
            DGAMMA_T(ISLIP) = DGAMMA_T_ALL(NOEL_PT+KM,ICRYS,ISLIP)
              ACCGAM_T(ISLIP) = ACCGAM_T_ALL(NOEL_PT+KM,ICRYS,ISLIP) 
69          CONTINUE
           
C
C           Calculate {\fp_tau}^{-1} for the crystal.
C            See equation 20 of Kalidindi et al., 1992.
C
          DO 75 I = 1,3
          DO 75 J = 1,3
            AUX1(I,J) = ZERO
            DO 70 ISLIP = 1,NSLIP
              AUX1(I,J) = AUX1(I,J) +
     +          DGAMMA_T(ISLIP)*SMAT_0_ALL(NOEL_PT+KM,ICRYS,ISLIP,I,J)
70            CONTINUE
            AUX2(I,J) = AIDENT(I,J) - AUX1(I,J)
75          CONTINUE

          FP_TAU_INV(1,1) = FP_T_INV(1,1)*AUX2(1,1) +
     +                        FP_T_INV(1,2)*AUX2(2,1) +
     +                        FP_T_INV(1,3)*AUX2(3,1)

          FP_TAU_INV(1,2) = FP_T_INV(1,1)*AUX2(1,2) +
     +                        FP_T_INV(1,2)*AUX2(2,2) +
     +                        FP_T_INV(1,3)*AUX2(3,2)

          FP_TAU_INV(1,3) = FP_T_INV(1,1)*AUX2(1,3) +
     +                        FP_T_INV(1,2)*AUX2(2,3) +
     +                        FP_T_INV(1,3)*AUX2(3,3)

          FP_TAU_INV(2,1) = FP_T_INV(2,1)*AUX2(1,1) +
     +                        FP_T_INV(2,2)*AUX2(2,1) +
     +                        FP_T_INV(2,3)*AUX2(3,1)

          FP_TAU_INV(2,2) = FP_T_INV(2,1)*AUX2(1,2) +
     +                        FP_T_INV(2,2)*AUX2(2,2) +
     +                        FP_T_INV(2,3)*AUX2(3,2)

          FP_TAU_INV(2,3) = FP_T_INV(2,1)*AUX2(1,3) +
     +                        FP_T_INV(2,2)*AUX2(2,3) +
     +                        FP_T_INV(2,3)*AUX2(3,3)

          FP_TAU_INV(3,1) = FP_T_INV(3,1)*AUX2(1,1) +
     +                        FP_T_INV(3,2)*AUX2(2,1) +
     +                        FP_T_INV(3,3)*AUX2(3,1)

          FP_TAU_INV(3,2) = FP_T_INV(3,1)*AUX2(1,2) +
     +                        FP_T_INV(3,2)*AUX2(2,2) +
     +                        FP_T_INV(3,3)*AUX2(3,2)

          FP_TAU_INV(3,3) = FP_T_INV(3,1)*AUX2(1,3) +
     +                        FP_T_INV(3,2)*AUX2(2,3) +
     +                        FP_T_INV(3,3)*AUX2(3,3)
C
C            Calculate the det of {\fp_tau}^{-1}
C
          DET_FP_TAU_INV = 
     +      FP_TAU_INV(1,1)*(FP_TAU_INV(2,2)*FP_TAU_INV(3,3)-
     +                         FP_TAU_INV(3,2)*FP_TAU_INV(2,3)) -
     +          FP_TAU_INV(2,1)*(FP_TAU_INV(1,2)*FP_TAU_INV(3,3)-
     +                       FP_TAU_INV(3,2)*FP_TAU_INV(1,3)) +
     +          FP_TAU_INV(3,1)*(FP_TAU_INV(1,2)*FP_TAU_INV(2,3)-
     +                       FP_TAU_INV(2,2)*FP_TAU_INV(1,3))
      
            
      
          IF (DET_FP_TAU_INV .LE. ZERO) THEN
c            WRITE(80,*) 'WARNING: MAIN:'
c            WRITE(80,*) 'TOTAL_TIME = ',TOTAL_TIME,' DTIME = ',DTIME
c            WRITE(80,*) 'INOEL = ',NOEL_PT+KM,' ICRYS = ',ICRYS
c            WRITE(80,*) 'WARNING: DET_FP_TAU_INV :-'
c            WRITE(80,*) 'WARNING: DET_FP_TAU_INV = ',DET_FP_TAU_INV
            CALL PRTMAT(FP_TAU_INV,3,3)
          ENDIF
            
C
C              Normalize {\fp_tau}^{-1} such that its determinant
C              is equal to unity.
C
            
          VAL = ONE/(DET_FP_TAU_INV**ONE_THIRD)
          FP_TAU_INV(1,1) = FP_TAU_INV(1,1)*VAL
          FP_TAU_INV(1,2) = FP_TAU_INV(1,2)*VAL
          FP_TAU_INV(1,3) = FP_TAU_INV(1,3)*VAL
          FP_TAU_INV(2,1) = FP_TAU_INV(2,1)*VAL
          FP_TAU_INV(2,2) = FP_TAU_INV(2,2)*VAL
          FP_TAU_INV(2,3) = FP_TAU_INV(2,3)*VAL
          FP_TAU_INV(3,1) = FP_TAU_INV(3,1)*VAL
          FP_TAU_INV(3,2) = FP_TAU_INV(3,2)*VAL
          FP_TAU_INV(3,3) = FP_TAU_INV(3,3)*VAL
C
C           Calculate the elastic deformation gradient
C           F^*(\tau) = F(\tau)*F^p^{-1}(\tau)
C          Equation 5 of Kalidindi
C

          FSTAR_TAU(1,1) = F_TAU(1,1)*FP_TAU_INV(1,1) +
     +                       F_TAU(1,2)*FP_TAU_INV(2,1) +
     +                       F_TAU(1,3)*FP_TAU_INV(3,1)

          FSTAR_TAU(1,2) = F_TAU(1,1)*FP_TAU_INV(1,2) +
     +                       F_TAU(1,2)*FP_TAU_INV(2,2) +
     +                       F_TAU(1,3)*FP_TAU_INV(3,2)

          FSTAR_TAU(1,3) = F_TAU(1,1)*FP_TAU_INV(1,3) +
     +                       F_TAU(1,2)*FP_TAU_INV(2,3) +
     +                       F_TAU(1,3)*FP_TAU_INV(3,3)

          FSTAR_TAU(2,1) = F_TAU(2,1)*FP_TAU_INV(1,1) +
     +                       F_TAU(2,2)*FP_TAU_INV(2,1) +
     +                       F_TAU(2,3)*FP_TAU_INV(3,1)

          FSTAR_TAU(2,2) = F_TAU(2,1)*FP_TAU_INV(1,2) +
     +                       F_TAU(2,2)*FP_TAU_INV(2,2) +
     +                       F_TAU(2,3)*FP_TAU_INV(3,2)

          FSTAR_TAU(2,3) = F_TAU(2,1)*FP_TAU_INV(1,3) +
     +                       F_TAU(2,2)*FP_TAU_INV(2,3) +
     +                       F_TAU(2,3)*FP_TAU_INV(3,3)

          FSTAR_TAU(3,1) = F_TAU(3,1)*FP_TAU_INV(1,1) +
     +                       F_TAU(3,2)*FP_TAU_INV(2,1) +
     +                       F_TAU(3,3)*FP_TAU_INV(3,1)

          FSTAR_TAU(3,2) = F_TAU(3,1)*FP_TAU_INV(1,2) +
     +                       F_TAU(3,2)*FP_TAU_INV(2,2) +
     +                       F_TAU(3,3)*FP_TAU_INV(3,2)

          FSTAR_TAU(3,3) = F_TAU(3,1)*FP_TAU_INV(1,3) +
     +                       F_TAU(3,2)*FP_TAU_INV(2,3) +
     +                       F_TAU(3,3)*FP_TAU_INV(3,3)
C
C              Calculate Det of F^{*}(\tau)
C
          DET_FSTAR_TAU = 
     +          FSTAR_TAU(1,1)*(FSTAR_TAU(2,2)*FSTAR_TAU(3,3)-
     +                        FSTAR_TAU(3,2)*FSTAR_TAU(2,3)) -
     +          FSTAR_TAU(2,1)*(FSTAR_TAU(1,2)*FSTAR_TAU(3,3)-
     +                      FSTAR_TAU(3,2)*FSTAR_TAU(1,3)) +
     +          FSTAR_TAU(3,1)*(FSTAR_TAU(1,2)*FSTAR_TAU(2,3)-
     +                      FSTAR_TAU(2,2)*FSTAR_TAU(1,3))

          IF (DET_FSTAR_TAU .LE. ZERO) THEN
            WRITE(80,*) 'WARNING: MAIN:'
            WRITE(80,*) 'WARNING: DET_FSTAR_TAU = ',DET_FSTAR_TAU
            WRITE(80,*) 'INOEL = ',NOEL_PT+KM,' ICRYS = ',ICRYS
            WRITE(80,*) 'WARNING: DET_FSTAR_TAU :-'
            CALL PRTMAT(FSTAR_TAU,3,3)
          ENDIF
            
      
          DET_FSTAR_TAU_INV = ONE/DET_FSTAR_TAU
C
C          Calculate the elastic strain {E^*}_{\tau}
C         Equation 3 of Kalidindi
C           minus eigenstrain, modified by Yin Zhang 27th Nov 2017
           
          ESTAR_TAU(1,1) = 
     +          ONE_HALF*(FSTAR_TAU(1,1)*FSTAR_TAU(1,1) +
     +                  FSTAR_TAU(2,1)*FSTAR_TAU(2,1) +
     +                  FSTAR_TAU(3,1)*FSTAR_TAU(3,1) - 
     +                  AIDENT(1,1))
C - EIGEN(NOEL_PT+KM,1)

          ESTAR_TAU(2,2) =
     +          ONE_HALF*(FSTAR_TAU(1,2)*FSTAR_TAU(1,2) +
     +                  FSTAR_TAU(2,2)*FSTAR_TAU(2,2) +
     +                    FSTAR_TAU(3,2)*FSTAR_TAU(3,2) - 
     +                  AIDENT(2,2))
C - EIGEN(NOEL_PT+KM,2)
     
          ESTAR_TAU(3,3) = 
     +          ONE_HALF*(FSTAR_TAU(1,3)*FSTAR_TAU(1,3) +
     +                   FSTAR_TAU(2,3)*FSTAR_TAU(2,3) +
     +                  FSTAR_TAU(3,3)*FSTAR_TAU(3,3) - 
     +                  AIDENT(3,3))
C - EIGEN(NOEL_PT+KM,3)
     
          ESTAR_TAU(1,2) = 
     +          ONE_HALF*(FSTAR_TAU(1,1)*FSTAR_TAU(1,2) +
     +                  FSTAR_TAU(2,1)*FSTAR_TAU(2,2) +
     +                  FSTAR_TAU(3,1)*FSTAR_TAU(3,2) - 
     +                  AIDENT(1,2))

          ESTAR_TAU(1,3) = 
     +          ONE_HALF*(FSTAR_TAU(1,1)*FSTAR_TAU(1,3) +
     +                  FSTAR_TAU(2,1)*FSTAR_TAU(2,3) +
     +                  FSTAR_TAU(3,1)*FSTAR_TAU(3,3) - 
     +                  AIDENT(1,3))

          ESTAR_TAU(2,3) = 
     +          ONE_HALF*(FSTAR_TAU(1,2)*FSTAR_TAU(1,3) +
     +                  FSTAR_TAU(2,2)*FSTAR_TAU(2,3) +
     +                  FSTAR_TAU(3,2)*FSTAR_TAU(3,3) - 
     +                  AIDENT(2,3))

          ESTAR_TAU(2,1) = ESTAR_TAU(1,2)
          ESTAR_TAU(3,1) = ESTAR_TAU(1,3)
          ESTAR_TAU(3,2) = ESTAR_TAU(2,3)


C          
C           Calculate the stress {T^*}_{\tau}
C            Equation 2 of Kalidindi
C
            TSTAR_TAU(1,1) = 
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,1,1,1,1)*ESTAR_TAU(1,1) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,1,1,1,2)*ESTAR_TAU(1,2) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,1,1,1,3)*ESTAR_TAU(1,3) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,1,1,2,1)*ESTAR_TAU(2,1) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,1,1,2,2)*ESTAR_TAU(2,2) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,1,1,2,3)*ESTAR_TAU(2,3) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,1,1,3,1)*ESTAR_TAU(3,1) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,1,1,3,2)*ESTAR_TAU(3,2) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,1,1,3,3)*ESTAR_TAU(3,3)

            TSTAR_TAU(2,2) = 
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,2,2,1,1)*ESTAR_TAU(1,1) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,2,2,1,2)*ESTAR_TAU(1,2) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,2,2,1,3)*ESTAR_TAU(1,3) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,2,2,2,1)*ESTAR_TAU(2,1) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,2,2,2,2)*ESTAR_TAU(2,2) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,2,2,2,3)*ESTAR_TAU(2,3) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,2,2,3,1)*ESTAR_TAU(3,1) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,2,2,3,2)*ESTAR_TAU(3,2) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,2,2,3,3)*ESTAR_TAU(3,3)

            TSTAR_TAU(3,3) = 
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,3,3,1,1)*ESTAR_TAU(1,1) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,3,3,1,2)*ESTAR_TAU(1,2) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,3,3,1,3)*ESTAR_TAU(1,3) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,3,3,2,1)*ESTAR_TAU(2,1) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,3,3,2,2)*ESTAR_TAU(2,2) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,3,3,2,3)*ESTAR_TAU(2,3) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,3,3,3,1)*ESTAR_TAU(3,1) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,3,3,3,2)*ESTAR_TAU(3,2) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,3,3,3,3)*ESTAR_TAU(3,3)
          
            TSTAR_TAU(1,2) = 
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,1,2,1,1)*ESTAR_TAU(1,1) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,1,2,1,2)*ESTAR_TAU(1,2) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,1,2,1,3)*ESTAR_TAU(1,3) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,1,2,2,1)*ESTAR_TAU(2,1) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,1,2,2,2)*ESTAR_TAU(2,2) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,1,2,2,3)*ESTAR_TAU(2,3) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,1,2,3,1)*ESTAR_TAU(3,1) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,1,2,3,2)*ESTAR_TAU(3,2) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,1,2,3,3)*ESTAR_TAU(3,3)

            TSTAR_TAU(1,3) = 
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,1,3,1,1)*ESTAR_TAU(1,1) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,1,3,1,2)*ESTAR_TAU(1,2) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,1,3,1,3)*ESTAR_TAU(1,3) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,1,3,2,1)*ESTAR_TAU(2,1) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,1,3,2,2)*ESTAR_TAU(2,2) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,1,3,2,3)*ESTAR_TAU(2,3) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,1,3,3,1)*ESTAR_TAU(3,1) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,1,3,3,2)*ESTAR_TAU(3,2) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,1,3,3,3)*ESTAR_TAU(3,3)

            TSTAR_TAU(2,3) =
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,2,3,1,1)*ESTAR_TAU(1,1) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,2,3,1,2)*ESTAR_TAU(1,2) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,2,3,1,3)*ESTAR_TAU(1,3) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,2,3,2,1)*ESTAR_TAU(2,1) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,2,3,2,2)*ESTAR_TAU(2,2) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,2,3,2,3)*ESTAR_TAU(2,3) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,2,3,3,1)*ESTAR_TAU(3,1) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,2,3,3,2)*ESTAR_TAU(3,2) +
     +      ELMAT_G_ALL(NOEL_PT+KM,ICRYS,2,3,3,3)*ESTAR_TAU(3,3)

            TSTAR_TAU(2,1) = TSTAR_TAU(1,2)
            TSTAR_TAU(3,1) = TSTAR_TAU(1,3)
            TSTAR_TAU(3,2) = TSTAR_TAU(2,3)

C          If(NOEL_PT+KM.eq.1) then
C           write(*,*) DGAMMA_T_ALL(NOEL_PT+KM,1,4)
C          endif
C
C             Calculate the Cauchy stress {T}_{\tau}
C             for each crystal
C             Use equation 4 of Kalidindi
C
          T_TAU(1,1) = 
     +          (FSTAR_TAU(1,1)*TSTAR_TAU(1,1)*FSTAR_TAU(1,1) +
     +           FSTAR_TAU(1,1)*TSTAR_TAU(1,2)*FSTAR_TAU(1,2) +
     +           FSTAR_TAU(1,1)*TSTAR_TAU(1,3)*FSTAR_TAU(1,3) +
     +           FSTAR_TAU(1,2)*TSTAR_TAU(2,1)*FSTAR_TAU(1,1) +
     +           FSTAR_TAU(1,2)*TSTAR_TAU(2,2)*FSTAR_TAU(1,2) +
     +           FSTAR_TAU(1,2)*TSTAR_TAU(2,3)*FSTAR_TAU(1,3) +
     +           FSTAR_TAU(1,3)*TSTAR_TAU(3,1)*FSTAR_TAU(1,1) +
     +           FSTAR_TAU(1,3)*TSTAR_TAU(3,2)*FSTAR_TAU(1,2) +
     +           FSTAR_TAU(1,3)*TSTAR_TAU(3,3)*FSTAR_TAU(1,3))*
     +        DET_FSTAR_TAU_INV
            
          T_TAU(2,2) = 
     +          (FSTAR_TAU(2,1)*TSTAR_TAU(1,1)*FSTAR_TAU(2,1) +
     +           FSTAR_TAU(2,1)*TSTAR_TAU(1,2)*FSTAR_TAU(2,2) +
     +           FSTAR_TAU(2,1)*TSTAR_TAU(1,3)*FSTAR_TAU(2,3) +
     +           FSTAR_TAU(2,2)*TSTAR_TAU(2,1)*FSTAR_TAU(2,1) +
     +           FSTAR_TAU(2,2)*TSTAR_TAU(2,2)*FSTAR_TAU(2,2) +
     +           FSTAR_TAU(2,2)*TSTAR_TAU(2,3)*FSTAR_TAU(2,3) +
     +           FSTAR_TAU(2,3)*TSTAR_TAU(3,1)*FSTAR_TAU(2,1) +
     +           FSTAR_TAU(2,3)*TSTAR_TAU(3,2)*FSTAR_TAU(2,2) +
     +           FSTAR_TAU(2,3)*TSTAR_TAU(3,3)*FSTAR_TAU(2,3))*
     +       DET_FSTAR_TAU_INV

          T_TAU(3,3) =
     +          (FSTAR_TAU(3,1)*TSTAR_TAU(1,1)*FSTAR_TAU(3,1) +
     +           FSTAR_TAU(3,1)*TSTAR_TAU(1,2)*FSTAR_TAU(3,2) +
     +           FSTAR_TAU(3,1)*TSTAR_TAU(1,3)*FSTAR_TAU(3,3) +
     +           FSTAR_TAU(3,2)*TSTAR_TAU(2,1)*FSTAR_TAU(3,1) +
     +           FSTAR_TAU(3,2)*TSTAR_TAU(2,2)*FSTAR_TAU(3,2) +
     +           FSTAR_TAU(3,2)*TSTAR_TAU(2,3)*FSTAR_TAU(3,3) +
     +           FSTAR_TAU(3,3)*TSTAR_TAU(3,1)*FSTAR_TAU(3,1) +
     +           FSTAR_TAU(3,3)*TSTAR_TAU(3,2)*FSTAR_TAU(3,2) +
     +           FSTAR_TAU(3,3)*TSTAR_TAU(3,3)*FSTAR_TAU(3,3))*
     +       DET_FSTAR_TAU_INV
     
          T_TAU(1,2) = 
     +          (FSTAR_TAU(1,1)*TSTAR_TAU(1,1)*FSTAR_TAU(2,1) +
     +           FSTAR_TAU(1,1)*TSTAR_TAU(1,2)*FSTAR_TAU(2,2) +
     +           FSTAR_TAU(1,1)*TSTAR_TAU(1,3)*FSTAR_TAU(2,3) +
     +           FSTAR_TAU(1,2)*TSTAR_TAU(2,1)*FSTAR_TAU(2,1) +
     +           FSTAR_TAU(1,2)*TSTAR_TAU(2,2)*FSTAR_TAU(2,2) +
     +           FSTAR_TAU(1,2)*TSTAR_TAU(2,3)*FSTAR_TAU(2,3) +
     +           FSTAR_TAU(1,3)*TSTAR_TAU(3,1)*FSTAR_TAU(2,1) +
     +           FSTAR_TAU(1,3)*TSTAR_TAU(3,2)*FSTAR_TAU(2,2) +
     +           FSTAR_TAU(1,3)*TSTAR_TAU(3,3)*FSTAR_TAU(2,3))*
     +       DET_FSTAR_TAU_INV

          T_TAU(1,3) = 
     +          (FSTAR_TAU(1,1)*TSTAR_TAU(1,1)*FSTAR_TAU(3,1) +
     +           FSTAR_TAU(1,1)*TSTAR_TAU(1,2)*FSTAR_TAU(3,2) +
     +           FSTAR_TAU(1,1)*TSTAR_TAU(1,3)*FSTAR_TAU(3,3) +
     +           FSTAR_TAU(1,2)*TSTAR_TAU(2,1)*FSTAR_TAU(3,1) +
     +           FSTAR_TAU(1,2)*TSTAR_TAU(2,2)*FSTAR_TAU(3,2) +
     +           FSTAR_TAU(1,2)*TSTAR_TAU(2,3)*FSTAR_TAU(3,3) +
     +           FSTAR_TAU(1,3)*TSTAR_TAU(3,1)*FSTAR_TAU(3,1) +
     +           FSTAR_TAU(1,3)*TSTAR_TAU(3,2)*FSTAR_TAU(3,2) +
     +           FSTAR_TAU(1,3)*TSTAR_TAU(3,3)*FSTAR_TAU(3,3))*
     +       DET_FSTAR_TAU_INV

          T_TAU(2,3) = 
     +          (FSTAR_TAU(2,1)*TSTAR_TAU(1,1)*FSTAR_TAU(3,1) +
     +           FSTAR_TAU(2,1)*TSTAR_TAU(1,2)*FSTAR_TAU(3,2) +
     +           FSTAR_TAU(2,1)*TSTAR_TAU(1,3)*FSTAR_TAU(3,3) +
     +           FSTAR_TAU(2,2)*TSTAR_TAU(2,1)*FSTAR_TAU(3,1) +
     +           FSTAR_TAU(2,2)*TSTAR_TAU(2,2)*FSTAR_TAU(3,2) +
     +           FSTAR_TAU(2,2)*TSTAR_TAU(2,3)*FSTAR_TAU(3,3) +
     +           FSTAR_TAU(2,3)*TSTAR_TAU(3,1)*FSTAR_TAU(3,1) +
     +           FSTAR_TAU(2,3)*TSTAR_TAU(3,2)*FSTAR_TAU(3,2) +
     +           FSTAR_TAU(2,3)*TSTAR_TAU(3,3)*FSTAR_TAU(3,3))*
     +       DET_FSTAR_TAU_INV
 
           T_TAU(2,1) = T_TAU(1,2)
           T_TAU(3,1) = T_TAU(1,3)
           T_TAU(3,2) = T_TAU(2,3)

C        IF(TOTAL_TIME.GE.10) THEN
C         Interval = 20.0
C         IF(STEP_TIME.LE.30) THEN
C            Interval = 1.0
C         END IF
C         IF(MOD(STEP_TIME,Interval).LE.DTIME) THEN
C          WRITE(100,152) ESTAR_TAU(1,1), ESTAR_TAU(2,2), ESTAR_TAU(3,3),
C     +    ESTAR_TAU(1,2), ESTAR_TAU(1,3), ESTAR_TAU(2,3), T_TAU(2,2),
C     +    STATE_OLD(KM,1)
C152       FORMAT(' ',8E12.4)
C         END IF
C        END IF
C
C          Calculate accumulated slip
C

          DO 315 ISLIP = 1,NSLIP
            ACCGAM_TAU(ISLIP) = ACCGAM_T(ISLIP) +
     +                            DGAMMA_T(ISLIP)
315         CONTINUE
           
C
C              Calculate s^{\alpha}_{\tau}
C

          IF (AH0 .EQ. ZERO) THEN 
C               Non-hardening case 
            DO 84 ISLIP = 1,NSLIP
              S_TAU(ISLIP) = S_T(ISLIP)
84            CONTINUE
          ELSE
C               Hardening case
C             write(*,*) 'DET_FP',DTIME
             DO 85 ISLIP = 1,NSLIP      
             IF(DABS(DGAMMA_T(ISLIP))/DTIME .LE. 1.E-10) THEN 
                 AH_T(ISLIP) = 0.
           ELSE
                 SA_T = S_T(ISLIP) 
             SFRAC = ONE - SA_T/SSAT

                IF (SFRAC .GE. ZERO) THEN
                    AH_T(ISLIP) = AH0*(DABS(SFRAC)**RHARD)
                ELSE
                    AH_T(ISLIP) = -AH0*(DABS(SFRAC)**RHARD)
             ENDIF              
      
           ENDIF             
85            CONTINUE
           
C
C                Calculate the hardening matrix h^{\alpha \beta}
C                at the beginning of the step
C
             DO 90 ISLIP = 1,NSLIP      
             DO 90 JSLIP = 1,NSLIP
              AHMAT_T(ISLIP,JSLIP) = QLAT(ISLIP,JSLIP)*AH_T(JSLIP)
90            CONTINUE

            DO 110 ISLIP = 1,NSLIP
              S_TAU(ISLIP) = S_T(ISLIP)
            DO 110 JSLIP = 1,NSLIP
              S_TAU(ISLIP) = S_TAU(ISLIP) + 
     +                               AHMAT_T(ISLIP,JSLIP)*
     +                         DABS(DGAMMA_T(JSLIP))
110            CONTINUE
          ENDIF
C   Read Back Stress
        Bij(1,1) = STATE_OLD(KM,10)
        Bij(2,2) = STATE_OLD(KM,11)
        Bij(3,3) = STATE_OLD(KM,12)
        Bij(1,2) = STATE_OLD(KM,13)
        Bij(1,3) = STATE_OLD(KM,14)
        Bij(2,3) = STATE_OLD(KM,15)
        Bij(2,1) = Bij(1,2)
        Bij(3,1) = Bij(1,3)
        Bij(3,2) = Bij(2,3)

C      Calculate the resolved shear stresses tau^\alpha(\tau)
C     at the end of the step.
C     Equation 11 of Kalidindi
C     Plus backstress term, modified by Yin Zhang, Bij = B*delta_ij
C       (Tij - Bij)*Sij 
          DO 120 ISLIP = 1,NSLIP
            RSS_TAU(ISLIP) = 
     +            TSTAR_TAU(1,1)*SMAT_0_ALL(NOEL_PT+KM,ICRYS,ISLIP,1,1) +
     +            TSTAR_TAU(1,2)*SMAT_0_ALL(NOEL_PT+KM,ICRYS,ISLIP,1,2) +
     +            TSTAR_TAU(1,3)*SMAT_0_ALL(NOEL_PT+KM,ICRYS,ISLIP,1,3) +
     +            TSTAR_TAU(2,1)*SMAT_0_ALL(NOEL_PT+KM,ICRYS,ISLIP,2,1) +
     +            TSTAR_TAU(2,2)*SMAT_0_ALL(NOEL_PT+KM,ICRYS,ISLIP,2,2) +
     +            TSTAR_TAU(2,3)*SMAT_0_ALL(NOEL_PT+KM,ICRYS,ISLIP,2,3) +
     +            TSTAR_TAU(3,1)*SMAT_0_ALL(NOEL_PT+KM,ICRYS,ISLIP,3,1) +
     +            TSTAR_TAU(3,2)*SMAT_0_ALL(NOEL_PT+KM,ICRYS,ISLIP,3,2) +
     +            TSTAR_TAU(3,3)*SMAT_0_ALL(NOEL_PT+KM,ICRYS,ISLIP,3,3)

          RSS_TAU_EFF(ISLIP) = RSS_TAU(ISLIP) -
     +       Bij(1,1)*SMAT_0_ALL(NOEL_PT+KM,ICRYS,ISLIP,1,1) -
     +       Bij(1,2)*SMAT_0_ALL(NOEL_PT+KM,ICRYS,ISLIP,1,2) -
     +       Bij(1,3)*SMAT_0_ALL(NOEL_PT+KM,ICRYS,ISLIP,1,3) -
     +       Bij(2,1)*SMAT_0_ALL(NOEL_PT+KM,ICRYS,ISLIP,2,1) -
     +       Bij(2,2)*SMAT_0_ALL(NOEL_PT+KM,ICRYS,ISLIP,2,2) -
     +       Bij(2,3)*SMAT_0_ALL(NOEL_PT+KM,ICRYS,ISLIP,2,3) -
     +       Bij(3,1)*SMAT_0_ALL(NOEL_PT+KM,ICRYS,ISLIP,3,1) -
     +       Bij(3,2)*SMAT_0_ALL(NOEL_PT+KM,ICRYS,ISLIP,3,2) -
     +       Bij(3,3)*SMAT_0_ALL(NOEL_PT+KM,ICRYS,ISLIP,3,3)




120          CONTINUE
C
C      Calculate the shear increments at the end
C     of the step \delta\gamma^\alpha(\tau). These will be stored
C     for use in the next increment
C

          DO 130 ISLIP = 1,NSLIP
              S_ALPHA = S_TAU(ISLIP) 
              TAU_ALPHA = DABS(RSS_TAU_EFF(ISLIP))

              IF ((TAU_ALPHA .EQ. ZERO) 
     +            .OR. (DTIME .EQ. ONE)) THEN
                DGAMMA_TAU(ISLIP) = ZERO
              GO TO 130
              ENDIF

             
             IF((DLOG10(DABS(TAU_ALPHA/S_ALPHA))/AM) .LT. 90.D0) THEN
                IF (RSS_TAU(ISLIP) .GE. ZERO) THEN
                    DGAMMA_TAU(ISLIP) =  DTIME*GDOT0*
     +                                  DABS(TAU_ALPHA/S_ALPHA)**(ONE/AM)              
                ELSE
                  DGAMMA_TAU(ISLIP) = -DTIME*GDOT0*
     +                                  DABS(TAU_ALPHA/S_ALPHA)**(ONE/AM)              
                ENDIF
              ELSE
                 WRITE(*,*)'TIME STEP TOO LARGE!! STOPPING!!'
                 STOP
              END IF
130          CONTINUE
C 
C                 Calculate the volume averaged Cauchy stress
C
          TBAR_TAU(1,1) = TBAR_TAU(1,1) +
     +                      T_TAU(1,1)*VOL_FRAC_ALL(NOEL_PT+KM,ICRYS) 
          TBAR_TAU(1,2) = TBAR_TAU(1,2) + 
     +                      T_TAU(1,2)*VOL_FRAC_ALL(NOEL_PT+KM,ICRYS) 
          TBAR_TAU(1,3) = TBAR_TAU(1,3) + 
     +                      T_TAU(1,3)*VOL_FRAC_ALL(NOEL_PT+KM,ICRYS) 
          TBAR_TAU(2,1) = TBAR_TAU(2,1) + 
     +                      T_TAU(2,1)*VOL_FRAC_ALL(NOEL_PT+KM,ICRYS) 
          TBAR_TAU(2,2) = TBAR_TAU(2,2) + 
     +                      T_TAU(2,2)*VOL_FRAC_ALL(NOEL_PT+KM,ICRYS) 
          TBAR_TAU(2,3) = TBAR_TAU(2,3) + 
     +                      T_TAU(2,3)*VOL_FRAC_ALL(NOEL_PT+KM,ICRYS) 
          TBAR_TAU(3,1) = TBAR_TAU(3,1) + 
     +                      T_TAU(3,1)*VOL_FRAC_ALL(NOEL_PT+KM,ICRYS) 
          TBAR_TAU(3,2) = TBAR_TAU(3,2) + 
     +                      T_TAU(3,2)*VOL_FRAC_ALL(NOEL_PT+KM,ICRYS) 
          TBAR_TAU(3,3) = TBAR_TAU(3,3) +
     +                      T_TAU(3,3)*VOL_FRAC_ALL(NOEL_PT+KM,ICRYS) 
C
C                Store  plastic deformation gradient
C       

          FP_T_INV_ALL(NOEL_PT+KM,ICRYS,1,1) = FP_TAU_INV(1,1)
          FP_T_INV_ALL(NOEL_PT+KM,ICRYS,1,2) = FP_TAU_INV(1,2)
          FP_T_INV_ALL(NOEL_PT+KM,ICRYS,1,3) = FP_TAU_INV(1,3)
          FP_T_INV_ALL(NOEL_PT+KM,ICRYS,2,1) = FP_TAU_INV(2,1)
          FP_T_INV_ALL(NOEL_PT+KM,ICRYS,2,2) = FP_TAU_INV(2,2)
          FP_T_INV_ALL(NOEL_PT+KM,ICRYS,2,3) = FP_TAU_INV(2,3)
          FP_T_INV_ALL(NOEL_PT+KM,ICRYS,3,1) = FP_TAU_INV(3,1)
          FP_T_INV_ALL(NOEL_PT+KM,ICRYS,3,2) = FP_TAU_INV(3,2)
          FP_T_INV_ALL(NOEL_PT+KM,ICRYS,3,3) = FP_TAU_INV(3,3)
C
C              Store the slip resistances, shear increments
C              and the accumulated sshears backstress, for all the slip systems
C

          DO 140 ISLIP = 1,NSLIP
            S_T_ALL(NOEL_PT+KM,ICRYS,ISLIP) = S_TAU(ISLIP)
            DGAMMA_T_ALL(NOEL_PT+KM,ICRYS,ISLIP) = DGAMMA_TAU(ISLIP)
          ACCGAM_T_ALL(NOEL_PT+KM,ICRYS,ISLIP) = ACCGAM_TAU(ISLIP)

          Bij(1,1) = Bij(1,1) +
     +    BH*C44*EXP(-3000*DABS(ACCGAM_TAU(ISLIP)))*DGAMMA_TAU(ISLIP)*
     +    SMAT_0_ALL(NOEL_PT+KM,ICRYS,ISLIP,1,1) 
          Bij(2,2) = Bij(2,2) +
     +    BH*C44*EXP(-3000*DABS(ACCGAM_TAU(ISLIP)))*DGAMMA_TAU(ISLIP)*
     +    SMAT_0_ALL(NOEL_PT+KM,ICRYS,ISLIP,2,2) 
          Bij(3,3) = Bij(3,3) +
     +    BH*C44*EXP(-3000*DABS(ACCGAM_TAU(ISLIP)))*DGAMMA_TAU(ISLIP)*
     +    SMAT_0_ALL(NOEL_PT+KM,ICRYS,ISLIP,3,3) 
          Bij(1,2) = Bij(1,2) +
     +    BH*C44*EXP(-3000*DABS(ACCGAM_TAU(ISLIP)))*DGAMMA_TAU(ISLIP)*
     +    (SMAT_0_ALL(NOEL_PT+KM,ICRYS,ISLIP,1,2) +
     +     SMAT_0_ALL(NOEL_PT+KM,ICRYS,ISLIP,2,1))/2
          Bij(1,3) = Bij(1,3) +
     +    BH*C44*EXP(-3000*DABS(ACCGAM_TAU(ISLIP)))*DGAMMA_TAU(ISLIP)*
     +    (SMAT_0_ALL(NOEL_PT+KM,ICRYS,ISLIP,1,3) + 
     +     SMAT_0_ALL(NOEL_PT+KM,ICRYS,ISLIP,3,1))/2
          Bij(2,3) = Bij(2,3) +
     +    BH*C44*EXP(-3000*DABS(ACCGAM_TAU(ISLIP)))*DGAMMA_TAU(ISLIP)*
     +    (SMAT_0_ALL(NOEL_PT+KM,ICRYS,ISLIP,2,3) +
     +     SMAT_0_ALL(NOEL_PT+KM,ICRYS,ISLIP,3,2))/2
140          CONTINUE

C
C           Calculate the plastic work increment
C
          DO 150 ISLIP = 1,NSLIP
            PLASTIC_WORK_INC = PLASTIC_WORK_INC + 
     +                           DABS(RSS_TAU(ISLIP)*
     +                           DGAMMA_TAU(ISLIP))*
     +                           VOL_FRAC_ALL(NOEL_PT+KM,ICRYS)
150          CONTINUE
C
C             Calculate the texture from the elastic rotations
C        
          IF ((TOTAL_TIME .GE. TEXCALTIME) .AND. (DTIME .NE. ONE)) THEN
          CALL POLAR(FSTAR_TAU,RSTAR_TAU,USTAR_TAU)
            
          QTEX(1,1) = RSTAR_TAU(1,1)*Q(1,1) +
     +                        RSTAR_TAU(1,2)*Q(2,1) +
     +                  RSTAR_TAU(1,3)*Q(3,1)
     
          QTEX(1,2) = RSTAR_TAU(1,1)*Q(1,2) + 
     +                      RSTAR_TAU(1,2)*Q(2,2) +
     +                  RSTAR_TAU(1,3)*Q(3,2)
     
          QTEX(1,3) = RSTAR_TAU(1,1)*Q(1,3) + 
     +                      RSTAR_TAU(1,2)*Q(2,3) +
     +                  RSTAR_TAU(1,3)*Q(3,3)
     
          QTEX(2,1) = RSTAR_TAU(2,1)*Q(1,1) + 
     +                      RSTAR_TAU(2,2)*Q(2,1) +
     +                  RSTAR_TAU(2,3)*Q(3,1)
     
          QTEX(2,2) = RSTAR_TAU(2,1)*Q(1,2) + 
     +                      RSTAR_TAU(2,2)*Q(2,2) +
     +                  RSTAR_TAU(2,3)*Q(3,2)
     
          QTEX(2,3) = RSTAR_TAU(2,1)*Q(1,3) + 
     +                      RSTAR_TAU(2,2)*Q(2,3) +
     +                  RSTAR_TAU(2,3)*Q(3,3)
     
          QTEX(3,1) = RSTAR_TAU(3,1)*Q(1,1) +
     +                      RSTAR_TAU(3,2)*Q(2,1) +
     +                  RSTAR_TAU(3,3)*Q(3,1)
     
          QTEX(3,2) = RSTAR_TAU(3,1)*Q(1,2) + 
     +                      RSTAR_TAU(3,2)*Q(2,2) +
     +                  RSTAR_TAU(3,3)*Q(3,2)
     
          QTEX(3,3) = RSTAR_TAU(3,1)*Q(1,3) + 
     +                      RSTAR_TAU(3,2)*Q(2,3) +
     +                  RSTAR_TAU(3,3)*Q(3,3)
C
C             Calculate the Euler angles TH, PH, OM for the
C             reoriented crystal
C          
          CALL EULANG(QTEX,TH,PH,OM,ICRYS,IEULERERR)
          IF (IEULERERR .EQ. 1) THEN
            WRITE(90,*) 'ERROR: MAIN:'
            WRITE(90,*) 'ERROR: Euler angles not calculated!!'
            WRITE(90,*) 'INOEL = ',NOEL_PT+KM,' ICRYS = ',ICRYS
             WRITE(90,*)'FSTAR_TAU:-'
              CALL PRTMAT(FSTAR_TAU,3,3)
             WRITE(90,*)'RSTAR_TAU:-'
              CALL PRTMAT(RSTAR_TAU,3,3)
          ELSE
            TH = TH*180./PI
            PH = PH*180./PI
            OM = OM*180./PI
            THK = TH
            PHIK = PH - 90.
            PSIK = 90. - OM
             WRITE(90,'(I5,2X,E11.4,4(2X,F10.4))')
     +            NOEL_PT+KM, TOTAL_TIME, 
     +                TH, PH, OM,
C                 WRITE(90,'(I5,2X,E11.4,4(2X,F10.4))')
C    +              NOEL_PT+KM, TOTAL_TIME, 
C    +              PSIK, THK, PHIK, 
     +              VOL_FRAC_ALL(NOEL_PT+KM,ICRYS)*VOLUME
          ENDIF
          
c            DO 160 ISLIP = 1,NSLIP
c            WRITE(100,'(I5,2X,E11.4,2X,I5,3(2X,E11.4))')
c     +              NOEL_PT+KM, TOTAL_TIME, ISLIP,
c     +              ACCGAM_TAU(ISLIP), DABS(RSS_TAU(ISLIP))/1.E06,
c     +              DABS(S_TAU(ISLIP))/1.E06
c160         CONTINUE          
        ENDIF

500        CONTINUE
C****************************************************************
C            End of loop over all crystals at an integration point
C*****************************************************************

C
C            Calculate the mean normal pressure -(1/3)tr{\bar T}
C

        PBAR = TBAR_TAU(1,1) + TBAR_TAU(2,2) + TBAR_TAU(3,3)
        PBAR = -ONE_THIRD*PBAR
C
C             Calculate the deviatoric part of {\bar T}
C
        TBARDEV_TAU(1,1) = TBAR_TAU(1,1) + PBAR*AIDENT(1,1)
        TBARDEV_TAU(1,2) = TBAR_TAU(1,2) + PBAR*AIDENT(1,2)
        TBARDEV_TAU(1,3) = TBAR_TAU(1,3) + PBAR*AIDENT(1,3)
        TBARDEV_TAU(2,1) = TBAR_TAU(2,1) + PBAR*AIDENT(2,1)
        TBARDEV_TAU(2,2) = TBAR_TAU(2,2) + PBAR*AIDENT(2,2)
        TBARDEV_TAU(2,3) = TBAR_TAU(2,3) + PBAR*AIDENT(2,3)
        TBARDEV_TAU(3,1) = TBAR_TAU(3,1) + PBAR*AIDENT(3,1)
        TBARDEV_TAU(3,2) = TBAR_TAU(3,2) + PBAR*AIDENT(3,2)
        TBARDEV_TAU(3,3) = TBAR_TAU(3,3) + PBAR*AIDENT(3,3)
C
C             Calculate the equivalent tensile stress
C

        EQSTRESS_TAU = TBARDEV_TAU(1,1)*TBARDEV_TAU(1,1) +
     +                   TBARDEV_TAU(1,2)*TBARDEV_TAU(1,2) +
     +                 TBARDEV_TAU(1,3)*TBARDEV_TAU(1,3) +
     +                 TBARDEV_TAU(2,1)*TBARDEV_TAU(2,1) +
     +                 TBARDEV_TAU(2,2)*TBARDEV_TAU(2,2) +
     +                 TBARDEV_TAU(2,3)*TBARDEV_TAU(2,3) +
     +                 TBARDEV_TAU(3,1)*TBARDEV_TAU(3,1) +
     +                 TBARDEV_TAU(3,2)*TBARDEV_TAU(3,2) +
     +                 TBARDEV_TAU(3,3)*TBARDEV_TAU(3,3)
        EQSTRESS_TAU = DSQRT(THREE_HALF*EQSTRESS_TAU)
C
C       Abaqus/Explicit uses the stress measure (R^T){\bar T}R
C       To calculate this, first  calculate the rotation 
C       R(\tau) from R = FU^{-1}
C                
        CALL MATINV(U_TAU,U_TAU_INV,DET_U_TAU_INV)

          R_TAU(1,1) = F_TAU(1,1)*U_TAU_INV(1,1) +
     +                 F_TAU(1,2)*U_TAU_INV(2,1) +
     +                 F_TAU(1,3)*U_TAU_INV(3,1)

          R_TAU(1,2) = F_TAU(1,1)*U_TAU_INV(1,2) +
     +                 F_TAU(1,2)*U_TAU_INV(2,2) +
     +                 F_TAU(1,3)*U_TAU_INV(3,2)

          R_TAU(1,3) = F_TAU(1,1)*U_TAU_INV(1,3) +
     +                 F_TAU(1,2)*U_TAU_INV(2,3) +
     +                 F_TAU(1,3)*U_TAU_INV(3,3)

          R_TAU(2,1) = F_TAU(2,1)*U_TAU_INV(1,1) +
     +                 F_TAU(2,2)*U_TAU_INV(2,1) +
     +                 F_TAU(2,3)*U_TAU_INV(3,1)

          R_TAU(2,2) = F_TAU(2,1)*U_TAU_INV(1,2) +
     +                 F_TAU(2,2)*U_TAU_INV(2,2) +
     +                 F_TAU(2,3)*U_TAU_INV(3,2)

          R_TAU(2,3) = F_TAU(2,1)*U_TAU_INV(1,3) +
     +                 F_TAU(2,2)*U_TAU_INV(2,3) +
     +                 F_TAU(2,3)*U_TAU_INV(3,3)

          R_TAU(3,1) = F_TAU(3,1)*U_TAU_INV(1,1) +
     +                 F_TAU(3,2)*U_TAU_INV(2,1) +
     +                 F_TAU(3,3)*U_TAU_INV(3,1) 
                     
          R_TAU(3,2) = F_TAU(3,1)*U_TAU_INV(1,2) +
     +                 F_TAU(3,2)*U_TAU_INV(2,2) +
     +                 F_TAU(3,3)*U_TAU_INV(3,2)

          R_TAU(3,3) = F_TAU(3,1)*U_TAU_INV(1,3) +
     +                 F_TAU(3,2)*U_TAU_INV(2,3) +
     +                 F_TAU(3,3)*U_TAU_INV(3,3)
C
C               Calculate the stress measure (R^T){\bar T}R
C
          STRESS(1,1) = R_TAU(1,1)*TBAR_TAU(1,1)*R_TAU(1,1) +
     +                  R_TAU(1,1)*TBAR_TAU(1,2)*R_TAU(2,1) +
     +                  R_TAU(1,1)*TBAR_TAU(1,3)*R_TAU(3,1) +
     +                  R_TAU(2,1)*TBAR_TAU(2,1)*R_TAU(1,1) +
     +                  R_TAU(2,1)*TBAR_TAU(2,2)*R_TAU(2,1) +
     +                  R_TAU(2,1)*TBAR_TAU(2,3)*R_TAU(3,1) +
     +                  R_TAU(3,1)*TBAR_TAU(3,1)*R_TAU(1,1) +
     +                  R_TAU(3,1)*TBAR_TAU(3,2)*R_TAU(2,1) +
     +                  R_TAU(3,1)*TBAR_TAU(3,3)*R_TAU(3,1)

          STRESS(2,2) = R_TAU(1,2)*TBAR_TAU(1,1)*R_TAU(1,2) +
     +                  R_TAU(1,2)*TBAR_TAU(1,2)*R_TAU(2,2) +
     +                  R_TAU(1,2)*TBAR_TAU(1,3)*R_TAU(3,2) +
     +                  R_TAU(2,2)*TBAR_TAU(2,1)*R_TAU(1,2) +
     +                  R_TAU(2,2)*TBAR_TAU(2,2)*R_TAU(2,2) +
     +                  R_TAU(2,2)*TBAR_TAU(2,3)*R_TAU(3,2) +
     +                  R_TAU(3,2)*TBAR_TAU(3,1)*R_TAU(1,2) +
     +                  R_TAU(3,2)*TBAR_TAU(3,2)*R_TAU(2,2) +
     +                  R_TAU(3,2)*TBAR_TAU(3,3)*R_TAU(3,2)

          STRESS(3,3) = R_TAU(1,3)*TBAR_TAU(1,1)*R_TAU(1,3) +
     +                  R_TAU(1,3)*TBAR_TAU(1,2)*R_TAU(2,3) +
     +                  R_TAU(1,3)*TBAR_TAU(1,3)*R_TAU(3,3) +
     +                  R_TAU(2,3)*TBAR_TAU(2,1)*R_TAU(1,3) +
     +                  R_TAU(2,3)*TBAR_TAU(2,2)*R_TAU(2,3) +
     +                  R_TAU(2,3)*TBAR_TAU(2,3)*R_TAU(3,3) +
     +                  R_TAU(3,3)*TBAR_TAU(3,1)*R_TAU(1,3) +
     +                  R_TAU(3,3)*TBAR_TAU(3,2)*R_TAU(2,3) +
     +                  R_TAU(3,3)*TBAR_TAU(3,3)*R_TAU(3,3)
     
          STRESS(1,2) = R_TAU(1,1)*TBAR_TAU(1,1)*R_TAU(1,2) +
     +                  R_TAU(1,1)*TBAR_TAU(1,2)*R_TAU(2,2) +
     +                  R_TAU(1,1)*TBAR_TAU(1,3)*R_TAU(3,2) +
     +                  R_TAU(2,1)*TBAR_TAU(2,1)*R_TAU(1,2) +
     +                  R_TAU(2,1)*TBAR_TAU(2,2)*R_TAU(2,2) +
     +                  R_TAU(2,1)*TBAR_TAU(2,3)*R_TAU(3,2) +
     +                  R_TAU(3,1)*TBAR_TAU(3,1)*R_TAU(1,2) +
     +                  R_TAU(3,1)*TBAR_TAU(3,2)*R_TAU(2,2) +
     +                  R_TAU(3,1)*TBAR_TAU(3,3)*R_TAU(3,2)

          STRESS(1,3) = R_TAU(1,1)*TBAR_TAU(1,1)*R_TAU(1,3) +
     +                  R_TAU(1,1)*TBAR_TAU(1,2)*R_TAU(2,3) +
     +                  R_TAU(1,1)*TBAR_TAU(1,3)*R_TAU(3,3) +
     +                  R_TAU(2,1)*TBAR_TAU(2,1)*R_TAU(1,3) +
     +                  R_TAU(2,1)*TBAR_TAU(2,2)*R_TAU(2,3) +
     +                  R_TAU(2,1)*TBAR_TAU(2,3)*R_TAU(3,3) +
     +                  R_TAU(3,1)*TBAR_TAU(3,1)*R_TAU(1,3) +
     +                  R_TAU(3,1)*TBAR_TAU(3,2)*R_TAU(2,3) +
     +                  R_TAU(3,1)*TBAR_TAU(3,3)*R_TAU(3,3)

          STRESS(2,3) = R_TAU(1,2)*TBAR_TAU(1,1)*R_TAU(1,3) +
     +                  R_TAU(1,2)*TBAR_TAU(1,2)*R_TAU(2,3) +
     +                  R_TAU(1,2)*TBAR_TAU(1,3)*R_TAU(3,3) +
     +                  R_TAU(2,2)*TBAR_TAU(2,1)*R_TAU(1,3) +
     +                  R_TAU(2,2)*TBAR_TAU(2,2)*R_TAU(2,3) +
     +                  R_TAU(2,2)*TBAR_TAU(2,3)*R_TAU(3,3) +
     +                  R_TAU(3,2)*TBAR_TAU(3,1)*R_TAU(1,3) +
     +                  R_TAU(3,2)*TBAR_TAU(3,2)*R_TAU(2,3) +
     +                  R_TAU(3,2)*TBAR_TAU(3,3)*R_TAU(3,3)
     
          STRESS(2,1) = STRESS(1,2)
          STRESS(3,1) = STRESS(1,3)
          STRESS(3,2) = STRESS(2,3)
      
          
C
C      Calculate the plastic work increment, plastic strain increment,
C      total plastic dissipation

        IF (PLASTIC_WORK_INC .EQ. ZERO) THEN
          DEQPS = ZERO
        ELSE
          DEQPS = PLASTIC_WORK_INC/EQSTRESS_TAU
        ENDIF
        ENER_INELAS_NEW(KM) = ENER_INELAS_OLD(KM) +
     +                        PLASTIC_WORK_INC/DENSITY(KM)

C
C         Update the stress and state variables
C

        DO 170 I = 1,NDIR
          STRESS_NEW(KM,I) = STRESS(I,I)
170        CONTINUE
        IF (NSHR .NE. 0) THEN
          STRESS_NEW(KM,NDIR+1) = STRESS(1,2)
          IF (NSHR .NE. 1) THEN
            STRESS_NEW(KM,NDIR+2) = STRESS(2,3)
            IF (NSHR .NE. 2) THEN
            STRESS_NEW(KM,NDIR+3) = STRESS(1,3)
            ENDIF 
          ENDIF 
        ENDIF 
        
        IF (DTIME .NE. ONE) THEN
          STATE_NEW(KM,1) = STATE_OLD(KM,1) + DEQPS
          STATE_NEW(KM,2) = DEQPS/DTIME
          STATE_NEW(KM,3) = EQSTRESS_TAU
          
            STATE_NEW(KM,4) = ESTAR_TAU(1,1)
            STATE_NEW(KM,5) = ESTAR_TAU(2,2)
            STATE_NEW(KM,6) = ESTAR_TAU(3,3)
            STATE_NEW(KM,7) = ESTAR_TAU(1,2)
            STATE_NEW(KM,8) = ESTAR_TAU(1,3)
            STATE_NEW(KM,9) = ESTAR_TAU(2,3)
            STATE_NEW(KM,10) = Bij(1,1)
            STATE_NEW(KM,11) = Bij(2,2)
            STATE_NEW(KM,12) = Bij(3,3)
            STATE_NEW(KM,13) = Bij(1,2) 
            STATE_NEW(KM,14) = Bij(1,3)
            STATE_NEW(KM,15) = Bij(2,3)
            STATE_NEW(KM,16) = RSS_TAU(1)
            STATE_NEW(KM,17) = RSS_TAU(2)
            STATE_NEW(KM,18) = RSS_TAU(3)

C	  DEBUG
C     Debug session
C           if(NOEL_PT+KM==1) then
C			 write(*,*) C11,C12,C44,GDOT0,AM,S0,AH0,SSAT,RHARD,
C     +       QL,TEXCALTIME,TEXFREQ
C 
C           endif
C           if(NOEL_PT+KM==6) then
C			 write(*,*) C11,C12,C44,GDOT0,AM,S0,AH0,SSAT,RHARD,
C     +       QL,TEXCALTIME,TEXFREQ
C 
C           endif

        ENDIF
C        
C          Update the specific internal energy
C
          STRESS_POWER = 0.
        DO 180 I = 1,NDIR
            STRESS_POWER = STRESS_POWER +
     +      ONE_HALF*((STRESS_OLD(KM,I)+STRESS_NEW(KM,I))*
     +            STRAIN_INC(KM,I))
180       CONTINUE
          IF (NSHR .NE. 0) THEN
            STRESS_POWER = STRESS_POWER +
     +      ONE_HALF*((STRESS_OLD(KM,NDIR+1)+STRESS_NEW(KM,NDIR+1))*
     +            STRAIN_INC(KM,NDIR+1))
            IF (NSHR .NE. 1) THEN
              STRESS_POWER = STRESS_POWER +
     +        ONE_HALF*((STRESS_OLD(KM,NDIR+2)+STRESS_NEW(KM,NDIR+2))*
     +              STRAIN_INC(KM,NDIR+2))
              IF (NSHR .NE. 2) THEN
                STRESS_POWER = STRESS_POWER +
     +          ONE_HALF*((STRESS_OLD(KM,NDIR+3)+STRESS_NEW(KM,NDIR+3))*
     +                STRAIN_INC(KM,NDIR+3))
              ENDIF
            ENDIF
          ENDIF

          ENER_INTERN_NEW(KM) = ENER_INTERN_OLD(KM) +
     +                          STRESS_POWER/DENSITY(KM)
      
1000      CONTINUE
C***********************************************************************
C          End of loop over elements in an NBLOCK. Update the element
C          number pointer.
C***********************************************************************
      NOEL_PT = NOEL_PT + NBLOCK
      
        
        IF ((TOTAL_TIME .GE. TEXCALTIME) .AND. (NOEL_PT .EQ. NOEL)) THEN
          TEXCALTIME = TEXCALTIME + TEXFREQ
        ENDIF
      
      RETURN
      END




C**************************************************************************
C    Utility Subroutines
C****************************************************************************
      SUBROUTINE MATINV(A,A_INV,DET_A_INV)

C      This subroutine calculates the inverse of a {3 x 3} matrix and the
C      determinant of the inverse
C****************************************************************************      

      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION A(3,3), A_INV(3,3)
      
      PARAMETER(ZERO=0., ONE=1.)

      DET_A = A(1,1)*(A(2,2)*A(3,3) - A(3,2)*A(2,3)) -
     +              A(2,1)*(A(1,2)*A(3,3) - A(3,2)*A(1,3)) +
     +              A(3,1)*(A(1,2)*A(2,3) - A(2,2)*A(1,3))

      IF (DET_A .LE. ZERO) THEN
        WRITE(80,*) 'WARNING: SUBROUTINE MATINV:'
        WRITE(80,*) 'WARNING: DET of MAT is zero/negative!!'
      ENDIF

      DET_A_INV = ONE/DET_A

      A_INV(1,1) = DET_A_INV*(A(2,2)*A(3,3)-A(3,2)*A(2,3))
      A_INV(1,2) = DET_A_INV*(A(3,2)*A(1,3)-A(1,2)*A(3,3))
      A_INV(1,3) = DET_A_INV*(A(1,2)*A(2,3)-A(2,2)*A(1,3))
      A_INV(2,1) = DET_A_INV*(A(3,1)*A(2,3)-A(2,1)*A(3,3))
      A_INV(2,2) = DET_A_INV*(A(1,1)*A(3,3)-A(3,1)*A(1,3))
      A_INV(2,3) = DET_A_INV*(A(2,1)*A(1,3)-A(1,1)*A(2,3))
      A_INV(3,1) = DET_A_INV*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
      A_INV(3,2) = DET_A_INV*(A(3,1)*A(1,2)-A(1,1)*A(3,2))
      A_INV(3,3) = DET_A_INV*(A(1,1)*A(2,2)-A(2,1)*A(1,2))

      RETURN
      END
      
C****************************************************************************
      SUBROUTINE ROTMAT(TH,PH,OM,Q)
C****************************************************************************

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Q(3,3)

      STH = DSIN(TH)
      CTH = DCOS(TH)

      SPH = DSIN(PH)
      CPH = DCOS(PH)

      SOM = DSIN(OM)
      COM = DCOS(OM)

      Q(1,1) = CPH*COM - SPH*SOM*CTH
      Q(1,2) = SPH*COM + CPH*SOM*CTH
      Q(1,3) = SOM*STH

      Q(2,1) = -CPH*SOM - SPH*COM*CTH
      Q(2,2) = -SPH*SOM + CPH*COM*CTH
      Q(2,3) =  COM*STH

      Q(3,1) = SPH*STH
      Q(3,2) = -CPH*STH
      Q(3,3) = CTH

      RETURN
      END

C****************************************************************************
      SUBROUTINE PRTMAT(A,M,N)

C      Print the matrix A of dimension {m x n}
C****************************************************************************

      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION A(M,N)

      DO 10 I = 1,M
        WRITE(80,*) (A(I,J),J=1,N)
10      CONTINUE

      RETURN
      END

C****************************************************************************
      SUBROUTINE PRTMAT4(A,I,J,K,L)

C      Print the matrix A of dimension {i x j x k x l}
C****************************************************************************

      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION A(I,J,K,L)

      DO 10 II = 1,I
      DO 10 JJ = 1,J
      DO 10 KK = 1,K
      DO 10 LL = 1,L
        WRITE(80,*)'A(',II,',',JJ,',',KK,',',LL,') = ',A(II,JJ,KK,LL)
10      CONTINUE

      RETURN
      END

C****************************************************************************
      SUBROUTINE PRTVEC(A,L)

C      Print the vector A of dimension {l}
C****************************************************************************

      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION A(L)

      WRITE(80,*) (A(I), I=1,L)

      RETURN
      END

C****************************************************************************
      SUBROUTINE POLAR(F,R,U)

C      This subroutine performs the right polar decomposition 
C      [F] = [R][U] of the deformation gradient [F] into a 
C       rotation [R] and the right stretch tensor [U]
C****************************************************************************

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION F(3,3), C(3,3), R(3,3), U(3,3), U_INV(3,3),
     +              U_TEMP(3,3), U_EIGVAL(3), OMEGA(3), EIGVEC(3,3)

      PARAMETER(ZERO=0.D0)

C      Check if the determinant of [F] is greater than zero.
C      If not, then print diagnostic and stop.

      DET_F = F(1,1)*(F(2,2)*F(3,3) - F(3,2)*F(2,3)) -
     +              F(2,1)*(F(1,2)*F(3,3) - F(3,2)*F(1,3)) +
     +              F(3,1)*(F(1,2)*F(2,3) - F(2,2)*F(1,3))

        IF (DET_F .LE. ZERO) THEN
          WRITE(80,*) 'ERROR: SUBROUTINE POLAR:'
          WRITE(80,*) 'ERROR: DETF is negative/zero.....stopping!!'
          WRITE(80,*) 'DET_F = ',DET_F
        STOP
        ENDIF

C      Calculate the right Cauchy Green strain tensor [C]

      C(1,1) = F(1,1)*F(1,1) + F(2,1)*F(2,1) + F(3,1)*F(3,1)
      C(1,2) = F(1,1)*F(1,2) + F(2,1)*F(2,2) + F(3,1)*F(3,2)
      C(1,3) = F(1,1)*F(1,3) + F(2,1)*F(2,3) + F(3,1)*F(3,3)
      C(2,1) = F(1,2)*F(1,1) + F(2,2)*F(2,1) + F(3,2)*F(3,1)
      C(2,2) = F(1,2)*F(1,2) + F(2,2)*F(2,2) + F(3,2)*F(3,2)
      C(2,3) = F(1,2)*F(1,3) + F(2,2)*F(2,3) + F(3,2)*F(3,3)
      C(3,1) = F(1,3)*F(1,1) + F(2,3)*F(2,1) + F(3,3)*F(3,1)
      C(3,2) = F(1,3)*F(1,2) + F(2,3)*F(2,2) + F(3,3)*F(3,2)
      C(3,3) = F(1,3)*F(1,3) + F(2,3)*F(2,3) + F(3,3)*F(3,3)

C      Calculate the eigenvalues and eigenvectors of [C]

      CALL SPECTRAL(C,OMEGA,EIGVEC)

C      Calculate the principal values of the stretch [U]

      U_EIGVAL(1) = DSQRT(OMEGA(1))
      U_EIGVAL(2) = DSQRT(OMEGA(2))
      U_EIGVAL(3) = DSQRT(OMEGA(3))

      U(1,1) = U_EIGVAL(1)
      U(1,2) = ZERO
      U(1,3) = ZERO
      U(2,1) = ZERO
      U(2,2) = U_EIGVAL(2)
      U(2,3) = ZERO
      U(3,1) = ZERO
      U(3,2) = ZERO
      U(3,3) = U_EIGVAL(3)

C      Calculate the complete tensor [U]

      U_TEMP(1,1) = EIGVEC(1,1)*U(1,1)*EIGVEC(1,1) +
     +                    EIGVEC(1,1)*U(1,2)*EIGVEC(1,2) +
     +                    EIGVEC(1,1)*U(1,3)*EIGVEC(1,3) +
     +                    EIGVEC(1,2)*U(2,1)*EIGVEC(1,1) +
     +                    EIGVEC(1,2)*U(2,2)*EIGVEC(1,2) +
     +                    EIGVEC(1,2)*U(2,3)*EIGVEC(1,3) +
     +                    EIGVEC(1,3)*U(3,1)*EIGVEC(1,1) +
     +                    EIGVEC(1,3)*U(3,2)*EIGVEC(1,2) +
     +                    EIGVEC(1,3)*U(3,3)*EIGVEC(1,3)

      U_TEMP(2,2) = EIGVEC(2,1)*U(1,1)*EIGVEC(2,1) +
     +                    EIGVEC(2,1)*U(1,2)*EIGVEC(2,2) +
     +                    EIGVEC(2,1)*U(1,3)*EIGVEC(2,3) +
     +                    EIGVEC(2,2)*U(2,1)*EIGVEC(2,1) +
     +                    EIGVEC(2,2)*U(2,2)*EIGVEC(2,2) +
     +                    EIGVEC(2,2)*U(2,3)*EIGVEC(2,3) +
     +                    EIGVEC(2,3)*U(3,1)*EIGVEC(2,1) +
     +                    EIGVEC(2,3)*U(3,2)*EIGVEC(2,2) +
     +                    EIGVEC(2,3)*U(3,3)*EIGVEC(2,3)

      U_TEMP(3,3) = EIGVEC(3,1)*U(1,1)*EIGVEC(3,1) +
     +                    EIGVEC(3,1)*U(1,2)*EIGVEC(3,2) +
     +                    EIGVEC(3,1)*U(1,3)*EIGVEC(3,3) +
     +                    EIGVEC(3,2)*U(2,1)*EIGVEC(3,1) +
     +                    EIGVEC(3,2)*U(2,2)*EIGVEC(3,2) +
     +                    EIGVEC(3,2)*U(2,3)*EIGVEC(3,3) +
     +                    EIGVEC(3,3)*U(3,1)*EIGVEC(3,1) +
     +                    EIGVEC(3,3)*U(3,2)*EIGVEC(3,2) +
     +                    EIGVEC(3,3)*U(3,3)*EIGVEC(3,3)

      U_TEMP(1,2) = EIGVEC(1,1)*U(1,1)*EIGVEC(2,1) +
     +                    EIGVEC(1,1)*U(1,2)*EIGVEC(2,2) +
     +                    EIGVEC(1,1)*U(1,3)*EIGVEC(2,3) +
     +                    EIGVEC(1,2)*U(2,1)*EIGVEC(2,1) +
     +                    EIGVEC(1,2)*U(2,2)*EIGVEC(2,2) +
     +                    EIGVEC(1,2)*U(2,3)*EIGVEC(2,3) +
     +                    EIGVEC(1,3)*U(3,1)*EIGVEC(2,1) +
     +                    EIGVEC(1,3)*U(3,2)*EIGVEC(2,2) +
     +                    EIGVEC(1,3)*U(3,3)*EIGVEC(2,3)

      U_TEMP(2,3) = EIGVEC(2,1)*U(1,1)*EIGVEC(3,1) +
     +                    EIGVEC(2,1)*U(1,2)*EIGVEC(3,2) +
     +                    EIGVEC(2,1)*U(1,3)*EIGVEC(3,3) +
     +                    EIGVEC(2,2)*U(2,1)*EIGVEC(3,1) +
     +                    EIGVEC(2,2)*U(2,2)*EIGVEC(3,2) +
     +                    EIGVEC(2,2)*U(2,3)*EIGVEC(3,3) +
     +                    EIGVEC(2,3)*U(3,1)*EIGVEC(3,1) +
     +                    EIGVEC(2,3)*U(3,2)*EIGVEC(3,2) +
     +                    EIGVEC(2,3)*U(3,3)*EIGVEC(3,3)

      U_TEMP(1,3) = EIGVEC(1,1)*U(1,1)*EIGVEC(3,1) +
     +                    EIGVEC(1,1)*U(1,2)*EIGVEC(3,2) +
     +                    EIGVEC(1,1)*U(1,3)*EIGVEC(3,3) +
     +                    EIGVEC(1,2)*U(2,1)*EIGVEC(3,1) +
     +                    EIGVEC(1,2)*U(2,2)*EIGVEC(3,2) +
     +                    EIGVEC(1,2)*U(2,3)*EIGVEC(3,3) +
     +                    EIGVEC(1,3)*U(3,1)*EIGVEC(3,1) +
     +                    EIGVEC(1,3)*U(3,2)*EIGVEC(3,2) +
     +                    EIGVEC(1,3)*U(3,3)*EIGVEC(3,3)

      U_TEMP(2,1) = U_TEMP(1,2)
      U_TEMP(3,1) = U_TEMP(1,3)
      U_TEMP(3,2) = U_TEMP(2,3)

      U(1,1) = U_TEMP(1,1)
      U(1,2) = U_TEMP(1,2)
      U(1,3) = U_TEMP(1,3)
      U(2,1) = U_TEMP(2,1)
      U(2,2) = U_TEMP(2,2)
      U(2,3) = U_TEMP(2,3)
      U(3,1) = U_TEMP(3,1)
      U(3,2) = U_TEMP(3,2)
      U(3,3) = U_TEMP(3,3)

C      Calculate the inverse of [U]

      CALL MATINV(U,U_INV,DET_U_NIV)

C      Calculate the rotation [R]

      R(1,1) = F(1,1)*U_INV(1,1) + F(1,2)*U_INV(2,1) +
     +             F(1,3)*U_INV(3,1)
      R(1,2) = F(1,1)*U_INV(1,2) + F(1,2)*U_INV(2,2) +
     +             F(1,3)*U_INV(3,2)
      R(1,3) = F(1,1)*U_INV(1,3) + F(1,2)*U_INV(2,3) +
     +             F(1,3)*U_INV(3,3)
      R(2,1) = F(2,1)*U_INV(1,1) + F(2,2)*U_INV(2,1) +
     +             F(2,3)*U_INV(3,1)
      R(2,2) = F(2,1)*U_INV(1,2) + F(2,2)*U_INV(2,2) +
     +             F(2,3)*U_INV(3,2)
      R(2,3) = F(2,1)*U_INV(1,3) + F(2,2)*U_INV(2,3) +
     +             F(2,3)*U_INV(3,3)
      R(3,1) = F(3,1)*U_INV(1,1) + F(3,2)*U_INV(2,1) +
     +             F(3,3)*U_INV(3,1)
      R(3,2) = F(3,1)*U_INV(1,2) + F(3,2)*U_INV(2,2) +
     +             F(3,3)*U_INV(3,2)
      R(3,3) = F(3,1)*U_INV(1,3) + F(3,2)*U_INV(2,3) +
     +             F(3,3)*U_INV(3,3)

      RETURN
      END

C**********************************************************************
      SUBROUTINE SPECTRAL(A,D,V)

C      This subroutine calculates the spectral decomposition
C      of a symmetric 3 x 3 matrix
C
C      The eigenvalues and eigenvextors of a symmetric 3 x 3 matrix 
C      [A] are calculated; the 3 eigenvalues are stored in
C      ascending order in the vector {D} and the eigenvectors in the
C      corresponding columns of the matrix [V]
C**********************************************************************

      IMPLICIT REAL*8 (A-H,O-Z)

      PARAMETER(NP=3)

      DIMENSION D(NP),V(NP,NP)
      DIMENSION A(3,3),E(NP,NP)


      DO 10 I = 1,3
          DO 20 J= 1,3
            E(I,J) = A(I,J)
20        CONTINUE
10      CONTINUE

      CALL JACOBI(E,3,NP,D,V,NROT)
      CALL EIGSRT(D,V,3,NP)

      RETURN
      END

C**********************************************************************
      SUBROUTINE JACOBI(A,N,NP,D,V,NROT)

C      This subroutine computes all eigenvalues and eigenvectors of
C      a real symmetric matrix [A], which is of size N x N, stored in
C      a physical NP x NP array. On output, the elements of [A] above
C      the diagonal are destroyed, but the diagonal and sub-diagonal
C      are unchanged and give full information about the original
C      symmetric matrix. Vector {D} resturns the eigenvalues of [A]
C      in its first N elements. [V] is a matrix with the same logical
C      and physical dimensions as [A] whose columns contain, an output,
C      the normalized eigenvectors of [A]. NROT returns the number of
C      Jacobi rotations which were required

C      Note: This subroutine is taken from "Numerical Recipes", 
C            Page 346.
C**********************************************************************

      IMPLICIT REAL*8 (A-H,O-Z)

      PARAMETER (NMAX =100, ZERO=0., ONE=1., TWO=2.)

      DIMENSION A(NP,NP),D(NP),V(NP,NP),B(NMAX),Z(NMAX)

C      Initialize [V] to the identity matrix

      DO 12 IP = 1,N      
        DO 11 IQ = 1,N
          V(IP,IQ) = ZERO
11        CONTINUE
          V(IP,IP) = ONE
12      CONTINUE

C      Inititalize {B} and {D} to the diagonal of [A], and {Z} to
C      zero. The vector {Z} will accumulate terms of the form T*A_pq
C      as in equation (11.1.14)

      DO 13 IP = 1,N
        B(IP) = A(IP,IP)
        D(IP) = B(IP)
        Z(IP) = ZERO
13      CONTINUE
C
      NROT = 0
      DO 24 I = 1,50

C      Sum off-diagonal elements

          SM = 0.D0
          DO 15 IP = 1, N-1
            DO 14 IQ = IP + 1, N
            SM = SM + DABS ( A(IP,IQ ))
14          CONTINUE
15        CONTINUE

C      If SUM = ZERO, then return. This is the normal return
C      which relies on quadratic convergence to machine underflow

          IF ( SM .EQ. ZERO) RETURN

C      In the first three sweeps carry out the PQ rotation only if
C      |A_pq| > TRESH, where TRESH is some threshold value, 
C      see equation (11.1.25); thereafter TRESH = ZERO

          IF ( I .LT. 4) THEN
            TRESH = 0.2D0*SM/N**2
          ELSE
            TRESH = ZERO
          ENDIF
C
          DO 22 IP = 1, N-1
            DO 21 IQ = IP+1,N
              G = 100.D0*DABS(A(IP,IQ))

C      After four sweeps, skip the rotation if the off-diagonal
C      element is small.

            IF ((I .GT. 4) .AND. (DABS(D(IP))+G .EQ. DABS(D(IP)))
     +            .AND. ( DABS(D(IQ))+G .EQ. DABS(D(IQ)))) THEN
                A(IP,IQ) = ZERO
              ELSE IF ( DABS(A(IP,IQ)) .GT. TRESH) THEN
                H = D(IQ) - D(IP)
                IF (DABS(H)+G .EQ. DABS(H)) THEN

C      T = 1./(2.*THETA), equation (11.1.10)

                T =A(IP,IQ)/H
              ELSE
                THETA = 0.5D0*H/A(IP,IQ)
                T = ONE/(DABS(THETA)+DSQRT(ONE+THETA**2))
                IF (THETA .LT. ZERO) T = -T
              ENDIF
              C = ONE/DSQRT(ONE + T**2)
              S = T*C
              TAU = S/(ONE + C)
              H = T*A(IP,IQ)
              Z(IP) = Z(IP) - H
              Z(IQ) = Z(IQ) + H
              D(IP) = D(IP) - H
              D(IQ) = D(IQ) + H
              A(IP,IQ) = ZERO

C      Case of rotations 1 <= J < P
                        
              DO 16 J = 1, IP-1
                G = A(J,IP)
                H = A(J,IQ)
                A(J,IP) = G - S*(H + G*TAU)
                A(J,IQ) = H + S*(G - H*TAU)
16              CONTINUE

C      Case of rotations P < J < Q

              DO 17 J = IP+1, IQ-1
                G = A(IP,J)
                H = A(J,IQ)
                A(IP,J) = G - S*(H + G*TAU)
                A(J,IQ) = H + S*(G - H*TAU)
17              CONTINUE

C      Case of rotations Q < J <= N

              DO 18 J = IQ+1, N
                  G = A(IP,J)
                H = A(IQ,J)
                A(IP,J) = G - S*(H + G*TAU)
                A(IQ,J) = H + S*(G - H*TAU)
18              CONTINUE
              DO 19 J = 1,N
                G = V(J,IP)
                H = V(J,IQ)
                V(J,IP) = G - S*(H + G*TAU)
                V(J,IQ) = H + S*(G - H*TAU)
19              CONTINUE
              NROT = NROT + 1
              ENDIF
21          CONTINUE
22        CONTINUE

C      Update {D} with the sum of T*A_pq, and reinitialize {Z}

        DO 23 IP = 1, N
          B(IP) = B(IP) + Z(IP)
          D(IP) = B(IP)
          Z(IP) = ZERO
23        CONTINUE
24      CONTINUE

C      If the algorithm has reached this stage, then there are
C      too many sweeps, print diagnostic and stop.

      WRITE (80,*) 'ERROR: SUBROUTINE JACOBI:'
      WRITE (80,*) 'ERROR: 50 iterations in Jacobi..not converged!!'

      RETURN
      END

C**********************************************************************
      SUBROUTINE EIGSRT(D,V,N,NP)

C      Given the eigenvalues {D} and eigenvectors [V] as output from
C      the subroutine JACOBI, this subroutine sorts the eigenvalues
C      into its ascending order, and rearranges the columns of [V]
C      accordingly.

C      Note: This subroutine is taken from "Numerical Recipes", 
C      Page 348.
C**********************************************************************

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION D(NP),V(NP,NP)

      DO 10 I = 1,N-1
        K = I
        P = D(I)
        DO 20 J = I+1,N
          IF (D(J) .GE. P) THEN
            K = J
            P = D(J)
          END IF
20        CONTINUE
        IF (K .NE. I) THEN
          D(K) = D(I)
          D(I) = P
          DO 30 J = 1,N
            P = V(J,I)
            V(J,I) = V(J,K)
            V(J,K) = P
30          CONTINUE
          ENDIF
10      CONTINUE

      RETURN
      END

C**********************************************************************
      SUBROUTINE EULANG(Q,TH,PH,OM,ICRYS,IEULERERR)
C*******************************************************************

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Q(3,3)

      PARAMETER(ZERO=0., ONE=1., TWO=2.)

        PI = 4.D0*DATAN(ONE)
      ICHK = 0
        IEULERERR = 0
      IF (DABS(Q(3,3))-ONE .GT. 1.0D-6) THEN
          WRITE(90,*) 'Q',((Q(I,J),J=1,3),I=1,3)
        WRITE(90,*) 'Q(3,3) > ONE'
        IEULERERR = 1
        RETURN
      ENDIF
      DO 10 I = 1,3
      DO 10 J = 1,3
        IF (DABS(Q(I,J)) .LT. 1.0D-6) Q(I,J) = ZERO
10      CONTINUE
        IF (DABS(DABS(Q(3,3))-ONE) .LT. 1.0D-6) THEN
          CALL EULCHECK2(Q,TH,PH,OM,ICHK)
        IF (ICHK .NE. 1) GO TO 20
        RETURN
      ENDIF
      TH = DACOS(Q(3,3))
      STH = DSIN(TH)
      OM = DATAN2(Q(1,3)/STH,Q(2,3)/STH)
      PH = DATAN2(Q(3,1)/STH,-Q(3,2)/STH)
      CALL EULCHECK1(Q,TH,PH,OM,ICHK)
      IF (ICHK .EQ. 1) RETURN
      TH = TWO*PI - TH
      STH = DSIN(TH)
      OM = DATAN2(Q(1,3)/STH,Q(2,3)/STH)
      PH = DATAN2(Q(3,1)/STH,-Q(3,2)/STH)
      CALL EULCHECK1(Q,TH,PH,OM,ICHK)
20      IF (ICHK .NE. 1) THEN
          WRITE(90,*) 'ERROR: SUBROUTINE EULANG:'
          WRITE(90,*) 'ERROR: Failed to find Euler angles!!'
        WRITE(90,*) 'ICRYS = ',ICRYS
        WRITE(90,*) 'Q',((Q(J,K),K=1,3),J=1,3)
        IEULERERR = 1
        RETURN
      ENDIF

      RETURN
      END

C*******************************************************************
        SUBROUTINE EULCHECK1(Q,TH,PH,OM,ICHK)
C*******************************************************************

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Q(3,3)

      TOL = 1.0D-3

      A =  DCOS(PH)*DCOS(OM) - DSIN(PH)*DSIN(OM)*DCOS(TH)
      B = -DSIN(OM)*DCOS(PH) - DCOS(OM)*DSIN(PH)*DCOS(TH)
      C =  DCOS(OM)*DSIN(PH) + DSIN(OM)*DCOS(PH)*DCOS(TH)
      D = -DSIN(PH)*DSIN(OM) + DCOS(PH)*DCOS(OM)*DCOS(TH)

        IF (DABS(A-Q(1,1)) .LT. TOL 
     +          .AND. DABS(B-Q(2,1)) .LT. TOL
     +      .AND. DABS(C-Q(1,2)) .LT. TOL 
     +      .AND. DABS(D-Q(2,2)) .LT. TOL) ICHK = 1

      RETURN
      END

C*******************************************************************
        SUBROUTINE EULCHECK2(Q,TH,PH,OM,ICHK)
C*******************************************************************

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Q(3,3)

      PARAMETER(ZERO=0., ONE=1., TWO=2.)

      TOL=1.0D-3
      Q(3,3) = ONE*Q(3,3)/DABS(Q(3,3))
      TH = DACOS(Q(3,3))
      IF (DABS(Q(1,3)) .GT. TOL) RETURN
      IF (DABS(Q(2,3)) .GT. TOL) RETURN
      IF (DABS(Q(3,1)) .GT. TOL) RETURN
      IF (DABS(Q(3,2)) .GT. TOL) RETURN
      IF (Q(3,3) .EQ. ONE) THEN
        IF (DABS(Q(1,1) - Q(2,2)) .GT. TOL) RETURN
        IF (DABS(Q(1,2) + Q(2,1)) .GT. TOL) RETURN
      ELSE
        IF (DABS(Q(1,1) + Q(2,2)) .GT. TOL) RETURN
        IF (DABS(Q(1,2) - Q(2,1)) .GT. TOL) RETURN
      ENDIF
      PH = DATAN2(Q(1,2),Q(1,1))
      OM = ZERO
      ICHK = 1

      RETURN
      END

C*******************************************************************
