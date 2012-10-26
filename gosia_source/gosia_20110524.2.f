c     Modified rachel version based on 20081208.10.  
c     This version uses PRT, option 9,I to select additional output:
c     I = 11  print excitation amplitudes & probabilities during the integration over omega
c     1<=I<=6 print the adiabatic exponential.  I is the multipolarity used to look up the adiab. exp.
c     I < 0   print the electric collision functions for lambda = abs(I)
c     All other output and calculations are identical to 20081208.10
c
C                            GOSIA 20110524 (64-bit)
C
C                            24 May 2011 Update
C
C       http://www.pas.rochester.edu/~cline/Gosia/index.html
C
C       Gosia was developed in 1980 by T. Czosnyka, D. Cline and C.Y. Wu at
C       the University of Rochester, Rochester, NY, USA.  
C
C       The coding of the Gosia suite of codes was maintained by Tomasz
C       Czosnyka from 1980 until his untimely death in 2006. Since 2007 Nigel
C       Warr has upgraded and maintained the coding of Gosia and Gosia2. 
C
C       Responsibility for the management and development of the Gosia suite
C       of codes was assumed by the Steering Committee of the Gosia Users Group
C       in April 2008.
C
C       The Gosia Steering Committee comprises:
C          Douglas Cline,(University of Rochester)
C          Adam Hayes,  (University fo Rochester)
C          Pawel Napiorkowski, (Warsaw University)
C          Nigel Warr,  (University of Cologne)  
C
C       Valuable to Gosia contributions were made by:
C          L. Hasselgren (Uppsala)
C          A.B. Hayes (Rochester)
C          R. Ibbotson (Rochester)
C          A.E. Kavka (Uppsala and Rochester)
C          B. Kotlinski (Warsaw and Rochester)
C          J. Srebrny  (Warsaw)
C          Eh. Vogt (Munchen and Rochester)
C   
C       References and Credits
C          T. Czosnyka, D. Cline and C. Y. Wu,
C          Bull. Am. Phys. Soc. 28, 745 (1983)
C          Internal laboratory report UR/NSRL 308/1986
C 
C          Some concepts used come from the 1978 Winther, de Boer code 
C          COULEX and from the deexcitation code CEGRY developed by Cline 
C          and coworkers at Rochester.  However, the parts taken from 
C          both codes are in most cases completely rewritten, so the 
C          similarity of variable and routine names may be misleading.
C
C       Resources
C          It is recommended that users check the GOSIA website at 
C          Rochester periodically for important updates to the code and 
C          the manual, as well as sample input and output files and other 
C          information. Chapter 11 of this manual provides novice users 
C          with instructions, tutorials, etc.
C
C          http://www.pas.rochester.edu/~cline/Gosia/index.html
C
C          If you need additional information, please contact:
C
C          Prof. Douglas Cline
C          Department of Physics and Astronomy 
C          University of Rochester
C          Rochester, NY 14627, U.S.A.           Phone (585)275-4934
C          Cline@pas.rochester.edu
C          http://www.pas.rochester.edu/~cline/
C
C       Compiling Gosia
C          Gosia compiles on most 64-bit systems with GNU g77, using
C          the default compiler settings.  Previous versions of the
C          Gosia code did not explicitly specify 64-bit precision and
C          were intended to be compiled by the user at the highest
C          machine precision available.  The current availability of
C          64-bit machines and the accuracy problems which may arise
C          when relying on 32-bit precision led to the decision to make
C          this code explicitly 64-bit.  Modifying the code to run
C          with 32-bit precision is discouraged.
C
C       CHRONOLOGY OF MAJOR CHANGES:
C
C          (7 November 2011, N. Warr) gosia-20110524.2
C            - Fix bug in lambda=6 mu=3 collision functions
C          (10 October 2011, N. Warr) gosia-20110524.1
C            - Trap floating point rounding errors in OP,INTI
C          (24 May 2011, A. Hayes) gosia-20110524
C            - Add support for Rachel
C          (24 June 2010, N. Warr) gosia-20081208.10
C            - Fix a bug which uses the wrong sign for the relativistic
C              correction in inverse kinematics
C          (22 February 2010, N. Warr) gosia-20081208.9
C            - Fix a bug in the function TRINT which caused a discontinuity
C          (18 September 2009, N. Warr) gosia-20081208.8
C            - Increased dimensions of ZETA array to allow up to 999 matrix
C              elements (correctly this time, I hope!)
C          (14 September 2009, N. Warr) gosia-20081208.7
C            - Increased dimensions of ZETA array to allow up to 999 matrix
C              elements
C            - Explicit initialisation of variables that were generating
C              warnings on latest gfortran compiler
C          (16 August 2009, N. Warr) gosia-20081208.6
C          Bug fix
C            - Increased dimensions of variables in common VLIN to 101
C          (20 July 2009, P. Napiorkowski) gosia-20081208.5
C          Bug fix
C            - Integration over PIN diodes was incorrect
C          (2 April 2009, N. Warr) gosia-20081208.4
C          Bug fix
C            - E1 polarization was incorrect for projectile excitation
C          (1 February 2009, N. Warr) gosia-20081208.3
C          Bug fix
C            - Preserve IKIN flag in OP,INTI
C          (27 January 2009, N. Warr) gosia-20081208.2
C          Bug fix
C            - Fix OP,INTI for target excitations
C          (16 January 2009, N. Warr) gosia-20081208.1
C          Bug fix
C            - Give an error if limit in omega steps is exceeded
C          (8 December 2008, N. Warr) gosia-20081208
C            - New option OP,INTI which works like OP,INTG but the theta
C              values for the meshpoints are given as angles for the recoiling
C              target nuclei in the laboratory frame. This should help in
C              inverse kinematics cases.
C          (30 June 2008, N. Warr) gosia-20080630
C            - Reordering of all variable declarations so variables in 
C              common blocks are together.
C            - Increased size of arrays for levels from 75 to 100.
C            - Change in many format statements to accomodate 100 levels
C               (I2 -> I3).
C            - Increased size of arrays for substates from 600 to 1200.
C            - Restructure function EFFIX to allow more flexible selection
C              of efficiency calbration type.
C            - Added Radware efficiency calibration (Pawel J. Napiorkowski).
C            - Added EFF, option to CONT to make it easier to select
C              efficiency method.
C            - Bugfix in gremlin efficiency method: a sign was wrong in
C              Woods-Saxon term.
C            - Bugfix in initialisation of conversion coefficients for
C              interpolation method (only part of array was initialized).
C            - Add CONTINUE statements to prevent GOTO the ENDDO of a loop.
C              This is deprecated and will be removed in future.
C            - Call either LAGRAN or SPLNER to do interpolation according
C              to ISPL variable.
C            - Approximation in TRINT  if Arg is very large (then ratio is
C              one) as in pawel.
C            - Replaced OPEN(NAME=xxx) with OPEN(FILE=xxx) as the former is
C              an extension to the language, while the latter is standard.
C            - Bugfix in BRICC: we need the absolute value of the
C              difference in energy levels not the signed difference.
C            - Use JZB to indicate which unit to read from as in gosia2.
C            - Use IUNIT3 to indicate which unit is TAPE3 as in gosia2.
C            - Use irix variable to select units 7 and 12 as in gosia2.
C            - Add dummy variable Op2 to ANGULA as in gosia2.
C            - In OPENF we refuse to open units 25 and 26 unless reading
C              from unit 5 as in gosia2 (however, gosia always reads from unit
C              5, so this does nothing).
C            - Bugfix for dimension of bten, which should have been 1600
C              and was different in different places.
C            - Use DCMPLX instead of CMPLX for portability.
C            - Use DIMAG instead of IMAG for portability.
C            - Reorder DATA statement in NEWCNV for portability.
C            - Use LAGRAN for interpolation if there are less than three
C              data points even if ISPL = 1.
C          (26 June 2008, N. Warr) gosia-20080519.1
C          Bug fixes        
C            - CC array is initialised properly (all 50 elements, not just
C               first 20).
C            - Various GOTOs which go to the ENDDO statement of a loop are
C               changed to use a CONTINUE statement just before the ENDDO.
C               This is apparently no longer permitted in Fortran 2003, so
C               some compilers are already forbidding it.
C            - Arrays are passed with dimension '*' now, so they correctly
C               inherit the
C               dimension from the calling function.
C            - cpo and cpo1 arrays in CONV have been extended from 51 to
C               101.
C            - xx and yy arrays in EFFIX have been extended from 51 to 101.
C            - SAVE arh added in LAGRAN.
C            - The sign of the Woods-Saxon part of the gremlin efficiency
C               calibration has been corrected.
C            - NAME= has been replaced with FILE= in open statements in
C               BRIC as this is an extension, not part of the Fortran standard.
C            - The energy of the gamma is now the absolute value of the
C               difference in energy of the levels, so we don't have negative
C               energies!
C
C          (19 May 2008, N. Warr) gosia-20080519
C            Incorporated the internal conversion coefficient program BrIcc
C            as OP,BRIC
C            Provided option to make CONV use the file generated by OP,BRIC
C            Increased the number of energy levels from 75 to 100.
C            Reordered the declarations for each common block to facilitate
C            building singl program file.
C            Added variable to Gosia to reduce differences between Gosia and
C            Gosia2.
C  
C          (8 May 2008, N. Warr) gosia-20080508
C            Increased number of matrix elements from 500 to 999
C            Added OP,SELE to incorporate the separate program SELECT into
C            Gosia.
C
C          (7 May 2008, N. Warr) gosia-20080507
C            Corrected constants in SEQ
C            Added spline switch SPL to CONT
C
C          (7 May 2008, N. Warr) gosia-20080507
C            Corrected constants in SEQ
C            Added SPL option to CONT
C
C          (18 April 2008, N. Warr) gosia-20080418
C            Increased the size of varios arrays for interpolation  
C    
C          (July 2007, N. Warr) Changes to the Input Format:
C            Tapes 1 and 2 in the pre-2007 versions have been reassigned 
C            to tapes 11 and 12, respectively, in order to make switching 
C            between Gosia and Gosia2 (the mutual excitation code) easier.
C            This change affects only OP,FILE (optional--see below) and 
C            the default naming of files.  For example, on a typical Unix-like
C            system, the file formerly called "fort.2" will now be called 
C            "fort.12" and will contain the best set of matrix elements 
C            written at the end of the minimization.
C     
C          (July 2007, N. Warr) Bugs Fixed  
C            The routine DJMM relies on the values of DJM being
C            preserved between repeated calls to DJMM, which occurred
C            automatically on many older systems (e.g. DEC Alpha, VAX).
C            On some newer machines, DJM was effectively randomized
C            between calls, causing unpredictable errors including
C            negative values of chi-squared.  This was fixed by adding
C            the command "SAVE DJM" to the routine DJMM.
C            The routine ALLOC now handles error conditions gracefully,
C            and execution is halted in the event of a fatal error.
C            The WRN,X. switch in the CONT section of OP,GOSI and 
C            OP,COUL was unintentionally disabled in the 2007 version.
C            It is restored in the present update.
C
C          (July 2007, N. Warr) Explicit 64-bit Precision Upgrade 
C            All routines and functions including the main routine now
C            have "IMPLICIT NONE" declared, and all variables are
C            explicitly defined as either REAL*8, COMPLEX*16, or
C            INTEGER*4.  Numerical constants have been changed as
C            necessary to double precision.  Archaic functions have
C            been updated (below, "Archaic Functions"), in part to
C            preserve 64-bit precision during type-conversions.
C            (Precision in type conversion may be limited in some cases
C            by the compiler.)
C
C          (July 2007, N. Warr) Structure and Standards  
C            Sections of the code have been restructured using Spag 
C            (Polyhedron Software) under the academic license.  This
C            included unraveling of loops and goto statements, and
C            indenting the source code to make loops and if statements
C            more obvious.  The initialization in the main routine has
C            been slightly restructured, mainly to make it similar to
C            the 2007 version of Gosia2.  Other sections have been
C            restructured for clarity, without altering their function
C            (e.g. WTHREJ).
C
C          (July 2007, N. Warr) Common Blocks 
C            The common blocks ME2D, CCC, KIN, COEX, CAUX0, and LCZP
C            were re-ordered so that the 64-bit real variables come
C            before the 32-bit integer variables in order to
C            eliminate alignment problems.  Several unused common
C            blocks were removed from routines.
C
C          (July 2007, N. Warr) Archaic Functions
C            All instances of the following archaic functions have
C            been replaced by their modern counterparts.
C
C                Archaic    Replacement        Archaic    Replacement
C                IFIX       INT                MIN0       MIN  
C                FLOAT      REAL               AMIN1      MIN  
C                IABS       ABS                ALOG10     LOG10 
C                MAX0       MAX                ALOG       LOG     
C                AMAX1      MAX                
C         
C          (June 2006, T. Czosnyka) - The size of the array of
C            internal conversion coefficients (CC) has been increased
C            to 50.
C
C          (Nov 2000, T. Czosnyka) - A Jaeri efficiency calibration
C            has been added.
C
C          (2000) - A FITEFF efficiency calibration has been added
C            with credit to P. Olbratowski, P. Napiorkowski.
C
C          (July 1997, T. Czosnyka) - Known matrix elements of all
C            multipolarities may now be entered as data in OP,YIEL.
C            Note that this necessitates adding the multipole order
C            LAMBDA as the first field in the new input format:
C            LAMBDA, NS1, NS2, ME, DME 
C            where LAMBDA=7 and 8 represent M1 and M2, respectively.
C
C          (September 1996, T. Czosnyka) - The PIN diode particle
C            detector option has been added.  See the entry for
C            "PIN,X."  under the sub-option CONT in the Gosia manual.
C
C          (May 1995, T. Czosnyka) - Added a matrix element generator
C            "OP,THEO" following the "general structure of matrix
C            elements" as given in Bohr & Mottelson vol. II.  Refer to
C            the Gosia manual.
C
C          (April 1991, T. Czosnyka) - The OP,RAW function has been
C            added.  OP,RAW handles non-efficiency-corrected spectra
C            and allows the definition of Ge detector "clusters."  Up
C            to 20 clusters can be defined.  This effectively increases
C            the number of physical Ge detectors to 200, while the
C            number of data sets (i.e. single detectors + cluster
C            detectors) is still limited to 32.
C
C          (April 1991, T. Czosnyka) - Output is now written on unit
C            22 to avoid mixing it with system messages on some
C            systems.
C
C          (November 1990, T. Czosnyka) - The level scheme data
C            arrays have been increased to the following sizes:
C            number of levels   = 75
C            gamma-ray yields   = 32 x 1500 
C            magnetic substates = 600 
C            matrix elements    = 500
C
C
C          (April 1990, T. Czosnyka) - The dominant multipolarity
C            switch is now ignored by the code and does not need to be
C            set.  Full Q-maps are now calculated for electric
C            multipole orders E1 through E4.  The electric matrix
C            elements up to multipole order E6 may be entered and fit.
C            The Xi and Zeta function ranges are now calculated for
C            each multipolarity individually.
C
C          (1990, Eh. Vogt, T. Czosnyka) - OP,FILE has been added,
C            giving the user the option of specifying descriptive names
C            of the input and output files in the Gosia input, rather
C            than using the Fortran default names fort.1, fort.2, etc.
C            Refer to the Gosia website for sample input files that use
C            OP,FILE.
C
C          (March 1989, T. Czosnyka) - The code has been updated to
C            allow input of data from 32 Ge detectors.  [As of the 2007
C            version, this means a total of 32 X 1500 data points.]
C
C          (1980, T. Czosnyka, D. Cline, C.Y. Wu) - Original version.
C
C---------------------------------------------------------------------------
C PROGRAM GOSIA
C
C Calls: ADHOC, ALLOC, ANGULA, ARCCOS, ARCTG, CMLAB, COORD, DECAY, DJMM,
C        EFFIX, ELMT, FAKP, FHIP, FTBM, INTG, INVKIN, KLOPOT, KONTUR, LAGRAN,
C        LOAD, MINI, MIXR, MIXUP, OPENF, PATH, PRELM, PTICC, QFIT, READY,
C        SETIN, SIMIN, SNAKE, SPLNER, STING, TACOS, TAPMA, TEMB, TENS, WSIXJ,
C        WTHREJ
C
C Uses global variables:
C      ABC    - absorption coefficients
C      ACCA   - accuracy
C      ACCUR  - accuracy required
C      AGELI  - angles of Ge detectors
C      AKAVKA - efficiency curve parameters
C      ARM    - excitation amplitudes of substates.
C      AVJI   - average J (N.B. here it is G(1))
C      B      - table of factorials
C      BEQ    - identifier for angle for rotations
C      BETAR  - recoil beta
C      CAT    - substates of levels (n_level, J, m)
C      CC     - conversion coefficients
C      CNOR   - normalization factors
C      CORF   - internal correction factors
C      DEVD   -
C      DEVU   -
C      DIPOL  - E1 polarization parameter
C      DIX    - Ge parameters (inner & outer radius, length, distance)
C      DLOCK  - limit derivative below which matrix element is fixed if LOCKS=1
C      DS     - integrated rutherford cross-section
C      DSE    - rutherford cross section at given energy integrated over angles
C      DSG    - differential gamma-ray yield at meshpoints
C      DSIGS  - dsigma for each experiment
C      DYEX   - error on experimental yield
C      EAMX   - known matrix elements and their errors
C      ELM    - matrix elements
C      ELMH   -
C      ELML   - lower limit on matrix elements
C      ELMU   - upper limit on matrix elements
C      EMMA   - Controls number of magnetic substates in full coulex calc.
C      EN     - energy of level
C      EP     - bombarding energy
C      ERR    - error flag
C      EXPO   - adiabatic exponential
C      FIEL   - K (N.B. here it is G(6))
C      FIEX   - phi range of particle detector
C      GAMMA  - Gamma (N.B. here it is G(2))
C      GFAC   - g (N.B. here it is G(5))
C      GRAD   - partial derivative of chi squared wrt. each matrix element
C      HLM    - matrix elements before minimisation
C      HLMLM  - old value of matrix element or chi squared
C      IAMX   - index of matrix element for known matrix element
C      IAMY   - level indices of pair of levels for which matrix element is known
C      IAX    - axial symmetry flag
C      IBYP   - flag to indicate whether we calculate <\alpha_k>
C      ICLUST - cluster number for each experiment and detector
C      ICS    - read internal correction factors flag (OP,CONT switch CRF,)
C      IDIVE  - number of subdivisions
C      IDRN   - index of normalising transition for yields
C      IEXP   - experiment number
C      IFAC   - spin/parity phase factor
C      IFBFL  - calculate derivatives with forward-backward method
C      IFMO   - include correction to angular distance for finite recoil distance.
C      ILE    - yield number for each detector
C      IMIN   -
C      INHB   - inhibit error flag (LERF) setting in POMNOZ
C      INNR   - independent normalisation switch (see OP,CONT INR,)
C      INTERV - default accuracy check parameter for Adams-Moulton (see OP,CONT:INT)
C      INTR   - flag to swap chisqr and log(chisqr)
C      IP     - table of prime numbers
C      IPRM   - various flags to control output
C      IPS1   - terminate after calculating and storing internal correction factors
C      IRAWEX - flag to indicate raw uncorrected yield
C      ISEX   -
C      ISKIN  - kinematic flag (0,1)
C      ISMAX  - number of substates used
C      ISO    - isotropic flag
C      ITMA   - identify detectors according to OP,GDET
C      ITS    - create tape 18 file (OP,CONT switch SEL,)
C      ITTE   - thick target experiment flag
C      IUNIT3 - unit for TAPE3
C      IVAR   - indicates a limit or correlation is set
C      IWF    - warning flag
C      IY     - index for yields
C      IZ     - Z of investigated nucleus
C      IZ1    - Z of non-investigated nucleus
C      JENTR  - flag set to 0 normally, 1 in OP,ERRO
C      JSKIP  - Experiments to skip during minimisation.
C      JZB    - unit to read from
C      KFERR  - error flag for minimization
C      KSEQ   - index of level
C      KVAR   -
C      LAMAX  - number of multipolarities to calculate
C      LAMBDA - list of multipolarities to calculate
C      LASTCL - index of last detector in cluster
C      LDNUM  - number of matrix elements with each multipolarity populating each level
C      LEAD   - pair of levels involved in each matrix element
C      LIFCT  - index for lifetimes
C      LMAX   - ground-state spin + 1
C      LMAXE  - maximum multipolarity needed for calculation
C      LNORM  - normalisation constant control
C      LNY    - use logs to calculate chi squared
C      LOCKF  - flag to fix matrix elements with most significant derivative
C      LOCKS  - lock flag. If LOCKS=1, fix at first stage of minimization
C      LP1    - maximum number of experiments (50)
C      LP10   - maximum number of substates (1200)
C      LP11   - LP8 - 1 (2800)
C      LP12   - number of steps of omega (365)
C      LP13   - LP9 + 1 (47901)
C      LP14   - maximum space for collision functions (4900)
C      LP2    - maximum number of matrix elements (1500)
C      LP3    - maximum number of levels (100)
C      LP4    - maximum number of yields (1500)
C      LP6    - maximum number of gamma detectors (32)
C      LP7    - start of collision functions (45100)
C      LP8    - (104)
C      LP9    - length of ZETA - 2100 (47900)
C      MAGA   - number of magnetic substates in approximate calculation
C      MAGEXC - flag: 0 means no magnetic excitations, 1 means with mag. exc.
C      MEMAX  - number of matrix elements
C      MEMX6  - number of matrix elements with E1...6 multipolarity
C      MULTI  - number of matrix elements having given multipolarity
C      NAMX   - number of known matrix elements
C      NANG   - number of gamma-ray detectors for each experiment
C      NBRA   - number of branching ratios
C      NCM    - calculate kinematics assuming this state for final state (default = 2)
C      NDIM   - maximum number of levels
C      NDST   - number of data sets
C      NEXPT  - number of experiments
C      NLIFT  - number of lifetimes
C      NLOCK  - number of elemnts to fix if LOCKF=1
C      NMAX   - number of levels
C      NMAX1  - number of levels with decays
C      NYLDE  - number of yields
C      ODL    - results of OP,GDET calculation
C      PARX   - [for maps]
C      PARXM  - [for maps]
C      POWER  - x (N.B. here it is G(7))
C      QAPR   - approximate Coulomb amplitudes
C      SA     - ratio of matrix elements for correlated elements
C      SE     - seed for random number generator of OP,RAND
C      SGW    - number of standard deviations to generate warning (see control option WRN,X)
C      SPIN   - spin of level
C      SUBCH1 - partial chisqr
C      SUBCH2 - partial chisqr
C      SUMCL  - sum of yields for clusters
C      TAU    - lifetime in picoseconds
C      THICK  - thickness of each absorber type
C      TIMEC  - Tau_C (N.B. here it is G(4))
C      TIMEL  - lifetimes and their errors
C      TLBDG  - theta of particle detector in lab frame (in degrees)
C      TREP   - theta of recoiling nucleus (in radians)
C      UPL    - upper limits for all gamma detectors
C      VINF   - speed of projectile at infinty
C      XA     - A of investigated nucleus
C      XA1    - A of non-investigated nucleus
C      XI     - xi coupling coefficients
C      XIR    - [for maps]
C      XLAMB  - Lambda* (N.B. here it is G(3))
C      XV     - energy meshpoints (sometimes theta meshpoints) where we calculate exact Coulex
C      YEXP   - experimental yields
C      YGN    - gamma yield calculated without correction to angular distribution from finite recoil distance
C      YGP    - gamma yield calculated with correction to angular distribution from finite recoil distance
C      YNRM   - relative normalisation for gamma detectors
C      YV     - scattering angle meshpoints where we calculate exact Coulex
C      ZETA   - various coefficients
C      ZPOL   - dipole term (GDR excitation)
C      ZV     - energy meshpoints

      PROGRAM GOSIA
      IMPLICIT NONE
      REAL*8 acof , ap , ARCCOS , ARCTG , arg , ax , bcof , be2 , 
     &       be2a , be2b , be2c
      REAL*8 bk , bl , bm , bmx , bten , bu , ccc , 
     &       ccd , cf , chilo , chiok , chis0 , chisl , chisq , chiss , 
     &       cnst
      REAL*8 cocos , conu , d , decen , dedx , dsd , dsig , dst
      REAL*8 dsx , dsxm , effi , eh1 , elmi , ELMT , emhl1 , emn , emx ,
     &       enb
      REAL*8 eng , enh , esd , esp , ess , 
     &       fi0 , fi1 , fic , fiex1 , figl , fipo1 , fm , gth
      REAL*8 hen , het , p , pfi , 
     &       ph1 , ph2 , pi , po1 , po2 , polm , pop1 , pr , pv
      REAL*8 q1 , q2 , qc , qfac , qr , qui , r , r1 , r2 , r3 , r4 , 
     &       rem , remax , rl , rlr , rm , rx , ry
      REAL*8 rz , s , s11 , s12 , s21 , s22 , sbe , sf , sh , sh1 , 
     &       sh2 , SIMIN , slim
      REAL*8 summm , sz1 , sz2 , TACOS , tau1 , tau2 , test , 
     &       tetrc , tfac , thc , title , tmn , tmx , todfi
      REAL*8 tta , tth , tting , ttttt , txx , u , 
     &       val , waga , wph , wpi , WSIXJ , wth , wthh , 
     &       WTHREJ
      REAL*8 xep , xi1 , xi2 , xk1 , xk2 , xl1 , xlevb , 
     &       xlk , xm1 , xm2 , xm3 , xtest , xw , xx , xxi , 
     &       ycorr
      REAL*8 yy , yyd1 , yydd , yyy , zmir , zp , zz
      REAL*8 ttttx ! Only gosia1 and pawel
      INTEGER*4 i , i122 , iapx , ib , ibaf , icg , icll , ict , ictl , 
     &          id , ideff , idf
      INTEGER*4 idr , iecd , ient , ifbp , ifc , ifm , ifwd , 
     &          ig1 , ig2 , ih1 , ih2 , ihlm , ihuj , ii , ij
      INTEGER*4 ija0 , ijaja , ijan , ijk , ijx , ile1 , ilevls , 
     &          ilx , im , imode , in1 , in2 , inclus , ind , 
     &          ind1 , ind2 , indx
      INTEGER*4 inko , inm1 , inm2 , inn , inpo , intend , intvh , 
     &          inva , inx1 , iobl , iocc , iopri , iosr , ipd , iph
      INTEGER*4 ipine , ipinf , ipo1 , ipo2 , ipo3 , ipp , iprc , 
     &          ipri , irea , irep , irfix , irix , isip , iske , iskf
      INTEGER*4 isko , iskok , isoh , ispa , ispb , itno , 
     &          itp , iuy , iva , iva1 , ivarh , ivari , ivrh
      INTEGER*4 ixj , ixl , ixm , iyr , izcap , j , ja , 
     &          jan , jan1 , jb , jb1 , jb2 , jd , jde , jdy , je
      INTEGER*4 jex , jexp , jfi , jfre , jgd , jgl , jgl1 , jgr , jgs ,
     &          jj , jj1 , jjjj , jjlx , jjx , jk , jkloo , jktt , jl , 
     &          jmm , jmpin
      INTEGER*4 jp , jphd , jpin , jrls , js , jt , jtp , jyi , jyi1 , 
     &          jyi2 , jyv , jz , k , kb , kclust , kerf , kex
      INTEGER*4 kh , kh1 , kh2 , kk , kk1 , kk2 , kkk , kl , kloop , 
     &          kmat , kq , ktt , kuku , l , la , la1 , lam , lamd
      INTEGER*4 lamh , lb , lck1 , lck2 , levl , lex , lexp , 
     &          lfagg , lfini , lh1 , lh2 , liscl , lkj
      INTEGER*4 lkj1 , ll , lli , lll , lmax1 , lmaxh , locat , 
     &          loct , lp0 , lpin
      INTEGER*4 ltrn , ltrn1 , ltrn2 , lu , lx , lxd , magh , MEM
      INTEGER*4 memax1 , memh , memx4 , mend , mexl , 
     &          mfla , mlt , mm , mpin , ms , n , na , na1 , naa , 
     &          nallow
      INTEGER*4 naxfl , nb1 , nb2 , nbands , nch , ndima , ndum , 
     &          ne , nf , nfd , nfdd , 
     &          nfi , nflr , nft , nged
      INTEGER*4 ngpr , ni , nksi , nl , nmaxh , nmemx , nnl , 
     &          nogeli , npce , npce1 , npct , npct1 , 
     &          npt , nptl , nptx , ns1
      INTEGER*4 ns2 , ntap , ntt , numcl , nval , nz
      INTEGER*4 iskin_protect
      CHARACTER*4 oph , op1 , opcja , op2
      CHARACTER*1 prp
      DIMENSION ihlm(32) , esp(20) , dedx(20) , bten(1600) , ! bten dimension = 16 * maxlevels
     &          fiex1(100,100,2) , title(20) , pfi(101) , zmir(6,2,50) ,
     &          iecd(50) , wpi(100,2) , tau1(10) , eng(10) , 
     &          tau2(10,7) , xl1(7) , qui(8,10) , cf(8,2) , 
     &          ivarh(1500) , liscl(200) , dsxm(100,100,100) , 
     &          levl(50) , xlevb(50,2) , bm(8,20,20,3) , mlt(1500) , 
     &          ivari(1500) , jpin(50) , ideff(50) , iskin_protect(50)
      INTEGER*4 ICLUST , LASTCL , IRAWEX
      REAL*8 SUMCL      
      COMMON /CLUST / ICLUST(50,200) , LASTCL(50,20) , SUMCL(20,1500) , 
     &                IRAWEX(50)
      INTEGER*4 NDST
      COMMON /CCCDS / NDST(50)
      INTEGER*4 INHB
      COMMON /INHI  / INHB
      REAL*8 BEQ
      COMMON /IDENT / BEQ
      REAL*8 ABC, AKAVKA, THICK
      COMMON /EFCAL / ABC(8,10) , AKAVKA(9,200) , THICK(200,7)
      REAL*8 TETACM, TREP, DSIGS
      COMMON /TCM   / TETACM(50) , TREP(50) , DSIGS(50)
      REAL*8 BETAR
      COMMON /BREC  / BETAR(50)
      COMPLEX*16 EXPO
      COMMON /ADBXI / EXPO(1500)
      REAL*8 DIX, ODL
      COMMON /DIMX  / DIX(4) , ODL(200)
      REAL*8 DELTA, ENDEC, ENZ
      INTEGER*4 ITMA
      COMMON /TRA   / DELTA(1500,3) , ENDEC(1500) , ITMA(50,200) , 
     &                ENZ(200)
      REAL*8 CNOR
      INTEGER*4 INNR
      COMMON /CINIT / CNOR(32,100) , INNR
      REAL*8 SE
      COMMON /XRA   / SE
      REAL*8 HLM
      COMMON /HHH   / HLM(1500)
      REAL*8 VACDP, QCEN, DQ, XNOR, AKS
      INTEGER*4 IBYP
      COMMON /VAC   / VACDP(3,100) , QCEN , DQ , XNOR , AKS(6,100) ,
     &                IBYP
      REAL*8 EAMX
      INTEGER*4 NAMX, IAMX, IAMY
      COMMON /ME2D  / EAMX(100,2) , NAMX , IAMX(100) , IAMY(100,2)
      REAL*8 TIMEL
      INTEGER*4 LIFCT
      COMMON /LIFE1 / LIFCT(50) , TIMEL(2,50)
      REAL*8 DEVD, DEVU
      COMMON /DFTB  / DEVD(1500) , DEVU(1500)
      INTEGER*4 KFERR
      COMMON /ERRAN / KFERR
      INTEGER*4 LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &          LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      INTEGER*4 ISKIN
      COMMON /SECK  / ISKIN(50)
      REAL*8 XV, YV, ZV, DSG, DSE, DS
      COMMON /VLIN  / XV(101) , YV(101) , ZV(101) , DSG(101) ,
     &                DSE(101) , DS
      REAL*8 GRAD , HLMLM , ELMH
      COMMON /DUMM  / GRAD(1500) , HLMLM(1500) , ELMH(1500)
      REAL*8 BRAT
      INTEGER*4 IBRC , NBRA
      COMMON /BRNCH / BRAT(50,2) , IBRC(2,50) , NBRA
      REAL*8 YEXP, CORF , DYEX , UPL , YNRM
      INTEGER*4 IY , NYLDE , IDRN , ILE
      COMMON /YEXPT / YEXP(32,1500) , IY(1500,32) , CORF(1500,32) , 
     &                DYEX(32,1500) , NYLDE(50,32) , UPL(32,50) , 
     &                YNRM(32,50) , IDRN , ILE(32)
      REAL*8 YGN , YGP
      INTEGER*4 IFMO
      COMMON /YTEOR / YGN(1500) , YGP(1500) , IFMO
      REAL*8 TAU
      INTEGER*4 KSEQ
      COMMON /LEV   / TAU(100) , KSEQ(1500,4)
      REAL*8 PARX, PARXM, XIR
      COMMON /MAP   / PARX(50,12,5) , PARXM(50,4,10,6) , XIR(6,50)
      REAL*8 EG, CC, AGELI, Q
      INTEGER*4 NICC, NANG, ISPL
      COMMON /CCC   / EG(50) , CC(50,5) , AGELI(50,200,2) , Q(3,200,8) ,
     &                NICC , NANG(200) , ISPL
      REAL*8 G(7) , AVJI , GAMMA , XLAMB , TIMEC , GFAC , FIEL , POWER
      COMMON /GGG   / AVJI , GAMMA , XLAMB , TIMEC , GFAC , FIEL , POWER
      EQUIVALENCE(AVJI,G)
      COMPLEX*16 ARM
      COMMON /AZ    / ARM(1200,7)
      REAL*8 EPS, EROOT, FIEX
      INTEGER*4 IEXP, IAXS
      COMMON /KIN   / EPS(50) , EROOT(50) , FIEX(50,2) , IEXP , IAXS(50)
      REAL*8 XI
      COMMON /CXI   / XI(1500)
      INTEGER*4 LAMDA , LEAD , LDNUM , LAMMAX , MULTI
      COMMON /CLCOM / LAMDA(8) , LEAD(2,1500) , LDNUM(8,100) , LAMMAX , 
     &                MULTI(8)
      REAL*8 EN , SPIN , ACCUR , DIPOL , ZPOL , ACCA
      INTEGER*4 ISO
      COMMON /COEX  / EN(100) , SPIN(100) , ACCUR , DIPOL , ZPOL , 
     &                ACCA , ISO
      INTEGER*4 IMIN , LNORM
      COMMON /MINNI / IMIN , LNORM(50)
      REAL*8 XA , XA1 , EP , TLBDG , VINF
      INTEGER*4 NEXPT , IZ , IZ1
      COMMON /CX    / NEXPT , IZ , XA , IZ1(50) , XA1(50) , EP(50) , 
     &                TLBDG(50) , VINF(50)
      INTEGER*4 MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR
      COMMON /CEXC  / MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR(1500)
      INTEGER*4 IPRM
      COMMON /PRT   / IPRM(20)
      REAL*8 ZETA
      INTEGER*4 LZETA
      COMMON /CCOUP / ZETA(155600) , LZETA(8)
      REAL*8 B
      COMMON /CB    / B(20)
      INTEGER*4 LMAX
      COMMON /CLM   / LMAX
      INTEGER*4 IFAC
      COMMON /CLCOM0/ IFAC(100)
      REAL*8 CAT
      INTEGER*4 ISMAX
      COMMON /CLCOM8/ CAT(1200,3) , ISMAX
      LOGICAL ERR
      COMMON /CLCOM9/ ERR
      REAL*8 ELM , ELMU , ELML , SA
      COMMON /COMME / ELM(1500) , ELMU(1500) , ELML(1500) , SA(1500)
      INTEGER*4 NMAX , NDIM , NMAX1
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      INTEGER*4 INTERV
      COMMON /CEXC9 / INTERV(50)
      REAL*8 EMMA
      INTEGER*4 NCM
      COMMON /CAUX0 / EMMA(100) , NCM
      INTEGER*4 IPATH , MAGA
      COMMON /PTH   / IPATH(100) , MAGA(100)
      REAL*8 QAPR
      INTEGER*4 IAPR , ISEX
      COMMON /APRCAT/ QAPR(1500,2,7) , IAPR(1500,2) , ISEX(100)
      REAL*8 SGW , SUBCH1 , SUBCH2
      INTEGER*4 IWF
      COMMON /WARN  / SGW , SUBCH1 , SUBCH2 , IWF
      INTEGER*4 ITTE
      COMMON /THTAR / ITTE(50)
      REAL*8 DLOCK
      INTEGER*4 LOCKF , NLOCK , IFBFL , LOCKS
      COMMON /FIT   / LOCKF , NLOCK , IFBFL , LOCKS , DLOCK
      INTEGER*4 LERF , IDIVE
      COMMON /APRX  / LERF , IDIVE(50,2)
      INTEGER*4 JSKIP
      COMMON /SKP   / JSKIP(50)
      INTEGER*4 ITS
      COMMON /TRB   / ITS
      INTEGER*4 KVAR
      COMMON /SEL   / KVAR(1500)
      INTEGER*4 JENTR , ICS
      COMMON /ERCAL / JENTR , ICS
      INTEGER*4 LNY , INTR , IPS1
      COMMON /LOGY  / LNY , INTR , IPS1
      REAL*8 PILOG
      INTEGER*4 IP , IPI , KF
      COMMON /FAKUL / IP(26) , IPI(26) , KF(101,26) , PILOG(26)
      INTEGER*4 NLIFT
      COMMON /LIFE  / NLIFT
      INTEGER*4 IUNIT3 , JZB
      COMMON /SWITCH/ JZB , IUNIT3
      DATA (eng(k),k=1,10)/.05 , .06 , .08 , .1 , .15 , .2 , .3 , .5 , 
     &      1. , 1.5/
C     Absorption coefficients in units of 1/cm for Ge
      DATA (tau1(k),k=1,10)/17.656 , 10.726 , 5.076 , 2.931 , 1.3065 , 
     &      .8828 , .5959 , .4357 , .3041 , .2472/
C     Absorption coefficients in units of 1/cm for Al, C, Fe, Cu, Ag/Cd/Sn, Ta
C     and Pb at the energies 0.05, 0.06, 0.08, 0.1, 0.15, 0.2, 0.3, 0.5, 1, 1.5
C     MeV
      DATA (tau2(k,1),k=1,10)/.9883 , .7473 , .5442 , .4592 , .3718 , 
     &      .3302 , .2814 , .2278 , .1657 , .1350/
      DATA (tau2(k,2),k=1,10)/1.014 , .7443 , .5195 , .4261 , .3362 , 
     &      .2967 , .2518 , .2038 , .1479 , .1204/
      DATA (tau2(k,3),k=1,10)/15.167 , 9.405 , 4.652 , 2.889 , 1.525 , 
     &      1.135 , .8643 , .6592 , .4703 , .3830/
      DATA (tau2(k,4),k=1,10)/23.184 , 14.182 , 6.777 , 4.059 , 1.970 , 
     &      1.384 , .9936 , .7473 , .5274 , .4297/
      DATA (tau2(k,5),k=1,10)/84.351 , 51.445 , 23.822 , 13.070 , 
     &      4.774 , 2.605 , 1.339 , .7925 , .5005 , .4032/
      DATA (tau2(k,6),k=1,10)/93.364 , 58.559 , 125.96 , 70.713 , 
     &      25.302 , 12.541 , 5.193 , 2.215 , 1.077 , .8176/
      DATA (tau2(k,7),k=1,10)/89.809 , 56.338 , 27.009 , 62.966 , 
     &      22.933 , 11.334 , 4.540 , 1.813 , .8020 , .5900/
      DATA q1/0./,q2/0./,iph/0/
      DATA cnst/0./,sh1/0./,irfix/0/,jfre/0/ ! Only gosia1 and pawel

C     Initialize prime numbers
      IP(1) = 2
      IP(2) = 3
      IP(3) = 5
      IP(4) = 7
      IP(5) = 11
      IP(6) = 13
      IP(7) = 17
      IP(8) = 19
      IP(9) = 23
      IP(10) = 29
      IP(11) = 31
      IP(12) = 37
      IP(13) = 41
      IP(14) = 43
      IP(15) = 47
      IP(16) = 53
      IP(17) = 59
      IP(18) = 61
      IP(19) = 67
      IP(20) = 71
      IP(21) = 73
      IP(22) = 79
      IP(23) = 83
      IP(24) = 89
      IP(25) = 97
      IP(26) = 101

C     Initialize pointers
      lp0 = 155600 ! Size of ZETA array
      LP1 = 50 ! Maximum number of experiments
      LP2 = 1500 ! Maximum number of matrix elements
      LP3 = 100 ! Maximum number of levels
      LP4 = 1500
      LP6 = 32 ! Maximum number of gamma detectors
      LP7 = lp0 - 4900 ! Start of collision coefficients in ZETA
      LP8 = LP3*28 + 1
      LP9 = lp0 - LP3*28
      LP10 = 1200 ! Maximum number of substates
      LP11 = LP8 - 1
      LP12 = 365 ! Maximum number of steps of omega (dimension of ADB, SH, CH)
      LP13 = LP9 + 1
      LP14 = 4900 ! Maximum number of collision coefficients

      JZB = 5

C     Initialize normalization to 1.
      DO i = 1 , LP3 ! LP3 = 100 (maximum number of levels)
         DO j = 1 , LP6 ! LP6 = 32 (maximum number of gamma detectors)
            CNOR(j,i) = 1.
         ENDDO
      ENDDO

      IUNIT3 = 3 ! Is 33 in gosia2
      IBYP = 0
      INHB = 0
      BEQ = -983872.
      ipinf = 0
      iyr = 0
      pi = 3.141592654
      INNR = 0
      itno = 0
      chisq = 0.
      chilo = 0.
      IWF = 1 ! Turn on warnings
      ifm = 0 ! Fast minimisation switch off by default
      IPS1 = 11
      ifwd = -1
      INTR = 0
      LNY = 0
      JENTR = 0 ! Flag to indicate we are not in OP,ERRO
      ICS = 0
      ISPL = 0 ! Flag to indicate we should use LAGRAN not SPLNER

      DO i = 1 , LP1 ! LP1 = 50 (maximum number of experiments)
         ideff(i) = 0
         jpin(i) = 0
         iecd(i) = 0
      ENDDO
      txx = 0.
      SGW = 3.
      SUBCH1 = 0.
      SUBCH2 = 0.
      ITS = 0 ! Create tape 18 flag
      iosr = 0
      LOCKS = 0
      DLOCK = 1.1
      kerf = 0
      IFBFL = 0
      NLOCK = 0
      LOCKF = 0
      DO i = 1 , LP4 ! LP4 = 1500
         DO j = 1 , LP6 ! LP6 = 32 (maximum number of gamma detectors)
            CORF(i,j) = 1.
         ENDDO
      ENDDO
      DO i = 1 , 20
         IPRM(i) = 1
      ENDDO
      DO i = 1 , 50
         DO j = 1 , 5
            CC(i,j) = 0.
         ENDDO
      ENDDO
      IPRM(4) = -2
      IPRM(5) = 11111
      IPRM(6) = 11111
      IPRM(7) = 0
      IPRM(16) = 0
      IPRM(17) = 0
      IPRM(18) = 0
      IPRM(19) = 0
      IPRM(20) = 0
      DO i = 1 , LP1 ! LP1 = 50 (maximum number of experiments)
         DO j = 1 , 5
            IF ( j.NE.5 ) THEN
               DO k = 1 , 10
                  DO kuku = 1 , 6
                     PARXM(i,j,k,kuku) = 0.
                  ENDDO
               ENDDO
            ENDIF
            DO k = 1 , 12
               PARX(i,k,j) = 0.
            ENDDO
         ENDDO
      ENDDO
      DO k = 1 , LP1 ! LP1 = 50 (maximum number of experiments)
         IDIVE(k,1) = 1
         IDIVE(k,2) = 1
         DO iuy = 1 , 6
            XIR(iuy,k) = 0.
         ENDDO
      ENDDO
      iobl = 0
      lfagg = 0
      izcap = 12800
      KFERR = 0
      NDIM = LP3 ! LP3 = 100 (maximum number of levels)
      ISO = 1
      B(1) = 1.
      DO i = 2 , 20
         B(i) = B(i-1)*(i-1)
      ENDDO
      LMAXE = 0
      CALL FAKP
      CALL FHIP
      NCM = 2 ! Default final state for kinematics calculation (OP,CONT NCM,)
      DO ijx = 1 , LP1 ! LP1 = 50 (maximum number of experiments)
         INTERV(ijx) = 1
      ENDDO
      la = 0
      ipo3 = 1
      indx = 0
      ACCUR = .00001
      icg = 1
      ient = 1
      jphd = 1 ! Print header flag
      DIPOL = 0.005
      MAGEXC = 0 ! Initially flag that we don't need magnetic excitations
      LAMMAX = 0
      DO lam = 1 , 8
         DO lexp = 1 , LP3 ! LP3 = 100 (maximum number of levels)
            LDNUM(lam,lexp) = 0
         ENDDO
         MULTI(lam) = 0
         LAMDA(lam) = 0
      ENDDO
      DO j = 1 , LP2 ! LP2 = 1500 (maximum number of matrix elements)
         EXPO(j) = (1.,0.)
         KVAR(j) = 1
         ELM(j) = 0.
      ENDDO
      DO j = 1 , LP1 ! LP1 = 50 (maximum number of experiments)
         JSKIP(j) = 1
         ISKIN(j) = 0
      ENDDO
      DO j = 1 , LP3 ! LP3 = 100 (maximum number of levels)
         ISEX(j) = 1111
      ENDDO
      ISEX(1) = 0
      ACCA = .00001
      oph = '    '
      nmemx = LP2 + 9 ! LP2 = 1500 (maximum number of matrix elements)
      IEXP = 1
      IMIN = 0
      i122 = 0
      DO j = 1 , LP2 ! LP2 = 1500 (maximum number of matrix elements)
         DO k = 1 , 2
            DO l = 1 , 7
               QAPR(j,k,l) = 0.
            ENDDO
         ENDDO
      ENDDO
      ERR = .FALSE.
      opcja = '    '
      intend = 0 ! End of initialization

C.............................................................................
C     Start reading input file.
 100  READ (JZB,99001) op1 , op2
99001 FORMAT (1A3,1A4)
      
      IF ( op1.EQ.'OP, ' ) THEN
         IF ( op2.EQ.'GOSI' ) oph = op2
         IF ( op2.EQ.'GOSI' ) opcja = op2

C        Treat OP,FILE (attach files to fortran units)
         IF ( op2.EQ.'FILE' ) THEN
            CALL OPENF
            GOTO 100 ! End of OP,FILE - back to input loop
         ENDIF

C        Print header         
         IF ( jphd.EQ.1 ) WRITE (22,99002)
99002    FORMAT ('1'/1X,125('*')/1X,125('*')/1X,50('*'),25X,50('*')/1X,
     &           50('*'),10X,'GOSIA',10X,50('*')/1X,50('*'),25X,50('*')
     &           /1X,125('*')/1X,125('*')////)
         IF ( jphd.EQ.1 ) WRITE (22,99003)
99003    FORMAT (1X/20X,'ROCHESTER COULOMB EXCITATION DATA ANALYSIS ',
     &           'CODE BY T.CZOSNYKA,D.CLINE AND C.Y.WU'/50X,
     &           'VERS. 20110524.2   NOV 2011'//////)
         jphd = 0 ! Set print header flag to zero, so we don't repeat header

C        Handle OP,GDET (germanium detectors)
         IF ( op2.EQ.'GDET' ) THEN
            nl = 7
            READ (JZB,*) nfdd ! number of physical detectors

            nfd = ABS(nfdd) ! Negative value means graded absorber
            IF ( nfdd.LE.0 ) THEN
               REWIND 8
               DO i = 1 , nl
                  WRITE (8,*) (tau2(l,i),l=1,10)
               ENDDO
               WRITE (8,*) (eng(l),l=1,10)
            ENDIF

C           Write file for gamma-ray energy dependence of Ge solid-angle
C           attenuation coefficients
            REWIND 9
            WRITE (9,*) nfd
            DO i = 1 , nfd ! For each detector
               READ (JZB,*) (DIX(k),k=1,4) ! radius of core, outer radius, length, distance
               READ (JZB,*) (xl1(k),k=1,nl) ! thicknesses of 7 kinds of absorber
               IF ( DIX(1).LE.0. ) DIX(1) = .01
               WRITE (9,*) DIX(4) ! length
               IF ( nfdd.LE.0 ) WRITE (8,*) (xl1(k),k=1,nl)
               ind = 1
               IF ( xl1(5).GT.0. ) ind = 3
               IF ( xl1(6).GT.0. ) ind = 4
               IF ( xl1(7).GT.0. ) ind = 5
               WRITE (9,*) eng(ind) ! First energy
               CALL QFIT(qui,tau1,tau2,eng,xl1,cf,nl,ind)
               WRITE (22,99004) i
99004          FORMAT (10X,'DETECTOR',1X,1I2)
               DO k = 1 , 8
                  WRITE (22,99005) k , cf(k,1) , cf(k,2)
99005             FORMAT (1X,//5X,'K=',1I1,2X,'C1=',1E14.6,2X,'C2=',
     &                    1E14.6/5X,'ENERGY(MEV)',5X,'FITTED QK',5X,
     &                    'CALC.QK',5X,'PC.DIFF.'/)
                  WRITE (9,*) cf(k,1) , cf(k,2) , qui(k,ind)
                  DO l = 1 , 10
                     arg = (eng(l)-eng(ind))**2
                     qc = (qui(k,ind)*cf(k,2)+cf(k,1)*arg)/(cf(k,2)+arg)
                     WRITE (22,99006) eng(l) , qc , qui(k,l) , 
     &                                100.*(qc-qui(k,l))/qui(k,l)
99006                FORMAT (8X,1F4.2,6X,1F9.4,5X,1F9.4,3X,1E10.2)
                  ENDDO
               ENDDO
            ENDDO
            GOTO 100 ! End of OP,GDET - back to input loop

C        Treat OP,RAND (randomise matrix elements)
         ELSEIF ( op2.EQ.'RAND' ) THEN
            READ (JZB,*) SE ! Seed for random number generator
            CALL MIXUP
            WRITE (22,99007)
99007       FORMAT (1X///5X,'MATRIX ELEMENTS RANDOMIZED...'///)
            CALL PRELM(2)
            GOTO 100 ! End of OP,RAND - back to input loop

C        Treat OP,TROU (troubleshooting)
         ELSEIF ( op2.EQ.'TROU' ) THEN
            ITS = 1 ! Create tape 18 flag
            READ (JZB,*) kmat , rlr
            GOTO 100 ! End of OP,TROU - back to input loop

C        Treat OP,REST (restart)
         ELSEIF ( op2.EQ.'REST' ) THEN
            irix = 12
            REWIND irix
            memax1 = MEMAX + 1
            DO lkj = 1 , MEMAX
               READ (irix,*) ELM(lkj)
            ENDDO
            DO lkj = 1 , memax1
               READ (JZB,*) lkj1 , xlk
               IF ( lkj1.EQ.0 ) GOTO 120
               ELM(lkj1) = xlk
            ENDDO
 120        WRITE (22,99008)
99008       FORMAT (1X///5X,'*****',2X,
     &              'RESTART-MATRIX ELEMENTS OVERWRITTEN',2X,'*****'///)
            DO kk = 1 , MEMAX
               la = mlt(kk)
               IF ( ivari(kk).GE.10000 ) THEN
                  kk1 = ivari(kk)/10000
                  kk2 = ivari(kk) - 10000*kk1
                  la1 = la
                  IF ( kk2.GE.100 ) THEN
                     la1 = kk2/100
                     kk2 = kk2 - 100*la1
                  ENDIF
                  inx1 = MEM(kk1,kk2,la1)
C      ELML(KK)=ELML(INX1)*ELM(KK)/ELM(INX1)
C      ELMU(KK)=ELMU(INX1)*ELM(KK)/ELM(INX1)
                  SA(kk) = ELM(kk)/ELM(inx1)
                  IVAR(kk) = 1000 + inx1
                  IF ( ELMU(kk).LE.ELML(kk) ) THEN
                     elmi = ELMU(kk)
                     ELMU(kk) = ELML(kk)
                     ELML(kk) = elmi
                  ENDIF
               ENDIF
            ENDDO
            CALL PRELM(2) ! Parameter is 4 in gosia2
            GOTO 100 ! End of OP,REST - back to input loop

C     Treat OP,SELE
         ELSEIF ( op2.EQ.'SELE' ) THEN
            CALL SELECT
            GOTO 2000 ! End of execution

C     Treat OP,BRIC
         ELSEIF ( op2.EQ.'BRIC' ) THEN
            CALL BRICC
            GOTO 100 ! End of OP,BRIC - back to input loop

C        Treat other options
         ELSE

C           Treat OP,RE,A (release A)
            IF ( op2.EQ.'RE,A' ) GOTO 900
           
C           Treat OP,RE,F (release F)
            IF ( op2.EQ.'RE,F' ) GOTO 900

C           Treat OP,ERRO (calculate errors)
            IF ( op2.EQ.'ERRO' ) THEN
               READ (JZB,*) idf , ms , mend , irep , ifc , remax
               rem = LOG(remax)
               LOCKS = 0
               LOCKF = 0
               JENTR = 1 ! Flag to indicate we are in OP,ERRO
               sh = 1.
               ifbp = 0
               inpo = 1
               inko = 1
               IF ( iosr.NE.0 .AND. idf.NE.0 ) THEN
                  inn = 0
                  ij = MULTI(1)
                  IF ( ij.NE.0 ) THEN
                     DO ij = 1 , NMAX
                        lxd = LDNUM(1,ij)
                        IF ( lxd.NE.0 ) THEN
                           DO ijk = 1 , lxd
                              inn = inn + 1
                           ENDDO
                        ENDIF
                     ENDDO
                     inpo = inn + 1
                  ENDIF
                  DO ij = 1 , NMAX
                     lxd = LDNUM(2,ij)
                     IF ( lxd.NE.0 ) THEN
                        DO ijk = 1 , lxd
                           inn = inn + 1
                        ENDDO
                     ENDIF
                  ENDDO
                  inko = inn
                  IF ( irep.NE.2 ) THEN
                     WRITE (IUNIT3,*) NMAX , MEMAX , inpo , inko
                     DO inn = 1 , NMAX
                        WRITE (IUNIT3,*) inn , SPIN(inn) , EN(inn)
                     ENDDO
                     DO inn = 1 , MEMAX
                        WRITE (IUNIT3,*) inn , LEAD(1,inn) , LEAD(2,inn)
                     ENDDO
                     DO inn = 1 , MEMAX
                        WRITE (IUNIT3,*) inn , ELM(inn)
                     ENDDO
                  ENDIF ! IF ( irep.NE.2 )
               ENDIF ! IF ( iosr.NE.0 .AND. idf.NE.0 )
               IF ( irep.NE.0 ) THEN
                  REWIND 15
                  READ (15,*) (DEVD(kh1),DEVU(kh1),kh1=1,MEMAX)
               ELSE
                  DO kh1 = 1 , MEMAX
                     DEVD(kh1) = ELML(kh1) - ELM(kh1)
                     DEVU(kh1) = ELMU(kh1) - ELM(kh1)
                  ENDDO
               ENDIF
               IF ( IMIN.EQ.0 ) CALL CMLAB(0,dsig,ttttt)
               IF ( ERR ) GOTO 2000 ! Error
               IF ( IMIN.NE.0 ) GOTO 400
               GOTO 1300 ! End of OP,ERRO

C           Treat OP,RE,C (release C)
            ELSEIF ( op2.EQ.'RE,C' ) THEN
               jfre = 1
               irfix = 0
               GOTO 1000 ! End of OP,RE,C

C           Treat OP,TITL (title)
            ELSEIF ( op2.EQ.'TITL' ) THEN
               READ (JZB,99009) (title(k),k=1,20)
99009          FORMAT (20A4)
               WRITE (22,99010) (title(k),k=1,20)
99010          FORMAT (10X,20A4/10X,100('-'))
               GOTO 100 ! End of OP,TITL - back to input loop

            ELSE

C              Treat OP,GOSI
               IF ( op2.EQ.'GOSI' ) GOTO 200

C              Treat OP,COUL
               IF ( op2.EQ.'COUL' ) GOTO 200

C              Treat OP,EXIT
               IF ( op2.EQ.'EXIT' ) THEN
                  GOTO 430 ! End of OP,EXIT

C              Treat OP,MINI
               ELSEIF ( op2.EQ.'MINI' ) THEN
                  READ (JZB,*) imode , nptl , chiok , conu , xtest , 
     &                 LOCKF , NLOCK , IFBFL , LOCKS , DLOCK
                  op2 = opcja
                  IMIN = IMIN + 1
                  IF ( IMIN.NE.1 ) GOTO 1400
                  GOTO 1200 ! End of OP,MINI

C              Treat OP,THEO
               ELSEIF ( op2.EQ.'THEO' ) THEN
                  irix = 12
                  REWIND (irix)
                  ibaf = 1
                  DO jb = 1 , LP1 ! LP1 = 50 (maximum number of experiments)
                     DO lb = 1 , 2
                        xlevb(jb,lb) = 0
                     ENDDO
                  ENDDO
                  READ (JZB,*) nbands ! Number of bands
                  IF ( nbands.LE.0 ) ibaf = 0
                  nbands = ABS(nbands)
                  DO nl = 1 , 8
                     DO jb = 1 , nbands
                        DO jl = 1 , nbands
                           DO kl = 1 , 3
                              bm(nl,jb,jl,kl) = 0.
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
                  DO jb = 1 , nbands
                     READ (JZB,*) bk , ilevls ! K of band, number of levels in band
                     READ (JZB,*) (levl(ib),ib=1,ilevls) ! Level list for band
                     DO kb = 1 , ilevls
                        inva = levl(kb)
                        xlevb(inva,2) = bk
                        xlevb(inva,1) = DBLE(jb)
                     ENDDO
                  ENDDO
                  DO nl = 1 , 8
                     READ (JZB,*) nnl ! Multipolarity
 126                 IF ( nnl.LE.0 ) GOTO 130
                     READ (JZB,*) jb1 , jb2 ! band indices
                     IF ( jb1.NE.0 ) THEN
                        READ (JZB,*) (bm(nnl,jb1,jb2,j),j=1,3) ! intrinsic moments
                        DO j = 1 , 3
                           bm(nnl,jb2,jb1,j) = bm(nnl,jb1,jb2,j)
                        ENDDO
                        GOTO 126
                     ENDIF
                  ENDDO
 130              DO kb = 1 , MEMAX
                     IF ( ibaf.NE.0 ) THEN
                        ind1 = LEAD(1,kb)
                        ind2 = LEAD(2,kb)
                        xi1 = SPIN(ind1)
                        xi2 = SPIN(ind2)
                        lamd = mlt(kb)
                        nb1 = INT(xlevb(ind1,1)+.1)
                        nb2 = INT(xlevb(ind2,1)+.1)
                        xk1 = xlevb(ind1,2)
                        xk2 = xlevb(ind2,2)
                        xm1 = bm(lamd,nb1,nb2,1)
                        xm2 = bm(lamd,nb1,nb2,2)
                        xm3 = bm(lamd,nb1,nb2,3)
                        ELM(kb) = ELMT(xi1,xi2,lamd,nb1,nb2,xk1,xk2,xm1,
     &                            xm2,xm3)
                        IF ( ABS(ELM(kb)).LT.1E-6 ) ELM(kb) = 1.E-6
                        irix = 12
                        WRITE (irix,*) ELM(kb)
                     ENDIF
                  ENDDO
                  GOTO 100 ! End of OP,THEO - back to input loop

C              Treat OP,YIEL
               ELSEIF ( op2.EQ.'YIEL' ) THEN
                  CALL ADHOC(oph,idr,nfd,ntap,iyr)
                  GOTO 100 ! End of OP,YIEL - back to input loop

C              Treat OP,INTG
               ELSEIF ( op2.EQ.'INTG' ) THEN
                  REWIND 14
                  lfagg = 1
                  IF ( SPIN(1).LT..25 ) ISO = 0
                  DO lx = 1 , NEXPT ! For each experiment
                     lpin = 1
                     IF ( ipinf.NE.0 ) THEN
                        IF ( jpin(lx).NE.0 ) lpin = jpin(lx)
                     ENDIF
                     IEXP = lx
                     tth = TLBDG(lx)
                     enh = EP(lx)
                     DO mpin = 1 , lpin ! For each pin diode
                        IF ( iecd(lx).EQ.1 ) THEN ! Circular detector
                           READ (JZB,*) ne , ntt , emn , emx , wth , 
     &                          wph , wthh
                           mfla = 1
                           CALL COORD(wth,wph,wthh,ntt,0,pfi,wpi,tth,lx,
     &                                tmn,tmx)
                        ELSE
                           READ (JZB,*) ne , ntt , emn , emx , tmn , tmx
                           mfla = 0
                           IF ( ntt.LT.0 ) mfla = 1
                        ENDIF
                        ntt = ABS(ntt)
                        jan = NANG(lx)
                        jan1 = NDST(lx)
                        IF ( IRAWEX(lx).EQ.0 ) jan1 = jan
                        IF ( iecd(lx).EQ.1 ) THEN ! Circular detector
                           WRITE (14,*) ne , ntt , emn , emx , tmn , 
     &                                  tmx , jan1 , wth , wph , wthh
                        ELSE
                           WRITE (14,*) ne , ntt , emn , emx , tmn , 
     &                                  tmx , jan1 , tmx , tmx , tmx
                        ENDIF
                        READ (JZB,*) (XV(i),i=1,ne)
                        IF ( iecd(lx).NE.1 ) READ (JZB,*)
     &                       (YV(i),i=1,ntt)
                        IF ( tth.LT.0. ) ELMH(2*lx-1) = YV(1)
                        IF ( tth.LT.0. ) ELMH(2*lx) = YV(ntt)
                        DO kloop = 1 , ne ! For each energy meshpoint
                           enb = XV(kloop)
                           EP(lx) = enb
                           DO ktt = 1 , ntt
                              tta = SIGN(YV(ktt),tth)
                              IF ( IAXS(lx).NE.0 ) THEN ! If not axial symmetry
                                 IF ( iecd(lx).NE.1 ) THEN
                                    IF ( kloop.EQ.1 ) THEN
                                       READ (JZB,*) nfi ! Number of phi ranges
                                       READ (JZB,*) 
     &                                    (fiex1(ktt,jfi,1),fiex1(ktt,
     &                                    jfi,2),jfi=1,nfi)
                                       IF ( tth.LT.0. ) THEN
                                         DO jfi = 1 , nfi ! For each phi angle
                                         fiex1(ktt,jfi,1)
     &                                      = fiex1(ktt,jfi,1) + 180.
                                         fiex1(ktt,jfi,2)
     &                                      = fiex1(ktt,jfi,2) + 180.
                                         ENDDO
                                       ENDIF
                                    ENDIF
                                 ENDIF
                              ENDIF ! If not axial symmetry
                              TLBDG(lx) = tta
                              IF ( kloop.EQ.1 ) THEN
                                 IF ( iecd(lx).NE.0 ) THEN
                                    nfi = 1
                                    fiex1(ktt,1,1) = wpi(ktt,1) ! Lower phi limit
                                    fiex1(ktt,1,2) = wpi(ktt,2) ! Upper phi limit
                                 ENDIF
                              ENDIF
                              CALL CMLAB(lx,dsig,tetrc)
                              IF ( ERR ) GOTO 2000 ! Error
                              tting = TLBDG(lx)
                              IF ( ERR ) GOTO 1900 ! Troubleshoot
                              CALL LOAD(lx,1,1,0.D0,jj)
                              CALL ALLOC(ACCUR)
                              CALL SNAKE(lx,ZPOL)
                              CALL SETIN
                              DO j = 1 , LMAX ! For each spin up to ground-state spin + 1
                                 polm = DBLE(j-1) - SPIN(1)
                                 CALL LOAD(lx,2,1,polm,jj)
                                 CALL STING(jj)
                                 CALL PATH(jj)
                                 CALL INTG(IEXP)
                                 CALL TENB(j,bten,LMAX)
                              ENDDO
                              CALL TENS(bten)
                              CALL DECAY(ccd,0,ccc)
                              DO j = 1 , LP2 ! LP2 = 1500 (maximum number of matrix elements)
                                 DO ijan = 1 , 20
                                    SUMCL(ijan,j) = 0.
                                 ENDDO
                              ENDDO
                              ija0 = 0
                              DO ijan = 1 , jan ! For each detector angle
                                 IF ( IAXS(lx).EQ.0 ) nfi = 1
                                 DO jyi = 1 , idr
                                    GRAD(jyi) = 0.
                                 ENDDO
                                 todfi = 0.
                                 DO jfi = 1 , nfi ! For each phi angle
                                    fi0 = fiex1(ktt,jfi,1)/57.2957795
                                    fi1 = fiex1(ktt,jfi,2)/57.2957795
                                    gth = AGELI(IEXP,ijan,1)
                                    fm = (fi0+fi1)/2.
                                    figl = AGELI(IEXP,ijan,2)
                                    CALL ANGULA(YGN,idr,1,fi0,fi1,tetrc,
     &                                 gth,figl,ijan,op2)
                                    IF ( IFMO.NE.0 ) THEN ! If correction due to recoil
                                       id = ITMA(IEXP,ijan) ! Get detector identity
                                       d = ODL(id) ! Get result of OP,GDET calculation
                                       rx = d*SIN(gth)*COS(figl-fm)
     &                                    - .25*SIN(tetrc)*COS(fm)
                                       ry = d*SIN(gth)*SIN(figl-fm)
     &                                    - .25*SIN(tetrc)*SIN(fm)
                                       rz = d*COS(gth) - .25*COS(tetrc)
                                       rl = SQRT(rx*rx+ry*ry+rz*rz)
                                       sf = d*d/rl/rl
                                       thc = TACOS(rz/rl)
                                       fic = ATAN2(ry,rx)
                                       CALL ANGULA(YGP,idr,1,fi0,fi1,
     &                                    tetrc,thc,fic,ijan,op2)
                                       DO ixl = 1 , idr ! For each decay
                                         ixm = KSEQ(ixl,3)
                                         tfac = TAU(ixm)
                                         YGN(ixl) = YGN(ixl)
     &                                      + .01199182*tfac*BETAR(IEXP)
     &                                      *(sf*YGP(ixl)-YGN(ixl))
                                       ENDDO ! Loop on decays
                                    ENDIF ! If correction due to recoil
                                    IF ( IRAWEX(lx).NE.0 ) THEN
                                       ipd = ITMA(lx,ijan) ! Get identity of detector
                                       DO jyi = 1 , idr ! For each decay
                                         ni = KSEQ(jyi,3)
                                         nf = KSEQ(jyi,4)
                                         decen = EN(ni) - EN(nf)
                                         cocos = SIN(tetrc)*SIN(gth)
     &                                      *COS(fm-figl) + COS(tetrc)
     &                                      *COS(gth)
                                         decen = decen*(1.+BETAR(lx)
     &                                      *cocos)
                                         CALL EFFIX(ipd,decen,effi)
                                         YGN(jyi) = YGN(jyi)*effi
                                       ENDDO
                                       inclus = ICLUST(lx,ijan) ! Cluster number for detector ijan
                                       IF ( inclus.NE.0 ) THEN
                                         DO jyi = 1 , idr ! For each decay
                                         SUMCL(inclus,jyi)
     &                                      = SUMCL(inclus,jyi)
     &                                      + YGN(jyi)
                                         ENDDO
                                         IF ( ijan.NE.LASTCL(lx,inclus)
     &                                      ) GOTO 132 ! If it is not the last detector in the cluster
                                         DO jyi = 1 , idr ! For each decay
                                         YGN(jyi) = SUMCL(inclus,jyi)
                                         ENDDO
                                       ENDIF
                                    ENDIF
                                    IF ( jfi.EQ.1 ) ija0 = ija0 + 1
                                    DO jyi = 1 , idr ! For each decay
                                       GRAD(jyi) = GRAD(jyi) + YGN(jyi)
                                    ENDDO ! Loop on decays jyi
                                    todfi = todfi + ABS(fi1-fi0)
                                 ENDDO ! For each phi angle jfi
                                 IF ( IAXS(lx).EQ.0 ) todfi = 6.283185
                                 ax = 1.
                                 IF ( mfla.EQ.1 ) ax = 1./todfi
                                 dsx = dsig
                                 IF ( mfla.NE.1 ) dsx = dsig*todfi
                                 dsxm(mpin,kloop,ktt) = dsx
                                 WRITE (17,*) lx , mpin , kloop , ktt , 
     &                                  dsx
                                 WRITE (14,*) lx , enb , tting , ija0 , 
     &                                  dsx , 
     &                                  (GRAD(jyi)*dsig*ax,jyi=1,idr)
                                 IF ( IPRM(11).EQ.1 ) THEN
                                    WRITE (22,99048) lx , ija0 , enb , 
     &                                 tta
                                    IF ( tta.LT.0. ) WRITE (22,99017)
     &                                 tting
99017                               FORMAT (5X,
     &                             'RESPECTIVE TARGET SCATTERING ANGLE='
     &                             ,1F7.3,1X,'DEG'/)
                                    DO jyi = 1 , idr
                                       ni = KSEQ(jyi,3)
                                       nf = KSEQ(jyi,4)
                                       WRITE (22,99049) ni , nf , 
     &                                    SPIN(ni) , SPIN(nf) , 
     &                                    GRAD(jyi)*dsig*ax , GRAD(jyi)
     &                                    /GRAD(IDRN)
                                    ENDDO ! Loop on decays jyi
                                 ENDIF ! If printout of yields at meshpoints
 132                             CONTINUE
                              ENDDO ! Loop on detector angles ijan
                           ENDDO ! Loop on theta angles ktt
                        ENDDO ! Loop on energy meshpoints kloop
                     ENDDO ! Loop on pin diodes mpin
                      
                     EP(lx) = enh
                     TLBDG(lx) = tth
                  ENDDO ! Loop on experiments lx
                  REWIND 14
                  REWIND 15
                  iske = 0
                  DO na = 1 , LP6 ! LP6 = 32 (maximum number of gamma detectors)
                     ILE(na) = 1
                  ENDDO
                  ilx = 0
C                 We have now performed the full coulex calculation at each of the
C                 meshpoints, so now we start the integration
                  DO lx = 1 , NEXPT ! Loop over experiments
C                    Read tape 17
                     REWIND 17
                     DO ijaja = 1 , 300000
                        READ (17,*,END=134) jjlx , jmpin , jkloo , 
     &                        jktt , dsx
                        IF ( jjlx.EQ.lx ) dsxm(jmpin,jkloo,jktt) = dsx
                     ENDDO
 134                 na = NANG(lx)
                     IF ( lx.NE.1 ) THEN
                        DO na1 = 1 , LP6 ! LP6 = 32 (maximum number of gamma detectors)
                           ILE(na1) = ILE(na1) + NYLDE(lx-1,na1)
                        ENDDO
                     ENDIF
                     READ (JZB,*) nptx ! Number of meshpoints for stopping powers
                     IF ( nptx.NE.0 ) THEN
                        READ (JZB,*) (esp(i),i=1,nptx) ! Energy
                        READ (JZB,*) (dedx(i),i=1,nptx) ! Stopping power
                        npt = nptx
                     ENDIF
                     READ (JZB,*) npce , npct
                     mfla = 0
                     IF ( npct.LT.0 ) mfla = 1
                     IF ( iecd(lx).EQ.1 ) mfla = 1
                     npct = ABS(npct)
                     IF ( npct.GT.100 )
     &                  STOP 'ABS(NI2) is limited to 100!'
                     npce = npce + MOD(npce,2)
                     npct = npct + MOD(npct,2)
                     mpin = 1
                     IF ( ipinf.NE.0 ) THEN
                        IF ( jpin(lx).NE.0 ) mpin = jpin(lx)
                     ENDIF
                     dst = 0.
                     DO lpin = 1 , mpin ! Loop over pin diodes
                        ilx = ilx + 1
                        IF ( ilx.NE.1 )
     &                       CALL TAPMA(lx,iske,isko,iskf,nflr,idr,0,
     &                       nft,enb)
                        READ (14,*) ne , ntt , emn , emx , tmn , tmx , 
     &                              jan , wth , wph , wthh
                        iocc = (ne+ntt)*idr
                        IF ( iocc.GT.izcap ) GOTO 1800
                        hen = (emx-emn)/npce
                        npce1 = npce + 1
                        het = (tmx-tmn)/npct ! Step in theta in degrees
                        npct1 = npct + 1
                        IF ( iecd(lx).EQ.1 ) ! Circular detector
     &                       CALL COORD(wth,wph,wthh,npct1,1,pfi,wpi,
     &                       TLBDG(lx),lx,tmn,tmx)
                        IF ( iecd(lx).NE.1 ) THEN
                           IF ( mfla.EQ.1 ) READ (JZB,*)
     &                          (pfi(j),j=1,npct1)
                        ENDIF
                        het = het/57.2957795 ! Step in theta in radians
                        
C                       Interpolate stopping power for each of the energies
C                       that we need. esp is an array of energies and dedx is
C                       an array containing the stopping powers at those
C                       energies. Function is unweighted sqrt. The energies
C                       are not the energies we gave for the meshpoints, but
C                       the range over which we integrate the bombarding energy
C                       with the number of steps specified.
                        DO j = 1 , npce1
                           xx = (j-1)*hen + emn
                           IF ( ISPL.EQ.0 )
     &                        CALL LAGRAN(esp,dedx,npt,1,xx,yy,3,1)
                           IF ( ISPL.EQ.1 )
     &                        CALL SPLNER(esp,dedx,npt,xx,yy,3)
                           HLMLM(j) = 1./yy
                        ENDDO
                         
C                       Now we calculate for all the mesh points. 
                        naa = NDST(lx)
                        IF ( IRAWEX(lx).EQ.0 ) naa = NANG(lx)
                        iskf = naa - 1
                        DO ja = 1 , naa ! Loop over detector angles
                           icll = 3 ! Weighting mode
                           DO je = 1 , ne ! ne = number of energy mesh points
                              lu = ILE(ja)
                              isko = (je-1)*naa*ntt + ja - 1
                              CALL TAPMA(lx,iske,isko,iskf,ntt,idr,1,
     &                           nft,enb)
                              IF ( nft.EQ.1 ) GOTO 1900 ! Troubleshoot
                              DO jd = 1 , idr ! For each decay
                                 DO jtp = 1 , ntt ! ntt = number of theta meshpoints
                                    IF ( jd.EQ.1 .AND. ja.EQ.1 )
     &                                 DSG(jtp) = dsxm(lpin,je,jtp)
                                    jyv = (jtp-1)*idr + jd
                                    YV(jtp) = ZETA(jyv) ! Point yield
                                 ENDDO ! Loop on theta meshpoints jtp
                                 DO jt = 1 , npct1 ! number of equal divisions in theta for interpolation
                                    xx = (jt-1)*het + tmn/57.2957795
                                    IF ( ISPL.EQ.0 )
     &                                 CALL LAGRAN(XV,YV,ntt,jt,xx,yy,2,
     &                                 icll) ! interpolate point yield at theta = xx
                                    IF ( ISPL.EQ.1 )
     &                                 CALL SPLNER(XV,YV,ntt,xx,yy,2) ! interpolate point yield at theta = xx
                                    IF ( ISPL.EQ.0 )
     &                                 CALL LAGRAN(XV,DSG,ntt,jt,xx,zz,
     &                                 2,icll) ! interpolate gamma yield at theta = xx
                                    IF ( ISPL.EQ.1 )
     &                                 CALL SPLNER(XV,DSG,ntt,xx,zz,
     &                                 2) ! interpolate gamma yield at theta = xx
                                    IF ( mfla.EQ.1 ) yy = yy*pfi(jt)
     &                                 /57.2957795
                                    IF ( yy.LE.0. ) yy = 1.E-15
                                    IF ( mfla.EQ.1 ) zz = zz*pfi(jt)
     &                                 /57.2957795
                                    XI(jt) = yy*SIN(xx) ! yy = integral of point yields over phi
                                    IF ( jd.EQ.1 .AND. ja.EQ.1 ) HLM(jt)
     &                                 = zz*SIN(xx) ! zz = integral over phi of Rutherford cross section
                                 ENDDO ! Loop on equal theta divisions jt
                                 icll = 4
                                 locat = ntt*idr + (je-1)*idr + jd
C                                Integrate point yields over theta using Simpson's rule
                                 ZETA(locat) = SIMIN(npct1,het,XI)
C                                If it is first decay and angle, integrate Rutherford cross section over theta
                                 IF ( jd.EQ.1 .AND. ja.EQ.1 ) DSE(je)
     &                                = SIMIN(npct1,het,HLM)
                                 ZV(je) = enb
                              ENDDO ! Loop on decays jd
                           ENDDO ! Loop on energy meshpoints je

C    Interpolation over energy:
C    The array ZV contains the energies of the meshpoints and the elements of the YV
C    array are set to the angle-integrated yield for each decay at the corresponding
C    energy, while DSE contains the Rutherford cross section for those energies. Since
C    the energies of the meshpoints are not necessarily equally spaced, we need to
C    interpolate to a set of equally spaced energies separated by "hen" starting from
C    "emn". To get the contribution from each energy, dE = 1 / (stopping power). Note
C    that we only evaluate the Rutherford cross section for the first decay and first
C    angle, since it is the same for all.

                           icll = 3
                           DO jd = 1 , idr ! For each decay
                              DO jtp = 1 , ne ! For each energy meshpoint
                                 jyv = (jtp-1)*idr + jd + ntt*idr
                                 YV(jtp) = ZETA(jyv)
                              ENDDO ! Loop on energy meshpoints jtp
                              DO jt = 1 , npce1 ! npce1 is number of equal energy steps
                                 xx = (jt-1)*hen + emn

C                                Interpolate the angle-integrated yield for this energy
                                 IF ( ISPL.EQ.0 )
     &                                CALL LAGRAN(ZV,YV,ne,jt,xx,yy,2,
     &                                icll)
                                 IF ( ISPL.EQ.1 )
     &                                CALL SPLNER(ZV,YV,ne,xx,yy,2)

C                                Interpolate Rutherford cross-section for this energy
                                 IF ( jd.EQ.1 .AND. ja.EQ.1 .AND. ! Only for first decay and angle
     &                                ISPL.EQ.0 )
     &                                CALL LAGRAN(ZV,DSE,ne,jt,xx,zz,2,
     &                                icll) ! Interpolate for this energy
                                 IF ( jd.EQ.1 .AND. ja.EQ.1 .AND.
     &                                ISPL.EQ.1 )
     &                                CALL SPLNER(ZV,DSE,ne,xx,zz,2) ! Interpolate for this energy
                                 IF ( jd.EQ.1 .AND. ja.EQ.1 ) HLM(jt)
     &                             = zz*HLMLM(jt) ! HLMLM = 1 / stopping power
                                 XI(jt) = yy*HLMLM(jt)
                              ENDDO ! Loop on equal energy steps

C   So now after this loop, we have XI containing the angle-integrated yield times dE for 
C   a set of equally spaced energies, so we use Simpson's rule to integrate them and store
C   in GRAD(jd). The first time, we also have in HLM a set of Rutherford cross-sections for
C   equally spaced energies, which we integrate in the same way.
                              icll = 4
                              IF ( jd.EQ.1 .AND. ja.EQ.1 )
     &                             DS = SIMIN(npce1,hen,HLM) ! integrate
                              GRAD(jd) = SIMIN(npce1,hen,XI)
                           ENDDO ! Loop over decays jd

                           IF ( ja.EQ.1 ) dst = dst + DS
                           IF ( ja.EQ.1 ) WRITE (22,99018) DS , lx
99018                      FORMAT (1X/////5X,
     &                            'INTEGRATED RUTHERFORD CROSS SECTION='
     &                            ,1E9.4,2X,'FOR EXP.',1I2///)

                           WRITE (22,99019) lx , ja , emn , emx , tmn , 
     &                            tmx
99019                      FORMAT (1X,//50X,'INTEGRATED YIELDS'//5X,
     &                             'EXPERIMENT ',1I2,2X,'DETECTOR ',
     &                             1I2/5X,'ENERGY RANGE ',1F8.3,'---',
     &                             1F8.3,1X,'MEV',3X,
     &                             'SCATTERING ANGLE RANGE ',1F7.3,
     &                             '---',1F7.3,1X,'DEG'//5X,'NI',5X,
     &                             'NF',5X,'II',5X,'IF',5X,'YIELD',5X,
     &                             'NORMALIZED YIELD'/)
                           DO jd = 1 , idr
                              WRITE (15,*) GRAD(jd)
                           ENDDO
                           DO jd = 1 , idr
                              ni = KSEQ(jd,3)
                              nf = KSEQ(jd,4)
                              WRITE (22,99049) ni , nf , SPIN(ni) , 
     &                               SPIN(nf) , GRAD(jd) , GRAD(jd)
     &                               /GRAD(IDRN) ! IDRN is the normalising transition
                           ENDDO
                        ENDDO ! Loop over detector angles ja

                        IF ( iecd(lx).EQ.1 ) THEN ! Circular detector
                           IF ( jpin(lx).EQ.0 ) THEN
                              CALL COORD(wth,wph,wthh,1,2,pfi,wpi,
     &                           TLBDG(lx),lx,txx,txx)
                              WRITE (22,99020) FIEX(lx,1)*57.2957795 , 
     &                               FIEX(lx,2)*57.2957795 , lx
99020                         FORMAT (//5X,
     &                          'WARNING: THE PHI ANGLE WAS REPLACED BY'
     &                          ,1X,F8.3,1X,'TO',F8.3,3X,
     &                          'FOR EXPERIMENT',2X,I3)
                              IF ( TLBDG(lx).LT.0 ) THEN
                                 FIEX(lx,1) = FIEX(lx,1) + 3.14159265
                                 FIEX(lx,2) = FIEX(lx,2) + 3.14159265
                              ENDIF ! If theta_lab < 0
                           ENDIF ! If no pin diodes
                        ENDIF ! If circular detector
                        iske = iske + ne*ntt*naa
                     ENDDO ! Loop over pin diodes
                     IF ( mpin.GT.1 ) WRITE (22,99021) dst , lx
99021                FORMAT (1x//2x,
     &                      'Total integrated Rutherford cross section='
     &                      ,1E8.3,' for exp. ',1I2/)
                  ENDDO
                  REWIND 17 ! Added PJN (17Jul2009)
                  IF ( ipinf.NE.0 ) THEN
                     ngpr = 0
                     DO lx = 1 , NEXPT ! For each experiment
                        nged = NDST(lx)
                        IF ( IRAWEX(lx).EQ.0 ) nged = NANG(lx)
                        IF ( lx.NE.1 ) ngpr = ngpr + idr*jpin(lx-1)
     &                       *NDST(lx-1)
                        lpin = jpin(lx)
                        IF ( lpin.EQ.0 ) lpin = 1
                        DO jgd = 1 , nged ! For each angle or dataset
                           DO jd = 1 , idr
                              GRAD(jd) = 0.
                           ENDDO
                           DO mpin = 1 , lpin ! For each pin diode
                              REWIND 15
                              ndum = ngpr + (jgd-1)*idr + (mpin-1)
     &                          *nged*idr ! Was jgd instead of nged (PJN 17Jul2009)
                              IF ( ndum.NE.0 ) THEN
                                 DO jd = 1 , ndum
                                    READ (15,*) xx
                                 ENDDO
                              ENDIF
                              DO jd = 1 , idr ! For each decay
                                 READ (15,*) xx
                                 GRAD(jd) = GRAD(jd) + xx
                              ENDDO ! Loop on decays jd
                           ENDDO ! Loop on pin diodes mpin
                           WRITE (17,*) (GRAD(jd),jd=1,idr)
                        ENDDO ! Loop on angle or dataset jgd
                     ENDDO ! Loop on experiment lx
                     REWIND 15
                     REWIND 17
                     DO lx = 1 , NEXPT ! For each experiment
                        nged = NDST(lx)
                        IF ( IRAWEX(lx).EQ.0 ) nged = NANG(lx)
                        DO ija0 = 1 , nged ! For each angle or dataset
                           READ (17,*) (GRAD(jdy),jdy=1,idr)
                           DO jd = 1 , idr ! For each decay
                              WRITE (15,*) GRAD(jd)
                           ENDDO ! Loop on decays jd
                        ENDDO ! Loop on angle or dataset ija0
                     ENDDO ! Loop on experiments lx
                  ENDIF
                  GOTO 100 ! End of OP,INTG - back to input loop

C              Treat OP,INTI
               ELSEIF ( op2.EQ.'INTI' ) THEN
                  DO lx = 1 , NEXPT ! For each experiment store original ISKIN
                     iskin_protect(lx) = ISKIN(lx)
                  ENDDO
                  REWIND 14
                  lfagg = 1
                  IF ( SPIN(1).LT..25 ) ISO = 0
                  DO lx = 1 , NEXPT ! For each experiment
                     lpin = 1
                     IF ( ipinf.NE.0 ) THEN
                        IF ( jpin(lx).NE.0 ) lpin = jpin(lx)
                     ENDIF
                     IEXP = lx
                     tth = TLBDG(lx)
                     enh = EP(lx)
                     DO mpin = 1 , lpin ! For each pin diode
                        IF ( iecd(lx).EQ.1 ) THEN ! Circular detector
                           READ (JZB,*) ne , ntt , emn , emx , wth , 
     &                          wph , wthh
                           mfla = 1
                           CALL COORD(wth,wph,wthh,ntt,0,pfi,wpi,tth,lx,
     &                                tmn,tmx)
                        ELSE
                           READ (JZB,*) ne , ntt , emn , emx , tmn , tmx
                           mfla = 0
                           IF ( ntt.LT.0 ) mfla = 1
                        ENDIF
                        ntt = ABS(ntt)
                        jan = NANG(lx)
                        jan1 = NDST(lx)
                        IF ( IRAWEX(lx).EQ.0 ) jan1 = jan
                        IF ( iecd(lx).EQ.1 ) THEN ! Circular detector
                           WRITE (14,*) ne , ntt , emn , emx , tmn , 
     &                                  tmx , jan1 , wth , wph , wthh
                        ELSE
                           WRITE (14,*) ne , ntt , emn , emx , tmn , 
     &                                  tmx , jan1 , tmx , tmx , tmx
                        ENDIF
                        READ (JZB,*) (XV(i),i=1,ne)
                        IF ( iecd(lx).NE.1 ) READ (JZB,*)
     &                       (YV(i),i=1,ntt)
                        IF ( tth.LT.0. ) ELMH(2*lx-1) = YV(1)
                        IF ( tth.LT.0. ) ELMH(2*lx) = YV(ntt)
                        DO kloop = 1 , ne ! For each energy meshpoint
                           enb = XV(kloop)
                           EP(lx) = enb
                           DO ktt = 1 , ntt
                              tta = YV(ktt)
                              IF ( tth.LT.0 )
     &                           CALL INVKIN(EP(lx),EN(NCM),IZ1(lx),
     &                                       XA,XA1(lx),YV(ktt),tta,
     &                                       1,ISKIN(lx))
                              tta = SIGN(tta, tth)
                              IF ( IAXS(lx).NE.0 ) THEN ! If not axial symmetry
                                 IF ( iecd(lx).NE.1 ) THEN
                                    IF ( kloop.EQ.1 ) THEN
                                       READ (JZB,*) nfi ! Number of phi ranges
                                       READ (JZB,*) 
     &                                    (fiex1(ktt,jfi,1),fiex1(ktt,
     &                                    jfi,2),jfi=1,nfi)
                                       IF ( tth.LT.0. ) THEN
                                         DO jfi = 1 , nfi ! For each phi angle
                                         fiex1(ktt,jfi,1)
     &                                      = fiex1(ktt,jfi,1) + 180.
                                         fiex1(ktt,jfi,2)
     &                                      = fiex1(ktt,jfi,2) + 180.
                                         ENDDO
                                       ENDIF
                                    ENDIF
                                 ENDIF
                              ENDIF ! If not axial symmetry
                              TLBDG(lx) = tta
                              IF ( kloop.EQ.1 ) THEN
                                 IF ( iecd(lx).NE.0 ) THEN
                                    nfi = 1
                                    fiex1(ktt,1,1) = wpi(ktt,1) ! Lower phi limit
                                    fiex1(ktt,1,2) = wpi(ktt,2) ! Upper phi limit
                                 ENDIF
                              ENDIF
                              CALL CMLAB(lx,dsig,tetrc)
                              IF ( ERR ) GOTO 2000 ! Error
                              tting = TLBDG(lx)
                              IF ( ERR ) GOTO 1900 ! Troubleshoot
                              CALL LOAD(lx,1,1,0.D0,jj)
                              CALL ALLOC(ACCUR)
                              CALL SNAKE(lx,ZPOL)
                              CALL SETIN
                              DO j = 1 , LMAX ! For each spin up to ground-state spin + 1
                                 polm = DBLE(j-1) - SPIN(1)
                                 CALL LOAD(lx,2,1,polm,jj)
                                 CALL STING(jj)
                                 CALL PATH(jj)
                                 CALL INTG(IEXP)
                                 CALL TENB(j,bten,LMAX)
                              ENDDO
                              CALL TENS(bten)
                              CALL DECAY(ccd,0,ccc)
                              DO j = 1 , LP2 ! LP2 = 1500 (maximum number of matrix elements)
                                 DO ijan = 1 , 20
                                    SUMCL(ijan,j) = 0.
                                 ENDDO
                              ENDDO
                              ija0 = 0
                              DO ijan = 1 , jan ! For each detector angle
                                 IF ( IAXS(lx).EQ.0 ) nfi = 1
                                 DO jyi = 1 , idr
                                    GRAD(jyi) = 0.
                                 ENDDO
                                 todfi = 0.
                                 DO jfi = 1 , nfi ! For each phi angle
                                    fi0 = fiex1(ktt,jfi,1)/57.2957795
                                    fi1 = fiex1(ktt,jfi,2)/57.2957795
                                    gth = AGELI(IEXP,ijan,1)
                                    fm = (fi0+fi1)/2.
                                    figl = AGELI(IEXP,ijan,2)
                                    CALL ANGULA(YGN,idr,1,fi0,fi1,tetrc,
     &                                 gth,figl,ijan,op2)
                                    IF ( IFMO.NE.0 ) THEN ! If correction due to recoil
                                       id = ITMA(IEXP,ijan) ! Get detector identity
                                       d = ODL(id) ! Get result of OP,GDET calculation
                                       rx = d*SIN(gth)*COS(figl-fm)
     &                                    - .25*SIN(tetrc)*COS(fm)
                                       ry = d*SIN(gth)*SIN(figl-fm)
     &                                    - .25*SIN(tetrc)*SIN(fm)
                                       rz = d*COS(gth) - .25*COS(tetrc)
                                       rl = SQRT(rx*rx+ry*ry+rz*rz)
                                       sf = d*d/rl/rl
                                       thc = TACOS(rz/rl)
                                       fic = ATAN2(ry,rx)
                                       CALL ANGULA(YGP,idr,1,fi0,fi1,
     &                                    tetrc,thc,fic,ijan,op2)
                                       DO ixl = 1 , idr ! For each decay
                                         ixm = KSEQ(ixl,3)
                                         tfac = TAU(ixm)
                                         YGN(ixl) = YGN(ixl)
     &                                      + .01199182*tfac*BETAR(IEXP)
     &                                      *(sf*YGP(ixl)-YGN(ixl))
                                       ENDDO ! Loop on decays
                                    ENDIF ! If correction due to recoil
                                    IF ( IRAWEX(lx).NE.0 ) THEN
                                       ipd = ITMA(lx,ijan) ! Get identity of detector
                                       DO jyi = 1 , idr ! For each decay
                                         ni = KSEQ(jyi,3)
                                         nf = KSEQ(jyi,4)
                                         decen = EN(ni) - EN(nf)
                                         cocos = SIN(tetrc)*SIN(gth)
     &                                      *COS(fm-figl) + COS(tetrc)
     &                                      *COS(gth)
                                         decen = decen*(1.+BETAR(lx)
     &                                      *cocos)
                                         CALL EFFIX(ipd,decen,effi)
                                         YGN(jyi) = YGN(jyi)*effi
                                       ENDDO
                                       inclus = ICLUST(lx,ijan) ! Cluster number for detector ijan
                                       IF ( inclus.NE.0 ) THEN
                                         DO jyi = 1 , idr ! For each decay
                                         SUMCL(inclus,jyi)
     &                                      = SUMCL(inclus,jyi)
     &                                      + YGN(jyi)
                                         ENDDO
                                         IF ( ijan.NE.LASTCL(lx,inclus)
     &                                      ) GOTO 432 ! If it is not the last detector in the cluster
                                         DO jyi = 1 , idr ! For each decay
                                         YGN(jyi) = SUMCL(inclus,jyi)
                                         ENDDO
                                       ENDIF
                                    ENDIF
                                    IF ( jfi.EQ.1 ) ija0 = ija0 + 1
                                    DO jyi = 1 , idr ! For each decay
                                       GRAD(jyi) = GRAD(jyi) + YGN(jyi)
                                    ENDDO ! Loop on decays jyi
                                    todfi = todfi + ABS(fi1-fi0)
                                 ENDDO ! For each phi angle jfi
                                 IF ( IAXS(lx).EQ.0 ) todfi = 6.283185
                                 ax = 1.
                                 IF ( mfla.EQ.1 ) ax = 1./todfi
                                 dsx = dsig
                                 IF ( mfla.NE.1 ) dsx = dsig*todfi
                                 dsxm(mpin,kloop,ktt) = dsx
                                 WRITE (17,*) lx , mpin , kloop , ktt , 
     &                                  dsx
                                 WRITE (14,*) lx , enb , tting , ija0 , 
     &                                  dsx , 
     &                                  (GRAD(jyi)*dsig*ax,jyi=1,idr)
                                 IF ( IPRM(11).EQ.1 ) THEN
                                    WRITE (22,99048) lx , ija0 , enb , 
     &                                 tta
                                    IF ( tta.LT.0. ) WRITE (22,99017)
     &                                 tting
                                    DO jyi = 1 , idr
                                       ni = KSEQ(jyi,3)
                                       nf = KSEQ(jyi,4)
                                       WRITE (22,99049) ni , nf , 
     &                                    SPIN(ni) , SPIN(nf) , 
     &                                    GRAD(jyi)*dsig*ax , GRAD(jyi)
     &                                    /GRAD(IDRN)
                                    ENDDO ! Loop on decays jyi
                                 ENDIF ! If printout of yields at meshpoints
 432                             CONTINUE
                              ENDDO ! Loop on detector angles ijan
                           ENDDO ! Loop on theta angles ktt
                        ENDDO ! Loop on energy meshpoints kloop
                     ENDDO ! Loop on pin diodes mpin
                      
                     EP(lx) = enh
                     TLBDG(lx) = tth
                  ENDDO ! Loop on experiments lx
                  REWIND 14
                  REWIND 15
                  iske = 0
                  DO na = 1 , LP6 ! LP6 = 32 (maximum number of gamma detectors)
                     ILE(na) = 1
                  ENDDO
                  ilx = 0
C                 We have now performed the full coulex calculation at each of the
C                 meshpoints, so now we start the integration
                  DO lx = 1 , NEXPT ! Loop over experiments
C                    Read tape 17
                     REWIND 17
                     DO ijaja = 1 , 300000
                        READ (17,*,END=434) jjlx , jmpin , jkloo , 
     &                        jktt , dsx
                        IF ( jjlx.EQ.lx ) dsxm(jmpin,jkloo,jktt) = dsx
                     ENDDO
 434                 na = NANG(lx)
                     IF ( lx.NE.1 ) THEN
                        DO na1 = 1 , LP6 ! LP6 = 32 (maximum number of gamma detectors)
                           ILE(na1) = ILE(na1) + NYLDE(lx-1,na1)
                        ENDDO
                     ENDIF
                     READ (JZB,*) nptx ! Number of meshpoints for stopping powers
                     IF ( nptx.NE.0 ) THEN
                        READ (JZB,*) (esp(i),i=1,nptx) ! Energy
                        READ (JZB,*) (dedx(i),i=1,nptx) ! Stopping power
                        npt = nptx
                     ENDIF
                     READ (JZB,*) npce , npct
                     mfla = 0
                     IF ( npct.LT.0 ) mfla = 1
                     IF ( iecd(lx).EQ.1 ) mfla = 1
                     npct = ABS(npct)
                     IF ( npct.GT.100 )
     &                  STOP 'ABS(NI2) is limited to 100!'
                     npce = npce + MOD(npce,2)
                     npct = npct + MOD(npct,2)
                     mpin = 1
                     IF ( ipinf.NE.0 ) THEN
                        IF ( jpin(lx).NE.0 ) mpin = jpin(lx)
                     ENDIF
                     dst = 0.
                     DO lpin = 1 , mpin ! Loop over pin diodes
                        ilx = ilx + 1
                        IF ( ilx.NE.1 )
     &                       CALL TAPMA(lx,iske,isko,iskf,nflr,idr,0,
     &                       nft,enb)
                        READ (14,*) ne , ntt , emn , emx , tmn , tmx , 
     &                              jan , wth , wph , wthh
                        iocc = (ne+ntt)*idr
                        IF ( iocc.GT.izcap ) GOTO 1800
                        hen = (emx-emn)/npce
                        npce1 = npce + 1
                        het = (tmx-tmn)/npct ! Step in theta in degrees
                        npct1 = npct + 1
                        IF ( iecd(lx).EQ.1 ) ! Circular detector
     &                       CALL COORD(wth,wph,wthh,npct1,1,pfi,wpi,
     &                       TLBDG(lx),lx,tmn,tmx)
                        IF ( iecd(lx).NE.1 ) THEN
                           IF ( mfla.EQ.1 ) READ (JZB,*)
     &                          (pfi(j),j=1,npct1)
                        ENDIF
                        het = het/57.2957795 ! Step in theta in radians
                        
C                       Interpolate stopping power for each of the energies
C                       that we need. esp is an array of energies and dedx is
C                       an array containing the stopping powers at those
C                       energies. Function is unweighted sqrt. The energies
C                       are not the energies we gave for the meshpoints, but
C                       the range over which we integrate the bombarding energy
C                       with the number of steps specified.
                        DO j = 1 , npce1
                           xx = (j-1)*hen + emn
                           IF ( ISPL.EQ.0 )
     &                        CALL LAGRAN(esp,dedx,npt,1,xx,yy,3,1)
                           IF ( ISPL.EQ.1 )
     &                        CALL SPLNER(esp,dedx,npt,xx,yy,3)
                           HLMLM(j) = 1./yy
                        ENDDO
                         
C                       Now we calculate for all the mesh points. 
                        naa = NDST(lx)
                        IF ( IRAWEX(lx).EQ.0 ) naa = NANG(lx)
                        iskf = naa - 1
                        DO ja = 1 , naa ! Loop over detector angles
                           icll = 3 ! Weighting mode
                           DO je = 1 , ne ! ne = number of energy mesh points
                              lu = ILE(ja)
                              isko = (je-1)*naa*ntt + ja - 1
                              CALL TAPMA(lx,iske,isko,iskf,ntt,idr,1,
     &                           nft,enb)
                              IF ( nft.EQ.1 ) GOTO 1900 ! Troubleshoot
                              DO jd = 1 , idr ! For each decay
                                 DO jtp = 1 , ntt ! ntt = number of theta meshpoints
                                    IF ( jd.EQ.1 .AND. ja.EQ.1 )
     &                                 DSG(jtp) = dsxm(lpin,je,jtp)
                                    jyv = (jtp-1)*idr + jd
                                    YV(jtp) = ZETA(jyv) ! Point yield
                                 ENDDO ! Loop on theta meshpoints jtp
                                 DO jt = 1 , npct1 ! number of equal divisions in theta for interpolation
                                    xx = (jt-1)*het + tmn/57.2957795
                                    IF ( ISPL.EQ.0 )
     &                                 CALL LAGRAN(XV,YV,ntt,jt,xx,yy,2,
     &                                 icll) ! interpolate point yield at theta = xx
                                    IF ( ISPL.EQ.1 )
     &                                 CALL SPLNER(XV,YV,ntt,xx,yy,2) ! interpolate point yield at theta = xx
                                    IF ( ISPL.EQ.0 )
     &                                 CALL LAGRAN(XV,DSG,ntt,jt,xx,zz,
     &                                 2,icll) ! interpolate gamma yield at theta = xx
                                    IF ( ISPL.EQ.1 )
     &                                 CALL SPLNER(XV,DSG,ntt,xx,zz,
     &                                 2) ! interpolate gamma yield at theta = xx
                                    IF ( mfla.EQ.1 ) yy = yy*pfi(jt)
     &                                 /57.2957795
                                    IF ( yy.LE.0. ) yy = 1.E-15
                                    IF ( mfla.EQ.1 ) zz = zz*pfi(jt)
     &                                 /57.2957795
                                    XI(jt) = yy*SIN(xx) ! yy = integral of point yields over phi
                                    IF ( jd.EQ.1 .AND. ja.EQ.1 ) HLM(jt)
     &                                 = zz*SIN(xx) ! zz = integral over phi of Rutherford cross section
                                 ENDDO ! Loop on equal theta divisions jt
                                 icll = 4
                                 locat = ntt*idr + (je-1)*idr + jd
C                                Integrate point yields over theta using Simpson's rule
                                 ZETA(locat) = SIMIN(npct1,het,XI)
C                                If it is first decay and angle, integrate Rutherford cross section over theta
                                 IF ( jd.EQ.1 .AND. ja.EQ.1 ) DSE(je)
     &                                = SIMIN(npct1,het,HLM)
                                 ZV(je) = enb
                              ENDDO ! Loop on decays jd
                           ENDDO ! Loop on energy meshpoints je

C    Interpolation over energy:
C    The array ZV contains the energies of the meshpoints and the elements of the YV
C    array are set to the angle-integrated yield for each decay at the corresponding
C    energy, while DSE contains the Rutherford cross section for those energies. Since
C    the energies of the meshpoints are not necessarily equally spaced, we need to
C    interpolate to a set of equally spaced energies separated by "hen" starting from
C    "emn". To get the contribution from each energy, dE = 1 / (stopping power). Note
C    that we only evaluate the Rutherford cross section for the first decay and first
C    angle, since it is the same for all.

                           icll = 3
                           DO jd = 1 , idr ! For each decay
                              DO jtp = 1 , ne ! For each energy meshpoint
                                 jyv = (jtp-1)*idr + jd + ntt*idr
                                 YV(jtp) = ZETA(jyv)
                              ENDDO ! Loop on energy meshpoints jtp
                              DO jt = 1 , npce1 ! npce1 is number of equal energy steps
                                 xx = (jt-1)*hen + emn

C                                Interpolate the angle-integrated yield for this energy
                                 IF ( ISPL.EQ.0 )
     &                                CALL LAGRAN(ZV,YV,ne,jt,xx,yy,2,
     &                                icll)
                                 IF ( ISPL.EQ.1 )
     &                                CALL SPLNER(ZV,YV,ne,xx,yy,2)

C                                Interpolate Rutherford cross-section for this energy
                                 IF ( jd.EQ.1 .AND. ja.EQ.1 .AND. ! Only for first decay and angle
     &                                ISPL.EQ.0 )
     &                                CALL LAGRAN(ZV,DSE,ne,jt,xx,zz,2,
     &                                icll) ! Interpolate for this energy
                                 IF ( jd.EQ.1 .AND. ja.EQ.1 .AND.
     &                                ISPL.EQ.1 )
     &                                CALL SPLNER(ZV,DSE,ne,xx,zz,2) ! Interpolate for this energy
                                 IF ( jd.EQ.1 .AND. ja.EQ.1 ) HLM(jt)
     &                             = zz*HLMLM(jt) ! HLMLM = 1 / stopping power
                                 XI(jt) = yy*HLMLM(jt)
                              ENDDO ! Loop on equal energy steps

C   So now after this loop, we have XI containing the angle-integrated yield times dE for 
C   a set of equally spaced energies, so we use Simpson's rule to integrate them and store
C   in GRAD(jd). The first time, we also have in HLM a set of Rutherford cross-sections for
C   equally spaced energies, which we integrate in the same way.
                              icll = 4
                              IF ( jd.EQ.1 .AND. ja.EQ.1 )
     &                             DS = SIMIN(npce1,hen,HLM) ! integrate
                              GRAD(jd) = SIMIN(npce1,hen,XI)
                           ENDDO ! Loop over decays jd

                           IF ( ja.EQ.1 ) dst = dst + DS
                           IF ( ja.EQ.1 ) WRITE (22,99018) DS , lx

                           WRITE (22,99019) lx , ja , emn , emx , tmn , 
     &                            tmx
                           DO jd = 1 , idr
                              WRITE (15,*) GRAD(jd)
                           ENDDO
                           DO jd = 1 , idr
                              ni = KSEQ(jd,3)
                              nf = KSEQ(jd,4)
                              WRITE (22,99049) ni , nf , SPIN(ni) , 
     &                               SPIN(nf) , GRAD(jd) , GRAD(jd)
     &                               /GRAD(IDRN) ! IDRN is the normalising transition
                           ENDDO
                        ENDDO ! Loop over detector angles ja

                        IF ( iecd(lx).EQ.1 ) THEN ! Circular detector
                           IF ( jpin(lx).EQ.0 ) THEN
                              CALL COORD(wth,wph,wthh,1,2,pfi,wpi,
     &                           TLBDG(lx),lx,txx,txx)
                              WRITE (22,99020) FIEX(lx,1)*57.2957795 , 
     &                               FIEX(lx,2)*57.2957795 , lx
                              IF ( TLBDG(lx).LT.0 ) THEN
                                 FIEX(lx,1) = FIEX(lx,1) + 3.14159265
                                 FIEX(lx,2) = FIEX(lx,2) + 3.14159265
                              ENDIF ! If theta_lab < 0
                           ENDIF ! If no pin diodes
                        ENDIF ! If circular detector
                        iske = iske + ne*ntt*naa
                     ENDDO ! Loop over pin diodes
                     IF ( mpin.GT.1 ) WRITE (22,99021) dst , lx
                  ENDDO
                  REWIND 17 ! Added PJN (17Jul2009)
                  IF ( ipinf.NE.0 ) THEN
                     ngpr = 0
                     DO lx = 1 , NEXPT ! For each experiment
                        nged = NDST(lx)
                        IF ( IRAWEX(lx).EQ.0 ) nged = NANG(lx)
                        IF ( lx.NE.1 ) ngpr = ngpr + idr*jpin(lx-1)
     &                       *NDST(lx-1)
                        lpin = jpin(lx)
                        IF ( lpin.EQ.0 ) lpin = 1
                        DO jgd = 1 , nged ! For each angle or dataset
                           DO jd = 1 , idr
                              GRAD(jd) = 0.
                           ENDDO
                           DO mpin = 1 , lpin ! For each pin diode
                              REWIND 15
                              ndum = ngpr + (jgd-1)*idr + (mpin-1)
     &                          *nged*idr ! Was jgd instead of nged (PJN 17Jul2009)
                              IF ( ndum.NE.0 ) THEN
                                 DO jd = 1 , ndum
                                    READ (15,*) xx
                                 ENDDO
                              ENDIF
                              DO jd = 1 , idr ! For each decay
                                 READ (15,*) xx
                                 GRAD(jd) = GRAD(jd) + xx
                              ENDDO ! Loop on decays jd
                           ENDDO ! Loop on pin diodes mpin
                           WRITE (17,*) (GRAD(jd),jd=1,idr)
                        ENDDO ! Loop on angle or dataset jgd
                     ENDDO ! Loop on experiment lx
                     REWIND 15
                     REWIND 17
                     DO lx = 1 , NEXPT ! For each experiment
                        nged = NDST(lx)
                        IF ( IRAWEX(lx).EQ.0 ) nged = NANG(lx)
                        DO ija0 = 1 , nged ! For each angle or dataset
                           READ (17,*) (GRAD(jdy),jdy=1,idr)
                           DO jd = 1 , idr ! For each decay
                              WRITE (15,*) GRAD(jd)
                           ENDDO ! Loop on decays jd
                        ENDDO ! Loop on angle or dataset ija0
                     ENDDO ! Loop on experiments lx
                  ENDIF
                  DO lx = 1 , NEXPT ! For each experiment restore original ISKIN
                     ISKIN(lx) = iskin_protect(lx)
                  ENDDO
                  GOTO 100 ! End of OP,INTI - back to input loop

C              Treat OP,CORR
               ELSEIF ( op2.EQ.'CORR' ) THEN
                  CALL READY(idr,ntap,0)
                  REWIND 3
                  REWIND 15
                  REWIND 4
                  GOTO 1200 ! End of OP,CORR
               ELSE

C                 Treat OP,POIN
                  IF ( op2.EQ.'POIN' ) GOTO 1200

C                 Treat OP,MAP
                  IF ( op2.EQ.'MAP ' ) iobl = 1

C                 Treat OP,STAR
                  IF ( op2.EQ.'STAR' ) GOTO 1200

C                 Treat OP,SIXJ
                  IF ( op2.EQ.'SIXJ' ) THEN
                     DO k = 1 , 2
                        l = 4*k
                        DO j = 1 , 80
                           ixj = j - 1
                           DO ms = 1 , 5
                              mend = 2*(ms-3) + ixj
                              WRITE (14,*) WSIXJ(l,4,4,ixj,mend,ixj-4) ,
     &                               WSIXJ(l,4,4,ixj,mend,ixj-2) , 
     &                               WSIXJ(l,4,4,ixj,mend,ixj) , 
     &                               WSIXJ(l,4,4,ixj,mend,ixj+2) , 
     &                               WSIXJ(l,4,4,ixj,mend,ixj+4)
                           ENDDO
                        ENDDO
                     ENDDO
                     GOTO 2000 ! End of OP,SIXJ - normal end of execution

C                 Treat OP,RAW (raw uncorrected gamma yields)
                  ELSEIF ( op2.EQ.'RAW ' ) THEN
C                    Read absorber coefficients from unit 8
                     REWIND 8
                     DO l = 1 , 8
                        READ (8,*) (ABC(l,j),j=1,10) ! Absorption coefficients
                        DO j = 1 , 10
                           ABC(l,j) = LOG(ABC(l,j))
                        ENDDO
                     ENDDO
                     DO l = 1 , nfd
                        READ (8,*) (THICK(l,j),j=1,7) ! thickness of absorbers
                     ENDDO
                     DO l = 1 , LP1 ! LP1 = 50 (maximum number of experiments)
                        DO j = 1 , 200
                           ICLUST(l,j) = 0
                        ENDDO
                        DO j = 1 , 20
                           LASTCL(l,j) = 0
                        ENDDO
                        IRAWEX(l) = 0
                     ENDDO

C                    Read input from standard input
                     DO l = 1 , LP1 ! LP1 = 50 (maximum number of experiments)
                        READ (JZB,*) mexl ! experiment number
                        IF ( mexl.EQ.0 ) GOTO 100 ! Back to input loop
                        IRAWEX(mexl) = 1
                        n = NANG(mexl)
                        DO j = 1 , n
                           jj = ITMA(mexl,j) ! Get identity of detector
                           READ (JZB,*) (AKAVKA(k,jj),k=1,8) ! efficiency curve parameters
                           AKAVKA(9,jj) = ideff(mexl)
                        ENDDO
                        READ (JZB,*) kclust ! number of clusters
                        IF ( kclust.NE.0 ) THEN
                           DO j = 1 , kclust
                              READ (JZB,*) numcl ! Number of detectors for this cluster
                              READ (JZB,*) (liscl(k),k=1,numcl) ! Indices of logical detectors
                              LASTCL(l,j) = liscl(numcl) ! Index of last detector in cluster
                              DO k = 1 , numcl
                                 kk = liscl(k)
                                 ICLUST(l,kk) = j ! Set cluster number
                              ENDDO
                           ENDDO
                        ENDIF
                     ENDDO
                     GOTO 100 ! End of OP,RAW - back to input loop

C                 Treat OP,MAP
                  ELSEIF ( op2.EQ.'MAP ' ) THEN
                     GOTO 1200 ! End of OP,MAP 
                  ENDIF ! IF ( op2.EQ.'SIXJ' )
               ENDIF
            ENDIF
         ENDIF
      ENDIF ! End of if (op1.eq."OP, ") if statement

      WRITE (22,99022) op1 , op2
99022 FORMAT (5X,'UNRECOGNIZED OPTION',1X,1A3,1A4)
      GOTO 2000 ! Normal end of execution

C     Treat suboptions of OP,COUL and OP,GOSI
 200  READ (JZB,99023) op1 ! Read the suboption
99023 FORMAT (1A4)
      IF ( op1.EQ.'    ' ) GOTO 100 ! Back to input loop

C     Treat suboption LEVE (levels)
      IF ( op1.EQ.'LEVE' ) THEN
         NMAX = 0
         IF ( ABS(IPRM(1)).EQ.1 ) WRITE (22,99024)
99024    FORMAT (1X/40X,'LEVELS',//5X,'INDEX',5X,'PARITY',9X,'SPIN',11X,
     &           'ENERGY(MEV)')
         ndima = NDIM + 1
         DO k = 1 , ndima
            READ (JZB,*) ipo1 , ipo2 , po2 , po1 ! level number, parity, spin, energy
            IF ( ipo1.EQ.0 ) GOTO 200
            IF ( ipo1.EQ.1 .AND. ABS(po2).LT.1.E-6 ) ISO = 0
            NMAX = NMAX + 1
            SPIN(ipo1) = po2
            IF ( k.EQ.1 ) iph = ipo2
            iprc = ipo2 - iph
            IF ( iprc.NE.0 ) iprc = 1
            IFAC(ipo1) = (-1)**(iprc-INT(po2-SPIN(1)))
            EN(ipo1) = po1
            prp = '+'
            IF ( ipo2.EQ.-1 ) prp = '-'
            IF ( ABS(IPRM(1)).EQ.1 ) WRITE (22,99025) ipo1 , prp , 
     &           SPIN(ipo1) , EN(ipo1)
99025       FORMAT (5X,1I3,11X,1A1,10X,1F4.1,8X,1F10.4)
         ENDDO

C     Treat suboption ME (matrix elements)
      ELSEIF ( op1.EQ.'ME  ' ) THEN
         DO k = 1 , nmemx
            IF ( op2.EQ.'GOSI' ) THEN
               READ (JZB,*) ipo1 , ipo2 , po1 , bl , bu ! lamda, 0, 0, 0, 0 OR ind1, ind2, me, lo, hi
               iopri = 2
               icg = 2
            ELSE
               iopri = 1
               READ (JZB,*) ipo1 , ipo2 , po1 ! lambda, 0, 0 OR ind1, ind2, me
            ENDIF
            IF ( ipo1.NE.0 ) THEN
               IF ( ipo2.EQ.0 ) THEN
                  IF ( ipo1.LE.la ) GOTO 1600 ! Error - wrong sequence of multipolarities
                  LAMMAX = LAMMAX + 1
                  LAMDA(LAMMAX) = ipo1
                  ipo3 = 0
                  IF ( indx.EQ.0 ) GOTO 220
               ELSE
                  MULTI(la) = MULTI(la) + 1
                  indx = indx + 1
                  IF ( ipo1.GT.ABS(ipo2) ) GOTO 1500 ! Error - M.E. does not belong to the upper triangle
                  IF ( ipo1.NE.ipo3 ) THEN
                     IF ( ipo1.LT.ipo3 ) GOTO 1700 ! Error - repeated appearance of the state
                     ipo3 = ipo1
                  ENDIF
                  ELM(indx) = po1
                  mlt(indx) = la
                  LEAD(1,indx) = ipo1
                  LEAD(2,indx) = ABS(ipo2)
                  LDNUM(la,ipo1) = LDNUM(la,ipo1) + 1
                  IF ( op2.EQ.'GOSI' ) THEN
                     IF ( ipo2.LT.0 ) THEN ! If negative, bl and bu are indices
                                           ! to which we fix this element
                        IVAR(indx) = 10000*INT(bl) + INT(bu)
                     ELSE                 ! Otherwise they are limits
                        ELMU(indx) = bu
                        ELML(indx) = bl
                        IF ( ABS(bl-bu).LT.1.E-6 ) THEN
                           IVAR(indx) = 0 ! Fixed
                        ELSE
                           IVAR(indx) = 2
                           IF ( la.GT.4 ) IVAR(indx) = 1
                        ENDIF
                     ENDIF
                     isip = ISEX(ipo1) + 1
                     ISEX(ABS(ipo2)) = MIN(isip,ISEX(ABS(ipo2)))
                  ENDIF
                  GOTO 250
               ENDIF
            ENDIF
            DO kk = 1 , indx
               IF ( ABS(ELM(kk)).LE.1.E-6 ) ELM(kk) = 1.E-6
               IF ( IVAR(kk).GE.10000 ) THEN ! Correlated
                  kk1 = IVAR(kk)/10000
                  kk2 = IVAR(kk) - 10000*kk1
                  la1 = la
                  IF ( kk2.GE.100 ) THEN
                     la1 = kk2/100
                     kk2 = kk2 - 100*la1
                  ENDIF
                  inx1 = MEM(kk1,kk2,la1)
                  ELML(kk) = ELML(inx1)*ELM(kk)/ELM(inx1)
                  ELMU(kk) = ELMU(inx1)*ELM(kk)/ELM(inx1)
                  SA(kk) = ELM(kk)/ELM(inx1)
                  ivari(kk) = IVAR(kk)
                  IVAR(kk) = 1000 + inx1
                  IF ( ELMU(kk).LE.ELML(kk) ) THEN
                     elmi = ELMU(kk)
                     ELMU(kk) = ELML(kk)
                     ELML(kk) = elmi
                  ENDIF
               ENDIF
            ENDDO
            IF ( ipo1.EQ.0 ) GOTO 300
 220        la = ipo1
            IF ( la.GT.LMAXE .AND. la.LE.6 ) LMAXE = la
 250        CONTINUE
         ENDDO
 300     MEMAX = indx
         IF ( la.GT.6 ) MAGEXC = 1 ! Flag that we need magnetic excitations
         memx4 = MULTI(1) + MULTI(2) + MULTI(3) + MULTI(4)
         MEMX6 = memx4 + MULTI(5) + MULTI(6)
         IF ( ABS(IPRM(1)).EQ.1 ) CALL PRELM(iopri)
         DO kh = 1 , NMAX
            IF ( ISEX(kh).EQ.1111 ) ISEX(kh) = 1
         ENDDO
         DO kh = 1 , MEMAX
            ivarh(kh) = IVAR(kh)
         ENDDO

C     Treat suboption CONT (control)
      ELSEIF ( op1.EQ.'CONT' ) THEN
 350     READ (JZB,99026) op1 , fipo1
99026    FORMAT (1A4,1F7.1)
         ipo1 = INT(fipo1)
         IF ( op1.EQ.'ACP,' ) ACCA = 10.**(-fipo1)
         IF ( op1.EQ.'SEL,' ) ITS = 2
         IF ( op1.EQ.'SMR,' ) iosr = 1
         IF ( op1.EQ.'SPL,' ) ISPL = ipo1
         IF ( op1.EQ.'EFF,' ) THEN
            DO jjx = 1 , ipo1
               READ (JZB,*) ipo2 , ijx
               ideff(ipo2) = ijx
            ENDDO
         ENDIF
         IF ( op1.EQ.'FMI,' ) ifm = 1
         IF ( op1.EQ.'TEN,' ) itno = 1
         IF ( op1.EQ.'NCM,' ) NCM = ipo1
         IF ( op1.EQ.'WRN,' ) SGW = fipo1
         IF ( op1.EQ.'INT,' ) THEN
            DO jjx = 1 , ipo1
               READ (JZB,*) ipo2 , ijx
               INTERV(ipo2) = ijx
            ENDDO
         ELSE
            IF ( op1.EQ.'VAC,' ) THEN
               DO jjx = 1 , 7
                  READ (JZB,*) ijx , val
                  IF ( ijx.EQ.0 ) GOTO 350
                  G(ijx) = val
               ENDDO
            ELSE
               IF ( op1.EQ.'DIP,' ) DIPOL = 0.001*fipo1
               IF ( op1.EQ.'ACC,' ) ACCUR = 10.**(-fipo1)
               IF ( op1.EQ.'PRT,' ) THEN
                  DO jjx = 1 , 20
                     READ (JZB,*) inm1 , inm2
                     IF ( inm1.EQ.0 ) GOTO 350
                     IPRM(inm1) = inm2
                  ENDDO
                  GOTO 350
               ELSEIF ( op1.NE.'FIX,' ) THEN
                  IF ( op1.EQ.'SKP,' ) THEN
                     DO jjx = 1 , ipo1
                        READ (JZB,*) ijx
                        JSKIP(ijx) = 0
                     ENDDO
                     GOTO 350
                  ELSE
                     IF ( op1.EQ.'CRF,' ) ICS = 1
                     IF ( op1.EQ.'LCK,' ) THEN
 352                    READ (JZB,*) lck1 , lck2
                        IF ( lck1.EQ.0 ) GOTO 350
                        DO jjx = lck1 , lck2
                           ivarh(jjx) = 0
                           IVAR(jjx) = 0
                        ENDDO
                        GOTO 352
                     ELSE
                        IF ( op1.EQ.'INR,' ) INNR = 1
                        IF ( op1.EQ.'CRD,' ) THEN
                           DO jjx = 1 , ipo1
                              READ (JZB,*) ipo2
                              iecd(ipo2) = 1
                           ENDDO
                           GOTO 350
                        ELSE
                           IF ( op1.EQ.'CCF,' ) IPS1 = ipo1
                           IF ( op1.EQ.'PIN,' ) THEN
                              ipine = ipo1
                              ipinf = 1
                              DO ipp = 1 , ipine
                                 READ (JZB,*) ig1 , ig2
                                 jpin(ig1) = ig2
                              ENDDO
                              GOTO 350
                           ELSE
                              IF ( op1.NE.'END,' ) GOTO 350
                              GOTO 200
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
            READ (JZB,*) nallow
            DO jjx = 1 , nallow
               READ (JZB,*) ijk
               IVAR(ijk) = -IVAR(ijk)
            ENDDO
            DO jjx = 1 , MEMAX
               IF ( IVAR(jjx).GE.0 ) THEN
                  IF ( IVAR(jjx).LE.999 ) IVAR(jjx) = 0
               ENDIF
            ENDDO
            DO jjx = 1 , MEMAX
               IF ( IVAR(jjx).LT.0 ) IVAR(jjx) = -IVAR(jjx)
               ivarh(jjx) = IVAR(jjx)
            ENDDO
         ENDIF
         GOTO 350 ! Back to beginning of CONT loop

C     Treat suboption EXPT
      ELSEIF ( op1.EQ.'EXPT' ) THEN
         READ (JZB,*) NEXPT , IZ , XA
         G(1) = 3.             ! AVJI
         G(2) = .02            ! GAMMA
         G(3) = .0345          ! XLAMB
         G(4) = 3.5            ! TIMEC
         G(5) = DBLE(IZ)/XA    ! GFAC
         G(6) = 6.E-06         ! FIEL
         G(7) = .6             ! POWER
         DO k = 1 , NEXPT ! Zn, An, E_p, THETA_lab, M_c, M_A, IAX, phi1, phi2, ikin, ln
            READ (JZB,*) IZ1(k) , XA1(k) , EP(k) , TLBDG(k) , EMMA(k) ,
     &           MAGA(k) , IAXS(k) , fi0 , fi1 , ISKIN(k) , LNORM(k)
            ITTE(k) = 0
            IF ( XA1(k).LT.0. ) ITTE(k) = 1
            XA1(k) = ABS(XA1(k))
            FIEX(k,1) = fi0/57.2957795 ! Convert to radians
            FIEX(k,2) = fi1/57.2957795
            IF ( TLBDG(k).LT.0. ) THEN
               FIEX(k,1) = FIEX(k,1) + 3.14159265
               FIEX(k,2) = FIEX(k,2) + 3.14159265
            ENDIF
         ENDDO

C     Else we don't recognize the suboption
      ELSE
         WRITE (22,99027) op1
99027    FORMAT (5X,'UNRECOGNIZED SUBOPTION',1X,1A4)
         GOTO 2000 ! Normal end of execution
      ENDIF
      GOTO 200 ! Get next suboption

C     Handle OP,ERRO      
 400  IF ( ICS.EQ.1 ) THEN
         REWIND 11
         DO kh1 = 1 , LP4
            READ (11) (CORF(kh1,kh2),kh2=1,LP6) ! LP6 = 32 (maximum number of gamma detectors)
         ENDDO
      ELSE
         CALL FTBM(0,chiss,idr,0,chilo,bten)
         REWIND 11
         DO kh1 = 1 , LP4
            WRITE (11) (CORF(kh1,kh2),kh2=1,LP6) ! LP6 = 32 (maximum number of gamma detectors)
         ENDDO
      ENDIF

      CALL FTBM(3,chiss,idr,1,chilo,bten)
      chis0 = chiss
      WRITE (22,99028) chis0
99028 FORMAT (1X///10X,'***** CENTRAL CHISQ=',1E12.4,1X,'*****'//)
      INHB = 1
      chisl = chiss
      DO kh = 1 , MEMAX
         HLM(kh) = ELM(kh)
      ENDDO
      IF ( idf.EQ.1 ) THEN
         IFBFL = 1
         IF ( irep.NE.2 ) GOTO 700
         IF ( iosr.EQ.0 ) GOTO 700
         REWIND IUNIT3
         READ (IUNIT3,*) ll , mm , kk , inn
         DO inn = 1 , ll
            READ (IUNIT3,*) mm , yyy , zz
         ENDDO
         DO inn = 1 , MEMAX
            READ (IUNIT3,*) mm , ll , kk
         ENDDO
         DO inn = 1 , MEMAX
            READ (IUNIT3,*) mm , yyy
         ENDDO
 450     READ (IUNIT3,*) mm , ll
         IF ( mm.EQ.0 ) THEN
            BACKSPACE IUNIT3
            GOTO 700
         ELSE
            READ (IUNIT3,*) kk , ll , yyy
            READ (IUNIT3,*) (SA(mm),mm=1,MEMAX)
            GOTO 450
         ENDIF
      ELSE
         naxfl = 0
         IF ( ms.EQ.0 ) mend = MEMAX
         IF ( ms.EQ.0 ) ms = 1
         DO kh = ms , mend ! Loop over matrix elements
            DO ij = 1 , 2
               pv = (ELMU(kh)-ELML(kh))/100.
               IF ( ij.NE.1 .OR. (ELM(kh)-ELML(kh)).GE.pv ) THEN
                  IF ( ij.NE.2 .OR. (ELMU(kh)-ELM(kh)).GE.pv ) THEN
                     DO kh1 = 1 , MEMAX
                        SA(kh1) = 0.
                     ENDDO
                     IF ( IVAR(kh).EQ.0 ) GOTO 500
                     SA(kh) = 1.*(-1)**ij
                     kh1 = kh
                     CALL KONTUR(idr,chis0,chisl,ifbp,-1,kh1,sh,bten,
     &                           rem)
                     ELM(kh) = HLM(kh)
                  ENDIF
               ENDIF
            ENDDO
            REWIND 15
            WRITE (15,*) (DEVD(ij),DEVU(ij),ij=1,MEMAX)
 500        CONTINUE
         ENDDO
      ENDIF
 600  IF ( ifbp.EQ.1 ) THEN
         REWIND 17
         DO lkj = 1 , MEMAX
            READ (17,*) ELM(lkj)
         ENDDO
         WRITE (22,99029)
99029    FORMAT (1X///20X,'*** BEST POINT FOUND (TAPE17) ***'///)
         CALL PRELM(3)
      ENDIF
      IF ( naxfl.EQ.0 ) WRITE (22,99051)
      IF ( naxfl.NE.0 ) WRITE (22,99050)
      WRITE (22,99030)
99030 FORMAT (40X,'ESTIMATED ERRORS'//5X,'INDEX',5X,'NI',5X,'NF',5X,
     &        'ME AND ERRORS'//)
      DO kh1 = 1 , MEMAX
         IF ( IVAR(kh1).NE.0 .AND. IVAR(kh1).LE.999 ) THEN
            WRITE (22,99031) kh1 , LEAD(1,kh1) , LEAD(2,kh1) , HLM(kh1)
     &                       , DEVD(kh1) , DEVU(kh1) , DEVD(kh1)
     &                       *100./ABS(HLM(kh1)) , DEVU(kh1)
     &                       *100./ABS(HLM(kh1))
99031       FORMAT (6X,1I3,5X,1I3,4X,1I3,5X,1F9.5,2X,'(',1F9.5,' ,',
     &              1F9.5,')','......',1F7.1,' ,',1F7.1,1X,'PC')
         ENDIF
      ENDDO
      IF ( naxfl.NE.0 ) WRITE (22,99050)
      IF ( naxfl.EQ.0 ) WRITE (22,99051)
      WRITE (22,99032)
99032 FORMAT (40X,'ESTIMATED ERRORS',//5X,'INDEX',5X,'NI',5X,'NF',5X,
     &        'B(E,ML)(OR QUADRUPOLE MOMENT)',' AND ERRORS'//)

      DO kh2 = 1 , MEMAX
         IF ( IVAR(kh2).NE.0 .AND. IVAR(kh2).LE.999 ) THEN
            ispa = LEAD(2,kh2)
            IF ( LEAD(1,kh2).NE.LEAD(2,kh2) ) THEN
               sbe = 2.*SPIN(ispa) + 1.
               be2 = HLM(kh2)*HLM(kh2)/sbe
               be2a = HLM(kh2) + DEVD(kh2)
               be2b = HLM(kh2) + DEVU(kh2)
               be2c = be2b
               IF ( ABS(be2a).GT.ABS(be2b) ) be2b = be2a
               IF ( ABS(be2a-be2c).LT.1.E-6 ) be2a = be2c
               IF ( be2a/HLM(kh2).LE.0. .OR. be2b/HLM(kh2).LE.0. )
     &              be2a = 0.
               be2a = be2a**2/sbe
               be2b = be2b**2/sbe
               WRITE (22,99052) kh2 , LEAD(2,kh2) , LEAD(1,kh2) , be2 , 
     &                          be2a - be2 , be2b - be2
            ELSE
               ispb = INT(SPIN(ispa))*2
               qfac = 3.170662*WTHREJ(ispb,4,ispb,-ispb,0,ispb)
               WRITE (22,99052) kh2 , LEAD(2,kh2) , LEAD(1,kh2) , 
     &                          HLM(kh2)*qfac , DEVD(kh2)*qfac , 
     &                          DEVU(kh2)*qfac
            ENDIF
         ENDIF
      ENDDO
      GOTO 2000 ! Normal end of execution

 700  irea = 0
      IF ( ms.LT.0 ) irea = 1
      IF ( ms.EQ.0 ) mend = MEMAX
      IF ( ms.EQ.0 ) ms = 1
 800  naxfl = 1
      IF ( irea.EQ.1 ) READ (JZB,*) ms , mend
      IF ( ms.NE.0 ) THEN
         DO kh = ms , mend ! For matrix elements
            IF ( ifc.NE.1 ) THEN
               REWIND 18
               DO kh1 = 1 , kh
                  READ (18,*) (KVAR(jyi),jyi=1,MEMAX)
               ENDDO
               DO kh1 = 1 , MEMAX ! For each matrix element
                  ivrh = IVAR(kh1)
                  IF ( KVAR(kh1).EQ.0 ) IVAR(kh1) = 0
                  KVAR(kh1) = ivrh
               ENDDO
            ENDIF
            DO ij = 1 , 2
               sh = DEVU(kh)
               IF ( ij.EQ.1 ) sh = DEVD(kh)
               IF ( ABS(sh).LT.1.E-6 ) sh = (-1)**ij*ABS(HLM(kh))/10.
               ELM(kh) = HLM(kh) + 1.5*sh
               mm = 0
               DO kh1 = 1 , MEMAX ! For each matrix element
                  IF ( ifc.EQ.1 ) KVAR(kh1) = IVAR(kh1)
                  mm = mm + IVAR(kh1)
               ENDDO
               IF ( mm.EQ.0 ) WRITE (22,99033) kh
99033          FORMAT (10X,'ME=',1I3,5X,'NO FREE MATRIX ELEMENTS')
               IF ( mm.NE.0 ) THEN
                  KFERR = 1
                  IF ( iosr.EQ.1 ) WRITE (IUNIT3,*) kh , kh ! For sigma program
                  IF ( iosr.EQ.1 ) WRITE (IUNIT3,*) kh , ij , ELM(kh)
                  LOCKS = 1
                  DLOCK = .05
                  CALL MINI(chiss,-1.D0,2,.0001D0,1000,idr,100000.D0,0,
     &                      iosr,kh,bten)
                  DO kh1 = 1 , MEMAX ! For each matrix element
                     SA(kh1) = (ELM(kh1)-HLM(kh1))/ABS(sh)
                  ENDDO
                  CALL KONTUR(idr,chis0,chisl,ifbp,inpo,kh,sh,bten,rem)
               ENDIF
               DO kh1 = 1 , MEMAX ! For each matrix element
                  IF ( ifc.EQ.1 ) IVAR(kh1) = KVAR(kh1)
                  ELM(kh1) = HLM(kh1)
               ENDDO
            ENDDO
            IF ( ifc.NE.1 ) THEN
               DO kh1 = 1 , MEMAX ! For each matrix element
                  IVAR(kh1) = KVAR(kh1)
               ENDDO
            ENDIF
            REWIND 15
            WRITE (15,*) (DEVD(kh1),DEVU(kh1),kh1=1,MEMAX)
         ENDDO ! Loop on matrix elements kh
         IF ( irea.EQ.1 ) GOTO 800
      ENDIF
      IF ( iosr.NE.0 ) THEN
         im = 0
         WRITE (IUNIT3,*) im , im
      ENDIF
      GOTO 600

 900  jfre = 0
      irfix = 0
      IF ( op2.EQ.'RE,F' ) irfix = 1
 1000 DO jrls = 1 , MEMAX ! For each matrix element
         IF ( IVAR(jrls).NE.0 .OR. irfix.NE.1 ) THEN
            IF ( IVAR(jrls).GT.999 ) THEN
               IF ( jfre.EQ.1 ) GOTO 1100
            ENDIF
            IVAR(jrls) = 2
            ELML(jrls) = -ABS(ELML(jrls))
            ELMU(jrls) = ABS(ELMU(jrls))
            IF ( jrls.GT.MEMX6 ) IVAR(jrls) = 1
         ENDIF
 1100    CONTINUE
      ENDDO ! For each matrix element jrls
      DO jrls = 1 , MEMAX
         ivarh(jrls) = IVAR(jrls)
      ENDDO
      GOTO 100 ! Back to input loop

 1200 CALL CMLAB(0,dsig,ttttt) ! Options MAP, STAR, POINT, MINI etc.
      IF ( ERR ) GOTO 2000 ! Error
      IF ( op2.EQ.'POIN' ) READ (JZB,*) ifwd , slim
      ient = 1
      icg = 1
      IF ( SPIN(1).LT.1.E-6 ) ISO = 0
      IF ( iobl.LT.1 ) THEN
         IF ( op2.NE.'GOSI' ) THEN
            iapx = 0
            DO ii = 1 , LP6 ! LP6 = 32 (maximum number of gamma detectors)
               ILE(ii) = 1
            ENDDO
            nch = 0
            DO jexp = 1 , NEXPT ! For each experiment
               IEXP = jexp
               ttttt = TREP(IEXP)
               dsig = DSIGS(IEXP)
               IF ( op2.NE.'STAR' ) THEN
                  jmm = IEXP
                  IF ( IEXP.NE.1 ) THEN
                     DO lli = 1 , LP6 ! LP6 = 32 (maximum number of gamma detectors)
                        ILE(lli) = ILE(lli) + NYLDE(IEXP-1,lli)
                     ENDDO
                  ENDIF
               ENDIF
               fi0 = FIEX(IEXP,1) ! Lower phi limit
               fi1 = FIEX(IEXP,2) ! Upper phi limit
               CALL LOAD(IEXP,1,icg,0.D0,jj)
               CALL ALLOC(ACCUR)
               CALL SNAKE(IEXP,ZPOL)
               CALL SETIN
               DO j = 1 , LMAX ! For each spin up to ground-state spin + 1
                  polm = DBLE(j-1) - SPIN(1)
                  CALL LOAD(IEXP,2,icg,polm,jj)
                  CALL STING(jj)
                  CALL PATH(jj)
                  CALL INTG(IEXP)
                  CALL TENB(j,bten,LMAX)
                  pr = 0.
                  IF ( op2.EQ.'STAR' .OR. IPRM(19).EQ.1 )
     &                 WRITE (22,99034) (DBLE(j)-1.-SPIN(1)) , IEXP
99034             FORMAT (1X//40X,'EXCITATION AMPLITUDES'//10X,'M=',
     &                    1F5.1,5X,'EXPERIMENT',1X,1I2//5X,'LEVEL',2X,
     &                    'SPIN',2X,'M',5X,'REAL AMPLITUDE',2X,
     &                    'IMAGINARY AMPLITUDE'//)
                  DO k = 1 , ISMAX ! For substates
                     pr = pr + DBLE(ARM(k,5))**2 + DIMAG(ARM(k,5))**2
                     IF ( op2.EQ.'STAR' .OR. IPRM(19).EQ.1 )
     &                    WRITE (22,99035) INT(CAT(k,1)) , CAT(k,2) , 
     &                    CAT(k,3) , DBLE(ARM(k,5)) , DIMAG(ARM(k,5))
99035                FORMAT (7X,1I2,3X,1F4.1,2X,1F5.1,2X,1E14.6,2X,
     &                       1E14.6)
                  ENDDO ! Loop on substates k
                  IF ( op2.EQ.'STAR' .OR. IPRM(19).EQ.1 )
     &                 WRITE (22,99036) pr
99036             FORMAT (1X/5X,'SUM OF PROBABILITIES=',1E14.6)
               ENDDO ! Loop over spins j
               CALL TENS(bten)
               IF ( itno.NE.0 ) THEN ! write statistical tensors on tape 17
                  DO k = 2 , NMAX
                     WRITE (17,*) k
                     DO kk = 1 , 4
                        in1 = (k-1)*28 + 1 + (kk-1)*7
                        in2 = in1 + 2*kk - 2
                        WRITE (17,*) (ZETA(kkk),kkk=in1,in2)
                     ENDDO
                  ENDDO
               ENDIF
               summm = 0.
               DO jgl = 2 , NMAX
                  loct = (jgl-1)*28 + 1
                  summm = summm + ZETA(loct)
               ENDDO
               pop1 = 1. - summm
               jgl = 1
               IF ( op2.EQ.'STAR' .OR. IPRM(19).EQ.1 ) WRITE (22,99053)
     &              jgl , pop1
               DO jgl = 2 , NMAX
                  loct = (jgl-1)*28 + 1
                  IF ( op2.EQ.'STAR' .OR. IPRM(19).EQ.1 )
     &                 WRITE (22,99053) jgl , ZETA(loct)
               ENDDO
               IF ( op2.NE.'STAR' ) THEN
                  CALL DECAY(ccd,0,ccc)
                  nogeli = NANG(IEXP) ! Number of detector angles for expt
                  jgl1 = 0
                  DO js = 1 , LP2 ! LP2 = 1500 (maximum number of matrix elements)
                     DO jgl = 1 , 20
                        SUMCL(jgl,js) = 0.
                     ENDDO
                  ENDDO
                  DO jgl = 1 , nogeli ! For each detector angle
                     IF ( IRAWEX(IEXP).NE.0 ) THEN
                        IF ( op2.EQ.'POIN' .AND. IPRM(20).EQ.1 )
     &                       WRITE (23,99037) IEXP , jgl , EP(IEXP) , 
     &                       TLBDG(IEXP)
99037                   FORMAT (1x//50x,'CALCULATED YIELDS'//5x,
     &                          'EXPERIMENT ',1I2,2x,'DETECTOR ',1I2/5x,
     &                          'ENERGY ',1F10.3,1x,'MEV',2x,'THETA ',
     &                          1F7.3,1x,'DEG'//5x,'NI',5x,'NF',5x,'II',
     &                          5x,'IF',5x,'E(MeV)',5x,'EFFICIENCY'/)
                     ENDIF
                     gth = AGELI(IEXP,jgl,1)
                     figl = AGELI(IEXP,jgl,2)
                     fm = (fi0+fi1)/2.
                     CALL ANGULA(YGN,idr,1,fi0,fi1,ttttt,gth,figl,jgl,
     &                 op2)
                     IF ( IFMO.NE.0 ) THEN
                        id = ITMA(IEXP,jgl) ! Get identity of detector
                        d = ODL(id) ! Get results of OP,GDET for that detector
                        rx = d*SIN(gth)*COS(figl-fm) - .25*SIN(ttttt)
     &                       *COS(fm)
                        ry = d*SIN(gth)*SIN(figl-fm) - .25*SIN(ttttt)
     &                       *SIN(fm)
                        rz = d*COS(gth) - .25*COS(ttttt)
                        rl = SQRT(rx*rx+ry*ry+rz*rz)
                        thc = TACOS(rz/rl)
                        sf = d*d/rl/rl
                        fic = ATAN2(ry,rx)
                        CALL ANGULA(YGP,idr,1,fi0,fi1,ttttt,thc,fic,jgl,
     &                    op2)
                        DO ixl = 1 , idr
                           ixm = KSEQ(ixl,3)
                           tfac = TAU(ixm)
                           YGN(ixl) = YGN(ixl)
     &                                + .01199182*tfac*BETAR(IEXP)
     &                                *(sf*YGP(ixl)-YGN(ixl))
                        ENDDO
                     ENDIF
                     IF ( IRAWEX(IEXP).NE.0 ) THEN
                        ipd = ITMA(IEXP,jgl) ! Get identity of detector
                        DO jyi = 1 , idr ! For each decay
                           ni = KSEQ(jyi,3)
                           nf = KSEQ(jyi,4)
                           decen = EN(ni) - EN(nf)
                           cocos = SIN(ttttt)*SIN(gth)*COS(fm-figl)
     &                             + COS(ttttt)*COS(gth)
                           decen = decen*(1.+BETAR(IEXP)*cocos)
                           CALL EFFIX(ipd,decen,effi)
                           IF ( op2.EQ.'POIN' .AND. IPRM(20).EQ.1 )
     &                          WRITE (23,99049) ni , nf , SPIN(ni) , 
     &                                 SPIN(nf) , decen , effi
                           YGN(jyi) = YGN(jyi)*effi
                        ENDDO
                        inclus = ICLUST(IEXP,jgl) ! Cluster number for detector jgl
                        IF ( inclus.NE.0 ) THEN
                           DO jyi = 1 , idr ! For each decay
                              SUMCL(inclus,jyi) = SUMCL(inclus,jyi)
     &                           + YGN(jyi)
                           ENDDO
                           IF ( jgl.NE.LASTCL(IEXP,inclus) ) GOTO 1205 ! If it is not the last detector in the cluster
                           DO jyi = 1 , idr ! For each decay
                              YGN(jyi) = SUMCL(inclus,jyi)
                           ENDDO
                        ENDIF
                     ENDIF
                     jgl1 = jgl1 + 1
                     lu = ILE(jgl1)
                     IF ( op2.EQ.'POIN' .OR. IPRM(11).EQ.1 )
     &                    WRITE (22,99048) IEXP , jgl1 , EP(IEXP) , 
     &                    TLBDG(IEXP)
                     jmm = 0
C---- this bit removed in gosia2 start
                     ttttx = TLBDG(IEXP)/57.2957795
                     YGN(IDRN) = YGN(IDRN)*dsig*SIN(ttttx)
                     DO jyi = 1 , idr
                        IF ( jyi.NE.IDRN ) YGN(jyi) = YGN(jyi)
     &                       *dsig*SIN(ttttx)
                     ENDDO
C---- this bit removed in gosia2 end
                     DO jyi = 1 , idr
                        ni = KSEQ(jyi,3)
                        nf = KSEQ(jyi,4)
                        IF ( op2.EQ.'POIN' .OR. IPRM(11).EQ.1 )
     &                       WRITE (22,99049) ni , nf , SPIN(ni) , 
     &                       SPIN(nf) , YGN(jyi) , YGN(jyi)/YGN(IDRN)
                        IF ( ifwd.EQ.1 ) THEN
                           IF ( (YGN(jyi)/YGN(IDRN)).GE.slim ) THEN
                              IF ( jgl1.EQ.1 ) sh1 = YGN(IDRN)
                              jmm = jmm + 1
                              CORF(jmm,1) = DBLE(ni)
                              CORF(jmm,2) = DBLE(nf)
                              CORF(jmm,3) = YGN(jyi)/sh1 ! Not divided by sh1 in gosia2
                              IF ( YGN(jyi).GE.YGN(IDRN) ) CORF(jmm,4)
     &                             = CORF(jmm,3)/20.
                              IF ( YGN(jyi).LT.YGN(IDRN) ) CORF(jmm,4)
     &                             = CORF(jmm,3)
     &                             *(.05+.2*(1.-YGN(jyi)/YGN(IDRN)))
                           ENDIF
                        ENDIF
                        IF ( op2.EQ.'CORR' ) THEN
                           READ (15,*) yydd
                           nch = nch + 1
                           jjjj = IY(lu,jgl1)/1000
                           jyi1 = IY(lu,jgl1) - jjjj*1000
                           IF ( IY(lu,jgl1).EQ.jyi .OR. jjjj.EQ.jyi .OR.
     &                          jyi1.EQ.jyi ) THEN
                              IF ( IY(lu,jgl1).GE.1000 ) THEN
                                 jyi2 = jyi1 - jjjj
                                 IF ( jyi2.LE.0 ) GOTO 1202
                                 DO ihuj = 1 , jyi2
                                    READ (15,*) yyd1
                                 ENDDO
                                 yydd = yydd + yyd1
                                 YGN(jyi) = YGN(jyi) + YGN(jyi1)
                                 REWIND 15
                                 DO ihuj = 1 , nch
                                    READ (15,*) yyd1
                                 ENDDO
                              ENDIF
                              IF ( IEXP.EQ.1 .AND. lu.EQ.NYLDE(1,1)
     &                             .AND. jgl1.EQ.1 )
     &                             cnst = yydd/YGN(jyi)
                              CORF(lu,jgl1) = YEXP(jgl1,lu)
                              YEXP(jgl1,lu) = YEXP(jgl1,lu)
     &                           /yydd*YGN(jyi)
                              DYEX(jgl1,lu) = DYEX(jgl1,lu)
     &                           /yydd*YGN(jyi)
                              lu = lu + 1
                           ENDIF
                        ENDIF
 1202                   CONTINUE
                     ENDDO
                     IF ( ifwd.EQ.1 ) THEN
                        xw = 1.
                        WRITE (4,*) IEXP , jgl1 , ABS(IZ1(IEXP)) , 
     &                              ABS(XA1(IEXP)) , ABS(EP(IEXP)) , 
     &                              jmm , xw
                        DO jyi = 1 , jmm
                           WRITE (4,*) INT(CORF(jyi,1)) , 
     &                                 INT(CORF(jyi,2)) , CORF(jyi,3) , 
     &                                 CORF(jyi,4)
                        ENDDO
                     ENDIF
 1205                CONTINUE
                  ENDDO ! Loop on detector angles jgl
                  IF ( op2.EQ.'CORR' ) THEN
                     jgl1 = 0
                     DO jgl = 1 , nogeli ! For each detector
                        IF ( IRAWEX(jexp).NE.0 ) THEN
                           inclus = ICLUST(jexp,jgl) ! Cluster number for detector jgl
                           IF ( inclus.NE.0 ) THEN
                              IF ( jgl.NE.LASTCL(jexp,inclus) ) ! If detector is not the last in the cluster
     &                             GOTO 1206
                           ENDIF
                        ENDIF
                        jgl1 = jgl1 + 1
                        READ (3,*) ne , na , zp , ap , xep , nval , waga
                        WRITE (4,*) ne , na , zp , ap , EP(IEXP) , 
     &                              nval , waga
                        WRITE (22,99038) IEXP , jgl1
99038                   FORMAT (///10X,'EXPERIMENT',1X,I2,8X,'DETECTOR',
     &                          1X,I2,//9X,'NI',5X,'NF',5X,'YEXP',8X,
     &                          'YCOR',8X,'COR.F'/)
                        ile1 = ILE(jgl1)
                        DO itp = 1 , nval
                           READ (3,*) ns1 , ns2 , fiex1(1,1,1) , 
     &                                fiex1(1,1,2)
                           ltrn = IY(ile1+itp-1,jgl1)
                           IF ( ltrn.LT.1000 ) THEN
                              ns1 = KSEQ(ltrn,3)
                              ns2 = KSEQ(ltrn,4)
                           ELSE
                              ltrn1 = ltrn/1000
                              ns1 = KSEQ(ltrn1,3)*100
                              ns2 = KSEQ(ltrn1,4)*100
                              ltrn2 = ltrn - ltrn1*1000
                              ns1 = ns1 + KSEQ(ltrn2,3)
                              ns2 = ns2 + KSEQ(ltrn2,4)
                           ENDIF
                           ycorr = YEXP(jgl1,ile1+itp-1)*cnst ! Not multiplied by cnst in gosia2
                           WRITE (4,*) ns1 , ns2 , ycorr , 
     &                                 DYEX(jgl1,ile1+itp-1)*cnst ! Not multiplied by cnst in gosia2
                           WRITE (22,99039) ns1 , ns2 , 
     &                            CORF(ile1+itp-1,jgl1) , ycorr , 
     &                            ycorr/CORF(ile1+itp-1,jgl1)
99039                      FORMAT (5X,I4,5X,I4,3X,E8.3,4X,E8.3,4X,E8.3)
                        ENDDO ! Loop over itp
 1206                   CONTINUE
                     ENDDO ! Loop over jgl
                  ENDIF ! if ( op2.EQ. 'CORR')
               ENDIF
            ENDDO ! Loop over jexp
            IF ( op2.EQ.'STAR' ) oph = op2
            IF ( op2.NE.'STAR' ) THEN
               IF ( op2.EQ.'CORR' ) THEN
                  ntap = 4
                  CALL READY(idr,ntap,ipri)
                  REWIND ntap
               ENDIF
            ENDIF
            GOTO 100 ! Back to input loop
         ENDIF ! if (op2 .NE. 'GOSI') if statement
      ENDIF ! if ( iobl.LT.1 ) if statement

 1300 IF ( iobl.GE.1 ) THEN ! OP,ERRO
         ient = 1
         icg = 2
         nmaxh = NMAX
         lmax1 = LMAX
         sh1 = SPIN(1) ! Save ground-state spin
         sh2 = SPIN(2) ! Save spin of first excited state
         ih1 = IFAC(1)
         ih2 = IFAC(2)
         magh = MAGEXC
         lmaxh = LMAXE
         isoh = ISO
         ISO = 0
         eh1 = ELM(1)
         lh1 = LEAD(1,1)
         lh2 = LEAD(2,1)
         lamh = LAMMAX
         memh = MEMAX
         DO kh = 1 , 8 ! For each multipolarity
            ihlm(kh) = MULTI(kh)
            ihlm(kh+24) = LDNUM(kh,2)
            ihlm(kh+8) = LAMDA(kh)
            ihlm(kh+16) = LDNUM(kh,1)
         ENDDO
         DO jexp = 1 , NEXPT ! For each experiment
            IEXP = jexp
            intvh = INTERV(IEXP)
            DO jgs = 1 , MEMAX
               DO jgr = 1 , 7
                  QAPR(jgs,1,jgr) = 0.
               ENDDO
            ENDDO
            DO iuy = 1 , 6
               XIR(iuy,IEXP) = 0.
            ENDDO
            emhl1 = EMMA(IEXP)
            EMMA(IEXP) = DBLE(MAGA(IEXP))
            jde = 2
            IF ( MAGA(IEXP).EQ.0 ) jde = 1
            DO iuy = 1 , 6
               zmir(iuy,1,IEXP) = 0.
               zmir(iuy,2,IEXP) = 0.
            ENDDO
            CALL LOAD(IEXP,1,2,0.D0,jj)
            DO jgs = 1 , LMAX ! For each spin up to ground-state spin + 1
               polm = DBLE(jgs-1) - SPIN(1)
               CALL LOAD(IEXP,3,2,polm,jj)
               CALL PATH(jj)
               CALL LOAD(IEXP,2,2,polm,jj)
               ictl = 1
               DO kk = 1 , 6
                  ll = ihlm(kk)
                  IF ( ll.NE.0 ) THEN
                     lfini = ll + ictl - 1
                     ict = ictl
                     DO lll = ict , lfini
                        ictl = ictl + 1
                        IF ( jgs.EQ.1 ) XIR(kk,IEXP)
     &                       = MAX(XIR(kk,IEXP),ABS(XI(lll)))
                        r1 = ABS(QAPR(lll,1,1))
                        r2 = ABS(QAPR(lll,1,4))
                        r3 = ABS(QAPR(lll,1,7))
                        rm = MAX(r1,r2,r3)
                        bmx = MAX(ABS(ELMU(lll)),ABS(ELML(lll)))
                        zmir(kk,2,IEXP)
     &                     = MAX(zmir(kk,2,IEXP),rm*bmx/ABS(ELM(lll)),
     &                     rm)
                        r1 = ABS(QAPR(lll,1,2))
                        r2 = ABS(QAPR(lll,1,3))
                        r3 = ABS(QAPR(lll,1,5))
                        r4 = ABS(QAPR(lll,1,6))
                        rm = MAX(r1,r2,r3,r4)
                        zmir(kk,1,IEXP)
     &                     = MAX(zmir(kk,1,IEXP),rm*bmx/ABS(ELM(lll)),
     &                     rm)
                     ENDDO
                     IF ( zmir(kk,1,IEXP).LT..5 ) zmir(kk,1,IEXP) = .5
                     IF ( zmir(kk,2,IEXP).LT..5 ) zmir(kk,2,IEXP) = .5
                  ENDIF
               ENDDO
            ENDDO
            DO kk = 1 , 6
               XIR(kk,IEXP) = XIR(kk,IEXP)*1.01
               DO kh = 1 , 8
                  MULTI(kh) = 0
                  LAMDA(kh) = 0
                  LDNUM(kh,2) = 0
                  LDNUM(kh,1) = 0
               ENDDO
               NMAX = 2
               ELM(1) = 1.
               LEAD(1,1) = 1
               LEAD(2,1) = 2
               SPIN(1) = 0.
               IFAC(1) = 1
               LAMMAX = 1
               MEMAX = 1
               MAGEXC = 0
               kkk = 0
               icg = 1
               IF ( ihlm(kk).NE.0 ) THEN
                  MULTI(kk) = 1
                  LAMDA(1) = kk
                  SPIN(2) = DBLE(kk)
                  IFAC(2) = 1
                  LDNUM(kk,1) = 1
                  icg = 1
                  CALL LOAD(IEXP,1,icg,0.D0,jj)
                  CALL LOAD(IEXP,2,icg,0.D0,jj)
                  CALL PATH(1)
                  sz1 = MIN(zmir(kk,1,IEXP),10.D0)
                  sz2 = zmir(kk,2,IEXP)/50.
                  acof = 2.4009604E-3/zmir(kk,2,IEXP)
                  bcof = 8.163265E-4
                  DO jd = 1 , jde
                     nksi = 5
                     IF ( jd.EQ.2 ) nksi = 10
                     IF ( MAGA(IEXP).EQ.0 ) nksi = 10
                     DO jk = 1 , 3
                        ZETA(jk) = 0.
                     ENDDO
                     nz = 50
                     IF ( jd.EQ.1 .AND. MAGA(IEXP).NE.0 ) nz = 1
                     DO jk = 1 , nksi
                        XI(1) = XIR(kk,IEXP)*(jk-1)/(nksi-1)
                        IF ( jk.EQ.1 ) XI(1) = .02
                        s11 = 0.
                        s21 = 0.
                        s12 = 0.
                        s22 = 0.
                        ph1 = 0.
                        ph2 = 0.
                        DO jz = 1 , nz
                           ZETA(jd) = sz2*jz
                           IF ( jd.EQ.1 .AND. MAGA(IEXP).NE.0 ) ZETA(jd)
     &                          = sz1
                           IF ( ZETA(jd).LT..1 ) INTERV(IEXP) = 1000
                           IF ( ZETA(jd).GE..1 ) INTERV(IEXP) = intvh
                           CALL ALLOC(ACCUR)
                           CALL SNAKE(IEXP,ZPOL)
                           CALL SETIN
                           CALL STING(1)
                           IF ( kk.GT.2 ) THEN
                              ARM(1,5) = (.9999999,0.)
                              ARM(2,5) = (1.2E-6,0.)
                              ARM(1,6) = (.9999998,0.)
                              ARM(2,6) = (.9E-6,0.)
                              DO kh = 1 , 4
                                 ARM(1,kh) = (-1.E-6,0.)
                                 ARM(2,kh) = (1.E-6,0.)
                              ENDDO
                           ENDIF
                           CALL INTG(IEXP)
                           jp = 2
                           IF ( MAGA(IEXP).NE.0 .AND. jd.EQ.2 ) jp = 3
                           p = DBLE(ARM(1,5))
                           r = DIMAG(ARM(1,5))
                           qr = DBLE(ARM(jp,5))
                           s = DIMAG(ARM(jp,5))
                           test = p*p + r*r + qr*qr + s*s
                           p = p/SQRT(test)
                           s = ABS(r/s)
                           IF ( jk.EQ.1 ) THEN
                              IF ( MAGA(IEXP).EQ.0 ) THEN
                                 q1 = 0.
                                 GOTO 1302
                              ELSEIF ( jd.EQ.2 .OR. MAGA(IEXP).EQ.0 )
     &                                 THEN
                                 q1 = 0.
                                 GOTO 1302
                              ENDIF
                           ENDIF
                           q1 = ARCTG(s,ph1,pi)
                           ph1 = q1
 1302                      IF ( jk.EQ.1 ) THEN
                              IF ( jd.EQ.1 .AND. MAGA(IEXP).NE.0 ) THEN
                                 q2 = 0.
                                 GOTO 1304
                              ENDIF
                           ENDIF
                           q2 = ARCCOS(p,ph2,pi)
                           ph2 = q2
 1304                      q1 = q1/ZETA(jd)/2.
                           q2 = q2/ZETA(jd)
                           IF ( jd.EQ.1 .AND. MAGA(IEXP).NE.0 ) q2 = -q2
                           IF ( jd.NE.1 .OR. MAGA(IEXP).EQ.0 ) THEN
                              s11 = s11 + q1
                              s12 = s12 + q1*jz
                              s21 = s21 + q2
                              s22 = s22 + jz*q2
                           ENDIF
                        ENDDO
                        IF ( jd.EQ.1 .AND. MAGA(IEXP).NE.0 ) THEN
                           PARX(IEXP,2*kk-1,jk) = q1
                           PARX(IEXP,2*kk,jk) = q2
                        ELSE
                           PARXM(IEXP,1,jk,kk) = acof*(2.*s12-51.*s11)
                           PARXM(IEXP,2,jk,kk) = bcof*(101.*s11-3.*s12)
                           PARXM(IEXP,3,jk,kk) = acof*(2.*s22-51.*s21)
                           PARXM(IEXP,4,jk,kk) = bcof*(101.*s21-3.*s22)
                        ENDIF
                     ENDDO ! Loop over jk
                  ENDDO ! Loop over jd
               ENDIF
            ENDDO ! Loop over kk
            EMMA(IEXP) = emhl1
            NMAX = nmaxh
            SPIN(1) = sh1 ! Restore ground-state spin
            SPIN(2) = sh2 ! Restore spin of first excited state
            IFAC(1) = ih1
            IFAC(2) = ih2
            MAGEXC = magh
            ISO = isoh
            ELM(1) = eh1
            LEAD(1,1) = lh1
            LEAD(2,1) = lh2
            LAMMAX = lamh
            MEMAX = memh
            DO kh = 1 , 8 ! For each multipolarity
               LDNUM(kh,2) = ihlm(kh+24)
               MULTI(kh) = ihlm(kh)
               LAMDA(kh) = ihlm(kh+8)
               LDNUM(kh,1) = ihlm(kh+16)
            ENDDO
            INTERV(IEXP) = intvh
         ENDDO ! Loop over experiments jexp

         irix = 7
         REWIND irix
         DO iuy = 1 , 6
            WRITE (irix,*) (XIR(iuy,jj),jj=1,NEXPT)
            WRITE (irix,*) (zmir(iuy,1,jj),zmir(iuy,2,jj),jj=1,NEXPT)
         ENDDO
         DO jj = 1 , NEXPT ! For each experiment
            DO jk = 1 , 4
               DO kuku = 1 , 6
                  WRITE (irix,*) (PARXM(jj,jk,jl,kuku),jl=1,10)
               ENDDO
            ENDDO
            DO jk = 1 , 12
               WRITE (irix,*) (PARX(jj,jk,jl),jl=1,5)
            ENDDO
         ENDDO
         DO jj = 1 , 2
            DO jj1 = 1 , LP1 ! LP1 = 50 (maximum number of experiments)
               IDIVE(jj1,jj) = 1
            ENDDO
         ENDDO
      ELSE ! iobl .lt. 1
         irix = 7
         REWIND irix
         DO iuy = 1 , 6
            READ (irix,*) (XIR(iuy,jj),jj=1,NEXPT)
            READ (irix,*) (zmir(iuy,1,jj),zmir(iuy,2,jj),jj=1,NEXPT)
         ENDDO
         DO jj = 1 , NEXPT ! For each experiment
            DO jk = 1 , 4
               DO kuku = 1 , 6
                  READ (irix,*) (PARXM(jj,jk,jl,kuku),jl=1,10)
               ENDDO
            ENDDO
            DO jk = 1 , 12
               READ (irix,*) (PARX(jj,jk,jl),jl=1,5)
            ENDDO
         ENDDO
         DO jgs = 1 , MEMAX ! For each matrix element
            DO jgr = 1 , 7
               QAPR(jgs,1,jgr) = 0.
            ENDDO
         ENDDO
      ENDIF

C     Handle map
      IF ( IPRM(12).NE.0 ) THEN ! gosia2 has additional .OR. op2 .eq 'MAP '
         IPRM(12) = 0
         DO jex = 1 , NEXPT
            DO lex = 1 , 6
               IF ( MULTI(lex).NE.0 ) THEN
                  WRITE (22,99040) jex , XIR(lex,jex)
99040             FORMAT (1X//30X,'EXPERIMENT',1X,1I2,10X,'MAX.XI=',
     &                    1F6.4)
                  WRITE (22,99041) lex , zmir(lex,2,jex)
99041             FORMAT (1X/30X,'E',1I1,8X,'MI=0',5X,'MAX.ZETA=',
     &                    1F6.3//)
                  WRITE (22,99054)
                  DO kex = 1 , 10
                     xxi = XIR(lex,jex)*(kex-1)/9.
                     WRITE (22,99055) xxi , 
     &                                (PARXM(jex,ilx,kex,lex),ilx=1,4)
                  ENDDO
                  IF ( MAGA(jex).NE.0 ) THEN
                     WRITE (22,99042) lex , zmir(lex,1,jex)
99042                FORMAT (1X//30X,'E',1I1,8X,'MI=+/-1',5X,
     &                       'MAX.ZETA=',1F6.3//)
                     WRITE (22,99054)
                     DO kex = 1 , 5
                        xxi = XIR(lex,jex)*(kex-1)/4.
                        u = 0.
                        WRITE (22,99055) xxi , u , PARX(jex,2*lex-1,kex)
     &                         , u , PARX(jex,2*lex,kex)
                     ENDDO ! Loop on kex
                  ENDIF ! if maga(jex).ne.0
               ENDIF ! if multi(lex).ne.0
            ENDDO ! Loop on lex
         ENDDO ! Loop on jex
      ENDIF ! IPRM(12).ne.0
      IF ( op2.NE.'GOSI' .AND. op2.NE.'ERRO' ) GOTO 100 ! Back to input loop
      IF ( op2.EQ.'ERRO' ) GOTO 400

 1400 DO kh1 = 1 , MEMAX
         HLM(kh1) = ELM(kh1)
      ENDDO
      lfagg = 0
      DO kh1 = 1 , MEMAX
         IVAR(kh1) = ivarh(kh1)
      ENDDO
      CALL MINI(chisq,chiok,nptl,conu,imode,idr,xtest,0,0,0,bten)
      IF ( IPS1.EQ.0 ) GOTO 2000 ! Normal end of execution
      IMIN = IMIN + 1
      DO iva = 1 , LP1 ! LP1 = 50 (maximum number of experiments)
         JSKIP(iva) = 1
      ENDDO
      irix = 12
      REWIND irix
      DO lkj = 1 , MEMAX
         WRITE (irix,*) ELM(lkj)
      ENDDO
      IF ( ifm.EQ.1 ) CALL PRELM(3) ! ifm = fast minimisation switch
      IF ( ifm.NE.1 ) GOTO 100 ! Back to input loop
      GOTO 2000 ! Normal end of execution

C.............................................................................
 1500 WRITE (22,99043)
99043 FORMAT (5X,'ERROR-M.E. DOES NOT BELONG TO THE UPPER TRIANGLE')
      GOTO 1900 ! Troubleshoot

C.............................................................................
 1600 WRITE (22,99044)
99044 FORMAT (5X,'ERROR-WRONG SEQUENCE OF MULTIPOLARITIES')
      GOTO 1900 ! Troubleshoot

C.............................................................................
 1700 WRITE (22,99045)
99045 FORMAT (5X,'ERROR-REPEATED APPEARANCE OF THE STATE')
      GOTO 1900 ! Troubleshoot

C.............................................................................
 1800 WRITE (22,99046)
99046 FORMAT (1X///10X,'ERROR-INSUFFICIENT SPACE FOR E-THETA INTEGR ',
     &        'ATION')
      GOTO 1900 ! Troubleshoot

C.............................................................................
C     Troubleshooting
 1900 IF ( ITS.NE.0 ) THEN
         iva = 0
         WRITE (18,*) iva , iva , iva , chisq
         IF ( ITS.NE.2 ) THEN
            WRITE (15,*) iva , chisq , chisq , chisq , chisq
            CALL KLOPOT(kmat,rlr) ! Troubleshooting
         ENDIF
      ENDIF

C     End of execution
 2000 WRITE (22,99047)
99047 FORMAT (15X,'********* END OF EXECUTION **********')
      STOP

C     Handle OP,EXIT
 430  IF ( IPRM(18).NE.0 ) CALL PTICC(idr)
      IF ( oph.EQ.'GOSI' ) THEN
         IF ( lfagg.NE.1 ) THEN
            IF ( IMIN.NE.0 ) THEN
               IF ( IPRM(4).EQ.-1 ) IPRM(4) = 111111
               iskok = IPRM(7) + IPRM(8) + IPRM(13) + IPRM(14)
               IF ( iskok.NE.0 .OR. IPRM(4).NE.111111 ) THEN
                  IF ( iskok.NE.0 ) THEN
                     IF ( IPRM(7).EQ.1 ) IPRM(7) = -1
                     IF ( IPRM(8).EQ.1 ) IPRM(8) = -1
                     IF ( IPRM(3).EQ.1 .AND. NBRA.NE.0 ) IPRM(3) = -1
                     IF ( IPRM(13).EQ.1 ) IPRM(13) = -1
                     IF ( IPRM(14).EQ.1 ) IPRM(14) = -1
                  ENDIF
                  CALL MINI(chisq,chiok,+1,conu,2000,idr,xtest,2,0,0,
     &                      bten)
               ENDIF
            ENDIF
            CALL MIXR(iva,1,chisq,chilo)
            IF ( IPRM(15).NE.0 .AND. KFERR.NE.1 .AND. iyr.NE.0 ) THEN
               WRITE (22,99011)
99011          FORMAT (1X//20X,'CALCULATED LIFETIMES'//5X,'LEVEL',5X,
     &                 'LIFETIME(PSEC)',5X,'EXP',8X,'ERROR'/)
               DO iva = 2 , NMAX
                  DO iva1 = 1 , NLIFT
                     IF ( LIFCT(iva1).EQ.iva ) GOTO 122
                  ENDDO
                  WRITE (22,99012) iva , TAU(iva)
99012             FORMAT (6X,1I3,7X,1E10.4)
                  GOTO 124
 122              WRITE (22,99013) iva , TAU(iva) , TIMEL(1,iva1) , 
     &                   TIMEL(2,iva1)
99013             FORMAT (6X,1I3,7X,1E10.4,5X,1E10.4,4X,1E10.4)
 124              IF ( iva.EQ.NMAX ) THEN
                     IF ( NAMX.GE.1 ) THEN
                        WRITE (22,99014)
99014                   FORMAT (5x,//,
     &         'CALCULATED AND EXPERIMENTAL MATRIX ELEMENTS'
     &         ,//)
                        WRITE (22,99015)
99015                   FORMAT (5x,'NI ','NF ',' EXP. ME   ',
     &                     'CURRENT ME','   SIGMA')
                        DO kq = 1 , NAMX
                           ni = IAMY(kq,1)
                           nf = IAMY(kq,2)
                           ind = IAMX(kq)
                           ess = ELM(ind)
                           esd = EAMX(kq,1)
                           dsd = EAMX(kq,2)
                           WRITE (22,99016) ni , nf , esd , ess , 
     &                        (ess-esd)/dsd
99016                      FORMAT (4x,1I3,1x,1I3,1x,1F9.4,1x,1F9.4,1x,
     &                              1F9.4)
                        ENDDO
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
            IF ( IMIN.NE.0 ) CALL PRELM(3)
         ENDIF
      ENDIF
      GOTO 1900 ! End of OP,EXIT - troubleshoot

99048 FORMAT (1X//50X,'CALCULATED YIELDS'//5X,'EXPERIMENT ',1I2,2X,
     &        'DETECTOR ',1I2/5X,'ENERGY ',1F10.3,1X,'MEV',2X,'THETA ',
     &        1F7.3,1X,'DEG'//5X,'NI',5X,'NF',5X,'II',5X,'IF',5X,
     &        'YIELD',5X,'NORMALIZED YIELD'/)
99049 FORMAT (4X,1I3,4X,1I3,3X,1F4.1,3X,1F4.1,3X,1E11.5,3X,1E11.5)
99050 FORMAT (1X///44X,'OVERALL')
99051 FORMAT (1X///43X,'DIAGONAL')
99052 FORMAT (6X,1I3,5X,1I3,4X,1I3,5X,1F10.5,2X,'(',1F10.5,' ,',1F10.5,
     &        ')')
99053 FORMAT (2X,'LEVEL',1X,1I2,10X,'POPULATION',1X,1E14.6)
99054 FORMAT (5X,'XI',13X,'Q1',22X,'Q2'///13X,'SLOPE',2X,'INTERCEPT',7X,
     &        'SLOPE',5X,'INTERCEPT'//)
99055 FORMAT (2X,1F6.4,3X,1E8.2,2X,1E8.2,6X,1E8.2,2X,1E8.2)
      END
 
C----------------------------------------------------------------------
C FUNCTION ARCOS
C
C Called by: GOSIA
C Calls:     TACOS
C
C Purpose: calculates an arccosine in a particular range
C
C Formal parameters:
C      A      - argument
C      F      - range
C      Pi     - Pi must be set to 3.14159... before calling ARCCOS
C
C Return value:
C      arccosine(A) within range of F
 
      REAL*8 FUNCTION ARCCOS(A,F,Pi)
      IMPLICIT NONE
      REAL*8 A , an , F , Pi , q , qa , qap , TACOS
      INTEGER*4 ie , j , k

      q = TACOS(A)
      qa = q
      qap = q
      IF ( q.LE.F ) THEN
         DO j = 1 , 20
            an = 2*j*Pi
            DO k = 1 , 2
               qap = qa
               ie = (-1)**k
               qa = an + ie*q
               IF ( qa.GT.F ) GOTO 100
            ENDDO
         ENDDO
      ENDIF
 100  ARCCOS = qa
      IF ( (qa-F).GT.Pi/2. ) ARCCOS = qap
      END
 
C----------------------------------------------------------------------
C FUNCTION ARCTG
C
C Called by: GOSIA
C
C Purpose: calculates an arctangent in a particular range
C
C Formal parameters:
C      A      - argument
C      F      - range
C      Pi     - Pi must be set to 3.14159... before calling ARCTG
C
C Return value:
C      arctangent(A) within range of F
 
      REAL*8 FUNCTION ARCTG(A,F,Pi)
      IMPLICIT NONE
      REAL*8 A , an , F , Pi , q , qa , qap
      INTEGER*4 ie , j , k

      q = ATAN(A)
      qa = q
      qap = q
      IF ( q.LE.F ) THEN
         DO j = 1 , 40
            an = j*Pi
            DO k = 1 , 2
               qap = qa
               ie = (-1)**k
               qa = an + ie*q
               IF ( qa.GT.F ) GOTO 100
            ENDDO
         ENDDO
      ENDIF
 100  ARCTG = qa
      IF ( (qa-F).GT.Pi/4. ) ARCTG = qap
      END
 
C----------------------------------------------------------------------
C SUBROUTINE LOAD
C
C Called by: FTBM, GOSIA
C Calls:     LSLOOP
C
C Purpose: calculates various parameters, xi, psi which are stored in
C variables XI (common CXI) and PSI (common PCOM).
C
C Uses global variables:
C      CAT    - substates of levels (n_level, J, m)
C      DIPOL  - E1 polarization parameter
C      EMMA   - Controls number of magnetic substates in full coulex calc.
C      EN     - energy of level
C      EP     - bombarding energy
C      ERR    - error flag
C      IPATH  - index of substate in level with same m as substate Irld
C      ISHA   - is half-integer spin
C      ISMAX  - number of substates used
C      IZ     - Z of investigated nucleus
C      IZ1    - Z of not-investated nucleus
C      LAMDA  - list of multipolarities to calculate
C      LAMMAX - number of multipolarities to calculate
C      LDNUM  - number of matrix elements with each multipolarity populating level
C      LEAD   - pair of levels involved in each matrix element
C      LMAX   - ground-state spin + 1
C      LMAXE  - maximum multipolarity needed for calculation
C      LP7    - maximum number of zeta coefficients (45100)
C      LP10   - maximum number of magnetic substates (1200)
C      LZETA  - index into ZETA array for zeta for a given multipolarity
C      MAGA   - number of magnetic substates in approximate calculation
C      MAGEXC - flag: 0 means no magnetic excitations, 1 means with mag. exc.
C      MEMAX  - number of matrix elements
C      NMAX   - number of levels
C      NSTART - index in CAT of first substate associated with a level
C      NSTOP  - index in CAT of last substate associated with a level
C      PSI    - psi coefficients
C      QAPR   - approximate Coulomb amplitudes
C      SPIN   - spin of level
C      VINF   - speed of projectile at infinity
C      XA     - A of investigated nucleus
C      XA1    - A of not-investated nucleus
C      XI     - xi coupling constants
C      ZPOL   - dipole term (GDR excitation)
C
C Formal parameters:
C      Iexp   - Number of experiment
C      Ient   - Flag : 1, 2, 3 (read only)
C      Icg    - Flag : 1 = full coulex, 2 = approximate coulex (read only)
C      Polm   - (read only)
C      Joj    - index of substate (write only)
 
      SUBROUTINE LOAD(Iexp,Ient,Icg,Polm,Joj)
      IMPLICIT NONE
      REAL*8 a1 , a2 , aaz2 , aaz3 , aazz , ah , cpsi , dep , eta , 
     &       etan , Polm , pp1 , pp2
      REAL*8 ppp , rlam , ssqrt , szet , wrt , wrtm , z1 , z2 , zet , 
     &       zsqa
      INTEGER*4 i , i1 , i2 , i3 , Icg , Ient , Iexp , ir , is , ispi , 
     &          ispo
      INTEGER*4 jj , jjj , Joj , la , lam , lam1 , ld , m , m1
      INTEGER*4 m2 , mstop , n , n2 , n3 , nn , nz
      INTEGER*4 LAMDA , LEAD , LDNUM , LAMMAX , MULTI
      COMMON /CLCOM / LAMDA(8) , LEAD(2,1500) , LDNUM(8,100) , LAMMAX , 
     &                MULTI(8)
      INTEGER*4 LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &          LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      REAL*8 EN , SPIN , ACCUR , DIPOL , ZPOL , ACCA
      INTEGER*4 ISO
      COMMON /COEX  / EN(100) , SPIN(100) , ACCUR , DIPOL , ZPOL , 
     &                ACCA , ISO
      INTEGER*4 ISHA
      COMMON /PSPIN / ISHA
      INTEGER*4 MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR
      COMMON /CEXC  / MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR(1500)
      REAL*8 PSI
      COMMON /PCOM  / PSI(1500)
      REAL*8 ZETA
      INTEGER*4 LZETA
      COMMON /CCOUP / ZETA(155600) , LZETA(8)
      INTEGER*4 LMAX
      COMMON /CLM   / LMAX
      REAL*8 CAT
      INTEGER*4 ISMAX
      COMMON /CLCOM8/ CAT(1200,3) , ISMAX
      LOGICAL ERR
      COMMON /CLCOM9/ ERR
      INTEGER*4 NMAX , NDIM , NMAX1
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      REAL*8 XA , XA1 , EP , TLBDG , VINF
      INTEGER*4 NEXPT , IZ , IZ1
      COMMON /CX    / NEXPT , IZ , XA , IZ1(50) , XA1(50) , EP(50) , 
     &                TLBDG(50) , VINF(50)
      INTEGER*4 NSTART , NSTOP
      COMMON /CEXC0 / NSTART(101) , NSTOP(100)
      REAL*8 XI
      COMMON /CXI   / XI(1500)
      REAL*8 EMMA
      INTEGER*4 NCM
      COMMON /CAUX0 / EMMA(100) , NCM
      REAL*8 QAPR
      INTEGER*4 IAPR , ISEX
      COMMON /APRCAT/ QAPR(1500,2,7) , IAPR(1500,2) , ISEX(100)
      INTEGER*4 IPATH , MAGA
      COMMON /PTH   / IPATH(100) , MAGA(100)
      DIMENSION etan(100) , cpsi(8)
      
      LMAX = INT(SPIN(1)+1.1) ! ground-state spin + 1
      IF ( Ient.EQ.1 ) THEN
         ISHA = 0
         ispi = INT(SPIN(1)+.51)
         ispo = INT(SPIN(1)+.49)
         IF ( ispi.NE.ispo ) ISHA = 1 ! Half-integer spin
         z1 = DBLE(ABS(IZ1(Iexp)))
         z2 = DBLE(IZ)
         a1 = XA1(Iexp)
         a2 = XA
         IF ( IZ1(Iexp).GT.0 ) THEN
            ZPOL = DIPOL*EP(Iexp)*a2/(z2*z2*(1.+a1/a2))
         ELSE
            ZPOL = DIPOL*EP(Iexp)*a2/(z2*z2*(1.+a2/a1))
            ah = a1
            a1 = a2
            a2 = ah
         ENDIF

C        Calculate xi and store it in XI in common CXI
C        The value 6.349770 is 197.33/1.44*sqrt(2/931.49).
C        i.e. hbar c / e^2 * sqrt(2 / amu).
         eta = z1*z2*SQRT(a1/EP(Iexp))/6.349770
         DO m = 1 , NMAX
            dep = (1.0+a1/a2)*EN(m)
            zet = dep/EP(Iexp)
            szet = SQRT(1.0-zet)
            etan(m) = eta/szet
         ENDDO
         DO n = 1 , MEMAX
            i1 = LEAD(1,n)
            i2 = LEAD(2,n)
            XI(n) = etan(i1) - etan(i2)
         ENDDO

C        Calculate C_\lambda \over (s Z_1 Z_2)^\lambda 
         aazz = 1./(1.+a1/a2)/z1/z2
         cpsi(1) = 5.169286*aazz

C        Electric excitations up to LMAXE
         IF ( LMAXE.NE.1 ) THEN
            aaz2 = aazz*aazz
            cpsi(2) = 14.359366*aaz2
            IF ( LMAXE.NE.2 ) THEN
               aaz3 = aazz*aaz2
               cpsi(3) = 56.982577*aaz3
               IF ( LMAXE.NE.3 ) THEN
                  aazz = aaz2*aaz2
                  cpsi(4) = 263.812653*aazz
                  IF ( LMAXE.NE.4 ) THEN
                     aaz2 = aaz3*aaz2
                     cpsi(5) = 1332.409500*aaz2
                     IF ( LMAXE.NE.5 ) THEN
                        aazz = aaz3*aaz3
                        cpsi(6) = 7117.691577*aazz
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF

C        Magnetic excitations
         IF ( MAGEXC.NE.0 ) THEN ! If magnetic excitations are required
            aazz = VINF(Iexp)/95.0981942
            cpsi(7) = aazz*cpsi(1)
            IF ( LAMMAX.NE.8 ) cpsi(8) = aazz*cpsi(2)
         ENDIF

C        Calculate psi and store in PSI in common PCOM
         zsqa = z1*SQRT(a1)
         i3 = 1
         ppp = 1. + a1/a2
         DO i1 = 1 , LAMMAX ! For each calculated multipolarity
            lam = LAMDA(i1)
            lam1 = lam
            IF ( lam.GT.6 ) lam1 = lam - 6
            DO n2 = 1 , NMAX ! For each level
               nn = LDNUM(lam,n2) ! Number of matrix elements connected to this level by this multipolarity
               IF ( nn.NE.0 ) THEN
                  n3 = LEAD(1,i3)
                  pp1 = EP(Iexp) - ppp*EN(n3)
                  DO m1 = 1 , nn ! For each matrix element connected to level with this multipolarity
                     m2 = LEAD(2,i3)
                     i2 = i3
                     i3 = i3 + 1
                     pp2 = EP(Iexp) - ppp*EN(m2)
                     PSI(i2) = cpsi(lam)*zsqa*(pp1*pp2)
     &                         **((2.*DBLE(lam1)-1.)/4.)
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
         IF ( Ient.EQ.1 ) RETURN
      ENDIF

C     Initialise NSTART and NSTOP arrays
      DO n = 1 , NMAX ! For each level
         NSTART(n) = 0
         NSTOP(n) = 0
      ENDDO

      is = 1
      NSTART(1) = 1
      DO n = 1 , NMAX ! For each level
         wrt = Polm - EMMA(Iexp)
         wrtm = Polm + EMMA(Iexp)
         IF ( Icg.EQ.2 ) wrt = Polm - DBLE(MAGA(Iexp))
         IF ( Icg.EQ.2 ) wrtm = Polm + DBLE(MAGA(Iexp))
         IF ( wrtm.LT.-SPIN(n) ) THEN
            NSTART(n) = 0 ! Level has no interesting substates
         ELSE
            IF ( ABS(wrt).GT.SPIN(n) ) wrt = -SPIN(n)
            IF ( wrtm.GT.SPIN(n) ) wrtm = SPIN(n)
            mstop = INT(wrtm-wrt+1.01)
            DO i = 1 , mstop ! For each substate
               CAT(is,1) = n               ! Number of level
               CAT(is,2) = SPIN(n)         ! Spin of level
               CAT(is,3) = wrt + DBLE(i-1) ! m quantum number of substate
               IF ( n.EQ.1 .AND. ABS(CAT(is,3)-Polm).LT.1.E-6 ) Joj = is
               is = is + 1
            ENDDO ! Loop on substates i
         ENDIF
         NSTART(n+1) = is ! First substate of level n+1
         NSTOP(n) = is - 1 ! Last substate of level n
      ENDDO ! Loop on levels n

      ISMAX = is - 1 ! ISMAX is the number of substates used
      IF ( ISMAX.LE.LP10 ) THEN ! LP10 is max. number of substates (1200)
         IF ( Ient.EQ.3 ) RETURN
         nz = 0
         DO jj = 1 , 7
            DO jjj = 1 , MEMAX ! For each matrix element
               QAPR(jjj,1,jj) = 0.
               QAPR(jjj,2,jj) = 0.
            ENDDO
         ENDDO

C        Initialise pointers to ZETA array
         DO i = 1 , 8
            LZETA(i) = 0
         ENDDO

         DO i1 = 1 , LAMMAX ! For each multipolarity
            lam = LAMDA(i1)
            IF ( Icg.NE.2 .OR. lam.LE.6 ) THEN
               la = lam
               IF ( lam.GT.6 ) lam = lam - 6
               rlam = DBLE(lam)
               ssqrt = SQRT(2.*rlam+1.)
               LZETA(la) = nz
               ir = 0
 10            ir = ir + 1
               IF ( ir.LE.ISMAX ) THEN
                  n = CAT(ir,1) ! number of level for substate ir
                  IF ( Icg.NE.1 ) THEN
                     IF ( MAGA(Iexp).EQ.0 .AND. ir.NE.IPATH(n) ) GOTO 10
                     IF ( ABS(ir-IPATH(n)).GT.1 ) GOTO 10
                  ENDIF
                  ld = LDNUM(la,n) ! Number of matrix elements connected to this level by this multipolarity
                  IF ( ld.EQ.0 ) THEN
                     ir = ir + NSTOP(n) - NSTART(n)
                  ELSE
                     CALL LSLOOP(ir,n,nz,ld,lam,la,ssqrt,Icg,Iexp)
                  ENDIF
                  GOTO 10
               ENDIF
            ENDIF
         ENDDO ! Loop over multipolarity

C        Make sure we haven't reached the limit of zeta coefficients          
         IF ( nz.GT.LP7 ) THEN ! LP7 is 45100
            WRITE (22,99001) LP7
99001       FORMAT (1x,
     &              'ERROR - NUMBER OF ELEMENTS IN ZETA ARRAY EXCEEDS',
     &              'ZEMAX',5X,'(ZEMAX =',I6,')')
         ELSE
            RETURN
         ENDIF
      ELSE
         WRITE (22,99002) LP10
99002    FORMAT (' ERROR-ISMAX EXCEEDS MAGMAX',5X,'(MAGMAX =',I4,')')
      ENDIF

      ERR = .TRUE. ! Set error flag
      END
 
C----------------------------------------------------------------------
C SUBROUTINE LSLOOP
C
C Called by: LOAD
C Calls:     CODE7, LEADF, WTHREJ
C
C Purpose: calculates the coupling parameter zeta and stores it in the
C          ZETA array starting at the beginning of this array (note that
C          this array has other things in it as well as zeta).
C
C Uses global variables:
C      CAT    - substates of levels (n_level, J, m)
C      ELM    - matrix elements
C      IAPR   - index of initial and final levels for each matrix element
C      IFAC   - spin/parity phase factor
C      ISO    - isotropic flag
C      LP7    - maximum number of zeta coefficients (45100)
C      MAGA   - number of magnetic substates in approximate calculation
C      NSTART - index in CAT of first substate associated with a level
C      NSTOP  - index in CAT of last substate associated with a level
C      PSI    - psi coefficients
C      QAPR   - approximate Coulomb amplitudes
C      SPIN   - spin of level
C      ZETA   - zeta coupling coefficients
C
C Formal parameters:
C      Ir     - index of first substate of level
C      N      - index of level
C      Nz     - index into ZETA array for this multipolarity
C      Ld     - number of matrix elements with this multipolarity
C      Lam    - lambda
C      La     - 1...6 for E1...6 or 7,8 for M1,2
C      Ssqrt  - sqrt(2 * lambda + 1)
C      Icg    - 1 = full coulex, 2 = approximate coulex (read only)
C      Iexp   - experiment number
C
C \zeta_{kn}^{(\lambda n)} = \sqrt{2 \lambda + 1} *
C                            (-1)^{I_n - M_n} *
C                            \three_j{I_n \lambda I_k -M_n \mu M_k} *
C                            \psi_{kn}
C
C For the evaluation of the 3-j symbol, ins = 2 I_n, lam2 = 2 \lambda,
C inr = 2 I_k, jg1 = -2 M_n, jg2 = 2 * \mu, jrmir = 2 * M_k. Note that the
C parameters to WTHREJ are all doubled, so that this routine can cope with
C half-integers.
 
      SUBROUTINE LSLOOP(Ir,N,Nz,Ld,Lam,La,Ssqrt,Icg,Iexp)
      IMPLICIT NONE
      REAL*8 phz , rmir , rmis , Ssqrt , WTHREJ
      INTEGER*4 i2 , i3 , Icg , Iexp , iiex , indx , inqa , inr , 
     &          ins , Ir , is , is1 , is2 , ismin
      INTEGER*4 isplus , jg1 , jg2 , jrmir , La , Lam , lam2 , Ld , 
     &          LEADF , m , MEM , mrange , mt , N , Nz
      REAL*8 EN , SPIN , ACCUR , DIPOL , ZPOL , ACCA
      INTEGER*4 ISO
      COMMON /COEX  / EN(100) , SPIN(100) , ACCUR , DIPOL , ZPOL , 
     &                ACCA , ISO
      REAL*8 PSI
      COMMON /PCOM  / PSI(1500)
      REAL*8 ZETA
      INTEGER*4 LZETA
      COMMON /CCOUP / ZETA(155600) , LZETA(8)
      REAL*8 CAT
      INTEGER*4 ISMAX
      COMMON /CLCOM8/ CAT(1200,3) , ISMAX
      INTEGER*4 NSTART , NSTOP
      COMMON /CEXC0 / NSTART(101) , NSTOP(100)
      REAL*8 QAPR
      INTEGER*4 IAPR , ISEX
      COMMON /APRCAT/ QAPR(1500,2,7) , IAPR(1500,2) , ISEX(100)
      INTEGER*4 IPATH , MAGA
      COMMON /PTH   / IPATH(100) , MAGA(100)
      INTEGER*4 LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &          LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      REAL*8 ELM , ELMU , ELML , SA
      COMMON /COMME / ELM(1500) , ELMU(1500) , ELML(1500) , SA(1500)
      INTEGER*4 IFAC
      COMMON /CLCOM0/ IFAC(100)
      
      lam2 = 2*Lam
      inr = CAT(Ir,2)*2. ! 2 * Spin of substate Ir
      rmir = CAT(Ir,3)   ! m quantum number of substate Ir
      jrmir = 2.*rmir
      DO i2 = 1 , Ld
         m = LEADF(N,i2,La) ! Index of final level
         indx = MEM(N,m,La) ! Index of matrix element
         IAPR(indx,1) = N   ! Index of initial level
         IAPR(indx,2) = m   ! Index of final level
         ismin = 0
         ins = SPIN(m)*2.
         is1 = NSTART(m) ! Index of first substate of level m
         IF ( is1.NE.0 ) THEN
            isplus = INT(rmir-CAT(is1,3)) - Lam
            IF ( isplus.LT.0 ) THEN
               ismin = isplus
               isplus = 0
            ENDIF
            is2 = is1 + isplus - 1
            mrange = 2*Lam + 1 + ismin
            IF ( is2+mrange.GT.NSTOP(m) ) mrange = NSTOP(m) - is2
            IF ( mrange.GT.0 ) THEN
               DO i3 = 1 , mrange
                  is = is2 + i3
                  rmis = CAT(is,3) ! m quantum number of substate is
                  IF ( ISO.NE.0 .OR. rmis.LE..1 .OR. rmir.LE..1 ) THEN
                     jg1 = -rmis*2.
                     jg2 = (rmis-rmir)*2.
                     IF ( Icg.NE.2 .OR. ABS(jg2).LE.2*MAGA(Iexp) ) THEN
                        IF ( La.LE.6 .OR. jg2.NE.0 ) THEN
                           Nz = Nz + 1
                           IF ( Nz.LE.LP7 ) THEN
                              iiex = (ins+jg1)/2
                              phz = (-1.0)**iiex
                              ZETA(Nz) = phz*PSI(indx) ! This is really zeta
     &                           *Ssqrt*WTHREJ(ins,lam2,inr,jg1,jg2,
     &                           jrmir)
                              IF ( Icg.NE.1 ) THEN
                                 mt = CAT(is,1) ! level number of substate is
                                 CALL CODE7(Ir,is,N,mt,inqa,indx)
                                 IF ( ABS(ELM(indx)).LT.1.E-6 )
     &                                ELM(indx) = 1.E-6
                                 IF ( inqa.NE.-1 ) THEN
                                    QAPR(indx,1,inqa) = ZETA(Nz)
     &                                 *ELM(indx)
                                    IF ( ISO.EQ.0 .AND. inqa.EQ.1 )
     &                                 QAPR(indx,1,7) = QAPR(indx,1,1)
     &                                 *IFAC(m)
                                 ENDIF
                              ENDIF
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDIF ! If isotropic or rmis < 1 or rmir < 1
               ENDDO ! Loop on substates
            ENDIF ! If range of substates is greater than 0
         ENDIF ! If there are substates
      ENDDO ! Loop on matrix elements
      END
 
C----------------------------------------------------------------------
C FUNCTION LEADF
C
C Called by: LAIAMP, LSLOOP, NEWLV, SEQ
C
C Uses global variables:
C      LDNUM  - number of matrix elements with each multipolarity populating levels
C      LEAD   - pair of levels involved in each matrix element
C      MULTI  - number of matrix elements with a given multipolarity
C
C Formal parameters:
C      N1     - index of initial level
C      N2     - index of matrix element for given level and multipolarity
C      N3     - multipolarity
C
C Purpose: calculate the level number for the final level associated with the
C      matrix element index N2, initial level index N1 and multipolarity N3.
 
      INTEGER*4 FUNCTION LEADF(N1,N2,N3)
      IMPLICIT NONE
      INTEGER*4 k , lsum , N1 , n1m , N2 , N3 , n3m
      INTEGER*4 LAMDA , LEAD , LDNUM , LAMMAX , MULTI
      COMMON /CLCOM / LAMDA(8) , LEAD(2,1500) , LDNUM(8,100) , LAMMAX , 
     &                MULTI(8)
      lsum = 0
      n3m = N3 - 1
      IF ( n3m.NE.0 ) THEN
         DO k = 1 , n3m ! Loop over multipolarities lower than one required
            lsum = lsum + MULTI(k)
         ENDDO
      ENDIF

C     lsum now points to start of the multipolarity N3
       
      n1m = N1 - 1
      IF ( n1m.NE.0 ) THEN
         DO k = 1 , n1m ! Loop over levels below the selected one
            lsum = lsum + LDNUM(N3,k)
         ENDDO
      ENDIF

C     lsum now points to start of level N1 for multipolarity N3
      
      n1m = lsum + N2

C     n1m now points to the appropriate matrix element N2 for initial level N1
C     and multipolarity N3
      
      LEADF = LEAD(2,n1m) ! Get the final level for this matrix element
      END
 
C----------------------------------------------------------------------
C FUNCTION MEM
C
C Called by: SEQ, NEWLV
C
C Purpose: calculates an index to a matrix element given two level indices
C      and the multipolarity.
C
C Uses global variables:
C      LDNUM  - number of matrix elements with each multipolarity populating level
C      LEAD   - pair of levels involved in each matrix element
C      MULTI  - number of matrix elements having a given multipolarity
C
C Formal parameters:
C      N1     - level number for first level
C      N2     - level number for second level
C      N3     - multipolarity
C
C Return value:
C      Index of matrix element
 
      INTEGER*4 FUNCTION MEM(N1,N2,N3)
      IMPLICIT NONE
      INTEGER*4 k , msum , N1 , n1m , N2 , N3 , n3m
      INTEGER*4 LAMDA , LEAD , LDNUM , LAMMAX , MULTI
      COMMON /CLCOM / LAMDA(8) , LEAD(2,1500) , LDNUM(8,100) , LAMMAX , 
     &                MULTI(8)
      msum = 0
      IF ( N3.NE.1 ) THEN
         n3m = N3 - 1
         DO k = 1 , n3m ! For each multipolarity up to one below the one we want
            msum = msum + MULTI(k) ! Add the number of matrix elements for that multipolarity
         ENDDO
      ENDIF

C     msum is now an index to the start of the matrix elements for the chosen multipolarity

      n1m = N1 - 1
      IF ( n1m.NE.0 ) THEN
         DO k = 1 , n1m ! For each level up to one below the one we want
            msum = msum + LDNUM(N3,k) ! Add the number of matrix elements for that level and multipolarity
         ENDDO
      ENDIF

C     msum is now an index to the start of the matrix elements for the appropriate multipolarity and level

      n1m = msum + 1
      n3m = n1m + LDNUM(N3,N1)
      DO k = n1m , n3m ! Loop over matrix elements associated with that level and multipolarity
         msum = msum + 1
         IF ( LEAD(2,k).EQ.N2 ) GOTO 100 ! If it is the right one goto 100
      ENDDO

 100  MEM = msum ! MEM is now the index to the matrix element we want
      END
 
C----------------------------------------------------------------------
C SUBROUTINE CMLAB
C
C Called by: GOSIA
C Calls:     TASIN
C
C Purpose: calculate for center of mass frame
C
C Uses global variables:
C      BETAR  - recoil beta
C      DSIGS  - dsigma for each experiment
C      EN     - level energies
C      EP     - bombarding energy
C      EPS    - epsilon
C      EROOT  - sqrt(epsilon^2 - 1)
C      ERR    - error flag
C      IPRM   - printing flags (see suboption PRT of OP,CONT)
C      ISKIN  - kinematic flag (0,1)
C      IZ     - Z of investigated nucleus
C      IZ1    - Z of not-investated nucleus
C      NCM    - calculate kinematics assuming this state for final state (default = 2)
C      NEXPT  - number of experiments
C      NMAX   - number of level energies
C      TETACM - theta of particle detector in center of mass frame
C      TLBDG  - theta of particle detector in lab frame (in degrees)
C      VINF   - speed of projectile at infinity
C      XA     - A of investigated nucleus
C      XA1    - A of not-investated nucleus
C      TREP   - theta of recoiling nucleus (in radians)
C
C Formal parameters:
C      Ii     - experiment number (or zero for all experiments)
C      Dsig   - dsigma
C      Tetrn  - theta of recoiling nucleus in lab frame (in radians)

      SUBROUTINE CMLAB(Ii,Dsig,Tetrn)
      IMPLICIT NONE
      REAL*8 a1 , a2 , ared , d2a , dista , dists , Dsig , emax , epmin
      REAL*8 r3 , TASIN , tau , taup , tcmdg , tcmrad , Tetrn , tlbrad ,
     &       tmxdg , z1 , z2 , zcmdg , zcmrad , zlbrad
      INTEGER*4 iflaa , Ii , lexp , lexp0 , lexp1 , n
      LOGICAL ERR
      COMMON /CLCOM9/ ERR
      INTEGER*4 ISKIN
      COMMON /SECK  / ISKIN(50)
      INTEGER*4 IPRM
      COMMON /PRT   / IPRM(20)
      REAL*8 TETACM, TREP, DSIGS
      COMMON /TCM   / TETACM(50) , TREP(50) , DSIGS(50)
      REAL*8 BETAR
      COMMON /BREC  / BETAR(50)
      REAL*8 EMMA
      INTEGER*4 NCM
      COMMON /CAUX0 / EMMA(100) , NCM
      INTEGER*4 NMAX , NDIM , NMAX1
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      REAL*8 XA , XA1 , EP , TLBDG , VINF
      INTEGER*4 NEXPT , IZ , IZ1
      COMMON /CX    / NEXPT , IZ , XA , IZ1(50) , XA1(50) , EP(50) , 
     &                TLBDG(50) , VINF(50)
      REAL*8 EPS, EROOT, FIEX
      INTEGER*4 IEXP, IAXS
      COMMON /KIN   / EPS(50) , EROOT(50) , FIEX(50,2) , IEXP , IAXS(50)
      REAL*8 EN , SPIN , ACCUR , DIPOL , ZPOL , ACCA
      INTEGER*4 ISO
      COMMON /COEX  / EN(100) , SPIN(100) , ACCUR , DIPOL , ZPOL , 
     &                ACCA , ISO
      DATA r3/0./

      lexp0 = 1
      lexp1 = NEXPT
      IF ( Ii.NE.0 ) lexp0 = Ii
      IF ( Ii.NE.0 ) lexp1 = Ii
      DO lexp = lexp0 , lexp1 ! For each experiment
         iflaa = 0
         IF ( TLBDG(lexp).LT.0 ) iflaa = 1
         IF ( IPRM(1).EQ.1 ) THEN
            IF ( Ii.EQ.0 .AND. IPRM(10).EQ.1 ) WRITE (22,99001) lexp
99001       FORMAT (1X,///10X,'** EXPERIMENT',1X,1I2,1X,'**'//)
         ENDIF
         TLBDG(lexp) = ABS(TLBDG(lexp))
         a1 = XA1(lexp)
         IF ( IZ1(lexp).LT.0 ) a1 = XA
         a2 = XA
         IF ( IZ1(lexp).LT.0 ) a2 = XA1(lexp)
         z1 = DBLE(ABS(IZ1(lexp)))
         z2 = DBLE(IZ)
         IF ( IPRM(1).EQ.1 ) THEN
            IF ( IZ1(lexp).LT.0 .AND. (Ii.EQ.0 .AND. IPRM(10).EQ.1) )
     &           WRITE (22,99002) IZ , XA , ABS(IZ1(lexp)) , XA1(lexp)
99002       FORMAT (5X,'PROJECTILE EXCITATION OF(',1I3,',',1F7.3,
     &              ') ON(',1I3,',',1F7.3,')')
            IF ( IZ1(lexp).GT.0 .AND. (Ii.EQ.0 .AND. IPRM(10).EQ.1) )
     &           WRITE (22,99003) IZ , XA , IZ1(lexp) , XA1(lexp)
99003       FORMAT (5X,'TARGET EXCITATION OF(',1I3,',',1F7.3,') BY(',
     &              1I3,',',1F7.3,')')
         ENDIF
C
C        dists is Cline's estimate of the maximum safe bombarding energy
         dists = 1.44*(a1+a2)*z1*z2/((a1**.33333+a2**.33333)*1.25+5.)/a2
C        dista is 0.05 * distance of closest approach for head-on collisions
         dista = 0.0719949*(1.0+a1/a2)*z1*z2/EP(lexp)
C        d2a is the distance of closest approach for head-on collisions in fm
C        q^2/4/pi/epsilon_0 * (1+a1/a2) * Z1 * Z2 / Ep. For Ep in MeV and d2a
C        in fm, q^2/4/pi/epsilon_0 = 1.44
         d2a = 20.0*dista ! = 1.44 * (1.0+a1/a2)*z1*z2/EP(lexp)
C        VINF is the initial velocity of the incoming projectile (at infinity)
C        VINF = sqrt(2 * EP / 931.494028 * A1) : 931.494028 = 1 AMU
         VINF(lexp) = 0.0463365*SQRT(EP(lexp)/a1)

C        If IPRM(1) we want extra printout
         IF ( IPRM(1).EQ.1 ) THEN
            IF ( Ii.EQ.0 .AND. IPRM(10).EQ.1 ) WRITE (22,99004) EP(lexp)
     &           , VINF(lexp)
99004       FORMAT (5X,'ENERGY',1X,1F10.3,1X,'MEV',5X,'BETA',1X,1E14.6)
            IF ( EP(lexp).GT.dists .AND. (Ii.EQ.0 .AND. IPRM(10).EQ.1) )
     &           WRITE (22,99005) (EP(lexp)/dists-1.)*100.
99005       FORMAT (5X,'***** ','BE CAREFUL-ACCORDING',
     &              ' TO D.CLINE BOMBARDING ENERGY',1X,1F6.2,1X,'PC',1X,
     &              ' TOO HIGH FOR HEAD-ON COLLISIONS! *****')
            IF ( Ii.EQ.0 .AND. IPRM(10).EQ.1 ) WRITE (22,99006) d2a
99006       FORMAT (5X,
     &             'DISTANCE OF CLOSEST APPROACH FOR HEAD-ON COLLISIONS'
     &             ,1X,1F10.4,1X,'FM')
         ENDIF

C        Final kinetic energy \v{E} = E_P - \Delta E (1 + m_P / m_T)
C        Here we set ared = (1 + m_P / m_T)
C        The maximum excitation energy corresponds to \v{E} = 0, so
C        \DeltaE = E_P \over {1 + m_P / m_T) = E_P/ared
C        We check that there are no states defined which are higher than this.

         tlbrad = TLBDG(lexp)/57.2957795 ! Theta of detector to radians
         ared = 1.0 + a1/a2 ! reduced mass
         emax = EP(lexp)/ared ! Maximum excitation energy
         DO n = 1 , NMAX ! For each level
            IF ( EN(n).GT.emax ) GOTO 50 ! Give error if energy of state too high
         ENDDO

C        Gosia calculates assuming the kinematics for all states are approximately
C        those corresponding to the state NCM (by default NCM = 2 : the first excited
C        state). So for this energy we calculate the \v{E} and store it in epmin.
C        We also calculate tau defined as (a1/a2)*sqrt(E_P / \v{E}) for this value
C        of \v{E}.
C        A value of tau less than 1 corresponds to normal kinematics, so the full
C        range of theta in the centre of mass system corresponds to the full range
C        in the lab system. However, for tau greater than 1 (i.e. inverse kinematics)
C        there are two possible values for the lab angle for a given centre of mass
C        angle and there is a maximum lab angle, which can be attained: tmxdg given
C        by SIN(tmxdg) = 1 / tau.
         epmin = EP(lexp) - EN(NCM)*ared
         taup = SQRT(EP(lexp)/epmin)
         tau = taup*a1/a2
         IF ( tau.LE.1.0 ) GOTO 100 ! No limit on scattering angle
         tmxdg = TASIN(1.0/tau)*57.2957795 ! Maximum lab angle in degrees
         IF ( tmxdg.GE.TLBDG(lexp) ) GOTO 100 ! Within limit of scattering angle

         WRITE (22,99007) tmxdg , lexp
99007    FORMAT (1X,'ERROR- MAXIMUM SCATTERING ANGLE IS ',F7.2,
     &           ' DEGREES',' FOR EXPERIMENT ',1I2)
         GOTO 200 ! Error

 50      WRITE (22,99008) emax , lexp
99008    FORMAT (1X,'ERROR- MAXIMUM EXCITATION ENERGY IS ',F8.4,' MEV',
     &           ' FOR EXPERIMENT ',1I2)
         GOTO 200 ! Error

C        Calculate centre of mass angle
 100     tcmrad = tlbrad + TASIN(tau*SIN(tlbrad)) ! In radians
         tcmdg = tcmrad*57.2957795 ! and in degrees

C        In inverse kinematics, for a given lab angle, there are two solutions
C        for the centre of mass angle.
         IF ( tau.GT.1.0 ) THEN ! Inverse kinematics
            IF ( IPRM(1).EQ.1 ) THEN
               IF ( Ii.EQ.0 .AND. IPRM(10).EQ.1 ) WRITE (22,99009)
     &              tcmdg , lexp
99009          FORMAT (5X,'SECOND POSSIBLE CM SCATTERING ANGLE IS',F7.2,
     &                 ' DEGREES FOR EXPERIMENT ',1I2)
            ENDIF
            IF ( ISKIN(lexp).NE.1 ) THEN ! If ISKIN is set, take the second solution
               tcmdg = 180. + 2.*TLBDG(lexp) - tcmdg
               tcmrad = tcmdg/57.2957795
            ENDIF
         ENDIF

C        EPS is "epsilon" the eccentricity parameter.
         EPS(lexp) = 1./SIN(tcmrad/2.)
         TETACM(lexp) = tcmrad
         IF ( IPRM(1).EQ.1 ) THEN
            IF ( Ii.EQ.0 .AND. IPRM(10).EQ.1 ) WRITE (22,99010) tcmdg , 
     &           EPS(lexp)
99010       FORMAT (5X,'CM SCATTERING ANGLE',1X,1F10.3,1X,'DEG',5X,
     &              'EPSILON',1X,1F10.4)
         ENDIF

C        If Z1 is negative, we are interested in target excitations, but if it
C        is positive, we want the projectile excitation, so calculate the lab
C        recoil energy of appropriate particle and store it in BETAR (we will
C        convert this to beta of the recoil later)
         IF ( IZ1(lexp).GT.0 ) BETAR(lexp) = a1*a2/(a1+a2)
     &        **2*(1.+taup*taup-2.*taup*COS(tcmrad))*epmin
         IF ( IZ1(lexp).LT.0 ) BETAR(lexp) = (a2/(a1+a2))
     &        **2*(1.+tau*tau+2.*tau*COS(tcmrad))*epmin

C        More additional printout
         IF ( IPRM(1).EQ.1 ) THEN
            IF ( Ii.EQ.0 .AND. IPRM(10).EQ.1 ) WRITE (22,99011)
     &           BETAR(lexp)
99011       FORMAT (5X,'RECOIL ENERGY(MEV)',2X,1F10.4)
         ENDIF

C        This is the beta of the recoiling particle of interest (target or projectile
C        depending on sign of Z1, which is used as a flag)
         BETAR(lexp) = .0463365*SQRT(BETAR(lexp)/XA) ! 0.0463365 = sqrt(2/931.494028)
         IF ( IPRM(1).EQ.1 ) THEN
            IF ( Ii.EQ.0 .AND. IPRM(10).EQ.1 ) WRITE (22,99012)
     &           BETAR(lexp)
99012       FORMAT (5X,'RECOIL BETA',2X,1E14.6)
            IF ( Ii.EQ.0 .AND. IPRM(10).EQ.1 ) WRITE (22,99013) EP(lexp)
     &           /(dists*.5*(1.+EPS(lexp)))
99013       FORMAT (5X,'BOMBARDING ENERGY=',1F10.3,1X,
     &              'OF SAFE BOMBARDING ENERGY AT THIS ANGLE')
         ENDIF

C        iflaa = 0 when projectile detected, = 1 when target detected
C        r3 is the Jacobian dOmega/domega
         IF ( iflaa.NE.1 ) THEN ! Projectile detected
            IF ( ABS(tcmdg-180.).LT.1.E-5 ) THEN
               r3 = (1.-tau)**2
            ELSE
               r3 = SIN(tlbrad)/SIN(tcmrad)
               r3 = r3*r3*ABS(COS(tcmrad-tlbrad))
               r3 = 1./r3
            ENDIF
         ENDIF

C        Calculate the values for the target. In the centre of mass system, the
C        target and projectile angles differ by 180 degrees
         zcmdg = 180. - tcmdg ! Target angle in degrees in cm system
         zcmrad = zcmdg/57.2957795 ! and in radians
         zlbrad = ATAN(SIN(zcmrad)/(COS(zcmrad)+taup)) ! target theta in lab (radians)

C        iflaa = 0 when projectile detected, = 1 when target detected
C        r3 is the Jacobian dOmega/domega
         IF ( iflaa.NE.0 ) THEN ! Target detected, but theta is for projectile!
            IF ( ABS(tcmdg-180.).LT.1.E-5 ) THEN
               r3 = (1.+taup)**2
               TLBDG(lexp) = 0.
            ELSE
               r3 = SIN(zlbrad)/SIN(zcmrad)
               r3 = r3*r3
               r3 = r3*ABS(COS(zcmrad-zlbrad))
               r3 = 1./r3
               TLBDG(lexp) = zlbrad*57.2955795
            ENDIF
         ENDIF

C        Now calculate dsigma
         Dsig = 250.*r3*SQRT(EP(lexp)/(EP(lexp)-ared*EN(NCM)))
     &          *dista*dista*(EPS(lexp))**4
         EROOT(lexp) = SQRT(EPS(lexp)*EPS(lexp)-1.)
         DSIGS(lexp) = Dsig
         Tetrn = zlbrad
         IF ( IZ1(lexp).LT.0. ) Tetrn = tlbrad
         TREP(lexp) = Tetrn
      ENDDO ! Loop over experiments lexp

      IPRM(10) = 0 ! Turn off printing so we don't write things twice
      RETURN

C     An error has occured, so set error flag and return
 200  ERR = .TRUE. ! Set error flag
      END
 
C----------------------------------------------------------------------
C SUBROUTINE QE
C
C Called by: SNAKE
C
C Purpose: calculate Qe values
C
C Formal parameters:
C      C     - cosh(omega) + epsilon
C      D     - sqrt(epsilon^2 - 1) * sinh(omega)
C      B2    - B^2 = (epsilon * cosh(omega) + 1)^2
C      C2    - C^2
C      D2    - D^2
C      B4    - B^4
C      B6    - B^6
C      D3    - D^3
C      B8    - B^8
C      C4    - C^4
C      D4    - D^4
C      B10   - B^10
C      D5    - D^5
C      B12   - B^12
C      D6    - D^6
C      Lmda  - lambda
C      Pol   - E1 polarisation factor to allow for GDR excitation
C      Cq    - array where the results are returned
C
C We used different formulae depending on lambda (see the table of electric
C collision functions in the gosia manual).
C
C Note that we multiply by Pol, which is the E1 dipole correction factor for
C the GDR in the case of the quadrupole. This is correct. It isn't really a
C quadrupole, but the operator has the same form as the E2 operator, so it
C is easiest to consider the effect as a correction to this operator.
C
C Lmda = lambda (1 = E1, 2 = E2... 6 = E6)
 
      SUBROUTINE QE(C,D,B2,C2,D2,B4,B6,D3,B8,C4,D4,B10,D5,B12,D6,Lmda,
     &              Pol,Cq)
      IMPLICIT NONE
      REAL*8 B10 , B12 , B2 , B4 , B6 , B8 , C , C2 , C4 , Cq , D , D2 ,
     &       D3 , D4 , D5 , D6 , Pol
      INTEGER*4 Lmda
      DIMENSION Cq(7)

      IF ( Lmda.EQ.2 ) THEN ! E2
         Cq(1) = 0.75*(2.0*C2-D2)/B4*Pol
         Cq(2) = -1.83711730*C*D/B4*Pol
         Cq(3) = -0.91855865*D2/B4*Pol
         RETURN
      ELSEIF ( Lmda.EQ.3 ) THEN ! E3
         Cq(1) = 1.875*C*(2.0*C2-3.0*D2)/B6
         Cq(2) = -1.62379763*(4.0*C2-D2)*D/B6
         Cq(3) = -5.13489890*C*D2/B6
         Cq(4) = 2.09631373*D3/B6
         RETURN
      ELSEIF ( Lmda.EQ.4 ) THEN ! E4
         Cq(1) = 1.09375000*(8.0*C4-24.0*C2*D2+3.0*D4)/B8
         Cq(2) = -4.89139867*C*(4.0*C2-3.0*D2)*D/B8
         Cq(3) = -3.45874113*(6.0*C2-D2)*D2/B8
         Cq(4) = 12.9414244*C*D3/B8
         Cq(5) = 4.57548440*D4/B8
         RETURN
      ELSEIF ( Lmda.EQ.5 ) THEN ! E5
         Cq(1) = 1.230468*C*(-14.*C2*(9.*D2+B2)+30.*B4)/B10
         Cq(2) = -1.347911*D*(35.*C2*(-3.*D2+B2)+5.*B4)/B10
         Cq(3) = -35.662372*D2*C*(-3.*D2+2.*B2)/B10
         Cq(4) = 7.279552*D3*(9.*C2-B2)/B10
         Cq(5) = 30.884521*D4*C/B10
         Cq(6) = -9.766543*D5/B10
         RETURN
      ELSEIF ( Lmda.EQ.6 ) THEN ! E6
         Cq(1) = 2.707031*(21.*C2*(-C2*(11.*D2+4.*B2)+5.*B4)-5.*B6)/B12
         Cq(2) = -17.543567*D*C*(3.*C2*(-11.*D2+B2)+5.*B4)/B12
         Cq(3) = -13.869408*D2*(3.*C2*(-11.*D2+5.*B2)+B4)/B12
         Cq(4) = 27.738815*D3*C*(-11.*D2+8.*B2)/B12
         Cq(5) = 15.193177*D4*(11.*C2-B2)/B12
         Cq(6) = -71.262308*D5*C/B12
         Cq(7) = -20.571656*D6/B12
         GOTO 99999
      ENDIF ! E1
      Cq(1) = 0.5*C/B2
      Cq(2) = -0.35355339*D/B2
99999 END
 
C----------------------------------------------------------------------
C SUBROUTINE QM
C
C Called by: SNAKE
C
C Purpose: calculate Qm values
C
C Formal parameters:
C      C     - cosh(omega) + epsilon
C      D     - sqrt(epsilon^2 - 1) * sinh(omega)
C      B2    - B^2 = (epsilon * cosh(omega) + 1)^2
C      B4    - B^4
C      Ert   - sqrt(epsilon^2 -1)
C      Lmda  - lambda
C      Cq    - array where the results are returned
C
C We used different formulae depending on lambda (see the table of magnetic
C collision functions in the gosia manual).
C
C Lmda = lambda + 6. i.e. Lmda = 7 is M1, Lmda = 8 is M2.
 
 
      SUBROUTINE QM(C,D,B2,B4,Ert,Lmda,Cq)
      IMPLICIT NONE
      REAL*8 B2 , B4 , C , Cq , D , Ert
      INTEGER*4 Lmda
      DIMENSION Cq(7)

      IF ( Lmda.EQ.8 ) THEN
         Cq(1) = -.9185586536*C*Ert/B4
         Cq(2) = -Cq(1)*D/C
         GOTO 99999
      ENDIF
      Cq(1) = -.3535533905*Ert/B2
99999 END
 
C----------------------------------------------------------------------
C SUBROUTINE SNAKE
C
C Called by: FTBM, GOSIA
C Calls:     QE, QM, QRANGE
C
C Purpose: evaluate and store the dimensionless collision functions Qe and Qm.
C
C Uses global variables:
C      CH     - table of cosh values
C      EPS    - epsilon
C      EROOT  - sqrt(epsilon^2 -1)
C      LOCQ   - location of collision function in ZETA array
C      LP7    - start of collision functions in ZETA (45100)
C      SH     - table of sinh values
C      ZETA   - various coefficients (here the collision functions)
C
C Formal parameters:
C      Nexp   - experiment number
C      Zpol   - dipole term (GDR excitation)
C
C The function QE is used to calculate Qe and QM to calculate Qm, but first
C we call QRANGE to determine the range over which we need to calculate them.
C
C The results are stored in the ZETA array, but not starting from the
C beginning, which is where zeta itself is written, but from ZETA(LP7).
C
C LOCQ (in ALLC) is used as an index to these values.
C
C EROOT is set in CMLAB to \sqrt(\epsilon^2 - 1).
C
C Note that when we call QE and QM that lmda = 1...6 for E1...6 and 7,8 for
C M1, M2.
 
      SUBROUTINE SNAKE(Nexp,Zpol)
      IMPLICIT NONE
      REAL*8 b10 , b12 , b2 , b4 , b6 , b8 , c , c2 , c4 , c6 , 
     &       chi , cq , d , d2 , d3 , d4 , d5 , d6
      REAL*8 ert , pol , shi , Zpol
      INTEGER*4 ibm , icm , icnt , idm , irl , j , k , lloc , lmd , 
     &          lmda , mimx , Nexp , nind , nlm
      DIMENSION lloc(8) , cq(7) , irl(8)
      REAL*8 EPS, EROOT, FIEX
      INTEGER*4 IEXP, IAXS
      COMMON /KIN   / EPS(50) , EROOT(50) , FIEX(50,2) , IEXP , IAXS(50)
      REAL*8 ZETA
      INTEGER*4 LZETA
      COMMON /CCOUP / ZETA(155600) , LZETA(8)
      INTEGER*4 LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &          LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      INTEGER*4 LOCQ
      COMMON /ALLC  / LOCQ(8,7)
      REAL*8 CH , SH
      COMMON /HIPER / SH(365) , CH(365)
      
      icnt = 0
 100  icnt = icnt + 1

C     Calculate range over which we will want Qe and Qm
      CALL QRANGE(icnt,nlm,lloc,ibm,icm,idm,irl)
      IF ( nlm.EQ.0 ) RETURN

C     Calculate some parameters, which we will pass to QE or QM
      chi = CH(icnt) ! \cosh(\omega)
      shi = SH(icnt) ! \sinh(\omega)
      b2 = EPS(Nexp)*chi + 1.
      pol = 1. - Zpol/b2 ! E1 polarisation term
      b2 = b2*b2 ! b^2 = (\epsilon \cosh(\omega) + 1)^2
      IF ( ibm.NE.2 ) THEN
         b4 = b2*b2
         IF ( ibm.NE.4 ) THEN
            b6 = b4*b2
            IF ( ibm.NE.6 ) THEN
               b8 = b4*b4
               IF ( ibm.NE.8 ) THEN
                  b10 = b6*b4
                  IF ( ibm.NE.10 ) b12 = b6*b6
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      IF ( icm.NE.0 ) THEN
         c = chi + EPS(Nexp) ! c = \cosh(\omega) + \epsilon
         IF ( icm.NE.1 ) THEN
            c2 = c*c
            IF ( icm.NE.2 ) THEN
               c4 = c2*c2
               IF ( icm.NE.4 ) c6 = c2*c4
            ENDIF
         ENDIF
      ENDIF
      IF ( idm.NE.0 ) THEN
         d = EROOT(Nexp)*shi ! d = \sinh(\omega) * \sqrt(epsilon^2 - 1)
         IF ( idm.NE.1 ) THEN
            d2 = d*d
            IF ( idm.NE.2 ) THEN
               d3 = d*d2
               IF ( idm.NE.3 ) THEN
                  d4 = d2*d2
                  IF ( idm.NE.4 ) THEN
                     d5 = d3*d2
                     IF ( idm.NE.5 ) d6 = d3*d3
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      DO j = 1 , nlm
         lmda = lloc(j)
         IF ( lmda.GT.6 ) THEN
            lmd = lmda
            lmda = lmda - 6
            ert = EROOT(Nexp)
            CALL QM(c,d,b2,b4,ert,lmda,cq)
            mimx = lmda
            DO k = 1 , mimx
               nind = LOCQ(lmd,k) + icnt
               ZETA(nind+LP7) = cq(k) ! These are the collision functions
            ENDDO
         ELSE
            CALL QE(c,d,b2,c2,d2,b4,b6,d3,b8,c4,d4,b10,d5,b12,d6,lmda,
     &              pol,cq)
            mimx = lmda + 1
            DO k = 1 , mimx
               nind = LOCQ(lmda,k) + icnt
               ZETA(nind+LP7) = cq(k) ! These are the collision functions
            ENDDO
         ENDIF
      ENDDO
      IF ( icnt.EQ.365) THEN
        STOP 'Sorry, I can only do 365 steps in omega. You need more!'
      ENDIF
      GOTO 100
      END
 
C----------------------------------------------------------------------
C SUBROUTINE SPITQ
C
C Called by: INTG
C Calls:
C
C Purpose: write out collision functiosn to a file
C
C Uses global variables:
C
C Formal parameters:
C      Ixpt   - experiment number
C      Mult   - multipolarity

      SUBROUTINE SPITQ(Ixpt,Mult)
      IMPLICIT NONE
      REAL*8 ZETA
      INTEGER*4 LZETA
      COMMON /CCOUP / ZETA(155600) , LZETA(8)
      REAL*8 EPS, EROOT, FIEX
      INTEGER*4 IEXP, IAXS
      COMMON /KIN   / EPS(50) , EROOT(50) , FIEX(50,2) , IEXP , IAXS(50)
      INTEGER*4 LOCQ
      COMMON /ALLC  / LOCQ(8,7)
      INTEGER*4 IRA , MAXLA
      COMMON /RNG   / IRA(8) , MAXLA
      INTEGER*4 LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &          LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      INTEGER*4 mimmex
      INTEGER*4 Ixpt , Mult
      INTEGER*4 ibm , icm , icnt , idm , irl  , k , lloc , 
     &           nind , nlm
      DIMENSION lloc(8) , irl(8)
      REAL*8 collfunc , w0

c     Passed in the experiment number (Ixpt) and the multipolarity (Mult).
c     I am doing only for electric for now.
c     This routine is to read the electric collision functions from the ZETA array and print
c     them to the standard output file (23)  

C     Borrowed this line from subroutine below.
      w0 = IRA(MAXLA) * .03 + .03 ! Maximum omega to calculate for (steps of 0.03)

c     Look at how LAIAMP uses the ZETA array.
c     Are these the collision functions?

      WRITE(99,72072) Ixpt , Mult
72072 FORMAT('The *stored* collision functions ',
     &       'for experiment ' , i2 , ' and mult E' , i1)
      WRITE(99,71972)
71972 FORMAT('   |--------mu=0----------|--------mu=2---------',
     &       '|--------mu=2...')
      WRITE(99,71672)
71672 FORMAT('   |icnt nind  Qe         |icnt nind  Qe        ',
     &       '|icnt nind  Qe...')
      icnt = 0
 100  icnt = icnt + 1


      mimmex = Mult + 1
c     WRITE(99,71772) mimmex

C     With QRANGE I think I am checking that icnt is still in range to index the QE values
c     for this multipolarity and experiment.  There must be only one epsilon 
c     value per experiment.
      CALL QRANGE(icnt , nlm , lloc , ibm , icm , idm , irl)
c     WRITE(99,71772) mimmex
c71772 FORMAT(i4)
c     test to here

      IF ( nlm.EQ.0 ) RETURN

      DO k = 1 , mimmex
         nind = LOCQ(Mult,k) + icnt
         collfunc = ZETA(nind+LP7)  ! These are the collision functions
         WRITE(99,71572) icnt , nind , collfunc
71572    FORMAT(1X,I4,1X,I6,1X,D11.4,$)
      ENDDO
      WRITE(99,71872)
71872 FORMAT('')
      GOTO 100 

      END
 
C----------------------------------------------------------------------
C SUBROUTINE FHIP
C
C Called by: GOSIA
C
C Purpose: generates a table of the hyperbolic funcions sinh and cosh for
C later use. Note that these are in steps of \Delta\omega = 0.03. These are
C stored in the common block HIPER.
C
C Uses global variables:
C      CH     - table of cosh values
C      LP12   - number of steps of omega (365)
C      SH     - table of sinh values
C
C LP12 (from common MGN) is the number of values to calculate. This is set to
C 365  in GOSIA, which is the dimension of the arrays.
 
      SUBROUTINE FHIP
      IMPLICIT NONE
      REAL*8 er , ex , w
      INTEGER*4 j
      INTEGER*4 LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &          LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      REAL*8 CH , SH
      COMMON /HIPER / SH(365) , CH(365)
      
      w = -.03
      DO j = 1 , LP12
         w = w + .03
         ex = EXP(w)
         er = 1./ex
         SH(j) = (ex-er)/2.
         CH(j) = (ex+er)/2.
      ENDDO
      END
 
C----------------------------------------------------------------------
C SUBROUTINE ALLOC
C
C Called by: FTBM, GOSIA
C Calls:     RANGEL
C
C Purpose: to calculate and store the ranges of the integration over omega
C for each multipolarity.
C
C Uses global variables:
C      IRA    - range of omega for each multipolarity needed for accuracy
C      LOCQ   - location of collision function in ZETA array
C      LP14   - maximum length of space for collision functions (4900)
C
C Formal parameters:
C      Accur  - accuracy required
C
C We set up the LOCQ array, which indexes the start of the block of collision
C function coefficients for each omega. For a given multipolarity, lambda,
C there are lambda+1 collision functions, which have to be evaluated for a
C set of different omega values. We don't integrate over all possible omega
C values, but estimate how many values we need to achieve the accuracy Accur.
C The function RANGEL calculates how many we need for each multipolarity and
C stores it in IRA.
C
C Later (in SNAKE) we will store the values of the collision functions for
C 2 * IRA(1) values for E1, 3 * IRA(2) values for E2, 4 * IRA(3) values for
C E3... 3 * IRA(8) values for M2 in that order.
C
C We are limited to a maximum of LP14 (=4900) values in total.
 
      SUBROUTINE ALLOC(Accur)
      IMPLICIT NONE
      REAL*8 Accur
      INTEGER*4 iflag , j , k , k1 , load
      INTEGER*4 LOCQ
      COMMON /ALLC  / LOCQ(8,7)
      INTEGER*4 IRA , MAXLA
      COMMON /RNG   / IRA(8) , MAXLA
      INTEGER*4 LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &          LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
C     Call RANGEL to determine the range of the integration over omega, which
C     depends on the accuracy Accur.
      CALL RANGEL(Accur)

C     First zero all the elements
      load = 0
      iflag = 0
      DO j = 1 , 8
         DO k = 1 , 7
            LOCQ(j,k) = 0
         ENDDO
      ENDDO
      
C     Now store values for E1...E6
      DO k = 1 , 6
         k1 = k + 1
         DO j = 1 , k1
            LOCQ(k,j) = load
            load = load + IRA(k)
         ENDDO
      ENDDO
      
C     And for M1, M2
      DO k = 7 , 8
         k1 = k - 6
         DO j = 1 , k1
            LOCQ(k,j) = load
            load = load + IRA(k) ! IRA(k) is the number of omega values needed for requested accuracy
         ENDDO
      ENDDO

      IF ( load.LE.LP14 ) RETURN ! The Q-functions must fit in the last LP14 words of ZETA

      WRITE (22,99001)
99001 FORMAT (5X,'NO SPACE FOR Q FUNCTIONS TABULATION'//5X,
     &        'SORRY,JOB WILL BE BRUTALLY TERMINATED!')
      STOP 'JOB TERMINATED BY ALLOC' ! Added N. Warr Jul2007
      END
 
C----------------------------------------------------------------------
C SUBROUTINE RANGEL
C
C Called by: ALLOC
C
C Purpose: to determine the range of the integration over omega.
C
C Uses global variables:
C      ACC50  - accuracy required for integration
C      IRA    - range for omega for each multipolarity
C      MULTI  - number of matrix elements with each multipolarity populating levels
C
C Formal parameters:
C      Acc1   - the desired accuracy
C
C \omega_max >= \alpha_\lambda - {1 \over \lambda} \ln(a_c)
C where a_c is Acc1 here.
C
C The gosia documentation gives a table for \alpha_\lambda: E1 = -0.693,
C E2 = 0.203, E3 = 0.536, E4 = 0.716, E5 = 0.829, E6 = 0.962, M1 = 0.203,
C M2 = 0.536.
C
C Note that first we work out omega, but then we work out the appropriate
C index, knowing that we are always using steps of 0.03.

 
      SUBROUTINE RANGEL(Acc1)
      IMPLICIT NONE
      REAL*8 Acc1 , acl , w
      INTEGER*4 i
      REAL*8 ACC50
      COMMON /A50   / ACC50
      INTEGER*4 IRA , MAXLA
      COMMON /RNG   / IRA(8) , MAXLA
      INTEGER*4 LAMDA , LEAD , LDNUM , LAMMAX , MULTI
      COMMON /CLCOM / LAMDA(8) , LEAD(2,1500) , LDNUM(8,100) , LAMMAX , 
     &                MULTI(8)
      acl = -LOG(Acc1)
      ACC50 = Acc1/50.
      DO i = 1 , 8 ! Loop over multipolarity 1..6 = E1..6, 7,8 = M1,M2
         IF ( MULTI(i).NE.0 ) THEN
            IF ( i.EQ.2 .OR. i.EQ.7 ) THEN ! E2 or M1
               w = acl/2. + .203
            ELSEIF ( i.EQ.3 .OR. i.EQ.8 ) THEN ! E3 or M2
               w = acl/3. + .536
            ELSEIF ( i.EQ.4 ) THEN ! E4
               w = acl/4. + .716
            ELSEIF ( i.EQ.5 ) THEN ! E5
               w = acl/5. + .829
            ELSEIF ( i.EQ.6 ) THEN ! E6
               w = acl/6. + .962
            ELSE
               w = acl - .693 ! E1
            ENDIF
            w = w/.03        ! We step in steps of \Delta\omega = 0.03
            IRA(i) = INT(w+1.5)
         ELSE
            IRA(i) = 0
         ENDIF
      ENDDO
      IF ( IRA(7).NE.0 ) IRA(7) = IRA(7) + 1
      IF ( IRA(8).NE.0 ) IRA(8) = IRA(8) + 1
      END
 
C----------------------------------------------------------------------
C SUBROUTINE QRANGE
C
C Called by: SNAKE
C
C Purpose: determine the range for which we will need Qe and Qm values.
C
C Uses global variables:
C      IRA    - range to integrate over omega (readonly)
C      MAXLA  - multipolarity to calculate (writeonly)
C      MULTI  - number of matrix elements having given multipolarity (readonly)
C
C Formal parameters:
C      Icnt   - index of omega to calculate
C      Nlm    - returns the number of l,m values
C      Lloc   - array of l values to calculate
C      Ibm    -
C      Icm    -
C      Idm    -
C      Irl    - range to integrate over omega for each multipolarity
 
      SUBROUTINE QRANGE(Icnt,Nlm,Lloc,Ibm,Icm,Idm,Irl)
      IMPLICIT NONE
      INTEGER*4 Ibm , Icm , Icnt , Idm , Irl , is , k , ke , km , 
     &          l , ld , Lloc , ls
      INTEGER*4 nlend , Nlm
      DIMENSION Lloc(8) , Irl(8)
      INTEGER*4 IRA , MAXLA
      COMMON /RNG   / IRA(8) , MAXLA
      INTEGER*4 LAMDA , LEAD , LDNUM , LAMMAX , MULTI
      COMMON /CLCOM / LAMDA(8) , LEAD(2,1500) , LDNUM(8,100) , LAMMAX , 
     &                MULTI(8)
      IF ( Icnt.EQ.1 ) THEN
         Nlm = 0
         DO l = 1 , 8
            Lloc(l) = 0
            Irl(l) = 0
         ENDDO
         DO k = 1 , 6
            ke = 7 - k
            km = 13 - k
            IF ( km.LE.8 ) THEN
               IF ( MULTI(km).NE.0 ) THEN
                  Nlm = Nlm + 1
                  Lloc(Nlm) = km
                  Irl(Nlm) = IRA(km)
               ENDIF
            ENDIF
            IF ( MULTI(ke).NE.0 ) THEN
               Nlm = Nlm + 1
               Lloc(Nlm) = ke
               Irl(Nlm) = IRA(ke)
            ENDIF
         ENDDO
         nlend = INT((DBLE(Nlm)+1.1)/2.)
         DO k = 1 , nlend
            ke = Nlm - k + 1
            ls = Lloc(ke)
            is = Irl(ke)
            Lloc(ke) = Lloc(k)
            Irl(ke) = Irl(k)
            Lloc(k) = ls
            Irl(k) = is
         ENDDO
         l = 0
         DO k = 1 , 6
            IF ( MULTI(k).NE.0 ) l = k
         ENDDO
         Icm = MIN(4,l)
         Ibm = 2*l
         Idm = l
         l = 0
         DO k = 7 , 8
            ke = k - 6
            IF ( MULTI(k).NE.0 ) l = ke
         ENDDO
         Ibm = MAX(Ibm,2*l)
         Idm = MAX(Idm,l)
         IF ( Icm.EQ.1 .AND. l.GT.1 ) Icm = 2
         MAXLA = Lloc(1)
         RETURN
      ELSE
         IF ( Irl(Nlm).GE.Icnt ) RETURN
         ld = Lloc(Nlm)
         Lloc(Nlm) = 0
         Nlm = Nlm - 1
         IF ( Nlm.EQ.0 ) RETURN
         IF ( ld.GT.6 ) RETURN
         l = Lloc(Nlm)
         IF ( l.GT.6 ) l = l - 6
         Icm = MIN(2,l)
         Ibm = 2*l
         Idm = l
      ENDIF
      END
 
C----------------------------------------------------------------------
C SUBROUTINE AMPDER
C
C Called by: INTG
C Calls:     LAISUM, NEWLV
C
C Purpose: to calculate the derivatives of the amplitudes needed for the
C Adams-Moulton predictor-corrector method.
C
C Uses global variables:
C      ARM    - excitation amplitudes of substates.
C      CAT    - substates of levels (n_level, J, m)
C      ELM    - matrix elements
C      EXPO   - exponents of adiabatic term
C      ISG    - sign of omega
C      ISG1   - index of omega
C      ISMAX  - number of substates used
C      ISSTAR - first substate for given level
C      ISSTO  - last substate for given level
C      LAMDA  - list of multipolarities to calculate
C      LAMMAX - number of multipolarities to calculate
C      LAMR   - flag = 1 if we should calculate this multipolarity
C      LZETA  - index in ZETA to coupling coefficients for given multipolarity
C      MSTORE - index of final level number and index of matrix element
C      NMAX   - number of levels
C      NPT    - index in ADB array (this is omega / 0.03)
C      NSTART - index in CAT of first substate associated with a level
C      NSTOP  - index in CAT of last substate associated with a level
C
C Formal parameters:
C      I57    - switch which is either 5 or 7. This tells LAISUM to access either ARM(I,5) or ARM(I,7)

      SUBROUTINE AMPDER(I57)
      IMPLICIT NONE
      REAL*8 rsg
      INTEGER*4 i1 , I57 , ibg , iend , iflg , indx , ir , is2 , k , 
     &          lam , lax , ld , m , mm , n , nhold , nz
      COMPLEX*16 ARM
      COMMON /AZ    / ARM(1200,7)
      REAL*8 ELM , ELMU , ELML , SA
      COMMON /COMME / ELM(1500) , ELMU(1500) , ELML(1500) , SA(1500)
      INTEGER*4 LAMDA , LEAD , LDNUM , LAMMAX , MULTI
      COMMON /CLCOM / LAMDA(8) , LEAD(2,1500) , LDNUM(8,100) , LAMMAX , 
     &                MULTI(8)
      INTEGER*4 NMAX , NDIM , NMAX1
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      REAL*8 D2W
      INTEGER*4 NPT , NDIV , KDIV , LAMR , ISG , NSW , ISG1
      COMMON /CAUX  / NPT , NDIV , KDIV , LAMR(8) , ISG , D2W , NSW , 
     &                ISG1
      INTEGER*4 ISSTAR , ISSTO , MSTORE
      COMMON /PINT  / ISSTAR(101) , ISSTO(100) , MSTORE(2,100)
      COMPLEX*16 EXPO
      COMMON /ADBXI / EXPO(1500)
      REAL*8 ZETA
      INTEGER*4 LZETA
      COMMON /CCOUP / ZETA(155600) , LZETA(8)
      REAL*8 CAT
      INTEGER*4 ISMAX
      COMMON /CLCOM8/ CAT(1200,3) , ISMAX
      INTEGER*4 NSTART , NSTOP
      COMMON /CEXC0 / NSTART(101) , NSTOP(100)
C     Zero ARM(k,4) and ARM(k,6) for each substate used
      DO k = 1 , ISMAX ! ISMAX is number of substates used
         ARM(k,6) = (0.,0.)
         ARM(k,4) = (0.,0.)
      ENDDO

      ISG1 = ISG
      IF ( NPT.EQ.1 ) ISG1 = ABS(ISG1)
      rsg = DBLE(ISG)

      DO i1 = 1 , LAMMAX ! LAMMAX is number of multipolarities to calculate
         lam = LAMDA(i1) ! For each value of lambda, the user wants to calculate for
         lax = lam
         nz = LZETA(lam) ! Index into ZETA array for each multipolarity
         IF ( LAMR(lam).NE.0 ) THEN ! LAMR is flag to decide if we calculate for this multipolarity
            iflg = 1
            nhold = 1
 20         CALL NEWLV(nhold,ld,lam)
            IF ( ld.EQ.0 ) THEN ! If there are no decays
 30            nhold = nhold + 1
               IF ( NSTART(nhold).NE.0 ) GOTO 20
               GOTO 30
            ELSE
               ir = NSTART(nhold) - 1 ! Get first substate - 1 for this level
 40            ir = ir + 1 ! ir is a substate
               IF ( ir.LE.ISMAX ) THEN
                  n = CAT(ir,1) ! Level number of substate ir
                  IF ( n.NE.nhold ) THEN
                     DO mm = 1 , ld ! Loop over matrix elements
                        m = MSTORE(1,mm) ! Index of final level
                        IF ( m.NE.nhold ) THEN
                           indx = MSTORE(2,mm) ! Index of matrix element in ELM
                           ibg = ISSTAR(mm) ! First substate for this level
                           iend = ISSTO(mm) ! Last substate for this level
                           DO is2 = ibg , iend ! Loop over substates for level
                              ARM(is2,4) = ARM(is2,4) + ARM(is2,6)
     &                           *ELM(indx)/EXPO(indx)
                              ARM(is2,6) = (0.,0.)
                           ENDDO
                        ENDIF
                     ENDDO
 42                  CALL NEWLV(n,ld,lam)
                     IF ( ld.EQ.0 ) THEN ! if ld is zero, skip all the states for this level
                        ir = ir + NSTOP(n) - NSTART(n) + 1
                        n = n + 1
                        IF ( n.LE.NMAX ) GOTO 42
                        GOTO 100 ! IF this was the last level, loop back over lambda
                     ELSE
                        nhold = n
                     ENDIF
                  ENDIF ! If n .ne. nhold
                  CALL LAISUM(ir,n,rsg,lax,ld,nz,I57)
                  GOTO 40
               ENDIF ! If IR .le ISMAX
            ENDIF ! If LD .ne. 0
         ENDIF ! If LAMR(lam) .ne. 0
 100     CONTINUE
      ENDDO ! Loop over lambda
      END
 
C----------------------------------------------------------------------
C SUBROUTINE LAISUM
C
C Called by: AMPDER, STING
C Calls:     FAZA
C
C Purpose: evaluate the sum  over matrix elements.
C
C Uses global variables:
C      ARM    - excitation amplitudes of substates.
C      CAT    - substates of levels (n_level, J, m)
C      ELM    - matrix elements
C      EXPO   - adiabatic exponential
C      ISG    - sign of omega
C      ISG1   - index of omega
C      ISHA   - is half-integer spin
C      ISO    - isotropic flag
C      ISSTAR - first substate for given matrix element index
C      ISSTO  - last substate for given matrix element index
C      KDIV   - index for division
C      LOCQ   - location of collision functions in ZETA array
C      LP7    - start of collision functions in ZETA (45100)
C      MSTORE - index of final level number and index of matrix element
C      NDIV   - number of divisions
C      NPT    - index in ADB array (this is omega / 0.03)
C      NSTART - index in CAT of first substate associated with a level
C      NSTOP  - index in CAT of last substate associated with a level
C      ZETA   - various coefficients
C
C Formal parameters:
C      Ir     - index of substate
C      N      - index of level
C      Rsg    - sign of omega
C      Lam    - multipolarity
C      Ld     - number of levels connected to level N by this multipolarity Lam
C      Nz     - index into ZETA array for this multipolarity
C      I57    - switch which is either 5 or 7 so we access ARM(I,5) or ARM(I,7)
C
C   \sum_{lmn} \zeta^{lm}_{kn} . M^(1)_{kn} f_{lm}(\omega) a_n(\omega)
C where
C   f_{lm} = -i Q_{lm} exp(i \xi_{kn} (\epsilon \sinh(\omega) + \omega))
C and
C   M^(1)_{kn} = <k||E(M)\lambda||n>
C
C EXPO is exp(i * xi * sinh(w) + w) calculated in function EXPON.
C ARM are the excitation amplitudes of the substates.
C q is the Qe or Qm calculated by the functions QE and QM, respectively and
C stored in ZETA array in the function SNAKE.
C z is the coupling parameter zeta, calculated in the function LSLOOP.
      
      SUBROUTINE LAISUM(Ir,N,Rsg,Lam,Ld,Nz,I57)
      IMPLICIT NONE
      REAL*8 q , rmir , rmis , rmu , Rsg , z
      INTEGER*4 i2 , i3 , I57 , iii , indq , indx , Ir , irs , is , 
     &          is1 , is2 , ismin , isplus
      INTEGER*4 la , Lam , Ld , m , mrange , mua , N , Nz
      COMPLEX*16 FAZA , pamp , pamp1
      INTEGER*4 ISHA
      COMMON /PSPIN / ISHA
      COMPLEX*16 ARM
      COMMON /AZ    / ARM(1200,7)
      REAL*8 EN , SPIN , ACCUR , DIPOL , ZPOL , ACCA
      INTEGER*4 ISO
      COMMON /COEX  / EN(100) , SPIN(100) , ACCUR , DIPOL , ZPOL , 
     &                ACCA , ISO
      REAL*8 D2W
      INTEGER*4 NPT , NDIV , KDIV , LAMR , ISG , NSW , ISG1
      COMMON /CAUX  / NPT , NDIV , KDIV , LAMR(8) , ISG , D2W , NSW , 
     &                ISG1
      INTEGER*4 ISSTAR , ISSTO , MSTORE
      COMMON /PINT  / ISSTAR(101) , ISSTO(100) , MSTORE(2,100)
      COMPLEX*16 EXPO
      COMMON /ADBXI / EXPO(1500)
      REAL*8 ZETA
      INTEGER*4 LZETA
      COMMON /CCOUP / ZETA(155600) , LZETA(8)
      INTEGER*4 LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &          LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      REAL*8 CAT
      INTEGER*4 ISMAX
      COMMON /CLCOM8/ CAT(1200,3) , ISMAX
      REAL*8 ELM , ELMU , ELML , SA
      COMMON /COMME / ELM(1500) , ELMU(1500) , ELML(1500) , SA(1500)
      INTEGER*4 LOCQ
      COMMON /ALLC  / LOCQ(8,7)
      INTEGER*4 NSTART , NSTOP
      COMMON /CEXC0 / NSTART(101) , NSTOP(100)
      
      rmir = CAT(Ir,3) ! m quantum number of substate Ir
      iii = 0
      IF ( Lam.GT.6 ) iii = 1
      la = Lam
      IF ( Lam.GT.6 ) Lam = Lam - 6
      DO i2 = 1 , Ld
         pamp = (0.,0.)
         m = MSTORE(1,i2) ! Index of final level
         indx = MSTORE(2,i2) ! Index of matrix element in ELM
         ismin = 0
         is1 = NSTART(m) ! Index of first substate for level m
         IF ( is1.NE.0 ) THEN
            isplus = INT(rmir-CAT(is1,3)) - Lam
            IF ( isplus.LT.0 ) THEN
               ismin = isplus
               isplus = 0
            ENDIF
            is2 = is1 + isplus - 1
            mrange = 2*Lam + 1 + ismin
            IF ( is2+mrange.GT.NSTOP(m) ) mrange = NSTOP(m) - is2
            IF ( mrange.GT.0 ) THEN
               DO i3 = 1 , mrange
                  is = is2 + i3
                  rmis = CAT(is,3) ! m quantum number of substate is
                  IF ( ISO.NE.0 .OR. rmir.LE..1 .OR. rmis.LE..1 ) THEN
                     rmu = rmis - rmir
                     mua = ABS(rmu) + 1.1 ! delta-m + 1
C                    Only consider electromagnetic and delta-m .NE. 0 magnetic
C                    contribution
                     IF ( la.LE.6 .OR. mua.NE.1 ) THEN
                        indq = LOCQ(Lam,mua) + NPT ! Index to Q function
                        Nz = Nz + 1                ! Index to Zeta
                        z = ZETA(Nz)               ! Zeta
                        q = ZETA(indq+LP7)         ! Q-function
                        IF ( NDIV.NE.0 ) q = ZETA(indq+LP7) + DBLE(KDIV)
     &                       *(ZETA(indq+LP7+ISG1)-ZETA(indq+LP7))
     &                       /DBLE(NDIV)
                        pamp1 = FAZA(la,mua,rmu,Rsg)*q*z
                        IF ( ISO.NE.0 .OR. rmir.LE..1 ) THEN
                           pamp = pamp1*ARM(is,I57) + pamp
                           IF ( ISO.EQ.0 .AND. rmis.GT..1 ) GOTO 10
                        ENDIF
                        IF ( N.NE.m ) THEN ! Same level
                           irs = (-1)**(INT(rmir+rmis)-ISHA+iii) ! ISHA = 1 if half-integer spins, iii=0 for E, 1 for M
                           ARM(is,6) = ARM(is,6) + irs*pamp1*ARM(Ir,I57)
                           ISSTAR(i2) = MIN(is,ISSTAR(i2))
                           ISSTO(i2) = MAX(is,ISSTO(i2))
                        ENDIF
                     ENDIF
                  ENDIF ! If ISO.NE. 0 or either substate is spin 1
 10               CONTINUE
               ENDDO ! Loop on mrange
               IF ( N.EQ.m ) THEN ! N and m are level numbers, so if it is the same level, EXPO = 1
                  ARM(Ir,4) = ARM(Ir,4) + pamp*ELM(indx)
               ELSE
                  ARM(Ir,4) = ARM(Ir,4) + pamp*ELM(indx)*EXPO(indx)
               ENDIF
            ENDIF
         ENDIF
      ENDDO ! Loop on levels
      Lam = la
      END

C-----------------------------------------------------------------------
C FUNCTION EXPON
C
C Called by: NEWLV
C Calls:     TCEXP
C
C Purpose: calculates the exponential:
C       exp(i \xi_{kn} (\epsilon \sinh(\omega) + \omega))
C
C Uses global variables:
C      EXPO   - adiabatic exponential
C      ADB    - adiabatic function
C      XI     - xi coupling coefficients
C
C Formal parameters:
C      Inx    - index in XI array
C      Npt    - index in ADB array (this is omega / 0.03)
C      Isg    - sign of omega
C      Isg1   - index of omega
C      Ndiv   - number of divisions
C      Kdiv   - index for division
C
C Return value:
C      the exponential
C
C ci is sqrt(-1)
C XI (from common block CXI) are the XI coupling constants calculated in
C the function LOAD.
C ADB is the adiabatic parameters \epsilon \sinh(\omega) + \omega calculated
C in the function SETIN.

      COMPLEX*16 FUNCTION EXPON(Inx,Npt,Isg,Isg1,Ndiv,Kdiv)
      IMPLICIT NONE
      INTEGER*4 Inx , Isg , Isg1 , Kdiv , Ndiv , Npt
      COMPLEX*16 expo1 , ci , expox , TCEXP
      REAL*8 ADB
      COMMON /ADX   / ADB(365)
      REAL*8 XI
      COMMON /CXI   / XI(1500)
      DATA ci/(0.,1.)/
      
      expox = TCEXP(ci*XI(Inx)*ADB(Npt)*Isg)
      EXPON = expox
      IF ( Ndiv.NE.0 ) THEN
         expo1 = TCEXP(ci*XI(Inx)*ADB(Npt+Isg1)*Isg)
         EXPON = expox + DBLE(Kdiv)*(expo1-expox)/DBLE(Ndiv)
      ENDIF
      END
 
C----------------------------------------------------------------------
C FUNCTION FAZA
C
C Called by: LAISUM
C
C Purpose: calculate phase. QE and QM calculate only the magnitude of the
C          collision function for positive omega. Here we multiply by -i and
C          take into account whether the collision function is real or
C          imaginary. We also multiply by the sign of omega if the collision
C          function is asymmetric wrt. omega.
C
C Formal parameters:
C      La     - lambda 1...6 = E1...6, 7,8 = M1,M2
C      Mi     - mu + 1
C      Rmu    - mu
C      Rsg    - sign of omega
C Return value:
C      Phase

      COMPLEX*16 FUNCTION FAZA(La,Mi,Rmu,Rsg)
      IMPLICIT NONE
      INTEGER*4 ieven , La , Mi
      REAL*8 Rmu , Rsg
      COMPLEX*16 ci
      DATA ci/(0.,1.)/ ! sqrt(-1)

      IF ( La.GT.6 ) THEN ! M1, M2 multipolarity
         FAZA = -ci
         IF ( Rmu.LT.0. ) FAZA = -FAZA
         IF ( La.EQ.7 ) RETURN
         IF ( Mi.EQ.2 ) RETURN
         FAZA = DCMPLX(Rsg,0.D0)
         IF ( Rmu.LT.0. ) FAZA = -FAZA
         GOTO 99999
      ELSE                ! E1...6
         ieven = (-1)**Mi
         IF ( ieven.LE.0 ) THEN
            FAZA = -ci ! mu is even
            RETURN
         ENDIF
      ENDIF
      FAZA = DCMPLX(Rsg,0.D0) ! mu is odd, so sign changes with sign of omega
99999 END
 
C----------------------------------------------------------------------
C SUBROUTINE SETIN
C
C Called by FTBM, GOSIA
C
C Purpose: calculate the adiabatic parameter:
C \epsilon \sinh(\omega) + \omega
C
C Uses global variables:
C      ADB    - adiabatic function
C      EPS    - epsilon
C      IEXP   - experiment number
C      LP12   - number of steps of omega (365)
C      SH     - table of sinh values
C
C Note that it uses the tables of sinh calculated by FHIP (SH in common
C block HIPER) and that both the sinh table and this table of the adiabatic
C parameter are in steps of \Delta\omega = 0.03. The resulting table of
C adiabatic parameters are stored in ADB in common block ADX.
C
C LP12 (from common MGN) is the number of values to calculate. This is set to
C 365  in GOSIA, which is the dimension of the array.
 
      SUBROUTINE SETIN
      IMPLICIT NONE
      INTEGER*4 k
      INTEGER*4 LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &          LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      REAL*8 CH , SH
      COMMON /HIPER / SH(365) , CH(365)
      REAL*8 ADB
      COMMON /ADX   / ADB(365)
      REAL*8 EPS, EROOT, FIEX
      INTEGER*4 IEXP, IAXS
      COMMON /KIN   / EPS(50) , EROOT(50) , FIEX(50,2) , IEXP , IAXS(50)
      
      DO k = 1 , LP12
         ADB(k) = EPS(IEXP)*SH(k) + .03*(k-1)
      ENDDO
      END
 
C----------------------------------------------------------------------
C SUBROUTINE STING
C
C Called by: FTBM, GOSIA
C Calls:     LAIAMP, LAISUM, NEWLV
C
C Purpose: calculate and store excitation amplitudes,
C
C Uses global variables:
C      ARM    - excitation amplitudes of substates.
C      ELM    - matrix elements
C      EXPO   - adiabatic exponential
C      IFLG   - flag to determine whether to calculate exponential (so we don't calculate twice)
C      IRA    - limit of omega for integration for each multipolarity
C      ISG    - sign of omega
C      ISMAX  - number of substates used
C      ISSTAR - first substate for given level
C      ISSTO  - last substate for given level
C      KDIV   - index of division
C      LAMDA  - list of multipolarities to calculate
C      LAMMAX - number of multipolarities to calculate
C      LAMR   - flag = 1 if we should calculate this multipolarity
C      LDNUM  - number of matrix elements with each multipolarity populating levels
C      LZETA  - index in ZETA to coupling coefficients for a given multipolarity
C      MAXLA  - multipolarity to calculate here
C      MSTORE - index of final level number and index of matrix element
C      NDIV   - number of divisions
C      NPT    - index in ADB array (this is omega / 0.03)
C
C Formal parameters:
C      Irld   - index into ARM array
 
      SUBROUTINE STING(Irld)
      IMPLICIT NONE
      REAL*8 rsg , w0
      INTEGER*4 i , i57 , ibg , iend , indx , Irld , is2 , j , j1 , 
     &          jj , lam , ld , maxh , mm , n , nz
      INTEGER*4 LAMDA , LEAD , LDNUM , LAMMAX , MULTI
      COMMON /CLCOM / LAMDA(8) , LEAD(2,1500) , LDNUM(8,100) , LAMMAX , 
     &                MULTI(8)
      COMPLEX*16 ARM
      COMMON /AZ    / ARM(1200,7)
      REAL*8 ELM , ELMU , ELML , SA
      COMMON /COMME / ELM(1500) , ELMU(1500) , ELML(1500) , SA(1500)
      COMPLEX*16 EXPO
      COMMON /ADBXI / EXPO(1500)
      INTEGER*4 IFLG
      COMMON /FLA   / IFLG
      INTEGER*4 ISSTAR , ISSTO , MSTORE
      COMMON /PINT  / ISSTAR(101) , ISSTO(100) , MSTORE(2,100)
      REAL*8 ZETA
      INTEGER*4 LZETA
      COMMON /CCOUP / ZETA(155600) , LZETA(8)
      REAL*8 D2W
      INTEGER*4 NPT , NDIV , KDIV , LAMR , ISG , NSW , ISG1
      COMMON /CAUX  / NPT , NDIV , KDIV , LAMR(8) , ISG , D2W , NSW , 
     &                ISG1
      REAL*8 CAT
      INTEGER*4 ISMAX
      COMMON /CLCOM8/ CAT(1200,3) , ISMAX
      INTEGER*4 IRA , MAXLA
      COMMON /RNG   / IRA(8) , MAXLA
      maxh = MAXLA ! Save MAXLA, so we can restore it later
 100  ISG = -1
      n = 1
      rsg = -1.
      IFLG = 1
      w0 = IRA(MAXLA)*.03 + .03 ! Maximum omega to calculate for (steps of 0.03)
      
      DO j = 1 , ISMAX ! For substate used, zero ARM array
         DO jj = 1 , 6
            ARM(j,jj) = (0.,0.)
         ENDDO
      ENDDO
      ARM(Irld,5) = (1.,0.) ! Set y_n to 1

      DO j = 1 , 8
         LAMR(j) = 0 ! Initially mark that we shouldn't calculate any multipolarity
      ENDDO
      LAMR(MAXLA) = 1 ! Mark that we should calculate this multipolarity
      
      NPT = IRA(MAXLA) + 1 ! Number of omega values to calculate for this multipolarity
      
      IF ( MAXLA.EQ.7 .AND. IRA(2).NE.0 ) THEN ! Special case of M1
         LAMR(2) = 1
         NPT = NPT - 1
         w0 = w0 - .03
      ENDIF

      NDIV = 0
      KDIV = 0

      DO j = 1 , 4 ! Loop over terms for Adams-Moulton corrector-predictor
         NPT = NPT - 1
         DO j1 = 1 , LAMMAX ! Loop up to maximum multipolarity to calculate
            lam = LAMDA(j1) ! Get the multipolarity
            IF ( LAMR(lam).NE.0 ) THEN ! If this multipolarity should be calculated

C              Calculate and store exponentials in EXPO
               CALL NEWLV(n,ld,lam) ! n = 1, so for ground-state

               IF ( ld.NE.0 ) THEN ! If there are levels connected to g.s. by the right multipolarity
                  nz = LZETA(lam) ! Index into zeta array for this multipolarity
                  ld = LDNUM(lam,1) ! Number of levels connected to ground state for this multipolarity
                  i57 = 5 ! Use ARM(I,5) in LAISUM for excitation amplitudes
C                 Calculate sum over matrix elements
                  CALL LAISUM(Irld,n,rsg,lam,ld,nz,i57)
                  DO mm = 1 , ld ! Loop over levels
                     indx = MSTORE(2,mm) ! Index of matrix element in ELM
                     ibg = ISSTAR(mm) ! First substate for this level
                     iend = ISSTO(mm) ! Last substate for this level
                     DO is2 = ibg , iend ! Loop over substates
                        ARM(is2,4) = ARM(is2,4) + ARM(is2,6)*ELM(indx)
     &                               /EXPO(indx)
                        ARM(is2,6) = (0.,0.)
                     ENDDO ! Loop over substates
                  ENDDO ! Loop over matrix elements for ground state for this multipolarity
               ELSEIF ( j1.EQ.MAXLA ) THEN ! Else if it is the last multipolarity
                  IRA(MAXLA) = -IRA(MAXLA) ! Make IRA negative
                  DO jj = 1 , LAMMAX
                     lam = LAMDA(jj)
                     IF ( IRA(lam).GT.0 ) GOTO 105
                  ENDDO
 105              MAXLA = LAMDA(jj) ! Advance MAXLA to next multipolarity
                  GOTO 100 ! Back to start
               ENDIF

            ENDIF ! If we should calculate this multipolarity
         ENDDO ! Loop over multipolarities
         IF ( j.EQ.4 ) GOTO 200 ! We've set everything up, so finish
         
         DO i = 1 , ISMAX ! Shift terms up one
            ARM(i,j) = ARM(i,4)
            ARM(i,4) = (0.,0.)
         ENDDO

      ENDDO ! Loop over terms for Adams-Moulton corrector-predictor

C     Calculate amplitude
 200  CALL LAIAMP(Irld,w0)
      
      MAXLA = maxh ! Restore MAXLA

      DO jj = 1 , 8
         IRA(jj) = ABS(IRA(jj)) ! Make sure all the IRA are positive again
      ENDDO
      END
 
C----------------------------------------------------------------------
C SUBROUTINE LAIAMP
C
C Called by: STING
C Calls:     FAZA1, LEADF, STAMP, TCABS
C
C Purpose: calculate excitation amplitudes
C
C Uses global variables:
C      ARM    - excitation amplitudes of substates.
C      CAT    - substates of levels (n_level, J, m)
C      ELM    - matrix elements
C      EPS    - epsilon
C      EROOT  - sqrt(epsilon^2 - 1)
C      IEXP   - number of experiment
C      ISG    - sign of omega
C      LAMDA  - list of multipolarities to calculate
C      LAMMAX - number of multipolarities to calculate
C      LAMR   - flag = 1 if we should calculate this multipolarity
C      LDNUM  - number of matrix elements with each multipolarity populating levels
C      LZETA  - index in ZETA to coupling coefficients for given multipolarity
C      NSTART - index in CAT of first substate associated with a level
C      NSTOP  - index in CAT of last substate associated with a level
C      XI     - xi coupling coefficients
C      ZETA   - various coefficients
C
C Formal parameters:
C      Ir     - index of substate
C      W0     - omega limit
      
      SUBROUTINE LAIAMP(Ir,W0)
      IMPLICIT NONE
      REAL*8 epsi , errt , pm , ppp , rmir , rmis , rmu , TCABS , W0 , 
     &       xiv , z
      INTEGER*4 i1 , i2 , i3 , indx , Ir , is , is1 , is2 , ismin , 
     &          isplus , la , lam
      INTEGER*4 ld , LEADF , m , MEM , mrange , mua , nz
      COMPLEX*16 STAMP , dis , uhuj
      INTEGER*4 LAMDA , LEAD , LDNUM , LAMMAX , MULTI
      COMMON /CLCOM / LAMDA(8) , LEAD(2,1500) , LDNUM(8,100) , LAMMAX , 
     &                MULTI(8)
      COMPLEX*16 ARM
      COMMON /AZ    / ARM(1200,7)
      REAL*8 D2W
      INTEGER*4 NPT , NDIV , KDIV , LAMR , ISG , NSW , ISG1
      COMMON /CAUX  / NPT , NDIV , KDIV , LAMR(8) , ISG , D2W , NSW , 
     &                ISG1
      REAL*8 ZETA
      INTEGER*4 LZETA
      COMMON /CCOUP / ZETA(155600) , LZETA(8)
      REAL*8 CAT
      INTEGER*4 ISMAX
      COMMON /CLCOM8/ CAT(1200,3) , ISMAX
      REAL*8 ELM , ELMU , ELML , SA
      COMMON /COMME / ELM(1500) , ELMU(1500) , ELML(1500) , SA(1500)
      INTEGER*4 NSTART , NSTOP
      COMMON /CEXC0 / NSTART(101) , NSTOP(100)
      REAL*8 EPS, EROOT, FIEX
      INTEGER*4 IEXP, IAXS
      COMMON /KIN   / EPS(50) , EROOT(50) , FIEX(50,2) , IEXP , IAXS(50)
      REAL*8 XI
      COMMON /CXI   / XI(1500)
      ppp = 0.
      epsi = EPS(IEXP) ! epsilon
      errt = EROOT(IEXP) ! sqrt(epsilon^2 - 1)
      rmir = CAT(Ir,3) ! m quantum number of substate Ir
      
      DO i1 = 1 , LAMMAX ! Loop on multipolarity
         lam = LAMDA(i1) ! Get multipolarity
         nz = LZETA(lam) ! nz is an index into ZETA array for this multipolarity
         IF ( LAMR(lam).NE.0 ) THEN
            la = lam
            IF ( lam.GT.6 ) lam = lam - 6 ! la = 7,8 for M1,M2
            ld = LDNUM(la,1) ! Number of matrix elements with multipolarity la, connecting to ground state
            IF ( ld.NE.0 ) THEN
               DO i2 = 1 , ld ! Loop on matrix elements of that multipolarity connected to ground state
                  m = LEADF(1,i2,la) ! m is level index connected to ground state by element i2, mul. la
                  indx = MEM(1,m,la) ! indx is the index of the matrix element connecting this level to the ground state with this multipolarity
                  xiv = XI(indx) ! xi value
                  ismin = 0
                  is1 = NSTART(m) ! Index of first substate for level m
                  IF ( NSTART(m).NE.0 ) THEN
                     isplus = INT(rmir-CAT(is1,3)) - lam
                     IF ( isplus.LT.0 ) THEN
                        ismin = isplus
                        isplus = 0
                     ENDIF
                     is2 = is1 + isplus - 1
                     mrange = 2*lam + 1 + ismin
                     IF ( is2+mrange.GT.NSTOP(m) ) mrange = NSTOP(m)
     &                    - is2
                     IF ( mrange.GT.0 ) THEN ! If there are substates for level m
                        DO i3 = 1 , mrange
                           is = is2 + i3
                           nz = nz + 1
                           z = ZETA(nz) ! zeta coefficient
                           rmis = CAT(is,3) ! m quantum number of substate is
                           rmu = rmis - rmir
                           mua = ABS(rmu) + 1.1 ! delta-mu + 1

C                          Only consider electromagnetic and delta-mu .NE. 0 magnetic
C                          contribution
                           IF ( lam.LE.6 .OR. mua.NE.1 ) THEN
C                             calculate complex phase (dis)
                              CALL FAZA1(la,mua,rmir,rmis,dis,rmu)
                              pm = ELM(indx)*z ! Matrix element * zeta
C                             estimate amplitude
                              uhuj = STAMP(epsi,errt,xiv,.03D0,W0,lam,
     &                               mua)
                              ARM(is,5) = dis*pm*uhuj
                              ppp = ppp + TCABS(ARM(is,5))
     &                              *TCABS(ARM(is,5))
                           ENDIF

                        ENDDO ! Loop over substates
                     ENDIF ! If there are substates
                  ENDIF ! If there are substates for level m
               ENDDO ! Loop on matrix elements connected to ground state with multipolarity la
            ENDIF ! If there are matrix elements of this multipolarity connecting to the ground state
         ENDIF
      ENDDO ! Loop on lambda
      ARM(Ir,5) = DCMPLX(SQRT(1.-ppp),0.D0)
      END
 
C----------------------------------------------------------------------
C SUBROUTINE FAZA1
C
C Called by: LAIAMP
C
C Purpose: calculate complex phase. Note it is only called for E1...6, M1 with
C          mu non zero and not at all for M2.
C
C Formal parameters:
C      La     - lambda 1...6 = E1...6, 7,8 = M1,M2
C      Mi     - mu + 1
C      Rmir   - m of first state
C      Rmis   - m of second state
C      Dis    - complex phase
C      Rmu    - mu

      SUBROUTINE FAZA1(La,Mi,Rmir,Rmis,Dis,Rmu)
      IMPLICIT NONE
      INTEGER*4 ieven , irs , La , Mi
      REAL*8 Rmir , Rmis , Rmu
      COMPLEX*16 Dis , ci
      DATA ci/(0.,1.)/ ! sqrt(-1)

      irs = (-1)**INT(Rmir+Rmis)
      IF ( La.EQ.7 ) THEN ! M1
         Dis = -ci*irs
         IF ( Rmu.LT.0. ) Dis = -Dis
         GOTO 99999
      ELSE ! E1...6
         ieven = (-1)**Mi
         IF ( ieven.LE.0 ) THEN
            Dis = -ci*irs
            RETURN
         ENDIF
      ENDIF
      Dis = DCMPLX(-DBLE(irs),0.D0)
99999 END
 
C----------------------------------------------------------------------
C SUBROUTINE TRINT
C
C Called by: STAMP
C Calls:     POL4
C
C Purpose: calculate sine and cosine integrals (Si and Ci). Note, that we
C actually calculate pi/2-Si and -Ci. Note also that a constant added to these
C values doesn't make any difference, because STAMP always subtracts them
C pairwise from each other.
C
C Formal parameters:
C      Arg    - value of x for which to evaluate sine and cosine integrals
C      Si     - returned sine integral at that value (actually pi/2 - Si)
C      Ci     - returned cosine integral at that value (actually -Ci)
C
C For small x we use the series expansion. See Abramowitz and Stegun Handbook
C of Mathematical Functions with Formulas, Graphs and Mathematical Tables,
C National Bureau of Standards, 8th Ed. P232 Eqs. 5.2.14 and 5.2.16, except we
C calculate pi/2 - Si and -Ci:
C
C pi/2 - Si = pi/2 - x + x^3 / (3! * 3) - x^5 / (5! * 5) + x^7 / (7! * 7) + ...
C pi/2                       = 1.57079632679
C 1 / (3! * 3) = 1 / 18      = 0.05555555
C 1 / (5! * 5) = 1 / 600     = 1.666667E-3
C 1 / (7! * 7) = 1 / 35280   = 2.83446E-5
C
C -Ci = -Gamma - ln(x) + x^2 / (2! * 2) - x^4 / (4! * 4) + x^6 / (6! * 6) - ...
C Gamma        = Euler Gamma = 0.577215664902
C 1 / (2! * 2) = 1 / 4       = 0.25
C 1 / (4! * 4) = 1 / 96      = 0.0104166
C 1 / (6! * 6) = 1 / 4320    = 2.31481E-4
C 1 / (8! * 8) = 1 / 322560  = 3.10019E-6
C
C For large x we use the rational approximations. See Abramowitz and Stegun
C Handbook of Mathematical Functions with Formulas, Graphs and Mathematical
C Tables, National Bureau of Standards, 8th Ed. P233 Eqs. 5.2.38 and 5.2.39 to
C calculate the auxillary functions f and g and then use 5.2.8 and 5.2.9 to
C obtain the values of pi/2 - Si and -Ci from f and g.
      
      SUBROUTINE TRINT(Arg,Si,Ci)
      IMPLICIT NONE
      REAL*8 a , Arg , c , Ci , f , g , POL4 , s , Si

      a = Arg*Arg

C     If Arg is small, use the polynomial expansion. The coefficients are
C     evaluated from Abramowitz and Stegun 5.2.14 and 5.2.16 as shown above:
      IF ( Arg.LT.1. ) THEN
         Si = POL4(0.D0,2.83446712D-5,-1.66666667D-3,.055555555D0,-1.D0,
     &        a)
         Si = Si*Arg
         Si = Si + 1.57079632679D0 ! This is actually pi/2 - Si
         Ci = POL4(-3.100198413D-6,2.314814815D-4,-.0104166667D0,.25D0,
     &        0.D0,a)
         Ci = Ci - LOG(Arg) - 0.577215664902D0 ! This is actually -Ci
         GOTO 99999
      ENDIF

C     Otherwise use the expansion in terms of sine and cosine
      s = SIN(Arg)
      c = COS(Arg)

C     Here we use an approximation. If Arg is quite large, a is very large 
C     and the four polynomials are all huge. Moreover, the four polynomials 
C     are almost identical, so the ratios are unity. So in this case, 
C     f = 1./Arg and g=1./a is a good approximation.
      
      f = 1.
      g = 1.

C     From Abramowitz and Stegun 5.2.38 and 5.2.39 we have the following
C     relations for the auxillary functions f and g, using the coefficients
C     from that reference:
      IF ( a.LE.1.D+8 ) THEN
         f = POL4(1.D0,38.027246D0,265.187033D0,335.67732D0,38.102495D0,
     &       a)
         f = f/POL4(1.D0,40.021433D0,322.624911D0,570.23628D0,
     &       157.105423D0,a)
         g = POL4(1.D0,42.242855D0,302.757865D0,352.018498D0,
     &       21.821899D0,a)
         g = g/POL4(1.D0,48.196927D0,482.485984D0,1114.978885D0,
     &       449.690326D0,a)
      ENDIF
      
      f = f/Arg
      g = g/a

C     From Abramowitz and Stegun 5.2.8 and 5.2.9:      
      Si = f*c + g*s ! This is actually pi/2 - Si compared to Abramowitz and Stegun
      Ci = g*c - f*s ! This is actually -Ci compared to Abramowitz and Stegun
99999 END
 
C----------------------------------------------------------------------
C FUNCTION POL4
C
C Called by: TRINT
C
C Purpose: evaluate a 4th order polynomial
C
C Formal parameters:
C      C0     - coefficient of polynomial
C      C1     - coefficient of polynomial
C      C2     - coefficient of polynomial
C      C3     - coefficient of polynomial
C      C4     - coefficient of polynomial
C
C Evaluate C0 * A^4 + C1 * A^3 + C2 * A^2 + C3 * A + C4
 
      REAL*8 FUNCTION POL4(C0,C1,C2,C3,C4,A)
      IMPLICIT NONE
      REAL*8 A , C0 , C1 , C2 , C3 , C4

      POL4 = C4 + A*(C3+A*(C2+A*(C1+A*C0)))
      END
 
C----------------------------------------------------------------------
C FUNCTION STAMP
C
C Called by: LAIAMP
C Calls:     TRINT
C
C Purpose: Estimate amplitude
C
C Formal parameters:
C      Epsi   - epsilon for this experiment
C      Errt   - sqrt(epsilon^2 - 1) for this experiment
C      Xiv    - value of xi
C      Dw     - step in omega (0.03)
C      W0     - value of omega
C      Lmda   - lambda (1...2 for E1...2 and 7 for M1 - other multipolarites forbidden)
C      Mua    - mu + 1
C
C Note that the pre-factors in the fct values correspond to those of the
C collision functions.
C
C Return value:
C      Estimated amplitude
 
      COMPLEX*16 FUNCTION STAMP(Epsi,Errt,Xiv,Dw,W0,Lmda,Mua)
      IMPLICIT NONE
      REAL*8 a , axi , b , bic , bic2 , bis , bis2 , ca , cb , cia , 
     &       cib , cic , cis , Dw , dwi , Epsi , Errt , ex , exa , fct
      INTEGER*4 la , Lmda , mi , Mua
      REAL*8 sa , sb , sia , sib , W0 , Xiv
      DATA fct/0./

      mi = Mua - 1
      axi = ABS(Xiv) ! Absolute value of xi
      la = Lmda
      IF ( Lmda.EQ.7 ) la = 3

      IF ( axi.LT.1.E-5 ) THEN ! Small absolute values of xi
         a = -2.*W0
         IF ( la.EQ.3 ) a = -W0
         exa = EXP(a)
         dwi = 3*Dw
         cic = exa*(EXP(dwi)-1.)
         STAMP = DCMPLX(cic,0.D0)
         IF ( la.EQ.2 ) THEN
            IF ( mi.EQ.0 ) fct = 3.*(3.-Epsi*Epsi)/Epsi/Epsi/Epsi/Epsi
            IF ( mi.EQ.1 ) fct = 1.837117307*Errt/Epsi/Epsi/Epsi/Epsi ! 1.837117307 = 3/2 * sqrt(3/2)
            IF ( mi.EQ.2 ) fct = -3.674234613*Errt*Errt/Epsi/Epsi/ ! 3.674234613 = 3 * sqrt(3/2)
     &                           Epsi/Epsi
         ELSEIF ( la.EQ.3 ) THEN
            fct = -1.414213562*Errt/Epsi/Epsi ! 1.414213562 = sqrt(2)
         ELSE
            IF ( mi.EQ.0 ) fct = 1./Epsi/Epsi
            IF ( mi.EQ.1 ) fct = 1.414213562*Errt/Epsi/Epsi ! 1.414213562 = sqrt(2)
         ENDIF
      ELSE ! Larger absolute values of xi
         ex = EXP(W0)/2.
         b = axi*(Epsi*ex+W0)
         CALL TRINT(b,sib,cib)
         sb = SIN(b)/b
         cb = COS(b)/b
         bis = sb + cib
         bic = cb - sib
         bis2 = -sb/b
         bic2 = -cb/b
         dwi = -3.*Dw
         exa = EXP(dwi)
         a = axi*(Epsi*ex*exa+W0+dwi)
         sa = SIN(a)/a
         ca = COS(a)/a
         CALL TRINT(a,sia,cia)
         cis = sa + cia - bis
         cic = ca - sia - bic
         IF ( la.EQ.1 ) THEN
            STAMP = DCMPLX(cic,cis)
         ELSE
            dwi = (bic2-cis+ca/a)/2.
            exa = (bis2+cic+sa/a)/2.
            STAMP = DCMPLX(dwi,exa)
         ENDIF
         IF ( la.EQ.2 ) THEN
            IF ( mi.EQ.0 ) fct = .75*(3.-Epsi*Epsi)*axi*axi/Epsi/Epsi
            IF ( mi.EQ.1 ) fct = 1.837117307*Errt*axi*axi/Epsi/Epsi ! 1.837117307 = 3/2 * sqrt(3/2)
            IF ( mi.EQ.2 ) fct = -.9185586535*Errt*Errt*axi*axi/Epsi/ ! 0.9185586535 = 3/4 * sqrt(3/2)
     &                           Epsi
         ELSEIF ( la.EQ.3 ) THEN
           fct = -.3535533905*Errt*axi*axi ! 0.3535533907 = 1/4 * sqrt(2)
         ELSE
            IF ( mi.EQ.0 ) fct = .5*axi/Epsi
            IF ( mi.EQ.1 ) fct = .3535533907*Errt*axi/Epsi ! 0.3535533907 = 1/2 * sqrt(1/2)
         ENDIF
      ENDIF

      STAMP = STAMP*fct
      STAMP = CONJG(STAMP)
      END
 
C----------------------------------------------------------------------
C SUBROUTINE RESET
C
C Called by: INTG
C
C Purpose: to advance by one step. This means f(n-3) is set to the old value
C          of f(n-2), f(n-2) is set to the old value of f(n-1) and f(n-1) is
C          set to the old value of f(n).
C
C Uses global variables:
C      ARM    - excitation amplitudes of substates.
C      CAT    - substates of levels (n_level, J, m)
C      ISMAX  - number of substates used
C      NMAX   - number of levels
C      NSTART - index in CAT of first substate associated with a level
C
C Formal parameters:
C      Iso    - isotropic flag
 
      SUBROUTINE RESET(Iso)
      IMPLICIT NONE
      INTEGER*4 ir , Iso , j
      INTEGER*4 NSTART , NSTOP
      COMMON /CEXC0 / NSTART(101) , NSTOP(100)
      COMPLEX*16 ARM
      COMMON /AZ    / ARM(1200,7)
      REAL*8 CAT
      INTEGER*4 ISMAX
      COMMON /CLCOM8/ CAT(1200,3) , ISMAX
      INTEGER*4 NMAX , NDIM , NMAX1
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      
      IF ( Iso.EQ.0 ) THEN
         DO j = 1 , NMAX ! Loop over levels
            ir = NSTART(j) - 1 ! Index of first substate of level - 1
 20         ir = ir + 1
            ARM(ir,1) = ARM(ir,2)
            ARM(ir,2) = ARM(ir,3)
            ARM(ir,3) = ARM(ir,4)
            IF ( CAT(ir,3).LT.-.1 ) GOTO 20 ! m quantum number of substate ir
         ENDDO
         GOTO 99999
      ENDIF
       
      DO j = 1 , ISMAX ! Loop over substates
         ARM(j,1) = ARM(j,2)
         ARM(j,2) = ARM(j,3)
         ARM(j,3) = ARM(j,4)
      ENDDO
99999 END
 
C----------------------------------------------------------------------
C SUBROUTINE HALF
C
C Called by: INTG
C
C Purpose: to halve the step size for the integeration in INTG.
C
C Uses global variables:
C      ARM    - excitation amplitudes of substates.
C      CAT    - substates of levels (n_level, J, m)
C      ISMAX  - number of substates used
C      NMAX   - number of levels
C      NSTART - index in CAT of first substate associated with a level
C
C Formal parameters:
C      Iso    - isotropic flag
 
      SUBROUTINE HALF(Iso)
      IMPLICIT NONE
      INTEGER*4 ir , Iso , j
      COMPLEX*16 fpom
      INTEGER*4 NSTART , NSTOP
      COMMON /CEXC0 / NSTART(101) , NSTOP(100)
      COMPLEX*16 ARM
      COMMON /AZ    / ARM(1200,7)
      REAL*8 CAT
      INTEGER*4 ISMAX
      COMMON /CLCOM8/ CAT(1200,3) , ISMAX
      INTEGER*4 NMAX , NDIM , NMAX1
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      IF ( Iso.EQ.0 ) THEN
         DO j = 1 , NMAX ! Loop over levels
            ir = NSTART(j) - 1 ! Index of first substate of level - 1
 20         ir = ir + 1
            fpom = ARM(ir,3)
            ARM(ir,1) = -.0625*(ARM(ir,1)+ARM(ir,4))
     &                  + .5625*(ARM(ir,2)+ARM(ir,3))
            ARM(ir,3) = ARM(ir,3)*.75 + .375*ARM(ir,4) - ARM(ir,2)/8.
            ARM(ir,2) = fpom
            IF ( CAT(ir,3).LT.-.1 ) GOTO 20
         ENDDO
         GOTO 99999
      ENDIF
       
      DO j = 1 , ISMAX ! Loop over substates
         fpom = ARM(j,3)
         ARM(j,1) = -.0625*(ARM(j,4)+ARM(j,1))
     &              + .5625*(ARM(j,2)+ARM(j,3))
         ARM(j,3) = ARM(j,3)*.75 + .375*ARM(j,4) - ARM(j,2)/8.
         ARM(j,2) = fpom
      ENDDO
99999 END
 
C----------------------------------------------------------------------
C SUBROUTINE DOUBLE
C
C Called by: INTG
C
C Purpose: to double the step size for the integeration in INTG.
C
C Uses global variables:
C      ARM    - excitation amplitudes of substates.
C      CAT    - substates of levels (n_level, J, m)
C      ISMAX  - number of substates used
C      NMAX   - number of levels
C      NSTART - index in CAT of first substate associated with a level
C
C Formal parameters:
C      Iso    - isotropic flag
 
      SUBROUTINE DOUBLE(Iso)
      IMPLICIT NONE
      INTEGER*4 ir , Iso , j
      COMPLEX*16 fpom
      INTEGER*4 NSTART , NSTOP
      COMMON /CEXC0 / NSTART(101) , NSTOP(100)
      COMPLEX*16 ARM
      COMMON /AZ    / ARM(1200,7)
      REAL*8 CAT
      INTEGER*4 ISMAX
      COMMON /CLCOM8/ CAT(1200,3) , ISMAX
      INTEGER*4 NMAX , NDIM , NMAX1
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      
      IF ( Iso.EQ.0 ) THEN
         DO j = 1 , NMAX ! Loop over levels
            ir = NSTART(j) - 1 ! Index of first substate of level - 1
 20         ir = ir + 1
            fpom = ARM(ir,2)
            ARM(ir,2) = -8.*ARM(ir,3) + 6.*ARM(ir,2) + 3.*ARM(ir,4)
            ARM(ir,1) = -16.*ARM(ir,1) + 9.*ARM(ir,2) + 9.*fpom - 
     &                  ARM(ir,4)
            ARM(ir,3) = fpom
            IF ( CAT(ir,3).LT.-.1 ) GOTO 20
         ENDDO
         GOTO 99999
      ENDIF
       
      DO j = 1 , ISMAX ! Loop over substates
         fpom = ARM(j,2)
         ARM(j,2) = -8.*ARM(j,3) + 6.*ARM(j,2) + 3.*ARM(j,4)
         ARM(j,1) = -16.*ARM(j,1) + 9.*ARM(j,2) + 9.*fpom - ARM(j,4)
         ARM(j,3) = fpom
      ENDDO
99999 END
 
C----------------------------------------------------------------------
C SUBROUTINE PATH
C
C Called by: GOSIA
C
C Purpose: Calculate path for each level
C
C Uses global variables:
C      CAT    - substates of levels (n_level, J, m)
C      IPATH  - index of substate in level with same m as substate Irld
C      NMAX   - number of levels
C      NSTART - index in CAT of first substate associated with a level
C      NSTOP  - index in CAT of last substate associated with a level
C
C Formal parameters:
C      Irld   - index into ARM array
 
      SUBROUTINE PATH(Irld)
      IMPLICIT NONE
      REAL*8 spm , vl
      INTEGER*4 i , Irld , isp , ist , j
      INTEGER*4 NSTART , NSTOP
      COMMON /CEXC0 / NSTART(101) , NSTOP(100)
      INTEGER*4 IPATH , MAGA
      COMMON /PTH   / IPATH(100) , MAGA(100)
      INTEGER*4 NMAX , NDIM , NMAX1
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      REAL*8 CAT
      INTEGER*4 ISMAX
      COMMON /CLCOM8/ CAT(1200,3) , ISMAX
      spm = CAT(Irld,3) ! m quantum number for substate Irld
      DO i = 2 , NMAX ! For each level except ground state
         IPATH(i) = 0
         ist = NSTART(i) ! Index of first substate for level
         IF ( ist.NE.0 ) THEN ! If this is non-zero
            isp = NSTOP(i) ! Index of last substate for level
            DO j = ist , isp ! For each substate of level
               vl = CAT(j,3) ! m quantum number for substate j
               IF ( ABS(vl-spm).LT.1.E-6 ) GOTO 50 ! Jump if they have the same m
            ENDDO
         ENDIF
         GOTO 100
 50      IPATH(i) = j ! Store it
 100     CONTINUE
      ENDDO
      IPATH(1) = Irld ! Special case of ground state
      END
 
C----------------------------------------------------------------------
C SUBROUTINE INTG
C
C Called by: FTBM, GOSIA
C Calls:     AMPDER, DOUBLE, HALF, RESET
C
C Purpose: the main integration routine.
C
C Uses global variables:
C      ACC50  - accuracy required for integration
C      ARM    - excitation amplitudes of substates.
C      CAT    - substates of levels (n_level, J, m)
C      D2W    - step in omega (= 0.03)
C      IFAC   - spin/parity phase factor
C      IFLG   - flag to determine whether to calculate exponential (so we don't calculate twice)
C      INTERV - default accuracy check parameter (see OP,CONT:INT)
C      IPATH  - index of substate in level with same m as substate Irld
C      IRA    - limit of omega for integration for each multipolarity
C      ISG    - sign of omega
C      ISMAX  - number of substates used
C      ISO    - Isotropic flag
C      KDIV   - index for division
C      LAMR   - flag = 1 if we should calculate this multipolarity
C      MAXLA  - multipolarity to calculate
C      NDIV   - number of divisions
C      NMAX   - number of levels
C      NPT    - index in ADB array (this is omega / 0.03)
C      NSTART - index in CAT of first substate associated with a level
C      NSW    - step in omega
C
C Formal parameters:
C      Ien    - experiment number
C
C Note that if it finds that the step size for the integral is too small, it
C calls DOUBLE to increase it by a factor of two, or if it finds that the
C step size is too big, it decreases it by a factor of two by calling HALF.
C
C We use the the 4th order Adams-Moulton predictor-corrector method for
C solving an ordinary differential equation. We use an adaptive version, which
C can change the step size (increase or decrease) in order to get the desired
C accuracy.
C
C The predictor is given as:
C
C y(n+1)_p = y(n) + h/24 * {55*f(n) - 59*f(n-1) + 37*f(n-2) - 9*f(n-3)}
C
C and the corrector is:
C
C y(n+1)_c = y(n) + h/24 * {9*f_p(n+1) + 19*f(n) - 5*f(n-1) + f(n-2)}
C
C The error is |E(n+1)| ~ 19/270 * {y_p(n+1) - y_c(n+1)}
C
C In this function:
C                   D2W        = h
C                   ARM(ir, 1) = f(n-3)
C                   ARM(ir, 2) = f(n-2)
C                   ARM(ir, 3) = f(n-1)
C                   ARM(ir, 4) = f(n)
C                   ARM(ir, 5) = y(n) initially
C                   ARM(ir, 5) = y_c(n+1) finally
C                   ARM(ir, 6) is not used
C                   ARM(ir, 7) = y_p(n+1)
C
C The function RESET is called to advance n by one. i.e. f(n-3) is set to the
C old value of f(n-2), f(n-2) to the old value of f(n-1) and f(n-1) to the old
C value of f(n).

 
      SUBROUTINE INTG(Ien)
      IMPLICIT NONE
      REAL*8 f , rim , rl , srt
      INTEGER*4 i , i57 , Ien , ihold , intend , ir , ir1 , k , kast , 
     &          mir , n
      COMPLEX*16 hold
      REAL*8 EN , SPIN , ACCUR , DIPOL , ZPOL , ACCA
      INTEGER*4 ISO
      COMMON /COEX  / EN(100) , SPIN(100) , ACCUR , DIPOL , ZPOL , 
     &                ACCA , ISO
      COMPLEX*16 ARM
      COMMON /AZ    / ARM(1200,7)
      INTEGER*4 IRA , MAXLA
      COMMON /RNG   / IRA(8) , MAXLA
      REAL*8 ACC50
      COMMON /A50   / ACC50
      INTEGER*4 IFAC
      COMMON /CLCOM0/ IFAC(100)
      REAL*8 CAT
      INTEGER*4 ISMAX
      COMMON /CLCOM8/ CAT(1200,3) , ISMAX
      INTEGER*4 NMAX , NDIM , NMAX1
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      REAL*8 D2W
      INTEGER*4 NPT , NDIV , KDIV , LAMR , ISG , NSW , ISG1
      COMMON /CAUX  / NPT , NDIV , KDIV , LAMR(8) , ISG , D2W , NSW , 
     &                ISG1
      INTEGER*4 IFLG
      COMMON /FLA   / IFLG
      INTEGER*4 NSTART , NSTOP
      COMMON /CEXC0 / NSTART(101) , NSTOP(100)
      INTEGER*4 IPATH , MAGA
      COMMON /PTH   / IPATH(100) , MAGA(100)
      INTEGER*4 INTERV
      COMMON /CEXC9 / INTERV(50)
      
      real*8 adamtemp,TCABS                  ! Rachel modification
c      logical diderrcheck                    ! Rachel modification
      integer*4 MEM                          ! Rachel modification
      integer*4 indx                         ! Rachel modification
      COMPLEX*16 EXPO
      COMMON /ADBXI / EXPO(1500)
      INTEGER*4 IPRM
      COMMON /PRT   / IPRM(20)
c-------ADDITIONAL OUTPUT TO BE USED BY RACHEL.PY. MAR. 16 2011------------
      ! Rachel modification: print the collision functions to unit 99
      IF (IPRM(9).LT.0) CALL SPITQ(Ien,ABS(IPRM(9)))
      IF (IPRM(9).eq.11) THEN                ! If option was to print exc. amp. of substates
c     Write a blank line to mark the start of the experiment.
        WRITE(99,14617) 
14617   FORMAT("   ")
      END IF
c-------END OF ADDITIONAL RACHEL OUTPUT.-----------------------------------
      
      intend = INTERV(Ien) ! Default accuracy set by INT option of OP,CONT
      D2W = .03 ! We use steps of 0.03 in omega
      NSW = 1
      kast = 0
      NDIV = 0
      KDIV = 0
 100  IF ( (NPT+NSW).GT.IRA(MAXLA) .AND. ISG.GT.0 ) RETURN
      DO i = 1 , 8
         LAMR(i) = 0
         IF ( (NPT+NSW).LT.IRA(i) ) LAMR(i) = 1
      ENDDO
C     Predictor 
      IF ( ISO.EQ.0 ) THEN
         DO n = 1 , NMAX ! For each level
            ir = NSTART(n) - 1 ! First substate - 1
 120        ir = ir + 1
            ARM(ir,7) = ARM(ir,5)
     &                  + D2W/24.*(55.0*ARM(ir,4)-59.0*ARM(ir,3)
     &                  +37.0*ARM(ir,2)-9.0*ARM(ir,1))
            mir = CAT(ir,3) ! m quantum number of substate ir
            ir1 = ir - 2*mir
            ARM(ir1,7) = IFAC(n)*ARM(ir,7)
            IF ( DBLE(mir).LT.-0.1 ) GOTO 120
         ENDDO
      ELSE
         DO ir = 1 , ISMAX
            ARM(ir,7) = ARM(ir,5)
     &                  + D2W/24.*(55.0*ARM(ir,4)-59.0*ARM(ir,3)
     &                  +37.0*ARM(ir,2)-9.0*ARM(ir,1))
         ENDDO
      ENDIF
      NPT = NPT + NSW*ISG ! NPT loops over omega values, ISG is -1 at first then +1
      IF ( NPT.GT.0 ) THEN
         IF ( NDIV.EQ.0 ) GOTO 200
         KDIV = KDIV + 1
         IF ( KDIV.LT.NDIV ) GOTO 200
         KDIV = 0
         NPT = NPT + ISG
         IF ( NPT.GT.0 ) GOTO 200
      ENDIF
      NPT = -NPT + 2 ! We decreased omega to zero, so now start increasing
      ISG = 1
 200  CALL RESET(ISO)
      IFLG = 1
      i57 = 7 ! Tell LAISUM to use ARM(I,7) for excitation amplitudes

C     Calculate derivatives of amplitudes
      CALL AMPDER(i57)
      
C     Corrector
      IF ( ISO.EQ.0 ) THEN
         DO n = 1 , NMAX ! For each level
            ir = NSTART(n) - 1 ! First substate - 1
 220        ir = ir + 1
            ARM(ir,5) = ARM(ir,5)
     &                  + D2W/24.*(9.0*ARM(ir,4)+19.0*ARM(ir,3)
     &                  -5.0*ARM(ir,2)+ARM(ir,1))
            mir = CAT(ir,3) ! m quantum number of substate ir
            ir1 = ir - 2*mir
            ARM(ir1,5) = IFAC(n)*ARM(ir,5)
            IF ( DBLE(mir).LT.-0.1 ) GOTO 220
         ENDDO
      ELSE
         DO ir = 1 , ISMAX ! For each substate
            ARM(ir,5) = ARM(ir,5)
     &                  + D2W/24.*(9.0*ARM(ir,4)+19.0*ARM(ir,3)
     &                  -5.0*ARM(ir,2)+ARM(ir,1))
         ENDDO
      ENDIF
      kast = kast + 1
      IFLG = 0
      i57 = 5 ! Tell LAISUM to use ARM(I,5) for excitation amplitudes

C     Calculate derivatives of amplitudes
      CALL AMPDER(i57)
      IF ( (LAMR(2)+LAMR(3)).NE.0 ) THEN
         IF ( kast.GE.intend ) THEN
            kast = 0
            f = 0.
            DO k = 1 , NMAX ! For each level
               ihold = IPATH(k)
               IF ( ihold.NE.0 ) THEN
                  hold = ARM(ihold,5) - ARM(ihold,7)
                  rl = DBLE(hold)
                  rim = DIMAG(hold)
                  srt = rl*rl + rim*rim
                  f = MAX(f,srt)
               ENDIF
            ENDDO

C           Decide if we have appropriate accuracy (strictly it should be
C           f = SQRT(f)*19./270. but the difference is not all that large).
C
            f = SQRT(f)/14.
            IF ( f.GT.ACCUR .OR. f.LT.ACC50 ) THEN
               IF ( f.LT.ACC50 ) THEN
                  CALL DOUBLE(ISO) ! Double step size
                  D2W = 2.*D2W
                  NSW = 2*NSW
                  intend = (DBLE(intend)+.01)/2.
                  IF ( intend.EQ.0 ) intend = 1
                  IF ( NSW.LT.1 ) THEN
                     NDIV = (DBLE(NDIV)+.01)/2.
                     IF ( NDIV.LT.2 ) THEN
                        NDIV = 0
                        NSW = 1
                     ENDIF
                  ENDIF
               ELSE
                  CALL HALF(ISO) ! Halve step size
                  D2W = D2W/2.
                  NSW = (DBLE(NSW)+.01)/2.
                  intend = 2*intend
                  IF ( NSW.LT.1 ) THEN
                     NDIV = 2*NDIV
                     IF ( NDIV.EQ.0 ) NDIV = 2
                  ENDIF
               ENDIF
            ENDIF
             
         ENDIF ! if kast>=intend
      ENDIF

c-------ADDITIONAL OUTPUT TO BE USED BY RACHEL.PY. MAR. 16 2011------------
c     I am trying to output the probabilities and amplitudes at each step

      if(IPRM(9).eq.11) then                ! If option was to print exc. amp. of substates
        write(99,14619) NPT,D2W
14619   format(2X,I4,2x,f5.3,$)
c       if(diderrcheck) then 
c         write(99,14621) sqrt(f)/14.              ! print out the error term 
c14621     format(2x,E11.4,$)
c       else
c         write(99,14622)                    !  if no error (f) term calc'd
c14622     format(' none',$)                  ! error wasn't checked
c       end if
      end if
      if((IPRM(9).GT.0).and.(IPRM(9).le.6)) then ! if option was to print adiab exp
c       Note that it looks up the the terms by multipolarity
c       so I select it by lambda = IPRM(9)
        indx = MEM(1,2,IPRM(9))                  ! Index for matrix element from level N to level m with multipolarity La
        write(99,14699)DBLE(EXPO(indx)),DIMAG(EXPO(indx))
14699   format(2x,'2',2x,D11.4,2x,D11.4,$)
        write(99,14623)    ! close the line
      else if(IPRM(9).eq.11) then                 ! if option was to print excitation amplitudes of substates  
        do ir = 1, ismax
          adamtemp = TCABS(ARM(ir,i57))**2    ! probability(step)
          write(99,14618)ir,DBLE(ARM(ir,i57)),DIMAG(ARM(ir,i57)),
     &                   adamtemp
        enddo
14618   format(2x,I4,2x,D11.4,2x,D11.4,2x,D11.4,2x,$)
        write(99,14623)    ! close the line
      endif


14623 format('')
c-------END OF ADDITIONAL RACHEL OUTPUT.-----------------------------------


      GOTO 100
      END
 
C----------------------------------------------------------------------
C SUBROUTINE NEWLV
C
C Called by: AMPDER, STING
C Calls:     EXPON, LEADF, MEM
C
C Purpose: Setup a new level which can be excited from ground state. We store
C       ISSTAR, ISSTO and MSTORE for the level and calculate and store the
C       exponential: exp(i \xi_{kn} (\epsilon \sinh(\omega) + \omega))
C
C Uses global variables:
C      EXPO   - adiabatic exponential
C      IFLG   - flag to determine whether to calculate exponential (so we don't calculate twice)
C      ISG    - sign of omega
C      ISG1   - index of omega
C      ISSTAR - index of last substate for that level
C      ISSTO  - index of first substate for that level
C      KDIV   - index for division
C      LDNUM  - number of matrix elements with each multipolarity populating level
C      MSTORE - index of final level number and index of matrix element
C      NDIV   - number of divisions
C      NPT    - index in ADB array (this is omega / 0.03)
C      NSTART - index in CAT of first substate associated with a level
C      NSTOP  - index in CAT of last substate associated with a level
C
C Formal parameters:
C      N      - level number
C      Ld     - Number of matrix elements for level N multipolarity La
C      La     - multipolarity
C
C Note that the exponential is calculated by EXPON. This file does the
C storage part.
      
      SUBROUTINE NEWLV(N,Ld,La)
      IMPLICIT NONE
      INTEGER*4 i2 , indx , La , Ld , LEADF , m , MEM , N
      COMPLEX*16 EXPON
      INTEGER*4 LAMDA , LEAD , LDNUM , LAMMAX , MULTI
      COMMON /CLCOM / LAMDA(8) , LEAD(2,1500) , LDNUM(8,100) , LAMMAX , 
     &                MULTI(8)
      REAL*8 D2W
      INTEGER*4 NPT , NDIV , KDIV , LAMR , ISG , NSW , ISG1
      COMMON /CAUX  / NPT , NDIV , KDIV , LAMR(8) , ISG , D2W , NSW , 
     &                ISG1
      INTEGER*4 ISSTAR , ISSTO , MSTORE
      COMMON /PINT  / ISSTAR(101) , ISSTO(100) , MSTORE(2,100)
      COMPLEX*16 EXPO
      COMMON /ADBXI / EXPO(1500)
      INTEGER*4 IFLG
      COMMON /FLA   / IFLG
      INTEGER*4 NSTART , NSTOP
      COMMON /CEXC0 / NSTART(101) , NSTOP(100)
      Ld = LDNUM(La,N) ! Get number of levels connected to level N by multipolarity La
      IF ( Ld.EQ.0 ) RETURN ! Return if there aren't any

      DO i2 = 1 , Ld ! For each level
         m = LEADF(N,i2,La) ! Get the other level associated
         ISSTAR(i2) = NSTOP(m) ! Get the index of last substate for that level
         ISSTO(i2) = NSTART(m) ! Get the index of first substate for that level
         MSTORE(1,i2) = m ! Store the final level number
         indx = MEM(N,m,La) ! Index for matrix element from level N to level m with multipolarity La
         MSTORE(2,i2) = indx ! Store index of matrix element
         IF ( IFLG.NE.0 ) THEN
            IF ( m.NE.N ) EXPO(indx) = EXPON(indx,NPT,ISG,ISG1,NDIV,KDIV
     &                                 )
         ENDIF
      ENDDO
      END
 
C----------------------------------------------------------------------
C SUBROUTINE CODE7
C
C Called by: LSLOOP
C
C Purpose:
C
C Uses global variables:
C      IAPR   - index of initial and final levels for each matrix element
C      IPATH  - index of substate in level with same m as substate Irld
C
C Formal parameters:
C      Ir     - index of initial substate
C      Is     - index of final substate
C      N      - index of initial level
C      Mt     - index of final level
C      Inqa   - result of operation
C      Indx   - Index of matrix element
 
      SUBROUTINE CODE7(Ir,Is,N,Mt,Inqa,Indx)
      IMPLICIT NONE
      INTEGER*4 idm , idn , Indx , Inqa , Ir , Is , ism , Mt , N
      INTEGER*4 IPATH , MAGA
      COMMON /PTH   / IPATH(100) , MAGA(100)
      REAL*8 QAPR
      INTEGER*4 IAPR , ISEX
      COMMON /APRCAT/ QAPR(1500,2,7) , IAPR(1500,2) , ISEX(100)
      
      IAPR(Indx,1) = N  ! Index of initial level
      IAPR(Indx,2) = Mt ! Index of final level
      IF ( IPATH(N).EQ.0 .OR. IPATH(Mt).EQ.0 ) THEN
         Inqa = -1
         GOTO 99999
      ELSE
         idn = Ir - IPATH(N)
         idm = Is - IPATH(Mt)
         ism = idn + idm + 3
         IF ( ism.EQ.2 ) THEN
            Inqa = 2
            IF ( idn.GT.idm ) Inqa = 3
            RETURN
         ELSEIF ( ism.EQ.3 ) THEN
            Inqa = 4
            RETURN
         ELSEIF ( ism.EQ.4 ) THEN
            Inqa = 5
            IF ( idn.GT.idm ) Inqa = 6
            RETURN
         ELSEIF ( ism.NE.5 ) THEN
            Inqa = 1
            RETURN
         ENDIF
      ENDIF
      Inqa = 7
99999 END
 
C----------------------------------------------------------------------
C SUBROUTINE APRAM
C
C Called by: FTBM
C Calls:     NEWCAT, PODZIEL, POMNOZ
C
C Purpose: calculate approximate value of the Coulomb excitation amplitudes.
C
C Uses global parameters:
C      ARM    - excitation amplitudes of substates.
C      ELM    - matrix elements
C      IDIVE  - number of subdivisions
C      LERF   - error flag for expansion in POMNOZ
C      MAGA   - number of magnetic substates in approximate calculation
C      MEMX6  - number of matrix elements with E1...6 multipolarity
C      QAPR   - approximate Coulomb amplitudes
C
C Formal parameters:
C      Iexp   - experiment number
C      Inc    - flag: first time we call after LOAD, Inc=0, afterwards Inc=1
C      Indx   - index of matrix element
C      Irld   - index into ARM array
C      Acca   - accuracy required

      SUBROUTINE APRAM(Iexp,Inc,Indx,Irld,Acca)
      IMPLICIT NONE
      REAL*8 Acca , accah , uwa
      INTEGER*4 i1 , i56 , i7 , Iexp , img , Inc , Indx , Irld , itm , 
     &          j , jidim , jj , k , ktoto , l , l1 , l2 , l3 , m
      COMPLEX*16 ARM
      COMMON /AZ    / ARM(1200,7)
      REAL*8 QAPR
      INTEGER*4 IAPR , ISEX
      COMMON /APRCAT/ QAPR(1500,2,7) , IAPR(1500,2) , ISEX(100)
      INTEGER*4 IPATH , MAGA
      COMMON /PTH   / IPATH(100) , MAGA(100)
      INTEGER*4 MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR
      COMMON /CEXC  / MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR(1500)
      REAL*8 ELM , ELMU , ELML , SA
      COMMON /COMME / ELM(1500) , ELMU(1500) , ELML(1500) , SA(1500)
      INTEGER*4 LERF , IDIVE
      COMMON /APRX  / LERF , IDIVE(50,2)
      LERF = 0
      accah = Acca
 100  i7 = 7
      itm = -1
      img = 3
      i1 = 1
      IF ( MAGA(Iexp).EQ.0 ) THEN
         i7 = 4
         i1 = 4
         img = 1
      ENDIF
      IF ( Inc.EQ.0 ) GOTO 300
      IF ( LERF.EQ.0 ) CALL NEWCAT(Iexp,jidim)
      IF ( LERF.EQ.0 ) CALL PODZIEL(3,Iexp) ! Subdivide
      i56 = 5
      DO k = 1 , jidim
         ARM(k,2) = (0.,0.)
         ARM(k,5) = (0.,0.)
      ENDDO
      ARM(Irld+1,5) = (1.,0.)

 200  ktoto = 0
      LERF = 0

      l1 = IDIVE(Iexp,1)
      DO l3 = 1 , l1
         Acca = accah*l3/l1
         CALL POMNOZ(Acca,1,i56,ktoto,img,jidim) ! Expansion for L=1
         IF ( LERF.NE.0 ) THEN
            CALL PODZIEL(1,Iexp) ! Subdivide
            GOTO 100
         ENDIF
      ENDDO

      l2 = IDIVE(Iexp,2)
      DO l3 = 1 , l2
         Acca = accah + accah*l3/l2
         CALL POMNOZ(Acca,2,i56,ktoto,img,jidim) ! Expansion for L=2
         IF ( LERF.NE.0 ) THEN
            CALL PODZIEL(2,Iexp) ! Subdivide
            GOTO 100
         ENDIF
      ENDDO

      DO l = 1 , MEMX6 ! Matrix elements for E1...6
         DO m = i1 , i7
            QAPR(l,1,m) = -QAPR(l,1,m)
         ENDDO
      ENDDO

      DO l3 = 1 , l1
         Acca = accah*2. + accah*l3/l1
         CALL POMNOZ(Acca,1,i56,ktoto,img,jidim) ! Expansion for L=1
      ENDDO

      Acca = accah
      DO l = 1 , MEMX6 ! Matrix elements for E1...6
         DO m = i1 , i7
            QAPR(l,1,m) = -QAPR(l,1,m)
         ENDDO
      ENDDO

      IF ( Inc.NE.0 .OR. itm.NE.0 ) THEN
         IF ( Inc.EQ.0 ) THEN
            DO l = 1 , jidim
               ARM(l,6) = ARM(l,6) - ARM(l,7)
               ARM(l,6) = 50.*ARM(l,6)/ELM(Indx)
            ENDDO
            DO l = 1 , 2
               DO j = i1 , i7
                  QAPR(Indx,l,j) = QAPR(Indx,l,j)/.99
               ENDDO
            ENDDO
            DO jj = 2 , jidim
               ARM(jj-1,6) = ARM(jj,6)
            ENDDO
            GOTO 99999
         ELSE
            DO jj = 2 , jidim
               ARM(jj-1,5) = ARM(jj,5)
            ENDDO
            RETURN
         ENDIF
      ENDIF

C     Initialise (Inc = 0)
 300  itm = itm + 1
      i56 = itm + 6
      DO k = 1 , jidim
         ARM(k,i56) = (0.,0.)
      ENDDO

      ARM(Irld+1,i56) = (1.,0.)
      uwa = -itm*.0298019802 + 1.01
      DO l = 1 , 2
         DO j = i1 , i7
            QAPR(Indx,l,j) = QAPR(Indx,l,j)*uwa
         ENDDO
      ENDDO

      DO j = 1 , jidim
         ARM(j,2) = (0.,0.)
      ENDDO
      GOTO 200

99999 END
 
C----------------------------------------------------------------------
C SUBROUTINE NEWCAT
C
C Called by: APRAM
C Calls:     FXIS1, FXIS2
C
C Purpose: create a new catalog of matrix elements
C
C Uses global variables:
C      IAPR   - index of initial and final levels for each matrix element
C      MAGA   - number of magnetic substates in approximate calculation
C      MULTI  - number of matrix elements having given multipolarity
C      NMAX   - number of levels
C      PARX   -
C      PARXM  -
C      QAPR   - approximate Coulomb amplitudes
C      XI     - xi coupling coefficients
C      XIR    - [for maps]
C
C Formal parameters:
C     Iexp    - experiment number
C     Jidim   - 
      
      SUBROUTINE NEWCAT(Iexp,Jidim)
      IMPLICIT NONE
      REAL*8 a , b , FXIS1 , FXIS2 , q1 , q2 , wg , wl , xp , xx , zt
      INTEGER*4 Iexp , ist , istop , Jidim , k , kk , n , 
     &          ng , nl
      REAL*8 PARX, PARXM, XIR
      COMMON /MAP   / PARX(50,12,5) , PARXM(50,4,10,6) , XIR(6,50)
      REAL*8 XI
      COMMON /CXI   / XI(1500)
      INTEGER*4 LAMDA , LEAD , LDNUM , LAMMAX , MULTI
      COMMON /CLCOM / LAMDA(8) , LEAD(2,1500) , LDNUM(8,100) , LAMMAX , 
     &                MULTI(8)
      INTEGER*4 NMAX , NDIM , NMAX1
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      INTEGER*4 IPATH , MAGA
      COMMON /PTH   / IPATH(100) , MAGA(100)
      REAL*8 QAPR
      INTEGER*4 IAPR , ISEX
      COMMON /APRCAT/ QAPR(1500,2,7) , IAPR(1500,2) , ISEX(100)
      Jidim = NMAX + 1
      IF ( MAGA(Iexp).NE.0 ) Jidim = 3*NMAX + 1
      ist = 1
      DO kk = 1 , 6 ! For each multipolarity E1 to E6
         IF ( MULTI(kk).NE.0 ) THEN ! If there are matrix elements for this multipolarity
            istop = MULTI(kk) - 1 + ist ! Last matrix element for this multipolarity
            DO k = ist , istop ! For each matrix element for this multipolarity
               xx = ABS(XI(k))
               xx = xx/XIR(kk,Iexp)

               DO n = 1 , 7 , 3 ! For n in 1, 4 and 7
                  IF ( MAGA(Iexp).NE.0 .OR. n.EQ.4 ) THEN
                     zt = QAPR(k,1,n)
                     zt = ABS(zt)
                     xp = 9.*xx
                     nl = INT(xp) + 1
                     wg = xp - DBLE(nl-1)
                     ng = nl + 1
                     wl = DBLE(nl) - xp
                     a = wg*PARXM(Iexp,1,ng,kk) + wl*PARXM(Iexp,1,nl,kk)
                     b = wg*PARXM(Iexp,2,ng,kk) + wl*PARXM(Iexp,2,nl,kk)
                     q1 = a*zt + b
                     a = wg*PARXM(Iexp,3,ng,kk) + wl*PARXM(Iexp,3,nl,kk)
                     b = wg*PARXM(Iexp,4,ng,kk) + wl*PARXM(Iexp,4,nl,kk)
                     q2 = a*zt + b
                     QAPR(k,2,n) = QAPR(k,1,n)*q2*FXIS2(k,n) ! FXIS2 = 1
                     QAPR(k,1,n) = QAPR(k,1,n)*q1*FXIS1(k,n) ! FXIS1 = sign(xi)
                     IF ( IAPR(k,1).EQ.IAPR(k,2) ) THEN
                        QAPR(k,1,n) = 0.
                        QAPR(k,2,n) = QAPR(k,2,n)/2.
                     ENDIF
                  ENDIF
               ENDDO ! Loop over n

               IF ( MAGA(Iexp).NE.0 ) THEN
                  DO n = 2 , 6
                     IF ( n.NE.4 ) THEN ! For N in 2, 3, 5, 6
                        zt = QAPR(k,1,n)
                        zt = ABS(zt)
                        xp = 4.*xx
                        nl = INT(xp) + 1
                        wg = xp - DBLE(nl-1)
                        ng = nl + 1
                        wl = DBLE(nl) - xp
                        q1 = wg*PARX(Iexp,2*kk-1,ng)
     &                       + wl*PARX(Iexp,2*kk-1,nl)
                        q2 = wg*PARX(Iexp,2*kk,ng)
     &                       + wl*PARX(Iexp,2*kk,nl)
                        QAPR(k,2,n) = QAPR(k,1,n)*q2*FXIS2(k,n) ! FXIS2 = sign(xi)
                        QAPR(k,1,n) = QAPR(k,1,n)*q1*FXIS1(k,n) ! FXIS1 = 1
                     ENDIF
                  ENDDO ! Loop over n
               ENDIF
            ENDDO ! Loop over matrix elements for this mult. k
            ist = istop + 1
         ENDIF
      ENDDO ! Loop over electric multipolarities kk
      END
 
C----------------------------------------------------------------------
C SUBROUTINE POMNOZ
C
C Called by: APRAM
C Calls:     TCABS
C
C Purpose: perform the expansion to calculate the approximate Coulomb
C          amplitudes
C
C Note: pomnoz is Polish for multiply
C
C Uses global variables:
C      ARM    - excitation amplitudes of substates.
C      IAPR   - index of initial and final levels for each matrix element
C      INHB   - inhibit error flag setting (LERF)
C      IPATH  - index of substate in level with same m as substate Irld
C      ISEX   -
C      LERF   - error flag which is set here and used in APRAM
C      MEMX6  - number of matrix elements with E1...6 multipolarity
C      QAPR   - approximate Coulomb amplitudes
C
C Formal parameters:
C      Acca   - accuracy required
C      L
C      Iw
C      Img
C      Jidim
C      Ktoto  - number of iterations needed
 
      SUBROUTINE POMNOZ(Acca,L,Iw,Ktoto,Img,Jidim)
      IMPLICIT NONE
      REAL*8 Acca , sig , TCABS , test , u
      INTEGER*4 Img , Iw , Jidim , k , kk , Ktoto , L , m , mc , mc1 , 
     &          mw , mw1
      COMPLEX*16 ci
      INTEGER*4 INHB
      COMMON /INHI  / INHB
      REAL*8 QAPR
      INTEGER*4 IAPR , ISEX
      COMMON /APRCAT/ QAPR(1500,2,7) , IAPR(1500,2) , ISEX(100)
      INTEGER*4 IPATH , MAGA
      COMMON /PTH   / IPATH(100) , MAGA(100)
      INTEGER*4 MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR
      COMMON /CEXC  / MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR(1500)
      COMPLEX*16 ARM
      COMMON /AZ    / ARM(1200,7)
      INTEGER*4 LERF , IDIVE
      COMMON /APRX  / LERF , IDIVE(50,2)
      DATA ci/(0.,-1.)/ ! -sqrt(-1)

      sig = 1.
      IF ( L.NE.2 ) sig = -1.
      DO kk = 1 , Jidim
         ARM(kk,1) = ARM(kk,Iw)
      ENDDO

      DO k = 1 , 100 ! Perform up to 100 iterations
         Ktoto = Ktoto + 1
         DO m = 1 , MEMX6 ! Matrix elements for E1...6
            mw1 = IAPR(m,1)
            mc1 = IAPR(m,2)
            IF ( IPATH(mw1).NE.0 .AND. IPATH(mc1).NE.0 ) THEN
               mw = IPATH(mw1) + 1
               mc = IPATH(mc1) + 1
               IF ( Ktoto.GE.ISEX(mc1) ) THEN
                  IF ( Img.EQ.1 ) THEN
                     ARM(mw,2) = ARM(mw,2) + QAPR(m,L,4)*ARM(mc,1)
                     ARM(mc,2) = ARM(mc,2) + sig*QAPR(m,L,4)*ARM(mw,1)
                  ELSE
                     ARM(mw,2) = ARM(mw,2) + QAPR(m,L,4)*ARM(mc,1)
                     ARM(mc,2) = ARM(mc,2) + sig*QAPR(m,L,4)*ARM(mw,1)
                     ARM(mw-1,2) = ARM(mw-1,2) + QAPR(m,L,2)*ARM(mc,1)
                     ARM(mc,2) = ARM(mc,2) + sig*QAPR(m,L,2)*ARM(mw-1,1)
                     ARM(mw-1,2) = ARM(mw-1,2) + QAPR(m,L,1)*ARM(mc-1,1)
                     ARM(mc-1,2) = ARM(mc-1,2) + sig*QAPR(m,L,1)
     &                             *ARM(mw-1,1)
                     ARM(mw,2) = ARM(mw,2) + QAPR(m,L,3)*ARM(mc-1,1)
                     ARM(mc-1,2) = ARM(mc-1,2) + sig*QAPR(m,L,3)
     &                             *ARM(mw,1)
                     ARM(mw,2) = ARM(mw,2) + QAPR(m,L,5)*ARM(mc+1,1)
                     ARM(mc,2) = ARM(mc,2) + sig*QAPR(m,L,6)*ARM(mw+1,1)
                     ARM(mc+1,2) = ARM(mc+1,2) + sig*QAPR(m,L,5)
     &                             *ARM(mw,1)
                     ARM(mw+1,2) = ARM(mw+1,2) + QAPR(m,L,6)*ARM(mc,1)
                     ARM(mw+1,2) = ARM(mw+1,2) + QAPR(m,L,7)*ARM(mc+1,1)
                     ARM(mc+1,2) = ARM(mc+1,2) + sig*QAPR(m,L,7)
     &                             *ARM(mw+1,1)
                  ENDIF
               ENDIF
            ENDIF
         ENDDO

C        Calculate accuracy we have achieved
         test = 0.
         DO m = 1 , Jidim
            ARM(m,1) = ARM(m,2)/k
            ARM(m,2) = (0.,0.)
            IF ( L.NE.1 ) ARM(m,1) = ARM(m,1)*ci
            ARM(m,Iw) = ARM(m,Iw) + ARM(m,1)
            IF ( k.GT.5 ) THEN
               u = TCABS(ARM(m,Iw))
               test = test + u*u
            ENDIF
         ENDDO
C        Test to see if we have achieved required accuracy
         IF ( ABS(test-1.).LT.Acca ) GOTO 99999 ! Accuracy OK, so end
      ENDDO ! Iteration loop

      IF ( INHB.NE.1 ) LERF = 1 ! Accuracy not achieved, so set error flag
99999 END
 
C----------------------------------------------------------------------
C SUBROUTINE TENB
C
C Called by: FTBM, GOSIA
C Calls:     WTHREJ
C
C Purpose: calculate the state of polarization of the decaying level
C
C Uses global variables:
C      ARM    - excitation amplitudes of substates.
C      CAT    - substates of levels (n_level, J, m)
C      NMAX   - number of levels
C      NSTART - index in CAT of first substate associated with a level
C      NSTOP  - index in CAT of last substate associated with a level
C      SPIN   - spin of level
C
C Formal parameters:
C      Icl    - multipolarity
C      Bten   - result
C      Lmax   - maximum multipolarity to calculate for
C
C Note that the parameters to WTHREJ are all doubled, so that this routine
C can cope with half-integers.

      SUBROUTINE TENB(Icl,Bten,Lmax)
      IMPLICIT NONE
      REAL*8 Bten , ce , fc , si , WTHREJ , x
      INTEGER*4 i , Icl , iha , ila , ilg , ind , isi , ite , jm , 
     &          jmp , k , kk , kp , l , ll , Lmax , lp , m
      INTEGER*4 mm , mp , ms , msp
      DIMENSION Bten(*)
      REAL*8 EN , SPIN , ACCUR , DIPOL , ZPOL , ACCA
      INTEGER*4 ISO
      COMMON /COEX  / EN(100) , SPIN(100) , ACCUR , DIPOL , ZPOL , 
     &                ACCA , ISO
      REAL*8 CAT
      INTEGER*4 ISMAX
      COMMON /CLCOM8/ CAT(1200,3) , ISMAX
      INTEGER*4 NMAX , NDIM , NMAX1
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      INTEGER*4 NSTART , NSTOP
      COMMON /CEXC0 / NSTART(101) , NSTOP(100)
      COMPLEX*16 ARM
      COMMON /AZ    / ARM(1200,7)
      
      iha = (-1)**INT(2.*SPIN(1)+.01)
      IF ( Icl.EQ.1 ) THEN
         ms = 16*(NMAX-1)
         DO i = 1 , ms
            Bten(i) = 0.
         ENDDO
      ENDIF

      DO i = 2 , NMAX ! For each level except ground state
         ms = NSTART(i) ! First substate of level
         IF ( ms.NE.0 ) THEN
            msp = NSTOP(i) ! Last substate of level
            si = SPIN(i) ! Spin of level
            isi = INT(2.*si+.01)
            ce = SQRT(2.*si+1.)
            DO kp = 1 , 7 , 2
               k = kp - 1
               kk = 2*k
               IF ( isi.GE.k ) THEN
                  ila = -1
                  DO lp = 1 , kp
                     ila = -ila
                     l = lp - 1
                     ll = 2*l
                     ind = k*k/4 + lp + (i-2)*16
                     DO m = ms , msp
                        mm = m
                        mp = m + l
                        jm = INT(2.01*CAT(mm,3)) ! 2 * m quantum number of substate mm
                        IF ( mp.GT.NSTOP(i) ) GOTO 4
                        ilg = (-1)**INT(si-CAT(mp,3)) ! 2 * m quantum number of substate mp
                        jmp = -INT(2.01*CAT(mp,3))
                        fc = WTHREJ(isi,kk,isi,jmp,ll,jm)
                        ite = 1
 2                      IF ( ila.EQ.1 ) x = DBLE(ARM(mp,5))
     &                       *DBLE(ARM(mm,5)) + DIMAG(ARM(mp,5))
     &                       *DIMAG(ARM(mm,5))
                        IF ( ila.NE.1 ) x = DBLE(ARM(mp,5))
     &                       *DIMAG(ARM(mm,5)) - DBLE(ARM(mm,5))
     &                       *DIMAG(ARM(mp,5))
                        Bten(ind) = Bten(ind) + x*fc*ilg
                        IF ( ite.EQ.2 ) GOTO 6
 4                      IF ( iha.NE.1 .OR. Icl.NE.Lmax ) THEN
                           ite = 2
                           mp = mp - 2*l
                           IF ( mp.GE.NSTART(i) ) THEN
                              jmp = INT(2.01*CAT(mp,3)) ! 2 * m quantum number of substate mp
                              jm = -jm
                              fc = WTHREJ(isi,kk,isi,jmp,ll,jm)
                              ilg = (-1)**INT(si+CAT(mp,3))
                              GOTO 2
                           ENDIF
                        ENDIF
 6                      CONTINUE
                     ENDDO ! Loop over m
                     IF ( Icl.EQ.Lmax ) Bten(ind) = Bten(ind)
     &                    *ce/(2.*SPIN(1)+1.)
                  ENDDO ! Loop over lp
               ENDIF ! If isi.GE.k
            ENDDO ! Loop over kp
         ENDIF ! If ms.NE.0
      ENDDO ! Loop over level i
      END
 
C----------------------------------------------------------------------
C SUBROUTINE TENS
C
C Called by: FTBM, GOSIA
C Calls:     DJMM
C
C Purpose:
C
C Uses global variables:
C      IAXS   - axial symmetry flag (readonly)
C      IEXP   - experiment number (readonly)
C      NMAX   - number of levels (readonly)
C      TETACM - theta of particle detector in center of mass frame (readonly)
C      ZETA   - various coefficients (read/write)
C
C Formal parameter:
C      Bten   - 
 
      SUBROUTINE TENS(Bten)
      IMPLICIT NONE
      REAL*8 arg , Bten , DJMM
      INTEGER*4 i , ind , inz , iph , ix , k , k1 , kp , l , lp , 
     &          lpp , lx , lxx
      DIMENSION Bten(*)
      REAL*8 EPS, EROOT, FIEX
      INTEGER*4 IEXP, IAXS
      COMMON /KIN   / EPS(50) , EROOT(50) , FIEX(50,2) , IEXP , IAXS(50)
      REAL*8 TETACM, TREP, DSIGS
      COMMON /TCM   / TETACM(50) , TREP(50) , DSIGS(50)
      REAL*8 ZETA
      INTEGER*4 LZETA
      COMMON /CCOUP / ZETA(155600) , LZETA(8)
      INTEGER*4 NMAX , NDIM , NMAX1
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      ix = NMAX*28
      arg = 1.570796327 + TETACM(IEXP)/2.
      DO i = 1 , ix
         ZETA(i) = 0.
      ENDDO

      DO i = 2 , NMAX ! For each level except ground state
         DO kp = 1 , 7 , 2
            k = kp - 1
            k1 = INT(DBLE(k)/2.+.01)
            IF ( k.EQ.0 ) THEN
               ind = (i-2)*16 + 1
               inz = (i-1)*28 + 1
               ZETA(inz) = Bten(ind)
            ELSE
               DO lp = 1 , kp
                  IF ( IAXS(IEXP).NE.0 .OR. lp.EQ.1 ) THEN
                     inz = (i-1)*28 + k1*7 + lp
                     l = lp - 1
                     DO lpp = 1 , kp
                        ind = k*k/4 + lpp + (i-2)*16
                        lx = lpp - 1
                        lxx = lx
 2                      iph = (-1)**(l+INT(DBLE(lxx)/2.))
                        ZETA(inz) = ZETA(inz) + Bten(ind)
     &                              *iph*DJMM(arg,k,lx,l)
                        IF ( lpp.NE.1 ) THEN
                           IF ( lx.GE.0 ) THEN
                              lx = -lx
                              lxx = lx - 1
                              GOTO 2
                           ENDIF ! if lx .ge. 0
                        ENDIF ! if lpp .ne. 1
                     ENDDO ! Loop over lpp
                  ENDIF ! if iaxs .ne.0 .or. lp.eq.1
               ENDDO ! Loop over lp
            ENDIF ! If k .eq. 0
         ENDDO ! Loop over kp
      ENDDO ! Loop over level i
      END
 
C----------------------------------------------------------------------
C FUNCTION DJMM
C
C Called by: GOSIA, ROTATE, TENS
C
C Purpose: calculate the rotation functions D^k_{\xi \xi^\prime}
C
C Uses global variables:
C      B      - array of factorials
C      BEQ    - identifier for angle for rotations
C
C Formal parameters:
C      Beta   - rotation angle
C      K      - K
C      Kpp    - \xi
C      Kp     - \xi^\prime
C
C Return value:
C      Element of rotation matrix D
C
C Note that to be efficient, this function remembers values that it has
C previously calculated. For this to work, the variable djm has to be
C declared with a SAVE statement, otherwise the variable is an automatic one,
C which is created freshly each time the function is called.

 
      REAL*8 FUNCTION DJMM(Beta,K,Kpp,Kp)
      IMPLICIT NONE
      REAL*8 b1 , b2 , be , Beta , cb , ctb , djm , f , g , 
     &       sb , sk , ul
      INTEGER*4 iczy , ifla , ifza , ill , j , ja , jb , jc , jd , K , 
     &          Kp , Kpp , lca , loc , mas , mis
      DIMENSION djm(525) , iczy(525)
      REAL*8 BEQ
      COMMON /IDENT / BEQ
      REAL*8 B
      COMMON /CB    / B(20)
      SAVE djm , iczy ! Added N. Warr Jul2007
      
      ifza = 1
      IF ( Beta.LT.0. ) ifza = (-1)**(Kp+Kpp)
      sk = DBLE(K)
      ul = sk*((sk-1.)*(4.*sk+7)/6.+1.)
      lca = INT(ul+.1)

C     Calculate position in djm and iczy arrays
      loc = lca + (2*K+1)*Kp + Kpp + K + 1

      IF ( ABS(BEQ-ABS(Beta)).GT.1.E-6 ) THEN ! If beta doesn't match the identifier, initialise
         BEQ = ABS(Beta)
         DO ill = 1 , 525
            iczy(ill) = 0
         ENDDO
      ELSEIF ( iczy(loc).EQ.1 ) THEN ! We have already calculated it, so return
         DJMM = djm(loc)*ifza
         GOTO 99999
      ENDIF

C     We have to calculate it
      be = BEQ/2.
      cb = COS(be)
      sb = SIN(be)
      ifla = 0
      IF ( BEQ.GT..01 .AND. ABS(BEQ-6.2832).GT..01 ) ifla = 1
      IF ( ifla.NE.1 ) THEN
         IF ( Kp.EQ.Kpp ) THEN
            sb = 1.
         ELSE
            DJMM = 0.
            RETURN
         ENDIF
      ENDIF
       
      ctb = cb*cb/sb/sb
      ja = K + Kp + 1 ! K + \xi^\prime (+1 for indexing array)
      jb = K - Kp + 1 ! K - \xi^\prime (+1 for indexing array)
      jc = K + Kpp + 1 ! K + \xi (+1 for indexing array)
      jd = K - Kpp + 1 ! K - \xi (+1 for indexing array)
      b1 = B(ja)*B(jb)*B(jc)*B(jd) ! B array holds factorials

      ja = Kp + Kpp
      jb = 2*K - Kp - Kpp
      IF ( ABS(BEQ-3.141592654).LT..01 .AND. ja.LT.0 ) ifla = 3
      IF ( ifla.EQ.3 ) cb = 1.
      f = (-1)**(K-Kp)*(cb**ja)*(sb**jb)*SQRT(b1)
      mis = 0
      IF ( ja.LT.0 ) mis = -ja
      mas = K - Kpp
      IF ( Kpp.LT.Kp ) mas = K - Kp
      ja = Kp + Kpp + mis + 1
      jb = K - Kpp - mis + 1
      jc = K - Kp - mis + 1
      jd = mis + 1
      b2 = B(ja)*B(jb)*B(jc)*B(jd)
      IF ( ifla.NE.3 ) THEN
         g = (-ctb)**mis/b2
         DJMM = g
         ja = mis + 1
         IF ( mas.GE.ja ) THEN
            DO j = ja , mas
               g = -g*ctb*(K-Kpp-j+1)*(K-Kp-j+1)/(Kp+Kpp+j)/j
               DJMM = DJMM + g
            ENDDO
         ENDIF
         IF ( ifla.EQ.0 ) DJMM = g
         DJMM = DJMM*f*ifza
         djm(loc) = DJMM/ifza
         iczy(loc) = 1
         RETURN
      ENDIF
      DJMM = f*ifza/((-sb*sb)**mis)/b2
      djm(loc) = DJMM/ifza
      iczy(loc) = 1
99999 END
 
C----------------------------------------------------------------------
C SUBROUTINE FTBM
C
C Called by: GOSIA, KONTUR, MINI
C Calls:     ALLOC, APRAM, BRANR, CEGRY, CHMEM, INTG, LOAD, MIXR, SETIN, SNAKE
C            STING, TENB, TENS
C
C Purpose: main routine to perform the calculation with a given set of matrix
C          elements.
C
C Uses global variables:
C      ACCA   - accuracy
C      ACCUR  - accuracy required
C      ARM    - excitation amplitudes of substates.
C      CAT    - substates of levels (n_level, J, m)
C      CHIS11 - chi squared
C      ELM    - matrix elements given by user
C      EMH    - matrix element
C      IAXS   - axial symmetry flag
C      IEXP   - experiment number
C      IGRD   -
C      ILE    - yield number for each detector
C      INM    - index of matrix element
C      INTR   - flag to swap chisqr and log(chisqr)
C      IPRM   - printing flags (see suboption PRT of OP,CONT)
C      ISMAX  - number of substates used
C      ISO    - isotropic flag
C      ITAK2  -
C      IY     - index of experimental yields
C      JSKIP  - Experiments to skip during minimisation.
C      KSEQ   - index into ELM for pair of levels, and into EN or SPIN
C      LFL    -
C      LFL1   -
C      LFL2   -
C      LMAX   - ground-state spin + 1
C      LP3    - maximum number of levels (100)
C      LP6    - maximum number of Ge detectors (32)
C      LP8    - (104)
C      LP9    - last 2100 words of ZETA array (47900)
C      LP10   - maximum number of magnetic substates (1200)
C      LP11   - LP8 - 1 (2800)
C      LP13   - LP9 + 1 (47901)
C      LP14   - maximum space for collision functions (4900)
C      MEMAX  - number of matrix elements
C      MEMX6  - number of matrix elements with E1...6 multipolarity
C      NANG   - number of gamma-ray detectors for each experiment
C      NDIM   - maximum number of levels (100)
C      NEXPT  - number of experiments
C      NLIFT  - number of lifetimes
C      NMAX   - number of levels
C      NSTART - index in CAT of first substate associated with a level
C      NSTOP  - index in CAT of last substate associated with a level
C      NWR    - number of datapoints used in fit
C      NYLDE  - number of yields
C      SPIN   - spin of level
C      ZETA   - the coupling constants
C      ZPOL   - dipole term (GDR excitation)
C
C Formal parameters:
C      Icll   -
C      Chisq  - chi square
C      Idr    - number of decays
C      Ncall  -
C      Chilo  - chi square of logs
C      Bten   -

      SUBROUTINE FTBM(Icll,Chisq,Idr,Ncall,Chilo,Bten)
      IMPLICIT NONE
      REAL*8 aval , Bten , Chilo , chis1 , chish , Chisq , chisx , 
     &       chx , fc , fx , polm , pr , prop , val , wz
      INTEGER*4 i1 , i11 , iapx , Icll , idec , Idr , iflg , ii , ile1 ,
     &          ile2 , ile3 , ilin , indx , inko
      INTEGER*4 inp , inpo , inpx , inzz , inzzz , issp , itemp , ixx , 
     &          izzz
      INTEGER*4 j , jj , jjgg , jjj , jk , jkl , jm , jmf , jmt , jmte ,
     &          jpp , jpz , jy , k , karm , kk , kk6 , kkx , kmt
      INTEGER*4 knm , kx , larm , lcc , lcou , licz , lix , llx , lm , 
     &          lmh , loc , loch , loct
      INTEGER*4 lp , lpit , lput , lpx , lpxd , ls , lst
      INTEGER*4 luu , lx , Ncall , nlin , nowr , npoz , nrest , nwyr
      DIMENSION jmte(6) , prop(6) , Bten(*)
      REAL*8 XA , XA1 , EP , TLBDG , VINF
      INTEGER*4 NEXPT , IZ , IZ1
      COMMON /CX    / NEXPT , IZ , XA , IZ1(50) , XA1(50) , EP(50) , 
     &                TLBDG(50) , VINF(50)
      INTEGER*4 NSTART , NSTOP
      COMMON /CEXC0 / NSTART(101) , NSTOP(100)
      REAL*8 EG, CC, AGELI, Q
      INTEGER*4 NICC, NANG, ISPL
      COMMON /CCC   / EG(50) , CC(50,5) , AGELI(50,200,2) , Q(3,200,8) ,
     &                NICC , NANG(200) , ISPL
      INTEGER*4 NWR
      COMMON /ILEWY / NWR
      REAL*8 CHIS11
      COMMON /CH1T  / CHIS11
      INTEGER*4 IGRD
      COMMON /IGRAD / IGRD
      REAL*8 EMH
      INTEGER*4 INM , LFL1 , LFL2 , LFL
      COMMON /LCZP  / EMH , INM , LFL1 , LFL2 , LFL
      INTEGER*4 ITAK2
      COMMON /UWAGA / ITAK2
      REAL*8 TAU
      INTEGER*4 KSEQ
      COMMON /LEV   / TAU(100) , KSEQ(1500,4)
      REAL*8 ZETA
      INTEGER*4 LZETA
      COMMON /CCOUP / ZETA(155600) , LZETA(8)
      REAL*8 EPS, EROOT, FIEX
      INTEGER*4 IEXP, IAXS
      COMMON /KIN   / EPS(50) , EROOT(50) , FIEX(50,2) , IEXP , IAXS(50)
      REAL*8 YEXP, CORF , DYEX , UPL , YNRM
      INTEGER*4 IY , NYLDE , IDRN , ILE
      COMMON /YEXPT / YEXP(32,1500) , IY(1500,32) , CORF(1500,32) , 
     &                DYEX(32,1500) , NYLDE(50,32) , UPL(32,50) , 
     &                YNRM(32,50) , IDRN , ILE(32)
      REAL*8 ELM , ELMU , ELML , SA
      COMMON /COMME / ELM(1500) , ELMU(1500) , ELML(1500) , SA(1500)
      INTEGER*4 LMAX
      COMMON /CLM   / LMAX
      REAL*8 EN , SPIN , ACCUR , DIPOL , ZPOL , ACCA
      INTEGER*4 ISO
      COMMON /COEX  / EN(100) , SPIN(100) , ACCUR , DIPOL , ZPOL , 
     &                ACCA , ISO
      REAL*8 CAT
      INTEGER*4 ISMAX
      COMMON /CLCOM8/ CAT(1200,3) , ISMAX
      COMPLEX*16 ARM
      COMMON /AZ    / ARM(1200,7)
      INTEGER*4 LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &          LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      INTEGER*4 NMAX , NDIM , NMAX1
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      INTEGER*4 IPATH , MAGA
      COMMON /PTH   / IPATH(100) , MAGA(100)
      INTEGER*4 IPRM
      COMMON /PRT   / IPRM(20)
      INTEGER*4 MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR
      COMMON /CEXC  / MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR(1500)
      INTEGER*4 JSKIP
      COMMON /SKP   / JSKIP(50)
      INTEGER*4 NLIFT
      COMMON /LIFE  / NLIFT
      INTEGER*4 LNY , INTR , IPS1
      COMMON /LOGY  / LNY , INTR , IPS1
      DATA pr/0./,lmh/0/,loc/0/,loch/0/

      issp = 0
      Chilo = 0.
      fx = 2.*SPIN(1) + 1.
      Chisq = 0.
      LFL = 0
      chis1 = 0.
      ixx = NDIM*MEMAX + LP11 ! LP11 is 2800
      IF ( ixx.GT.LP7 ) THEN
         STOP 'Too many matrix elements for the ZETA array'
      ENDIF

      DO i1 = 1 , ixx
         ZETA(i1) = 0.
      ENDDO

      DO ii = 1 , LP6 ! LP6 is 32
         ILE(ii) = 1
      ENDDO

      itemp = 0
      NWR = 0
      iapx = 1

      DO jkl = 1 , NEXPT ! For each experiment
         IEXP = jkl
         IGRD = 0
         LFL2 = 1
         IF ( ITAK2.EQ.-1 ) THEN
            DO larm = 1 , 4
               DO karm = 1 , LP10 ! LP10 is 1200
                  ARM(karm,larm) = (0.,0.)
               ENDDO
            ENDDO
         ENDIF
         iflg = 0
         IF ( IEXP.NE.1 ) THEN
            kk = NANG(IEXP) ! Number of detector angles
            DO jjj = 1 , LP6 ! LP6 is 32
               ILE(jjj) = ILE(jjj) + NYLDE(IEXP-1,jjj)
            ENDDO
         ENDIF
         lp = 3
         IF ( JSKIP(jkl).EQ.0 ) GOTO 200
         IF ( MAGA(IEXP).EQ.0 ) lp = 1
         IF ( Ncall.EQ.0 ) GOTO 150
         IF ( Icll.EQ.4 ) GOTO 100
 50      loch = LP3*(MEMAX-1) + NMAX + LP11 ! LP3 is 100, LP11 is 2800
         DO k = 1 , loch
            ZETA(k) = 0.
         ENDDO
         CALL LOAD(IEXP,1,2,0.D0,jj)
         DO k = 1 , LMAX ! For each multipolarity up to ground-state spin + 1
            fc = 2.
            IF ( k.EQ.LMAX ) fc = 1.
            IF ( DBLE(INT(SPIN(1))).LT.SPIN(1) ) fc = 2.
            loc = 0
            polm = DBLE(k-1) - SPIN(1) ! Multipolarity - ground-state spin
            CALL LOAD(IEXP,3,2,polm,jj) ! Calculate parameters
            CALL PATH(jj) ! Find path
            CALL LOAD(IEXP,2,2,polm,jj) ! Calculate parameters
            CALL APRAM(IEXP,1,1,jj,ACCA) ! Calculate excitation amplitudes
            IF ( Ncall.NE.0 ) THEN
               IF ( Icll.NE.3 ) THEN
                  DO indx = 1 , MEMX6 ! Loop over E1...6 matrix elements
                     CALL APRAM(IEXP,0,indx,jj,ACCA) ! Calculate excitation amplitudes
                     kx = 0
                     DO i11 = 1 , NMAX ! Loop over levels
                        IF ( NSTART(i11).NE.0 ) THEN
                           loc = LP3*(indx-1) + i11 + LP11
                           jpp = INT(2.*SPIN(i11)+1.)
                           lpx = MIN(lp,jpp)
                           IF ( ISO.NE.0 ) lpx = NSTOP(i11)
     &                          - NSTART(i11) + 1
                           DO lpxd = 1 , lpx ! Loop over substates for level
                              kx = kx + 1
                              ZETA(loc) = ZETA(loc) + fc*DBLE(ARM(kx,5))
     &                           *DBLE(ARM(kx,6))
     &                           /fx + fc*DIMAG(ARM(kx,5))
     &                           *DIMAG(ARM(kx,6))/fx
                           ENDDO ! Loop on lpxd
                        ENDIF ! IF ( NSTART(i11).NE.0 )
                     ENDDO ! Loop over levels
                  ENDDO ! Loop on E1...6 matrix elements
               ENDIF ! IF ( Icll.NE.3 )
            ENDIF ! IF ( Ncall.NE.0 )
            CALL TENB(k,Bten,LMAX)
         ENDDO ! Loop on multipolarity k

         IF ( loc.NE.0 ) THEN
            REWIND 14
            WRITE (14,*) (ZETA(i11),i11=LP8,loch)
         ENDIF
         CALL TENS(Bten)
         IF ( Ncall.EQ.0 ) GOTO 200
         IF ( Icll.GE.2 ) GOTO 200
         llx = 28*NMAX
         DO lx = 1 , llx
            ZETA(LP9+lx) = ZETA(lx) ! LP9 is 47900
         ENDDO
         IF ( Icll.NE.1 ) GOTO 200
 100     iapx = 0
         issp = 1
         CALL LOAD(IEXP,1,1,0.D0,jj) ! Calculate parameters
         CALL ALLOC(ACCUR)           ! Calculate ranges
         CALL SNAKE(IEXP,ZPOL)       ! Calculate collision functions
         CALL SETIN                  ! Calculate adiabatic parameters
         DO k = 1 , LMAX
            polm = DBLE(k-1) - SPIN(1)
            CALL LOAD(IEXP,2,1,polm,kk)
            IF ( IPRM(7).EQ.-1 ) WRITE (22,99001) polm , IEXP
99001       FORMAT (1X//40X,'EXCITATION AMPLITUDES'//10X,'M=',1F4.1,5X,
     &              'EXPERIMENT',1X,1I2//5X,'LEVEL',2X,'SPIN',2X,'M',5X,
     &              'REAL AMPLITUDE',2X,'IMAGINARY AMPLITUDE'//)
            CALL STING(kk) ! Calculate excitation amplitudes
            CALL PATH(kk)
            CALL INTG(IEXP) ! Integrate
            CALL TENB(k,Bten,LMAX)
            IF ( IPRM(7).EQ.-1 ) THEN
               DO j = 1 , ISMAX
                  WRITE (22,99002) INT(CAT(j,1)) , CAT(j,2) , CAT(j,3) ,
     &                             DBLE(ARM(j,5)) , DIMAG(ARM(j,5))
99002             FORMAT (7X,1I2,3X,1F4.1,2X,1F4.1,2X,1E14.6,2X,1E14.6)
               ENDDO
            ENDIF ! IF ( IPRM(7).EQ.-1 )
         ENDDO ! Loop on k
         CALL TENS(Bten)
         IF ( IPRM(7).EQ.-1 ) THEN
            DO jjgg = 2 , NMAX
               loct = (jjgg-1)*28 + 1
               WRITE (22,99003) jjgg , ZETA(loct)
99003          FORMAT (2X,'LEVEL',1X,1I2,10X,'POPULATION',1X,1E14.6)
            ENDDO
         ENDIF
         GOTO 200
 150     IF ( iflg.EQ.1 ) THEN
            itemp = 1
            iflg = 2
            GOTO 50
         ELSE
            IF ( iflg.EQ.2 ) GOTO 300
            itemp = 2
            iflg = 1
            GOTO 100
         ENDIF
 200     CALL CEGRY(Chisq,itemp,Chilo,Idr,nwyr,Icll,issp,0)
         issp = 0
         IF ( Ncall.EQ.0 .AND. JSKIP(jkl).NE.0 ) THEN
            IF ( Ncall.EQ.0 ) GOTO 150
            GOTO 200
         ELSE
            NWR = NWR + nwyr
            IF ( Icll.LE.2 .AND. JSKIP(jkl).NE.0 ) THEN
               IF ( IEXP.EQ.1 ) chish = CHIS11
               IF ( Icll.EQ.1 ) chis1 = CHIS11
               IF ( Icll.EQ.0 ) chis1 = Chisq
               LFL2 = 0
               IGRD = 1
               IF ( ITAK2.EQ.-1 ) LFL = 1
               REWIND 14
               READ (14,*) (ZETA(i11),i11=LP8,loch)
               DO larm = 1 , 4
                  DO karm = 1 , LP10
                     ARM(karm,larm) = (0.,0.)
                  ENDDO
               ENDDO
               chisx = 0.
               llx = 28*NMAX
               DO lix = 1 , llx
                  ZETA(LP9+lix) = ZETA(lix) ! LP9 is 47900
               ENDDO
               CALL CEGRY(chisx,itemp,Chilo,Idr,nwyr,0,0,1)
               DO knm = 1 , MEMAX ! Loop over matrix elements
                  INM = knm
                  chisx = 0.
                  EMH = ELM(INM)
                  ELM(INM) = 1.05*EMH
                  lcc = LP3*(INM-1) + LP11
                  DO lst = 2 , NMAX ! For all states except ground state
                     wz = ZETA(lst+lcc)
                     inpx = (lst-1)*28
                     DO jy = 1 , 4
                        inp = inpx + (jy-1)*7
                        IF ( jy.EQ.1 ) pr = ZETA(LP13+inp) + 1.E-12
                        jmf = 2*jy - 1
                        IF ( IAXS(IEXP).EQ.0 ) jmf = 1
                        DO jm = 1 , jmf
                           inp = inp + 1
                           ZETA(inp) = ZETA(inp+LP9)*(1.+.1*EMH*wz/pr)
                        ENDDO
                     ENDDO
                  ENDDO
                  CALL CEGRY(chisx,itemp,Chilo,Idr,nwyr,0,0,0)
                  ELM(INM) = EMH
               ENDDO
               IF ( ITAK2.EQ.-1 .AND. LFL1.NE.0 ) THEN
                  IF ( IPRM(17).NE.0 ) THEN
                     kmt = ABS(IPRM(17))
                     WRITE (22,99004) IEXP
99004                FORMAT (1X///20X,'EXPERIMENT',11X,1I2,5X,
     &                       'D(LOG(P))/D(LOG(ME)) MAP'/20X,52('-')///)
                     nlin = (NMAX-2)/6 + 1
                     nrest = NMAX - 1 - 6*(nlin-1)
                     DO ilin = 1 , nlin
                        npoz = 6
                        IF ( ilin.EQ.nlin ) npoz = nrest
                        inpo = (ilin-1)*6 + 2
                        inko = inpo + npoz - 1
                        lpit = 0
                        DO lm = inpo , inko
                           lpit = lpit + 1
                           jmte(lpit) = lm
                        ENDDO
                        WRITE (22,99005) (jmte(lm),lm=1,lpit)
99005                   FORMAT (5X,'LEVEL',6(8X,1I2,9X))
                        WRITE (22,99006)
     &                         (ZETA(LP13+(jpz-1)*28),jpz=inpo,inko)
99006                   FORMAT (1X,'EXC.PROB.',6(5X,1E10.4,4X))
                        DO jmt = 1 , kmt
                           lput = 0
                           DO ls = inpo , inko
                              lput = lput + 1
                              prop(lput) = 0.
                              DO lm = 1 , MEMX6
                                 inzz = ls + LP3*(lm-1) + LP11
                                 inzzz = LP13 + (ls-1)*28
                                 IF ( ABS(ZETA(inzzz)).LT.1.E-20 )
     &                                ZETA(inzzz) = 1.E-20
                                 val = 2.*ELM(lm)*ZETA(inzz)/ZETA(inzzz)
                                 aval = ABS(val)
                                 IF ( aval.GT.ABS(prop(lput)) ) THEN
                                    prop(lput) = val
                                    lmh = lm
                                    jmte(lput) = lm
                                 ENDIF
                              ENDDO
                              izzz = (lmh-1)*LP3 + LP11 + ls
                              ZETA(izzz) = 0.
                           ENDDO
                           WRITE (22,99007)
     &                            (jmte(lcou),prop(lcou),lcou=1,npoz)
99007                      FORMAT (10X,6(2X,'(',1X,1I3,1X,1E8.2,')',2X))
                        ENDDO
                     ENDDO
                     REWIND 14
                     READ (14,*) (ZETA(i11),i11=LP8,loch)
                     IF ( IPRM(17).LT.0 ) GOTO 300
                  ENDIF
                  LFL = 0
                  WRITE (22,99008) IEXP
99008             FORMAT (10X,'EXPERIMENT',1X,1I2/10X,
     &                    'D(LOG(Y)/D(LOG(ME))',//)
                  ile1 = ILE(1) + NYLDE(IEXP,1) - 1
                  ile3 = ILE(1)
                  licz = 0
                  DO ile2 = ile3 , ile1 ! For each experimental yield
                     licz = licz + 1
                     idec = IY(ile2,1) ! Decay number
                     IF ( idec.GT.1000 ) idec = idec/1000
                     luu = 6*licz - 5
                     jk = (luu-1)/LP10 + 1
                     kk = luu - LP10*(jk-1)
                     kk6 = kk + 5
                     WRITE (22,99009) KSEQ(idec,3) , KSEQ(idec,4) , ! Level numbers
     &                                (INT(DBLE(ARM(kkx,jk))),
     &                                DIMAG(ARM(kkx,jk)),kkx=kk,kk6)
99009                FORMAT (2X,1I2,'--',1I2,5X,
     &                       6('(',1I3,2X,1E8.2,')',3X))
                  ENDDO ! Loop on ile2
               ENDIF ! IF ( ITAK2.EQ.-1 .AND. LFL1.NE.0 )
            ENDIF ! IF ( Icll.LE.2 .AND. JSKIP(jkl).NE.0 )
         ENDIF ! ELSE of IF ( Ncall.EQ.0 .AND. JSKIP(jkl).NE.0 )
 300     CONTINUE
      ENDDO ! Loop on experiments

      IF ( ITAK2.EQ.-1 .AND. Icll.LT.2 ) ITAK2 = 0
      IF ( Ncall.NE.0 ) THEN
         IF ( Icll.LE.2 ) THEN
            IF ( Icll.EQ.1 ) CALL CEGRY(Chisq,itemp,Chilo,Idr,nowr,7,
     &                                  issp,0)
         ENDIF
         CALL BRANR(Chisq,NWR,Chilo) ! Branching ratios
         CALL MIXR(NWR,0,Chisq,Chilo) ! Mixing ratios
         CALL CHMEM(NWR,Chisq,Chilo) ! Compare matrix elements
         NWR = NWR + NLIFT
         Chisq = Chisq/NWR
         IF ( INTR.NE.0 ) THEN
            chx = Chisq ! Swap chisqr and log(chisqr)
            Chisq = Chilo
            Chilo = chx
         ENDIF
      ENDIF
      END
 
C----------------------------------------------------------------------
C SUBROUTINE MINI
C
C Called by: GOSIA
C Calls:     FTBM, LIMITS
C
C Purpose: perform the minimization
C
C Uses global variables:
C      CHIS11 - chi squared
C      CORF   - internal correction factors
C      DEVD   -
C      DEVU   -
C      DLOCK  - limit of derivative below which matrix element fixed if LOCKS=1
C      ELM    - matrix elements
C      ELMH   - matrix element
C      GRAD   - partial derivative of chi squared wrt. matrix element
C      HLMLM  - old value of matrix element or chi squared
C      ICS    - read internal correction factors from file rather than recalculating
C      IFBFL  - calculate derivatives with forward-backward method
C      INTR   - flag to swap chisqr and log(chisqr)
C      IPRM   - printing flags (see suboption PRT of OP,CONT)
C      IPS1   - terminate after calculating and writing correction factors
C      ITAK2  -
C      IUNIT3 - unit for TAPE3
C      IVAR   - fixed, correlated or free flag
C      JENTR  - flag set to 0 normally, 1 in OP,ERRO
C      KFERR  - error flag for minimization
C      KVAR   -
C      LFL1   -
C      LNY    - use logs to calculate chi squared
C      LOCKF  - fix those with most significat derivative
C      LOCKS  - lock flag. if LOCKS=1, fix at first stage of minimisation
C      LP4    - 1500
C      LP6    - 32
C      MEMAX  - number of matrix elements
C      NLOCK  - number of matrix elements to lock
C      NWR    - number of datapoints used in fit
C      SA     - ratio of elements for correlated elements
C
C Formal parameters:
C      Chisq  - chi squared of minimization
C      Chiok  - desired chi squared
C      Nptl   - number of iterations allowed
C      Conv   - parameter for convergence test
C      Imode  - mode of minimization
C      Idr    -
C      Xtest  -
C      Ips    -
C      Is     - generate input for sigma program flag
C      Jjh    -
C      Bten   -
C
C FTBM does the main calculation and LIMITS makes sure the matrix elements
C don't go outside the limits specified by the user.

      SUBROUTINE MINI(Chisq,Chiok,Nptl,Conv,Imode,Idr,Xtest,Ips,Is,Jjh,
     &                Bten)
      IMPLICIT NONE
      REAL*8 a , a0 , a1 , b , Bten , c , ccd , chd , chil , chilo , 
     &       Chiok , chirf , chis12 , chis13 , chisf , chisp , Chisq , 
     &       chiss , chl
      REAL*8 chx , cmax , Conv , crit , dl , dm , f1 , f2 , flt
      REAL*8 gradp , ht , p , q , rfk , sel , shl , sumg1 , 
     &       sumg2 , sumht , uxa , xkat , Xtest
      INTEGER*4 i , icl1 , icl2 , icount , Idr , iht , iin , Imode , 
     &          indx1 , inmx , ino , ipas , ipm
      INTEGER*4 Ips , Is , istec , itf , j , jcoup , jcp , jin , 
     &          Jjh , jjj , jlin , jnm , jpr , jsa , jst
      INTEGER*4 kh2 , kkk , l , lnm , metf , mvfl , ncall , nlinn , 
     &          noflg , Nptl
      DIMENSION ipm(10) , Bten(*) , gradp(1500)
      REAL*8 GRAD , HLMLM , ELMH
      COMMON /DUMM  / GRAD(1500) , HLMLM(1500) , ELMH(1500)
      INTEGER*4 NWR
      COMMON /ILEWY / NWR
      REAL*8 CHIS11
      COMMON /CH1T  / CHIS11
      INTEGER*4 LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &          LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      INTEGER*4 ITAK2
      COMMON /UWAGA / ITAK2
      REAL*8 YEXP, CORF , DYEX , UPL , YNRM
      INTEGER*4 IY , NYLDE , IDRN , ILE
      COMMON /YEXPT / YEXP(32,1500) , IY(1500,32) , CORF(1500,32) , 
     &                DYEX(32,1500) , NYLDE(50,32) , UPL(32,50) , 
     &                YNRM(32,50) , IDRN , ILE(32)
      REAL*8 DEVD, DEVU
      COMMON /DFTB  / DEVD(1500) , DEVU(1500)
      INTEGER*4 IPRM
      COMMON /PRT   / IPRM(20)
      REAL*8 EMH
      INTEGER*4 INM , LFL1 , LFL2 , LFL
      COMMON /LCZP  / EMH , INM , LFL1 , LFL2 , LFL
      INTEGER*4 MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR
      COMMON /CEXC  / MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR(1500)
      REAL*8 ELM , ELMU , ELML , SA
      COMMON /COMME / ELM(1500) , ELMU(1500) , ELML(1500) , SA(1500)
      INTEGER*4 KVAR
      COMMON /SEL   / KVAR(1500)
      REAL*8 DLOCK
      INTEGER*4 LOCKF , NLOCK , IFBFL , LOCKS
      COMMON /FIT   / LOCKF , NLOCK , IFBFL , LOCKS , DLOCK
      INTEGER*4 KFERR
      COMMON /ERRAN / KFERR
      INTEGER*4 LNY , INTR , IPS1
      COMMON /LOGY  / LNY , INTR , IPS1
      INTEGER*4 JENTR , ICS
      COMMON /ERCAL / JENTR , ICS
      INTEGER*4 IUNIT3 , JZB
      COMMON /SWITCH/ JZB , IUNIT3
      DATA chirf/0./,dm/0./,sumg2/0./

C     Initialise gradp to zero for each matrix element
      DO i = 1 , MEMAX
         gradp(i) = 0.
      ENDDO

C     Initialise some parameters to zero
      icount = 0
      lnm = 0
      LNY = 0
      INTR = 0
      metf = 0
      LFL1 = 0
      ncall = 0
      ITAK2 = 0

C     Handle the different modes
C     Imode = IJKL, where
C        I=1 => fast approximation to calculate chi squared and its partial derivatives
C        I=2 => full Coulomb excitation formalism, but derivatives with fast approximation
C
C        J=0 => steepest descent minimization
C        J=1 => gradient minimization
C
C        K=0 => absolute changes of matrix elements
C        K=1 => relative changes
C
C        L=0 => yields, branching ratios used to calculate chi squared
C        L=1 => logs used to claculate chi squared
      IF ( Imode.LT.2000 ) THEN ! Fast approximation for chi squared and derivatives
         icl1 = 0
         icl2 = 3
         IF ( Imode.GE.1100 ) metf = 1
         IF ( (Imode-1000-100*metf).GE.10 ) lnm = 1
         IF ( (Imode-1000-100*metf-10*lnm).EQ.1 ) LNY = 1 ! Use logs
         IF ( JENTR.EQ.1 ) GOTO 200 ! If we are in OP,ERRO, jump
         IF ( ICS.NE.0 ) THEN ! Read correction factors from file, rather than recalculating
            REWIND 11
            DO jnm = 1 , LP4 ! LP4 is 1500
               READ (11) (CORF(jnm,kh2),kh2=1,LP6) ! LP6 is 32
            ENDDO
            ICS = 0
            GOTO 200
         ENDIF
      ELSE ! Full Coulomb excitation formalism for chi squared, fast approx for derivatives
         icl1 = 1
         IF ( Imode.GE.2100 ) metf = 1
         IF ( (Imode-2000-100*metf).GE.10 ) lnm = 1
         IF ( (Imode-2000-100*metf-10*lnm).EQ.1 ) LNY = 1 ! Use logs
         icl2 = 4
         IF ( Ips.NE.0 ) THEN
            IF ( Ips.EQ.1 ) THEN
               IF ( IPRM(4).EQ.-1 ) ITAK2 = -2
            ELSE
               IF ( IPRM(4).LT.0 ) ITAK2 = -2
            ENDIF
            icl1 = 4
            IF ( ITAK2.EQ.-2 ) icl1 = 1
            IF ( icl1.EQ.4 ) GOTO 200
         ENDIF
      ENDIF

C     Call FTBM to perform a single calculation
 100  CALL FTBM(0,chiss,Idr,0,chl,Bten)

C     Write correction factors
      REWIND 11
      DO jnm = 1 , LP4
         WRITE (11) (CORF(jnm,kh2),kh2=1,LP6)
      ENDDO
       
      IF ( IPS1.EQ.0 ) RETURN ! If IPS1 = 0, terminate after writing correction factors
       
 200  noflg = 0
      ncall = 1
 300  sumht = 0.
      IF ( LNY.EQ.1 ) INTR = 1
      LFL1 = 1
      ITAK2 = ITAK2 + 1

      icount = icount + 1
      IF ( icount.GT.Nptl ) THEN
         IF ( KFERR.EQ.1 ) RETURN
         IF ( Ips.EQ.0 ) WRITE (22,99001) Nptl
99001    FORMAT (5X,'MINIMIZATION STOPPED-NUMBER OF STEPS NPTL=',1I5,1X,
     &           'EXCEEDED')
         IF ( Ips.EQ.0 ) WRITE (22,99010) chil
         INTR = 0
         RETURN
      ELSE
         IF ( ITAK2.EQ.IPRM(4) ) ITAK2 = -1
         IF ( ITAK2.EQ.-1 ) THEN
            IF ( KFERR.NE.1 ) THEN
               CALL FTBM(3,chd,Idr,1,chl,Bten)
               CHIS11 = chd*NWR
               CALL FTBM(icl1,Chisq,Idr,ncall,chilo,Bten)
            ENDIF
         ENDIF
         IF ( Ips.EQ.1 ) RETURN
         IF ( icl1.EQ.1 ) CALL FTBM(4,Chisq,Idr,ncall,chilo,Bten)
         IF ( IPRM(8).EQ.-1 .OR. IPRM(13).EQ.-1 ) THEN
            IF ( IPRM(8).EQ.-1 ) IPRM(8) = -2
            IF ( IPRM(13).EQ.-1 ) IPRM(13) = -2
            CALL FTBM(4,ccd,Idr,ncall,chl,Bten)
            IF ( Ips.EQ.2 ) RETURN
         ENDIF
         CALL FTBM(3,chis12,Idr,ncall,chilo,Bten)
         IF ( icl1.EQ.0 ) Chisq = chis12
         uxa = Chisq
         IF ( INTR.EQ.1 ) uxa = chilo
         ipas = 0
         IF ( uxa.LT.Chiok ) Chisq = uxa
         IF ( uxa.LT.Chiok ) GOTO 600
 350     ino = 1
         IF ( metf.EQ.1 ) ipas = ipas + 1
         IF ( IFBFL.EQ.1 ) ino = 2 ! IFBFL = 1 means use forward-backward method
         DO jjj = 1 , ino
            DO jnm = 1 , MEMAX
               GRAD(jnm) = 0.
               IF ( IVAR(jnm).EQ.1 .OR. IVAR(jnm).EQ.2 ) THEN
                  DO jcoup = 1 , MEMAX
                     ELMH(jcoup) = ELM(jcoup)
                  ENDDO
                  DO jcoup = 1 , MEMAX
                     IF ( jnm.NE.jcoup ) THEN
                        IF ( IVAR(jcoup).LT.1000 ) GOTO 355
                        jcp = IVAR(jcoup) - 1000
                        IF ( jcp.NE.jnm ) GOTO 355
                        IF ( IVAR(jnm).EQ.0 ) GOTO 355
                     ENDIF
                     flt = 1.01
                     IF ( jjj.EQ.2 ) flt = .99
                     ELM(jcoup) = ELMH(jcoup)*flt
 355                 CONTINUE
                  ENDDO
                  CALL FTBM(3,chis13,Idr,ncall,chx,Bten)
                  IF ( jjj.EQ.1 ) HLMLM(jnm) = chis13
                  IF ( IFBFL.NE.1 .OR. jjj.NE.1 ) THEN
                     IF ( jjj.EQ.2 ) chis12 = chis13
                     GRAD(jnm) = 100.*(HLMLM(jnm)-chis12)/ELMH(jnm)
                     IF ( IFBFL.EQ.1 ) GRAD(jnm) = GRAD(jnm)/2. ! Forward-backward
                     IF ( lnm.EQ.1 ) GRAD(jnm) = GRAD(jnm)
     &                    *ABS(ELMH(jnm))
                  ENDIF
                  DO jcoup = 1 , MEMAX
                     ELM(jcoup) = ELMH(jcoup)
                  ENDDO
               ENDIF ! If IVAR is 1 or 2
            ENDDO
         ENDDO
         IF ( KFERR.EQ.1 ) THEN
            GRAD(Jjh) = 0.
            IF ( Is.EQ.1 .AND. icount.EQ.1 ) WRITE (IUNIT3,*) ! For sigma program
     &           (NWR*GRAD(jnm),jnm=1,MEMAX)
         ENDIF
         IF ( metf.EQ.1 .AND. ipas.EQ.2 ) THEN
            DO jnm = 1 , MEMAX
               ELM(jnm) = DEVU(jnm)
            ENDDO
            shl = dm/20./sumg2
            sumg1 = 0.
            DO jnm = 1 , MEMAX
               GRAD(jnm) = (DEVD(jnm)*sumg2-GRAD(jnm))/shl
               sumg1 = sumg1 + GRAD(jnm)*GRAD(jnm)
            ENDDO
            sumg1 = SQRT(sumg1)
            p = 0.
            DO jnm = 1 , MEMAX
               GRAD(jnm) = GRAD(jnm)/sumg1
               DEVU(jnm) = ELM(jnm)
               sel = dm*GRAD(jnm)/100.
               IF ( lnm.EQ.1 ) sel = sel*ABS(DEVU(jnm))
               p = p + DEVD(jnm)*GRAD(jnm)
               ELM(jnm) = ELM(jnm) + sel
            ENDDO
            CALL FTBM(3,chis13,Idr,ncall,chx,Bten)
            shl = dm/100.
            DO jnm = 1 , MEMAX
               sel = dm*GRAD(jnm)/50.
               IF ( lnm.EQ.1 ) sel = sel*ABS(DEVU(jnm))
               ELM(jnm) = ELM(jnm) - sel
            ENDDO
            CALL FTBM(3,chis12,Idr,ncall,chx,Bten)
            q = (chis12+chis13-2.*Chisq)/shl/shl
            a0 = q*sumg2/sumg1 - p
            a1 = p*p - 1.
            sumg1 = SQRT(a0*a0+a1*a1+2.*a0*a1*p)
            DO jnm = 1 , MEMAX
               ELM(jnm) = DEVU(jnm)
               GRAD(jnm) = (GRAD(jnm)*a1+DEVD(jnm)*a0)/sumg1
            ENDDO
         ELSE
            sumg2 = 0.
            DO jnm = 1 , MEMAX
               IF ( IVAR(jnm).EQ.1 .OR. IVAR(jnm).EQ.2 ) sumg2 = sumg2 +
     &              GRAD(jnm)*GRAD(jnm)
            ENDDO
            IF ( sumg2.LT.1.E-10 ) GOTO 700
            sumg2 = SQRT(sumg2)
            DO jnm = 1 , MEMAX
               GRAD(jnm) = GRAD(jnm)/sumg2
            ENDDO
            IF ( metf.NE.0 ) THEN
               dm = 0.
               DO jnm = 1 , MEMAX
                  IF ( IVAR(jnm).EQ.2 .OR. IVAR(jnm).EQ.1 ) dm = dm + 
     &                 ELM(jnm)*ELM(jnm)*GRAD(jnm)*GRAD(jnm)
               ENDDO
               dm = SQRT(dm)
               DO jnm = 1 , MEMAX
                  DEVD(jnm) = GRAD(jnm)
                  DEVU(jnm) = ELM(jnm)
                  sel = dm*GRAD(jnm)/20.
                  IF ( lnm.EQ.1 ) sel = sel*ABS(ELM(jnm))
                  ELM(jnm) = ELM(jnm) - sel
               ENDDO
               IF ( IFBFL.EQ.0 ) CALL FTBM(3,chis12,Idr,ncall,chx,Bten)
               GOTO 350
            ENDIF
         ENDIF
         LFL1 = 0
         IF ( lnm.NE.0 ) THEN
            DO jnm = 1 , MEMAX
               GRAD(jnm) = GRAD(jnm)*ABS(ELM(jnm))
            ENDDO
         ENDIF
         sumg1 = 0.
         DO jnm = 1 , MEMAX
            sumg1 = sumg1 + GRAD(jnm)*GRAD(jnm)
         ENDDO
         sumg1 = SQRT(sumg1)
         DO jnm = 1 , MEMAX
            GRAD(jnm) = GRAD(jnm)/sumg1
         ENDDO
         IF ( LNY.EQ.1 ) Chisq = chilo
         IF ( noflg.EQ.0 ) chirf = Chisq
         noflg = 1
         chil = Chisq
         IF ( KFERR.NE.1 ) THEN
            IF ( MOD(icount,IPRM(5)).EQ.0 .OR. icount.EQ.1 )
     &           WRITE (22,99010) Chisq
            WRITE (*,99010) Chisq
            IF ( MOD(icount,IPRM(6)).EQ.0 ) THEN
               WRITE (22,99002)
99002          FORMAT (20X,'GRADIENT'//)
               nlinn = MEMAX/10 + 1
               DO jlin = 1 , nlinn
                  jsa = (jlin-1)*10 + 1
                  DO jin = 1 , 10
                     ipm(jin) = jsa + jin - 1
                  ENDDO
                  jst = MIN(jsa+9,MEMAX)
                  jpr = MIN(10,MEMAX-jsa+1)
                  WRITE (22,99003) (ipm(jin),jin=1,jpr)
99003             FORMAT (5X,10(5X,1I3,4X))
                  WRITE (22,99004) (GRAD(jin),jin=jsa,jst)
99004             FORMAT (5X,10(1X,1E10.4,1X)/)
               ENDDO
            ENDIF
         ENDIF
         IF ( chil.LT.Chiok ) GOTO 600 ! We've achieved desired chi square
         DO l = 1 , MEMAX
            HLMLM(l) = ELM(l)
         ENDDO
         DO l = 1 , MEMAX
            IF ( ABS(GRAD(l)).LE.DLOCK .AND. LOCKS.EQ.1 .AND. 
     &           icount.EQ.1 .AND. IVAR(l).LE.999 .AND. IVAR(l).NE.0 )
     &           THEN
               IF ( KFERR.NE.1 ) KVAR(l) = 0
               IF ( KFERR.NE.1 ) WRITE (22,99005) l , GRAD(l)
99005          FORMAT (1X,'MATRIX ELEMENT',1X,1I3,1X,'LOCKED',3X,
     &                 'DERIVATIVE=',1E14.6)
               IVAR(l) = 0
            ENDIF
         ENDDO
         istec = 0
      ENDIF
       
 400  DO j = 1 , MEMAX
         ELMH(j) = ELM(j)
      ENDDO

C     Find steepest gradient
      istec = istec + 1
      cmax = 0.
      INTR = 0
      inmx = 1
      DO iht = 1 , MEMAX
         IF ( ABS(GRAD(iht)).GT.cmax ) THEN
            cmax = ABS(GRAD(iht))
            inmx = iht
         ENDIF
      ENDDO
       
      ht = .01*ABS(ELM(inmx))/cmax
      mvfl = 0
      IF ( icount.NE.1 .AND. istec.EQ.1 ) THEN
         xkat = 0.
         DO j = 1 , MEMAX
            xkat = xkat + GRAD(j)*gradp(j)
         ENDDO
         DO j = 1 , MEMAX
            gradp(j) = GRAD(j)
         ENDDO
         IF ( xkat.GE..8 ) THEN
            a = 0.
            DO j = 1 , MEMAX
               IF ( IVAR(j).NE.0 .AND. IVAR(j).LE.999 ) THEN
                  a = MAX(a,ABS(GRAD(j)))
                  IF ( ABS(a-ABS(GRAD(j))).LT.1.E-9 ) iin = j
               ENDIF
            ENDDO
            WRITE (22,99011) iin
            IVAR(iin) = 0
            GRAD(iin) = 0.
            gradp(iin) = 0.
         ENDIF
      ENDIF
       
 500  DO j = 1 , MEMAX
         ELM(j) = ELMH(j) - ht*GRAD(j)
      ENDDO

      DO j = 1 , MEMAX
         IF ( IVAR(j).GE.1000 ) THEN  ! For correlated elements
            indx1 = IVAR(j) - 1000    ! Index of element to which it is correlated
            ELM(j) = ELM(indx1)*SA(j) ! SA is the ratio we require
         ENDIF
      ENDDO

      IF ( mvfl.EQ.0 ) THEN
         CALL FTBM(icl2,chisp,Idr,ncall,chilo,Bten)
         DO j = 1 , MEMAX
            ELM(j) = 2.*ELMH(j) - ELM(j)
         ENDDO
         CALL FTBM(icl2,chisf,Idr,ncall,chilo,Bten)
         c = (chisp+chisf-2.*chil)/ht/ht
         b = (chisp-chisf)/ht/2.
         dl = b*b - 2.*c*chil
         IF ( dl.GT.0. ) THEN
            f1 = chil
            f2 = b
         ELSE
            f1 = b
            f2 = c
         ENDIF
         mvfl = 1
         IF ( ABS(f2).LT.1.E-10 ) THEN
            ht = 1.
         ELSE
            ht = -f1/f2
         ENDIF
         GOTO 500
      ELSE
         CALL LIMITS
         CALL FTBM(icl2,Chisq,Idr,ncall,chilo,Bten)
         IF ( Chisq.GE.chil ) THEN
            ht = ht/2.
            IF ( ABS(ht).GE.Conv ) GOTO 500
         ELSE
            chil = Chisq
            sumht = sumht + ht
            IF ( ABS(ht/sumht).GE..01 ) GOTO 400
         ENDIF
         crit = 0.
         DO jjj = 1 , MEMAX
            crit = crit + (ELM(jjj)-HLMLM(jjj))**2
         ENDDO
         crit = SQRT(crit)
         IF ( crit.LT.Conv ) GOTO 700
         IF ( Chisq.GE.Chiok ) THEN
            rfk = chirf/Chisq
            IF ( rfk.LE.Xtest .OR. icount.GE.Nptl ) GOTO 300
            GOTO 100
         ENDIF
      ENDIF

C     Required chi square achieved       
 600  chil = Chisq
      IF ( Ips.EQ.0 ) WRITE (22,99006) icount
99006 FORMAT (5X,'AT STEP',1X,1I5,1X,'CHISQ CRITERION FULFILLED')
      IF ( Ips.EQ.0 ) WRITE (22,99010) chil
      RETURN
      
 700  IF ( LOCKF.EQ.0 ) THEN ! Terminate if convergence satisfied
         IF ( Chisq.GE.chil ) THEN
            DO jjj = 1 , MEMAX
               ELM(jjj) = ELMH(jjj)
            ENDDO
         ENDIF
         IF ( KFERR.EQ.1 ) RETURN
         IF ( Ips.EQ.0 ) WRITE (22,99007) icount , crit
99007    FORMAT (5X,'AT STEP',1X,1I5,'CONVERGENCE ACHIEVED(',1E14.6,')')
         IF ( Ips.EQ.0 ) WRITE (22,99010) MIN(chil,Chisq)
      ELSE ! Fix most significant chi squared derivatives
         DO kkk = 1 , NLOCK ! NLOCK is number of derivatives to fix
            a = 0.
            iin = 1
            DO jjj = 1 , MEMAX
               IF ( IVAR(jjj).NE.0 .AND. IVAR(jjj).LE.999 ) THEN
                  a = MAX(a,ABS(GRAD(jjj)))
                  IF ( ABS(a-ABS(GRAD(jjj))).LT.1.E-9 ) iin = jjj
               ENDIF
            ENDDO
            IVAR(iin) = 0
            WRITE (22,99011) iin
         ENDDO
         itf = 0
         DO jjj = 1 , MEMAX
            IF ( IVAR(jjj).LE.999 ) THEN
               IF ( IVAR(jjj).NE.0 ) itf = itf + 1
            ENDIF
         ENDDO
         IF ( itf.EQ.1 ) THEN
            metf = 0
            WRITE (22,99008)
99008       FORMAT (2x,'Warning - only one matrix element free',//2x,
     &              'Mode reset to single gradient, execution continues'
     &              ,/)
         ENDIF
         IF ( itf.NE.0 ) GOTO 300
         WRITE (22,99009)
99009    FORMAT (1X/////5X,'*****',2X,'ALL MATRIX ELEMENTS LOCKED!',2X,
     &           '*****'/////)
      ENDIF
      INTR = 0
      RETURN
       
99010 FORMAT (5X,'*** CHISQ=',1E14.6,1X,'***')
99011 FORMAT (1X/5X,'MATRIX ELEMENT',1X,1I3,1X,'LOCKED!')
      END
 
C----------------------------------------------------------------------
C SUBROUTINE CEGRY
C
C Called by: FTBM
C Calls:     ANGULA, DECAY, EFFIX, SIXEL, TACOS
C
C Purpose: calculate the gamma-ray deexcitation.
C
C Uses global variables:
C      AGELI  - angles of the Ge detectors
C      BETAR  - recoil beta
C      CNOR   - normalization factors
C      CORF   - internal correction factors
C      DEV    -
C      DYEX   - error on experimental yield
C      EMH    -
C      ENDEC  - energy difference for each matrix element
C      FIEX   - phi range of particle detector
C      ICLUST - cluster number for each experiment and detector
C      IDRN   - index of normalising transition for yields
C      IEXP   - number of experiment
C      IFMO   - include correction to angular distance for finite recoil distance.
C      IGRD   -
C      ILE    - yield number for each detector
C      IMIN   -
C      INM    - index of matrix element
C      INNR   - independent normalisation switch (see OP,CONT INR,)
C      IPRM   - printing flags (see suboption PRT of OP,CONT)
C      IRAWEX -
C      ITMA   - identify detectors according to OP,GDET
C      ITS    - create tape 18 file (OP,CONT switch SEL,)
C      IWF    - warning flag
C      IY     - index for yields
C      JSKIP  - Experiments to skip during minimisation.
C      KSEQ   - index into ELM for pair of levels, and into EN or SPIN
C      KVAR   -
C      LASTCL - index of last detector in cluster
C      LFL    -
C      LNORM  - normalization constant control
C      LP2    - maximum number of matrix elements (1500)
C      LP6    - maximum number of Ge detectors 32
C      LP10   - maximum number of magnetic substates 1200
C      NANG   - number of gamma-ray detectors for each experiment
C      NDST   - number of data sets
C      NEXPT  - number of experiments
C      NLIFT  - number of lifetimes
C      NMAX   - number of levels
C      NYLDE  - number of yields
C      ODL    - results of OP,GDET calculation
C      SGW    - number of standard deviations to generate warning (see control option WRN,X)
C      SPIN   - spin of level
C      SUBCH1 - partial chisqr
C      SUBCH2 - partial chisqr
C      SUMCL  - sum of yields for clusters
C      TAU    - lifetime in picoseconds
C      TREP   - theta of recoiling nucleus (in radians)
C      UPL    - upper limits for all gamma detectors
C      VACDP  - G_k for each level
C      YEXP   - experimental yield
C      YGN    - gamma yield calculated without correction to angular distribution from finite recoil distance
C      YGP    - gamma yield calculated with correction to angular distribution from finite recoil distance
C      YNRM   - relative normalization factors for gamma detectors
C
C Formal parameters:
C      Chisq  - chi squared
C      Itemp  -
C      Chilo  - chi squared of logs
C      Idr    - number of decays
C      Nwyr   - number of data points contributing to chi squared
C      Icall  -
C      Issp   -
C      Iredv  -
 
      SUBROUTINE CEGRY(Chisq,Itemp,Chilo,Idr,Nwyr,Icall,Issp,Iredv)
      IMPLICIT NONE
      REAL*8 ccc , ccd , Chilo , Chisq , cnr , cocos , d , decen
      REAL*8 dl , effi , fi0 , fi1 , fic , figl , fm , g
      REAL*8 gth , part , partl , rik , rl , rx , ry , 
     &       rys , rz , sf , sgm , sum3 , sumpr , TACOS
      REAL*8 tetrc , tfac , thc , wf
      INTEGER*4 iabc , Icall , id , idc , Idr , ifdu , ifxd , ii , 
     &          ile2 , inclus
      INTEGER*4 ipd , Iredv , Issp , Itemp , iva , iw , ixl , 
     &          ixm , iyex , jj
      INTEGER*4 jj1 , jk , jpc , k , k9 , kc , kj , kk , l , l1 , 
     &          lic , licz , ll1
      INTEGER*4 lth , lu , luu , na , nf , nf1 , ni , ni1 , Nwyr
      CHARACTER*4 wupl , war
      CHARACTER*4 op2
      DIMENSION part(32,50,2) , lic(32) , lth(1500) , cnr(32,50) , 
     &          partl(32,50,2)
      INTEGER*4 ICLUST , LASTCL , IRAWEX
      REAL*8 SUMCL      
      COMMON /CLUST / ICLUST(50,200) , LASTCL(50,20) , SUMCL(20,1500) , 
     &                IRAWEX(50)
      REAL*8 DEV
      COMMON /ODCH  / DEV(1500)
      INTEGER*4 NMAX , NDIM , NMAX1
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      REAL*8 DELTA, ENDEC, ENZ
      INTEGER*4 ITMA
      COMMON /TRA   / DELTA(1500,3) , ENDEC(1500) , ITMA(50,200) , 
     &                ENZ(200)
      REAL*8 BETAR
      COMMON /BREC  / BETAR(50)
      REAL*8 DIX, ODL
      COMMON /DIMX  / DIX(4) , ODL(200)
      REAL*8 VACDP, QCEN, DQ, XNOR, AKS
      INTEGER*4 IBYP
      COMMON /VAC   / VACDP(3,100) , QCEN , DQ , XNOR , AKS(6,100) ,
     &                IBYP
      REAL*8 CNOR
      INTEGER*4 INNR
      COMMON /CINIT / CNOR(32,100) , INNR
      INTEGER*4 IPRM
      COMMON /PRT   / IPRM(20)
      INTEGER*4 NLIFT
      COMMON /LIFE  / NLIFT
      REAL*8 TAU
      INTEGER*4 KSEQ
      COMMON /LEV   / TAU(100) , KSEQ(1500,4)
      INTEGER*4 IGRD
      COMMON /IGRAD / IGRD
      REAL*8 XA , XA1 , EP , TLBDG , VINF
      INTEGER*4 NEXPT , IZ , IZ1
      COMMON /CX    / NEXPT , IZ , XA , IZ1(50) , XA1(50) , EP(50) , 
     &                TLBDG(50) , VINF(50)
      INTEGER*4 IMIN , LNORM
      COMMON /MINNI / IMIN , LNORM(50)
      REAL*8 EMH
      INTEGER*4 INM , LFL1 , LFL2 , LFL
      COMMON /LCZP  / EMH , INM , LFL1 , LFL2 , LFL
      REAL*8 YGN , YGP
      INTEGER*4 IFMO
      COMMON /YTEOR / YGN(1500) , YGP(1500) , IFMO
      INTEGER*4 KVAR
      COMMON /SEL   / KVAR(1500)
      INTEGER*4 LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &          LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      REAL*8 EG, CC, AGELI, Q
      INTEGER*4 NICC, NANG, ISPL
      COMMON /CCC   / EG(50) , CC(50,5) , AGELI(50,200,2) , Q(3,200,8) ,
     &                NICC , NANG(200) , ISPL
      REAL*8 YEXP, CORF , DYEX , UPL , YNRM
      INTEGER*4 IY , NYLDE , IDRN , ILE
      COMMON /YEXPT / YEXP(32,1500) , IY(1500,32) , CORF(1500,32) , 
     &                DYEX(32,1500) , NYLDE(50,32) , UPL(32,50) , 
     &                YNRM(32,50) , IDRN , ILE(32)
      REAL*8 EPS, EROOT, FIEX
      INTEGER*4 IEXP, IAXS
      COMMON /KIN   / EPS(50) , EROOT(50) , FIEX(50,2) , IEXP , IAXS(50)
      REAL*8 SGW , SUBCH1 , SUBCH2
      INTEGER*4 IWF
      COMMON /WARN  / SGW , SUBCH1 , SUBCH2 , IWF
      REAL*8 EN , SPIN , ACCUR , DIPOL , ZPOL , ACCA
      INTEGER*4 ISO
      COMMON /COEX  / EN(100) , SPIN(100) , ACCUR , DIPOL , ZPOL , 
     &                ACCA , ISO
      INTEGER*4 JSKIP
      COMMON /SKP   / JSKIP(50)
      INTEGER*4 ITS
      COMMON /TRB   / ITS
      INTEGER*4 NDST
      COMMON /CCCDS / NDST(50)
      REAL*8 TETACM, TREP, DSIGS
      COMMON /TCM   / TETACM(50) , TREP(50) , DSIGS(50)
      DATA sum3/0./,sumpr/0./

      op2 = '    '
      ifxd = 0
      tetrc = TREP(IEXP) ! Theta of recoiling nucleus

C     If the user set print flag 13 to +1, it is set to -1 by OP,EXIT and then
C     if it is -1, it is set to -2 in MINI, which is called from there, which
C     in turn calls FTBM, which calls this function. In other words, this
C     routine is called with IPRM(13) set to -2 if the user sets IPRM(13) to 1
C     with CONT:PRT, and then does OP,EXIT

      IF ( Icall.EQ.4 .AND. IPRM(13).EQ.-2 ) THEN
         IPRM(13) = 0
         WRITE (22,99001)
99001    FORMAT (1X//20X,'NORMALIZATION CONSTANTS'//2X,'EXPERIMENT',5X,
     &           'DETECTORS(1-32)')
         DO jpc = 1 , NEXPT
            k = NDST(jpc)
            WRITE (22,99012) jpc , (CNOR(l,jpc),l=1,k)
         ENDDO
         WRITE (22,99002)
99002    FORMAT (1X//20X,'RECOMMENDED RELATIVE GE(LI) EFFICIENCIES'//2X,
     &           'EXPERIMENT')
         DO jpc = 1 , NEXPT
            IF ( ABS(cnr(1,jpc)).LT.1.E-9 ) cnr(1,jpc) = 1.
            k = NDST(jpc)
            WRITE (22,99012) jpc , (cnr(l,jpc)/cnr(1,jpc),l=1,k)
         ENDDO ! Loop on experiments
      ENDIF ! if Icall.EQ.4 .AND. IPRM(13).EQ.-2

      DO jpc = 1 , LP6 ! LP6 is 32
         lic(jpc) = 0
      ENDDO

      IF ( Icall.NE.7 ) THEN
         IF ( Itemp.EQ.0 ) THEN
            Nwyr = 0
            IF ( IGRD.NE.1 ) THEN
               IF ( IEXP.EQ.1 ) sumpr = 0.
               IF ( IEXP.EQ.1 ) sum3 = 0.
               DO jj = 1 , LP6 ! LP6 is 32
                  DO jk = 1 , 2
                     partl(jj,IEXP,jk) = 0.
                     part(jj,IEXP,jk) = 0.
                  ENDDO
               ENDDO
            ENDIF

            CALL DECAY(Chisq,NLIFT,Chilo)

            IF ( Icall.EQ.4 .AND. IPRM(14).EQ.-1 ) THEN
               IF ( IEXP.EQ.NEXPT ) IPRM(14) = 0
               WRITE (22,99003)
99003          FORMAT (1X//20X,'VACUUM DEPOLARIZATION COEFFICIENTS '//)
               WRITE (22,99004) IEXP
99004          FORMAT (5X,'EXPERIMENT',1X,1I2/5X,'LEVEL',10X,'G2',10X,
     &                 'G4',10X,'G6'/)
               DO iva = 2 , NMAX
                  WRITE (22,99005) iva , (VACDP(ii,iva),ii=1,3)
99005             FORMAT (7X,1I2,9X,3(1F6.4,6X))
               ENDDO
            ENDIF

            fi0 = FIEX(IEXP,1) ! Lower phi limit
            fi1 = FIEX(IEXP,2) ! Upper phi limit
            na = NANG(IEXP) ! Number of detector angles

            DO k = 1 , LP2 ! LP2 is 1500
               DO k9 = 1 , 20
                  SUMCL(k9,k) = 0.
               ENDDO
            ENDDO

            k9 = 0
            DO k = 1 , na ! For each detector angle
               gth = AGELI(IEXP,k,1) ! theta
               figl = AGELI(IEXP,k,2) ! phi
               ifxd = 0
               fm = (fi0+fi1)/2.
               IF ( Icall.EQ.4 ) ifxd = 1
               CALL ANGULA(YGN,Idr,ifxd,fi0,fi1,tetrc,gth,figl,k,op2)

C              Correct for finite recoil
               IF ( IFMO.NE.0 ) THEN
                  id = ITMA(IEXP,k) ! Get identity for detector
                  d = ODL(id) ! Results of OP,GDET for this detector
                  rx = d*SIN(gth)*COS(figl-fm) - .25*SIN(tetrc)*COS(fm)
                  ry = d*SIN(gth)*SIN(figl-fm) - .25*SIN(tetrc)*SIN(fm)
                  rz = d*COS(gth) - .25*COS(tetrc)
                  rl = SQRT(rx*rx+ry*ry+rz*rz)
                  sf = d*d/rl/rl
                  thc = TACOS(rz/rl)
                  fic = ATAN2(ry,rx)
                  CALL ANGULA(YGP,Idr,ifxd,fi0,fi1,tetrc,thc,fic,k,op2)
                  DO ixl = 1 , Idr ! For each decay
                     ixm = KSEQ(ixl,3) ! Initial level of ixl'th decay
                     tfac = TAU(ixm) ! Get lifetime
                     YGN(ixl) = YGN(ixl) + .01199182*tfac*BETAR(IEXP)
     &                          *(sf*YGP(ixl)-YGN(ixl))
                  ENDDO ! Loop on decays ixl
               ENDIF ! If correction for finite recoil

               IF ( IRAWEX(IEXP).NE.0 ) THEN
                  ipd = ITMA(IEXP,k) ! Get identity for detector
                  DO l = 1 , Idr
                     decen = ENDEC(l)
                     cocos = SIN(tetrc)*SIN(gth)*COS(fm-figl)
     &                       + COS(tetrc)*COS(gth)
                     decen = decen*(1.+BETAR(IEXP)*cocos)
                     CALL EFFIX(ipd,decen,effi)
                     YGN(l) = YGN(l)*effi
                  ENDDO
                  inclus = ICLUST(IEXP,k) ! Cluster number for detector k
                  IF ( inclus.NE.0 ) THEN
                     DO l = 1 , Idr ! For each decay
                        SUMCL(inclus,l) = SUMCL(inclus,l) + YGN(l)
                     ENDDO
                     IF ( k.NE.LASTCL(IEXP,inclus) ) GOTO 20 ! If it is not the last detector in the cluster
                     DO l = 1 , Idr ! For each decay
                        YGN(l) = SUMCL(inclus,l)
                     ENDDO
                  ENDIF
               ENDIF
               k9 = k9 + 1 ! Increment detector number
               IF ( Icall.EQ.4 .AND. IPRM(8).EQ.-2 ) THEN
                  WRITE (22,99006) IEXP , k9
99006             FORMAT (1X//5X,
     &                 'CALCULATED AND EXPERIMENTAL YIELDS   EXPERIMENT'
     &                 ,1X,1I2,1X,'DETECTOR',1X,1I2//6X,'NI',5X,'NF',7X,
     &                 'II',8X,'IF',9X,'ENERGY(MEV)',6X,'YCAL',8X,
     &                 'YEXP',7X,'PC. DIFF.',2X,'(YE-YC)/SIGMA')
               ENDIF
               lu = ILE(k9) ! Yield number for detector k9
               DO iabc = 1 , LP2 ! LP2 = 1500
                  lth(iabc) = 0
               ENDDO
               DO l = 1 , Idr ! For each decay
                  ni = KSEQ(l,3) ! Intial level of l'th decay
                  nf = KSEQ(l,4) ! Final level of l'th decay
                  IF ( l.EQ.IY(lu,k9) .OR. l.EQ.(IY(lu,k9)/1000) ) THEN
                     ifdu = 0
                     lic(k9) = lic(k9) + 1
                     licz = lic(k9)
                     Nwyr = Nwyr + 1
                     wf = CORF(lu,k9)
                     IF ( Icall.EQ.4 ) wf = 1.
                     IF ( Icall.EQ.1 .AND. Issp.EQ.1 ) wf = 1.
                     IF ( IY(lu,k9).GE.1000 ) THEN
                        ifdu = 1
                        l1 = IY(lu,k9)/1000
                        l1 = IY(lu,k9) - 1000*l1
                        YGN(l) = YGN(l) + YGN(l1)
                        lth(l1) = 1
                        IF ( Icall.EQ.4 .AND. IPRM(8).EQ.-2 ) THEN
                           war = '    '
                           sgm = (YEXP(k9,lu)-YGN(l)*CNOR(k9,IEXP))
     &                           /DYEX(k9,lu)
                           IF ( ABS(sgm).GE.SGW ) war = '*?!*'
                           ni1 = KSEQ(l1,3) ! Initial level of l1'th decay
                           nf1 = KSEQ(l1,4) ! Final level of l1'th decay
                           WRITE (22,99007) ni , ni1 , nf , nf1 , 
     &                            SPIN(ni) , SPIN(ni1) , SPIN(nf) , 
     &                            SPIN(nf1) , ENDEC(l) , ENDEC(l1) , 
     &                            YGN(l)*CNOR(k9,IEXP) , YEXP(k9,lu) , 
     &                            100.*(YEXP(k9,lu)-YGN(l)*CNOR(k9,IEXP)
     &                            )/YEXP(k9,lu) , sgm , war
99007                      FORMAT (4X,1I2,'+',1I2,'--',1I2,'+',1I2,3X,
     &                             1F4.1,'+',1F4.1,'--',1F4.1,'+',1F4.1,
     &                             3X,1F6.4,'+',1F6.4,2X,1E9.4,6X,1E9.4,
     &                             3X,1F6.1,5X,1F4.1,10X,1A4)
                           SUBCH1 = SUBCH1 + sgm*sgm
                        ENDIF
                     ENDIF
                     ry = YGN(l)*wf*CNOR(k9,IEXP) - YEXP(k9,lu)
                     IF ( ifdu.NE.1 ) THEN
                        IF ( Icall.EQ.4 .AND. IPRM(8).EQ.-2 ) THEN
                           war = '    '
                           sgm = (YEXP(k9,lu)-YGN(l)*CNOR(k9,IEXP))
     &                           /DYEX(k9,lu)
                           IF ( ABS(sgm).GE.SGW ) war = '*?!*'
                           WRITE (22,99013) ni , nf , SPIN(ni) , 
     &                            SPIN(nf) , ENDEC(l) , YGN(l)
     &                            *CNOR(k9,IEXP) , YEXP(k9,lu) , 
     &                            100.*(YEXP(k9,lu)-YGN(l)*CNOR(k9,IEXP)
     &                            )/YEXP(k9,lu) , sgm , war
                           SUBCH1 = SUBCH1 + sgm*sgm
                        ENDIF
                     ENDIF
                     rys = ry*ry
                     IF ( IGRD.EQ.1 ) Chisq = Chisq + rys/DYEX(k9,lu)
     &                    /DYEX(k9,lu)
                     IF ( k9.EQ.1 .AND. Iredv.EQ.1 ) DEV(licz) = ry
                     IF ( Iredv.NE.1 ) THEN
                        IF ( LFL.EQ.1 ) THEN
                           IF ( k9.EQ.1 ) THEN
                              luu = 6*licz - 5
                              jk = (luu-1)/LP10 + 1 ! LP10 is 1200
                              kk = luu - LP10*(jk-1)
                              rik = DEV(licz) + YEXP(k9,lu)
                              sgm = -DEV(licz)/DYEX(k9,lu)
                              IF ( ITS.EQ.1 .AND. KVAR(INM).NE.0 )
     &                             WRITE (17,*) ni , nf , sgm , YGN(l)
     &                             *CNOR(k9,IEXP)/DYEX(k9,lu)
                              IF ( ITS.EQ.1 .AND. INM.EQ.1 )
     &                             WRITE (15,*) IEXP , rik/CNOR(1,IEXP)
     &                             , CNOR(1,IEXP) , DYEX(k9,lu) , 
     &                             YEXP(k9,lu)
                              CALL SIXEL(rik,ry,EMH,jk,kk,INM,licz)
                           ENDIF
                        ENDIF
                     ENDIF
                     IF ( IGRD.NE.1 ) THEN
                        IF ( JSKIP(IEXP).NE.0 ) THEN
                           dl = DYEX(k9,lu)*DYEX(k9,lu)
                           part(k9,IEXP,1) = part(k9,IEXP,1) + YGN(l)
     &                        *YGN(l)*wf*wf/dl
                           part(k9,IEXP,2) = part(k9,IEXP,2) - 2.*YGN(l)
     &                        *wf*YEXP(k9,lu)/dl
                           sumpr = sumpr + YEXP(k9,lu)*YEXP(k9,lu)/dl
                           partl(k9,IEXP,1) = partl(k9,IEXP,1)
     &                        + YEXP(k9,lu)*YEXP(k9,lu)/dl
                           partl(k9,IEXP,2) = partl(k9,IEXP,2)
     &                        + LOG(wf*YGN(l)/YEXP(k9,lu))*YEXP(k9,lu)
     &                        *YEXP(k9,lu)/dl
                           sum3 = sum3 + YEXP(k9,lu)*YEXP(k9,lu)
     &                            *LOG(wf*YGN(l)/YEXP(k9,lu))**2/dl
                        ENDIF
                     ENDIF
                     lu = lu + 1
                  ELSE
                     IF ( JSKIP(IEXP).EQ.0 ) YGN(IDRN) = 1.E+10
                     ry = YGN(l)/YGN(IDRN)
                     IF ( Icall.EQ.4 .AND. IPRM(8).EQ.-2 ) THEN
                        wupl = '    '
                        IF ( ry.GT.UPL(k9,IEXP) .AND. lth(l).EQ.0 )
     &                       wupl = 'UPL!'
                        IF ( IPRM(16).NE.0 .OR. wupl.NE.'    ' ) THEN
                           IF ( wupl.EQ.'    ' ) WRITE (22,99008) ni , 
     &                          nf , SPIN(ni) , SPIN(nf) , ENDEC(l) , 
     &                          YGN(l)*CNOR(k9,IEXP) , wupl
99008                      FORMAT (6X,1I2,5X,1I2,7X,1F4.1,6X,1F4.1,9X,
     &                             1F6.4,6X,1E9.4,10X,1A4)
                           IF ( wupl.NE.'    ' ) THEN
                              sgm = (ry-UPL(k9,IEXP))/UPL(k9,IEXP)
                              WRITE (22,99013) ni , nf , SPIN(ni) , 
     &                               SPIN(nf) , ENDEC(l) , YGN(l)
     &                               *CNOR(k9,IEXP) , UPL(k9,IEXP)
     &                               *CNOR(k9,IEXP)*YGN(IDRN) , 
     &                               100.*(1.-YGN(l)/UPL(k9,IEXP)
     &                               /YGN(IDRN)) , sgm , wupl
                              SUBCH1 = SUBCH1 + sgm*sgm
                           ENDIF
                        ENDIF
                     ENDIF
                     IF ( ry.GE.UPL(k9,IEXP) .AND. lth(l).NE.1 ) THEN
                        Chisq = Chisq + (ry-UPL(k9,IEXP))
     &                          *(ry-UPL(k9,IEXP))/UPL(k9,IEXP)
     &                          /UPL(k9,IEXP)
                        Chilo = Chilo + LOG(ry/UPL(k9,IEXP))**2
                        IF ( IWF.NE.0 ) THEN ! If warning flag is set
                           WRITE (22,99009) IEXP , ni , nf , 
     &                            ry/UPL(k9,IEXP)
99009                      FORMAT (5X,'WARNING-EXP.',1I2,2X,'TRANS. ',
     &                             1I2,'--',1I2,5X,
     &                             'EXCEEDS UPPER LIMIT (RATIO=',1E14.6,
     &                             ')')
                        ENDIF
                     ENDIF
                  ENDIF
               ENDDO ! Loop on decays l
               IF ( IEXP.EQ.NEXPT ) IWF = 0 ! Turn off warnings now
               IF ( Icall.EQ.4 .AND. IPRM(8).EQ.-2 ) THEN
                  WRITE (22,99010) SUBCH1 - SUBCH2
99010             FORMAT (1X/50X,'CHISQ SUBTOTAL = ',E14.6)
                  SUBCH2 = SUBCH1
               ENDIF
 20            CONTINUE
            ENDDO ! Loop on detector angles k

            IF ( IGRD.EQ.1 ) RETURN
            IF ( IEXP.NE.NEXPT ) RETURN
            IF ( Icall.EQ.1 ) RETURN
         ELSE
            ifxd = 1
            IF ( Itemp.NE.2 ) ifxd = 0
            Nwyr = 1
            CALL DECAY(ccd,0,ccc)
            fi0 = FIEX(IEXP,1)
            fi1 = FIEX(IEXP,2)
            na = NANG(IEXP)
            DO k = 1 , LP2 ! LP2 is 1500
               DO kj = 1 , 20
                  SUMCL(kj,k) = 0
               ENDDO
            ENDDO
            k9 = 0
            DO k = 1 , na
               gth = AGELI(IEXP,k,1)
               figl = AGELI(IEXP,k,2)
               fm = (fi0+fi1)/2.

               CALL ANGULA(YGN,Idr,ifxd,fi0,fi1,tetrc,gth,figl,k,op2)

C              Correct for finite recoil
               IF ( IFMO.NE.0 ) THEN
                  id = ITMA(IEXP,k) ! Get identity for detector
                  d = ODL(id) ! Get results of OP,GDET for detector
                  rx = d*SIN(gth)*COS(figl-fm) - .25*SIN(tetrc)*COS(fm)
                  ry = d*SIN(gth)*SIN(figl-fm) - .25*SIN(tetrc)*SIN(fm)
                  rz = d*COS(gth) - .25*COS(tetrc)
                  rl = SQRT(rx*rx+ry*ry+rz*rz)
                  sf = d*d/rl/rl
                  thc = TACOS(rz/rl)
                  fic = ATAN2(ry,rx)
                  CALL ANGULA(YGP,Idr,ifxd,fi0,fi1,tetrc,thc,fic,k,op2)
                  DO ixl = 1 , Idr
                     ixm = KSEQ(ixl,3) ! Initial level of ixl'th decay
                     tfac = TAU(ixm)
                     IF ( tfac.GT.1.E+4 ) GOTO 25
                     YGN(ixl) = YGN(ixl) + .01199182*tfac*BETAR(IEXP)
     &                          *(sf*YGP(ixl)-YGN(ixl))
                  ENDDO
 25               IFMO = 0
                  WRITE (22,99011)
99011             FORMAT (1X,/,2X,'DURING THE MINIMIZATION',1X,
     &    'IT WAS NECESSARY TO SWITCH OFF THE TIME-OF-FLIGHT CORRECTION'
     &    )
               ENDIF ! if correction for finite recoil

               IF ( IRAWEX(IEXP).NE.0 ) THEN
                  ipd = ITMA(IEXP,k) ! Get identity of detector
                  DO l = 1 , Idr
                     decen = ENDEC(l)
                     cocos = SIN(tetrc)*SIN(gth)*COS(fm-figl)
     &                       + COS(tetrc)*COS(gth)
                     decen = decen*(1.+BETAR(IEXP)*cocos)
                     CALL EFFIX(ipd,decen,effi)
                     YGN(l) = YGN(l)*effi
                  ENDDO
                  inclus = ICLUST(IEXP,k) ! Cluster number for detector k
                  IF ( inclus.NE.0 ) THEN
                     DO l = 1 , Idr ! For each decay
                        SUMCL(inclus,l) = SUMCL(inclus,l) + YGN(l)
                     ENDDO
                     IF ( k.NE.LASTCL(IEXP,inclus) ) GOTO 40 ! If it is not the last detector in the cluster
                     DO l = 1 , Idr ! For each decay
                        YGN(l) = SUMCL(inclus,l)
                     ENDDO
                  ENDIF
               ENDIF
               k9 = k9 + 1
               iyex = NYLDE(IEXP,k9) + ILE(k9) - 1
               ile2 = ILE(k9)
               DO l = ile2 , iyex
                  IF ( JSKIP(IEXP).NE.0 ) THEN
                     idc = IY(l,k9)
                     IF ( idc.GE.1000 ) THEN
                        idc = idc/1000
                        ll1 = IY(l,k9) - idc*1000
                        YGN(idc) = YGN(idc) + YGN(ll1)
                     ENDIF
                     IF ( Itemp.EQ.1 ) THEN
                        CORF(l,k9) = CORF(l,k9)/(YGN(idc)+1.E-24)
                     ELSE
                        CORF(l,k9) = YGN(idc)
                        IF ( IMIN.LE.1 .AND. l.EQ.iyex ) CNOR(k9,IEXP)
     &                       = YEXP(k9,l)/YGN(idc)
                     ENDIF
                  ENDIF
               ENDDO ! Loop on l
 40            CONTINUE
            ENDDO ! Loop on k
            RETURN
         ENDIF ! if Itemp.EQ.0
      ENDIF ! if Icall.NE.7

C     Sort out normalisation coefficients
      DO jj = 1 , NEXPT ! For each experiment
         IF ( JSKIP(jj).NE.0 ) THEN
            kc = NDST(jj) ! Number of datasets
            DO jk = 1 , kc ! For each dataset
               cnr(jk,jj) = -.5*part(jk,jj,2)/part(jk,jj,1)
               IF ( INNR.NE.0 ) CNOR(jk,jj) = cnr(jk,jj)
            ENDDO ! Loop on datasets

C           If we want a common normalisation, sort it out here
            IF ( INNR.NE.1 ) THEN
               d = 0.
               g = 0.
               DO jj1 = jj , NEXPT ! For each experiment
                  IF ( LNORM(jj1).EQ.jj ) THEN
                     k = NDST(jj1) ! Number of datasets
                     DO jk = 1 , k ! For each dataset
                        d = d + YNRM(jk,jj1)*part(jk,jj1,1)*YNRM(jk,jj1)
                        g = g - .5*YNRM(jk,jj1)*part(jk,jj1,2)
                     ENDDO ! Loop on datasets
                  ENDIF ! IF ( LNORM(jj1).EQ.jj )
               ENDDO ! Loop on experiment
               IF ( LNORM(jj).EQ.jj ) THEN ! If this is the normalisation transition
                  CNOR(1,jj) = g*YNRM(1,jj)/d
                  k = NDST(jj) ! Number of datasets
                  IF ( k.NE.1 ) THEN
                     DO jk = 2 , k ! For each dataset
                        CNOR(jk,jj) = YNRM(jk,jj)*CNOR(1,jj)/YNRM(1,jj)
                     ENDDO ! Loop on jk
                  ENDIF ! IF ( k.NE.1 )
               ENDIF ! IF ( LNORM(jj).EQ.jj )
            ENDIF ! IF ( INNR.NE.1 )
             
         ENDIF ! IF ( JSKIP(jj).NE.0 )
      ENDDO ! Loop on experiment

C     If there is a common normalisation, normalise to it
      IF ( INNR.NE.1 ) THEN
         DO jj = 1 , NEXPT ! For each experiment
            IF ( LNORM(jj).NE.jj ) THEN
               iw = LNORM(jj) ! Get index of normalisation transition
               k = NDST(jj) ! Get number of datasets
               DO jk = 1 , k ! For each dataset
                  CNOR(jk,jj) = CNOR(1,iw)*YNRM(jk,jj)/YNRM(1,iw)
               ENDDO ! Loop on datasets
            ENDIF ! IF ( LNORM(jj).NE.jj )
         ENDDO ! Loop on experiment
      ENDIF ! IF ( INNR.NE.1 )

C     Calculate chi squared
      IF ( Icall.EQ.7 ) Chisq = 0.
      DO jj = 1 , NEXPT
         k = NDST(jj)
         DO jk = 1 , k
            Chilo = Chilo + partl(jk,jj,1)*LOG(CNOR(jk,jj))
     &              **2 + partl(jk,jj,2)*2.*LOG(CNOR(jk,jj))
            Chisq = Chisq + CNOR(jk,jj)*CNOR(jk,jj)*part(jk,jj,1)
     &              + CNOR(jk,jj)*part(jk,jj,2)
         ENDDO ! Loop on datasets
      ENDDO ! Loop on experiment

      Chisq = Chisq + sumpr
      Chilo = Chilo + sum3
      RETURN

99012 FORMAT (1X,1I2,2X,32(1E10.4,1X))
99013 FORMAT (6X,1I2,5X,1I2,7X,1F4.1,6X,1F4.1,9X,1F6.4,6X,1E9.4,6X,
     &        1E9.4,3X,1F6.1,5X,1F4.1,10X,1A4)
      END
 
C----------------------------------------------------------------------
C SUBROUTINE FAKP
C
C Called by: GOSIA
C Calls:     PRIM
C
C Purpose: calculate log of primes and the factoring of primes
C
C Uses global variables:
C      IP     - table of primes
C      IPI    - number of time a number is divisible by each prime in IP
C      KF     - sum of factors of primes
C      PILOG  - table of natural logs of primes
 
      SUBROUTINE FAKP
      IMPLICIT NONE
      INTEGER*4 i , k , l
      REAL*8 x
      REAL*8 PILOG
      INTEGER*4 IP , IPI , KF
      COMMON /FAKUL / IP(26) , IPI(26) , KF(101,26) , PILOG(26)
C     Calculate log of primes
      DO i = 1 , 26
         x = DBLE(IP(i))
         PILOG(i) = LOG(x)
      ENDDO

C     Initialise KF       
      DO l = 1 , 26
         KF(1,l) = 0
         KF(2,l) = 0
      ENDDO

C     Calculate factors of numbers
      DO k = 3 , 101
         CALL PRIM(k-1) ! Puts factors in IPI array
         DO i = 1 , 26
            KF(k,i) = KF(k-1,i) + IPI(i) ! IPI are number of times prime is a factor
         ENDDO
      ENDDO
      END
 
C----------------------------------------------------------------------
C SUBROUTINE PRIM
C
C Called by: FAKP
C
C Purpose: Find out how many times each prime divides a number N
C
C Uses global variables:
C      IP     - table of primes
C      IPI    - multipliers for each prime
C
C Formal parameters:
C      N      - number N
 
      SUBROUTINE PRIM(N)
      IMPLICIT NONE
      INTEGER*4 i , N , nni , nnk
      REAL*8 PILOG
      INTEGER*4 IP , IPI , KF
      COMMON /FAKUL / IP(26) , IPI(26) , KF(101,26) , PILOG(26)
      nnk = N
      DO i = 1 , 26
         nni = nnk
         IPI(i) = 0
 50      nni = nni/IP(i)
         IF ( IP(i)*nni.EQ.nnk ) THEN
            IPI(i) = IPI(i) + 1
            nnk = nni
            GOTO 50
         ENDIF
      ENDDO
      END
 
C----------------------------------------------------------------------
C SUBROUTINE SEQ
C
C Called by: ADHOC
C Calls:     CONV, F, GF, LEADF, MEM
C
C Purpose: in order to calculate the yields, we need to start with the highest
C level and calculate its yield, so as to work out the feeding for the lower
C levels and take this into account, gradually working our way down to the
C ground state.
C
C Uses global variables:
C      DELTA  - \delta_\lambda: index 1 = electric^2, 2 = magnetic^2, 3 = cross term
C      EN     - energy of level
C      ENDEC  - energy difference for each matrix element
C      FP     - F coefficient * DELTA^2
C      GKP    - Gk * DELTA^2
C      KLEC   - number of decays for each level
C      KSEQ   - indices for each decay (level1, level2, matrix element, multipolarity + 10)
C      LDNUM  - number of matrix elements with each multipolarity populating levels
C      LP2    - maximum number of matrix elements (1500)
C      LP3    - maximum number of levels (100)
C      MULTI  - number of matrix elements having a given multipolarity
C      NMAX   - number of levels
C      NMAX1  - number of levels with decays
C      SPIN   - spin of level
C      TAU    - normally lifetime in picoseconds (here it is used for energies, however)
C
C Formal parameters:
C      Idr    - returns number of items in KSEQ array.
C
C We store the order in the KSEQ array of common block LEV.
C
C Note that in the code, a multipolarity 1 = E1, 2 = E2 ... 6 = E6, 7 = M1,
C 8 = M2.
 
      SUBROUTINE SEQ(Idr)
      IMPLICIT NONE
      REAL*8 CONV , ega , egs , emax , F , GF , spinf , spini , twoi
      INTEGER*4 idecay , Idr , indx , inx , inx1 , ir , is , istr1 , 
     &          istr2 , j , js , jsave , k , kpa , l , la , la1 , 
     &          ld , LEADF
      INTEGER*4 m , m1 , m6 , MEM , mk , mule , mulm , n , n1 , nob
      INTEGER*4 NMAX , NDIM , NMAX1
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      REAL*8 DELTA, ENDEC, ENZ
      INTEGER*4 ITMA
      COMMON /TRA   / DELTA(1500,3) , ENDEC(1500) , ITMA(50,200) , 
     &                ENZ(200)
      INTEGER*4 LAMDA , LEAD , LDNUM , LAMMAX , MULTI
      COMMON /CLCOM / LAMDA(8) , LEAD(2,1500) , LDNUM(8,100) , LAMMAX , 
     &                MULTI(8)
      INTEGER*4 LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &          LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      REAL*8 EN , SPIN , ACCUR , DIPOL , ZPOL , ACCA
      INTEGER*4 ISO
      COMMON /COEX  / EN(100) , SPIN(100) , ACCUR , DIPOL , ZPOL , 
     &                ACCA , ISO
      REAL*8 TAU
      INTEGER*4 KSEQ
      COMMON /LEV   / TAU(100) , KSEQ(1500,4)
      REAL*8 FP , GKP
      INTEGER*4 KLEC
      COMMON /CATLF / FP(4,1500,3) , GKP(4,1500,2) , KLEC(100)
      DATA jsave/0/
      
      m6 = 0
      DO l = 1 , 6
         m6 = m6 + MULTI(l)
      ENDDO

      idecay = 0
      Idr = 0

      DO l = 1 , LP3 ! LP3 = 100 (number of levels)
         KLEC(l) = 0 ! Initialise KLEC to zero
      ENDDO

      DO k = 1 , LP2 ! LP2 = 1500 (number of matrix elements)
         DO j = 1 , 3
            DO l = 1 , 4
               FP(l,k,j) = 0.
               IF ( j.NE.3 ) GKP(l,k,j) = 0.
            ENDDO
            DELTA(k,j) = 0.
         ENDDO
      ENDDO

C     Store the energies in TAU array
      DO n = 1 , NMAX
         TAU(n) = EN(n)
      ENDDO

      DO n = 1 , NMAX ! Loop on levels
C        Find level with highest energy
         emax = 0.
         DO j = 1 , NMAX ! Loop on levels
            IF ( TAU(j).GE.emax ) THEN
               emax = TAU(j)
               jsave = j
            ENDIF
         ENDDO
         DO is = 1 , NMAX ! Loop on levels
            DO la = 1 , 8 ! Loop on multipolarities
               IF ( la.LE.3 .OR. la.EQ.7 .OR. la.EQ.8 ) THEN ! E3, M1, M2
                  ld = LDNUM(la,is) ! Number of levels connected to this one with this multipolarity
                  IF ( ld.NE.0 ) THEN
                     DO ir = 1 , ld ! For each level ir connected to level is with multipolarity la
                        m = LEADF(is,ir,la)
                        IF ( m.EQ.jsave .OR. is.EQ.jsave ) THEN
                           IF ( is.NE.jsave .OR. EN(m).LT.EN(is) ) THEN
                              IF ( m.NE.jsave .OR. EN(is).LT.EN(m) )
     &                             THEN
                                 indx = MEM(is,m,la) ! Matrix element from level is to level m with multipolarity la
                                 idecay = idecay + 1
                                 KSEQ(idecay,1) = m       ! Level
                                 KSEQ(idecay,2) = is      ! Level
                                 KSEQ(idecay,3) = indx    ! Matrix element
                                 KSEQ(idecay,4) = la + 10 ! Multipolarity + 10
                                 IF ( EN(m).LE.EN(is) ) THEN ! If the levels are degenerate, swap order
                                    KSEQ(idecay,1) = is
                                    KSEQ(idecay,2) = m
                                 ENDIF
                              ENDIF
                           ENDIF
                        ENDIF
                     ENDDO ! Loop on levels ir
                  ENDIF
               ENDIF
            ENDDO ! Loop on multipolarity la
         ENDDO ! Loop on levels is
         TAU(jsave) = -1.
      ENDDO ! Loop on levels n

C     Now for each decay, calculate transition amplitudes for each
C     multipolarity
      DO l = 1 , idecay ! For each decay
         istr1 = 0
         IF ( KSEQ(l,4).LT.10 ) GOTO 200 ! KSEQ(l,4) is 10 + multipolarity
         istr2 = 0
         n = KSEQ(l,1) ! Initial level
         m = KSEQ(l,2) ! Final level
         inx = KSEQ(l,3) ! Matrix element
         la = KSEQ(l,4) - 10 ! Multipolarity
         ega = EN(n) - EN(m)    ! ega = E_\gamma
         twoi = 1./SQRT(2.*SPIN(n)+1.)
         spini = SPIN(n) + .001
         spinf = SPIN(m) + .001
         egs = SQRT(ega)*twoi   ! egs = \sqrt{E_\gamma \over 2 I_1 + 1}
         js = l + 1
         la1 = 0
         inx1 = 0
         DO j = js , idecay ! For each decay
            IF ( KSEQ(j,4).GE.10 ) THEN ! KSEQ(j,4) is 10 + multipolarity
               n1 = KSEQ(j,1) ! Initial level
               m1 = KSEQ(j,2) ! Final level
               IF ( n1.EQ.n .AND. m1.EQ.m ) THEN ! Decays involving the same pair of levels
                  inx1 = KSEQ(j,3) ! Matrix element
                  la1 = KSEQ(j,4) - 10 ! Multipolarity
                  KSEQ(j,4) = KSEQ(j,4) - 10 ! Subtract ten to indicate we have handled this one
               ENDIF
            ENDIF
         ENDDO ! Loop on decays j
         KSEQ(l,4) = KSEQ(l,4) - 10 ! Subtract ten to indicate we have handled this one
         Idr = Idr + 1
         mule = 0
         mulm = 0
         nob = 1
 50      IF ( la.LE.3 ) THEN
            IF ( la.EQ.1 ) THEN
               DELTA(Idr,1) = 398.77393*ega*egs ! E1
               mule = 1
               istr1 = 1 ! In array CC and N parameter of CONV -> E1
            ELSEIF ( la.EQ.2 ) THEN
               DELTA(Idr,1) = 3.5002636*egs*ega*ega ! E2
               mule = 2
               istr1 = 2 ! In array CC and N parameter of CONV -> E2
            ELSEIF ( la.EQ.3 ) THEN
               DELTA(Idr,1) = 0.023891302*ega*ega*ega*egs ! E3
               mule = 3
               istr1 = 3 ! In array CC and N parameter of CONV -> E3
            ELSE
               GOTO 100
            ENDIF
            GOTO 150
         ENDIF
 100     la = la - 6
         IF ( la.EQ.2 ) THEN
            DELTA(Idr,2) = 0.036806836*ega*ega*egs ! M2
            mulm = 2
            istr2 = 5 ! In array CC and N parameter of CONV -> M2
         ELSE
            DELTA(Idr,2) = 4.1932861*ega*egs ! M1
            mulm = 1
            istr2 = 4 ! In array CC and N parameter of CONV -> M1
         ENDIF
 150     IF ( nob.NE.2 ) THEN
            IF ( mule.NE.1 ) THEN
               nob = nob + 1
               IF ( la.GT.3 ) inx1 = inx
               IF ( la1.NE.0 ) THEN
                  la = la1
                  GOTO 50
               ENDIF
            ENDIF
            inx1 = 0
         ENDIF
         DELTA(Idr,3) = DELTA(Idr,1)*DELTA(Idr,2)
         DELTA(Idr,1) = DELTA(Idr,1)*DELTA(Idr,1)
         DELTA(Idr,2) = DELTA(Idr,2)*DELTA(Idr,2)
         KSEQ(Idr,1) = inx
         KSEQ(Idr,2) = inx1
         KSEQ(Idr,3) = n
         KSEQ(Idr,4) = m
         IF ( inx.GT.m6 ) THEN
            KSEQ(Idr,2) = inx
            KSEQ(Idr,1) = 0
         ENDIF
         ENDEC(Idr) = EN(n) - EN(m) ! Energy difference between levels
         DO mk = 1 , 7 , 2
            kpa = mk/2 + 1
            k = mk - 1
            IF ( mule.GE.3 .OR. k.NE.6 ) THEN
               GKP(kpa,Idr,1) = GF(k,spini,spinf,mule)*DELTA(Idr,1)
     &                          *(1.+CONV(ega,istr1))
               GKP(kpa,Idr,2) = GF(k,spini,spinf,mulm)*DELTA(Idr,2)
     &                          *(1.+CONV(ega,istr2))
               FP(kpa,Idr,1) = F(k,spini,spinf,mule,mule)*DELTA(Idr,1)
               FP(kpa,Idr,3) = F(k,spini,spinf,mulm,mule)*DELTA(Idr,3)
               FP(kpa,Idr,2) = F(k,spini,spinf,mulm,mulm)*DELTA(Idr,2)
            ENDIF
         ENDDO ! Loop on mk
         DELTA(Idr,1) = DELTA(Idr,1)*(1.+CONV(ega,istr1))
         DELTA(Idr,2) = DELTA(Idr,2)*(CONV(ega,istr2)+1.)
         KLEC(n) = KLEC(n) + 1 ! Increment KLEC for initial level
 200     CONTINUE
      ENDDO ! Loop on decays l

      NMAX1 = 0
      DO n = 1 , NMAX ! For each level count those which have decays
         IF ( KLEC(n).NE.0 ) NMAX1 = NMAX1 + 1
      ENDDO
      END
 
C----------------------------------------------------------------------
C FUNCTION GF
C
C Called by: SEQ
C Calls:     WSIXJ
C
C Purpose: calculate the H_k coefficients to modify the statistical tensors
C to take feeding due to multiple excitation into account.
C
C Formal parameters:
C      K      - K
C      Sji    - initial spin
C      Sjf    - final spin
C      L      - L
C
C Return value:
C      H_k coefficient
C
C Note that the parameters to WSIXJ need to be doubled, so that this function
C can handle half-integers.

      REAL*8 FUNCTION GF(K,Sji,Sjf,L)
      IMPLICIT NONE
      INTEGER*4 i , ix , jfz , jiz , K , kz , L , lz
      REAL*8 phase , Sjf , Sji , WSIXJ
      
      GF = 0.
      IF ( L.EQ.0 ) RETURN
      ix = INT(Sji+Sjf+.0001)
      i = ix + L + K
      phase = 1.
      IF ( i/2*2.NE.i ) phase = -1.
      kz = K*2
      jiz = Sji*2
      jfz = Sjf*2
      lz = L*2
      GF = phase*SQRT((jiz+1.)*(jfz+1.))*WSIXJ(jiz,jiz,kz,jfz,jfz,lz)
      END
 
C----------------------------------------------------------------------
C FUNCTION F
C
C Called by: SEQ
C Calls:     WSIXJ, WTHREJ
C
C Purpose: evaluates the F coefficients.
C
C Formal coefficients:
C      K      - K
C      Sji    - initial spin
C      Sjf    - final spin
C      L1     - lambda
C      L2     - lambda'
C
C Return value:
C      F-coefficient
C
C We evaluate:
C F_k(\lambda \lambda^\prime I_2 I1) = (-1)^{I_1 + I_2 -l} *
C        \sqrt{(2 k + 1) (2 I_1 + 1) (2 \lambda + 1) (2 \lambda^\prime + 1) *
C        \threej{\lambda \lambda^\prime k 1 -1 0} *
C        \sixj{\lambda \lambda^\prime k I_1 I_1 I_2}
C
C Here \lambda = L1, \lambda^\prime = L2, I_1 = Sji, I_2 = Sjf, k = K
C
C Note that the code actually evaluates:
C \sixj{I_1 I_1 k \lambda^\prime \lambda I_2} which is equal to
C \sixj{\lambda \lambda^\prime k I_1 I_1 I_2} by the symmetry rules for 6-j
C symbols.
C
C Note also that both WTHREJ and WSIXJ need to have parameters which are twice
C the values to calculate, so that they can handle half-integers correctly.

      REAL*8 FUNCTION F(K,Sji,Sjf,L1,L2)
      IMPLICIT NONE
      INTEGER*4 ix , jfz , jiz , K , kz , l , L1 , l1z , L2 , l2z
      REAL*8 phase , Sjf , Sji , WSIXJ , WTHREJ
      
      F = 0.
      IF ( (L1*L2).EQ.0 ) RETURN
      ix = INT(Sji+Sjf+.0001)
      l = ix - 1
      phase = 1.
      IF ( l/2*2.NE.l ) phase = -1.
      kz = K*2
      jiz = Sji*2
      jfz = Sjf*2
      l1z = L1*2
      l2z = L2*2
      F = phase*SQRT((l1z+1.)*(l2z+1.)*(jiz+1.)*(kz+1.))
     &    *WTHREJ(l1z,l2z,kz,2,-2,0)*WSIXJ(jiz,jiz,kz,l2z,l1z,jfz)
      END
 
C----------------------------------------------------------------------
C FUNCTION CONV
C
C Called by: BRANR, PTICC, SEQ
C Calls:     LAGRAN, NEWCNV, SPLNER
C
C Purpose: calculate the conversion coefficient at a particular energy by
C interpolating over the values provided by the user.
C
C Uses global variables:
C      CC     - conversion coefficients
C      EG     - energies for conversion coefficients
C      NICC   - number of conversion coefficients
C
C Formal parameters:
C      Ega    - gamma energy
C      N      - multipolarity N=1,2,3 = E1,2,3 and N=4,5 = M1,2 (not as elsewhere!)
C
C Return value:
C      conversion coefficient interpolated to energy Ega

      REAL*8 FUNCTION CONV(Ega,N)
      IMPLICIT NONE
      REAL*8 cpo , cpo1 , cv , Ega , NEWCNV
      INTEGER*4 j , N , n1 , nen
      DIMENSION cpo(101) , cpo1(101)
      REAL*8 EG, CC, AGELI, Q
      INTEGER*4 NICC, NANG, ISPL
      COMMON /CCC   / EG(50) , CC(50,5) , AGELI(50,200,2) , Q(3,200,8) ,
     &                NICC , NANG(200) , ISPL
C     If the number of conversion coefficients entered by the user is negative
C     then use read the conversion coefficients from a file on unit 29.
      IF ( NICC.LE.0 ) THEN
         CONV=NEWCNV(Ega,N)
         RETURN
      ENDIF

      IF ( N.EQ.0 ) THEN ! If no multipolarity defined
         CONV = 0.0
      ELSEIF ( ABS(CC(1,N)).LT.1.E-9 ) THEN ! If no conversion coefficients given for this multipolarity
         CONV = 0.0
      ELSE
         nen = 4
         DO j = 1 , NICC ! Loop over coefficients provided by user
            IF ( Ega.LE.EG(j) ) GOTO 50
         ENDDO
 50      n1 = j - 2
         IF ( n1.LT.1 ) n1 = 1
         IF ( (j+1).GT.NICC ) n1 = n1 - 1
         IF ( NICC.LE.4 ) THEN
            n1 = 1
            nen = NICC
         ENDIF
         DO j = 1 , nen
            cpo(j) = CC(n1+j-1,N)
            cpo1(j) = EG(n1+j-1)
         ENDDO
C        Interpolate 
         IF ( ISPL.EQ. 0 ) CALL LAGRAN(cpo1,cpo,4,1,Ega,cv,2,1)
         IF ( ISPL.EQ. 1 ) CALL SPLNER(cpo1,cpo,4,Ega,cv,2)
         CONV = cv
         RETURN
      ENDIF
      END
 
C----------------------------------------------------------------------
C FUNCTION WTHREJ
C
C Called by: ELMT, F, GOSIA, LSLOOP, TENB
C
C Purpose: evaluates a Wigner 3-j symbol.
C
C Uses global variables:
C      IP     - table of prime numbers
C      KF     - sum of factors of primes
C      PILOG  - table of natural logs of primes
C
C Formal parameters:
C      J1     - twice the value of J1
C      J2     - twice the value of J2
C      J3     - twice the value of J3
C      M1     - twice the value of M1
C      M2     - twice the value of M2
C      M3     - twice the value of M3
C
C Return value:
C      The value of the 3-j symbol
C
C Note that the values of the parameters are doubled, so that this function
C can handle half-integers. In other words if you want to evaluate
C \threej(J1 J2 J3 M1 M2 M3) you need to use call the function as:
C WTHREJ(2 * J1, 2 * J2, 2 * J3, 2 * M1, 2 * M2, 2 * M3).
 
      REAL*8 FUNCTION WTHREJ(J1,J2,J3,M1,M2,M3)
      IMPLICIT NONE
      INTEGER*4 iz , iza , izb , izc , izd , ize , izexp , 
     &          izf , izmax , izmin , J1 , J2 , J3 , jabc , jabm , 
     &          jbma , jj1 , jj2
      INTEGER*4 jj3 , jjha , jjhb , jjhc , jjhd , jlp , jma , jmax , 
     &          jmb , jmc , jmd , jme , jmf , jta , jtb , jtc , jvo , 
     &          jvora , M1
      INTEGER*4 M2 , M3 , mm1 , mm2 , mm3 , n , nmax
      REAL*8 qsumlo , sumlo , vorz , wthrep , zuthre
      DIMENSION jvora(26)
      REAL*8 PILOG
      INTEGER*4 IP , IPI , KF
      COMMON /FAKUL / IP(26) , IPI(26) , KF(101,26) , PILOG(26)
      
      wthrep = 0.E+00
      jjha = (J1+J2-J3)/2 + 1
      jjhb = (J1-J2+J3)/2 + 1
      jjhc = (-J1+J2+J3)/2 + 1

      IF ( (jjha.LT.1) .OR. (jjhb.LT.1) .OR. (jjhc.LT.1) .OR. 
     &     ((M1+M2+M3).NE.0) ) THEN
         WTHREJ = wthrep
         GOTO 99999
      ENDIF

      jjhd = (J1+J2+J3+4)/2
      jmax = MAX(J1,J2,J3)
      IF ( jmax.NE.J1 ) THEN
         IF ( jmax.EQ.J2 ) THEN
            jj1 = J3
            jj2 = J1
            jj3 = J2
            mm1 = M3
            mm2 = M1
            mm3 = M2
            GOTO 100
         ELSEIF ( jmax.EQ.J3 ) THEN
            jj1 = J1
            jj2 = J2
            jj3 = J3
            mm1 = M1
            mm2 = M2
            mm3 = M3
            GOTO 100
         ENDIF
      ENDIF
      jj1 = J2
      jj2 = J3
      jj3 = J1
      mm1 = M2
      mm2 = M3
      mm3 = M1

 100  jma = (jj1+mm1)/2
      jmb = (jj1-mm1)/2
      jmc = (jj2+mm2)/2
      jmd = (jj2-mm2)/2
      jme = (jj3+mm3)/2
      jmf = (jj3-mm3)/2
      jabc = (jj1+jj2-jj3)/2
      jabm = (jj2-jj3-mm1)/2
      jbma = (jj1+mm2-jj3)/2
      izmin = MAX(jabm,jbma,0)
      izmax = MIN(jabc,jmb,jmc)
      nmax = MAX(jjhd,izmax+1)
      DO n = 1 , 26
         IF ( IP(n).GE.nmax ) GOTO 200
      ENDDO
      WTHREJ = wthrep
      GOTO 99999

 200  DO jlp = 1 , n
         jta = KF(jjha,jlp) + KF(jjhb,jlp) + KF(jjhc,jlp) - KF(jjhd,jlp)
         jtb = KF(jma+1,jlp) + KF(jmb+1,jlp) + KF(jmc+1,jlp)
         jtc = KF(jmd+1,jlp) + KF(jme+1,jlp) + KF(jmf+1,jlp)
         jvora(jlp) = jta + jtb + jtc
      ENDDO

      vorz = -1.E+00
      IF ( 2*(izmin/2).EQ.izmin ) vorz = +1.E+00
      IF ( izmin.LE.izmax ) THEN
         DO iz = izmin , izmax
            qsumlo = 0.E+00
            iza = iz + 1
            izb = jabc + 1 - iz
            izc = jmb + 1 - iz
            izd = jmc + 1 - iz
            ize = iz - jabm + 1
            izf = iz - jbma + 1
            DO jlp = 1 , n
               izexp = jvora(jlp) - 2*KF(iza,jlp) - 2*KF(izb,jlp)
     &                 - 2*KF(izc,jlp) - 2*KF(izd,jlp) - 2*KF(ize,jlp)
     &                 - 2*KF(izf,jlp)
               sumlo = izexp
               qsumlo = qsumlo + sumlo*PILOG(jlp)*(.5E+00)
            ENDDO
            zuthre = vorz*EXP(qsumlo)
            wthrep = wthrep + zuthre
            vorz = -vorz
         ENDDO
         jvo = jj1 - jj2 - mm3
         IF ( 4*(jvo/4).NE.jvo ) wthrep = -wthrep
      ENDIF
      WTHREJ = wthrep
99999 END
 
C----------------------------------------------------------------------
C FUNCTION WSIXJ
C
C Called by: F, GF, GKK, GOSIA
C
C Purpose: evaluates a Wigner 6-j symbol.
C
C Uses global variables:
C      IP     - table of prime numbers
C      KF     - sum of factors of primes
C      PILOG  - table of natural logs of primes
C
C Formal parameters:
C      J1     - twice the value of J1
C      J2     - twice the value of J2
C      J3     - twice the value of J3
C      L1     - twice the value of L1
C      L2     - twice the value of L2
C      L3     - twice the value of L3
C
C Return value:
C      The value of the 6-j symbol
C
C Note that the values of the parameters are doubled, so that this function
C can handle half-integers. In other words if you want to evaluate
C \sixj(J1 J2 J3 L1 L2 L3) you need to use call the function as:
C WSIXJ(2 * J1, 2 * J2, 2 * J3, 2 * L1, 2 * L2, 2 * L3).
 
      REAL*8 FUNCTION WSIXJ(J1,J2,J3,L1,L2,L3)
      IMPLICIT NONE
      INTEGER*4 irj , irl , isa , isb , isc , isumfa , iva , ivb , 
     &          ivc , ivd , ivorfa , iz , iza , izb , izc , izd , 
     &          ize , izf
      INTEGER*4 izg , izh , izmax , izmin , J1 , J2 , J3 , kqa , 
     &          kqb , kqc , kqd , kra , krb , krc , krd , ksa , ksb , 
     &          ksc , ksd
      INTEGER*4 kta , ktb , ktc , ktd , kua , kub , kuc , L1 , L2 , L3 ,
     &          n , nmax
      REAL*8 qsumfa , qsumlo , sumlo , vorz , wsixp , zusix
      DIMENSION isumfa(26) , ivorfa(26)
      REAL*8 PILOG
      INTEGER*4 IP , IPI , KF
      COMMON /FAKUL / IP(26) , IPI(26) , KF(101,26) , PILOG(26)
      
      wsixp = 0.E+00
      IF ( ((J1+J2-J3).GE.0) .AND. ((J1-J2+J3).GE.0) .AND. 
     &     ((-J1+J2+J3).GE.0) ) THEN
         IF ( ((J1+L2-L3).GE.0) .AND. ((J1-L2+L3).GE.0) .AND. 
     &        ((-J1+L2+L3).GE.0) ) THEN
            IF ( ((L1+J2-L3).GE.0) .AND. ((L1-J2+L3).GE.0) .AND. 
     &           ((-L1+J2+L3).GE.0) ) THEN
               IF ( ((L1+L2-J3).GE.0) .AND. ((L1-L2+J3).GE.0) .AND. 
     &              ((-L1+L2+J3).GE.0) ) THEN
                  kqa = (J1+J2-J3)/2
                  kqb = (J1-J2+J3)/2
                  kqc = (J2+J3-J1)/2
                  kqd = (J1+J2+J3)/2
                  kra = (J1+L2-L3)/2
                  krb = (J1-L2+L3)/2
                  krc = (L2+L3-J1)/2
                  krd = (J1+L2+L3)/2
                  ksa = (L1+J2-L3)/2
                  ksb = (L1-J2+L3)/2
                  ksc = (J2+L3-L1)/2
                  ksd = (L1+J2+L3)/2
                  kta = (L1+L2-J3)/2
                  ktb = (L1-L2+J3)/2
                  ktc = (L2+J3-L1)/2
                  ktd = (L1+L2+J3)/2
                  izmin = MAX(kqd,krd,ksd,ktd)
                  kua = kqa + kta + J3
                  kub = ksc + ktc + L1
                  kuc = krb + ktb + L2
                  izmax = MIN(kua,kub,kuc)
                  IF ( izmin.LE.izmax ) THEN
                     nmax = MAX(izmax+2,kqd+2,krd+2,ksd+2,ktd+2)
                     DO n = 1 , 26
                        IF ( IP(n).GE.nmax ) GOTO 5
                     ENDDO
                  ENDIF
                  GOTO 100
 5                vorz = -1.E+00
                  IF ( 2*(izmin/2).EQ.izmin ) vorz = +1.E+00
                  DO irl = 1 , n
                     iva = KF(kqa+1,irl) + KF(kqb+1,irl) + KF(kqc+1,irl)
     &                     - KF(kqd+2,irl)
                     ivb = KF(kra+1,irl) + KF(krb+1,irl) + KF(krc+1,irl)
     &                     - KF(krd+2,irl)
                     ivc = KF(ksa+1,irl) + KF(ksb+1,irl) + KF(ksc+1,irl)
     &                     - KF(ksd+2,irl)
                     ivd = KF(kta+1,irl) + KF(ktb+1,irl) + KF(ktc+1,irl)
     &                     - KF(ktd+2,irl)
                     ivorfa(irl) = iva + ivb + ivc + ivd
                  ENDDO
                  DO iz = izmin , izmax
                     sumlo = 0.E+00
                     iza = iz + 2
                     izb = iz - kqd + 1
                     izc = iz - krd + 1
                     izd = iz - ksd + 1
                     ize = iz - ktd + 1
                     izf = kua - iz + 1
                     izg = kub - iz + 1
                     izh = kuc - iz + 1
                     DO irj = 1 , n
                        isa = 2*KF(iza,irj) - 2*KF(izb,irj)
     &                        - 2*KF(izc,irj)
                        isb = -2*KF(izd,irj) - 2*KF(ize,irj)
     &                        - 2*KF(izf,irj)
                        isc = ivorfa(irj) - 2*KF(izg,irj)
     &                        - 2*KF(izh,irj)
                        isumfa(irj) = isa + isb + isc
                        qsumfa = isumfa(irj)
                        sumlo = sumlo + qsumfa*PILOG(irj)
                     ENDDO
                     qsumlo = (.5E+00)*sumlo
                     zusix = EXP(qsumlo)*vorz
                     wsixp = wsixp + zusix
                     vorz = -vorz
                  ENDDO
               ENDIF
            ENDIF
         ENDIF
      ENDIF
 100  WSIXJ = wsixp
      END
 
C----------------------------------------------------------------------
C SUBROUTINE LAGRAN
C
C Called by: CONV, EFFIX, GOSIA
C Calls:     FUNC, FUNC1
C
C Purpose: perform a Lagrangian interpolation
C
C Formal parameters:
C      X      - x-coordinate of input data
C      Y      - y-coordinate of input data
C      Ndata  - number of data points
C      Ipc    - index for storing results
C      Xx     - value for which to interpolate
C      Yy     - result of interpolation
C      Iscal  - mode: 1 = linear, 2 = exponential, 3 = square root
C      Irc    - weighting mode
C
C Note that the effect of FUNC and FUNC1 depends on Iscal:
C Iscal = 1   FUNC(y) = y        FUNC1(y) = y
C Iscal = 2   FUNC(y) = ln(y)    FUNC1(y) = exp(y)
C Iscal = 3   FUNC(y) = sqrt(y)  FUNC1(y) = y^2

      SUBROUTINE LAGRAN(X,Y,Ndata,Ipc,Xx,Yy,Iscal,Irc)
      IMPLICIT NONE
      REAL*8 arh , FUNC , FUNC1 , t , w , X , Xx , Y , y1 , Yy
      INTEGER*4 i , Ipc , Irc , Iscal , j , Ndata
      DIMENSION X(*) , Y(*) , w(101) , arh(101,101)
      SAVE arh
      
      IF ( Irc.EQ.2 ) THEN
      ELSEIF ( Irc.EQ.3 ) THEN
         DO i = 1 , Ndata
            t = 1.
            DO j = 1 , Ndata
               IF ( i.NE.j ) t = t*(Xx-X(j))/(X(i)-X(j))
            ENDDO
            arh(Ipc,i) = t
         ENDDO
         GOTO 100
      ELSEIF ( Irc.EQ.4 ) THEN
         GOTO 100
      ELSE
         DO i = 1 , Ndata
            w(i) = 1.
            DO j = 1 , Ndata
               IF ( i.NE.j ) w(i) = w(i)*(Xx-X(j))/(X(i)-X(j))
            ENDDO
         ENDDO
      ENDIF

      Yy = 0.
      DO j = 1 , Ndata
         y1 = Y(j)
         Yy = Yy + w(j)*FUNC(y1,Iscal)
      ENDDO
      Yy = FUNC1(Yy,Iscal)
      RETURN

 100  Yy = 0.
      DO j = 1 , Ndata
         y1 = Y(j)
         Yy = Yy + arh(Ipc,j)*FUNC(y1,Iscal)
      ENDDO
      Yy = FUNC1(Yy,Iscal)
      END
 
C----------------------------------------------------------------------
C FUNCTION FUNC
C
C Called by: LAGRAN
C
C Purpose: evaluates f(y) = y, f(y) = log_e(y) or f(y) = sqrt(y), depending
C on the flag I
C
C Formal parameters:
C      Y      - argument to evaluate
C      I      - mode: 1 = linear, 2 = log, 3 = square root
C
C Return value:
C      returns either y, log(y) or sqrt(y) depending on mode
C
C Note that this is the inverse of FUNC1
      
      REAL*8 FUNCTION FUNC(Y,I)
      IMPLICIT NONE
      INTEGER*4 I
      REAL*8 Y
      
      IF ( I.EQ.2 ) THEN
         IF ( Y.LT.1.E-12 ) Y = 1.E-12
         FUNC = LOG(Y)
         RETURN
      ELSEIF ( I.EQ.3 ) THEN
         FUNC = SQRT(Y)
         GOTO 99999
      ENDIF
      FUNC = Y
99999 END
 
C----------------------------------------------------------------------
C FUNCTION FUNC1
C
C Called by: LAGRAN
C
C Purpose: evaluates f(y) = y, f(y) = e^y or f(y) = y^2, depending on the
C flag I
C
C Formal parameters:
C      Y      - argument to evaluate
C      I      - mode: 1 = linear, 2 = exp, 3 = square
C
C Return value:
C      returns either y, exp(y) or y^2 depending on mode
C
C Note that this is the inverse of FUNC

      REAL*8 FUNCTION FUNC1(Y,I)
      IMPLICIT NONE
      INTEGER*4 I
      REAL*8 Y
      
      IF ( I.EQ.2 ) THEN
         FUNC1 = EXP(Y)
         RETURN
      ELSEIF ( I.EQ.3 ) THEN
         FUNC1 = Y*Y
         GOTO 99999
      ENDIF
      FUNC1 = Y
99999 END
 
C----------------------------------------------------------------------
C SUBROUTINE GKVAC
C
C Called from: DECAY
C Calls:       GKK
C
C Purpose: calculate the nuclear deorientation and store the results in the
C GKI array of common GVAC
C
C Uses global variables:
C      BETAR  - recoil beta
C      GKI    - G_k for a single level
C      IEXP   - experiment number
C      ITTE   - thick target experiment flag
C      SPIN   - spin of level
C      TAU    - lifetime in picoseconds
C      VACDP  - G_k for each level
C      XLAMB  - Lambda*       (this is G(3) in GOSIA)
C
C Formal parameters:
C      Il     - level index
 
      SUBROUTINE GKVAC(Il)
      IMPLICIT NONE
      REAL*8 beta , sp , time
      INTEGER*4 i , Il
      REAL*8 TAU
      INTEGER*4 KSEQ
      COMMON /LEV   / TAU(100) , KSEQ(1500,4)
      REAL*8 BETAR
      COMMON /BREC  / BETAR(50)
      REAL*8 G(7) , AVJI , GAMMA , XLAMB , TIMEC , GFAC , FIEL , POWER
      COMMON /GGG   / AVJI , GAMMA , XLAMB , TIMEC , GFAC , FIEL , POWER
      EQUIVALENCE(AVJI,G)
      REAL*8 XA , XA1 , EP , TLBDG , VINF
      INTEGER*4 NEXPT , IZ , IZ1
      COMMON /CX    / NEXPT , IZ , XA , IZ1(50) , XA1(50) , EP(50) , 
     &                TLBDG(50) , VINF(50)
      REAL*8 GKI , SUM
      COMMON /GVAC  / GKI(3) , SUM(3)
      REAL*8 EPS, EROOT, FIEX
      INTEGER*4 IEXP, IAXS
      COMMON /KIN   / EPS(50) , EROOT(50) , FIEX(50,2) , IEXP , IAXS(50)
      REAL*8 VACDP, QCEN, DQ, XNOR, AKS
      INTEGER*4 IBYP
      COMMON /VAC   / VACDP(3,100) , QCEN , DQ , XNOR , AKS(6,100) ,
     &                IBYP
      REAL*8 EN , SPIN , ACCUR , DIPOL , ZPOL , ACCA
      INTEGER*4 ISO
      COMMON /COEX  / EN(100) , SPIN(100) , ACCUR , DIPOL , ZPOL , 
     &                ACCA , ISO
      INTEGER*4 ITTE
      COMMON /THTAR / ITTE(50)
      IF ( ABS(XLAMB).GE.1.E-9 ) THEN
         IF ( ITTE(IEXP).EQ.0 ) THEN
            sp = SPIN(Il) ! Spin of level
            beta = BETAR(IEXP)
            time = TAU(Il) ! lifetime of level
            CALL GKK(IZ,beta,sp,time,Il)
            VACDP(1,Il) = GKI(1)
            VACDP(2,Il) = GKI(2)
            VACDP(3,Il) = GKI(3)
            GOTO 99999
         ENDIF
      ENDIF
      DO i = 1 , 3
         VACDP(i,Il) = 1.
      ENDDO
99999 END
 
C----------------------------------------------------------------------
C SUBROUTINE GKK
C
C Called by: GKVAC
C Calls:     ATS, WSIXJ, XSTATIC
C
C Purpose: calculate time-dependent deorientation coefficients
C
C Uses global variables:
C      AKS    - <\alpha_k> values
C      AVJI   - average J  (this is G(1) in GOSIA)
C      DQ     - width of gaussian distribution
C      FIEL   - K          (this is G(6) in GOSIA)
C      GAMMA  - Gamma      (this is G(2) in GOSIA)
C      GFAC   - g          (this is G(5) in GOSIA)
C      GKI    - G_k for a single level
C      IBYP   - flag to indicate whether we calculate <\alpha_k>
C      POWER  - x          (this is G(7) in GOSIA)
C      QCEN   - center of gaussian distribution
C      SUM    - sum over 6-j symbol squared
C      TIMEC  - Tau_C      (this is G(4) in GOSIA)
C      XLAMB  - Lambda*    (this is G(3) in GOSIA)
C      XNOR   - normalisation factor
C
C Formal parameters:
C      Iz     - Z of nucleus
C      Beta   - v/c
C      Spin   - spin of state
C      Time   - lifetime of state
C      Il     - index into AKS array
C
C We start by calling XSTATIC to calculate the static part. This calculates
C QCEN (the centre of the gaussian charge state distribution), DQ (the
C gaussian width of this distribution) and XNOR (the normalization parameter
C such that the sum over probabilities is one).
C
C We calculate:
C <a_k> = \sum_l p(J_1) \sum_F {(2 F + 1)^2 \over 2 J_1 + 1} *
C                              {\sixj{F F k I I J_1}}^2.
C
C We include a correction to take into account the effect of nuclear lifetimes
C which are comparable to the mean time between random reorientations \tau_c.
C
C Note that certain values have defaults:
C AVJI = 3, GAMMA = 0.02, XLAMB = 0.0345, TIMEC = 3.5, GFAC = Z/A,
C FIEL = 6E-6 and POWER = 0.6, which are set in GOSIA, where they are treated
C as an array called G in the order of the values in the GGG common block.
C However, the user may change them using the VAC suboption of the CONT option
C of OP,COUL or OP,GOSI.
C
C Note that WSIXJ requires all its parameters to be doubled, so it can handle
C half-integers properly.
C
C The function ATS is used to determine the ground-state spin for a given
C element.

      SUBROUTINE GKK(Iz,Beta,Spin,Time,Il)
      IMPLICIT NONE
      REAL*8 alp , ATS , Beta , ccf , down , dwc , f , hmean , rk , sm ,
     &       Spin , Time , up , upc , valmi , w2 , wrt , WSIXJ , wsp , 
     &       xji , xlam
      INTEGER*4 i , if2 , ifq , Il , imean , inq , irk2 , ispin2 , 
     &          ixji2 , Iz , j , k , k1 , k2 , l , m , ncoup , nz
      REAL*8 GKI , SUM
      COMMON /GVAC  / GKI(3) , SUM(3)
      REAL*8 VACDP, QCEN, DQ, XNOR, AKS
      INTEGER*4 IBYP
      COMMON /VAC   / VACDP(3,100) , QCEN , DQ , XNOR , AKS(6,100) ,
     &                IBYP
      REAL*8 G(7) , AVJI , GAMMA , XLAMB , TIMEC , GFAC , FIEL , POWER
      COMMON /GGG   / AVJI , GAMMA , XLAMB , TIMEC , GFAC , FIEL , POWER
      EQUIVALENCE(AVJI,G)
      IF ( IBYP.NE.1 ) THEN
         imean = 0
         CALL XSTATIC(Iz,inq,ifq,Beta) ! inq and ifq are range of integral
         l = 0
         DO i = 1 , 6
            AKS(i,Il) = 0.
         ENDDO
 50      IF ( imean.EQ.1 ) inq = 1
         IF ( imean.EQ.1 ) ifq = 1

         DO j = inq , ifq
            l = l + 1
            nz = Iz - j
            xji = ATS(nz) ! Ground-state spin of atom
            sm = Spin
            IF ( imean.EQ.1 ) xji = AVJI
            IF ( Spin.GT.xji ) sm = xji
            ncoup = INT(2.*sm+.5) + 1
            SUM(1) = 0.
            SUM(2) = 0.
            SUM(3) = 0.
            valmi = Spin - xji
            IF ( valmi.LT.0. ) valmi = -valmi
            DO m = 1 , ncoup
               f = valmi + DBLE(m) - 1.
               DO k = 1 , 3
                  rk = 2.*DBLE(k)
                  if2 = f*2. + 0.0001
                  irk2 = rk*2. + 0.0001
                  ispin2 = Spin*2. + 0.0001
                  ixji2 = xji*2. + 0.0001
                  SUM(k) = SUM(k)
     &                     + ((2.*f+1.)*WSIXJ(if2,if2,irk2,ispin2,
     &                     ispin2,ixji2))**2/(2.*xji+1.)
               ENDDO
            ENDDO
            IF ( imean.NE.1 ) THEN
               DO k = 1 , 3
                  k1 = 2*k - 1
                  AKS(k1,Il) = AKS(k1,Il) + SUM(k)
     &                         *EXP(-((QCEN-DBLE(j))/DQ)**2/2.)/XNOR
               ENDDO
               IF ( imean.EQ.0 ) GOTO 100
            ENDIF
            DO k = 1 , 3
               k1 = 2*k
               AKS(k1,Il) = AKS(k1,Il) + SUM(k)
            ENDDO
 100        CONTINUE
         ENDDO ! Loop on j
         imean = imean + 1
         IF ( imean.EQ.1 ) GOTO 50
      ENDIF

      hmean = FIEL*Iz*(Beta**POWER) ! Mean magnetic field in fluctuating state
      wsp = 4789.*GFAC*hmean/AVJI ! 4789 is the nuclear magneton
      wsp = wsp*TIMEC
      wsp = wsp*wsp*AVJI*(AVJI+1.)/3.
      DO k = 1 , 3
         k2 = 2*k
         k1 = 2*k - 1
         wrt = wsp*k2*(k2+1)
         w2 = wrt
         wrt = -wrt/(1.-AKS(k2,Il))
         xlam = (1.-AKS(k2,Il))*(1.-EXP(wrt))/TIMEC
         up = (GAMMA*Time*AKS(k1,Il)+1.)/(Time*GAMMA+1.)
         up = up*XLAMB*Time + 1.       ! numerator
         down = Time*(xlam+XLAMB) + 1. ! denominator = r
         GKI(k) = up/down
         alp = 9.*xlam*xlam + 8.*xlam*TIMEC*(w2-xlam*xlam)
         alp = SQRT(alp) - 3.*xlam
         alp = alp/4./xlam/TIMEC                      ! alp is p
         upc = xlam*Time*(down-2.*alp*alp*Time*TIMEC) ! numerator
         dwc = (down+alp*Time)*(down+2.*alp*Time)     ! denominator
         ccf = 1. + upc/dwc                           ! ccf is correction factor
         GKI(k) = GKI(k)*ccf
      ENDDO
      END
 
C----------------------------------------------------------------------
C SUBROUTINE XSTATIC
C
C Called by: GKK
C
C Purpose: calculate the static part of the deorientation coefficients G. It
C is assumed to be a gaussian distribution.
C
C Uses global variables:
C      DQ     - width of gaussian distribution
C      QCEN   - center of gaussian distribution
C      XNOR   - normalisation factor
C
C Formal parameters:
C      Iz     - Z of nucleus
C      Ido    - lower limit for integral over gaussian to 3 sigma
C      Iup    - upper limit for integral over gaussian to 3 sigma
C      Beta   - beta
C
C We use: h = {1 \over {1 + (0.012008 \beta Z^0.45)^{5/3}}} and
C Q_0 = Z h^6 (QCEN here)
C \sigma_Q = \sqrt(Q_0 (1 - h)) (DQ here)
C We also calculate the normalization factor, needed to ensure that the sum
C is unity.
C The value 0.012008 is v'/c, where v' is taken from Nikolaev and Dmitriev,
C Phys. Lett. 82A, 277, to be 3.6 x 10^6 m/s. The 0.45 is the coefficient
C alpha from the same paper. The power of 5/3 is 1/k = 1/0.6 from the same
C paper.

      SUBROUTINE XSTATIC(Iz,Ido,Iup,Beta)
      IMPLICIT NONE
      REAL*8 Beta , h
      INTEGER*4 Ido , Iup , Iz , lq
      REAL*8 VACDP, QCEN, DQ, XNOR, AKS
      INTEGER*4 IBYP
      COMMON /VAC   / VACDP(3,100) , QCEN , DQ , XNOR , AKS(6,100) ,
     &                IBYP
      
      h = 1./(1.+(Iz**.45*.012008/Beta)**1.666667)
      QCEN = Iz*h**.6
      DQ = SQRT(QCEN*(1.-h))/2.
      
      Iup = INT(QCEN+3.*DQ+.5)
      Ido = INT(QCEN-3.*DQ-.5)
      IF ( Iup.GT.Iz ) Iup = Iz
      IF ( Ido.LT.1 ) Ido = 1
      XNOR = 0.
      DO lq = Ido , Iup
         XNOR = XNOR + EXP(-((QCEN-DBLE(lq))/DQ)**2/2.)
      ENDDO
      END
 
C----------------------------------------------------------------------
C FUNCTION ATS
C
C Called by: GKK
C
C Purpose: determine the atomic ground-state spin
C
C Formal parameters:
C      N      - number of electrons (Z - charge state)
C
C Return value:
C      truncation point
 
      REAL*8 FUNCTION ATS(N)
      IMPLICIT NONE
      INTEGER*4 m , N
      REAL*8 x , xm
      
      IF ( N.LE.0 .OR. N.GT.96 ) THEN
         ATS = 0.
         RETURN
      ELSE
         x = N/2. + 1
         m = N/2 + 1
         xm = DBLE(m)
         IF ( ABS(x-xm).GE.1.E-9 ) THEN
            IF ( m.EQ.1 .OR. m.EQ.2 .OR. m.EQ.3 .OR. m.EQ.6 .OR. 
     &           m.EQ.7 .OR. m.EQ.10 .OR. m.EQ.15 .OR. m.EQ.16 .OR. 
     &           m.EQ.19 .OR. m.EQ.24 .OR. m.EQ.25 .OR. m.EQ.28 .OR. 
     &           m.EQ.31 .OR. m.EQ.35 .OR. m.EQ.37 .OR. m.EQ.40 .OR. 
     &           m.EQ.41 .OR. m.EQ.44 ) THEN
               ATS = .5
               RETURN
            ELSEIF ( m.EQ.4 .OR. m.EQ.5 .OR. m.EQ.8 .OR. m.EQ.9 .OR. 
     &               m.EQ.11 .OR. m.EQ.17 .OR. m.EQ.18 .OR. m.EQ.20 .OR.
     &               m.EQ.26 .OR. m.EQ.27 .OR. m.EQ.36 .OR. m.EQ.42 .OR.
     &               m.EQ.43 .OR. m.EQ.45 ) THEN
               ATS = 1.5
               RETURN
            ELSEIF ( m.EQ.12 .OR. m.EQ.14 .OR. m.EQ.21 .OR. m.EQ.23 .OR.
     &               m.EQ.32 .OR. m.EQ.39 ) THEN
               ATS = 2.5
               RETURN
            ELSEIF ( m.EQ.13 .OR. m.EQ.22 .OR. m.EQ.38 ) THEN
               ATS = 4.5
               RETURN
            ELSEIF ( m.EQ.29 .OR. m.EQ.30 .OR. m.EQ.48 ) THEN
               ATS = 3.5
               RETURN
            ELSEIF ( m.EQ.33 ) THEN
               ATS = 7.5
               RETURN
            ELSEIF ( m.EQ.34 ) THEN
               ATS = 6.5
               GOTO 99999
            ELSEIF ( m.EQ.46 .OR. m.EQ.47 ) THEN
               ATS = 5.5
               RETURN
            ENDIF
         ENDIF
         m = m - 1
         IF ( m.EQ.4 .OR. m.EQ.8 .OR. m.EQ.17 .OR. m.EQ.26 .OR. 
     &        m.EQ.28 .OR. m.EQ.30 .OR. m.EQ.32 .OR. m.EQ.42 .OR. 
     &        m.EQ.45 .OR. m.EQ.48 ) THEN
            ATS = 2.
            RETURN
         ELSEIF ( m.EQ.10 .OR. m.EQ.36 ) THEN
         ELSEIF ( m.EQ.12 .OR. m.EQ.21 .OR. m.EQ.37 ) THEN
            ATS = 3.
            RETURN
         ELSEIF ( m.EQ.13 .OR. m.EQ.22 .OR. m.EQ.29 .OR. m.EQ.31 .OR. 
     &            m.EQ.34 .OR. m.EQ.38 .OR. m.EQ.47 ) THEN
            ATS = 4.
            RETURN
         ELSEIF ( m.EQ.33 ) THEN
            ATS = 8.
            RETURN
         ELSEIF ( m.EQ.46 ) THEN
            ATS = 6.
            RETURN
         ELSE
            ATS = 0.
            RETURN
         ENDIF
      ENDIF
      ATS = 1.
99999 END
 
C----------------------------------------------------------------------
C SUBROUTINE YLM
C
C Called by: ANGULA
C
C Purpose: evaluate the even spherical harmonics.
C
C Uses global variables:
C      IEXP   - experiment number
C      IAXS   - axial symmetry flag
C
C Formal parameters:
C      Theta  - theta for which to evaluate (read only)
C      Ylmr   - return value for that theta (write only)
C
C Ylmr(l,m) = 1 / \sqrt{4 \pi} Y_{2l}^{m - 1}
C
C Note the factor of 1 / \sqrt{4 \pi} compared to the orthonormal spherical
C harmonics.
C
C 0.0889703179  = sqrt(5) / (8 pi)
C 0.0298415518  = 3 / (32 pi)
C 0.0179325408  = sqrt(13) / (64 pi)
C 0.1089659406  = sqrt(30) / (16 pi)
C -0.2179318812 = -1 * sqrt(30) / (8 pi)
C 0.1248361677  = 3 * sqrt(70) / (64 pi)
C -0.3530900028 = -3 * sqrt(140) / (32 pi)
C 0.0943672726  = 3 * sqrt(10) / (32 pi)
C -0.1334554768 = -3 * sqrt(20) / (32 pi)
C etc.
      
      SUBROUTINE YLM(Theta,Ylmr)
      IMPLICIT NONE
      REAL*8 ct , ctsq , st , Theta , Ylmr
      INTEGER*4 i , j , l , lf , m
      REAL*8 EPS, EROOT, FIEX
      INTEGER*4 IEXP, IAXS
      COMMON /KIN   / EPS(50) , EROOT(50) , FIEX(50,2) , IEXP , IAXS(50)
      DIMENSION Ylmr(9,9) , st(7)
      
      ct = COS(Theta)
      ctsq = ct*ct
      IF ( IAXS(IEXP).EQ.0 ) THEN ! If axially symmetric
         Ylmr(1,1) = .0889703179*(3.*ctsq-1.)
         Ylmr(2,1) = .0298415518*((35.*ctsq-30.)*ctsq+3.)
         Ylmr(3,1) = .0179325408*(((231.*ctsq-315.)*ctsq+105.)*ctsq-5.)
         GOTO 99999
      ENDIF
      st(1) = SIN(Theta)
      DO i = 2 , 7
         j = i - 1
         st(i) = st(j)*st(1)
      ENDDO
      Ylmr(1,3) = .1089659406
      Ylmr(1,2) = -.2179318812*ct
      Ylmr(1,1) = .0889703179*(3.*ctsq-1.)
      Ylmr(2,5) = .1248361677
      Ylmr(2,4) = -.3530900028*ct
      Ylmr(2,3) = .0943672726*(7.*ctsq-1.)
      Ylmr(2,2) = -.1334554768*ct*(7.*ctsq-3.)
      Ylmr(2,1) = .0298415518*((35.*ctsq-30.)*ctsq+3.)
      Ylmr(3,7) = .1362755124
      Ylmr(3,6) = -.4720722226*ct
      Ylmr(3,5) = .100646136*(11.*ctsq-1.)
      Ylmr(3,4) = -.1837538634*ct*(11.*ctsq-3.)
      Ylmr(3,3) = .0918769316*((33.*ctsq-18.)*ctsq+1.)
      Ylmr(3,2) = -.1162161475*ct*((33.*ctsq-30.)*ctsq+5.)
      Ylmr(3,1) = .0179325408*(((231.*ctsq-315.)*ctsq+105.)*ctsq-5.)
      DO l = 1 , 3
         lf = 2*l + 1
         DO m = 2 , lf
            Ylmr(l,m) = Ylmr(l,m)*st(m-1)
         ENDDO
      ENDDO
99999 END
 
C----------------------------------------------------------------------
C SUBROUTINE DECAY
C
C Called by: CEGRY, GOSIA
C Calls:     GKVAC
C
C Purpose: Calculate the gamma decay following excitation.
C
C Uses global variables:
C      DELLA  - products of matrix elements: e1^2, e2^2, e1*e2
C      DELTA  - \delta_\lambda: index 1 = electric^2, 2 = magnetic^2, 3 = cross term
C      GKP    - Gk * DELTA^2
C      IAXS   - axial symmetry flag
C      IBYP   - flag to indicate whether we calculate <\alpha_k>
C      IEXP   - experiment number
C      KLEC   - number of decays for each level
C      KSEQ   - index into ELM for pair of levels, and into EN or SPIN
C      LIFCT  - index of level for lifetimes
C      NMAX   - number of levels
C      NMAX1  - number of levels with decays
C      TAU    - lifetime in picoseconds
C      TIMEL  - lifetimes and their errors
C      VACDP  - G_k for each level
C      ZETA   - various coefficients
C
C Formal parameters:
C      Chisq  - chi squared
C      Nlift  - number of lifetimes
C      Chilo  - chi squared of logs
      
      SUBROUTINE DECAY(Chisq,Nlift,Chilo)
      IMPLICIT NONE
      REAL*8 bsum , Chilo , Chisq , df , el1 , emt , emt1 , gk , vcd
      INTEGER*4 i , ibra , idr , idrh , ifn , il , inx , inx1 , iu , 
     &          j , jlt , k , kl , kq , l , l1 , lc1 , lc2 , n1 , n2
      INTEGER*4 Nlift
      REAL*8 DELTA, ENDEC, ENZ
      INTEGER*4 ITMA
      COMMON /TRA   / DELTA(1500,3) , ENDEC(1500) , ITMA(50,200) , 
     &                ENZ(200)
      REAL*8 TIMEL
      INTEGER*4 LIFCT
      COMMON /LIFE1 / LIFCT(50) , TIMEL(2,50)
      REAL*8 VACDP, QCEN, DQ, XNOR, AKS
      INTEGER*4 IBYP
      COMMON /VAC   / VACDP(3,100) , QCEN , DQ , XNOR , AKS(6,100) ,
     &                IBYP
      REAL*8 ZETA
      INTEGER*4 LZETA
      COMMON /CCOUP / ZETA(155600) , LZETA(8)
      REAL*8 TAU
      INTEGER*4 KSEQ
      COMMON /LEV   / TAU(100) , KSEQ(1500,4)
      INTEGER*4 NMAX , NDIM , NMAX1
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      REAL*8 ELM , ELMU , ELML , SA
      COMMON /COMME / ELM(1500) , ELMU(1500) , ELML(1500) , SA(1500)
      REAL*8 EPS, EROOT, FIEX
      INTEGER*4 IEXP, IAXS
      COMMON /KIN   / EPS(50) , EROOT(50) , FIEX(50,2) , IEXP , IAXS(50)
      REAL*8 FP , GKP
      INTEGER*4 KLEC
      COMMON /CATLF / FP(4,1500,3) , GKP(4,1500,2) , KLEC(100)
      REAL*8 DELLA
      COMMON /LCDL  / DELLA(1500,3)
      DIMENSION gk(4)
      DATA emt1/0./

      idr = 1
      DO il = 1 , NMAX1 ! For each level with decays
         l = KSEQ(idr,3) ! Initial level of idr'th decay
         n1 = 28*(l-1)
         ibra = KLEC(l) ! Number of decays from level l
         bsum = 0.
         idrh = idr
         DO j = 1 , ibra ! For each decay from level l
            inx = KSEQ(idr,1) ! Index to matrix element of idr'th decay
            inx1 = KSEQ(idr,2) ! Index 2 of idr'th decay
            el1 = 0.
            IF ( inx.NE.0 ) el1 = ELM(inx)
            emt = el1*el1
            DELLA(idr,1) = emt
            IF ( inx1.NE.0 ) emt1 = ELM(inx1)*ELM(inx1)
            bsum = bsum + DELTA(idr,1)*emt
            IF ( inx1.NE.0 ) THEN
               DELLA(idr,3) = el1*ELM(inx1)
               DELLA(idr,2) = emt1
               bsum = bsum + DELTA(idr,2)*emt1
            ENDIF
            idr = idr + 1
         ENDDO ! Loop on j

         idr = idrh
         TAU(l) = 1./bsum
         CALL GKVAC(l) ! Evaluate G_k

         DO j = 1 , ibra ! For each decay from level l
            l1 = KSEQ(idr,4) ! Final energy of idr'th decay
            n2 = 28*(l1-1)
            inx1 = KSEQ(idr,2) ! Index 2 of idr'th decay
            DO i = 1 , 4
               gk(i) = GKP(i,idr,1)*DELLA(idr,1)
            ENDDO
            IF ( inx1.NE.0 ) THEN
               DO i = 1 , 4
                  gk(i) = gk(i) + GKP(i,idr,2)*DELLA(idr,2)
               ENDDO
            ENDIF
            DO i = 1 , 4
               vcd = 1.
               IF ( i.NE.1 ) vcd = VACDP(i-1,l)
               gk(i) = gk(i)*TAU(l)
               ifn = 2*i - 1
               iu = (i-1)*7
               IF ( IAXS(IEXP).EQ.0 ) ifn = 1
               DO kq = 1 , ifn
                  lc1 = n1 + iu + kq
                  lc2 = n2 + iu + kq
                  ZETA(lc2) = ZETA(lc2) + gk(i)*vcd*ZETA(lc1)
               ENDDO
            ENDDO
            idr = idr + 1
         ENDDO ! Loop on j
      ENDDO ! Loop on l

      IBYP = 1 ! Set flag to indicate we have calculated <\alpha_k>
      IF ( Nlift.NE.0 .AND. IEXP.EQ.1 ) THEN
         DO jlt = 1 , Nlift ! For each lifetime
            kl = LIFCT(jlt) ! Get level for this lifetime
            df = (TAU(kl)-TIMEL(1,jlt))/TIMEL(2,jlt) ! TIMEL(1,X) is lifetime and TIMEL(2,X) is the error
            Chilo = Chilo + (LOG(TAU(kl)/TIMEL(1,jlt))*TIMEL(1,jlt)
     &              /TIMEL(2,jlt))**2 ! Log chisqr
            Chisq = Chisq + df*df ! Chisqr
         ENDDO
      ENDIF

      DO l = 2 , NMAX ! For each level except the ground state
         IF ( KLEC(l).NE.0 ) THEN ! If there are decays from this level
            n1 = 28*(l-1)
            DO j = 1 , 4
               vcd = 1.
               IF ( j.NE.1 ) vcd = VACDP(j-1,l) ! G_k for each level
               ifn = 2*j - 1
               iu = (j-1)*7
               DO k = 1 , ifn
                  lc1 = n1 + iu + k
                  ZETA(lc1) = ZETA(lc1)*vcd
               ENDDO
            ENDDO
         ENDIF
      ENDDO ! Loop on levels
      END
 
C----------------------------------------------------------------------
C SUBROUTINE ANGULA
C
C Called by: GOSIA, CEGRY
C Calls:     FIINT, FIINT1, RECOIL, YLM, YLM1
C
C Purpose: calculate angular distribution of emitted gamma rays
C
C Uses global variables:
C      BETAR  - recoil beta
C      DELLA  - products of matrix elements: e1^2, e2^2, e1*e2
C      ENDEC  - energy difference for each matrix element
C      ENZ    - something to do with the absorption
C      FP     - F coefficient * DELTA^2
C      IAXS   - axial symmetry flag
C      IEXP   - experiment number
C      ITMA   - identify detectors according to OP,GDET
C      ITTE   - thick target experiment flag
C      KSEQ   - index into ELM for pair of levels, and into EN or SPIN
C      Q      - solid angle attenuation coefficients
C      TAU    - lifetime in picoseconds
C      ZETA   - various coefficients
C
C Formal parameters:
C      Ygn    - Gamma-ray yield
C      Idr    - number of decays
C      Iful   - flag to select full basis or not
C      Fi0    - phi_0
C      Fi1    - phi_1
C      Trec   - Theta of recoiling nucleus
C      Gth    - Theta of gamma
C      Figl   - Phi of gamma
C      Ngl    - detector number
C      Op2    - The part after the OP, for the option we are processing
      
      SUBROUTINE ANGULA(Ygn,Idr,Iful,Fi0,Fi1,Trec,Gth,Figl,Ngl,Op2)
      IMPLICIT NONE
      REAL*8 alab , arg , at , attl , bt , f , Fi0 , fi01 , Fi1 ,
     &       fi11 , Figl , Gth , qv , sm , Trec , trec2 , Ygn , ylmr
      INTEGER*4 Idr , ifn , Iful , ig , il , inat , inx1 , 
     &          ipd , is , iu , ixs , j , ji , jj , jm , k
      INTEGER*4 kq , l , lf , lf1 , mind , Ngl , nlv
      CHARACTER*4 Op2
      DIMENSION f(4) , ylmr(9,9) , at(28) , alab(9,9) , attl(9,9) , 
     &          Ygn(*)
      REAL*8 ZETA
      INTEGER*4 LZETA
      COMMON /CCOUP / ZETA(155600) , LZETA(8)
      REAL*8 DELTA, ENDEC, ENZ
      INTEGER*4 ITMA
      COMMON /TRA   / DELTA(1500,3) , ENDEC(1500) , ITMA(50,200) , 
     &                ENZ(200)
      REAL*8 TAU
      INTEGER*4 KSEQ
      COMMON /LEV   / TAU(100) , KSEQ(1500,4)
      REAL*8 EG, CC, AGELI, Q
      INTEGER*4 NICC, NANG, ISPL
      COMMON /CCC   / EG(50) , CC(50,5) , AGELI(50,200,2) , Q(3,200,8) ,
     &                NICC , NANG(200) , ISPL
      REAL*8 EPS, EROOT, FIEX
      INTEGER*4 IEXP, IAXS
      COMMON /KIN   / EPS(50) , EROOT(50) , FIEX(50,2) , IEXP , IAXS(50)
      REAL*8 DELLA
      COMMON /LCDL  / DELLA(1500,3)
      REAL*8 FP , GKP
      INTEGER*4 KLEC
      COMMON /CATLF / FP(4,1500,3) , GKP(4,1500,2) , KLEC(100)
      REAL*8 BETAR
      COMMON /BREC  / BETAR(50)
      INTEGER*4 ITTE
      COMMON /THTAR / ITTE(50)
      REAL*8 XA , XA1 , EP , TLBDG , VINF
      INTEGER*4 NEXPT , IZ , IZ1
      COMMON /CX    / NEXPT , IZ , XA , IZ1(50) , XA1(50) , EP(50) , 
     &                TLBDG(50) , VINF(50)
      
      DO l = 1 , Idr ! For each decay

         nlv = KSEQ(l,3) ! Level number of l'th decay
         il = (nlv-1)*28
         inx1 = KSEQ(l,2) ! Index of l'th decay

         DO j = 1 , 4
            f(j) = FP(j,l,1)*DELLA(l,1)
         ENDDO

         IF ( inx1.NE.0 ) THEN
            DO j = 1 , 4
               f(j) = f(j) + 2.*FP(j,l,3)*DELLA(l,3) + FP(j,l,2)
     &                *DELLA(l,2)
            ENDDO
         ENDIF

         DO j = 1 , 4
            f(j) = f(j)*TAU(nlv)
            iu = (j-1)*7
            ifn = 2*j - 1
            IF ( IAXS(IEXP).EQ.0 ) ifn = 1
            DO kq = 1 , ifn
               is = iu + kq
               ig = is + il
               at(is) = ZETA(ig)*f(j)
            ENDDO
         ENDDO

         IF ( Iful.EQ.1 ) THEN
            DO j = 1 , 9
               DO k = 1 , 9
                  alab(j,k) = 0.
                  attl(j,k) = 0.
               ENDDO
            ENDDO
            DO j = 1 , 4
               lf = 2*j - 1
               lf1 = lf
               IF ( IAXS(IEXP).EQ.0 ) lf1 = 1
               DO k = 1 , lf1
                  inat = (j-1)*7 + k
                  alab(lf,k) = at(inat)
               ENDDO
            ENDDO
            bt = BETAR(IEXP) ! Get beta
            trec2 = SIGN(Trec, DBLE(IZ1(IEXP)))
            IF ( ITTE(IEXP).NE.1 ) CALL RECOIL(alab,attl,bt,trec2) ! Relativistic correction
            IF ( l.EQ.1 ) CALL YLM1(Gth,ylmr)
            ixs = IAXS(IEXP) ! Get axial symmetry flag
            fi01 = Fi0 - Figl ! Get lower phi limit
            fi11 = Fi1 - Figl ! Get upper phi limit
            CALL FIINT1(fi01,fi11,alab,ixs) ! Integrate over phi in lab frame
            Ygn(l) = alab(1,1)*.0795774715 ! 0.0795774715 = 1 / (4 pi)
            DO j = 2 , 9
               sm = ylmr(j,1)*alab(j,1)
               IF ( IAXS(IEXP).NE.0 ) THEN
                  DO k = 2 , j
                     sm = sm + 2.*ylmr(j,k)*alab(j,k)
                  ENDDO
               ENDIF
               ipd = ITMA(IEXP,Ngl) ! Detector ID
               arg = (ENDEC(l)-ENZ(ipd))**2
               qv = (Q(3,ipd,j-1)*Q(2,ipd,j-1)+Q(1,ipd,j-1)*arg)
     &              /(Q(2,ipd,j-1)+arg)
               Ygn(l) = Ygn(l) + sm*qv
            ENDDO
         ELSE
            ixs = IAXS(IEXP) ! Get axial symmetry flag
            fi01 = Fi0 - Figl ! Get lower phi limit
            fi11 = Fi1 - Figl ! Get upper phi limit
            CALL FIINT(fi01,fi11,at,ixs) ! Integrate over phi in recoiling nucleus frame, result in at
            IF ( l.EQ.1 ) CALL YLM(Gth,ylmr)
            Ygn(l) = at(1)*.0795774715 ! 0.0795774715 = 1 / (4 pi)
            DO jj = 1 , 3
               ji = jj*7 + 1
               sm = ylmr(jj,1)*at(ji)
               IF ( IAXS(IEXP).NE.0 ) THEN
                  mind = 2*jj + 1
                  DO jm = 2 , mind
                     ji = ji + 1
                     sm = ylmr(jj,jm)*at(ji)*2. + sm
                  ENDDO
               ENDIF
               ipd = ITMA(IEXP,Ngl) ! Detector ID
               arg = (ENDEC(l)-ENZ(ipd))**2
               qv = (Q(3,ipd,2*jj)*Q(2,ipd,2*jj)+Q(1,ipd,2*jj)*arg)
     &              /(Q(2,ipd,2*jj)+arg) ! solid angle attenuation coefficients
               Ygn(l) = Ygn(l) + sm*qv
            ENDDO
         ENDIF
      ENDDO ! Loop over decays

      IF ( Op2.EQ.'INTG' .OR. Op2.EQ.'INTI' ) RETURN

C     In gosia2, we multiply by dsig*SIN(ttx) here, but not in gosia
      END
 
C----------------------------------------------------------------------
C SUBROUTINE READY
C
C Called by: ADHOC
C Calls:     SZEREG
C
C Purpose: To read experimental yields from a file.
C
C Uses global variables:
C      DYEX   - error on experimental yield
C      IY     - index of experimental yields
C      KSEQ   - index into ELM for pair of levels, and into EN or SPIN
C      LP6    - 32
C      NDST   - number of data sets
C      NEXPT  - number of experiments
C      NYLDE  - number of yields
C      YEXP   - experimental yield
C
C Formal parameters:
C      Idr    - number of decays
C      Ntap   - unit for yield file
C      Ipri   - printing flag (Ipri=1 gives additional output)
C
C NTAP is the unit number of the file from which we should read the
C experimental yields
 
      SUBROUTINE READY(Idr,Ntap,Ipri)
      IMPLICIT NONE
      REAL*8 ap , u , w , waga , xep , zp
      INTEGER*4 idc , idc1 , idcx , Idr , ii , Ipri , 
     &          iytot , iytt , j , k , kk , kkl , lbg
      INTEGER*4 lxp , nanx , nde , nde1 , ne , ns1
      INTEGER*4 ns2 , ns3 , ns4 , nsxh , nsyh , Ntap , nval
      DIMENSION iytot(32)
      REAL*8 YEXP, CORF , DYEX , UPL , YNRM
      INTEGER*4 IY , NYLDE , IDRN , ILE
      COMMON /YEXPT / YEXP(32,1500) , IY(1500,32) , CORF(1500,32) , 
     &                DYEX(32,1500) , NYLDE(50,32) , UPL(32,50) , 
     &                YNRM(32,50) , IDRN , ILE(32)
      REAL*8 XA , XA1 , EP , TLBDG , VINF
      INTEGER*4 NEXPT , IZ , IZ1
      COMMON /CX    / NEXPT , IZ , XA , IZ1(50) , XA1(50) , EP(50) , 
     &                TLBDG(50) , VINF(50)
      REAL*8 TAU
      INTEGER*4 KSEQ
      COMMON /LEV   / TAU(100) , KSEQ(1500,4)
      INTEGER*4 LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &          LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      INTEGER*4 NDST
      COMMON /CCCDS / NDST(50)
      
C     Rewind yield file
      REWIND Ntap

      DO k = 1 , LP6 ! LP6 = 32
         iytot(k) = 0
      ENDDO
      IF ( Ipri.EQ.1 ) WRITE (22,99001)
99001 FORMAT (5X/47X,'REPRINT OF EXPERIMENTAL DATA TO BE FITTED'//)
      DO lxp = 1 , NEXPT ! For each experiment
         DO kkl = 1 , LP6
            NYLDE(lxp,kkl) = 0
         ENDDO
         ii = NDST(lxp) ! Number of datasets
         DO kk = 1 , ii ! iexp, ng, zp, ag, ep, nd, wt
            READ (Ntap,*) ne , nanx , zp , ap , xep , nval , waga
            IF ( Ipri.EQ.1 ) WRITE (22,99002) ne , zp , ap , xep , 
     &                              NDST(ne) , waga
99002       FORMAT (1X,///5X,'EXPERIMENT',1X,1I2/2X,'PROJECTILE',1X,'(',
     &              1F4.0,',',1F4.0,')',1X,1F7.3,1X,'MEV',1X,'---',1I1,
     &              1X,'GE(LI) DETECTOR(S)',2X,'WEIGHT=',1E8.2/20X,
     &              '** EXPERIME','NTAL YIELDS **')
            IF ( Ipri.EQ.1 ) WRITE (22,99003)
99003       FORMAT (4X,'DECAY',1X,'IS',2X,'IF',1(9X,'YIELD+/-ERROR',9X)
     &              /)
            DO j = 1 , nval
               READ (Ntap,*) ns1 , ns2 , u , w ! ii, if, y, dy
               nsxh = ns1
               nsyh = ns2
               IF ( ns1.GE.100 ) THEN
                  ns1 = ns1/100
                  ns2 = ns2/100
               ENDIF
               DO nde = 1 , Idr
                  IF ( ns1.EQ.KSEQ(nde,3) .AND. ns2.EQ.KSEQ(nde,4) )
     &                 GOTO 10
               ENDDO
               IF ( Ipri.EQ.1 ) WRITE (22,99005) ns1 , ns2
               GOTO 40
 10            idc = nde
               iytot(kk) = iytot(kk) + 1
               idc1 = 0
               IF ( nsxh.GE.100 ) THEN
                  ns3 = nsxh - 100*ns1
                  ns4 = nsyh - 100*ns2
                  DO nde1 = 1 , Idr
                     IF ( ns3.EQ.KSEQ(nde1,3) .AND. ns4.EQ.KSEQ(nde1,4)
     &                    ) GOTO 20
                  ENDDO
                  IF ( Ipri.EQ.1 ) WRITE (22,99005) ns3 , ns4
               ENDIF
               GOTO 30
 20            idcx = idc*1000 + nde1
               IF ( idc.GT.nde1 ) idcx = nde1*1000 + idc
               idc = idcx
 30            idc1 = idc
               IF ( idc1.GT.1000 ) idc1 = idc/1000
               IF ( Ipri.EQ.1 ) WRITE (22,99004) idc , KSEQ(idc1,3) , 
     &                                 KSEQ(idc1,4) , u , w
99004          FORMAT (2X,1I6,1X,1I3,1X,1I3,1(1E14.6,3X,1E14.6))
               iytt = iytot(kk)
               YEXP(kk,iytt) = u
               DYEX(kk,iytt) = w/(SQRT(waga)+1.E-4)
               IY(iytt,kk) = idc
 40            CONTINUE
            ENDDO
            iytt = iytot(kk)
            lbg = iytt - nval + 1
            CALL SZEREG(lbg,iytt,kk)
            NYLDE(lxp,kk) = nval
         ENDDO ! For each dataset kk
      ENDDO ! Loop on experiments lxp
99005 FORMAT (1X///5X,'ERROR-NO MATRIX ELEMENT BETWEEN STATES',1X,1I3,
     &        ' AND ',1I3,/10X,'THIS TRANSITION IGNORED',//)
      END
 
C----------------------------------------------------------------------
C SUBROUTINE BRANR
C
C Called by: FTBM
C Calls:     CONV
C
C Purpose: calculate the theoretical branching ratios and compare to the
C experimental ones.
C
C Uses global variables:
C      BRAT   - branching ratio and its error
C      DELTA  - \delta_\lambda: index 1 = electric^2, 2 = magnetic^2, 3 = cross term
C      ELM    - matrix elements
C      EN     - level energy
C      IBRC   - index branching ratios
C      IPRM   - printing flags (see suboption PRT of OP,CONT)
C      KSEQ   - index into ELM for pair of levels, and into EN or SPIN
C      MULTI  - number of matrix elements with given multipolarity
C      NBRA   - number of branching ratios
C
C Formal parameters:
C      Chisq  - chi squared
C      Nwyr   - number of data points contributing to chi squared
C      Chilo  - chi squared of logs

      SUBROUTINE BRANR(Chisq,Nwyr,Chilo)
      IMPLICIT NONE
      REAL*8 ch1 , ch2 , Chilo , Chisq , CONV , eng1 , eng2 , u
      INTEGER*4 i1 , i2 , iflg , iout , itt , j1 , j2 , 
     &          k , lab1 , lab2 , mul2
      INTEGER*4 n1 , n2 , Nwyr
      REAL*8 EN , SPIN , ACCUR , DIPOL , ZPOL , ACCA
      INTEGER*4 ISO
      COMMON /COEX  / EN(100) , SPIN(100) , ACCUR , DIPOL , ZPOL , 
     &                ACCA , ISO
      INTEGER*4 LAMDA , LEAD , LDNUM , LAMMAX , MULTI
      COMMON /CLCOM / LAMDA(8) , LEAD(2,1500) , LDNUM(8,100) , LAMMAX , 
     &                MULTI(8)
      REAL*8 BRAT
      INTEGER*4 IBRC , NBRA
      COMMON /BRNCH / BRAT(50,2) , IBRC(2,50) , NBRA
      REAL*8 DELTA, ENDEC, ENZ
      INTEGER*4 ITMA
      COMMON /TRA   / DELTA(1500,3) , ENDEC(1500) , ITMA(50,200) , 
     &                ENZ(200)
      INTEGER*4 IPRM
      COMMON /PRT   / IPRM(20)
      REAL*8 ELM , ELMU , ELML , SA
      COMMON /COMME / ELM(1500) , ELMU(1500) , ELML(1500) , SA(1500)
      REAL*8 TAU
      INTEGER*4 KSEQ
      COMMON /LEV   / TAU(100) , KSEQ(1500,4)
C     If no branching ratios were defined, return doing nothing
      IF ( NBRA.EQ.0 ) RETURN

C     If printing option is on, print something      
      IF ( IPRM(3).EQ.-1 ) WRITE (22,99001)
99001 FORMAT (1X,///10X,'EXP. AND CALCULATED BRANCHING RATIOS',//5X,
     &        'NS1',5X,'NF1',5X,'NS2',5X,'NF2',5X,'RATIO(1:2)',9X,
     &        'ERROR',7X,'CALC.RATIO',5X,'(EXP-CAL)/ERROR',//)

      Nwyr = Nwyr + NBRA ! Add to number of datapoints contributing to chisqr
      mul2 = MULTI(1) + MULTI(2)
      DO k = 1 , NBRA ! For each branching ratio
         ch1 = 0.
         ch2 = 0.
         iflg = 1
         itt = 1
         iout = 0
         n1 = IBRC(1,k) ! 1st matrix element
         n2 = IBRC(2,k) ! 2nd matrix element
         i1 = KSEQ(n1,1) ! Index of n1'th level
         i2 = KSEQ(n2,1) ! Index of n2'th level
         eng1 = EN(KSEQ(n1,3)) - EN(KSEQ(n1,4)) ! Energy of gamma 1
         eng2 = EN(KSEQ(n2,3)) - EN(KSEQ(n2,4)) ! Energy of gamma 2
         IF ( i1.NE.0 ) THEN
            IF ( i1.LE.MULTI(1) ) lab1 = 1
            IF ( i1.GT.MULTI(1) .AND. i1.LE.mul2 ) lab1 = 2
            IF ( i1.GT.mul2 ) lab1 = 3
         ENDIF
         IF ( i2.NE.0 ) THEN
            IF ( i2.LE.MULTI(1) ) lab2 = 1
            IF ( i2.GT.MULTI(1) .AND. i2.LE.mul2 ) lab2 = 2
            IF ( i2.GT.mul2 ) lab2 = 3
         ENDIF
         IF ( i1.NE.0 ) ch1 = ELM(i1)*ELM(i1)*DELTA(n1,1)
     &                        /(1.+CONV(eng1,lab1))
         IF ( i2.NE.0 ) ch2 = ELM(i2)*ELM(i2)*DELTA(n2,1)
     &                        /(1.+CONV(eng2,lab2))
         j1 = KSEQ(n1,2)
         IF ( j1.NE.0 ) THEN
            iflg = iflg + 1
            lab1 = lab1 + 2
            ch1 = ch1 + ELM(j1)*ELM(j1)*DELTA(n1,2)/(1.+CONV(eng1,lab1))
         ENDIF
         j2 = KSEQ(n2,2)
         IF ( j2.NE.0 ) THEN
            iflg = iflg + 1
            lab2 = lab2 + 2
            ch2 = ch2 + ELM(j2)*ELM(j2)*DELTA(n2,2)/(1.+CONV(eng2,lab2))
         ENDIF
         u = (ch1/ch2-BRAT(k,1))/BRAT(k,2)
         Chisq = Chisq + u*u
         Chilo = Chilo + (BRAT(k,1)*LOG(ch1/ch2/BRAT(k,1))/BRAT(k,2))**2
         IF ( IPRM(3).EQ.-1 ) WRITE (22,99002) KSEQ(n1,3) , KSEQ(n1,4) ,
     &                               KSEQ(n2,3) , KSEQ(n2,4) , BRAT(k,1)
     &                               , BRAT(k,2) , ch1/ch2 , -u
99002    FORMAT (5X,3(1I2,6X),1I2,5X,3(1F10.5,5X),5X,1F4.1)
      ENDDO ! Loop on branching ratios
       
      IF ( IPRM(3).EQ.-1 ) IPRM(3) = 0 ! Turn off printing option, so we don't print twice
      END
 
C----------------------------------------------------------------------
C SUBROUTINE LIMITS
C
C Called by: KONTUR, MINI
C
C Purpose: to constrain the matrix elements to within the limits specified by
C the user.
C
C Uses global variables:
C      ELM    - matrix elements
C      ELML   - lower limit on matrix elements
C      ELMU   - upper limit on matrix elements
C      IVAR   - indicates a limit or correlation is set
C      MEMAX  - number of matrix elements
 
      SUBROUTINE LIMITS
      IMPLICIT NONE
      INTEGER*4 j
      INTEGER*4 MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR
      COMMON /CEXC  / MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR(1500)
      REAL*8 ELM , ELMU , ELML , SA
      COMMON /COMME / ELM(1500) , ELMU(1500) , ELML(1500) , SA(1500)
      DO j = 1 , MEMAX ! Loop over matrix elements
         IF ( IVAR(j).NE.0 ) THEN ! If not fixed
            IF ( ELM(j).GT.ELMU(j) .OR. ELM(j).LT.ELML(j) ) THEN
               IF ( ELM(j).GT.ELMU(j) ) THEN
                  ELM(j) = ELMU(j)
                  WRITE (22,99001) j , ELM(j)
               ELSE
                  ELM(j) = ELML(j)
                  WRITE (22,99001) j , ELM(j)
               ENDIF
            ENDIF
         ENDIF
      ENDDO
99001 FORMAT (2X,'Warning - matrix element ',1I3,' reset to ',1F10.6)
      END
 
C----------------------------------------------------------------------
C SUBROUTINE SZEREG
C
C Called by: READY
C
C Purpose: sort out the decay sequence
C
C Note: szereg is Polish for "series"
C
C Uses global variables:
C      DYEX   - error on experimental yield
C      IY     - index for yields
C      YEXP   - experimental yield
C
C Formal parameters:
C      Lst    - first yield in set
C      Ls     - last yield in set
C      L      - number of dataset
 
      SUBROUTINE SZEREG(Lst,Ls,L)
      IMPLICIT NONE
      REAL*8 dyh , yh
      INTEGER*4 ia , ib , ih , inx , k , L , Ls , lsp , Lst , lst1
      REAL*8 YEXP, CORF , DYEX , UPL , YNRM
      INTEGER*4 IY , NYLDE , IDRN , ILE
      COMMON /YEXPT / YEXP(32,1500) , IY(1500,32) , CORF(1500,32) , 
     &                DYEX(32,1500) , NYLDE(50,32) , UPL(32,50) , 
     &                YNRM(32,50) , IDRN , ILE(32)
      IF ( Lst.EQ.Ls ) RETURN
      lst1 = Lst
      lsp = Ls - 1

 100  ia = IY(lst1,L)
      IF ( ia.GT.1000 ) ia = ia/1000

      inx = lst1
      DO k = lst1 , lsp
         ib = IY(k+1,L)
         IF ( ib.GT.1000 ) ib = ib/1000
         ia = MIN(ia,ib)
         IF ( ia.EQ.ib ) inx = k + 1
      ENDDO

C     Swap them
      IF ( inx.NE.lst1 ) THEN
         ih = IY(lst1,L)
         IY(lst1,L) = IY(inx,L)
         IY(inx,L) = ih
         yh = YEXP(L,lst1)
         dyh = DYEX(L,lst1)
         YEXP(L,lst1) = YEXP(L,inx)
         DYEX(L,lst1) = DYEX(L,inx)
         YEXP(L,inx) = yh
         DYEX(L,inx) = dyh
      ENDIF

      lst1 = lst1 + 1
      IF ( lst1.GT.lsp ) RETURN
      GOTO 100
      END
 
C----------------------------------------------------------------------
C SUBROUTINE SIXEL
C
C Called by: CEGRY
C
C Purpose:
C
C Uses global variables:
C      ARM    - excitation amplitudes of substates.
C      DEV    -
C      IEXP   - experiment number
C      ITS    - create tape 18 file (OP,CONT switch SEL,)
C      KVAR   -
C
C Formal parameters:
C      Rik    - DEV + YEXP
C      Rv     - difference between experimental and calculated yields
C      Em     - matrix element
C      Jk     -
C      Kk     -
C      Indx   - index of matrix element
C      Lu     -
 
      SUBROUTINE SIXEL(Rik,Rv,Em,Jk,Kk,Indx,Lu)
      IMPLICIT NONE
      REAL*8 a1 , al , al1 , c1 , c2 , Em , Rik , rn , Rv , rx
      INTEGER*4 Indx , j , j1 , Jk , Kk , kk6 , l , l1 , Lu
      COMPLEX*16 ARM
      COMMON /AZ    / ARM(1200,7)
      REAL*8 DEV
      COMMON /ODCH  / DEV(1500)
      REAL*8 EPS, EROOT, FIEX
      INTEGER*4 IEXP, IAXS
      COMMON /KIN   / EPS(50) , EROOT(50) , FIEX(50,2) , IEXP , IAXS(50)
      INTEGER*4 ITS
      COMMON /TRB   / ITS
      INTEGER*4 KVAR
      COMMON /SEL   / KVAR(1500)
      
      kk6 = Kk + 5
      rn = DEV(Lu)
      al = (Rv-rn)*20./Rik
      IF ( ITS.EQ.1 .AND. KVAR(Indx).NE.0 ) WRITE (18,*) Lu , Indx , 
     &     IEXP , al/Em
      al1 = ABS(al)
      IF ( ITS.EQ.2 ) WRITE (18,*) Lu , Indx , IEXP , al1
      IF ( al1.LE.ABS(DIMAG(ARM(kk6,Jk))) ) RETURN

      DO j = Kk , kk6
         a1 = ABS(DIMAG(ARM(j,Jk)))
         IF ( al1.GT.a1 ) THEN
            j1 = j + 1
            DO l = j1 , kk6
               l1 = kk6 + j1 - l
               c1 = DBLE(ARM(l1-1,Jk))
               c2 = DIMAG(ARM(l1-1,Jk))
               ARM(l1,Jk) = DCMPLX(c1,c2)
            ENDDO
            rx = DBLE(Indx)
            ARM(j,Jk) = DCMPLX(rx,al)
            GOTO 99999
         ENDIF
      ENDDO

99999 END
 
C----------------------------------------------------------------------
C SUBROUTINE PRELM
C
C Called by: GOSIA
C
C Purpose: print matrix elements
C
C Uses global variables:
C      ELM    - matrix elements
C      ELML   - lower limits on matrix elements
C      ELMU   - upper limits on matrix elements
C      HLM    - matrix elements before minimisation
C      IVAR   - indicates a limit or correlation is set
C      LDNUM  - number of matrix elements with each multipolarity populating levels
C      LEAD   - pair of levels involved in each matrix element
C      MULTI  - number of matrix elements having a given multipolarity
C      NMAX   - number of levels
C      SPIN   - spin of level
C
C Formal parameters:
C      Iop    - print flag (controls what is written to output).
 
      SUBROUTINE PRELM(Iop)
      IMPLICIT NONE
      REAL*8 b , pv , ste
      INTEGER*4 inx , Iop , isp , j , k , kk , l , m
      CHARACTER*3 wrn
      REAL*8 HLM
      COMMON /HHH   / HLM(1500)
      REAL*8 ELM , ELMU , ELML , SA
      COMMON /COMME / ELM(1500) , ELMU(1500) , ELML(1500) , SA(1500)
      INTEGER*4 MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR
      COMMON /CEXC  / MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR(1500)
      REAL*8 EN , SPIN , ACCUR , DIPOL , ZPOL , ACCA
      INTEGER*4 ISO
      COMMON /COEX  / EN(100) , SPIN(100) , ACCUR , DIPOL , ZPOL , 
     &                ACCA , ISO
      INTEGER*4 LAMDA , LEAD , LDNUM , LAMMAX , MULTI
      COMMON /CLCOM / LAMDA(8) , LEAD(2,1500) , LDNUM(8,100) , LAMMAX , 
     &                MULTI(8)
      INTEGER*4 NMAX , NDIM , NMAX1
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      inx = 0
      WRITE (22,99001)
99001 FORMAT (2X/40X,'MATRIX ELEMENTS',//)

      DO j = 1 , 8
         m = MULTI(j)
         IF ( m.NE.0 ) THEN
            WRITE (22,99002) j
99002       FORMAT (5X,'MULTIPOLARITY=',1I1)
            IF ( Iop.EQ.1 ) WRITE (22,99003)
99003       FORMAT (4X,'INDEX',3X,'NF',5X,'NS',10X,'ME')
            IF ( Iop.EQ.2 ) WRITE (22,99004)
99004       FORMAT (4X,'INDEX',3X,'NF',5X,'NS',10X,'ME',15X,'LIMITS')
            IF ( Iop.EQ.3 ) WRITE (22,99005)
99005       FORMAT (4X,'INDEX',3X,'NF',5X,'NS',10X,'ME',10X,'PC CHANGE',
     &              5X,'RED. TRANS. PROB.')

            DO k = 1 , NMAX
               l = LDNUM(j,k)
               IF ( l.NE.0 ) THEN
                  DO kk = 1 , l
                     inx = inx + 1

                     IF ( Iop.EQ.2 ) THEN ! Iop is 2
                        IF ( IVAR(inx).EQ.0 ) THEN ! Fixed
                           WRITE (22,99006) inx , LEAD(1,inx) , 
     &                            LEAD(2,inx) , ELM(inx)
                        ELSEIF ( IVAR(inx).GT.1000 ) THEN ! Correlation
                           WRITE (22,99007) inx , LEAD(1,inx) , 
     &                            LEAD(2,inx) , ELM(inx) , 
     &                            (IVAR(inx)-1000)
                        ELSE ! Limit
                           WRITE (22,99009) inx , LEAD(1,inx) , 
     &                            LEAD(2,inx) , ELM(inx) , ELML(inx) , 
     &                            ELMU(inx)
                        ENDIF

                     ELSEIF ( Iop.EQ.3 ) THEN ! Iop is 3
                        isp = LEAD(2,inx)
                        pv = (ELMU(inx)-ELML(inx))/100.
                        wrn = '   '
                        IF ( (ELM(inx)-ELML(inx)).LT.pv ) wrn = '*?*'
                        IF ( (ELMU(inx)-ELM(inx)).LT.pv ) wrn = '*?*'
                        ste = HLM(inx)
                        b = ELM(inx)*ELM(inx)/(2.*SPIN(isp)+1.)
                        IF ( LEAD(1,inx).EQ.LEAD(2,inx) ) b = 9999999.
                        WRITE (22,99009) inx , LEAD(1,inx) , LEAD(2,inx)
     &                         , ELM(inx) , 100.*(ELM(inx)-ste)/ste , 
     &                         b , wrn
                     ELSE ! Iop is 1
                        WRITE (22,99008) inx , LEAD(1,inx) , LEAD(2,inx)
     &                         , ELM(inx)
                     ENDIF

                  ENDDO ! Loop on kk
               ENDIF ! If l .ne. 0
            ENDDO ! Loop on k
         ENDIF ! If m .ne. 0
      ENDDO ! Loop on j

99006 FORMAT (5X,1I3,4X,1I3,4X,1I3,5X,1F10.5,5X,'FIXED')
99007 FORMAT (5X,1I3,4X,1I3,4X,1I3,5X,1F10.5,5X,'COUPLED TO',1X,1I3)
99008 FORMAT (5X,1I3,4X,1I3,4X,1I3,5X,1F10.5)
99009 FORMAT (5X,1I3,4X,1I3,4X,1I3,3(5X,1F10.5),1A3)
      END
 
C----------------------------------------------------------------------
C SUBROUTINE RECOIL
C
C Called by: ANGULA
C Calls:     ROTATE
C
C Purpose: correct for relativistic effects of recoiling nucleus.
C
C Formal parameters:
C      Alab   - matrix of (l,m) pairs in lab frame
C      Attl   - matrix of (l,m) pairs in rotated frame
C      Beta   - beta of recoil
C      Theta  - angle to rotate
C
C We transform into the frame of the recoiling nucleus, correct according to
C the method of Lesser and then rotate back to the laboratory frame.
 
      SUBROUTINE RECOIL(Alab,Attl,Beta,Theta)
      IMPLICIT NONE
      REAL*8 Alab , atemp , Attl , Beta , betasq , dum , hold , test , 
     &       Theta
      INTEGER*4 i , i1 , j , l , m
      DIMENSION Alab(9,9) , Attl(9,9) , atemp(16)
      
      hold = Alab(1,1)
      IF ( ABS(hold).LT.1.E-9 ) RETURN

C     Rotate into frame of recoiling nucleus
      CALL ROTATE(Alab,Attl,-Theta,7,2)

C     Correct for relativistic effects
      Attl(2,1) = (2./SQRT(15.))*(SQRT(5.)*Attl(1,1)-Attl(3,1))
      Attl(2,2) = -Attl(3,2)/SQRT(5.)
      Attl(4,1) = (4./SQRT(35.))*(3.*Attl(3,1)-SQRT(5.)*Attl(5,1))
      Attl(4,2) = (8.*SQRT(2.)*Attl(3,2)-5.*SQRT(3.)*Attl(5,2))
     &            /SQRT(35.)
      Attl(4,3) = (2./SQRT(7.))*(2.*Attl(3,3)-SQRT(3.)*Attl(5,3))
      Attl(4,4) = -Attl(5,4)
      Attl(6,1) = (10./SQRT(11.))*(Attl(5,1)-(3./SQRT(13.))*Attl(7,1))
      Attl(6,2) = (1./SQRT(11.))
     &            *(4.*SQRT(6.)*Attl(5,2)-5.*SQRT(35./13.)*Attl(7,2))
      Attl(6,3) = SQRT(4./11.)
     &            *(SQRT(21.)*Attl(5,3)-10.*SQRT(2./13.)*Attl(7,3))
      Attl(6,4) = SQRT(1./11.)*(8.*Attl(5,4)-15.*SQRT(3./13.)*Attl(7,4))
      Attl(6,5) = SQRT(4./11.)*(3.*Attl(5,5)-5.*SQRT(5./13.)*Attl(7,5))
      Attl(6,6) = -Attl(7,6)*SQRT(25./13.)
      Attl(8,1) = (56./SQRT(195.))*Attl(7,1)
      Attl(8,2) = (32./SQRT(65.))*Attl(7,2)
      Attl(8,3) = (8.*SQRT(3./13.))*Attl(7,3)
      Attl(8,4) = (16.*SQRT(2./39.))*Attl(7,4)
      Attl(8,5) = (8.*SQRT(11./65.))*Attl(7,5)
      Attl(8,6) = (16.*SQRT(2./65.))*Attl(7,6)
      Attl(8,7) = (8./SQRT(15.))*Attl(7,7)
      DO l = 2 , 8 , 2
         DO m = 1 , l
            Attl(l,m) = Beta*Attl(l,m)
         ENDDO
      ENDDO
      betasq = Beta*Beta
      IF ( betasq.GE.1.0E-10 ) THEN
         i1 = 0
         DO i = 1 , 7 , 2
            DO j = 1 , i
               i1 = i1 + 1
               atemp(i1) = Attl(i,j)
            ENDDO
         ENDDO
         dum = (2./5.)*SQRT(5.)*atemp(1) - (10./7.)*atemp(2) + (12./35.)
     &         *SQRT(5.)*atemp(5)
         Attl(3,1) = atemp(2) + betasq*dum
         dum = -(17./14.)*atemp(3) + (2./7.)*SQRT(6.)*atemp(6)
         Attl(3,2) = atemp(3) + betasq*dum
         dum = -(4./7.)*atemp(4) + (2./7.)*SQRT(3.)*atemp(7)
         Attl(3,3) = atemp(4) + betasq*dum
         dum = (8./7.)*SQRT(5.)*atemp(2) - (380./77.)*atemp(5)
     &         + (100./11.)*SQRT(1./13.)*atemp(10)
         Attl(5,1) = atemp(5) + betasq*dum
         dum = (20./21.)*SQRT(6.)*atemp(3) - (723./154.)*atemp(6)
     &         + (20./11.)*SQRT(70./39.)*atemp(11)
         Attl(5,2) = atemp(6) + betasq*dum
         dum = (20./21.)*SQRT(3.)*atemp(4) - (306./77.)*atemp(7)
     &         + (40./11.)*SQRT(14./39.)*atemp(12)
         Attl(5,3) = atemp(7) + betasq*dum
         dum = -(61./22.)*atemp(8) + (40./11.)*SQRT(3./13.)*atemp(13)
         Attl(5,4) = atemp(8) + betasq*dum
         dum = -(12./11.)*atemp(9) + (20./11.)*SQRT(5./13.)*atemp(14)
         Attl(5,5) = atemp(9) + betasq*dum
         dum = (210./11.)*SQRT(1./13.)*atemp(5) - (574./55.)*atemp(10)
         Attl(7,1) = atemp(10) + betasq*dum
         dum = (14./11.)*SQRT(210./13.)*atemp(6) - (1121./110.)
     &         *atemp(11)
         Attl(7,2) = atemp(11) + betasq*dum
         dum = (28./11.)*SQRT(42./13.)*atemp(7) - (104./11.)*atemp(12)
         Attl(7,3) = atemp(12) + betasq*dum
         dum = (84./11.)*SQRT(3./13.)*atemp(8) - (181./22.)*atemp(13)
         Attl(7,4) = atemp(13) + betasq*dum
         dum = (42./11.)*SQRT(5./13.)*atemp(9) - (358./55.)*atemp(14)
         Attl(7,5) = atemp(14) + betasq*dum
         Attl(7,6) = atemp(15)*(1.-(43./10.)*betasq)
         Attl(7,7) = atemp(16)*(1.-(8./5.)*betasq)
         Attl(9,1) = (672./5.)*SQRT(1./221.)*atemp(10)*betasq
         Attl(9,2) = (144./5.)*SQRT(21./221.)*atemp(11)*betasq
         Attl(9,3) = 36.*SQRT(12./221.)*atemp(12)*betasq
         Attl(9,4) = 24.*SQRT(22./221.)*atemp(13)*betasq
         Attl(9,5) = (144./5.)*SQRT(11./221.)*atemp(14)*betasq
         Attl(9,6) = (72./5.)*SQRT(2./17.)*atemp(15)*betasq
         Attl(9,7) = (24./5.)*SQRT(7./17.)*atemp(16)*betasq
      ENDIF

C     Rotate back into laboratory frame
      CALL ROTATE(Attl,Alab,Theta,9,1)
      test = ABS(1.0-Alab(1,1)/hold)
      IF ( test.GT.1.0E-07 ) THEN
         WRITE (22,99001) test
99001    FORMAT (' ERROR IN ROTATION',1X,1E10.3/)
      ENDIF
      END
 
C----------------------------------------------------------------------
C SUBROUTINE ROTATE
C
C Called by: RECOIL
C Calls:     DJMM
C
C Purpose: rotate the frame of reference
C
C Formal parameters:
c      Alab   - statistical tensor in lab frame
C      Attl   - statistical tensor rotated
C      Theta  - angle to rotate by
C      K2     - maximum dimension to rotate
C      Kd     - step over dimension
C
C We use:
C \rho_{k \xi} = (-1)^\xi *
C  \sum_\xi^\prime{i^\xi\prime d^k_{xi^\prime \xi}
C     ({\pi + \theta \over 2}) \rho_{k \xi^\prime}}

      SUBROUTINE ROTATE(Alab,Attl,Theta,K2,Kd)
      IMPLICIT NONE
      REAL*8 Alab , Attl , djarg , DJMM , dkkk , sum , Theta
      INTEGER*4 idj , idm , idmp , j , k , K2 , ka , kappa , kapri , Kd
      DIMENSION Alab(9,9) , Attl(9,9)
      
      IF ( ABS(Theta).GT..01 ) THEN
         djarg = Theta
         DO ka = 1 , K2 , Kd
            idj = ka - 1
            DO kappa = 1 , ka
               idmp = kappa - 1
               sum = 0.0
               DO kapri = 1 , ka
                  idm = kapri - 1
                  dkkk = DJMM(djarg,idj,idm,idmp)
                  sum = sum + dkkk*Alab(ka,kapri)
               ENDDO
               IF ( ka.NE.1 ) THEN
                  DO kapri = 2 , ka
                     idm = -kapri + 1
                     dkkk = DJMM(djarg,idj,idm,idmp)
                     sum = sum + dkkk*Alab(ka,kapri)*(-1.0)**(kapri-1)
                  ENDDO
               ENDIF
               Attl(ka,kappa) = sum
            ENDDO
         ENDDO
         GOTO 99999
      ENDIF
      DO j = 1 , 9
         DO k = 1 , 9
            Attl(j,k) = Alab(j,k)
         ENDDO
      ENDDO
99999 END
 
C----------------------------------------------------------------------
C SUBROUTINE YLM1
C
C Called by: ANGULA
C
C Purpose: evaluate the odd spherical harmonics.
C
C Formal parameters:
C      Theta  - theta for which to evaluate
C      Ylmr   - return value for that theta
C
C Ylmr(l,m) = 1 / \sqrt{4 \pi} Y_{l - 1}^{m - 1}
C
C Note the factor of 1 / \sqrt{4 \pi} compared to the orthonormal spherical
C harmonics.
C
C Note also that YLM1 and YLM have some values in common.
C e.g. YLM1(5,3) = YLM(2,3)
 
      SUBROUTINE YLM1(Theta,Ylmr)
      IMPLICIT NONE
      REAL*8 ct , ctsq , st , Theta , Ylmr
      INTEGER*4 i , j , l , m
      DIMENSION Ylmr(9,9) , st(9)
      
      ct = COS(Theta)
      ctsq = ct*ct
      st(1) = SIN(Theta)
      DO i = 2 , 9
         j = i - 1
         st(i) = st(j)*st(1)
      ENDDO
      DO l = 2 , 9
         DO m = 1 , 9
            Ylmr(l,m) = 0.0
         ENDDO
      ENDDO
      Ylmr(2,2) = -SQRT(6.)/2.
      Ylmr(2,1) = SQRT(3.)*ct
      Ylmr(3,3) = SQRT(30.)/4.
      Ylmr(3,2) = -(SQRT(30.)/2.)*ct
      Ylmr(3,1) = (SQRT(5.)/2.)*(3.*ctsq-1.)
      Ylmr(4,4) = -SQRT(35.)/4.
      Ylmr(4,3) = (SQRT(210.)/4.)*ct
      Ylmr(4,2) = -(SQRT(21.)/4.)*(5.*ctsq-1.)
      Ylmr(4,1) = (SQRT(7.)/2.)*ct*(5.*ctsq-3.)
      Ylmr(5,5) = 3.*SQRT(70.)/16.
      Ylmr(5,4) = -(3.*SQRT(35.)/4.)*ct
      Ylmr(5,3) = (3.*SQRT(10.)/8.)*(7.*ctsq-1.)
      Ylmr(5,2) = -(3.*SQRT(5.)/4.)*ct*(7.*ctsq-3.)
      Ylmr(5,1) = (3./8.)*((35.*ctsq-30.)*ctsq+3.)
      Ylmr(6,6) = -3.*SQRT(77.)/16.
      Ylmr(6,5) = (3.*SQRT(770.)/16.)*ct
      Ylmr(6,4) = -(SQRT(385.)/16.)*(9.*ctsq-1.)
      Ylmr(6,3) = (SQRT(2310.)/8.)*ct*(3.*ctsq-1.)
      Ylmr(6,2) = -(SQRT(330.)/16.)*((21.*ctsq-14.)*ctsq+1.)
      Ylmr(6,1) = (SQRT(11.)/8.)*ct*((63.*ctsq-70.)*ctsq+15.)
      Ylmr(7,7) = SQRT(3003.)/32.
      Ylmr(7,6) = -(3.*SQRT(1001.)/16.)*ct
      Ylmr(7,5) = (3.*SQRT(182.)/32.)*(11.*ctsq-1.)
      Ylmr(7,4) = -(SQRT(1365.)/16.)*ct*(11.*ctsq-3.)
      Ylmr(7,3) = (SQRT(1365.)/32.)*((33.*ctsq-18.)*ctsq+1.)
      Ylmr(7,2) = -(SQRT(546.)/16.)*ct*((33.*ctsq-30.)*ctsq+5.)
      Ylmr(7,1) = (SQRT(13.)/16.)*(((231.*ctsq-315.)*ctsq+105.)*ctsq-5.)
      Ylmr(8,8) = -3.*SQRT(1430.)/64.
      Ylmr(8,7) = (3.*SQRT(5005.)/32.)*ct
      Ylmr(8,6) = -(3.*SQRT(770.)/64.)*(13.*ctsq-1.)
      Ylmr(8,5) = (3.*SQRT(770.)/32.)*(13.*ctsq-3.)*ct
      Ylmr(8,4) = -(3.*SQRT(70.)/64.)*((143.*ctsq-66.)*ctsq+3.)
      Ylmr(8,3) = (3.*SQRT(35.)/32.)*((143.*ctsq-110.)*ctsq+15.)*ct
      Ylmr(8,2) = -(SQRT(210.)/64.)
     &            *(((429.*ctsq-495.)*ctsq+135.)*ctsq-5.)
      Ylmr(8,1) = (SQRT(15.)/16.)
     &            *(((429.*ctsq-693.)*ctsq+315.)*ctsq-35.)*ct
      Ylmr(9,9) = 3.*SQRT(24310.)/256.
      Ylmr(9,8) = -(3.*SQRT(24310.)/64.)*ct
      Ylmr(9,7) = (SQRT(7293.)/64.)*(15.*ctsq-1.)
      Ylmr(9,6) = -(3.*SQRT(34034.)/64.)*(5.*ctsq-1.)*ct
      Ylmr(9,5) = (3.*SQRT(2618.)/128.)*((65.*ctsq-26.)*ctsq+1.)
      Ylmr(9,4) = -(SQRT(39270.)/64.)*((39.*ctsq-26.)*ctsq+3.)*ct
      Ylmr(9,3) = (3.*SQRT(595.)/64.)
     &            *(((143.*ctsq-143.)*ctsq+33.)*ctsq-1.)
      Ylmr(9,2) = -(3.*SQRT(34.)/64.)
     &            *(((715.*ctsq-1001.)*ctsq+385.)*ctsq-35.)*ct
      Ylmr(9,1) = (SQRT(17.)/128.)
     &            *((((6435.*ctsq-12012.)*ctsq+6930.)*ctsq-1260.)
     &            *ctsq+35.)
      DO l = 2 , 9
         Ylmr(l,1) = Ylmr(l,1)*.0795774715 ! 0.0795774715 = 1 / (4 pi)
         DO m = 2 , l
            Ylmr(l,m) = Ylmr(l,m)*st(m-1)*.0795774715 ! 0.0795774715 = 1 / (4 pi)
         ENDDO
      ENDDO
      END
 
C----------------------------------------------------------------------
C SUBROUTINE FIINT
C
C Called by: ANGULA
C
C Purpose: integrate over phi in frame of recoiling nucleus
C
C Formal parameters:
C      Fi0    - phi_0
C      Fi1    - phi_1
C      At     - return value
C      Ixs    - axial symmetry flag

      SUBROUTINE FIINT(Fi0,Fi1,At,Ixs)
      IMPLICIT NONE
      REAL*8 At , Fi0 , Fi1 , wsp
      INTEGER*4 Ixs , j , jf , js , m , mm
      DIMENSION At(28)
      
      IF ( Ixs.NE.0 ) THEN
         DO m = 2 , 7
            js = m/2
            mm = m - 1
            wsp = (SIN(mm*Fi1)-SIN(mm*Fi0))/mm
            js = js*7 + m
            jf = m + 21
            DO j = js , jf , 7
               At(j) = At(j)*wsp
            ENDDO
         ENDDO
         wsp = Fi1 - Fi0
      ENDIF
      IF ( Ixs.EQ.0 ) wsp = 6.283185308 ! 6.283185308 = 2 * pi
      DO j = 1 , 4
         js = (j-1)*7 + 1
         At(js) = At(js)*wsp
      ENDDO
      END
 
C----------------------------------------------------------------------
C SUBROUTINE FIINT1
C
C Called by: ANGULA
C
C Purpose: integrate over phi in lab frame
C
C Formal parameters:
C      Fi0    - phi_0
C      Fi1    - phi_1
C      Alab   - return value
C      Ixs    - axial symmetry flag

      SUBROUTINE FIINT1(Fi0,Fi1,Alab,Ixs)
      IMPLICIT NONE
      REAL*8 Alab , Fi0 , Fi1 , wsp
      INTEGER*4 Ixs , j , m , mm
      DIMENSION Alab(9,9)
      
      IF ( Ixs.NE.0 ) THEN
         DO m = 2 , 9
            mm = m - 1
            wsp = (SIN(mm*Fi1)-SIN(mm*Fi0))/mm
            DO j = 1 , 9
               Alab(j,m) = Alab(j,m)*wsp
            ENDDO
         ENDDO
         wsp = Fi1 - Fi0
      ENDIF
      IF ( Ixs.EQ.0 ) wsp = 6.283185308 ! 6.283185308 = 2 * pi
      DO j = 1 , 9
         Alab(j,1) = Alab(j,1)*wsp
      ENDDO
      END
 
C----------------------------------------------------------------------
C SUBROUTINE TAPMA
C
C Called by: GOSIA
C
C Purpose: read parameters for sensitivity maps
C
C Uses global variables:
C      DS     - differential cross section
C      XV     - energy meshpoints (sometimes theta meshpoints) where we calculate exact Coulex
C      YGN    - gamma yield calculated without correction to angular distribution from finite recoil distance
C      ZETA   - various coefficients
C
C Formal parameters:
C      Lx     - experiment number
C      Iske   -
C      Isko   -
C      Iskf   -
C      Nflr   -
C      Idr    - number of decays
C      Nco    -
C      Nft    - error flag: 0 = no error, 1 = error
C      Enb    - energy of meshpoint read from file
C
C Note that unit 14 is used internally for the purpose of sensitivity
C maps.
 
      SUBROUTINE TAPMA(Lx,Iske,Isko,Iskf,Nflr,Idr,Nco,Nft,Enb)
      IMPLICIT NONE
      REAL*8 emn , emx , en0 , Enb , tmn , tmx , tta
      INTEGER*4 Idr , Iske , Iskf , Isko , j , jf , jj , js , k , 
     &          Lx , lx1 , na , Nco , ne , nfil , nfilt , Nflr , Nft
      INTEGER*4 ng , ng1 , ntt
      REAL*8 XV, YV, ZV, DSG, DSE, DS
      COMMON /VLIN  / XV(101) , YV(101) , ZV(101) , DSG(101) ,
     &                DSE(101) , DS
      REAL*8 ZETA
      INTEGER*4 LZETA
      COMMON /CCOUP / ZETA(155600) , LZETA(8)
      REAL*8 YGN , YGP
      INTEGER*4 IFMO
      COMMON /YTEOR / YGN(1500) , YGP(1500) , IFMO
      Nft = 0
      nfilt = 0
      REWIND 14

C     Skip over unwanted records
      IF ( Iske.NE.0 ) THEN
 50      READ (14,*) ne , ntt , emn , emx , tmn , tmx , na , tmx , tmx ,
     &               tmx
         nfil = ne*ntt*na
         nfilt = nfilt + nfil
         DO j = 1 , nfil
            READ (14,*) lx1 , Enb , tta , ng , DS , (YGN(k),k=1,Idr)
         ENDDO
         IF ( nfilt.NE.Iske ) GOTO 50
      ENDIF

      IF ( Nco.EQ.0 ) RETURN

C     Read record
      READ (14,*) ne , ntt , emn , emx , tmn , tmx , na , tmx , tmx , 
     &            tmx
      IF ( Isko.NE.0 ) THEN
         DO j = 1 , Isko
            READ (14,*) lx1 , Enb , tta , ng , DS , (YGN(k),k=1,Idr)
         ENDDO
      ENDIF

      DO j = 1 , Nflr
         js = (j-1)*Idr + 1
         jf = js + Idr - 1
         READ (14,*) lx1 , Enb , tta , ng1 , DS , (ZETA(k),k=js,jf)
         IF ( lx1.NE.Lx ) Nft = 1
         IF ( Nft.EQ.1 ) GOTO 100
         XV(j) = tta/57.2957795
         IF ( Iskf.NE.0 .AND. j.NE.Nflr ) THEN
            DO jj = 1 , Iskf
               READ (14,*) lx1 , en0 , tta , ng , DS , (YGN(k),k=1,Idr)
            ENDDO
         ENDIF
      ENDDO

      RETURN
 100  WRITE (22,99001)
99001 FORMAT (10X///10X,'TAPE READ ERROR'/10X,'JOB ABORTED')
      END
 
C----------------------------------------------------------------------
C FUNCTION SIMIN
C
C Called by: GOSIA
C
C Purpose: Integrate under a curve defined by an array.
C
C Formal parameters:
C      Np     - number of points in array Y
C      H      - step between points
C      Y      - array of points
C
C Return value:
C      Integral under the array
C
 
      REAL*8 FUNCTION SIMIN(Np,H,Y)
      IMPLICIT NONE
      REAL*8 ee , H , sm , Y
      INTEGER*4 ik , in , Np
      DIMENSION Y(101)

      IF ( Np.GE.3 ) THEN
         ik = Np - 2
         sm = Y(1) + Y(Np)
         DO in = 1 , ik
            ee = in/2.
            sm = sm + 2.*Y(in+1)/(1.+INT(ee)-ee)
         ENDDO
         SIMIN = sm*H/3.
         RETURN
      ELSEIF ( Np.EQ.1 ) THEN
         SIMIN = Y(1)
         GOTO 99999
      ENDIF
      SIMIN = (Y(1)+Y(2))*H/2.
99999 END
 
C----------------------------------------------------------------------
C SUBROUTINE MIXUP
C
C Called by: GOSIA
C Calls:     RNDM
C
C Purpose: set the matrix elements to random values, as a starting value.
C
C Uses global variables:
C      ELM    - matrix elements
C      ELML   - lower limit on matrix elements
C      ELMU   - upper limit on matrix elements
C      IVAR   - indicates a limit or correlation is set
C      MEMAX  - number of matrix elements
C      SA     - ratio of elements for correlated elements
C      SE     - seed for random number generator
C
C It is called when the user gives the option OP,RAND
C
C Note that if IVAR = 0, then the matrix element is fixed, so we don't do
C anything here. If it is >= 10000, this means it is correlated to another
C matrix element, so use the correlation to determine the new value, which
C may have been changed when we randomized.
 
      SUBROUTINE MIXUP
      IMPLICIT NONE
      REAL*8 RNDM
      INTEGER*4 k , k1
      REAL*8 ELM , ELMU , ELML , SA
      COMMON /COMME / ELM(1500) , ELMU(1500) , ELML(1500) , SA(1500)
      INTEGER*4 MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR
      COMMON /CEXC  / MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR(1500)
      REAL*8 SE
      COMMON /XRA   / SE
C     Randomize all that are not fixed or correlated
      DO k = 1 , MEMAX ! For each matrix element
         IF ( IVAR(k).NE.0 .AND. IVAR(k).LE.999 ) ! Not fixed or correlated
     &        ELM(k) = ELML(k) + RNDM(SE)*(ELMU(k)-ELML(k))
      ENDDO

C     Now adjust the correlated elements, since we may have changed the
C     element to which it is correlated.
      DO k = 1 , MEMAX
         IF ( IVAR(k).GE.999 ) THEN ! Correlated
            k1 = IVAR(k) - 1000 ! Index to which it is correlated
            IF ( ABS(ELMU(k1)).LT.1.E-9 ) THEN
               ELM(k) = 0.
            ELSE
               ELM(k) = ELM(k1)*SA(k) ! SA is the ratio we require
            ENDIF
         ENDIF
      ENDDO
      END
 
C----------------------------------------------------------------------
C FUNCTION FXIS1
C
C Called by: NEWCAT
C
C Purpose: return -1 * sign(xi) except for N = 2,3,5 and 6
C
C Uses global variables:
C      XI     - xi coupling coefficients
C
C Formal parameters:
C      I      - index into XI array
C      N      - 
C
C Return value:
C      sign of xi
      
      REAL*8 FUNCTION FXIS1(I,N)
      IMPLICIT NONE
      INTEGER*4 I , N
      REAL*8 XI
      COMMON /CXI   / XI(1500)
      IF ( N.EQ.2 .OR. N.EQ.3 .OR. N.EQ.5 .OR. N.EQ.6 ) THEN
         FXIS1 = 1.
         GOTO 99999
      ENDIF
      FXIS1 = -SIGN(1.D0,XI(I))
99999 END
 
C----------------------------------------------------------------------
C FUNCTION FXIS2
C
C Called by: NEWCAT
C
C Purpose: return -1 * sign(xi) for N = 2,3,5 and 6
C
C Uses global variables:
C      XI     - xi coupling coefficients
C
C Formal parameters:
C      I      - index into XI array
C      N      - 
C
C Return value:
C      sign of xi

      REAL*8 FUNCTION FXIS2(I,N)
      IMPLICIT NONE
      INTEGER*4 I , N
      REAL*8 XI
      COMMON /CXI   / XI(1500)
      IF ( N.EQ.2 .OR. N.EQ.3 .OR. N.EQ.5 .OR. N.EQ.6 ) THEN
         FXIS2 = -SIGN(1.D0,XI(I))
         GOTO 99999
      ENDIF
      FXIS2 = 1.
99999 END
 
C----------------------------------------------------------------------
C SUBROUTINE PODZIEL
C
C Called by: APRAM
C
C Purpose: subdivide matrix operators if the summation doesn't converge.
C
C Note: podziel is Polish for "split"
C
C Uses global variables:
C      IDIVE  - number of subdivisions
C      LP2    - maximum number of matrix elements (1500)
C      QAPR   - approximate Coulomb amplitudes
C
C We use the identity: exp(A) \bar{a} = exp(A/2) exp(A/2)\bar{a}.
C
C Formal parameters:
C      I      - flag (I=1,2,3) I=3 means initialise
C      J      - experiment number
 
      SUBROUTINE PODZIEL(I,J)
      IMPLICIT NONE
      INTEGER*4 I , J , k , l , l1 , l2
      INTEGER*4 LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &          LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      REAL*8 QAPR
      INTEGER*4 IAPR , ISEX
      COMMON /APRCAT/ QAPR(1500,2,7) , IAPR(1500,2) , ISEX(100)
      INTEGER*4 LERF , IDIVE
      COMMON /APRX  / LERF , IDIVE(50,2)
      DATA l1/0/

      IF ( I.NE.3 ) THEN
         IF ( I.EQ.1 ) THEN
            l1 = IDIVE(J,1)
            IDIVE(J,1) = l1 + 1
            GOTO 100
         ELSE
            l1 = IDIVE(J,2)
            IDIVE(J,2) = l1 + 1
         ENDIF
      ENDIF

      l2 = IDIVE(J,2)
      IF ( I.EQ.3 ) l1 = 1
      DO k = 1 , LP2 ! For each matrix element
         DO l = 1 , 7
            QAPR(k,2,l) = QAPR(k,2,l)*l1/l2
         ENDDO
      ENDDO

      IF ( I.EQ.2 ) WRITE (22,99001) J , IDIVE(J,1) , l2
      IF ( I.NE.3 ) RETURN
      
 100  l2 = IDIVE(J,1)
      IF ( I.EQ.3 ) l1 = 1
      DO k = 1 , LP2 ! For each matrix element
         DO l = 1 , 7
            QAPR(k,1,l) = QAPR(k,1,l)*l1/l2
         ENDDO
      ENDDO
       
      IF ( I.EQ.1 ) WRITE (22,99001) J , l2 , IDIVE(J,2)
      RETURN
      
99001 FORMAT (5X,'*****',1X,'EXP(A) EXPANSION FAILURE!',1X,'*****'/5X,
     &        'EXPERIMENT',1X,1I2,3X,'NEW SUBDIVISION',1X,'(',1I1,',',
     &        1I1,')')
      END
 
C----------------------------------------------------------------------
C SUBROUTINE KLOPOT
C
C Called by: GOSIA
C
C Purpose: trouble shooting (see OP,TROU)
C
C Note: Klopot (with a slash through the "l" is Polish for trouble
C
C Uses global variables:
C      CORF   - internal correction factors
C      ELM    - matrix elements
C      ELML   - lower limit on matrix elements
C      ELMU   - upper limit on matrix elements
C      EP     - bombarding energy
C      IY     - index of experimental yields
C      KVAR   -
C      LP2    - maximum number of matrix elements (1500)
C      MEMAX  - number of matrix elements
C      NEXPT  - number of experiments
C      TLBDG  - theta of particle detector in lab frame (in degrees)
C      ZETA   - various coefficients
C      VINF   - speed of projectile at infinity
C
C Formal parameters:
C      K      - number of experimental yields giving largest and positive
C               components of the derivative of chi squared.
C      Rlr    - print out if matrix element exceeds Rlr.

      SUBROUTINE KLOPOT(K,Rlr)
      IMPLICIT NONE
      REAL*8 a , al , al1 , b , c , ch , d , dy , e , g , g1 , g2 , 
     &       rl , Rlr , sgm , u , umm , ump , ux
      INTEGER*4 i , iex , iexh , iexp , indx , inh , ipf , j , jm , 
     &          jp , K , l , lc , ll , lngt , loc , lu , nf , ni , nm , 
     &          np
      REAL*8 ELM , ELMU , ELML , SA
      COMMON /COMME / ELM(1500) , ELMU(1500) , ELML(1500) , SA(1500)
      REAL*8 YEXP, CORF , DYEX , UPL , YNRM
      INTEGER*4 IY , NYLDE , IDRN , ILE
      COMMON /YEXPT / YEXP(32,1500) , IY(1500,32) , CORF(1500,32) , 
     &                DYEX(32,1500) , NYLDE(50,32) , UPL(32,50) , 
     &                YNRM(32,50) , IDRN , ILE(32)
      INTEGER*4 MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR
      COMMON /CEXC  / MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR(1500)
      INTEGER*4 LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &          LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      REAL*8 ZETA
      INTEGER*4 LZETA
      COMMON /CCOUP / ZETA(155600) , LZETA(8)
      REAL*8 XA , XA1 , EP , TLBDG , VINF
      INTEGER*4 NEXPT , IZ , IZ1
      COMMON /CX    / NEXPT , IZ , XA , IZ1(50) , XA1(50) , EP(50) , 
     &                TLBDG(50) , VINF(50)
      INTEGER*4 KVAR
      COMMON /SEL   / KVAR(1500)
      DATA jm/0/,jp/0/

      REWIND 14
      REWIND 18
      ipf = 1
      lngt = 0
      indx = 1
      REWIND 15
      REWIND 17

      DO i = 1 , MEMAX ! Zero all matrix elements and limits
         ELM(i) = 0.
         ELMU(i) = 0.
         ELML(i) = 0.
      ENDDO

      iexh = 1
 100  g = 0.
      d = 0.
 200  READ (15,*) iex , a , b , c , e
      IF ( iex.NE.iexh ) THEN
         EP(iexh) = g/d
         TLBDG(iexh) = g
         VINF(iexh) = d
         iexh = iex
         BACKSPACE 15
         IF ( iex.NE.0 ) GOTO 100
         REWIND 15
         iexp = 1
      ELSE
         g = g + e*a/c/c
         d = d + a*a/c/c
         lngt = lngt + 1
         GOTO 200
      ENDIF
 300  g1 = 0.
      g2 = 0.
      inh = indx
      iexh = iexp
 400  READ (18,*) lu , indx , iexp , al
      IF ( indx.NE.0 ) THEN
         READ (17,*) ni , nf , sgm , al1
         READ (15,*) iex , a , b , c , e
         IF ( iexp.NE.1 .AND. ipf.NE.1 ) THEN
 420        READ (15,*) iex , a , b , c , e
            IF ( iexp.NE.iex ) GOTO 420
            ipf = 1
         ENDIF
         IF ( indx.EQ.inh ) THEN
            dy = al*al1/b/c
            g1 = e*dy + g1
            g2 = -2.*dy*a + g2
            WRITE (14,*) indx , iexp , ni , nf , dy , a , e , c
            GOTO 400
         ENDIF
      ENDIF
      loc = (iexh-1)*LP2 + inh ! LP2 = 1500
      ipf = 0
      ZETA(loc) = (VINF(iexh)*g1+TLBDG(iexh)*g2)/VINF(iexh)/VINF(iexh)
      inh = indx
      REWIND 15
      BACKSPACE 17
      BACKSPACE 18
      IF ( indx.NE.0 ) GOTO 300
      WRITE (14,*) indx , iexp , ni , nf , dy , a , e , b
      REWIND 14
      REWIND 17
 500  READ (14,*) indx , iexp , ni , nf , dy , a , e , b
      IF ( indx.EQ.0 ) THEN
         WRITE (17,*) indx , iexp , sgm , ni , nf , u
         REWIND 17
         ll = 0
         ch = 0.
 550     READ (17,*) indx , iexp , sgm , ni , nf , u
         IF ( indx.EQ.0 ) THEN
            WRITE (22,99001)
99001       FORMAT (2X////40X,'TROUBLESHOOTING ROUT',
     &              'INE HAS BEEN ACTIVATED...'//5X,
     &              'LOCAL MINIMUM ANALYSIS FOLLOWS:'//)
            WRITE (22,99002) ch/ll
99002       FORMAT (2X//5X,'CHISQ FOR FIRST GE(LI)S ONLY ',
     &              'WITH INDEPENDENT NORMALIZATION=',1E12.4//5X,
     &              'NORM.CONSTANTS:'//)
            DO i = 1 , NEXPT
               WRITE (22,99003) i , EP(i)
99003          FORMAT (5X,'EXP.',1X,1I2,5X,'C=',1E14.6)
            ENDDO
            WRITE (22,99004)
99004       FORMAT (1X//5X,'M.E.',20X,'RL',20X,'STRENGTH',//)
            DO i = 1 , MEMAX
               IF ( KVAR(i).NE.0 ) THEN
                  rl = LOG10(ELMU(i)/ABS(ELM(i)))
                  IF ( rl.GE.Rlr ) ELML(i) = 1.
                  WRITE (22,99005) i , rl , ELMU(i)/lngt
99005             FORMAT (6X,1I3,18X,1F4.1,20X,1E7.2)
               ENDIF
            ENDDO
            WRITE (22,99006)
99006       FORMAT (2X////40X,'ANALYSIS OF SIGNIFICANT DEPENDENCES'//)
            DO i = 1 , MEMAX ! For each matrix element
               IF ( KVAR(i).NE.0 ) THEN
                  lc = 0
                  IF ( ELML(i).GE..5 ) THEN
                     REWIND 17
 552                 READ (17,*) indx , iexp , sgm , ni , nf , al
                     IF ( indx.EQ.0 ) THEN
                        np = 0
                        nm = 0
                        DO j = 1 , lc
                           u = CORF(j,1)*CORF(j,2)*2.
                           IF ( ABS(u)/ELMU(i).GE..05 ) THEN
                              IF ( u.LT.0. ) nm = nm + 1
                              IF ( u.GT.0. ) np = np + 1
                           ENDIF
                        ENDDO
                        WRITE (22,99007) i , np , nm
99007                   FORMAT (1X/5X,10('*'),5X,'M.E.',1X,1I3,5X,1I3,
     &                          1X,'POSITIVE COMPONENTS',20X,1I3,1X,
     &                          'NEGATIVE COMPONENTS'///30X,'POSITIVE',
     &                          52X,'NEGATIVE'//5X,'EXP',2X,
     &                          'TRANSITION',2X,'SIGMA',3X,'DERIVATIVE',
     &                          3X,'D(SIGMA**2)/D(ME)',4X,'I',1X,'EXP',
     &                          2X,'TRANSITION',2X,'SIGMA',3X,
     &                          'DERIVATIVE',3X,'D(SIGMA**2)/D(ME)')
                        DO l = 1 , K ! For each of the important contributions to chisqr
                           ump = 0.
                           umm = 0.
                           DO j = 1 , lc
                              u = 2.*CORF(j,1)*CORF(j,2)
                              IF ( u.LT.0. ) THEN
                                 IF ( u.LE.umm ) THEN
                                    umm = u
                                    jm = j
                                 ENDIF
                              ELSEIF ( u.GE.ump ) THEN
                                 ump = u
                                 jp = j
                              ENDIF
                           ENDDO
                           WRITE (22,99008) IY(jp,1) , IY(jp,2) , 
     &                            IY(jp,3) , CORF(jp,1) , CORF(jp,2) , 
     &                            ump , IY(jm,1) , IY(jm,2) , IY(jm,3) ,
     &                            CORF(jm,1) , CORF(jm,2) , umm
99008                      FORMAT (5X,1I3,2X,1I3,'--',1I3,4X,1F4.1,4X,
     &                             1E9.2,7X,1E9.2,9X,'I',1X,1I3,2X,1I3,
     &                             '--',1I3,4X,1F4.1,4X,1E9.2,7X,1E9.2)
                           CORF(jp,1) = 0.
                           CORF(jm,1) = 0.
                        ENDDO
                     ELSE
                        IF ( indx.EQ.i ) THEN
                           lc = lc + 1
                           IY(lc,1) = iexp
                           IY(lc,2) = ni
                           IY(lc,3) = nf
                           CORF(lc,1) = sgm
                           CORF(lc,2) = al
                        ENDIF
                        GOTO 552
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO ! Loop over matrix elements
            RETURN
         ELSE
            ll = ll + 1
            ch = ch + sgm*sgm
            ux = 2.*sgm*u
            ELM(indx) = ELM(indx) + ux
            ELMU(indx) = ELMU(indx) + ABS(ux)
            GOTO 550
         ENDIF
      ELSE
         loc = (iexp-1)*LP2 + indx
         sgm = (e-a*EP(iexp))/b
         u = dy*EP(iexp)*b + a*ZETA(loc)/b
         WRITE (17,*) indx , iexp , sgm , ni , nf , u
         GOTO 500
      ENDIF
      END
 
C----------------------------------------------------------------------
C SUBROUTINE MIXR
C
C Called by: FTBM, GOSIA
C
C Purpose: calculate theoretical mixing ratio and compare to experimental one.
C
C Uses global variables:
C      DMIX   - 0.8326 * gamma energy
C      DMIXE  - mixing ratio and its error
C      IMIX   - decay associated with known mixing ratio
C      KSEQ   - index into ELM for pair of levels, and into EN or SPIN
C      LNY    - use logs to calculate chi squared
C      NDL    - number of mixing ratios
C
C Formal parameters:
C      Nw     - number of data points used to calculate chi squared
C      Ipsw   - printing flag (0 means no print, 1 means print)
C      Chi    - chi squared
C      Chilo  - chi squared using logs
 
      SUBROUTINE MIXR(Nw,Ipsw,Chi,Chilo)
      IMPLICIT NONE
      REAL*8 Chi , Chilo , dl
      INTEGER*4 i , inx , inx1 , Ipsw , it , Nw
      REAL*8 TAU
      INTEGER*4 KSEQ
      COMMON /LEV   / TAU(100) , KSEQ(1500,4)
      REAL*8 ELM , ELMU , ELML , SA
      COMMON /COMME / ELM(1500) , ELMU(1500) , ELML(1500) , SA(1500)
      REAL*8 DMIXE , DMIX
      INTEGER*4 IMIX , NDL
      COMMON /MIXD  / DMIXE(20,2) , DMIX(20) , IMIX(20) , NDL
      INTEGER*4 LNY , INTR , IPS1
      COMMON /LOGY  / LNY , INTR , IPS1
      dl = 0.
      
      IF ( NDL.EQ.0 ) RETURN
      Nw = Nw + NDL

      DO i = 1 , NDL ! For each mixing ratio
         it = IMIX(i) ! Decay for this mixing ratio
         inx = KSEQ(it,1) ! Index 1 of it'th decay
         inx1 = KSEQ(it,2) ! Index 2 of it'th decay
         IF ( ABS(ELM(inx1)).LT.1.E-5 ) ELM(inx1) = 1.E-5
         dl = DMIX(i)*ELM(inx)/ELM(inx1)
         IF ( Ipsw.EQ.1 ) DMIX(i) = dl
         Chi = Chi + (dl-DMIXE(i,1))**2/DMIXE(i,2)/DMIXE(i,2)
         IF ( LNY.EQ.1 ) Chilo = Chilo + 
     &                           (DMIXE(i,1)*LOG(ABS(dl/DMIXE(i,1)))
     &                           /DMIXE(i,2))**2
      ENDDO ! Loop on mixing ratios i

      IF ( Ipsw.EQ.0 ) RETURN
      
      WRITE (22,99001)
99001 FORMAT (1X//10X,'E2/M1 MIXING RATIOS'/10X,'TRANSITION',10X,
     &        'EXP.DELTA',10X,'CALC.DELTA',10X,'SIGMA'/)

      DO i = 1 , NDL ! For each mixing ratio
         dl = (DMIX(i)-DMIXE(i,1))/DMIXE(i,2) ! Relative error
         it = IMIX(i) ! Matrix element for this mixing ratio
         WRITE (22,99002) KSEQ(it,3) , KSEQ(it,4) , DMIXE(i,1) , DMIX(i)
     &                    , dl ! KSEQs are level numbers
99002    FORMAT (9X,1I3,'---',1I3,13X,1F7.2,12X,1F7.2,13X,1F5.2)
      ENDDO ! Loop on mixing ratios i

      END
 
C----------------------------------------------------------------------
C SUBROUTINE COORD
C
C Called by: GOSIA
C Calls:     TACOS, TASIN
C
C Purpose: calculate geometry for circular detector
C
C Uses global variables:
C      FIEX   - phi range of particle detector
C      ISKIN  - kinematic flag
C      IZ1    - Z of non-investigated nucleus
C      XA     - A of investigated nucleus
C      XA1    - A of non-investigated nucleus
C      YV     - scattering angle meshpoints where we calculate exact Coulex
C
C Formal parameters:
C      Wth    - theta of centre of detector (degrees) - readonly
C      Wph    - phi of centre of detector (degrees) - readonly
C      Wthh   - half angle subtended (degrees) - readonly
C      Naa    - number of theta divisions - readonly
C      Ifw    - flag: 0 for meshpoints, 1 for subdivisions, 2 for pin diodes -readonly
C      Pfi    - phi range for each theta value - writeonly
C      Wpi    - phi range of detector - read/write
C      Wtlb   - angle of particle detector in theta (degrees) in lab frame - readonly
C      Lz     - experiment number - readonly
C      Tyy    - lower limit of theta (degrees) - read/write
C      Tzz    - upper limit of theta (degrees) - read/write
 
      SUBROUTINE COORD(Wth,Wph,Wthh,Naa,Ifw,Pfi,Wpi,Wtlb,Lz,Tyy,Tzz)
      IMPLICIT NONE
      REAL*8 ga , gi , Pfi , rade , rmass , TACOS , TASIN , thetb , 
     &       ttcm , Tyy , Tzz
      REAL*8 wpa , Wph , Wpi , ws , Wth , Wthh , Wtlb , xaa , xph , 
     &       xth , xthh , za , za1 , zb , zl
      INTEGER*4 i , Ifw , Lz , Naa
      DIMENSION Pfi(101) , Wpi(100,2)
      REAL*8 XV, YV, ZV, DSG, DSE, DS
      COMMON /VLIN  / XV(101) , YV(101) , ZV(101) , DSG(101) ,
     &                DSE(101) , DS
      REAL*8 EPS, EROOT, FIEX
      INTEGER*4 IEXP, IAXS
      COMMON /KIN   / EPS(50) , EROOT(50) , FIEX(50,2) , IEXP , IAXS(50)
      REAL*8 XA , XA1 , EP , TLBDG , VINF
      INTEGER*4 NEXPT , IZ , IZ1
      COMMON /CX    / NEXPT , IZ , XA , IZ1(50) , XA1(50) , EP(50) , 
     &                TLBDG(50) , VINF(50)
      INTEGER*4 ISKIN
      COMMON /SECK  / ISKIN(50)
      DATA rade/57.2957795/ ! 180 / pi
      DATA ws/0./

      IF ( Ifw.EQ.0 ) THEN ! For meshpoints
         Tyy = Wth - Wthh ! Lower limit of theta
         Tzz = Wth + Wthh ! Upper limit of theta
      ENDIF

C     Convert to radians
      xth = Wth/rade ! theta of centre of detector in radians
      xph = Wph/rade ! phi of centre of detector in radians
      xthh = Wthh/rade ! half angle subtended in radians

C     pre-calculate trigonometric functions
      zl = TAN(xthh)
      za = COS(xth)
      za1 = SIN(xth)
      zb = COS(xthh)

      rmass = XA1(Lz)/XA ! Mass ratio for this experiment
      IF ( IZ1(Lz).LT.0 ) rmass = 1./rmass

C     Calculate size of each division (ws)
      IF ( Ifw.NE.2 ) THEN ! Unless we are using the pin diode option
         ws = (Tzz-Tyy)/(Naa+1)
         IF ( Ifw.EQ.1 ) ws = (Tzz-Tyy)/(Naa-1)
      ENDIF

      DO i = 1 , Naa ! Loop over theta divisions
         IF ( Ifw.NE.2 ) THEN ! Not pin diode option
            IF ( Ifw.EQ.0 ) YV(i) = Tyy + i*ws ! theta value for this step in degrees
            xaa = (Tyy+ws*(i-1))/rade ! and in radians
            IF ( Ifw.EQ.1 .AND. (i.EQ.1 .OR. i.EQ.Naa) ) THEN
               Pfi(i) = 0.
               GOTO 100
            ELSE
               IF ( Ifw.EQ.0 ) xaa = YV(i)/rade
            ENDIF
         ELSE ! Pin diode option
            xaa = ABS(Wtlb)/rade ! Detector angle theta in lab frame in radians
            IF ( Wtlb.GT.0. ) GOTO 50
            IF ( IZ1(Lz).LT.0 ) THEN
               IF ( XA.LE.XA1(Lz) ) GOTO 20
            ELSEIF ( XA1(Lz).LE.XA ) THEN
               GOTO 20
            ENDIF
            IF ( ISKIN(Lz).EQ.0 ) THEN ! ISKIN = 0 means take lower CM angle
               ttcm = xaa - TASIN(rmass*SIN(xaa))
               xaa = ABS(ttcm)/2.
               GOTO 50
            ENDIF
 20         ttcm = xaa + TASIN(rmass*SIN(xaa)) ! Take higher CM angle
            xaa = (3.14159265-ttcm)/2.
         ENDIF ! End of pin diode option

 50      gi = (za-COS(xaa)/zb)/(zl*za1)
         ga = TACOS(gi)
         wpa = ATAN(zl*SIN(ga)/(za1+zl*COS(ga)*za))
         wpa = ABS(wpa)
         IF ( Ifw.EQ.2 ) THEN ! Pin diode option
            FIEX(Lz,1) = (xph-wpa) ! phi min
            FIEX(Lz,2) = (xph+wpa) ! phi max
         ELSEIF ( Ifw.EQ.1 ) THEN ! Interpolation option
            Pfi(i) = 2.*wpa*rade
         ELSE ! Meshpoint option
            Wpi(i,1) = (xph-wpa)*rade ! Lower phi limit
            Wpi(i,2) = (xph+wpa)*rade ! Upper phi limit
         ENDIF
 100     CONTINUE
      ENDDO ! Loop on theta divisions i

C     If a negative value of theta was specified for a meshpoint value,
C     we use the target angle
      IF ( Wtlb.LT.0. .AND. Ifw.EQ.0 ) THEN
         DO i = 1 , Naa ! For each theta division
            xaa = YV(i)/rade ! theta in radians
            thetb = ATAN(SIN(2.*xaa)/(rmass-COS(2.*xaa)))*rade
            IF ( thetb.LT.0. ) thetb = 180. + thetb
            YV(i) = -1.*thetb
            Wpi(i,1) = Wpi(i,1) + 180.
            Wpi(i,2) = Wpi(i,2) + 180.
         ENDDO
      ENDIF
      END
 
C----------------------------------------------------------------------
C SUBROUTINE CHMEM
C
C Called by: FTBM
C
C Purpose: compare fitted matrix elements with known ones and calculate the
C          effect on the chi squared
C
C Uses global variables:
C      ELM    - matrix elements
C      EAMX   - known matrix elements and their error
C      NAMX   - number of known matrix elements
C      IAMX   - index of matrix element for known matrix element
C      IAMY   - level indices of pair of levels for which matrix element is known
 
      SUBROUTINE CHMEM(Nw,Chi,Chilo)
      IMPLICIT NONE
      REAL*8 Chi , Chilo , di
      INTEGER*4 ia , ib , Nw
      REAL*8 EAMX
      INTEGER*4 NAMX, IAMX, IAMY
      COMMON /ME2D  / EAMX(100,2) , NAMX , IAMX(100) , IAMY(100,2)
      REAL*8 ELM , ELMU , ELML , SA
      COMMON /COMME / ELM(1500) , ELMU(1500) , ELML(1500) , SA(1500)
      IF ( NAMX.EQ.0 ) RETURN
      Nw = Nw + NAMX
      DO ia = 1 , NAMX
         ib = IAMX(ia)
         IF ( IAMY(ia,1).NE.IAMY(ia,2) ) THEN
            di = (ELM(ib)-EAMX(ia,1))/EAMX(ia,2)
            Chilo = Chilo + 
     &              (LOG(ABS(ELM(ib)/EAMX(ia,1)))*ABS(EAMX(ia,1))
     &              /EAMX(ia,2))**2
            Chi = Chi + di*di
         ELSE
            di = (ELM(ib)-EAMX(ia,1))/EAMX(ia,2)
            Chilo = Chilo + 
     &              (LOG(ABS(ELM(ib)/EAMX(ia,1)))*ABS(EAMX(ia,1))
     &              /EAMX(ia,2))**2
            Chi = Chi + di*di
         ENDIF
      ENDDO
      END
 
C----------------------------------------------------------------------
C SUBROUTINE PTICC
C
C Called by: GOSIA
C Calls:     CONV
C
C Purpose: print the conversion coefficients
C
C Uses global variables:
C      EN     - energy of level
C      KSEQ   - index into ELM for pair of levels, and into EN or SPIN
C      MULTI  - number of matrix elements having a given multipolarity
C      SPIN   - spin of level
C
C Formal parameters:
C      Idr    - number of decays

      SUBROUTINE PTICC(Idr)
      IMPLICIT NONE
      REAL*8 cone1 , cone2 , conm1 , CONV , enet
      INTEGER*4 Idr , iinx , l , nf , ni
      INTEGER*4 LAMDA , LEAD , LDNUM , LAMMAX , MULTI
      COMMON /CLCOM / LAMDA(8) , LEAD(2,1500) , LDNUM(8,100) , LAMMAX , 
     &                MULTI(8)
      REAL*8 EN , SPIN , ACCUR , DIPOL , ZPOL , ACCA
      INTEGER*4 ISO
      COMMON /COEX  / EN(100) , SPIN(100) , ACCUR , DIPOL , ZPOL , 
     &                ACCA , ISO
      REAL*8 TAU
      INTEGER*4 KSEQ
      COMMON /LEV   / TAU(100) , KSEQ(1500,4)
      WRITE (22,99001)
99001 FORMAT (1X//20X,'CALCULATED INTERNAL CONVERSION ',
     &        'COEFFICIENTS FOR E1,E2 AND M1'//5X,'NI',5X,'NF',7X,'II',
     &        8X,'IF',9X,'ENERGY(MEV)',6X,'ICC(E1)',8X,'ICC(E2)',8X,
     &        'ICC(M1)')
      DO l = 1 , Idr
         iinx = KSEQ(l,1) ! Index of l'th decay
         ni = KSEQ(l,3) ! Initial level of l'th decay
         nf = KSEQ(l,4) ! Final level of l'th decay
         enet = EN(ni) - EN(nf)
         cone2 = CONV(enet,2)
         IF ( ABS(SPIN(ni)-SPIN(nf)).GT.2. ) cone2 = 0.
         conm1 = 0.
         cone1 = 0.
         IF ( iinx.LE.MULTI(1) ) cone1 = CONV(enet,1)
         IF ( ABS(SPIN(ni)-SPIN(nf)).LT.2. ) conm1 = CONV(enet,4)
         WRITE (22,99002) ni , nf , SPIN(ni) , SPIN(nf) , enet , cone1 ,
     &                    cone2 , conm1
99002    FORMAT (4X,I3,4X,I3,7X,F4.1,6X,F4.1,9X,F6.4,8X,E9.4,6X,E9.4,6X,
     &           E9.4)
      ENDDO
      END
 
C----------------------------------------------------------------------
C FUNCTION RNDM
C
C Called by: MIXUP
C
C Purpose: Generate a pseudo-random number based on the seed Se
C
C Formal parameters:
C      Se     - seed for random number
C
C It is used to generate random matrix elements as a starting position,
C when OP,RAND is called. The parameter to OP,RAND is the seed here.
 
      REAL*8 FUNCTION RNDM(Se)
      IMPLICIT NONE
      REAL*8 ai , p , r , rxdm , Se , t , u
      INTEGER*4 i
      DATA t/0./
      SAVE t

      IF ( Se.GT.32000. ) Se = 100.*t + .511
      Se = Se*Se
      u = LOG10(Se)
      i = INT(u) + 1
      t = Se/(10.**i)
      r = SQRT(SQRT(SQRT(t)))
      p = SQRT(SQRT(SQRT(.1)))
      rxdm = (r-p)/(1.-p)
      rxdm = 10.*rxdm
      ai = DBLE(INT(rxdm))
      RNDM = rxdm - ai
      END
 
C----------------------------------------------------------------------
C SUBROUTINE KONTUR
C
C Called by: GOSIA
C Calls:     FTBM, LIMITS, RK4
C
C Purpose: scans the chi^2 hypersurface in OP,ERRO
C
C Uses global variables:
C      DEVU   -
C      ELM    - matrix elements
C      ELML   - lower limit on matrix elements
C      ELMU   - upper limit on matrix elements
C      HLM    - previous values of matrix elements
C      INTR   - flag to swap chisqr and log(chisqr)
C      IPS1   - terminate after calculating and writing correction factors
C      LNY    - use logs to calculate chi squared
C      MEMAX  - number of matrix elements
C      NWR    - number of datapoints used in fit
C      SA     - ratio of elements for correlated elements
C      XV     - energy meshpoints where we calculate exact Coulex
C      YV     - scattering angle meshpoints where we calculate exact Coulex
C
C Formal parameters:
C      Idr    - number of decays
C      Chis0  -
C      Chil   -
C      Ifbf   -
C      Inpo   -
C      Jj     - matrix element
C      Sh     -
C      Bten   -
C      Rem    - natural log of the largest value the computer can represent
 
      SUBROUTINE KONTUR(Idr,Chis0,Chil,Ifbf,Inpo,Jj,Sh,Bten,Rem)
      IMPLICIT NONE
      REAL*8 ac , Bten , c , Chil , chilo , Chis0 , chis1 , chis2 , d1 ,
     &       d2 , f , h
      REAL*8 Rem , RK4 , sajj , Sh , t , v , ww , x , y
      INTEGER*4 i , Idr , Ifbf , Inpo , itl , ix , j , Jj , l , m
      DIMENSION f(3) , Bten(*)
      REAL*8 XV, YV, ZV, DSG, DSE, DS
      COMMON /VLIN  / XV(101) , YV(101) , ZV(101) , DSG(101) ,
     &                DSE(101) , DS
      REAL*8 DEVD, DEVU
      COMMON /DFTB  / DEVD(1500) , DEVU(1500)
      REAL*8 ELM , ELMU , ELML , SA
      COMMON /COMME / ELM(1500) , ELMU(1500) , ELML(1500) , SA(1500)
      INTEGER*4 MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR
      COMMON /CEXC  / MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR(1500)
      REAL*8 HLM
      COMMON /HHH   / HLM(1500)
      INTEGER*4 NWR
      COMMON /ILEWY / NWR
      INTEGER*4 LNY , INTR , IPS1
      COMMON /LOGY  / LNY , INTR , IPS1
      LNY = 0
      h = .05*ABS(HLM(Jj))
      IF ( Inpo.NE.-1 ) h = ABS(Sh)
 100  INTR = 0
      sajj = ABS(SA(Jj)) ! ratio of matrix elements for correlation
      DO l = 1 , MEMAX ! For each matrix element
         ELM(l) = HLM(l)
         SA(l) = SA(l)/sajj
      ENDDO
      YV(1) = 0.
      XV(1) = HLM(Jj)
      f(3) = 1.
      i = 1
 200  itl = 0
      v = ELMU(Jj) - ELM(Jj)
      IF ( SA(Jj).LT.0. ) v = ELM(Jj) - ELML(Jj)
      IF ( h.GT.v ) itl = 1
      IF ( h.GT.v ) h = v
      i = i + 1
      f(1) = f(3)
      DO j = 1 , MEMAX
         ELM(j) = .5*h*SA(j) + ELM(j)
      ENDDO
      CALL LIMITS ! Constrain matrix elements within limits
      CALL FTBM(3,chis1,Idr,1,chilo,Bten)
      IF ( chis1.LE.Chis0 ) THEN
         IF ( Inpo.EQ.-1 ) WRITE (22,99003) Jj , ELM(Jj) , chis1
         IF ( chis1.LE.Chil .AND. Inpo.NE.-1 ) THEN
            Ifbf = 1
            ix = 1
            Chil = chis1
            WRITE (22,99004) Chil
            GOTO 500
         ENDIF
      ENDIF
 300  ww = .5*(Chis0-chis1)*NWR
      IF ( ww.GE.Rem ) GOTO 700
      f(2) = EXP(ww)
      IF ( i.EQ.2 .AND. f(2).LT..1 .AND. ABS(XV(1)-HLM(Jj)).LT.1E-9 )
     &     THEN
         h = h/2.
         GOTO 100
      ELSE
         DO j = 1 , MEMAX ! For each matrix element
            ELM(j) = ELM(j) + .5*SA(j)*h
         ENDDO
         v = ELM(Jj)
         CALL LIMITS ! Constrain matrix elements within limits
         IF ( ABS(v-ELM(Jj)).GT.1.E-6 ) itl = 1
         CALL FTBM(3,chis2,Idr,1,chilo,Bten)
         IF ( chis2.LE.Chis0 ) THEN
            IF ( Inpo.EQ.-1 ) WRITE (22,99003) Jj , ELM(Jj) , chis2
            IF ( chis2.LE.Chil .AND. Inpo.NE.-1 ) THEN
               Ifbf = 1
               ix = 2
               Chil = chis2
               WRITE (22,99004) Chil
               GOTO 500
            ENDIF
         ENDIF
      ENDIF
 400  ww = .5*(Chis0-chis2)*NWR
      IF ( ww.GT.Rem ) GOTO 700
      f(3) = EXP(ww)
      IF ( itl.EQ.1 ) WRITE (22,99001) Jj
99001 FORMAT (5X,'WARNING-ME(',1I3,')',5X,
     &        'INTEGRATION STOPPED AT THE LIMIT')
      IF ( i.EQ.2 ) THEN
         IF ( itl.NE.1 ) THEN
            IF ( f(3).LT..1 .AND. ABS(XV(1)-HLM(Jj)).LT.1.E-9 ) THEN
               h = h/2.
               GOTO 100
            ELSEIF ( f(1).LE.f(2) .OR. f(2).LE.f(3) ) THEN
               IF ( f(1).LT.f(2) .AND. f(2).GT.f(3) ) THEN
                  d1 = f(2) - f(1)
                  d2 = f(3) - f(1)
                  ac = (d2-4.*d1)*h/(d2-2.*d1)/4.
                  DO l = 1 , MEMAX
                     ELM(l) = (ELM(l)-h*SA(l)) + ac*SA(l)
                  ENDDO
                  CALL LIMITS
                  XV(1) = ELM(Jj)
                  i = 1
                  CALL FTBM(3,chis1,Idr,1,chilo,Bten)
                  ww = .5*(Chis0-chis1)*NWR
                  IF ( ww.GE.Rem ) GOTO 700
                  f(3) = EXP(ww)
                  GOTO 200
               ELSE
                  i = 1
                  XV(1) = ELM(Jj)
                  IF ( Inpo.EQ.-1 ) h = 2.*h
                  GOTO 200
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      y = YV(i-1)
      YV(i) = RK4(y,h,f)
      XV(i) = ELM(Jj)
      IF ( NWR*(chis2-Chis0).LT.2. .AND. Inpo.EQ.-1 ) h = 2.*h
      IF ( itl.EQ.1 ) GOTO 600
      IF ( f(3).GE.1.E-3 ) GOTO 200
      GOTO 600
 500  REWIND 17
      DO l = 1 , MEMAX ! For each matrix element
         WRITE (17,*) ELM(l)
      ENDDO
      IF ( ix.EQ.1 ) GOTO 300
      IF ( ix.NE.2 ) GOTO 200
      GOTO 400
 600  c = YV(i)
      m = 0
      DO l = 1 , i
         YV(l) = 1.00001 - YV(l)/c
         IF ( m.EQ.0 .AND. YV(l).LT..317 ) m = l
      ENDDO
      x = (XV(m)-XV(m-1))*(.317-YV(m))/(YV(m-1)-YV(m))
      t = XV(m) - x - HLM(Jj)
      IF ( t.GE.0. ) DEVU(Jj) = t
      IF ( t.LT.0. ) DEVD(Jj) = t
      RETURN
 700  WRITE (22,99002) Jj
99002 FORMAT (5X,'** WARNING **',/,2X,'ME=',1I3,2X,
     &     'TOO FAR FROM THE MINIMUM TO CARRY OUT THE ERROR ESTIMATION!'
     &     ,/)
99003 FORMAT (5X,'ELM(',1I3,')=',1F10.6,5X,'CHISQ=',1E12.4)
99004 FORMAT (10X,'BETTER POINT FOUND...MATRIX ELEMENTS WRITTEN ON 17',
     &        3X,'CHISQ=',1E12.4)
      END
 
C----------------------------------------------------------------------
C FUNCTION RK4
C
C Called by: KONTUR
C
C Purpose:
C
C Formal parameters:
C      Y      - 
C      H      -
C      F      - array of three coefficients
C
C Return value:
 
      REAL*8 FUNCTION RK4(Y,H,F)
      IMPLICIT NONE
      REAL*8 F , H , Y
      DIMENSION F(3)

      RK4 = Y + H*(F(1)+4.*F(2)+F(3))/6.
      END
 
C----------------------------------------------------------------------
C SUBROUTINE QFIT
C
C Called by: GOSIA
C Calls:     GAMATT
C
C Purpose: for OP,GDET, fit attenuation by absorbers
C
C Formal parameters:
C      Qui    - attenuation coefficients
C      Tau1   - absorption coefficients' table
C      Tau2   - absorption coefficients' table
C      Eng    - gamma energy
C      Xl1    - thickness of absorbers
C      Cf     - coefficients of fit
C      Nl     - number of types of absorber (7)
C      Ind    - type of absorber
C
C Note the absorbers are: Al, C, Fe, Cu, Ag/Cd/Sn, Ta and Pb, respectively.
      
      SUBROUTINE QFIT(Qui,Tau1,Tau2,Eng,Xl1,Cf,Nl,Ind)
      IMPLICIT NONE
      REAL*8 ca , cb , Cf , cm , cn , co , d , d1 , d2 , Eng , Qui , 
     &       Tau1 , Tau2 , Xl1
      INTEGER*4 Ind , ind1 , k , Nl
      DIMENSION Tau1(10) , Eng(10) , Tau2(10,7) , Xl1(7) , Qui(8,10) , 
     &          Cf(8,2)

      CALL GAMATT(Qui,Tau1,Tau2,Xl1,Nl)
      
      ind1 = 5
      IF ( Ind.EQ.4 ) ind1 = 6
      IF ( Ind.EQ.5 ) ind1 = 7
      DO k = 1 , 8
         co = Qui(k,Ind)
         cn = Qui(k,10)
         cm = Qui(k,ind1)
         ca = (Eng(ind1)-Eng(Ind))**2
         cb = (Eng(10)-Eng(Ind))**2
         d = ca*(co-cn) - cb*(co-cm)
         d1 = ca*cm*(co-cn) - cb*cn*(co-cm)
         d2 = ca*cb*(cn-cm)
         Cf(k,1) = d1/d
         Cf(k,2) = d2/d
      ENDDO
      END
 
C----------------------------------------------------------------------
C SUBROUTINE GAMATT
C
C Called by: QFIT
C Calls:     GCF
C
C Purpose: calculate gamma attenuation in absorbers
C
C Formal parameters:
C      Qui    - attenuation
C      Tau1   - table of absorption coefficients
C      Tau2   - table of absorption coefficients
C      Xl1    - thickness of each kind of absorber
C      Nl     - number of kinds of absorber, we can treat (7)
C
C Note the absorbers are: Al, C, Fe, Cu, Ag/Cd/Sn, Ta and Pb, respectively.
 
      SUBROUTINE GAMATT(Qui,Tau1,Tau2,Xl1,Nl)
      IMPLICIT NONE
      INTEGER*4 i , i1 , k , Nl
      REAL*8 q , Qui , tau , Tau1 , Tau2 , thing , thing1 , thing3 , Xl1
      DIMENSION Tau1(10) , Tau2(10,7) , Xl1(7) , thing3(10) , q(9) , 
     &          Qui(8,10)

C     Treat absorbers
      DO i = 1 , 10 ! Loop over energies
         i1 = 1
         thing3(i) = 0.
 50      thing1 = -Tau2(i,i1)*Xl1(i1) + thing3(i)
         i1 = i1 + 1
         thing3(i) = thing1
         IF ( i1.LE.Nl ) GOTO 50 ! Loop over Nl absorbers
      ENDDO

C     Treat germanium
      DO i = 1 , 10 ! Loop over energies
         tau = Tau1(i)
         thing = thing3(i)
         CALL GCF(tau,thing,q)
         DO k = 2 , 9
            Qui(k-1,i) = q(k)
         ENDDO
      ENDDO
      END
 
C----------------------------------------------------------------------
C SUBROUTINE GCF
C
C Called by: GAMATT
C
C Purpose: calculate detection probability (probability that gamma of a given
C          energy is absorbed in the Ge but not in one of the absorbers).
C
C Formal parameters:
C      Tau    - Absorption coefficient for Ge at this energy (input)
C      Thing  - Absorption coefficient for absorbers at this energy (input)
C      Q      - Attenuation coefficient (output)
 
      SUBROUTINE GCF(Tau,Thing,Q)
      IMPLICIT NONE
      REAL*8 A , b , D , dl , ev , ex , f , fint , od , Q , R , 
     &       Tau , Thing , XL , xm , yl , yu
      INTEGER*4 i , j , k , m
      REAL*8 DIX, ODL
      COMMON /DIMX  / DIX(4) , ODL(200)
      DIMENSION f(101) , b(4) , Q(9)

      A = DIX(1)
      R = DIX(2)
      XL = DIX(3)
      D = DIX(4)
      
      b(1) = ATAN2(A,D+XL)
      b(2) = ATAN2(A,D)
      b(3) = ATAN2(R,D+XL)
      b(4) = ATAN2(R,D)
      DO k = 1 , 9 ! Loop over order of Legendre polynomial order
         Q(k) = 0.0
         DO j = 1 , 3
            yl = b(j)
            yu = b(j+1)
            dl = (yu-yl)/100.
            DO m = 1 , 101
               xm = yl + dl*(m-1)
               IF ( j.EQ.2 ) THEN
                  ex = -Tau*XL/COS(xm)
               ELSEIF ( j.EQ.3 ) THEN
                  ex = Tau*(D*TAN(xm)-R)/SIN(xm)
               ELSE
                  ex = Tau*(A-(D+XL)*TAN(xm))/SIN(xm)
               ENDIF
               f(m) = SIN(xm)*(1-EXP(ex))*EXP(Thing/COS(xm))
               IF ( j.EQ.1 ) f(m) = f(m)*EXP(-Tau*(A/SIN(xm)-D/COS(xm)))
               IF ( k.EQ.1 ) THEN ! Legendre polynomials order k
               ELSEIF ( k.EQ.3 ) THEN
                  f(m) = f(m)*(1.5*COS(xm)**2-0.5)
               ELSEIF ( k.EQ.4 ) THEN
                  f(m) = f(m)*(2.5*COS(xm)**3-1.5*COS(xm))
               ELSEIF ( k.EQ.5 ) THEN
                  f(m) = f(m)*(4.375*COS(xm)**4-3.75*COS(xm)**2+.375)
               ELSEIF ( k.EQ.6 ) THEN
                  f(m) = f(m)*((63.*COS(xm)**5-70.*COS(xm)**3+15.)/8.)
               ELSEIF ( k.EQ.7 ) THEN
                  f(m) = f(m)
     &                   *((21.*COS(xm)**2*(11.*COS(xm)**4-15.*COS(xm)
     &                   **2+5.)-5.)/16.)
               ELSEIF ( k.EQ.8 ) THEN
                  f(m) = f(m)
     &                   *(429.*COS(xm)**7-693.*COS(xm)**5+315.*COS(xm)
     &                   **3-35.*COS(xm))/16.
               ELSEIF ( k.EQ.9 ) THEN
                  f(m) = f(m)
     &                   *(6435.*COS(xm)**8-12012.*COS(xm)**6+6930.*COS
     &                   (xm)**4-1260.*COS(xm)**2+35.)/128.
               ELSE
                  f(m) = f(m)*COS(xm)
               ENDIF
            ENDDO
            ev = 0.0
            od = 0.0
            DO m = 2 , 98 , 2
               ev = ev + f(m)
               od = od + f(m+1)
            ENDDO
            fint = dl/3.*(f(1)+4.*(ev+f(100))+2.*od+f(101))
            Q(k) = Q(k) + fint
         ENDDO
      ENDDO
      DO i = 1 , 8
         Q(i+1) = Q(i+1)/Q(1)
      ENDDO
      Q(1) = Q(1)/2.
      END
 
C----------------------------------------------------------------------
C FUNCTION TCEXP
C
C Called by: EXPON
C
C Purpose: evaluates a complex exponential
C
C Formal parameters:
C      Z      - argument of exponential (complex)
C
C Return value:
C      complex exponential of Z.
 
      COMPLEX*16 FUNCTION TCEXP(Z)
      IMPLICIT NONE
      REAL*8 a , b , c , d
      COMPLEX*16 Z
      
      a = DBLE(Z)
      b = DIMAG(Z)
      a = EXP(a)
      c = a*COS(b)
      d = a*SIN(b)
      TCEXP = DCMPLX(c,d)
      END
 
C----------------------------------------------------------------------
C FUNCTION TCABS
C
C Called by: laiamp, pomnoz
C
C Purpose: evaluates the absolute value of a complex number
C
C Formal parameters:
C      Z      - argument for abs
C
C Return value:
C      absolute value of complex number Z.
 
      REAL*8 FUNCTION TCABS(Z)
      IMPLICIT NONE
      REAL*8 a , b
      COMPLEX*16 Z
      
      a = DBLE(Z)
      b = DIMAG(Z)
      IF ( ABS(a).LT.1.E-16 ) a = 0.
      IF ( ABS(b).LT.1.E-16 ) b = 0.
      TCABS = SQRT(a*a+b*b)
      END
 
C----------------------------------------------------------------------
C FUNCTION TASIN
C
C Called by: CMLAB, COORD, TACOS
C
C Purpose: calculate an arcsine(x)
C
C Formal parameters:
C      X      - value for which we are to evaluate the arcsine
C
C Return value:
C      arcsine of X
C
C We take care of the special case of abs(x) = 1. Otherwise, we evaluate
C arctan(x / sqrt(1 - x^2).

      REAL*8 FUNCTION TASIN(X)
      IMPLICIT NONE
      REAL*8 dol , test , war , X
      
      test = ABS(X) - 1.
      IF ( ABS(test).LT.1.E-9 ) THEN
         TASIN = 1.570796327 ! 1.570796327 is pi / 2
         IF ( X.LT.0. ) TASIN = -1.570796327
         GOTO 99999
      ENDIF
      dol = SQRT(1.-X*X)
      war = X/dol
      TASIN = ATAN(war)
99999 END
 
C----------------------------------------------------------------------
C FUNCTION TACOS
C
C Called by: ARCCOS, CEGRY, COORD, GOSIA
C Calls:     TASIN
C
C Purpose: evaluate arccosine(x)
C
C Formal parameters:
C      X      - value for which we are to evaluate the arccosine
C
C Return value:
C      arccosine of X
C
C We use: arccos(x) = pi/2 - arcsin(x)
 
      REAL*8 FUNCTION TACOS(X)
      IMPLICIT NONE
      REAL*8 TASIN , X
      
      TACOS = 1.570796327 - TASIN(X) ! 1.570796327 = pi / 2
      END
 
C----------------------------------------------------------------------
C SUBROUTINE OPENF
C
C Called by: GOSIA
C
C Purpose: open files to specified units.
C
C Uses global variables:
C      JZB    - unit to read from
C
C The function reads three integers, the first of which is the unit to use for
C the open statement. The second is 1 if the file is required to exist already,
C 2 if it is required not to exist and 3 if it does not matter. The third is 1
C if the file is formatted and 2 if it is unformatted. A second line is read,
C which gives the name of the file to associate with that unit. If the unit is
C zero, the function returns. It keeps looping until a unit zero is reached.
 
      SUBROUTINE OPENF
      IMPLICIT NONE
      INTEGER*4 i , j , k
      CHARACTER name*60 , opt1*20 , opt2*20
      INTEGER*4 IUNIT3 , JZB
      COMMON /SWITCH/ JZB , IUNIT3
 100  READ (JZB,*) i , j , k ! unit, old/new/unknown, formatted/unformatted
      IF ( i.EQ.0 ) RETURN
      IF ( j.EQ.1 ) opt1 = 'OLD'
      IF ( j.EQ.2 ) opt1 = 'NEW'
      IF ( j.EQ.3 ) opt1 = 'UNKNOWN'
      IF ( k.EQ.1 ) opt2 = 'FORMATTED'
      IF ( k.EQ.2 ) opt2 = 'UNFORMATTED'
      READ (JZB,99001) name ! name of file
99001 FORMAT (A)

C     If it is for unit 25 or 26 and we are not reading from unit 5, ignore it
      IF ( JZB.NE.5 .AND. (i.EQ.25 .OR. i.EQ.26) ) GOTO 100 ! For gosia2

C     Now open the file
      OPEN (i,IOSTAT=k,FILE=name,STATUS=opt1,FORM=opt2)
      IF ( k.EQ.0 ) WRITE (6,99002) 'OPENED ' , name
99002 FORMAT (1X,2A)
      WRITE (6,99003) ' IO-num = ' , i , opt1 , opt2
99003 FORMAT (1X,A,I4,2(1x,A))
      IF ( k.EQ.0 ) GOTO 100
      WRITE (6,99004) 'PROBLEMS OPENING ' , name , k
99004 FORMAT (A,A,I6)
      END
 
C----------------------------------------------------------------------
C SUBROUTINE EFFIX
C
C Called by: CEGRY, GOSIA
C Calls:     LAGRAN, SPLNER
C
C Purpose: calculate the efficiency of the detector at a given energy.
C
C Uses global variables:
C      ABC    - absorption coefficients
C      AKAVKA - efficiency curve parameters
C      THICK  - thickness of each absorber type
C
C Formal parameters:
C      Ipd    - detector number
C      En     - gamma-ray energy
C      Effi   - efficiency
C
C Note that it uses LAGRAN or SPLNER according to the ISPL flag to
C interpolate between the data points given by the user.
C
C The efficiency curve parameters are those of GREMLIN plus an extra control
C flag:
C     AKAVKA(1) = a0
C     AKAVKA(2) = a1
C     AKAVKA(3) = a2
C     AKAVKA(4) = a3
C     AKAVKA(5) = f - for F-factor
C     AKAVKA(6) = N - for F-factor
C     AKAVKA(7) = b - for Woods-saxon factor
C     AKAVKA(8) = c - for Woods-saxon factor
C     AKAVKA(9) = control flag
C
C Efficiency parametrizations (control flag):
C     0  - Gremlin
C     1  - Jaeri
C     2  - Fiteff
C     3  - Leuven
C     4  - Radware
      
      SUBROUTINE EFFIX(Ipd,En,Effi)
      IMPLICIT NONE
      REAL*8 d , Effi , En , enl , pw , s , t , w , xx , yy
      INTEGER*4 i , Ipd , j , l , ll , n
      DIMENSION xx(101) , yy(101)
      REAL*8 ABC, AKAVKA, THICK
      COMMON /EFCAL / ABC(8,10) , AKAVKA(9,200) , THICK(200,7)
      REAL*8 EG, CC, AGELI, Q
      INTEGER*4 NICC, NANG, ISPL
      COMMON /CCC   / EG(50) , CC(50,5) , AGELI(50,200,2) , Q(3,200,8) ,
     &                NICC , NANG(200) , ISPL
      
      Effi = 1.E-6
      En = En + 1.E-24
      enl = LOG(En)
      DO i = 1 , 10
         ll = 11 - i
         j = ll
         IF ( enl.GE.ABC(8,ll) ) GOTO 100
         j = -1
      ENDDO
 100  IF ( j.EQ.-1 ) Effi = 1.E-10
      IF ( j.EQ.-1 ) RETURN
      IF ( j.EQ.1 .OR. j.EQ.10 ) THEN
         s = 0.
         DO l = 1 , 7
            IF ( ABS(THICK(Ipd,l)).GE.1.E-9 ) THEN
               t = EXP(ABC(l,j))
               d = THICK(Ipd,l)
               s = s + t*d
            ENDIF
         ENDDO
      ELSE
         IF ( j.EQ.9 ) THEN
            xx(1) = ABC(8,8)
            xx(2) = ABC(8,9)
            xx(3) = ABC(8,10)
         ELSE
            xx(1) = ABC(8,j)
            xx(2) = ABC(8,j+1)
            xx(3) = ABC(8,j+2)
         ENDIF
         s = 0.
         DO l = 1 , 7
            IF ( ABS(THICK(Ipd,l)).GE.1.E-9 ) THEN
               IF ( j.EQ.9 ) THEN
                  yy(1) = ABC(l,8)
                  yy(2) = ABC(l,9)
                  yy(3) = ABC(l,10)
               ELSE
                  yy(1) = ABC(l,j)
                  yy(2) = ABC(l,j+1)
                  yy(3) = ABC(l,j+2)
               ENDIF
               IF ( ISPL.EQ.0 ) CALL LAGRAN(xx,yy,3,0,enl,t,1,1)
               IF ( ISPL.EQ.1 ) CALL SPLNER(xx,yy,3,enl,t,1)
               s = s + EXP(t)*THICK(Ipd,l)
            ENDIF
         ENDDO
      ENDIF
      Effi = EXP(-s)

C     Branch according to type of calibration
      IF ( (AKAVKA(8,Ipd).LE.-999.) .OR. (AKAVKA(9,Ipd).EQ.3.) ) THEN
         GOTO 1003 ! Leuven
      ELSEIF ( AKAVKA(9,Ipd).EQ.4. ) THEN
         GOTO 1004 ! Radware
      ELSEIF ( (AKAVKA(5,Ipd).GT.0. .AND. AKAVKA(5,Ipd).LT.10.) .OR. 
     &         (AKAVKA(9,Ipd).EQ.2.) ) THEN
         GOTO 1002 ! Fiteff
      ELSEIF ( (AKAVKA(5,Ipd).LT.10.) .AND. (AKAVKA(9,Ipd).NE.1.) ) THEN
         GOTO 1000 ! Gremlin
      ENDIF
      GOTO 1001 ! Jaeri

C-----------------------------------------------------------------
C     GREMLIN efficiency calibration
 1000 w = LOG(20.*En) ! E0 = 50 keV, so w = LOG(En/E0) with En in MeV
      pw = AKAVKA(1,Ipd) + AKAVKA(2,Ipd)*w + AKAVKA(3,Ipd)
     &     *w*w + AKAVKA(4,Ipd)*w*w*w
      Effi = Effi*EXP(pw)
      IF ( ABS(AKAVKA(5,Ipd)).GE.1.E-9 ) THEN ! F-factor
         n = INT(AKAVKA(6,Ipd)+.1)
         pw = w**n
         w = AKAVKA(5,Ipd)/pw
         Effi = Effi*EXP(w)
      ENDIF
      IF ( ABS(AKAVKA(8,Ipd)).LT.1.E-9 ) RETURN
      w = (AKAVKA(7,Ipd)-1000.*En)/AKAVKA(8,Ipd) ! Woods-saxon factor
      pw = EXP(w)
      IF ( ABS(pw-1.).LT.1.E-6 ) WRITE (22,99001)
99001 FORMAT (5x,'***** CRASH - EFFIX *****')
      Effi = Effi/(1.+pw) ! Older versions of gosia have a minus sign here, which is wrong
                          ! because it is not what is done in gremlin (FITFUN) or the gosia manual
      RETURN

C-----------------------------------------------------------------
C     JAERI efficiency calibration - TC, Nov.2000
 1001 w = LOG(En/.511)
      Effi = EXP(AKAVKA(1,Ipd)+AKAVKA(2,Ipd)
     &       *w-EXP(AKAVKA(3,Ipd)+AKAVKA(4,Ipd)*w))
      RETURN

C-----------------------------------------------------------------
C     FITEFF efficiency calibration by P.Olbratowski use
C     PJN@2000
 1002 w = LOG(En/AKAVKA(5,Ipd))
      pw = AKAVKA(2,Ipd)*w
      IF ( En.LT.AKAVKA(5,Ipd) ) pw = pw + 
     &     w*w*(AKAVKA(3,Ipd)+w*AKAVKA(4,Ipd))
      Effi = Effi*EXP(pw)*AKAVKA(1,Ipd)
      RETURN

C-----------------------------------------------------------------
C     Leuven efficiency calibration
 1003 Effi = AKAVKA(1,Ipd)
      w = LOG(1000.*En)
      DO i = 1 , 6
         Effi = Effi + AKAVKA(i+1,Ipd)*w**i
      ENDDO
      Effi = EXP(Effi)
      RETURN

C-----------------------------------------------------------------
C     Radware efficiency calibration
C     PJN@2008
 1004 w = LOG(En/.1)
      Effi = (AKAVKA(2,Ipd)+(AKAVKA(3,Ipd)+AKAVKA(4,Ipd)*w)*w)
     &       **(-AKAVKA(8,Ipd))
      w = LOG(En)
      Effi = (AKAVKA(5,Ipd)+(AKAVKA(6,Ipd)+AKAVKA(7,Ipd)*w)*w)
     &       **(-AKAVKA(8,Ipd)) + Effi
      Effi = AKAVKA(1,Ipd)*EXP(Effi**(-1/AKAVKA(8,Ipd)))
      RETURN

      END
 
C----------------------------------------------------------------------
C SUBROUTINE ADHOC
C
C Called by: GOSIA
C Calls: READY, SEQ
C
C Purpose: to handle the OP,YIEL option.
C
C Uses global variables:
C      AGELI  - angles of the Ge detectors
C      BRAT   - branching ratio and its error
C      CC     - conversion coefficients for different energies and multipolarities
C      DMIX   - 0.8326 * gamma energy
C      DMIXE  - mixing ratio and its error
C      EAMX   - known matrix elements and their error
C      EG     - energies for conversion coefficients
C      EN     - energy of level
C      ENZ    - depends on absorber
C      IAMX   - index of matrix element for known matrix element
C      IAMY   - level indices of pair of levels for which matrix element is known
C      IBRC   - index branching ratios
C      IDRN   - index of normalising transition for yields
C      IFMO   - include correction to angular distance for finite recoil distance.
C      IMIX   - decay associated with known mixing ratio
C      IPRM   - printing flags (see suboption PRT of OP,CONT)
C      ITMA   - identify detectors according to OP,GDET
C      ITS    - create tape 18 file (OP,CONT switch SEL,)
C      IVAR   - indicates a limit or correlation is set
C      JZB    - unit to read from
C      KSEQ   - index into ELM for pair of levels, and into EN or SPIN
C      LIFCT  - index for lifetimes
C      MEMAX  - number of matrix elements
C      NAMX   - number of known matrix elements
C      NANG   - number of gamma-ray detectors for each experiment
C      NBRA   - number of branching ratios
C      NDL    - number of mixing ratios
C      NDST   - number of data sets
C      NEXPT  - number of experiments
C      NICC   - number of conversion coefficients
C      NLIFT  - number of lifetimes
C      NYLDE  - number of yields
C      ODL    - results of OP,GDET calculation
C      Q      - solid angle attenuation coefficients
C      TAU    - lifetime in picoseconds
C      TIMEL  - lifetimes and their errors
C      UPL    - upper limits for all gamma detectors
C      YNRM   - relative normalization factors for gamma detectors
C
C Formal parameters:
C      Oph    - this indicates the option (GOSI, STAR etc.)
C      Idr    - number of decays
C      Nfd    - number of physical detectors
C      Ntap   - unit of yield file
C      Iyr    - flag set here
C
C Here we parse the input of the OP,YIEL command and store the values.
 
      SUBROUTINE ADHOC(Oph,Idr,Nfd,Ntap,Iyr)
      IMPLICIT NONE
      REAL*8 wamx , wbra , wdl , wlf
      INTEGER*4 iax , Idr , iexp1 , ilft , iosr , ipri , isrt1 , iuf
      INTEGER*4 Iyr , jic , jicc , juf , lb , li , licc , llia , lxt , 
     &          MEM , n1 , n2 , ndas , ndtp , Nfd , nistr , ns1 , ns2 , 
     &          ns3 , ns4 , Ntap , nvare
      CHARACTER*4 Oph
      INTEGER*4 NDST
      COMMON /CCCDS / NDST(50)
      REAL*8 DIX, ODL
      COMMON /DIMX  / DIX(4) , ODL(200)
      REAL*8 DELTA, ENDEC, ENZ
      INTEGER*4 ITMA
      COMMON /TRA   / DELTA(1500,3) , ENDEC(1500) , ITMA(50,200) , 
     &                ENZ(200)
      INTEGER*4 NLIFT
      COMMON /LIFE  / NLIFT
      REAL*8 DMIXE , DMIX
      INTEGER*4 IMIX , NDL
      COMMON /MIXD  / DMIXE(20,2) , DMIX(20) , IMIX(20) , NDL
      REAL*8 EAMX
      INTEGER*4 NAMX, IAMX, IAMY
      COMMON /ME2D  / EAMX(100,2) , NAMX , IAMX(100) , IAMY(100,2)
      REAL*8 TIMEL
      INTEGER*4 LIFCT
      COMMON /LIFE1 / LIFCT(50) , TIMEL(2,50)
      REAL*8 BRAT
      INTEGER*4 IBRC , NBRA
      COMMON /BRNCH / BRAT(50,2) , IBRC(2,50) , NBRA
      REAL*8 YEXP, CORF , DYEX , UPL , YNRM
      INTEGER*4 IY , NYLDE , IDRN , ILE
      COMMON /YEXPT / YEXP(32,1500) , IY(1500,32) , CORF(1500,32) , 
     &                DYEX(32,1500) , NYLDE(50,32) , UPL(32,50) , 
     &                YNRM(32,50) , IDRN , ILE(32)
      REAL*8 YGN , YGP
      INTEGER*4 IFMO
      COMMON /YTEOR / YGN(1500) , YGP(1500) , IFMO
      REAL*8 TAU
      INTEGER*4 KSEQ
      COMMON /LEV   / TAU(100) , KSEQ(1500,4)
      REAL*8 EG, CC, AGELI, Q
      INTEGER*4 NICC, NANG, ISPL
      COMMON /CCC   / EG(50) , CC(50,5) , AGELI(50,200,2) , Q(3,200,8) ,
     &                NICC , NANG(200) , ISPL
      REAL*8 EN , SPIN , ACCUR , DIPOL , ZPOL , ACCA
      INTEGER*4 ISO
      COMMON /COEX  / EN(100) , SPIN(100) , ACCUR , DIPOL , ZPOL , 
     &                ACCA , ISO
      REAL*8 XA , XA1 , EP , TLBDG , VINF
      INTEGER*4 NEXPT , IZ , IZ1
      COMMON /CX    / NEXPT , IZ , XA , IZ1(50) , XA1(50) , EP(50) , 
     &                TLBDG(50) , VINF(50)
      INTEGER*4 MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR
      COMMON /CEXC  / MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR(1500)
      INTEGER*4 IPRM
      COMMON /PRT   / IPRM(20)
      INTEGER*4 ITS
      COMMON /TRB   / ITS
      INTEGER*4 IUNIT3 , JZB
      COMMON /SWITCH/ JZB , IUNIT3
      
C     Read OP,YIEL parameters
      iosr = 0
      READ (JZB,*) IFMO ! IFLAG
      READ (JZB,*) NICC , nistr ! N1, N2
      READ (JZB,*) (EG(jicc),jicc=1,ABS(NICC)) ! E1,E2...
      Iyr = 1
      DO jic = 1 , nistr
        READ (JZB,*) isrt1 ! I1
         IF ( isrt1.GT.6 ) isrt1 = isrt1 - 3
         READ (JZB,*) (CC(jicc,isrt1),jicc=1,ABS(NICC)) ! CC(I1,1)...CC(I1,N1)
      ENDDO
      READ (JZB,*) (NANG(jicc),jicc=1,NEXPT) ! NANG(I)...NANG(NEXPT)

C     Read file for gamma-ray energy dependence of Ge solid-angle attenuation
C     coefficients Q
      REWIND 9
      READ (9,*) Nfd
      DO jicc = 1 , Nfd
         READ (9,*) ODL(jicc) ! DIX(4) - distance from target to front of detector
         READ (9,*) ENZ(jicc) ! Depends on absorber
         DO isrt1 = 1 , 8
            READ (9,*) (Q(licc,jicc,isrt1),licc=1,3)
         ENDDO
      ENDDO

C     Read detector identities, theta and phi
      DO jic = 1 , NEXPT ! For each experiment
         juf = NANG(jic)
         IF ( juf.LT.0 ) THEN ! If NANG < 0 use previous values
            juf = ABS(juf) ! Number of detector angles
            DO jicc = 1 , juf ! For each detector angle
               AGELI(jic,jicc,1) = AGELI(jic-1,jicc,1) ! theta same as previous detector
               AGELI(jic,jicc,2) = AGELI(jic-1,jicc,2) ! phi same as previous detector
               ITMA(jic,jicc) = ITMA(jic-1,jicc)
            ENDDO
            IF ( Oph.NE.'GOSI' ) NANG(jic) = ABS(NANG(jic))
         ELSE
            READ (JZB,*) (ITMA(jic,jicc),jicc=1,juf) ! IP(1)...IP(NANG(I))
            READ (JZB,*) (AGELI(jic,jicc,1),jicc=1,juf) ! Theta Ge det
            READ (JZB,*) (AGELI(jic,jicc,2),jicc=1,juf) ! Phi Ge det
         ENDIF
      ENDDO ! Loop jic on experiments

C     Call SEQ to calculate "chronological" order of levels, so we can
C     account for feeding
      CALL SEQ(Idr)

C     Convert angles into radians
      DO jic = 1 , NEXPT ! For each experiment
         juf = NANG(jic)
         juf = ABS(juf) ! Number of detector angles
         DO jicc = 1 , juf ! For each detector angle
            DO lxt = 1 , 2 ! 1 is theta, 2 is phi
               AGELI(jic,jicc,lxt) = AGELI(jic,jicc,lxt)*.0174532925 ! 0.017452925 = pi / 180
            ENDDO
         ENDDO
      ENDDO ! Loop jic on experiments

C     Set normalising transition
      TAU(1) = 1.E+25 ! Initialise ground-state lifetime to 1E25 picoseconds
      READ (JZB,*) ns1 , ns2 ! NS1, NS2
      DO li = 1 , Idr ! Search through decays for right pair of levels
         IF ( KSEQ(li,3).EQ.ns1 .AND. KSEQ(li,4).EQ.ns2 ) GOTO 100
      ENDDO
 100  IDRN = li ! Index of normalising transition
      IF ( Oph.NE.'GOSI' ) RETURN

C     Read upper limits and relative normalisation factors
      DO li = 1 , NEXPT ! Loop on experiments
         juf = NANG(li)
         IF ( juf.LT.0 ) THEN ! If NANG < 0, use same as previous
            juf = ABS(juf)
            NANG(li) = juf ! Number of detector angles
            NDST(li) = NDST(li-1) ! Number of datasets same as previous
            DO jicc = 1 , juf ! For each detector angle
               UPL(jicc,li) = UPL(jicc,li-1) ! Upper limits same as previous
               YNRM(jicc,li) = YNRM(jicc,li-1) ! Relative normalisation same as previous
            ENDDO
         ELSE
            READ (JZB,*) NDST(li) ! NDST
            ndas = NDST(li)
            READ (JZB,*) (UPL(jicc,li),jicc=1,ndas) ! UPL1...N
            READ (JZB,*) (YNRM(jicc,li),jicc=1,ndas) ! YNRM1...N
         ENDIF
      ENDDO ! Loop li on experiments

C     Read file for experimental yields
      READ (JZB,*) Ntap ! NTAP
      IF ( Ntap.NE.0 ) THEN
         ipri = IPRM(2)
         CALL READY(Idr,Ntap,ipri) ! Read yields from unit Ntap
         ndtp = 0
         DO iexp1 = 1 , NEXPT ! Loop on experiments
            juf = NDST(iexp1) ! Number of datasets
            DO iuf = 1 , juf ! Loop on datasets
               ndtp = ndtp + NYLDE(iexp1,iuf)
            ENDDO
         ENDDO ! Loop iexp1 on experiments

C        Count free variables
         nvare = 0
         DO iexp1 = 1 , MEMAX ! For each matrix element
            IF ( IVAR(iexp1).EQ.1 .OR. IVAR(iexp1).EQ.2 )
     &           nvare = nvare + 1
         ENDDO
         WRITE (22,99001) ndtp , nvare
99001    FORMAT (1X//5X,1I4,1X,'EXPERIMENTAL YIELDS',10X,1I3,1X,
     &           'MATRIX ELEMENTS TO BE VARIED'///)
      ENDIF ! IF ( Ntap.NE.0 )

C     Read branching ratios
      READ (JZB,*) NBRA , wbra ! NBRA, WBRA
      IF ( ITS.EQ.2 ) THEN
         REWIND 18
         WRITE (18,*) MEMAX
      ENDIF
      IF ( NBRA.NE.0 ) THEN
         WRITE (22,99002)
99002    FORMAT (40X,'BRANCHING RATIOS',//5X,'NS1',5X,'NF1',5X,'NS2',5X,
     &           'NF2',5X,'RATIO(1:2)',9X,'ERROR')
         DO lb = 1 , NBRA ! I1,I2,I3,I4,B,DB repeated NBRA times
            READ (JZB,*) ns1 , ns2 , ns3 , ns4 , BRAT(lb,1) , BRAT(lb,2)
            BRAT(lb,2) = BRAT(lb,2)/(SQRT(wbra)+1.E-10) ! Relative error
            WRITE (22,99003) ns1 , ns2 , ns3 , ns4 , BRAT(lb,1) , 
     &                       BRAT(lb,2)
99003       FORMAT (4X,1I3,5X,1I3,5X,1I3,5X,1I3,5X,1F10.5,5X,1F10.5)
            DO li = 1 , Idr ! Search decays for these pairs of levels
               IF ( KSEQ(li,3).EQ.ns3 .AND. KSEQ(li,4).EQ.ns4 ) THEN
                  IBRC(2,lb) = li ! Decay index for first pair
               ELSEIF ( KSEQ(li,3).EQ.ns1 .AND. KSEQ(li,4).EQ.ns2 ) THEN
                  IBRC(1,lb) = li ! Decay index for second pair
               ENDIF
            ENDDO
            IF ( ITS.EQ.2 ) THEN
               n1 = IBRC(1,lb) ! Decay of first pair
               n2 = IBRC(2,lb) ! Decay of second pair
               WRITE (18,*) KSEQ(n1,1) , KSEQ(n2,1)
               WRITE (18,*) KSEQ(n1,1) , KSEQ(n2,2)
               WRITE (18,*) KSEQ(n1,1) , KSEQ(n1,2)
               WRITE (18,*) KSEQ(n2,1) , KSEQ(n1,2)
               WRITE (18,*) KSEQ(n2,1) , KSEQ(n2,2)
               IF ( KSEQ(n1,2).NE.0 .AND. KSEQ(n2,2).NE.0 ) WRITE (18,*)
     &              KSEQ(n1,2) , KSEQ(n2,2)
            ENDIF
         ENDDO
         WRITE (22,99004) wbra
99004    FORMAT (5X,'BRANCHING RATIOS ARE TAKEN WITH WEIGHT',2X,1E14.6)
      ENDIF

C     Read lifetimes
      READ (JZB,*) NLIFT , wlf ! NL, WL
      IF ( NLIFT.NE.0 ) THEN
         WRITE (22,99005)
99005    FORMAT (1X///30X,'LIFETIMES(PSEC)'///5X,'LEVEL',9X,'LIFETIME',
     &           5X,'ERROR'/)
         DO ilft = 1 , NLIFT ! INDEX, T, DT repeated NL times
            READ (JZB,*) LIFCT(ilft) , TIMEL(1,ilft) , TIMEL(2,ilft)
            TIMEL(2,ilft) = TIMEL(2,ilft)/(SQRT(wlf)+1.E-10) ! Relative error
            WRITE (22,99006) LIFCT(ilft) , TIMEL(1,ilft) , TIMEL(2,ilft)
99006       FORMAT (6X,1I3,6X,1F10.2,3X,1F10.2)
         ENDDO
         WRITE (22,99007) wlf
99007    FORMAT (1X/10X,'LIFETIMES ARE TAKEN WITH WEIGHT',2X,1E14.6)
      ENDIF

C     Read known mixing ratios
      READ (JZB,*) NDL , wdl ! NDL, WDL
      IF ( NDL.NE.0 ) THEN
         WRITE (22,99008)
99008    FORMAT (1X//20X,'EXPERIMENTAL E2/M1 MIXING RATIOS'///10X,
     &           'TRANSITION',12X,'DELTA',10X,'ERROR'/)
         DO li = 1 , NDL ! IS, IF, DELTA, ERROR repeated NDL times
            READ (JZB,*) ns1 , ns2 , DMIXE(li,1) , DMIXE(li,2)
            DMIXE(li,2) = DMIXE(li,2)/(SQRT(wdl)+1.E-10)
            WRITE (22,99012) ns1 , ns2 , DMIXE(li,1) , DMIXE(li,2)
            DO lb = 1 , Idr ! Search through decays for right pair of levels
               IF ( KSEQ(lb,3).EQ.ns1 .AND. KSEQ(lb,4).EQ.ns2 ) THEN
                  IMIX(li) = lb ! Decay index
                  DMIX(li) = .8326*(EN(ns1)-EN(ns2)) ! 0.8326 * energy of gamma
                  IF ( ITS.EQ.2 ) WRITE (18,*) KSEQ(lb,1) , KSEQ(lb,2)
               ENDIF
            ENDDO
         ENDDO
         WRITE (22,99009) wdl
99009    FORMAT (/10X,'E2/M1 MIXING RATIOS ARE TAKEN WITH WEIGHT',2X,
     &           1E14.6)
      ENDIF
      IF ( ITS.EQ.2 ) WRITE (18,*) iosr , iosr

C     Read known matrix elements
      READ (JZB,*) NAMX , wamx ! NAMX, WAMX
      IF ( NAMX.EQ.0 ) RETURN
      WRITE (22,99010)
99010 FORMAT (1X//30X,'EXPERIMENTAL MATRIX ELEMENT(S)'///10X,
     &        'TRANSITION',10X,'MAT.EL.',10X,'ERROR'/)

      DO iax = 1 , NAMX ! LAMBDA, INDEX1, INDEX2, ME, DME repeated NAMX times
         READ (JZB,*) llia , ns1 , ns2 , EAMX(iax,1) , EAMX(iax,2)
         IAMY(iax,1) = ns1 ! Level index
         IAMY(iax,2) = ns2 ! Level index
         EAMX(iax,2) = EAMX(iax,2)/(SQRT(wamx)+1.E-10) ! Relative error of ME
         WRITE (22,99012) ns1 , ns2 , EAMX(iax,1) , EAMX(iax,2)
         IAMX(iax) = MEM(ns1,ns2,llia) ! Index to matrix element
      ENDDO
      WRITE (22,99011) wamx
99011 FORMAT (/10X,' MATRIX ELEMENT(S) ARE TAKEN WITH WEIGHT',2X,1E14.6)

99012 FORMAT (9X,1I3,'---',1I3,13X,1F9.4,8X,1F9.4)
      END
 
C----------------------------------------------------------------------
C FUNCTION ELMT
C
C Called by: GOSIA
C Calls:     WTHREJ
C
C Purpose: collective model matrix elements (OP,THEO)
C
C Formal parameters:
C      Xi1    - spin of initial level
C      Xi2    - spin of final level
C      Lam    - multipolarity
C      Nb1    - band number of initial level
C      Nb2    - band number of final level
C      Xk1    - initial level
C      Xk2    - final level
C      Xm1    - intrinsic moment Q1
C      Xm2    - intrinsic moment Q2
C      Xm3    - intrinsic moment Q3
C
C Return value:
C      Collective model matrix element
C
C Note that the parameters to WTHREJ are doubled to allow it to handle half
C integers.
 
      REAL*8 FUNCTION ELMT(Xi1,Xi2,Lam,Nb1,Nb2,Xk1,Xk2,Xm1,Xm2,Xm3)
      IMPLICIT NONE
      REAL*8 addt , fac , fct , pha1 , pha2 , s1 , s2 , WTHREJ , Xi1 , 
     &       Xi2 , Xk1 , Xk2 , xlam , Xm1 , Xm2 , Xm3 , xn
      INTEGER*4 i1 , i2 , ipha , k1 , k2 , l , la , Lam , llam , n , 
     &          Nb1 , Nb2

      la = Lam
      IF ( la.GT.6 ) la = la - 6
      xlam = DBLE(la)
      i1 = INT(2.*Xi1)
      i2 = INT(2.*Xi2)
      llam = 2*la
      k1 = INT(2.*Xk1)
      k2 = INT(2.*Xk2)
      fac = SQRT(2.*Xi1+1.)*SQRT(2.*Xi2+1.)
C-----In-band matrix element

      IF ( Nb1.NE.Nb2 ) THEN
C-----Interband, K-allowed
C-----One K=0
         IF ( ABS(k1-k2).GE.llam ) THEN
C-----Forbidden and K1-K2=lambda, Mikhailov formula
            addt = 0.
            IF ( k1.EQ.1 ) addt = (-1.)**((i1+1)/2)*(i1+1)/2.*Xm3
            xn = ABS(Xk1-Xk2) - xlam
            n = INT(xn+.1)
            IF ( n.EQ.0 ) THEN
               fct = 1.
            ELSEIF ( n.EQ.1 ) THEN
               fct = SQRT((Xi1-Xk1)*(Xi1+Xk1+1.))
            ELSE
               s1 = Xi1 - Xk1
               s2 = Xi1 + Xk1 + 1.
               DO l = 1 , n
                  s1 = s1*(Xi1-Xk1-DBLE(l))
                  s2 = s2*(Xi1+Xk2+1.+DBLE(l))
               ENDDO
               fct = SQRT(s1*s2)
            ENDIF
            pha1 = (-1.)**INT((Xi1-xlam+Xk2)+.1)
            ELMT = fac*pha1*fct*WTHREJ(i1,llam,i2,k2-llam,llam,-k2)
     &             *(Xm1+Xm2*(Xi2*(Xi2+1.)-Xi1*(Xi1+1.))+addt)
         ELSEIF ( k1.NE.0 .AND. k2.NE.0 ) THEN
C-----Both K's non-zero
            pha1 = (-1.)**((i1-llam+k2)/2)
            pha2 = (-1.)**((i1+k1)/2)*pha1
            ELMT = fac*(pha1*WTHREJ(i1,llam,i2,k1,k2-k1,-k2)
     &             *Xm1+pha2*WTHREJ(i1,llam,i2,-k1,k1+k2,-k2)*Xm2)
            RETURN
         ELSE
            ipha = (i1-llam+k2)/2
            IF ( k2.EQ.0 ) ipha = ((i2-llam+k1)/2)
            pha1 = (-1.)**ipha
            ELMT = fac*pha1*WTHREJ(i1,llam,i2,0,k2,-k2)*Xm1
            IF ( k2.EQ.0 ) ELMT = fac*pha1*WTHREJ(i2,llam,i1,0,k1,-k1)
     &                            *Xm1
            IF ( k1.NE.0 .OR. k2.NE.0 ) ELMT = ELMT*SQRT(2.)
            RETURN
         ENDIF
C-----K=0
      ELSEIF ( k1.NE.0 ) THEN
C-----In band, K.ne.0
         pha1 = (-1.)**((i1-llam+k1)/2)
         pha2 = (-1.)**((k1+i1)/2+1)*pha1
         ELMT = fac*(pha1*WTHREJ(i1,llam,i2,k1,0,-k1)
     &          *Xm1+pha2*WTHREJ(i1,llam,i2,-k1,2*k1,-k1)*Xm2)
         RETURN
      ELSE
         ELMT = fac*WTHREJ(i1,llam,i2,0,0,0)*Xm1
         RETURN
      ENDIF
      END
 
C----------------------------------------------------------------------
C SUBROUTINE SELECT
C
C Called by: GOSIA
C
C Purpose: integrate the functionality of the program SELECT into gosia as
C          OP,SELE
C
C PJN April 2008
C
      SUBROUTINE SELECT
      IMPLICIT NONE
      REAL*8 a , al , am , y
      INTEGER*4 i , ie , iexp , indx , ixf , j , l , lm , lu , lum , 
     &        lx , memax
      DIMENSION lm(1500) , y(175,1500) , a(1500,1500)

      ixf = 0
      REWIND 18
      READ (18,*) memax
      DO i = 1 , memax
         DO j = 1 , memax
           a(i,j) = 0.d0
         ENDDO
      ENDDO
      DO i = 1 , 1000
         READ (18,*) l , j
         IF ( l.EQ.0 ) GOTO 100
         IF ( j.NE.0 ) THEN
            a(l,j) = -1.d0
            a(j,l) = -1.d0
         ENDIF
      ENDDO
 100  ie = 1
 200  lum = 0
      DO i = 1 , 175
         DO j = 1 , memax
            y(i,j) = 0.
         ENDDO
      ENDDO
      DO i = 1 , 10000
         READ (18,*) lu , indx , iexp , al
         IF ( iexp.NE.ie ) GOTO 300
         lum = MAX0(lum,lu)
         y(lu,indx) = al
      ENDDO
 300  BACKSPACE 18
      ie = iexp
      IF ( ie.EQ.0 ) ixf = 1
      DO i = 1 , memax
         DO j = i , memax
           IF ( a(i,j).NE.-1.d0 ) THEN
               DO l = 1 , lum
                  a(i,j) = a(i,j) + y(l,i)*y(l,j)
               ENDDO
               a(j,i) = a(i,j)
            ENDIF
         ENDDO
      ENDDO
      IF ( ixf.NE.1 ) GOTO 200
      DO i = 1 , memax
         am = 0.
         DO j = 1 , memax
            am = MAX(a(i,j),am)
         ENDDO
         IF ( am.EQ.0.d0 ) am = 1.
         DO j = 1 , memax
           IF ( a(i,j).NE.-1.d0 ) THEN
               a(i,j) = a(i,j)/am
               IF ( a(i,j).LT..1d0 ) THEN
                 a(i,j) = 0.d0
                  GOTO 350
               ENDIF
            ENDIF
            a(i,j) = 1.d0
 350        CONTINUE
         ENDDO
         WRITE (10,*) (INT(a(i,j)),j=1,memax)
      ENDDO
      DO i = 1 , memax
         WRITE (22,99001) i
99001    FORMAT (10X,'ME=',1I3,3X,'PREDICTED CORRELATION'/)
         lx = 0
         DO j = 1 , memax
           IF ( a(i,j).NE.0.d0 ) THEN
               lx = lx + 1
               lm(lx) = j
            ENDIF
         ENDDO
         WRITE (22,*) (lm(j),j=1,lx)
      ENDDO
      END
C----------------------------------------------------------------------
C SUBROUTINE BRICC
C
C Called by: GOSIA
C Calls:     CCLKUP
C
C Purpose: evaluate internal conversion coefficients using the BrIcc
C          database for each transition energy that gosia needs. The
C          results are stored in the file on unit 29, which is read
C          the first time CONV is called.
C
C Uses global variables:
C      EN     - energy of level
C      IZ     - Z of investigated nucleus
C      LEAD   - pair of levels involved in each matrix element
C      MEMAX  - number of matrix elements
C
      SUBROUTINE BRICC
      IMPLICIT NONE
      REAL*8 temp , egamma
      INTEGER*4 i , j , ngamma
      DIMENSION egamma(1500)
      REAL*8 EN , SPIN , ACCUR , DIPOL , ZPOL , ACCA
      INTEGER*4 ISO
      COMMON /COEX  / EN(100) , SPIN(100) , ACCUR , DIPOL , ZPOL , 
     &                ACCA , ISO
      REAL*8 XA , XA1 , EP , TLBDG , VINF
      INTEGER*4 NEXPT , IZ , IZ1
      COMMON /CX    / NEXPT , IZ , XA , IZ1(50) , XA1(50) , EP(50) , 
     &                TLBDG(50) , VINF(50)
      INTEGER*4 MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR
      COMMON /CEXC  / MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR(1500)
      INTEGER*4 LAMDA , LEAD , LDNUM , LAMMAX , MULTI
      COMMON /CLCOM / LAMDA(8) , LEAD(2,1500) , LDNUM(8,100) , LAMMAX , 
     &                MULTI(8)
      INTEGER*4 IUNIT3 , JZB
      COMMON /SWITCH/ JZB , IUNIT3
      INTEGER*4 n1, n2
      CHARACTER*1024 idx_name, icc_name
      REAL*8 mycc(5), CCLKUP

C     Read the names of the files
      read (JZB,'(A)') idx_name
      read (JZB,'(A)') icc_name

C     Write to output
      write(22,'(/,3A)') 'OP,BRIC interpolation of conversion ',
     &  'coefficients from the BrIcc database, which will be used by',
     &  ' gosia.'
      write(22,'(2A)') 'Please cite T. Kibedi et al. NIM A589 (2008)',
     &  ' 202-229 for the conversion coefficients!'
      write(22,'(2A)') 'Energy [MeV]   E1           E2           E3',
     &  '           M1           M2'

C     Make sure we are at start of file that we want to write
      rewind(29)
      
C     Open the BrIcc database files
      OPEN (UNIT=30,FILE=idx_name,ACCESS='direct',RECL=2048,ERR=999,
     &      STATUS='OLD')
      OPEN (UNIT=31,FILE=icc_name, ACCESS='direct',RECL=44,ERR=999,
     &      FORM='UNFORMATTED',STATUS='OLD')

      ngamma = 0
      DO i = 1 , MEMAX ! For each matrix element
         n1 = LEAD(2,i) ! Upper level
         n2 = LEAD(1,i) ! Lower level
         IF ( n1.EQ.n2 ) GOTO 100 ! Ignore diagonal matrix elements

         temp = ABS(EN(n1) - EN(n2)) ! Energy of transition

C        Now look to see if we have it already
         DO j = 1, ngamma
            IF ( ABS(temp - egamma(j)).LT.1E-6 ) GOTO 100
         ENDDO

C        We get here if we don't have it, so add it to the list
         ngamma = ngamma + 1
         egamma(ngamma) = temp
         mycc(1) = CCLKUP(IZ, temp * 1E3, 1)
         mycc(2) = CCLKUP(IZ, temp * 1E3, 2)
         mycc(3) = CCLKUP(IZ, temp * 1E3, 3)
         mycc(4) = CCLKUP(IZ, temp * 1E3, 6)
         mycc(5) = CCLKUP(IZ, temp * 1E3, 7)
         WRITE(22,'(F7.4,3X,1P,5E13.3)') temp, (mycc(j),j=1,5)
         WRITE(29,'(F7.4,3X,1P,5E13.3)') temp, (mycc(j),j=1,5)
 100     CONTINUE
      ENDDO

C     Close BrIcc database files
      CLOSE (30)
      CLOSE (31)
      RETURN

 999  STOP 'Unable to open BrIcc database files'
      END
C----------------------------------------------------------------------
C FUNCTION NEWCNV
C
C Called by: CONV
C
C Purpose: calculate the conversion coefficient at a particular energy by
C interpolating over the values provided by the user.
C
C Formal parameters:
C      Ega    - gamma energy
C      N      - multipolarity N=1,2,3 = E1,2,3 and N=4,5 = M1,2 (not as elsewhere!)
C
C Return value:
C      conversion coefficient interpolated to energy Ega

      REAL*8 FUNCTION NEWCNV(Ega,N)
      IMPLICIT NONE

      INTEGER*4 isfirst, i, j, N, nenergies
      REAL*8 energies(1500), bricc(1500, 5), Ega
      SAVE energies, bricc, isfirst, nenergies
      DATA isfirst/1/

C     The first time, we need to read the data
      IF ( isfirst.eq.1 ) THEN
        rewind(29)
        isfirst = 0
        DO nenergies = 1, 1500
          READ(29,*,END=100) energies(nenergies),
     &      (bricc(nenergies,j),j=1,5)
        ENDDO
      ENDIF

C     Check multipolarity is valid
 100  IF ( N.LT.1.OR.N.GT.5 ) THEN
         NEWCNV = 0.0
         RETURN
      ENDIF

C     Search for the energy in the list

      DO i = 1, nenergies
        IF (ABS(Ega - energies(i)) .LT. 1E-3) THEN
           NEWCNV = bricc(i,N)
           return
        ENDIF
      ENDDO

C     We get here if the energy isn't in the list, so stop with an error
C     message
      WRITE (*,'(A,F7.3,A)')
     & 'Unable to find conversion coefficients for ',
     &  Ega, ' MeV'
      STOP 'Missing conversion coefficients'

      END
C----------------------------------------------------------------------
C SUBROUTINE SPLNER
C
C Called by: CCLKUP, EFFIX , GOSIA
C Calls:     FUNC, FUNC1, LAGRAN, SPLINE, SPLINT
C
C Purpose: interpolates using a cubic spline
C
C Formal parameters:
C      X      - x-coordinate of input data
C      Yr     - y-coordinate of input data
C      N      - number of data points
C      Xx     - value for which to interpolate
C      Yy     - result of interpolation
C      Iscal  - mode: 1 = linear, 2 = logarithmic, 3 = square root
C
C Note that the effect of FUNC and FUNC1 depends on Iscal:
C Iscal = 1   FUNC(y) = y        FUNC1(y) = y
C Iscal = 2   FUNC(y) = ln(y)    FUNC1(y) = exp(y)
C Iscal = 3   FUNC(y) = sqrt(y)  FUNC1(y) = y^2
C
      SUBROUTINE SPLNER(X,Yr,N,Xx,Yy,Iscal)
      IMPLICIT NONE
      REAL*8 FUNC , FUNC1
      INTEGER*4 N
      REAL*8 X(*) , Yr(*) , w(1500)
      REAL*8 yp1 , ypn , y(1500) , ys
      INTEGER*4 i , Iscal
      REAL*8 Xx , Yy

C     We need at least three points, so if we don't have them, use the
C     Lagrangian method instead
      IF ( N.LE. 3 ) THEN
        CALL LAGRAN(X,Yr,N,1,Xx,Yy,Iscal,1)
        RETURN
      ENDIF 

C     Apply the scaling function
      DO i = 1 , N
        y(i) = FUNC(Yr(i),Iscal)
      ENDDO

C     Get the slope at each end
      yp1 = (y(2)-y(1))/(X(2)-X(1))
      ypn = (y(N)-y(N-1))/(X(N)-X(N-1))

C     Fit a spline
      CALL SPLINE(X,y,N,yp1,ypn,w)

C     Evaluate the spline at the desired point
      CALL SPLINT(X,y,w,N,Xx,ys)

C     Apply the inverse scaling function
      Yy = FUNC1(ys,Iscal)
      RETURN
      END
 
C----------------------------------------------------------------------
C SUBROUTINE SPLINE
C
C Called by: SPLNER
C
C Purpose: fit a spline to data points
C
C Formal parameters:
C      X      - x-coordinate of input data
C      Y      - y-coordinate of input data
C      N      - number of data points
C      Yp1    - slope of first two data points
C      Yyn    - slope of last two data points
C      Y2     - second derivative vector
C
      SUBROUTINE SPLINE(X,Y,N,Yp1,Ypn,Y2)
      IMPLICIT NONE
      INTEGER*4 N
      REAL*8 Yp1 , Ypn , X(*) , Y(*) , Y2(*)
      INTEGER*4 i , k
      REAL*8 p , qn , sig , un , u(1500)
 
      IF ( Yp1.GT..99E30 ) THEN
         Y2(1) = 0.
         u(1) = 0.
      ELSE
         Y2(1) = -.5
         u(1) = (3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-Yp1)
      ENDIF
      DO i = 2 , N - 1
         sig = (X(i)-X(i-1))/(X(i+1)-X(i-1))
         p = sig*Y2(i-1) + 2.
         Y2(i) = (sig-1.)/p
         u(i) = (6.*((Y(i+1)-Y(i))/(X(i+1)-X(i))-(Y(i)-Y(i-1))/(X(i)-X(i
     &          -1)))/(X(i+1)-X(i-1))-sig*u(i-1))/p
      ENDDO
      IF ( Ypn.GT..99E30 ) THEN
         qn = 0.
         un = 0.
      ELSE
         qn = .5
         un = (3./(X(N)-X(N-1)))*(Ypn-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N) = (un-qn*u(N-1))/(qn*Y2(N-1)+1.)
      DO k = N - 1 , 1 , -1
         Y2(k) = Y2(k)*Y2(k+1) + u(k)
      ENDDO
      END
C----------------------------------------------------------------------
C SUBROUTINE SPLINT
C
C Called by: SPLNER
C
C Purpose: evaluate spline
C
C Formal parameters:
C      Xa     - x-coordinate of input data
C      Ya     - y-coordinate of input data
C      Ya2    - second derivative vector
C      N      - number of data points
C      Xx     - value for which to evaluate spline
C      Yy     - result
C
      SUBROUTINE SPLINT(Xa,Ya,Y2a,N,Xx,Yy)
c
c xa,ya - tabulated function
c n - number of tabulated points
c y2a - second derivatives vector
c xx - interpolated point
c yy - intetpoalted value
c
      IMPLICIT NONE
      INTEGER*4 N
      REAL*8 Xx , Yy , Xa(*) , Y2a(*) , Ya(*)
      INTEGER*4 k , khi , klo
      REAL*8 a , b , h
      DATA k/0/
 
      klo = 1
      khi = N
 100  IF ( khi-klo.GT.1 ) THEN
         k = (khi+klo)/2
         IF ( Xa(k).GT.Xx ) THEN
            khi = k
         ELSE
            klo = k
         ENDIF
         GOTO 100
      ENDIF
      h = Xa(khi) - Xa(klo)
C
C  Linear extrapolation
C
      IF ( ABS(h).LT.1E-9 ) THEN
         IF ( Xx.LT.Xa(1) ) h = 1
         IF ( Xx.GT.Xa(N) ) k = N - 1
         a = (Ya(k)-Ya(k+1))/(Xa(k)-Xa(k+1))
         b = (Ya(k)+Ya(k+1)-a*(Xa(k)+Xa(k+1)))*.5
         Yy = a*Xx + b
c         PRINT * , 'splint' , Xx , Yy , 'extrapolation'
         RETURN
      ENDIF
      a = (Xa(khi)-Xx)/h
      b = (Xx-Xa(klo))/h
      Yy = a*Ya(klo) + b*Ya(khi) + ((a**3-a)*Y2a(klo)+(b**3-b)*Y2a(khi))
     &     *(h**2)/6.
      END
C FUNCTION CCLKUP
C
C Called by: BRICC
C Calls:     SPLNER
C
C Purpose: looks up a single conversion coefficient for given Z, energy and
C          multipolarity in the BrIcc database, interpolating for the right
C          energy using a cubic spline on a log-log scale.
C
C Formal parameters:
C      Myz     - Z of nucleus to look up
C      Egamma  - Gamma ray energy for which we want to interpolate
C      Myimult - Multipolarity (1=E1,2=E2,3=E3,6=M1,7=M2)
      REAL*8 FUNCTION CCLKUP(Myz,Egamma,Myimult)
 
      IMPLICIT NONE
 
      INTEGER*4 iz , ia , flag(37) , nrec(37) , rec1(37) , Myz
      INTEGER*4 ishell , irec , ienergy , imult , Myimult
      REAL*8 x(500) , y(500) , result , Egamma
      REAL*4 binding_energy(500) , energy(37,500) , cc(37,500,10)
      CHARACTER*8 file(37) , name(37)
      CHARACTER*4 element
 
C     Initialise
      CCLKUP = 0.0D0

C     We can't calculate above 6000 keV
      IF ( Egamma.GE.6000 ) RETURN
  
C     Read the index record for that Z
      READ (30,REC=Myz,ERR=100) iz , element , ia , 
     &                          (flag(ishell),binding_energy(ishell),
     &                          file(ishell),name(ishell),nrec(ishell),
     &                          rec1(ishell),ishell=1,37)

C     Ignore internal pair conversion below 1100 keV
      IF ( Egamma.LE.1100 ) flag(37) = 0
      
C     Now read the internal conversion coefficient data records for each shell
      DO ishell = 1 , 37 ! Only first 37 subshells

C        Flag that shell doesn't contribute if we are below binding energy
         IF ( Egamma.LE.binding_energy(ishell) ) flag(ishell) = 0

         IF ( flag(ishell).EQ.-1 ) THEN ! If subshell is present in record

C           Read all the records for this subshell
            ienergy = 1
            DO irec = rec1(ishell) , rec1(ishell) + nrec(ishell) - 1
               READ (31,REC=irec,ERR=100) energy(ishell,ienergy) , 
     &               (cc(ishell,ienergy,imult),imult=1,10)
               ienergy = ienergy + 1
            ENDDO ! Loop on records
             
            IF ( cc(ishell,1,Myimult).EQ.0.0D0 ) flag(ishell) = 0
          ENDIF ! If subshell is present in record
      ENDDO ! Loop on subshells
 
C     Interpolate for each subshell
      DO ishell = 1 , 37 ! Only first 37 subshells                       

        IF ( flag(ishell).EQ.-1 ) THEN ! If subshell is present in record
 
C           If the energy is less than the first data point use its ICC, or
C           if it is more than the last data point, otherwise we interpolate
          IF ( Egamma.LE.energy(ishell,1) ) THEN
               result = DBLE(CC(ishell,1,Myimult))
               WRITE (22,'(A,F7.4,3A)') 'Warning Egamma=',Egamma/1.D3,
     &           ' is in regime where solid state effects dominate',
     &           ' conversion coefficients for shell ',name(ishell)
            ELSEIF ( Egamma .GE. energy(ishell,nrec(ishell)) ) THEN
               result = DBLE(cc(ishell,nrec(ishell),Myimult))
               WRITE (22,'(A,F7.4,3A)') 'Warning Egamma=',Egamma/1.D3,
     &           ' exceeds range of conversion coefficients table',
     &           ' for shell ',name(ishell)
            ELSE
C             Set up for interpolation
              DO ienergy = 1 , nrec(ishell)
                 x(ienergy) = LOG(DBLE(energy(ishell,ienergy)))
                 y(ienergy) = DBLE(cc(ishell,ienergy,Myimult)) ! Log for this is done in FUNC
              ENDDO
 
C             Perform spline over data
              CALL SPLNER(x,y,nrec(ishell),LOG(Egamma),result,2)
            ENDIF

C           Add the conversion coefficients of each subshell
            CCLKUP = CCLKUP + result
            
         ENDIF ! If subshell is present in record
      ENDDO ! Loop on subshells
 
      RETURN
 
 100  WRITE (*,*) 'ERROR - No data found for this Z ', Myz
      CCLKUP = -1.0D0
      RETURN
      END
 
C----------------------------------------------------------------------
C SUBROUTINE INVKIN
C
C Called by: GOSIA
C
C Purpose: calculate the angle of the scattered projectile in the lab frame
C          when the user gave the angle of the recoiling target nucleus in
C          the lab frame. There are two solutions to this problem, so Iflag
C          = 1 selects the larger angle one and Iflag = 2 the smaller one.
C          Note that the smaller angle (Iflag = 2) corresponds to very low
C          energies of the recoiling target nucleus, which probably either
C          don't get out of the target or don't get detected. So Iflag = 2
C          is probably not very useful! Also, this routine calculates the
C          correct value of the kinematic flag IKIN.
C
C Formal parameters:
C      E_p     - Beam energy in MeV (readonly)
C      E_x     - energy of excited state to use for kinematic in MeV (readonly)
C      I_Z     - Projectile/target flag. -ve if projectile excitation
C      M_inv   - mass of investigated nuclei in AMU (readonly)
C      M_non   - mass of non-investigated nuclei in AMU (readonly)
C      Theta_t - theta of recoiling target nucleus in lab frame (readonly)
C      Theta_p - theta of scattered projectile in lab frame (writeonly)
C      Iflag   - flag to select one of two possible solutions (readonly)
C      Ikin    - kinematic flag (writeonly)
      
      SUBROUTINE INVKIN(E_p, E_x , I_z , M_inv , M_non , Theta_t ,
     &                  Theta_p , Iflag , Ikin)

      IMPLICIT NONE
      REAL*8 E_p , M_inv , M_non , Theta_t , Theta_p , E_x , M_p , M_t
      REAL*8 ared , epmin , t , x(2), y , thres , tau , taup , TASIN
      INTEGER*4 Iflag , Ikin , I_z

C     Sort out which is the projectile and which is the target
      
      IF ( I_z.LT.0 ) THEN
         M_p = M_inv ! Projectile is investigated
         M_t = M_non ! Target is non investigated
      ELSE
         M_p = M_non ! Projectile is non investigated
         M_t = M_inv ! Target is investigated
      ENDIF

C     Reduced mass
      
      ared = 1 + M_p / M_t
      
C     Excitation energy of inelastically scattered particle when state at
C     energy E_x is excited
      
      epmin = E_p - E_x * ared
      
C     Tau
      
      taup = sqrt(E_p / epmin)
      tau = taup * M_p / M_t
      
C     Calculate the two solutions
      
      y = tan(theta_t/57.2957795)
      t = taup * taup * y * y * y * y -
     &      (1 + y * y) * (taup * taup * y * y - 1)
      t = sqrt(t)
      x(1) = (taup * y * y + t) / (1 + y * y)
      x(1) = atan2(sqrt(1 - x(1) * x(1)), tau + x(1))
      x(2) = (taup * y * y - t) / (1 + y * y)
      x(2) = atan2(sqrt(1 - x(2) * x(2)), tau + x(2))
      
C     Select the solution we want according to the flag. Note that the
C     solution with the lower angle corresponds to target recoils which
C     are probably undetectable.
      
      IF ( Iflag.EQ.1 ) THEN
         theta_p = MAX(x(1),x(2))*57.2957795
      ELSE
         theta_p = MIN(x(1),x(2))*57.2957795
      ENDIF

C     Trap spurious values where tau * sin(theta_p/57.2957795) is greater
C     than unity due to floating point rounding errors
      if (tau * SIN(theta_p/57.2957795) .gt. 1.0d0) 
     &  theta_p = 57.2957795 * TASIN(1.d0/tau)
      
C     Calculate angle of scattered projectile in centre of mass frame, for
C     which the maximum laboratory scattering angle is reached.
      
      t = acos(-1./tau)
      
C     Now calculate the arctangent of the corresponding angle for the
C     recoiling target nuclei in the laboratory frame
      
      thres = sin(t)/(taup-cos(t))
      
C     So now, if y = tan(theta_t_lab) > thres, we are above the maximum and
C     need the larger value of theta_p_cm, so we set Ikin to 1. Otherwise we
C     are below the maximum and need the smaller value so we choose Ikin = 0.
      
      IF ( y.GT.thres ) THEN
         Ikin = 1
      ELSE
         Ikin = 0
      ENDIF
      
      END
