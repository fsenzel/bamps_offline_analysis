!--------------------------------------------
! provided by subversion
!--------------------------------------------
! $HeadURL: svn+ssh://senzel@th.physik.uni-frankfurt.de/home/bamps/svn/full/offlineAnalysis/trunk/src/fullpythia_event.f90 $
! $LastChangedDate: $
! $LastChangedRevision: -1 $
! $LastChangedBy: $
!--------------------------------------------
!---------------------------------------------

      SUBROUTINE PYTHIASHOWEREVENT(SHOWER_TYPE,FLAVOR,SQRTS,P0,SEED)
!...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP,COUNTER
      character*1 tab
!...EXTERNAL statement links PYDATA on most machines.
      EXTERNAL PYDATA

!...Commonblocks.
      LOGICAL CHARM_FOUND, BOTTOM_FOUND, TOP_FOUND
      DOUBLE PRECISION PA(7,50)
      INTEGER SHOWER_TYPE
      INTEGER FLAVOR
      DOUBLE PRECISION SQRTS
      DOUBLE PRECISION P0
      INTEGER(kind=8) SEED
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYDATR/MRPY(6)
      COMMON/BAMPS/PA

!.. initialize particle vector
      DO 100 Iqq=1,7,1
        DO 200 Jqq=1,100,1 
          PA(Iqq,Jqq) = 0.D0
  200   CONTINUE
  100 CONTINUE

!.. set seed in PYTHIA
      MRPY(1)=SEED

!.. minimum momentum transfer
      CKIN(3)=P0

!.. set heavy quark mass to standard value of BAMPS
      PMAS(4,1)=1.3D0
      PMAS(5,1)=4.6D0

!.. set processes depending on requested shower type
      SELECT CASE ( SHOWER_TYPE )
        !.. inclusive pQCD processes
        CASE( 2 )
          MSEL=1

        !.. only pQCD processes with photons
        CASE ( 3 )
          MSEL=0
          MSUB(14)=1 ! f_i + fbar_i -> g + gamma
          MSUB(29)=1 ! f_i + g -> f_i + gamma
          MSUB(115)=1 ! g + g -> g + gamma
!         MSUB(114)=1 ! g + g -> gamma + gamma
!         MSUB(18)=1 ! f_i + fbar_i -> gamma + gamma

        !.. only heavy quark pQCD proceses ( without gluon splitting )
        CASE ( 4, 5 )
	      MSEL = shower_type
!         MSUB(81)=1 ! f_i + fbar_i -> Q_k + Qbar_k
!         MSUB(82)=1 ! g + g -> Q_k + Qbar_k
!         MSUB(83)=1 ! q_i + f_j -> Q_k + f_l
!         MSUB(84)=1 ! g + gamma -> Q_k + Qbar_k
!         MSUB(85)=1 ! gamma + gamma -> F_k + Fbar_k

        CASE DEFAULT
          WRITE( *, * ) 'Unknown shower type in PYTHIA fortran code.'
          CALL ABORT
      END SELECT

!.. TUNE D6T
!     MSTP(5)=109

      MSTP(61)=0       ! master switch for initial-state QCD and QED radiation.
      MSTP(91)=0       ! primordial kt switched off

!.. tell PYTHIA to use LHAPDF
      MSTP(52)=2
!.. choose PDF 10042 -> CTEQ6l (LO fit/NLO alphas)
      MSTP(51)=10042

!.. fragmentation switches
      MSTP(71)=1       !  master switch for final-state QCD and QED radiation.
      MSTP(111)=0      !  master switch for fragmentation and decay, as obtained with a PYEXEC call.
      MSTP(81)=0       !  Multiple Interactions
!       MSTJ(21)=0       !  form of particle decays.

!.. intialize p+p collision
      CALL PYINIT('CMS','p','p',sqrtS) ! 2.76 p+p LHC collisions

!.. generate PYTHIA event
400   CALL PYEVNT ! Generate one PYTHIA event.

!.. debugging
!      CALL PYEDIT(5)

!.. check for heavy quarks
      CHARM_FOUND = .FALSE.
      BOTTOM_FOUND = .FALSE.
      TOP_FOUND = .FALSE.
      DO 3 Iqq=1, N, 1
        IF( ABS(K(Iqq,2)).EQ.4 ) CHARM_FOUND = .TRUE. ! charm quarks was found
        IF( ABS(K(Iqq,2)).EQ.5 ) BOTTOM_FOUND = .TRUE. ! bottom quarks was found
        IF( ABS(K(Iqq,2)).EQ.6 ) TOP_FOUND = .TRUE. ! top quarks was found
3     CONTINUE

!.. reject inappropriate events
      IF( TOP_FOUND ) GOTO 400 ! event has top quark -> reject
      SELECT CASE( SHOWER_TYPE )
        CASE( 2 )
          IF( FLAVOR .EQ. 101 .AND. ( CHARM_FOUND .OR. BOTTOM_FOUND ) ) THEN
            GOTO 400 ! heavy quark was found in light parton shower
          ELSE IF( ABS(FLAVOR) .EQ. 4 .AND. .NOT. CHARM_FOUND ) THEN
            GOTO 400 ! charm is needed but not found
          ELSE IF( ABS(FLAVOR) .EQ. 5 .AND. .NOT. BOTTOM_FOUND ) THEN
            GOTO 400 ! bottom is needed but not found
          END IF
        CASE( 3 )
          IF( CHARM_FOUND .OR. BOTTOM_FOUND .OR. TOP_FOUND ) GOTO 400 ! heavy quark was found in photon shower
        CASE( 4 )
          IF( .NOT. CHARM_FOUND ) GOTO 400 ! no charm in charm shower was found
        CASE( 5 )
          IF( .NOT. BOTTOM_FOUND ) GOTO 400 ! no bottom in bottom shower was found
      END SELECT

!.. save particles
!      CALL PYLIST(1)
      DO 2 Iqq=1, N, 1
        PA(1,Iqq) = K(Iqq,2)
        PA(2,Iqq) = K(Iqq,3)
        PA(3,Iqq) = P(Iqq,1)
        PA(4,Iqq) = P(Iqq,2)
        PA(5,Iqq) = P(Iqq,3)
        PA(6,Iqq) = P(Iqq,4)
        PA(7,Iqq) = P(Iqq,5)
2     CONTINUE

!.. debugging
!      CALL PYLIST(2)
      
      END
