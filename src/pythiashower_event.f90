!--------------------------------------------
! provided by subversion
!--------------------------------------------
! $HeadURL: svn+ssh://senzel@th.physik.uni-frankfurt.de/home/bamps/svn/full/offlineAnalysis/trunk/src/fullpythia_event.f90 $
! $LastChangedDate: $
! $LastChangedRevision: -1 $
! $LastChangedBy: $
!--------------------------------------------
!---------------------------------------------

      SUBROUTINE PYTHIASHOWEREVENT(ptMin,seed)
!...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP,COUNTER
      character*1 tab
!...EXTERNAL statement links PYDATA on most machines.
      EXTERNAL PYDATA

!...Commonblocks.
      LOGICAL ACCEP
      DOUBLE PRECISION ptMin
      DOUBLE PRECISION PA(7,100)
      INTEGER(kind=8) seed
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYDATR/MRPY(6)
      COMMON/BAMPS/PA

      DO 100 Iqq=1,7,1
        DO 200 Jqq=1,100,1 
          PA(Iqq,Jqq) = 0.D0
  200   CONTINUE
  100 CONTINUE

!.. set seed in PYTHIA
      MRPY(1)=seed

      CKIN(3)=ptMin

!.. TUNE D6T
!       MSTP(5)=109

      MSTP(61)=0       ! master switch for initial-state QCD and QED radiation.
      MSTP(91)=0       ! primordial kt switched off

!tell PYTHIA to use LHAPDF
      MSTP(52)=2
!choose PDF 10042 -> CTEQ6l (LO fit/NLO alphas)
      MSTP(51)=10042

!fragmentation switches
      MSTP(71)=1       !  master switch for final-state QCD and QED radiation.
      MSTP(111)=0      !  master switch for fragmentation and decay, as obtained with a PYEXEC call.
      MSTP(81)=0       !  Multiple Interactions
!       MSTJ(21)=0       !  form of particle decays.

      CALL PYINIT('CMS','p','p',2.76d+3) ! LHC collisions
      CALL PYEVNT ! Generate one PYTHIA event.

!       CALL PYEDIT(5)
     
      ACCEP = .TRUE.
      DO 3 Iqq=1, N, 1
      IF ((K(Iqq,2).EQ.5).OR.(K(Iqq,2).EQ.6).OR.(K(Iqq,2).EQ.-5).OR.(K(Iqq,2).EQ.-6)&
     .OR.(K(Iqq,2).EQ.4).OR.(K(Iqq,2).EQ.-4)) ACCEP = .FALSE.
3     CONTINUE

      IF (ACCEP) THEN
!       CALL PYLIST(1)
      DO 2 Iqq=1, N, 1
        PA(1,Iqq) = K(Iqq,2)
        PA(2,Iqq) = K(Iqq,3)
        PA(3,Iqq) = P(Iqq,1)
        PA(4,Iqq) = P(Iqq,2)
        PA(5,Iqq) = P(Iqq,3)
        PA(6,Iqq) = P(Iqq,4)
        PA(7,Iqq) = P(Iqq,5) 
!         WRITE(*,*) K(Iqq,2), K(Iqq,3), P(Iqq,1), P(Iqq,2), P(Iqq,3), P(Iqq,4), P(Iqq,5)
2     CONTINUE
      CALL PYLIST(2)
      
      END IF
      
      END
