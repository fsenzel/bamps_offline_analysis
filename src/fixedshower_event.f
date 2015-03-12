C---------------------------------------------
C provided by subversion
C---------------------------------------------
C $HeadURL: svn+ssh://senzel@th.physik.uni-frankfurt.de/home/bamps/svn/full/offlineAnalysis/trunk/src/pyshow_event.f $
C $LastChangedDate: $
C $LastChangedRevision: -1 $
C $LastChangedBy: $
C---------------------------------------------
C---------------------------------------------

C   SAMPLE MAIN FOR Q-PYTHIA.1.0.
C   BRANCHING THROUGH PYSHOW OF A JET OF 100 GEV ALONG THE Z AXIS.
C   (based on sample file of Q-PYTHIA by N. Armesto, L. Cunqueiro and 
C    C. A. Salgado, Eur. Phys. J. C63 (2009) 679 [arXiv:0907.1014 [hep-ph]].) 

      SUBROUTINE FIXEDSHOWEREVENT(px,py,pz,flav1,flav2,seed)

C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      INTEGER PYK,PYCHGE,PYCOMP
C...EXTERNAL statement links PYDATA on most machines.
      EXTERNAL PYDATA

C...Commonblocks.
      LOGICAL ACCEP
      DOUBLE PRECISION ejet
      DOUBLE PRECISION PA(6,100)
      INTEGER flav1,flav2
      DOUBLE PRECISION px,py,pz
      INTEGER(kind=8) seed
      
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYDATR/MRPY(6),RRPY(100)
      COMMON/BAMPS/PA

      DO 100 Iqq=1,6,1
        DO 200 Jqq=1,100,1 
          PA(Iqq,Jqq) = 0.D0
  200   CONTINUE
  100 CONTINUE

C... set seed in PYTHIA
      MRPY(1)=seed

c   switches (default refers to PYTHIA-6.4.18 defaults).
c
c   MSTJ(47) and MSTJ(50) are fixed this way because we are running a
c   single parton which does not come from any hard scattering.
c
c   MSTJ(41) is fixed this way because we are only interested in QCD radiation.
c
c   MSTJ(41) must NOT be 11 or 12, as then FSR may go through PYPTFS
c   (kt-ordered cascade) in which medium effects have not been introduced.
c
      MSTJ(41)=1 ! QCD radiation only (default=2)
      MSTJ(42)=2 ! 2: angular ordering (default), 1: no angular ordering
      MSTJ(44)=2 ! options to run alphas (default=2)
      MSTJ(47)=0 ! no correction back to hard scattering element (default=3)
      MSTJ(50)=0 ! no coherence in first branching (default=3)
      PARJ(82)=1.D0 ! GeV, cut off for parton showers (default=1 GeV)            
      
         N=2
         K(1,1)=2
         K(1,2)=flav1 ! gluon = 21, for initial quarks it would be 1,2,3,4,5,6
         K(1,3)=0
         K(1,4)=0
         K(1,5)=0
         P(1,1)=px
         P(1,2)=py
         P(1,3)=pz
         P(1,4)=sqrt(px**2+py**2+pz**2)
         P(1,5)=PYMASS(K(1,2))
         V(1,1)=0.D0
         V(1,2)=0.D0
         V(1,3)=0.D0
         V(1,4)=0.D0
         V(1,5)=0.D0
         K(2,1)=1
         K(2,2)=flav2 ! gluon = 21, for initial quarks it would be 1,2,3,4,5,6
         K(2,3)=0
         K(2,4)=0
         K(2,5)=0
         P(2,1)=-px
         P(2,2)=-py
         P(2,3)=pz
         P(2,4)=sqrt(px**2+py**2+pz**2)
         P(2,5)=PYMASS(K(2,2))
         V(2,1)=0.D0
         V(2,2)=0.D0
         V(2,3)=0.D0
         V(2,4)=0.D0
         V(2,5)=0.D0

c
C         call qpygin0 ! generate initial nucleon-nucleon position.
c

	 ejet = sqrt(px**2+py**2+pz**2)

C          DO 5 I=1, 10, 1
 4       CALL PYSHOW(1,2,2.D0*ejet)
C       CONTINUE

         CALL PYEDIT(5)
       
         ACCEP = .TRUE.
         DO 3 Iqq=1, N, 1
         IF ((K(Iqq,2).EQ.5).OR.(K(Iqq,2).EQ.6).OR.(K(Iqq,2).EQ.-5)
     &   .OR.(K(Iqq,2).EQ.-6).OR.(K(Iqq,2).EQ.4).OR.(K(Iqq,2).EQ.-4)) 
     &    ACCEP = .FALSE.
 3       CONTINUE

         IF (ACCEP) THEN
!         CALL PYLIST(1)
         
         DO 2 Iqq=1, N, 1
           PA(1,Iqq) = K(Iqq,2)  
           PA(2,Iqq) = P(Iqq,1)
           PA(3,Iqq) = P(Iqq,2)
           PA(4,Iqq) = P(Iqq,3)
           PA(5,Iqq) = P(Iqq,4)
           PA(6,Iqq) = P(Iqq,5) 
C           WRITE(*,*) K(Iqq,2), P(Iqq,1), P(Iqq,2), P(Iqq,3),
C     &      P(Iqq,4), P(Iqq,5), K(Iqq,1)
C     &      P(Iqq,4), P(Iqq,5), K(Iqq,1)
C           WRITE(*,*) K(Iqq,2)
 2       CONTINUE
         END IF

      END
