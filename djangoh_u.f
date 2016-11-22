Ch 08/08/05
Ch *  updated for LEPTO65
C
      SUBROUTINE HSUSER(ICALL,X,Y,Q2)
C...user analysis routine:
C   ICALL=1 initialization: before generating events
C        =2 analysis:after  event have been generated
C        =3 final:after all generations
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
C...declarations for heracles
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSNUME/ SIGTOT,SIGTRR,SIGG(20),SIGGRR(20),NEVENT,NEVE(20)
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSIKP/  S,T,U,SS,TS,US,DKP,DKPS,DKQ,DKQS
      COMMON /HSLABP/ EH,PH,EQH,PQH,ESH,PSH,COSEH,SINEH
      PARAMETER (NMXHEP=2000)
      COMMON /HEPEVT/ NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &                JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),
     &                PHEP(5,NMXHEP),VHKK(4,NMXHEP)
      COMMON /HSCHNN/ ICHNN
      COMMON /HSONLY/ IHSONL
      CHARACTER*80 OUTFILENAM
      COMMON /HSOUTF/ OUTFILENAM,ICH
      REAL RTIME
C...declarations for lepto65
      COMMON /DJPASS/ NTOT,NPASS,NFAILL
      COMMON /DJFAIL/ NFAILI(10)
      COMMON/LEPTOU/CUT(14),LST(40),PARL(30),XSCH,YSCH,W2SCH,Q2SCH,USCH
      REAL          CUT            ,PARL    ,XSCH,YSCH,W2SCH,Q2SCH,USCH
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      REAL                      P        ,V
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      REAL                    PARU               ,PARJ
      INTEGER MSTU,MSTJ
C...Number of events from sophia
      COMMON /SPPASS/ NSOPH,NSPOUT,NFAILP,NSPACC
C...Hbook declarations
      INTEGER NWPAWC,CHANNELS,LIMITS
      REAL HMEMORY
      double precision HBINW
      PARAMETER (NWPAWC=100000)
      PARAMETER (CHANNELS=50,LIMITS=CHANNELS+1)
      COMMON/PAWC/HMEMORY(NWPAWC)
      INTEGER NID(1:200)
      DATA NID/200*0/
      double precision BINW(1:200)
      DATA BINW/200*0./
C
      DIMENSION IFLCNT(-6:6)
      DIMENSION NFAILC(0:10)
      DIMENSION NMIS(0:12)
      DIMENSION PSUM(4)
      LOGICAL LFIRST
      DATA LFIRST /.TRUE./
      DATA NEVMOD/10000/
      SAVE IFLCNT
C
      IF(LFIRST) THEN
        LFIRST=.FALSE.
        IEIN=1
        DO 10 I=-6,6
   10   IFLCNT(I)=0
        DO 11 I=0,10
   11   NFAILC(I)=0
        DO 12 I=0,12
   12   NMIS(I)=0
        CALL TIMEX(RTIME)
        TIMINI=RTIME
        NEVMOD=10000
        N10CNT=0
        GSP=SP-MPRO2-MEI2
        OPEN(31,FILE=OUTFILENAM(1:ICH)//'_evt.dat',STATUS='NEW')
        OPEN(6,FILE='luevents.dat',STATUS='NEW')
	CLOSE(6)

      ENDIF

C-----------------------------------------------------------------------
C...Initialization of user action
      GOTO (100,200,300) ICALL

 100  CONTINUE
C-----------------------------------------------------------------------
C...initialize event
      DO 13 I=-6,6
   13  IFLCNT(I)=0
C...initialize hbook
      CALL HLIMIT(NWPAWC)
      Q2MNL=DLOG10(Q2MIN)
      Q2MXL=DLOG10(Q2MAX)
      CALL HBOOK1(11,'<D[S]>/<D LOG10 Q2>',
     & CHANNELS,SNGL(Q2MNL),SNGL(Q2MXL),0.)
      BINW(11)=HBINW(11)
      XMIN1=XMIN
      XMINL=DLOG10(XMIN1)
      XMAX1=XMAX
      XMAXL=DLOG10(XMAX1)
      CALL HBOOK1(12,'<D[S]>/<D LOG10 X>',
     & CHANNELS,SNGL(XMINL),SNGL(XMAXL),0.)
      BINW(12)=HBINW(12)
      YMIN1=YMIN
      YMAX1=YMAX
      CALL HBOOK1(13,'<D[S]>/<D Y>',
     & CHANNELS,SNGL(YMIN1),SNGL(YMAX1),0.)
      BINW(13)=HBINW(13)
      WMIN1=WMIN
      WMAX1=SQRT(SP)
      CALL HBOOK1(14,'<D[S]>/<D W>',
     & CHANNELS,SNGL(WMIN1),SNGL(WMAX1),0.)
      BINW(14)=HBINW(14)
      WMIN1=0.5D0
      WMAX1=1.2D0
      CALL HBOOK1(15,'<D[S]>/<D W>',
     & CHANNELS,SNGL(WMIN1),SNGL(WMAX1),0.)
      BINW(15)=HBINW(15)
      YMIN1=DLOG10(YMIN)
      YMAX1=DLOG10(YMAX)
      CALL HBOOK1(16,'<D[S]>/<D LOG10 Y>',
     & CHANNELS,SNGL(YMIN1),SNGL(YMAX1),0.)
      BINW(16)=HBINW(16)
      WMIN1=DLOG10(MPRO)
      WMAX1=DLOG10(SP)/2D0
      CALL HBOOK1(17,'<D[S]>/<D LOG10 W>',
     & CHANNELS,SNGL(WMIN1),SNGL(WMAX1),0.)
      BINW(17)=HBINW(17)
      PMISMN=-1D0
      PMISMX= 1D0
      CALL HBOOK1(21,'<D[S]>/<D Pt>',
     & CHANNELS,SNGL(PMISMN),SNGL(PMISMX),0.)
      BINW(21)=HBINW(21)
      PMISMN=-2D0
      PMISMX= 2D0
      CALL HBOOK1(22,'<D[S]>/<D Pl>',
     & CHANNELS,SNGL(PMISMN),SNGL(PMISMX),0.)
      BINW(22)=HBINW(22)
      CALL HBOOK1(23,'<D[S]>/<D E>',
     & CHANNELS,SNGL(PMISMN),SNGL(PMISMX),0.)
      BINW(23)=HBINW(23)
      CALL HBOOK2(24,'<D[S]>/<D E><D W>',
     & CHANNELS,SNGL(PMISMN),SNGL(PMISMX),
     & CHANNELS,SNGL(WMIN1),SNGL(WMAX1),0.)
      BINW(24)=HBINW(24)
      XNLIN=0D0
      XMLIN=0.21D0
      CALL HBOOK1(31,'<D[S]>/<D X>',
     & CHANNELS,SNGL(XNLIN),SNGL(XMLIN),0.)
      BINW(31)=HBINW(31)
      YNLIN=0.99D0
      YMLIN=1D0
      CALL HBOOK1(32,'<D[S]>/<D Y>',
     & CHANNELS,SNGL(YNLIN),SNGL(YMLIN),0.)
      BINW(32)=HBINW(32)
      ENLIN=26D0
      EMLIN=70D0
      CALL HBOOK1(33,'<D[S]>/<D E>',
     & CHANNELS,SNGL(ENLIN),SNGL(EMLIN),0.)
      BINW(33)=HBINW(33)
      CALL HBOOK1(41,'<D[S]>/<D Q2H>',
     & CHANNELS,1.0,5.0,0.)
      BINW(41)=HBINW(41)
      CALL HBOOK1(42,'<D[S]>/<D XH>',
     & CHANNELS,0.0,1.0,0.)
      BINW(42)=HBINW(42)
C...for more accurate statistics (large event numbers):
      CALL HBARX(0)
      RETURN

 200  CONTINUE
C-----------------------------------------------------------------------
C...user's action with generated event
      IF (ICHNN.EQ.6.OR.ICHNN.EQ.7.OR.ICHNN.EQ.8.OR.ICHNN.EQ.12) THEN
        DKPRO=EPRO*PHEP(4,3)+PPRO*PHEP(3,3)
        Q2HAD=-TS
        XHAD=Q2HAD/(Y*GSP-2D0*DKPRO)
        YHAD=Q2HAD/XHAD/GSP
      ELSE
        Q2HAD=Q2SCH
        XHAD=XSCH
        YHAD=YSCH
      ENDIF
      W2HAD=YHAD*(1D0-XHAD)*GSP+MPRO2
      WHAD=SQRT(W2HAD)
C...Q2-, x-, y- distributions
      DSIGHB=SIGTOT/DFLOAT(NEVENT)
      DSIGH1=DSIGHB/BINW(11)
      Q2L=DLOG10(Q2)
c      CALL HFF1(11,NID(11),SNGL(Q2L),SNGL(DSIGH1))
      CALL HFILL(11,SNGL(Q2L),0.,SNGL(DSIGH1))
      DSIGH1=DSIGHB/BINW(12)
      XLOG=DLOG10(X)
c      CALL HFF1(12,NID(12),SNGL(XLOG),SNGL(DSIGH1))
      CALL HFILL(12,SNGL(XLOG),0.,SNGL(DSIGH1))
      DSIGH1=DSIGHB/BINW(13)
c      CALL HFF1(13,NID(13),SNGL(Y),SNGL(DSIGH1))
      CALL HFILL(13,SNGL(Y),0.,SNGL(DSIGH1))
      DSIGH1=DSIGHB/BINW(14)
c      CALL HFF1(14,NID(14),SNGL(WHAD),SNGL(DSIGH1))
      CALL HFILL(14,SNGL(WHAD),0.,SNGL(DSIGH1))
      DSIGH1=DSIGHB/BINW(15)
c      CALL HFF1(15,NID(15),SNGL(WHAD),SNGL(DSIGH1))
      CALL HFILL(15,SNGL(WHAD),0.,SNGL(DSIGH1))
      DSIGH1=DSIGHB/BINW(16)
      YLOG=DLOG10(Y)
c      CALL HFF1(16,NID(16),SNGL(YLOG),SNGL(DSIGH1))
      CALL HFILL(16,SNGL(YLOG),0.,SNGL(DSIGH1))
      WLOG=DLOG10(WHAD)
      DSIGH1=DSIGHB/BINW(17)
c      CALL HFF1(17,NID(17),SNGL(WLOG),SNGL(DSIGH1))
      CALL HFILL(17,SNGL(WLOG),0.,SNGL(DSIGH1))
      DSIGH1=DSIGHB/BINW(31)
c      CALL HFF1(31,NID(31),SNGL(X),SNGL(DSIGH1))
      CALL HFILL(31,SNGL(X),0.,SNGL(DSIGH1))
      DSIGH1=DSIGHB/BINW(32)
c      CALL HFF1(32,NID(32),SNGL(Y),SNGL(DSIGH1))
      CALL HFILL(32,SNGL(Y),0.,SNGL(DSIGH1))
      DSIGH1=DSIGHB/BINW(33)
c      CALL HFF1(33,NID(33),SNGL(ESH),SNGL(DSIGH1))
      CALL HFILL(33,SNGL(ESH),0.,SNGL(DSIGH1))

      DSIGH1=DSIGHB/BINW(41)
      Q2H=DLOG10(Q2HAD)
c      CALL HFF1(41,NID(41),SNGL(Q2H),SNGL(DSIGH1))
      CALL HFILL(41,SNGL(Q2H),0.,SNGL(DSIGH1))
      DSIGH1=DSIGHB/BINW(42)
c      CALL HFF1(42,NID(42),SNGL(XHAD),SNGL(DSIGH1))
      CALL HFILL(42,SNGL(XHAD),0.,SNGL(DSIGH1))

C...count flavors
      IF (NSPACC.EQ.0) THEN
        DO 20 I=-6,6
   20   IF (LST(25).EQ.I) IFLCNT(I)=IFLCNT(I)+1
      ENDIF

C...count errors in DJGEVT
      DO 21 I=0,10
   21 IF (NFAILI(I).NE.0) NFAILC(I)=NFAILC(I)+1

      IF (IHSONL.EQ.0) THEN
C...Check energy-momentum conservation
      IF (LST(21).EQ.0) THEN
      CALL LUEDIT(1)
C...Print every 'nevmod'th event
      IF (MOD(NEVHEP,NEVMOD).EQ.0) THEN
        CALL LULIST(1)
        WRITE(6,*) ' LST(21) = ',LST(21),'    MSTU(24) = ',MSTU(24),
     .                                   '      NSPACC = ',NSPACC
      ENDIF
C
      DO 22 J=1,4
   22 PSUM(J)=0D0
      DO 23 I=1,N
      DO 23 J=1,4
   23 IF (K(I,1).GT.0.AND.K(I,1).LE.10) PSUM(J)=PSUM(J)+P(I,J)
      PSUMT=SQRT(PSUM(1)**2+PSUM(2)**2)
      IF (PSUMT.GT.0.1D0) NMIS(1)=NMIS(1)+1
      IF (PSUMT.GT.1.0D0) NMIS(2)=NMIS(2)+1
      IF (PSUMT.GT.2.0D0) NMIS(3)=NMIS(3)+1
      IF (PSUMT.GT.10.0D0) NMIS(4)=NMIS(4)+1
      IF (ABS(PSUM(3)-PELE+PPRO).GT.0.1D0) NMIS(5)=NMIS(5)+1
      IF (ABS(PSUM(3)-PELE+PPRO).GT.1.0D0) NMIS(6)=NMIS(6)+1
      IF (ABS(PSUM(3)-PELE+PPRO).GT.2.0D0) NMIS(7)=NMIS(7)+1
      IF (ABS(PSUM(3)-PELE+PPRO).GT.10.0D0) NMIS(8)=NMIS(8)+1
      IF (ABS(PSUM(4)-EELE-EPRO).GT.0.1D0) NMIS(9)=NMIS(9)+1
      IF (ABS(PSUM(4)-EELE-EPRO).GT.1.0D0) NMIS(10)=NMIS(10)+1
      IF (ABS(PSUM(4)-EELE-EPRO).GT.2.0D0) NMIS(11)=NMIS(11)+1
      IF (ABS(PSUM(4)-EELE-EPRO).GT.10.0D0) NMIS(12)=NMIS(12)+1
      PSUML=PSUM(3)-PELE+PPRO
      ESUM=PSUM(4)-EELE-EPRO
      DSIGH1=DSIGHB/BINW(21)
      CALL HFILL(21,SNGL(PSUMT),0.,SNGL(DSIGH1))
      DSIGH1=DSIGHB/BINW(22)
      CALL HFILL(22,SNGL(PSUML),0.,SNGL(DSIGH1))
      DSIGH1=DSIGHB/BINW(23)
      CALL HFILL(23,SNGL(ESUM),0.,SNGL(DSIGH1))
      DSIGH1=DSIGHB/BINW(24)
      CALL HFILL(24,SNGL(ESUM),SNGL(WHAD),SNGL(DSIGH1))
      IF ((PSUMT.GT.10D0.OR.ABS(PSUM(3)-PELE+PPRO).GT.10D0
     &   .OR.ABS(PSUM(4)-EELE-EPRO).GT.10D0).AND.N10CNT.LT.3) THEN
        N10CNT=N10CNT+1
        WRITE(6,*) ' '
        WRITE(6,*) ' '
        WRITE(6,*) ' '
        WRITE(6,*) ' DJUSER: E-P-mismatch > 10GeV: '
        WRITE(6,*) ' NTOT = ',NTOT,'    NEVHEP = ',NEVHEP
        WRITE(6,*) ' ICHNN = ',ICHNN
        WRITE(6,*) ' PSUM = ',PSUM(1),PSUM(2)
        WRITE(6,*) '        ',PSUM(3),PSUM(4)
        WRITE(6,*) ' LST(21) = ',LST(21),'    MSTU(24) = ',MSTU(24)
        CALL LULIST(1)
        WRITE(6,*) ' '
        WRITE(6,*) ' X = ',X
        WRITE(6,*) ' Y = ',Y
        WRITE(6,*) ' Q2 = ',Q2
        WRITE(6,*) ' XHAD = ',XHAD
        WRITE(6,*) ' YHAD = ',YHAD
        WRITE(6,*) ' Q2HAD = ',Q2HAD
        WRITE(6,*) ' WHAD = ',WHAD
        WMOM=(P(2,4)+P(3,4))**2-(P(2,1)+P(3,1))**2
     &      -(P(2,2)+P(3,2))**2-(P(2,3)+P(3,3))**2
        WMOM=SQRT(WMOM)
        WRITE(6,*) ' WHAD from momenta on LUJETS = ',WMOM
        WMOM=P(2,4)*P(3,4)-P(2,1)*P(3,1)-P(2,2)*P(3,2)-P(2,3)*P(3,3)
        WMOM=(PHEP(4,6)+PHEP(4,5)-PHEP(4,1)-PHEP(4,3))**2
     &      -(PHEP(1,6)+PHEP(1,5)-PHEP(1,1)-PHEP(1,3))**2
     &      -(PHEP(2,6)+PHEP(2,5)-PHEP(2,1)-PHEP(2,3))**2
     &      -(PHEP(3,6)+PHEP(3,5)-PHEP(3,1)-PHEP(3,3))**2
        WMOM=SQRT(WMOM)
        WRITE(6,*) ' WHAD from momenta on HEPEVT = ',WMOM
        IF (ICHNN.GT.2) THEN
          OMEGA=(PH*DKQ+PQH*DKP)/(PH*EQH+PQH*EH)
          DCTHGA=(DKP/OMEGA-MEF2/2D0/EH)/PH
          CTHGA=1D0-DCTHGA
          DKPRO=OMEGA*(EPRO+CTHGA*PPRO)
          WHSA=Y*(1D0-X)*(SP-MPRO2)+MPRO2-2D0*(DKP-DKPS+DKPRO)
          WHSA=SQRT(WHSA)
          WRITE(6,*) ' WHAD as in HSACPT = ',WHSA
        ENDIF
      ENDIF
      ENDIF
      ENDIF

C...average time per event
      CONTINUE
      IF(MOD(NEVHEP,NEVMOD).EQ.0) THEN
        WRITE(LUNOUT,2200) NEVHEP
        WRITE(LUNOUT,2201) ICHNN
        WRITE(LUNOUT,2202) 'leptonic: ',X,Y,Q2
        WRITE(LUNOUT,2202) 'hadronic: ',XHAD,YHAD,Q2HAD
        WRITE(LUNOUT,2202) 'scaled:   ',XSCH,YSCH,Q2SCH
        WRITE(LUNOUT,2203) PHEP(1,1),PHEP(2,1),PHEP(3,1),PHEP(4,1)
     &                    ,PHEP(5,1)
        WRITE(LUNOUT,2204) PHEP(1,2),PHEP(2,2),PHEP(3,2),PHEP(4,2)
     &                    ,PHEP(5,2)
        WRITE(LUNOUT,2205) PHEP(1,3),PHEP(2,3),PHEP(3,3),PHEP(4,3)
     &                    ,PHEP(5,3)
 2200 FORMAT(//,' ***** Event No ',I12)
 2201 FORMAT(/,' CHANNEL = ',I3,9x,'X',12x,'Y',12x,'Q2')
 2202 FORMAT(1x,A10,9x,2(E9.3,4x),E9.3)
 2203 FORMAT('-----------',
     &      9x,'Px',10x,'Py',10x,'Pz',10x,'E ',10x,'M'/,
     &       ' HS e-scat: ',5F12.4)
 2204 FORMAT(' HS q-scat: ',5F12.4)
 2205 FORMAT(' HS ga-rad: ',5F12.4)
 2206 FORMAT(/,' LST(21) = ',I3)
        WRITE(LUNOUT,2206) LST(21)
        CALL LULIST(1)
        CALL TIMEX(RTIME)
        TIMFIN=RTIME
        RTIME=TIMFIN-TIMINI
        TIMINI=TIMFIN
        TIMEVT=RTIME/NEVMOD
        WRITE(LUNOUT,2001) NEVHEP,TIMEVT
      ENDIF
C...write events to file
      WRITE(31,2301) ICHNN
C...Added by Nicolas for test purpose
      WRITE(31,2303) X,Y,Q2,XHAD,YHAD,Q2HAD
C... -
      WRITE(31,2302) IDHEP(1)
      WRITE(31,2303) PHEP(1,1),PHEP(2,1),PHEP(3,1),PHEP(4,1),PHEP(5,1)
      WRITE(31,2302) IDHEP(2)
      WRITE(31,2303) PHEP(1,2),PHEP(2,2),PHEP(3,2),PHEP(4,2),PHEP(5,2)
      WRITE(31,2302) IDHEP(3)
      WRITE(31,2303) PHEP(1,3),PHEP(2,3),PHEP(3,3),PHEP(4,3),PHEP(5,3)
C...Write LULIST result to file for each event
      OPEN(6,FILE='luevents.dat',POSITION='append',STATUS='OLD')
      CALL LULIST(1)
C... -
      WRITE(6,2301) ICHNN
      WRITE(6,2303) X,Y,Q2
      CLOSE(6)
C... -
 2301 FORMAT(I4)
 2302 FORMAT(I4)
 2303 FORMAT(5E16.8)
      RETURN

 300  CONTINUE
C-----------------------------------------------------------------------
C...final call, overall output
C...program performance
C...save all histograms
        CALL HRPUT(0,OUTFILENAM(1:ICH)//'_his.paw','N')

      IF (IHSONL.EQ.0) THEN
        WRITE(LUNOUT,3001) NTOT,NPASS,NSOPH,NFAILL,NFAILP
        IF (NTOT.NE.NPASS) WRITE(LUNOUT,3003) PARL(24)
C...flavor distribution; new flavor coding (ckc)
        WRITE(LUNOUT,3004) IFLCNT(1),IFLCNT(2),IFLCNT(3)
     &                    ,IFLCNT(4),IFLCNT(5),IFLCNT(6)
     &                    ,IFLCNT(-1),IFLCNT(-2),IFLCNT(-3)
     &                    ,IFLCNT(-4),IFLCNT(-5),IFLCNT(-6)
C...failures in hadronization
        WRITE(LUNOUT,3005) NFAILC(1),NFAILC(2),NFAILC(3)
     &                    ,NFAILC(4),NFAILC(5),NFAILC(6)
     &                    ,NFAILC(7),NFAILC(8),NFAILC(9)
     &                    ,NFAILC(10)
C...energy-momentum mismatch
        WRITE(LUNOUT,3006) NMIS(1),NMIS(2),NMIS(3),NMIS(4)
     &                    ,NMIS(5),NMIS(6),NMIS(7),NMIS(8),NMIS(9)
     &                    ,NMIS(10),NMIS(11),NMIS(12)
      ENDIF

      RETURN

2001  FORMAT(/,' AVERAGE TIME PER EVENT AFTER ',I10,' EVENTS: ',
     &            1PE11.2,' SEC')
C
3001  FORMAT(/,' ******************************************************'
     F,'************************',/
     F      ,' Program performance ',/
     F      ,1X,I12,' Events were accepted by HERACLES,',/
     F      ,1X,I12,' Events passed fragmentation in LEPTO',/
     F      ,1X,I12,' Events passed fragmentation in SOPHIA,',/
     F      ,1X,I12,' Events failed fragmentation in LEPTO',/
     F      ,1X,I12,' Events failed fragmentation in SOPHIA')
3003  FORMAT(' Cross section was corrected: ',/
     F      ,' Not all events could be hadronized. ',/
     F      ,' Total cross section is now    SIGTOT = ',E12.5,' pb',/)
3004  FORMAT(/,' Distribution of flavors (DJ):',/
     F        ,10X,'d',9X,'u',9X,'s',9X,'c',9X,'b',9X,'t',/,1X,6I10,/
     F        ,7X,'dbar',6X,'ubar',6X,'sbar',6X,'cbar',6X,'bbar',6X
     F        ,'tbar',/,1X,6I10,/)
3005  FORMAT(/,' Errors in DJGEVT: Hadronization failed:',/
     F        ,' Nfail(1) = ',I8,'     Nfail(2)  = ',I8,/
     F        ,' Nfail(3) = ',I8,'     Nfail(4)  = ',I8,/
     F        ,' Nfail(5) = ',I8,'     Nfail(6)  = ',I8,/
     F        ,' Nfail(7) = ',I8,'     Nfail(8)  = ',I8,/
     F        ,' Nfail(9) = ',I8,'     Nfail(10) = ',I8,/)
3006  FORMAT(/,' Energy-momentum mismatch: ',/
     F        ,' PmisT >  0.1: NMIS(1)  = ',I8,/
     F        ,' PmisT >  1.0: NMIS(2)  = ',I8,/
     F        ,' PmisT >  2.0: NMIS(3)  = ',I8,/
     F        ,' PmisT > 10.0: NMIS(4)  = ',I8,/
     F        ,' PmisL >  0.1: NMIS(5)  = ',I8,/
     F        ,' PmisL >  1.0: NMIS(6)  = ',I8,/
     F        ,' PmisL >  2.0: NMIS(7)  = ',I8,/
     F        ,' PmisL > 10.0: NMIS(8)  = ',I8,/
     F        ,'     E >  0.1: NMIS(9)  = ',I8,/
     F        ,'     E >  1.0: NMIS(10) = ',I8,/
     F        ,'     E >  2.0: NMIS(11) = ',I8,/
     F        ,'     E > 10.0: NMIS(12) = ',I8,/)

      END
C
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C...dummy routines
C...when ARIADNE is not linked: uncomment the following two routines
C
      SUBROUTINE ARINIT(MODE)
      CHARACTER MODE*(*)
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON/LEPTOU/CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
       WRITE(LUNOUT,101)
  101  FORMAT(/'*************************************************',
     &   /'WARNING: Ariadne is inactive because:',
     &   /'       1/ ARIADNE program is not linked or',
     &   /'       2/ dummy routines in DJUSER are active',
     &   /'         (should be commented out when ARIADNE is linked',
     &   /'  CASCADE flag set to 12, i.e. ME+PS option',
     &   /'*************************************************'/)
        LST(8)=12
      RETURN
      END
cck---------------------------------------------------------------------
      SUBROUTINE AREXEC
cck..dummy
      RETURN
      END
cck---------------------------------------------------------------------
      SUBROUTINE ARPRDA
cck..dummy
      RETURN
      END

************************************************************************

C FUNCTION RETURNS THE BIN WIDTH OF A HISTOGRAM
      FUNCTION HBINW(ID)
      IMPLICIT NONE
      double precision HBINW
      INTEGER ID
      CHARACTER*80 CHTITL
      INTEGER NX,NY,NWT,LOC
      REAL XMI,XMA,YMI,YMA
      CALL HGIVE(ID,CHTITL,NX,XMI,XMA,NY,YMI,YMA,NWT,LOC)
      IF (NY.EQ.0) THEN
         HBINW=(DBLE(XMA)-DBLE(XMI))/DBLE(NX)
      ELSE
         HBINW=(DBLE(XMA)-DBLE(XMI))*(DBLE(YMA)-DBLE(YMI))
     .        /(DBLE(NX)*DBLE(NY))
      ENDIF
      RETURN
      END
