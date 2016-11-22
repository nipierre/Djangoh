C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C         For test run

      subroutine tstrrr
      IMPLICIT REAL*8 (A-H,M,O-Z)
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HSFIJK/ F1(2,2),F2(2,2),F3(2,2)
      write(*,*) ' '
      write(*,*) ' '
      write(*,*) ' '
      write(*,*) ' '
      write(*,*) ' TSTRRR: '
      write(*,*) ' xmin,xmax = ',xmin,xmax
      write(*,*) ' ymin,ymax = ',ymin,ymax
      do 10 ix=1,11
      x=xmin+(xmax-xmin)/1d1*dfloat(ix-1)
      if (ix.eq.1) x=x+1d-3
      GSP=SP-MEI2-MPRO2
      YMAXX=X*(1D0-4D0*MEI2*MPRO2/GSP/GSP)
     *        /(X*(1D0+X*MPRO2/GSP)+MEI2/GSP)
      ymaxl=dmin1(ymax,ymaxx)
      yminx=q2min/x/(sp-mei2-mpro2)
      yminl=dmax1(ymin,yminx)
      write(*,*) ' '
      write(*,*) ' x = ',x
      write(6,*) '      y,             F2,             F1,',
     *'             FL (0) '
      do 20 iy=1,31
      y=ymaxl-(ymaxl-yminl)/(2d0**dfloat(iy))
      CALL HSDELO(X,Y)
      q2=x*y*sp
      call HSSTRF(X,Q2)
      fl=f2(1,1)-2d0*x*f1(1,1)
      write(6,200) y,f2(1,1),f1(1,1),fl
      Q2L=Q2
c      TTR=6D0
c      IF (Q2L.LT.TTR) THEN
c        Q2L=TTR
c        XP21=SNGL(X)*SNGL(SP-MPRO2-MEI2)
c        IF (Q2L.GE.XP21) Q2L=XP21
c      ENDIF
      F2EM=F2(1,1)
      call HSLUFL(X,Q2L,F2EM,FLC)
      write(6,202) q2,q2l,flc
      call HSSTR1(X,Q2)
      fl=f2(1,1)-2d0*x*f1(1,1)
      write(6,201) f2(1,1),f1(1,1),fl
 20   continue
 10   continue
 200  format(f14.10,3e16.6)
 201  format(14x,3e16.6)
 202  format(46x,2f8.4,e16.6)
      return
      end












