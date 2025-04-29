C     PROGRAM PENTRU GENERAREA PROGRAMULUI UNUI PENDUL DUBLU
$NOTRUNCATE
      INCLUDE 'FGRAPH.FI'
      INCLUDE 'FGRAPH.FD'
      LOGICAL testcolor
      CHARACTER*1 cont

      INTEGER*2 caz,Vprog2(8),vp(10)
      INTEGER*2 ip(40,10),iip(40)
      REAL l1,l2,m1,m2
      REAL M11,M21,M22,J1,J2,M212,A,B,Md121,Md22
      REAL*8 H,TMAX
      COMMON /GEN/ip,iip,fp(40)
      COMMON /PARAM/l1,l2,m1,m2
      COMMON /VECT/Vprog2
      COMMON /COEF1/M11,M21,M22,J1,J2,M212,A,B,Md121,Md22
      DATA Vprog2/2,3,3,3,4,4,1,2/
      DATA l1/1.0/,l2/0.5/,m1/2.0/,m2/1.0/
c
      if(testcolor() )then
C     INITIALIZARI
      write(*,'(''    NGEN=''\)')
      read(*,'(I5)')ngen
      n=40
      nred=n/2
      call gettim(ihr,imin,isec,i100th)
      call seed(i100th)
      call random(vrnd)
C     n POPULATIA
      CALL POPINIT(n)
      cont='G'
      H=1.e-3
      TMAX=20.0
      PE=0.05
C    CALCULUL COEFICIENTILOR COEF1 CONSTANTI PT. UN PENDUL ALES
      g=9.81
      M11=m1*l1*l1/4.
      M21=m2*l1*l1
      M22=m2*l2*l2/4.
      J1=m1*l1*l1/12.
      J2=m2*l2*l2/12.
      M212=m2*l1*l2/2.
      A=J1+M11+J2+M21+M22
      B=J2+M22
      Md121=0.5*g*l1*(m1+2*m2)
      Md22=0.5*g*m2*l2
      igen=0
      do while(igen.lt.ngen)
1     do i=1,n
c     parametri de simulare pt. un cromozom
c     ip(k,1)=1  caz=1
c     ip(k,1)=2  caz=0
      caz=ip(i,1)
      TiM=0.1*(ip(i,2)+5)
      do j=1,8
      Vprog2(j)=ip(i,j+2)
      end do
C     SIMULAREA
      CALL SIMUP(caz,H,TMAX,TiM,PE,KP,ThnRM)
      fp(i)=ThnRM
      end do
C     ORDONARE
      CALL ORDP(n)
c     AFISAREA ELITEI
      igen=igen+1
      write(*,'(''  GENERATIA '',I5)')igen
      do i=1,10
      write(*,'(10I4,I8)')(ip(i,j),j=1,10),iip(i)
      write(*,'(''  fp='',F7.3)')fp(i)
      end do
cc      write (*,*) 'DACA DORITI CONTINUARE APASATI TASTA G  !'
cc      read (*,'(A)') cont
cc      if(cont.eq.'G'.or.cont.eq.'g') then
cc      read(*,*)
C      REORDONAREA INDICILOR
      do i=1,n
      iip(i)=i
      end do
c     recombinarea & mutatia
      CALL RECOMB(nred)
c      CALL MUTATIE(nred)
      CALL MUTATIE(n)
cc      else
cc      STOP
cc      end if
cc      go to 1
      end do
      write(*,*)  '        PENTRU ANIMATIE APASATI ENTER !'
      read(*,*)
C     SIMULAREA CU ANIMATIE
      TMAX=50.0
      do i=1,10
      vp(i)=ip(1,i)
      end do
      CALL DSIMUP(caz,H,TMAX,TiM,PE,KP)
      read(*,*)
      else
      write(*,*) 'NU MERGE GRAFICA'
      end if
      status=setvideomode($DEFAULTMODE)
      STOP
      END
      SUBROUTINE DSIMUP(caz,H,TMAX,TiM,PE,KP)
      INTEGER RUNGE
      INTEGER*2 caz,stel1,stel2,stel2n,Vprog2(8)

      REAL l1,l2,m1,m2
      REAL*8 Y(2),F(2),T,H,TMAX,Th1d,Th2d,Th1op,Th2op
      REAL M11,M21,M22,J1,J2,M212,A,B,Md121,Md22,
     * C,MC212,MS212
      INTEGER*2 clr
      COMMON /PARAM/l1,l2,m1,m2
      COMMON /VECT/Vprog2

      COMMON /COEF1/M11,M21,M22,J1,J2,M212,A,B,Md121,Md22
      COMMON /COEF2/C,MC212,MS212
C    CALCULUL COEFICIENTILOR COEF1 CONSTANTI PT.
      g=9.81
      M11=m1*l1*l1/4.
      M21=m2*l1*l1
      M22=m2*l2*l2/4.
      J1=m1*l1*l1/12.
      J2=m2*l2*l2/12.
      M212=m2*l1*l2/2.
      A=J1+M11+J2+M21+M22
      B=J2+M22
      Md121=0.5*g*l1*(m1+2*m2)
      Md22=0.5*g*m2*l2
C     POZITIILE INITIALE
      if (caz.eq.1)then
      call Pozinit1(Th1,Th1P,Th2,stel1,stel2)
      Oo=0.0
      else
      call Pozinit2(Th1,Th1P,Th2,stel1,stel2)
      Oo=-0.1
      end if
C     DATE PT.INITIALIZARE CICLU DE INTEGRARE
      KP=0
      nR=0

      T=0.D0
      flsp=1
c     pregg
      CALL pregg()
      Th1d=Th1
      Th2d=Th2
      icont=0
      clr=15
      CALL dpendul(Th1d,Th2d,clr)
      Th1op=Th1d
      Th2op=Th2d
C     INCEPE RASNIREA
      do while(T.lt.TMAX )
      if(flsp.eq.1)then
C     eliminarea succesiunilor nepermise
        if(stel2.ne.2 .and.stel2.ne.4)then
        stel2n=Vprog2(stel1)
        if(stel2.eq.1)then
        if(stel2n.eq.4 .or. stel2n.eq.3)then
        stel2n=stel2
        end if
        else if(stel2n.eq.2 .or. stel2n.eq.1)then
        stel2n=stel2
        end if
      stel2=stel2n
      end if
      end if
      Th1o=Th1
C     CALCULUL COEFICIENTILOR COEF2
      C=A+2.*M212*cos(Th2)
      MC212=M212*cos(Th2)
      MS212=M212*sin(Th2)
C     POZITIE/PARAM.MISCARE EL2
      call par_el2s(stel2,PE,TiM,Th2,Th2P,Th2PP)
      KP=KP+1
      TE=KP*PE
      Y(1)=Th1
      Y(2)=Th1P
      do while(T.le.TE .or. K.eq.1)
      K=RUNGE(2,Y,F,T,H)
      if(K.eq.1)then
      call DINAMIC(F,Y,Th2,Th2P,Th2PP)
      end if
      end do
      Th1=Y(1)
      Th1P=Y(2)
c      write(*,111)T,Th1,Th2
c      write(*,'(2I6)') stel1,stel2
c     dpendul
      Th1d=Th1
      Th2d=Th2
      CALL timp()
      clr=0
      CALL dpendul(Th1op,Th2op,clr)
      clr=15
      CALL dpendul(Th1d,Th2d,clr)
      Th1op=Th1d
      Th2op=Th2d
C     SCHIMBARE STARE PG FUNCTIE DE Th1 POZITIA EL1
      flsp=0
          if(stel1.eq.7)then
            if(Th1o.le.Th1)then

            stel1=8
            flsp=1
            end if
          else if(stel1.eq.8)then
	    if(Th1.ge.Oo)then

            stel1=1
            flsp=1
            end if
          else if(stel1.eq.1)then
	    if(Th1.gt.Oo)then

            stel1=2
            flsp=1
            end if
	 else if(stel1.eq.2)then
	    if(Th1o.ge.Th1 .and. Th1o.gt.Oo)then

            stel1=3
            nR=nR+1
	    flsp=1
	    Th1R=Th1o
            end if
        else if(stel1.eq.3)then
            if(Th1o.gt.Th1)then

            stel1=4
            flsp=1
            end if
         else if(stel1.eq.4)then
	    if(Th1.LE.Oo)then
	    stel1=5
            flsp=1
            end if
          else if(stel1.eq.5)then
            if(Th1o.ge.Th1)then

            stel1=6
            flsp=1
            end if
          else if(stel1.eq.6)then
            if(Th1o.le.Th1)then

            stel1=7
	    flsp=1
	    Th1L=Th1o
	    Oo=0.5*(Th1R+Th1L)
            end if
         end if
      end do
      return
      end

      subroutine popinit(n)
      INTEGER*2 ip(40,10),iip(40)
      COMMON /GEN/ip,iip,fp(40)
      np2=n/2
      do i=1,np2
      ip(i,1)=1
      ip(i,2)=irnd_m(10)
      do j=3,10
      ip(i,j)=irnd_m(4)
      end do
      end do
      do i=np2+1,n
      ip(i,1)=2
      ip(i,2)=irnd_m(10)
      do j=3,10
      ip(i,j)=irnd_m(4)
      end do
      end do
      do i=1,n
      iip(i)=i
      end do
      return
      end

      subroutine ordp(n)
C     ORDONARE DUPA PERFORMANTA
      integer*2 ip(40,10),iip(40),k,itemp,itmp(40,10)
      COMMON /GEN/ip,iip,fp(40)
c num. cromozomi n
c     ordonare fp si iip in ordine descrescatoare fp
      lim=n-1
      do i=1,lim
      in=i+1
      do j=in,n
      if(fp(i).lt.fp(j))then
      temp=fp(j)
      fp(j)=fp(i)
      fp(i)=temp
      itemp=iip(j)
      iip(j)=iip(i)
      iip(i)=itemp
      end if
      end do
      end do
c     ordonare cromozomi
      do j=1,n
      k=iip(j)
      do i=1,10
      itmp(j,i)=ip(k,i)
      end do
      end do
c     copiere in forma de lucru
      do i=1,n
      do j=1,10
      ip(i,j)=itmp(i,j)
      end do
      end do
      return
      end
      SUBROUTINE RECOMB(m)
      INTEGER*2 ip(40,10),iip(40)
      COMMON /GEN/ip,iip,fp(40)
      j=1
      do while(j.lt.m)
      kc1=irnd_m(m)
      kc2=irnd_m(m)
      do while(kc1.eq.kc2)
      kc2=irnd_m(m)
      end do
      ks=irnd_m(10)
      mpj=m+j
      mpj1=mpj+1
      call copycrom(ip,kc1,mpj)
      call copycrom(ip,kc2,mpj1)
      call random(prnd)
      if(ks.eq.1)then
      call cross(ip,mpj,mpj1,1)
      else if(ks.eq.2)then
      call cross(ip,mpj,mpj1,2)
      else if(prnd.lt.0.5)then
c      call random(prnd)
c      if(prnd.lt.0.5)then
      call crosst(ip,mpj,mpj1,ks)
      else
      call cross11(ip,mpj,mpj1)
      end if
c      end if
      j=j+2
      end do
      return
      end
      subroutine copycrom(ip,kc,mp)
c     COPIAZA CROMOZOM
      INTEGER*2 ip(40,10)
      do j=1,10
      ip(mp,j)=ip(kc,j)
      end do
      return
      end

      subroutine cross(ip,ic1,ic2,ig)
c     incrucisare de o gena
      INTEGER*2 ip(40,10),itemp
      itemp=ip(ic1,ig)
      ip(ic1,ig)=ip(ic2,ig)
      ip(ic2,ig)=itemp
      return
      end
      subroutine crosst(ip,ic1,ic2,kt)
C     INCRUCISARE LA PUNCTUL DE TAIERE kt
      INTEGER*2 ip(40,10),itemp(8)
      do j=3,kt
      itemp(j)=ip(ic1,j)
      end do
      do j=3,kt
      ip(ic1,j)=ip(ic2,j)
      ip(ic2,j)=itemp(j)
      end do
      return
      end
      subroutine cross11(ip,ic1,ic2)
C     SCHIMBA GENELE IMPARE INTRE ELE
      INTEGER*2 ip(40,10),itemp
      do j=3,10,2
      itemp=ip(ic1,j)
      ip(ic1,j)=ip(ic2,j)
      ip(ic2,j)=itemp
      end do
      return
      end

      subroutine mutatie(n)
C     REALIZEAZA MUTATIA LA DOI CROMOZOMI
      INTEGER*2 ip(40,10),iip(40)
      COMMON /GEN/ip,iip,fp(40)
      do j=1,5
      id1=irnd_m(n-10)
      id2=irnd_m(8)
      id3=irnd_m(4)
      ip(id1+10,id2+2)=id3
      end do
      return
      end

      SUBROUTINE SIMUP(caz,H,TMAX,TiM,PE,KP,ThnRM)
      INTEGER RUNGE
      INTEGER*2 caz,stel1,stel2,stel2n,Vprog2(8)

      REAL l1,l2,m1,m2
      REAL*8 Y(2),F(2),T,H,TMAX
      REAL M11,M21,M22,J1,J2,M212,A,B,Md121,Md22,
     * C,MC212,MS212
      COMMON /PARAM/l1,l2,m1,m2
      COMMON /VECT/Vprog2
c      COMMON /REZ/F1(1500),F2(1500),F3(1)
      COMMON /COEF1/M11,M21,M22,J1,J2,M212,A,B,Md121,Md22
      COMMON /COEF2/C,MC212,MS212
C     POZITIILE INITIALE
      if (caz.eq.1)then
      call Pozinit1(Th1,Th1P,Th2,stel1,stel2)
      Oo=0
      else
      call Pozinit2(Th1,Th1P,Th2,stel1,stel2)
      Oo=-0.1
      end if
cc	stel2=Vprog2(stel1)
C     DATE PT.INITIALIZARE CICLU DE INTEGRARE

      KP=0
      nR=0
      T=0.D0
      flsp=1

C     INCEPE RASNIREA
c      write(*,111)T,Th1,Th2
c  111 FORMAT (' ',3(3X,F12.7))
c      F1(1)=Th1
c      F2(1)=Th2
      do while(T.lt.TMAX .and. nR.ne.8)
      ThnRM=0.0
      if(flsp.eq.1)then
C     eliminarea succesiunilor nepermise
        if(stel2.ne.2 .and.stel2.ne.4)then
        stel2n=Vprog2(stel1)
        if(stel2.eq.1)then
        if(stel2n.eq.4 .or. stel2n.eq.3)then
        stel2n=stel2
	end if
	 else if(stel2n.eq.2 .or. stel2n.eq.1)then
        stel2n=stel2
        end if
      stel2=stel2n
      end if
c      stel2=stel2n
      end if
      Th1o=Th1
C     CALCULUL COEFICIENTILOR COEF2
      C=A+2.*M212*cos(Th2)
      MC212=M212*cos(Th2)
      MS212=M212*sin(Th2)
C     POZITIE/PARAM.MISCARE EL2
      call par_el2s(stel2,PE,TiM,Th2,Th2P,Th2PP)
      KP=KP+1
      TE=KP*PE
      Y(1)=Th1
      Y(2)=Th1P
      do while(T.le.TE .or. K.eq.1)
      K=RUNGE(2,Y,F,T,H)
      if(K.eq.1)then
      call DINAMIC(F,Y,Th2,Th2P,Th2PP)
      end if
      end do
      Th1=Y(1)
      Th1P=Y(2)
c      write(*,111)T,Th1,Th2
c      write(*,'(2I6)') stel1,stel2
c      F1(KP+1)=Th1
c      F2(KP+1)=Th2
C     SCHIMBARE STARE PG FUNCTIE DE Th1 POZITIA EL1
      flsp=0
          if(stel1.eq.7)then
            if(Th1o.le.Th1)then

            stel1=8
            flsp=1
            end if
          else if(stel1.eq.8)then
	    if(Th1.ge.Oo)then

            stel1=1
            flsp=1
            end if
          else if(stel1.eq.1)then
	    if(Th1.gt.Oo)then

            stel1=2
            flsp=1
            end if
         else if(stel1.eq.2)then
	    if(Th1o.ge.Th1 .and. Th1o.gt.Oo)then

            stel1=3
            nR=nR+1
            ThnRM=Th1o
	    flsp=1
	    Th1R=Th1o
            end if
        else if(stel1.eq.3)then
            if(Th1o.gt.Th1)then

            stel1=4
            flsp=1
            end if
         else if(stel1.eq.4)then
	    if(Th1.LE.Oo)then

            stel1=5
            flsp=1
            end if
          else if(stel1.eq.5)then
            if(Th1o.ge.Th1)then

            stel1=6
            flsp=1
            end if
          else if(stel1.eq.6)then
            if(Th1o.le.Th1)then

            stel1=7
	    flsp=1
	    Th1L=Th1o
	    Oo=0.5*(Th1R+Th1L)
            end if
         end if
      if(nR.ne.8)ThnRM=0.0
      end do
      return
      end

      FUNCTION RUNGE (N,Y,F,X,H)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 Y,F,X,H
      INTEGER RUNGE
      DIMENSION PHI(5),SAVEY(5),Y(N),F(N)
      DATA M/0/
      M=M+1
      GO TO (1,2,3,4,5),M
    1 RUNGE =1
      RETURN
    2 DO 22 J=1,N
      SAVEY(J)=Y(J)
      PHI(J)=F(J)
   22 Y(J)=SAVEY(J)+0.5*H*F(J)
      X=X+0.5*H
      RUNGE=1
      RETURN
    3 DO 33 J=1,N
      PHI(J)=PHI(J)+2.0*F(J)
   33 Y(J)=SAVEY(J)+0.5*H*F(J)
      RUNGE=1
      RETURN
    4 DO 44 J=1,N
      PHI(J)=PHI(J)+2.0*F(J)
   44 Y(J)=SAVEY(J)+H*F(J)
      X=X+0.5*H
      RUNGE=1
      RETURN
    5 DO 55 J=1,N
   55 Y(J)=SAVEY(J)+(PHI(J)+F(J))*H/6.0
      M=0
      RUNGE=0
      RETURN
      END

      subroutine par_el2s(stel2,PE,TiM,Th2,Th2P,Th2PP)
C     excitatie sinusoidala
      INTEGER*2 stel2
      PI=3.1415926
      PIP2=PI/2.0
C     STATIONARE
      if(stel2.eq.1 .or. stel2.eq.3)then
      if(stel2.eq.1)then
      Ti=0.0
      Th2=PIP2
      Th2P=0.0
      Th2PP=0.0
      else
      Ti=0.0
      Th2=0.0
      Th2P=0.0
      Th2PP=0.0
      end if
      else
C     MISCARE
      if(stel2.eq.4)then
      Ti=Ti+PE
      FI=Ti*PI/TiM
      if(Ti.ge.TiM)then
      stel2=1
      Ti=0.0
      FI=PI
      end if
      FIP=PI/TiM
      else
C    ( stel2=2 )
      Ti=Ti+PE
      FI=PI-Ti*PI/TiM
      if(Ti.ge.TiM)then
      stel2=3
      Ti=0.0
      FI=0.0
      end if
      FIP=-PI/TiM
      end if
      Th2=PI*(1.0-cos(FI))/4.0
      Th2P=PI*FIP*sin(FI)/4.0
      Th2PP=PI*FIP*FIP*cos(FI)/4.0
      end if
      return
      end
      subroutine DINAMIC(F,Y,Th2,Th2P,Th2PP)
      REAL M11,M21,M22,J1,J2,M212,A,B,Md121,Md22,
     * C,MC212,MS212
      COMMON /COEF1/M11,M21,M22,J1,J2,M212,A,B,Md121,Md22
      COMMON /COEF2/C,MC212,MS212
      REAL*8 Y(1),F(1),Y1,YP
      Y1=Y(1)
      YP=Y(2)
      F(1)=YP
      F(2)=(MS212*(2.0*YP*Th2P+Th2P*Th2P)-(MC212+B)*Th2PP-
     *   Md121*dsin(Y1)-Md22*dsin(Y1+Th2))/C
      return
      end

      subroutine Pozinit1(Th1,Th1P,Th2,stel1,stel2)
      INTEGER*2 stel1,stel2
C     simplificata Th2=PIP2
      stel1=7
      stel2=1
      Th2=0.5*3.1415926
      Th1=-0.126
      Th1P=0.0
      return
      end

      subroutine Pozinit2(Th1,Th1P,Th2,stel1,stel2)
      INTEGER*2 stel1,stel2
C     simplificata Th2=0
      stel1=5
      stel2=3
      Th2=0.0
      Th1=-0.0
      Th1P=0.0
      return
      end

      function irnd_m(nn)
      integer*2 irnd
      irnd=1
      ptres=1.0/nn
      tres=ptres
      call random(rnd)
      do while(rnd.gt.tres .and. irnd.lt.nn)
      irnd=irnd+1
      tres=tres+ptres
      end do
      irnd_m=irnd
      return
      end
c      subroutine random_m(vrnd)
c      call gettim(ihr,imin,isec,i100th)
c      call seed(i100th)
c      it=i100th/10
c      if(it.lt.1)it=1
c      do i=1,it
c      call random(vrnd)
c      end do
c      return
c      end
      SUBROUTINE pregg()
      INCLUDE 'FGRAPH.FD'
      INTEGER*2 dummy,halfx,halfy,xdom,ydom,clr
cc      LOGICAL testcolor
cc      DOUBLE PRECISION wxy
cc      RECORD/rccoord/curpos
      RECORD/videoconfig/myscreen
      RECORD/wxycoord/wxy
      COMMON myscreen

      CALL clearscreen($GCLEARSCREEN)
      xdom=myscreen.numxpixels-1
      ydom=myscreen.numypixels-1
      halfx=xdom/2
      halfy=ydom/2
c  VIEWPORT
      CALL setviewport(0,0,xdom,ydom)
c   WINDOW
      dummy=setwindow(.TRUE.,-2.,-2.5,2.,0.5)
c    BORDER
      clr=15
      dummy=setcolor(clr)
      dummy=rectangle_w($GBORDER,-2.,0.5,2.,-2.5)
c   TRIUNGHI
      CALL moveto_w(0.,0.,wxy)
      dummy=lineto_w(0.1,0.1)
      dummy=lineto_w(0.1,0.1)
      dummy=lineto_w(-0.1,0.1)
      dummy=lineto_w(0.,0.)
      dummy=ellipse_w($GBORDER,-0.015,0.015,0.015,-0.015)
      return
      END

      SUBROUTINE dpendul(Th1,Th2,clr)
      INCLUDE  'FGRAPH.FD'
      INTEGER*2 dummy,clr
      DOUBLE PRECISION wx1,wy1,wx2,wy2,Th1,Th2
      REAL l1,l2
      RECORD/videoconfig/myscreen
      RECORD/wxycoord/wxy
      COMMON myscreen
      l1=1.
      l2=0.5
      P=3.141592/2.
      dummy=setcolor(clr)
      CALL moveto_w(0.,0.,wxy)
      wx1=l1*dcos(Th1)
      wy1=l1*dsin(Th1)
      wx2=l2*dcos(Th1+Th2)+wx1
      wy2=l2*dsin(Th1+Th2)+wy1
      dummy=lineto_w(wy1,-wx1)
      dummy=lineto_w(wy2,-wx2)
      return
      end
      SUBROUTINE timp()
      mult=1000000
      pi=3.141592
      DO i=1,mult
      r=pi*pi/pi*pi/100.
      end do
      return
      end
      LOGICAL FUNCTION testcolor()
      INCLUDE 'FGRAPH.FD'
      INTEGER*2 status
      RECORD/videoconfig/myscreen
      COMMON myscreen
      status =setvideomode($VRES16COLOR)
      CALL getvideoconfig(myscreen)
      testcolor=.TRUE.
      IF (status.EQ.0)testcolor=.FALSE.

      END
