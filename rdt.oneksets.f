c   integration for rotating rdt, general dudy and rot
c       with townsend dki/dt = -kndUn/dxi
c   initialize isotropic for wave number radius = ak1 (=1)
      parameter(idim=20000)
      real*8    pa(6),pas(6),con(3),pi,dth,theta,dpa(6),cont(3)
      real*8    dispt(6),pkk(6,6)
      real*8    phij(9,idim),ak(idim,3),dphij(9,idim)
      real*8    phijknkm(3,3,3,3),em(3,3,3,3),cijsn(3,3,6),cijwn(3,3,6)
      dimension tout(2),tol(9),ymin(9),mark(2)
      dimension stset(10),sset(10)
      dimension iip(6),jjp(6),np(3,3)
      real*8    del(3,3),rot,dudy
      common dudy(3,3),rot(3),enu

      data pi/3.14159265/
      data iip/1,2,3,1,1,2/,jjp/1,2,3,2,3,3/
      data del/1.,0.,0.,0.,1.,0.,0.,0.,1./
      data np/1,4,5,4,2,6,5,6,3/
      
      data ak1/1./
c
      read (5,*) jmax,kmax,sfix
      read (5,*) ((dudy(i,j),i=1,3),j=1,3)
      read (5,*) (rot(i),i=1,3)
      read (5,*) dtout,tmax,tolf,tolmin
      isetm=1
      if (sfix<0) then
        read (5,*) isetm, (stset(i),sset(i),i=1,isetm)
        sfix=sset(1)
        iset=1
      endif
      ss=0
      do 1 i=1,3
      do 1 j=1,3
         sij=0.5*(dudy(i,j)+dudy(j,i))
         ss=ss+sij*sij
c        write (6,*) i,j,sij,ss
 1       continue
      ss=sqrt(2.*ss)
      enu=0
      if (sfix .gt. 0)   enu=ss/(2.*sfix*ak1**2)
 101   format (9hjmax,kmax,2i4,5h sfix,f6.1,7h dUidxj,9f6.1,
     *    4h rot,3f6.1)
      write (6,101) jmax,kmax,sfix,((dudy(i,j),i=1,3),j=1,3)
     * , (rot(i),i=1,3)
      if (jmax*kmax*2.gt.idim) then
        write (6,*) 'redimension'
        stop
      endif
c-------------------------------------------
c     runge parameters
      do 151 n=1,9
      ymin(n)=tolmin
 151  tol(n)=tolf
      mark(1)=1
      mark(2)=2
      itmax=tmax/dtout+1
c------------------------------------------
c      set k array
      dth=pi/jmax
      dacos=1./kmax
c      write (6,*) 'dth',dth,dacos
      i=0
      do 7 k=1,kmax
      do 7 m=1,2
      do 7 j=1,jmax
      i=i+1
      ak(i,3)=dacos*(k-.5)
      if (m.eq.2) ak(i,3)=-ak(i,3)
      sina=sqrt(1.-ak(i,3)**2)
      theta=dth*(j-.5)
      ak(i,1)=sina*dcos(theta)
      ak(i,2)=sina*dsin(theta)
      do 8 n=1,3
 8       ak(i,n)=ak(i,n)*ak1
c      write (6,'(3g12.5)')ak(i,1),ak(i,2),ak(i,3)
 7    continue
      imax=i
c     phij 7,8,9 are kx,ky,kz
      do 9 i=1,imax
      do 9 n=1,3
 9       phij(n+6,i)=ak(i,n)
c------------------------------------------------
c      initialize phij isotropic turbulence
c     phij=3/(8 pi k^4)(delij k^2 - ki kj)Et(k) townsend 3.10.4
c     townsend Et = E11 where E integrate to K
      do 3 i=1,imax
       do 4 n=1,6
           phij(n,i)= del(iip(n),jjp(n))
     *                  -ak(i,iip(n))*ak(i,jjp(n))/ak1**2
 4      continue
c      write (6,'(i2,5g11.4)') i,akk,phij(1,i),phij(2,i),
c     *              phij(3,i),phij(4,i)
 3    continue
c------------------------------------------
c          integrate one dtout at a time
      t=0.
      do 152 itt=1,itmax
      if (itt.gt.1) then
        do 153 i=1,imax
        x=t
        tout(1)=x
        tout(2)=x+dtout
        call runge(9,x,phij(1,i),dphij(1,i),
     *                   tol,ymin,dtout*.01,tout,mark)
 153    continue
        t=x
      endif
c...................... print
c      continuity check
      do 114 n=1,3
 114   con(n)=0.
c      write (6,*) 't',t
      do 115 i=1,imax
      do 113 n=1,3
          cont(n)=0.
      do 113 j=1,3
 113     cont(n)=cont(n)+phij(j+6,i)*phij(np(n,j),i)
c      write (6,*) 'i,c',i,(cont(n),n=1,3)
c      write (6,*) 'p',(phij(n,i),n=1,9)
      do 115 n=1,3
 115    con(n)=con(n)+cont(n) 
c      sum and print
      do 111 n=1,6
        pas(n)=pa(n)
        pa(n)=0.
        dpa(n)=0.
        dispt(n)=0.
      do 112 m=1,6
 112     pkk(n,m)=0.
      do 111 i=1,imax
        akk2=phij(7,i)**2+phij(8,i)**2+phij(9,i)**2
        dispt(n)=dispt(n)+2.*akk2*phij(n,i)
        dpa(n)=dpa(n)+dphij(n,i)
        pa(n)=pa(n)+phij(n,i)
        ii=iip(n)
        jj=jjp(n)
      do 111 m=1,6
         pkk(n,m)=pkk(n,m)+phij(m,i)*phij(ii+6,i)*phij(jj+6,i)/akk2
 111     continue
      tkeold=tke
      tke=.5*(pa(1)+pa(2)+pa(3))
      eps=.5*enu*(dispt(1)+dispt(2)+dispt(3))
      b11=pa(1)/(2.*tke)-1./3. 
      b22=pa(2)/(2.*tke)-1./3.
      b33=pa(3)/ (2.*tke)-1./3.
      b12=pa(4)/(2.*tke)
      b13=pa(5)/(2.*tke)
      b23=pa(6)/(2.*tke)
c       convert pkk to 4 index array 
      do 116 i=1,3
      do 116 j=1,3
      do 116 n=1,3
      do 116 m=1,3
 116     phijknkm(i,j,n,m)=pkk(np(i,j),np(m,n))
c     sum to give em(i,j,n,m)  FRSM  Eq. 9.4.4 
      do 117 i=1,3
      do 117 j=1,3
      do 117 n=1,3
      do 117 m=1,3
 117     em(i,j,n,m)= (phijknkm(i,n,m,j)+ phijknkm(j,n,i,m))/tke
c      form coefficients of strain and vorticity in pressure-strain FRSM eq 9.5.1
      do 118 i=1,3
      do 118 j=1,3
      cijsn(i,j,1)=2*em(i,j,1,1)
      cijsn(i,j,2)=2*em(i,j,2,2)
      cijsn(i,j,3)=2*em(i,j,3,3)
      cijsn(i,j,4)=2*(em(i,j,1,2)+em(i,j,2,1))
      cijsn(i,j,5)=2*(em(i,j,1,3)+em(i,j,3,1))
      cijsn(i,j,6)=2*(em(i,j,2,3)+em(i,j,3,2))
      cijwn(i,j,4)=2*(em(i,j,1,2)-em(i,j,2,1))
      cijwn(i,j,5)=2*(em(i,j,1,3)-em(i,j,3,1))
      cijwn(i,j,6)=2*(em(i,j,2,3)-em(i,j,3,2))
 118  continue

      if (t.eq.0) then
        tke0=tke
        write (6,*)  88
        write (6,*) 't beta kdk0 s'
        write (6,*) 'con1 con2 con3  DkDt dkdt dkdtP'
        write (6,*) ' uu vv ww   uv  uw vw '
        write (6,*) ' duu dvv dww   duv  duw dvw '
        write (6,*) ' b11 b22 b33  b12 b13 b23'
        write (6,*) ' dispuu dispvv dispww dispuv dispuw dispvw '
        write (6,*) ' cs1111 cs1122 cs1133 cs1112 cs1113 cs1123',
     *   '  cw1112 cw1113 cw1123'
        write (6,*) ' cs2211 cs2222 cs2233 cs2212 cs2213 cs2223',
     *   '  cw2212 cw2213 cw2223'
        write (6,*) ' cs3311 cs3322 cs3333 cs3312 cs3313 cs3323',
     *   '  cw3312 cw3313 cw3323'
        write (6,*) ' cs1211 cs1222 cs1233 cs1212 cs1213 cs1223',
     *   '  cw1212 cw1213 cw1223'
        write (6,*) ' cs1311 cs1322 cs1333 cs1312 cs1313 cs1323',
     *   '  cw1312 cw1313 cw1323'
        write (6,*) ' cs2311 cs2322 cs2333 cs2312 cs2313 cs2323',
     *   '  cw2312 cw2313 cw2323'
        dkdt1=0.
        dkdt2=0.
        dkdy3=0.
      else
        dkdt1=(tke-tkeold)/dtout
        dkdt2=.5*(dpa(1)+dpa(2)+dpa(3))
        dkdt3=-2.*tke*(b11*dudy(1,1)+b22*dudy(2,2)+b33*dudy(3,3)
     *         +b12*(dudy(1,2)+dudy(2,1))+b13*(dudy(1,3)+dudy(3,1))
     *         +b23*(dudy(3,2)+dudy(2,3)) )
      endif
      ssd=ss
      if (eps.gt.0.) then
           ssd=ss*tke/eps
           if (isetm>1) then
              if (t*ss.gt.stset(iset+1)) iset=min0(isetm-1,iset+1)
              sfix=sset(iset)+(sset(iset+1)-sset(iset))*
     *             (t*ss-stset(iset))/(stset(iset+1)-stset(iset))
           endif
           enu=enu*ssd/sfix
      endif
      write (6,'(/6g12.5)') t,t*ss, tke/tke0,ssd
      write (6,'(6g12.5)')(con(n),n=1,3),dkdt1,dkdt2,dkdt3
      write (6,'(6g12.5)')(pa(n),n=1,6),(dpa(n),n=1,6)
      write (6,'(6g12.5)') b11,b22,b33,b12,b13,b23
      write (6,'(6g12.5)') (dispt(n)/tke,n=1,6)
      write (6,'(9f8.3)') ((cijsn(iip(n),jjp(n),m)
     *    ,m=1,6),(cijwn(iip(n),jjp(n),m),m=4,6),n=1,6)
 152  continue
      stop
      end
c------------------------------------------------------
      subroutine diffeq(nnn,t,p,dp)
      real*8 p(9),dp(9),dudy,rot,bk(3),e(3,3,3)
      dimension np(3,3),iip(6),jjp(6)
      common dudy(3,3),rot(3),enu

      data e/0.,0.,0., 0.,0.,-1., 0.,1.,0.,
     *       0.,0.,1., 0.,0.,0., -1.,0.,0.,
     *       0.,-1.,0., 1.,0.,0., 0.,0.,0./  
      data iip/1,2,3,1,1,2/,jjp/1,2,3,2,3,3/
      data np/1,4,5,4,2,6,5,6,3/
c---------------------------------------------------
      do 1 i=1,3
 1       bk(i)=p(i+6)
c    
      akk2=bk(1)*bk(1)+bk(2)*bk(2)+bk(3)*bk(3)
      do 12 n=1,6
      ii=iip(n)
      jj=jjp(n)
      dp(n)=-2*enu*akk2*p(n)
      do 12 m=1,3
      dp(n)=dp(n)-p(np(m,jj))*dudy(ii,m)-p(np(m,ii))*dudy(jj,m)
      do 12 k=1,3
      dp(n)=dp(n)+p(np(m,jj))*(-2.*rot(k)*e(ii,k,m)
     *      +2.*dudy(k,m)*bk(ii)*bk(k)/akk2)
      dp(n)=dp(n)+p(np(m,ii))*(-2.*rot(k)*e(jj,k,m)
     *      +2.*dudy(k,m)*bk(jj)*bk(k)/akk2)
      do 12 nn=1,3
         dp(n)=dp(n)+2.*p(np(m,jj))*e(nn,k,m)*rot(k)*bk(ii)*bk(nn)/akk2
     *        +2.*p(np(m,ii))*e(nn,k,m)*rot(k)*bk(jj)*bk(nn)/akk2
 12   continue
      do 13 m=1,3
         dp(m+6)=0.
      do 13 k=1,3
 13      dp(m+6)=dp(m+6)-bk(k)*dudy(k,m)
c      write (6,'(a2,7g11.4)') 'd2',t,bk(1),bk(2),bk(3),p(1),p(4),dp(1)
      return
      end
c--------------------------------------------------
      subroutine print(n,xout,y,dy,j)
      real*8 y(2),dy(2)
      return
      end
c---------------------------------------------------
C                                                RUNGE
C      FIRST ORDER DIFF. EQ. ROUTINE -- ADJUSTS STEP SIZE
      SUBROUTINE RUNGE (N,X,Y,dy,TOL,YMIN,H,XOUT,MARK)
      real*8 y(2),dy(2)
      dimension tol(2),ymin(2),xout(2),mark(2)
      real*8 YA(50),FA(50),FB(50),FC(50),YKEEP(50)
      dimension sub(50)
      if (n .gt. 50) then
       write (6,*) 'n=',n, 'too large redimension runge'
       stop
      endif

      KBTWN=1
      KBIG=1
      KLOW=1
      NCOUNT=15
      J=MARK(1)
      MAX=MARK(2)
 230   DO 250 I=1,N
 250   SUB(I)=TOL(I)/32.0
 10    IF (MAX-J)20,30,30
20     RETURN
C
 30    A=XOUT(J)-X
       B=ABS(2.E-6*X)
      IF (A+B) 40,35,35
 35    IF(A-B) 50,50,60
 40   J=J+1
      GO TO 10
 50   CALL PRINT (N,XOUT,Y,DY,J)
      J=J+1
      GO TO 10
 60    IF (A-1.5*H) 70 ,70,80
 70    H=A
      GO TO 1000
 80    IF (A-3.*H) 90,1000,1000
 90    H=.5*A
C
C          DO RUNGE-KUTTE-MERSON INTEGRATION
C
 1000 XA=X+H/3.
      XB=X+.5*H
      CALL DIFFEQ (N,X,Y,DY)
      X=X+H
      DO 1030 I=1,N
      YKEEP(I)=Y(I)
      FA(I)=H*DY(I)
 1030 YA(I)=Y(I)+FA(I)/3.
      CALL DIFFEQ (N,XA,YA,DY)
      DO 1040 I=1,N
 1040 YA(I)=Y(I)+FA(I)/6.+H*DY(I)/6.
      CALL DIFFEQ (N,XA,YA,DY)
      DO 1050 I=1,N
      FB(I)=H*DY(I)
 1050  YA(I)=Y(I)+.125*FA(I)+.375*FB(I)
      CALL DIFFEQ (N,XB,YA,DY)
      DO 1060 I=1,N
      FC(I)=H*DY(I)
 1060 YA(I)=Y(I)+.5*FA(I)-1.5*FB(I)+2.*FC(I)
      CALL DIFFEQ (N,X,YA,DY)
      DO 1130 I=1,N
      Y(I)=Y(I)+FA(I)/6.+2./3.*FC(I)+H*DY(I)/6.
 1061 U=Y(I)
      IF (ABS(U)-YMIN(I)) 1130,1130,1090
 1090 KLOW=2
      E=.2*ABS(U-YA(I))
      IF (E-ABS(TOL(I)*U)) 1110,1110,1100
 1100 KBIG=2
      GO TO 1130
 1110 IF (E-ABS(SUB(I)*U)) 1130,1120,1120
 1120 KBTWN=2
 1130  CONTINUE
      GO TO (100,1135),KLOW
 1135 GO TO (1180,1140),KBIG
 1140 NCOUNT=NCOUNT-1
      IF (NCOUNT) 1150,1150,1170
 1150 WRITE (6,1160) X,H
      WRITE (6,1165) (I,Y(I),DY(I),I=1,N)
      RETURN
 1160  FORMAT (' STEP SIZE HALVED 15 TIMES CONSECUTIVELY SINCE LAST PRIN
     1T ',/' PROGRAM TERMINATED AT X= ',E16.8,', H=',E16.8,//
     2 3H  I ,13X,4HY(I),16X,5HDY(I),//)
 1165 FORMAT (I3,7X,E16.8,4X,E16.8)
 1170 KBIG=1
      IF (H-B) 1176,1172,1172
 1172 X=X-H
      H=.5*H
      DO 1174 I=1,N
 1174 Y(I)=YKEEP(I)
      KBTWN=1
      KLOW=1
      GO TO 1000
 1176 M=15-NCOUNT
      WRITE (6,1178) M,X,H
      WRITE (6,1165) (I,Y(I),DY(I),I=1,N)
      RETURN
 1178 FORMAT (' STEP SIZE BECAME TOO SMALL FOR COMPUTER.'
     1 ,/' IT HAS BEEN HALVED',I2,' TIMES CONSECUTIVELY.',/' PROGRAM TER
     2MINATED AT X =',E16.8,'  H= ',E16.8,//3H  I,13X,4HY(I),16X,5HDY(I)
     3,//)
 1180 NCOUNT=15
      GO TO (1190,1200), KBTWN
 1190 H=2.*H
 1200 KBTWN=1
      KLOW=1
 100  GO TO 10
      END
