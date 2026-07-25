!!***this program to use data and interpolate them via myself
!for JAN!
program Datausage

     integer ,parameter::m=288,n=145
     real*4 aa(m,n),per(m,n),mean(m,n,10),Z(m,n,12),Zr(m,n,10),F(12),FE(10),FD(10)
     character(len=2):: hour(4),duration(4)
     character(len=3):: vlevr(1:10),vlev(1:12)
	 character(len=4):: year(1:4)
     character(len=len("HrXXXXXXXXX.dat")):: out_file
     character(len=len("HXXXXXXXXX.dat")):: in_file
     character(len=len("mHXXXXXXXXX.dat")):: outm_file
     character(len=len("pHXXXXXXXXX.dat"))::outp_file
     character(len=87):: location_r, location_i
     !character(2):: location2
     integer:: DD(3),irec,t,ll,i,j,ioerr,tmax(4),c

     DATA DD/31,31,28/
     DATA vlev/'000','925','850','700','600','500','400','300','250','200','150','100'/
     DATA vlevr/'000','900','800','700','600','500','400','300','200','100'/
     !DATA hour/'00','06','12','18'/
     DATA year/'2016','2017','2017','1617'/
     DATA duration/'12','01','02','wi'/
     DATA tmax/124,124,112,360/
     undef=1e+15
	
     location_r = '/home/pari/Documents/university/master/THESIS/MYWORK-summer2018/3.flux-code/data-10lev/'
     location_i = '/home/pari/Documents/university/master/THESIS/MYWORK-summer2018/2.(12to10lev)code/'

      open(90,file="Z.txt" ,STATUS='unknown',Action="readwrite",iostat=ioerr)
  			                  	if (ioerr /= 0) then
					                      !write(*,*) "ioerrtxt=", ioerr
				                     	stop
					               else
                              !write(*,*) "file is ready"
                                end if
     mean(:,:,:)=0
	  c=4
     DO while (c<5)

 	    irec=0
   		   DO t=1,tmax(c)
                    irec=irec+1
                    DO ll=1,12
				                  100 format(A,A,A,A,A)
				                  write(in_file,100) "H", vlev(ll), year(c), duration(c), ".dat"
				                  !write(*,*) in_file
				                  !ioerr=200
				                  !write(*,*) "ioerr",ioerr
                            	open(20,file= in_file, form="unformatted"&
                              ,STATUS='old' ,ACCESS='direct',Action="read",iostat=ioerr,RECL=4*m*n)
  			                  	if (ioerr /= 0) then
					                      !write(*,*) "ioerr20=", ioerr
				                     	!stop
					               end if

			                   write(*,*) "t= " ,t

                            read(20, rec=t)  aa
                            Z(:,:,ll)= aa(:,:)
                            close(20)
                     end do

          do j=1,n
               do i=1,m
                    do ll=1,12
!write(*,*) Z(i,j,ll), i,j, ll
                         F(13-ll)=Z(i,j,ll)
                    end do
                    call VSPN(F,FE)
                    do ll= 1,10
                         Zr(i,j,ll)=FE(11-ll)
!write(*,*) Zr(i,j,ll), i,j, ll
                    end do
               end do
          end do

          do ll=1,10
               write(out_file, 100) "Hr", vlevr(ll),year(c),duration(c), ".dat"
			   !write(*,*) out_file
		!!!write(location2,2000) location, duration(counter) ,'/'
               open(21,file= location_r//out_file, form="unformatted"&
                ,STATUS='unknown',ACCESS='direct',Action="readwrite",iostat=ieorr,RECL=4*m*n)  !out10lev
				!write(*,*) "ioerr", ioerr
               do j=1,n
                    do i=1,m
                         mean(i,j,ll)=mean(i,j,ll) + Zr(i,j,ll)
                         aa(i,j)=Zr(i,j,ll)
                    end do
               end do
               write (21, rec= irec)  aa
          end do
     end do   !! end of t time

     mean=mean/(tmax(c)*1.0)
     do ll=1,10
          write(outp_file, 100) "pH", vlevr(ll),year(c),duration(c), ".dat"
		    write(outm_file, 100) "mH", vlevr(ll),year(c),duration(c), ".dat"
          !open(22,file=location_r//outm_file , form="unformatted", &
          !STATUS='new',ACCESS='direct',Action="write",RECL=4*m*n)   !!file of mean
          !open(23,file=location_r//outp_file , form="unformatted", &
          !STATUS='new',ACCESS='direct',Action="write",RECL=4*m*n)   !!file of per
          irec=0
          do t=1,tmax(c)
               irec = irec + 1
               read (21, rec= t)  aa
               per(:,:)= aa(:,:) - mean(:,:,ll)
               !write (23, rec= irec)  per
do i=1,m
do j=1,n
if (t==100 .or. t==200) then
write(90,*) "t,i,j,k",t,i,j,ll, per(i,j) 
write(90,*) "t,i,j,k",t,i,j,ll, per(i,j) 
write(*,*) per(i,j) , i,j, ll
end if
end do
end do
          end do
          close(23)

          aa(:,:)= mean(:,:,ll)
          !write (22, rec=1)  aa
     end do
     close(22)
	c=c+1
	end do

end program

!*********VPSN subroutin
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE VSPN(F,FE)
	PARAMETER (N=12,NE=10)
      REAL X(N),F(N),D(N),WK(2*N),FE(NE),XE(NE),FD(NE)
      LOGICAL SPLINE
      DATA X /100.0,150.0,200.0,250.0,300.0,400.0,500.0,600.0,700.0,850,925.0,1000.0/
      NWK=2*N
	XE(1)=100.0  
      DO  5  I=2,NE
       XE(I)=XE(I-1)+100.0
 5    CONTINUE

      SPLINE=.FALSE.
      CALL PCHEZ(N,X,F,D,SPLINE,WK,NWK,IERR)
      IF(IERR.LT.0)THEN
         WRITE(*,*)' AN ERROR CALLING PCHEZ, IERR=',IERR
         STOP
      ENDIF

      CALL PCHEV(N,X,F,D,NE,XE,FE,FD,IERR)
      IF(IERR.NE.0)THEN
         WRITE(*,*) ' AN ERROR CALLING PCHEV, IERR=',IERR
         STOP
      ENDIF
      RETURN
      END
      SUBROUTINE PCHEZ(N,X,F,D,SPLINE,WK,LWK,IERR)
      INTEGER  N, LWK, IERR
      REAL  X(N), F(N), D(N), WK(LWK)
      LOGICAL SPLINE
      INTEGER IC(2), INCFD
      REAL  VC(2)
      DATA IC(1) /0/
      DATA IC(2) /0/
      DATA INCFD /1/
      IF ( SPLINE ) THEN
        CALL  PCHSP (IC, VC, N, X, F, D, INCFD, WK, LWK, IERR)
      ELSE
        CALL  PCHIM (N, X, F, D, INCFD, IERR)
      ENDIF
      RETURN
      END
      SUBROUTINE PCHIM(N,X,F,D,INCFD,IERR)
      INTEGER  N, INCFD, IERR
      REAL  X(N), F(INCFD,N), D(INCFD,N)
      INTEGER  I, NLESS1
      REAL  DEL1, DEL2, DMAX, DMIN, DRAT1, DRAT2, DSAVE, H1, H2, HSUM, HSUMT3, THREE, W1, W2, ZERO
      REAL  PCHST
      DATA  ZERO /0./,  THREE /3./
      IF ( N.LT.2 )  GO TO 5001
      IF ( INCFD.LT.1 )  GO TO 5002
      DO 1  I = 2, N
         IF ( X(I).LE.X(I-1) )  GO TO 5003
    1 CONTINUE
      IERR = 0
      NLESS1 = N - 1
      H1 = X(2) - X(1)
      DEL1 = (F(1,2) - F(1,1))/H1
      DSAVE = DEL1
      IF (NLESS1 .GT. 1)  GO TO 10
      D(1,1) = DEL1
      D(1,N) = DEL1
      GO TO 5000
   10 CONTINUE
      H2 = X(3) - X(2)
      DEL2 = (F(1,3) - F(1,2))/H2
      HSUM = H1 + H2
      W1 = (H1 + HSUM)/HSUM
      W2 = -H1/HSUM
      D(1,1) = W1*DEL1 + W2*DEL2
      IF ( PCHST(D(1,1),DEL1) .LE. ZERO)  THEN
         D(1,1) = ZERO
      ELSE IF ( PCHST(DEL1,DEL2) .LT. ZERO)  THEN
         DMAX = THREE*DEL1
         IF (ABS(D(1,1)) .GT. ABS(DMAX))  D(1,1) = DMAX
      ENDIF
      DO 50  I = 2, NLESS1
         IF (I .EQ. 2)  GO TO 40
         H1 = H2
         H2 = X(I+1) - X(I)
         HSUM = H1 + H2
         DEL1 = DEL2
         DEL2 = (F(1,I+1) - F(1,I))/H2
   40    CONTINUE
         D(1,I) = ZERO
         IF ( PCHST(DEL1,DEL2) )  42, 41, 45
   41    CONTINUE
         IF (DEL2 .EQ. ZERO)  GO TO 50
         IF ( PCHST(DSAVE,DEL2) .LT. ZERO)  IERR = IERR + 1
         DSAVE = DEL2
         GO TO 50
   42    CONTINUE
         IERR = IERR + 1
         DSAVE = DEL2
         GO TO 50
   45    CONTINUE
         HSUMT3 = HSUM+HSUM+HSUM
         W1 = (HSUM + H1)/HSUMT3
         W2 = (HSUM + H2)/HSUMT3
         DMAX = AMAX1( ABS(DEL1), ABS(DEL2) )
         DMIN = AMIN1( ABS(DEL1), ABS(DEL2) )
         DRAT1 = DEL1/DMAX
         DRAT2 = DEL2/DMAX
         D(1,I) = DMIN/(W1*DRAT1 + W2*DRAT2)
   50 CONTINUE
      W1 = -H2/HSUM
      W2 = (H2 + HSUM)/HSUM
      D(1,N) = W1*DEL1 + W2*DEL2
      IF ( PCHST(D(1,N),DEL2) .LE. ZERO)  THEN
         D(1,N) = ZERO
      ELSE IF ( PCHST(DEL1,DEL2) .LT. ZERO)  THEN
         DMAX = THREE*DEL2
         IF (ABS(D(1,N)) .GT. ABS(DMAX))  D(1,N) = DMAX
      ENDIF
 5000 CONTINUE
      RETURN
 5001 CONTINUE
      IERR = -1
      CALL XERROR ('PCHIM -- NUMBER OF DATA POINTS LESS THAN TWO', 44, IERR, 1)
      RETURN
 5002 CONTINUE
      IERR = -2
      CALL XERROR ('PCHIM -- INCREMENT LESS THAN ONE' , 32, IERR, 1)
      RETURN
 5003 CONTINUE
      IERR = -3
      CALL XERROR ('PCHIM -- X-ARRAY NOT STRICTLY INCREASING', 40, IERR, 1)
      RETURN
      END
      REAL FUNCTION PCHST(ARG1,ARG2)
      REAL  ARG1, ARG2
      REAL  ONE, ZERO
      DATA  ZERO /0./,  ONE /1./
      PCHST = SIGN(ONE,ARG1) * SIGN(ONE,ARG2)
      IF ((ARG1.EQ.ZERO) .OR. (ARG2.EQ.ZERO))  PCHST = ZERO
      RETURN
      END
      SUBROUTINE PCHSP(IC,VC,N,X,F,D,INCFD,WK,NWK,IERR)
      INTEGER  IC(2), N, INCFD, NWK, IERR
      REAL  VC(2), X(N), F(INCFD,N), D(INCFD,N), WK(2,N)
      INTEGER  IBEG, IEND, INDEX, J, NM1
      REAL  G, HALF, ONE, STEMP(3), THREE, TWO, XTEMP(4), ZERO
      REAL  PCHDF
      DATA  ZERO /0./,  HALF /0.5/,  ONE /1./,  TWO /2./,  THREE /3./
      IF ( N.LT.2 )  GO TO 5001
      IF ( INCFD.LT.1 )  GO TO 5002
      DO 1  J = 2, N
         IF ( X(J).LE.X(J-1) )  GO TO 5003
    1 CONTINUE
      IBEG = IC(1)
      IEND = IC(2)
      IERR = 0
      IF ( (IBEG.LT.0).OR.(IBEG.GT.4) )  IERR = IERR - 1
      IF ( (IEND.LT.0).OR.(IEND.GT.4) )  IERR = IERR - 2
      IF ( IERR.LT.0 )  GO TO 5004
      IF ( NWK .LT. 2*N )  GO TO 5007
      DO 5  J=2,N
         WK(1,J) = X(J) - X(J-1)
         WK(2,J) = (F(1,J) - F(1,J-1))/WK(1,J)
    5 CONTINUE
      IF ( IBEG.GT.N )  IBEG = 0
      IF ( IEND.GT.N )  IEND = 0
      IF ( (IBEG.EQ.1).OR.(IBEG.EQ.2) )  THEN
         D(1,1) = VC(1)
      ELSE IF (IBEG .GT. 2)  THEN
         DO 10  J = 1, IBEG
            INDEX = IBEG-J+1
            XTEMP(J) = X(INDEX)
            IF (J .LT. IBEG)  STEMP(J) = WK(2,INDEX)
   10    CONTINUE
         D(1,1) = PCHDF (IBEG, XTEMP, STEMP, IERR)
         IF (IERR .NE. 0)  GO TO 5009
         IBEG = 1
      ENDIF
      IF ( (IEND.EQ.1).OR.(IEND.EQ.2) )  THEN
         D(1,N) = VC(2)
      ELSE IF (IEND .GT. 2)  THEN
         DO 15  J = 1, IEND
            INDEX = N-IEND+J
            XTEMP(J) = X(INDEX)
            IF (J .LT. IEND)  STEMP(J) = WK(2,INDEX+1)
   15    CONTINUE
         D(1,N) = PCHDF (IEND, XTEMP, STEMP, IERR)
         IF (IERR .NE. 0)  GO TO 5009
         IEND = 1
      ENDIF
      IF (IBEG .EQ. 0)  THEN
         IF (N .EQ. 2)  THEN
            WK(2,1) = ONE
            WK(1,1) = ONE
            D(1,1) = TWO*WK(2,2)
         ELSE
            WK(2,1) = WK(1,3)
            WK(1,1) = WK(1,2) + WK(1,3)
            D(1,1) =((WK(1,2) + TWO*WK(1,1))*WK(2,2)*WK(1,3)  + WK(1,2)**2*WK(2,3)) / WK(1,1)
         ENDIF
      ELSE IF (IBEG .EQ. 1)  THEN
         WK(2,1) = ONE
         WK(1,1) = ZERO
      ELSE
         WK(2,1) = TWO
         WK(1,1) = ONE
         D(1,1) = THREE*WK(2,2) - HALF*WK(1,2)*D(1,1)
      ENDIF
      NM1 = N-1
      IF (NM1 .GT. 1)  THEN
         DO 20 J=2,NM1
            IF (WK(2,J-1) .EQ. ZERO)  GO TO 5008
            G = -WK(1,J+1)/WK(2,J-1)
            D(1,J) = G*D(1,J-1) + THREE*(WK(1,J)*WK(2,J+1) + WK(1,J+1)*WK(2,J))
            WK(2,J) = G*WK(1,J-1) + TWO*(WK(1,J) + WK(1,J+1))
   20    CONTINUE
      ENDIF
      IF (IEND .EQ. 1)  GO TO 30
      IF (IEND .EQ. 0)  THEN
         IF (N.EQ.2 .AND. IBEG.EQ.0)  THEN
            D(1,2) = WK(2,2)
            GO TO 30
         ELSE IF ((N.EQ.2) .OR. (N.EQ.3 .AND. IBEG.EQ.0))  THEN
            D(1,N) = TWO*WK(2,N)
            WK(2,N) = ONE
            IF (WK(2,N-1) .EQ. ZERO)  GO TO 5008
            G = -ONE/WK(2,N-1)
         ELSE
            G = WK(1,N-1) + WK(1,N)
            D(1,N) = ((WK(1,N)+TWO*G)*WK(2,N)*WK(1,N-1)&
 + WK(1,N)**2*(F(1,N-1)-F(1,N-2))/WK(1,N-1))/G
            IF (WK(2,N-1) .EQ. ZERO)  GO TO 5008
            G = -G/WK(2,N-1)
            WK(2,N) = WK(1,N-1)
         ENDIF
      ELSE
         D(1,N) = THREE*WK(2,N) + HALF*WK(1,N)*D(1,N)
         WK(2,N) = TWO
         IF (WK(2,N-1) .EQ. ZERO)  GO TO 5008
         G = -ONE/WK(2,N-1)
      ENDIF
      WK(2,N) = G*WK(1,N-1) + WK(2,N)
      IF (WK(2,N) .EQ. ZERO)   GO TO 5008
      D(1,N) = (G*D(1,N-1) + D(1,N))/WK(2,N)
   30 CONTINUE
      DO 40 J=NM1,1,-1
         IF (WK(2,J) .EQ. ZERO)  GO TO 5008
         D(1,J) = (D(1,J) - WK(1,J)*D(1,J+1))/WK(2,J)
   40 CONTINUE
      RETURN
 5001 CONTINUE
      IERR = -1
      CALL XERROR ('PCHSP -- NUMBER OF DATA POINTS LESS THAN TWO', 44, IERR, 1)
      RETURN
 5002 CONTINUE
      IERR = -2
      CALL XERROR ('PCHSP -- INCREMENT LESS THAN ONE', 32, IERR, 1)
      RETURN
 5003 CONTINUE
      IERR = -3
      CALL XERROR ('PCHSP -- X-ARRAY NOT STRICTLY INCREASING'  , 40, IERR, 1)
      RETURN
 5004 CONTINUE
      IERR = IERR - 3
      CALL XERROR ('PCHSP -- IC OUT OF RANGE' , 24, IERR, 1)
      RETURN
 5007 CONTINUE
      IERR = -7
      CALL XERROR ('PCHSP -- WORK ARRAY TOO SMALL' , 29, IERR, 1)
      RETURN
 5008 CONTINUE
      IERR = -8
      CALL XERROR ('PCHSP -- SINGULAR LINEAR SYSTEM' , 31, IERR, 1)
      RETURN
 5009 CONTINUE
      IERR = -9
      CALL XERROR ('PCHSP -- ERROR RETURN FROM PCHDF' , 32, IERR, 1)
      RETURN
      END
      REAL FUNCTION PCHDF(K,X,S,IERR)
      INTEGER  K, IERR
      REAL  X(K), S(K)
      INTEGER  I, J
      REAL  VALUE, ZERO
      DATA  ZERO /0./
      IF (K .LT. 3)  GO TO 5001
      DO 10  J = 2, K-1
         DO 9  I = 1, K-J
            S(I) = (S(I+1)-S(I))/(X(I+J)-X(I))
    9    CONTINUE
   10 CONTINUE
      VALUE = S(1)
      DO 20  I = 2, K-1
         VALUE = S(I) + VALUE*(X(K)-X(I))
   20 CONTINUE
      IERR = 0
      PCHDF = VALUE
      RETURN
 5001 CONTINUE
      IERR = -1
      CALL XERROR ('PCHDF -- K LESS THAN THREE' , 26, IERR, 1)
      PCHDF = ZERO
      RETURN
      END
      SUBROUTINE PCHEV(N,X,F,D,NVAL,XVAL,FVAL,DVAL,IERR)
      INTEGER  N, NVAL, IERR
      REAL  X(N), F(N), D(N), XVAL(NVAL), FVAL(NVAL), DVAL(NVAL)
      INTEGER INCFD
      LOGICAL SKIP
      DATA SKIP /.TRUE./
      DATA INCFD /1/
      CALL PCHFD(N,X,F,D,INCFD,SKIP,NVAL,XVAL,FVAL,DVAL,IERR)
 5000 CONTINUE
      RETURN
      END
      SUBROUTINE PCHFD(N,X,F,D,INCFD,SKIP,NE,XE,FE,DE,IERR)
      INTEGER  N, INCFD, NE, IERR
      REAL  X(N), F(INCFD,N), D(INCFD,N), XE(NE), FE(NE), DE(NE)
      LOGICAL  SKIP
      INTEGER  I, IERC, IR, J, JFIRST, NEXT(2), NJ
      IF (SKIP)  GO TO 5
      IF ( N.LT.2 )  GO TO 5001
      IF ( INCFD.LT.1 )  GO TO 5002
      DO 1  I = 2, N
         IF ( X(I).LE.X(I-1) )  GO TO 5003
    1 CONTINUE
    5 CONTINUE
      IF ( NE.LT.1 )  GO TO 5004
      IERR = 0
      SKIP = .TRUE.
      JFIRST = 1
      IR = 2
   10 CONTINUE
         IF (JFIRST .GT. NE)  GO TO 5000
         DO 20  J = JFIRST, NE
            IF (XE(J) .GE. X(IR))  GO TO 30
   20    CONTINUE
         J = NE + 1
         GO TO 40
   30    CONTINUE
         IF (IR .EQ. N)  J = NE + 1
   40    CONTINUE
         NJ = J - JFIRST
         IF (NJ .EQ. 0)  GO TO 50
        CALL CHFDV (X(IR-1),X(IR), F(1,IR-1),F(1,IR), D(1,IR-1),D(1,IR),&
 NJ, XE(JFIRST), FE(JFIRST), DE(JFIRST), NEXT, IERC)
         IF (IERC .LT. 0)  GO TO 5005
         IF (NEXT(2) .EQ. 0)  GO TO 42
            IF (IR .LT. N)  GO TO 41
               IERR = IERR + NEXT(2)
               GO TO 42
   41       CONTINUE
               GO TO 5005
   42    CONTINUE
         IF (NEXT(1) .EQ. 0)  GO TO 49
            IF (IR .GT. 2)  GO TO 43
               IERR = IERR + NEXT(1)
               GO TO 49
   43       CONTINUE
               DO 44  I = JFIRST, J-1
                  IF (XE(I) .LT. X(IR-1))  GO TO 45
   44          CONTINUE
               GO TO 5005
   45          CONTINUE
               J = I
               DO 46  I = 1, IR-1
                  IF (XE(J) .LT. X(I)) GO TO 47
   46          CONTINUE
   47          CONTINUE
               IR = MAX0(1, I-1)
   49    CONTINUE
         JFIRST = J
   50 CONTINUE
      IR = IR + 1
      IF (IR .LE. N)  GO TO 10
 5000 CONTINUE
      RETURN
 5001 CONTINUE
      IERR = -1
      CALL XERROR ('PCHFD -- NUMBER OF DATA POINTS LESS THAN TWO' , 44, IERR, 1)
      RETURN
 5002 CONTINUE
      IERR = -2
      CALL XERROR ('PCHFD -- INCREMENT LESS THAN ONE', 32, IERR, 1)
      RETURN
 5003 CONTINUE
      IERR = -3
      CALL XERROR ('PCHFD -- X-ARRAY NOT STRICTLY INCREASING', 40, IERR, 1)
      RETURN
 5004 CONTINUE
      IERR = -4
      CALL XERROR ('PCHFD -- NUMBER OF EVALUATION POINTS LESS THAN ONE', 50, IERR, 1)
      RETURN
 5005 CONTINUE
      IERR = -5
      CALL XERROR ('PCHFD -- ERROR RETURN FROM CHFDV -- FATAL' , 41, IERR, 2)
      RETURN
      END
      SUBROUTINE CHFDV(X1,X2,F1,F2,D1,D2,NE,XE,FE,DE,NEXT,IERR)
      INTEGER  NE, NEXT(2), IERR
      REAL  X1, X2, F1, F2, D1, D2, XE(NE), FE(NE), DE(NE)
      INTEGER  I
      REAL  C2, C2T2, C3, C3T3, DEL1, DEL2, DELTA, H, X, XMI, XMA, ZERO
      DATA  ZERO /0./
      IF (NE .LT. 1)  GO TO 5001
      H = X2 - X1
      IF (H .EQ. ZERO)  GO TO 5002
      IERR = 0
      NEXT(1) = 0
      NEXT(2) = 0
      XMI = AMIN1(ZERO, H)
      XMA = AMAX1(ZERO, H)
      DELTA = (F2 - F1)/H
      DEL1 = (D1 - DELTA)/H
      DEL2 = (D2 - DELTA)/H
      C2 = -(DEL1+DEL1 + DEL2)
      C2T2 = C2 + C2
      C3 = (DEL1 + DEL2)/H
      C3T3 = C3+C3+C3
      DO 500  I = 1, NE
         X = XE(I) - X1
         FE(I) = F1 + X*(D1 + X*(C2 + X*C3))
         DE(I) = D1 + X*(C2T2 + X*C3T3)
         IF ( X.LT.XMI )  NEXT(1) = NEXT(1) + 1
         IF ( X.GT.XMA )  NEXT(2) = NEXT(2) + 1
  500 CONTINUE
      RETURN
 5001 CONTINUE
      IERR = -1
      CALL XERROR ('CHFDV -- NUMBER OF EVALUATION POINTS LESS THAN ONE'  , 50, IERR, 1)
      RETURN
 5002 CONTINUE
      IERR = -2
      CALL XERROR ('CHFDV -- INTERVAL ENDPOINTS EQUAL' , 33, IERR, 1)
      RETURN
      END
      SUBROUTINE XERROR(MESSG,NMESSG,NERR,LEVEL)
      CHARACTER*100 MESSG
      CALL XERRWV(MESSG,NMESSG,NERR,LEVEL,0,0,0,0,0.,0.)
      RETURN
      END
      SUBROUTINE XERRWV(MESSG,NMESSG,NERR,LEVEL,NI,I1,I2,NR,R1,R2)
      CHARACTER*100 MESSG
      WRITE(*,*) MESSG
      IF(NI.EQ.2)THEN
        WRITE(*,*) I1,I2
      ELSEIF(NI.EQ.1) THEN
        WRITE(*,*) I1
      ENDIF
      IF(NR.EQ.2) THEN
        WRITE(*,*) R1,R2
      ELSEIF(NR.EQ.1) THEN
        WRITE(*,*) R1
      ENDIF
      IF(ABS(LEVEL).LT.2)RETURN
      STOP
      END
	    

