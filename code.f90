PROGRAM flux_summer2018
! revised by Moheb on 27/12/93 (18/3/2015) to include the computation of 
! area average of wave activity for the Mediterranean domain and its subdomains
! revised by Moheb on 10/1/94 (30/3/2015) to include the computation of
! area average of 3D flux divergence for the Mediterranean domain and its subdomains
! edited by me... 2018

!IMPLICIT NONE
PARAMETER (nlev=10)
REAL*8 rr,omega,two_omega,gg,pi,gas_const,spe_heat,half,two,zero,one,quart,five
REAL*8 f_input1(10),f_input2(10),f_input3(10),f_input4(10),f_input_total,f_input1w(10)
REAL*8 f_input2w(10),f_input3w(10),f_input4w(10),f_input_totalw
REAL*8 f_input1c(10),f_input2c(10),f_input3c(10),f_input4c(10),f_input_totalc,f_input1e(10)
REAL*8 f_input2e(10),f_input3e(10),f_input4e(10),f_input_totale        
REAL*8 f_input1_3D,f_input2_3D,f_input3_3D,f_input4_3D,f_input1w_3D,f_input2w_3D,f_input3w_3D,f_input4w_3D
REAL*8 f_input1c_3D,f_input2c_3D,f_input3c_3D,f_input4c_3D,f_input1e_3D,f_input2e_3D,f_input3e_3D,f_input4e_3D
REAL*8 eamean(nlev),eamean_w(nlev),eamean_c(nlev),eamean_e(nlev)
REAL*8 fdmean(nlev),fdmean_w(nlev),fdmean_c(nlev),fdmean_e(nlev)
REAL*8 sa_t,sa_w,sa_c,sa_e
REAL*8 eamean_3D,eamean_w_3D,eamean_c_3D,eamean_e_3D
REAL*8 fdmean_3D,fdmean_w_3D,fdmean_c_3D,fdmean_e_3D

INTEGER m,n,i,j,k,tt,c,ioerr,irec,nw
PARAMETER (m=288,n=145,nw=m/2)              !m is for long(i), n is for lat (j)
PARAMETER (rr=6371.0D3,omega=7.292D-5,two_omega=2.D0*omega)
PARAMETER (gg=9.81D0,pi=3.14159265358979d0)
PARAMETER (gas_const=287.05D0,spe_heat=1.0035D3)
PARAMETER (n_st=74,n_end=n-16,n_tot=n_end-n_st+1)     !! n_tot=56, n-end=129
PARAMETER (half=0.5D0,one=1.D0,two=2.D0,zero=0.D0,quart=0.25D0,three=3.D0,five=5.D0)
PARAMETER (thd=one/three,tthd=two/three,fth=one/five)

LOGICAL ave_check
! read variables
REAL*4 zz(m,n),mmzz(m,n)

! temporary variable
REAL*4 out(m,n),temp(n,m)
REAL*8 ptemp(n),temp1(m,n),temp2(m,n)
REAL*4 Ave2(m,n),Ave4(m,n),Ave6(m,n),Ave8(m,n),Ave10(m,n),Ave12(m,n),Ave14(m,n)
REAL*8 Ave1(m,n,10),Ave3(m,n,10),Ave5(m,n,10),Ave7(m,n,10),Ave9(m,n,10),Ave11(m,n,10),Ave13(m,n,10)

REAL*8 mzz(m,n,10),mte3(n,10),mte3f(n,10),dd(m,n) 
!
REAL*8 mqly(n,10),mqlyh(n,10)
!
REAL*8 test(m,n)
!
REAL*8 uu3(m,n,10),vv3(m,n,10),te3(m,n,10)
REAL*8 uu3x(m,n,10),vv3x(m,n,10),te3x(m,n,10)
REAL*8 psi(m,n,10),phi_zm(n,10),psi_zm(n,10),phi_r(10),u_zm(n,10)
!
REAL*8 mpql(m,n,10),mql(n,10),pql(m,n,10),&
       pqlza(n,10),chsi(m,n,10),pqlx(m,n,10)
!
REAL*8 eddy_act(m,n,10),hg_mql(n,10)
!
REAL*8 f1(m,n,10),f2(m,n,10),f3(m,n,10),div(m,n,10)
REAL*8 divx(m,n),divy(m,n),divp(m,n),divxyp(m,n,10)
REAL*8 divyp(m,n,10),HMdivyp(m,n,10),VMdivyp(m,n)
REAL*8 VMf1(m,n),VMf2(m,n),VMf3(m,n),VMeddy_act(m,n),VMdiv(m,n),VMdivxyp(m,n)
REAL*8 HMf1(m,n,10),HMf2(m,n,10),HMf3(m,n,10),HMeddy_act(m,n,10),HMdiv(m,n,10),HMdivxyp(m,n,10)
!
REAL*8 fc(n),one_cos_lat(n),Gama,cos_lat(n),tang(n),cos_hlat(n-1)
REAL*8 alanda(n),coef_chsi(n),phi_a(n),rr_cos_phi(n)
!
REAL*8 Teta_r(10),dTetar_dp(10)
!
REAL*8 delta_phi,delta_landa,two_delta_landa,two_delta_phi,phi
COMMON /Average_S/ Phi_a,delta_phi,delta_landa
!
REAL*8 delta_landa_sqi,delta_phi_sqi
REAL*8 sum,pressure,exner,mwk(n),quart_rr,ffm

!
!REAL*8 wkas(n,m)
REAL*8 phi_1,phi_2,phi_m,sum_area,weight
!
INTEGER:: level(10),tmax(4)
CHARACTER*11 infile1, infile2
CHARACTER*84 location2
CHARACTER*81 location
CHARACTER*95 locationi
!CHARACTER*9 location2
!CHARACTER*6 location
!CHARACTER*13 locationi
CHARACTER*3 lev(10)
CHARACTER*2 duration(4)
!CHARACTER*4 year(4)

character(len=len("Hd_XX.dat")):: out_H      !Hd_year, duration
character(len=len("Vd_XX.dat")):: out_V
character(len=len("EnstrophyXX.dat")):: enstrophy
character(len=len("FxypXX.dat")):: Fxyp
character(len=len("FxypXX.dat")):: Fyp
character(len=len("FxyXX.dat")):: Fxy
character(len=len("WA-XX.dat")):: WA

!
DATA lev/'100','200','300','400','500','600','700','800','900','000'/
DATA level/100,200,300,400,500,600,700,800,900,1000/
!DATA year/'2016','2017','2017','1617'/
DATA duration/'12','01','02','wi'/
!DATA tmax/124,124,112,360/
DATA tmax/2,3,4,5/
!
!***********************************************
Gama=Gas_const/spe_heat    !R/Cp
location='/home/pari/Documents/university/master/THESIS/MYWORK(summer2018)/3.flux-code/out-'
locationi='/home/pari/Documents/university/master/THESIS/MYWORK-summer2018/1.levels-gs/10lev-data/out-m-p/'
!location='./out-'
!locationi='./input-data/'
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
	
	c=4
	do while (c<5)       !! do for each duration 
!!-----------------------------------------------------------------
		delta_phi=(pi/180.D0)*1.25    
		delta_phi_sqi=one/(delta_phi**2)	
		delta_landa=delta_phi              
		delta_landa_sqi=one/(delta_landa**2)
		two_delta_landa=two*delta_landa
		two_delta_phi=two*delta_phi


        	j=1
        	phi=REAL(j-1)*delta_phi -half*pi
	    	phi_a(j)=phi
        	tang(j)=tan(phi)
        	cos_lat(j)=cos(phi)
        	one_cos_lat(j)=one/(rr*cos_lat(j))
        	fc(j)=two_omega*sin(phi)
	     	alanda(j)=delta_landa_sqi/((cos(phi))**2)
		coef_chsi(j)=cos_lat(j)*rr*delta_landa 
!
      		DO  j=2,n-1                        !! we work in our domain so we r here!!
        	phi=REAL(j-1)*delta_phi -half*pi
	     	phi_a(j)=phi
        	tang(j)=tan(phi)
		cos_lat(j)=cos(phi)
        	fc(j)=two_omega*sin(phi)
        	one_cos_lat(j)=one/(rr*cos_lat(j))
        	alanda(j)=delta_landa_sqi/(cos(phi)**2)
        	coef_chsi(j)=cos_lat(j)*rr*delta_landa
      		END DO

      		j=n
        	phi=REAL(j-1)*delta_phi -half*pi
	    	phi_a(j)=phi
        	tang(j)=undef
		cos_lat(j)=cos(phi)
        	one_cos_lat(j)=undef
        	fc(j)=two_omega*sin(phi) 
		alanda(j)=undef
		coef_chsi(j)=cos_lat(j)*rr*delta_landa
	              
        	quart_rr=quart*rr
!
      		DO  j=1,n-1                         !! why this loop?? for what??
        	phi=phi_a(j)+0.5*delta_phi
        	cos_hlat(j)=cos(phi)         !! why ????
      		ENDDO
!------------------------------------------------------------------
		2000 format (A,A,A)
		write(location2,2000) location, duration(c) ,'/'
		write(*,*) location2                                       !!-- duration(c), ".dat"
		1001 format(A,A,A)
		
         	OPEN(600,FILE="wave-activity.txt",STATUS='unknown'&
         	,action="readwrite",iostat=ioerr)  !txt output of enstrophy!
					if (ioerr /= 0) then
					write(*,*) "ioerrwa-txt=", ioerr
				   	!stop
					end if 

         	OPEN(660,FILE=location2//WA,FORM='unformatted'&
         	,STATUS='unknown',iostat=ioerr,ACCESS='direct',RECL=4*m*n_tot*8)
					if (ioerr /= 0) then
					write(*,*) "ioerrwaB=", ioerr
				   	!stop
					end if 


		DO  120 k=1,10 
			1000 format(A,A,A,A)
			write(infile1,1000) "mH", lev(k),  duration(c), ".dat"
			write(infile2,1000) "pH", lev(k), duration(c), ".dat"
			write(*,*) infile1, infile2
 
         	OPEN(10+k,FILE=locationi//infile2,FORM='unformatted'&
	       ,STATUS='old',iostat=ioerr,ACCESS='direct',RECL=4*m*n)  !! perturbation
					if (ioerr /= 0) then
					write(*,*) "ioerr-p=", ioerr
				   	!stop
					end if
!
         	OPEN(70+k,FILE=locationi//infile1,FORM='unformatted'&
	       ,STATUS='old',ACCESS='direct',iostat=ioerr,RECL=4*m*n)   !! mean file
					if (ioerr /= 0) then
					write(*,*) "ioerr-m=", ioerr
				   	!stop
					end if
!
         	READ(70+k,REC=1) mmzz              ! GH :Z
         	DO  j=1,n
           		DO  i=1,m
             			mzz(i,j,k)=gg*mmzz(i,j)
           		END DO
         	END DO
		write(*,*) "mGEO-potential(50,50)=", mzz(50,50,k), k  
		write(*,*) "mGEO-potential(60,60)=", mzz(60,60,k), k  
         	CLOSE(70+k)
!
		120  CONTINUE 
                 
! now mzz is the time mean geopotential
!******************************************************************
! computing the potential temperature from the hydrostatic relation
!Computing the reference state
! zonal averaging
     DO k=1,10
       DO  j=1,n
         sum=zero
         DO  i=1,m
    	  sum=sum+mzz(i,j,k)
         END DO
         phi_zm(j,k)=sum/dble(m)       !phi
	  !write(99,*) "phi_zm",j,",",k,"=", phi_zm(j,k)    !tik
       END DO
     END DO 
     
!
       sum_area=zero
       DO  j=n_st,n_end-1
         phi_1=phi_a(j)
         phi_2=phi_a(j+1)
         phi_m=half*(phi_1+phi_2)
         weight=cos(phi_m)
         sum_area=sum_area+weight
       END DO

! computing the latitudinal mean of phi_zm over the domain
     DO  k=1,10
       sum=zero
       DO  j=n_st,n_end-1
         phi_1=phi_a(j)
         phi_2=phi_a(j+1)
         phi_m=half*(phi_1+phi_2)
         weight=cos(phi_m)                 !! i cant understand the weight?? why?
         sum=sum+half*(phi_zm(j,k)+phi_zm(j+1,k))*weight
       ENDDO
       phi_r(k)=sum/sum_area
     END DO
!
! computing the mean Coriolis Parameter over the domain
         phi_1=phi_a(n_st)
         phi_2=phi_a(n_end)
         phi_m=half*(phi_1+phi_2)
         ffm=two_omega*sin(phi_m)
!
     DO k=1,10
       DO  j=n_st-2,n_end+2
         psi_zm(j,k)=phi_zm(j,k)/ffm
       END DO
     ENDDO  
 ! psi_zm is now the time and zonal mean geostrophic streamfunction

! computing the time and zonal mean zonal velocity 
     DO k=1,10
       DO  j=n_st-1,n_end+1
         u_zm(j,k)=-(1/rr)*(psi_zm(j+1,k)-psi_zm(j-1,k))/two_delta_phi
       END DO
       write(*,*) u_zm(40,k)
     ENDDO      

! computing the time and zonal mean potential temperature 
! Computing the reference potential temperature
!Teta_r(1)=(phi_r(2)-phi_r(1))/1.0D4
DO  k=1,9         
  ptemp(k)=(phi_r(k+1)-phi_r(k))/(1.0D4)
  pressure=DBLE(Level(k)+50)*1.D2         !!average of p over each 2 lev,k and k+1
  exner=(1.D5/pressure)**Gama
  ptemp(k)=-pressure*exner*ptemp(k)/Gas_const  
ENDDO
! computing d/dp of reference potential temperature
Teta_r(1)=ptemp(1)
DO  k=2,9                 
  dTetar_dp(k)=(ptemp(k)-ptemp(k-1))/(1.0D4)
  Teta_r(k)=0.5d0*(ptemp(k)+ptemp(k-1))
         write(*,*) 'teta-r',teta_r(k)
ENDDO
Teta_r(10)=ptemp(9)
!
!******************************************************************
! computing the horizontal gradient of mql, eddy activity,and the fluxes
! computing zonal mean potential temperature
   DO  j=n_st-1,n_end+1
     DO  k=1,9         
       ptemp(k)=(phi_zm(j,k+1)-phi_zm(j,k))/(1.0D4)
       pressure=DBLE(Level(k)+50)*1.D2
       exner=(1.D5/pressure)**Gama
       mTe3(j,k)=-pressure*exner*ptemp(k)/Gas_const        
     ENDDO
     mTe3f(j,1)=mTe3(j,1)       
     DO  k=2,9
	   mTe3f(j,k)=0.5d0*(mTe3(j,k)+mTe3(j,k-1))   
     ENDDO
     mTe3f(j,10)=mTe3(j,9)
   ENDDO
   
!
DO  k=2,9    
	            write(*,*) 'teta-m',mTe3f(100,k)       
! computing mean vertical vorticity
     	DO  j=n_st-1,n_end+1
       	mwk(j)=(psi_zm(j-1,k)-two*psi_zm(j,k)+psi_zm(j+1,k))*delta_phi_sqi     
       	mwk(j)=mwk(j)-tang(j)*(psi_zm(j+1,k)-psi_zm(j-1,k))/two_delta_phi            
       	mwk(j)=((1/rr)**2)*mwk(j)
	END DO 

! computing the mean stretching term, and mean PV   
!
     	DO  j=n_st-1,n_end+1
	     	vertical=(mTe3(j,k)-mTe3(j,k-1))*1.D-4
            	vertical2=(teta_r(k-1)-2.0*teta_r(k)+teta_r(k+1))*1.D-8
            	vertical2=vertical2/dtetar_dp(k)**2
            	vertical=vertical/dtetar_dp(k)-vertical2*mTe3f(j,k)
            	mql(j,k)=fc(j)+mwk(j)+ffm*vertical
      	END DO 
! here mql is the linearized mean potential vorticity    
!
     	DO  j=n_st,n_end+1
          	mqlyh(j,k)=(mql(j,k)-mql(j-1,k))/delta_phi        
     	ENDDO 
     	DO  j=n_st,n_end
	       mqly(j,k)=0.5d0*(mqlyh(j,k)+mqlyh(j+1,k)) 
     	END DO

END DO    ! end of k loop

!******************************************************************
	mpql(:,:,:)=zero    !!! ,,,eddy anstrophy,,,

	Ave1(:,:,:)=zero
	Ave3(:,:,:)=zero
    	Ave5(:,:,:)=zero
    	Ave7(:,:,:)=zero
	Ave9(:,:,:)=zero
    	Ave11(:,:,:)=zero
    	Ave13(:,:,:)=zero	          

	 Ave2(:,:)=zero
	 Ave4(:,:)=zero
	 Ave6(:,:)=zero
	 Ave8(:,:)=zero 
	 Ave10(:,:)=zero
	 Ave12(:,:)=zero	          
    	 Ave14(:,:)=zero
    
   Ave50=0.;  Ave51=0.
   Ave52=0.;  Ave53=0.
   Ave54=0.;  Ave55=0.
   Ave56=0.;  Ave57=0.
   Ave58=0.;  Ave59=0.
   Ave60=0.;  Ave61=0.
   Ave62=0.;  Ave63=0.
   Ave64=0.;  Ave65=0.
   Ave66=0.;  Ave67=0. 
   Ave68=0.;  Ave69=0. 
   Ave70=0.;  Ave71=0.
   Ave72=0.;  Ave73=0.
   Ave74=0.;  Ave75=0.
   Ave76=0.;  Ave77=0.
   
! loop over time

	DO  tt=1, tmax(c)           !! end tt loop at 980 line
!      
  	WRITE(*,*) "c=",c,"tt=",tt
 	irec=tt
 
 	f_input1=0.
  	f_input2=0.
 	f_input3=0.
  	f_input4=0.
  	f_input_total=0.
 	f_input1w=0.
 	f_input2w=0.
 	f_input3w=0.
  	f_input4w=0.
 	f_input_totalw=0.
 	f_input1c=0.
 	f_input2c=0.
 	f_input3c=0.
 	f_input4c=0.
 	f_input_totalc=0.
 	f_input1e=0.
 	f_input2e=0.
 	f_input3e=0.
 	f_input4e=0.
  	f_input_totale=0.
  	eamean=0.0
  	eamean_w=0.0
  	eamean_c=0.0
  	eamean_e=0.0
  	fdmean=0.0
  	fdmean_w=0.0
  	fdmean_c=0.0
  	fdmean_e=0.0
            
   		DO k=1,10

      			READ(10+k,REC=irec) zz     !! from perturbation file
     			DO  j=n_st-2,n_end+2
       		DO  i=1,m
            		psi(i,j,k)=gg*DBLE(zz(i,j))                     ! here psi is the geopotential

       		END DO
     			END DO
		       write(*,*) "Z'(k) =",k, zz(50,50)  
		       write(*,*) "Z'(k) =",k, zz(60,60)  
   		END DO

! computing perturbation potential temperature (teta'=Te3)
  		DO  k=1,9
       		pressure=DBLE(Level(k)+50)*1.D2
       		exner=(1.D5/pressure)**Gama
       		const_exner=-pressure*exner/Gas_const
       		DO  j=n_st,n_end
         		DO  i=1,m
           		Te3(i,j,k)=(psi(i,j,k+1)-psi(i,j,k))/1.0D4
		   		Te3(i,j,k)=const_exner*Te3(i,j,k)
         		ENDDO
       		ENDDO
			write(*,*) 'teta-p,40,100,k',mTe3f(40,100,k) 
  		ENDDO
! computing the geostrophic streamfunction'
   		DO k=1,10
    			DO  j=n_st-2,n_end+2
       			DO  i=1,m
         			psi(i,j,k)=psi(i,j,k)/ffm
       			ENDDO
      			ENDDO
   		END DO
!

 		DO  k=2,9        !!!!end at line  ...
!  computing perturbation zonal and meridional velocities
         		DO  j=n_st,n_end
		   	DO  i=1,m 
		     		uu3(i,j,k)=-(1/rr)*(psi(i,j+1,k)-psi(i,j-1,k))/two_delta_phi 
			 	test(i,j)=half*(Te3(i,j,k-1)+Te3(i,j,k))       
		   	END DO
		   		vv3(1,j,k)=one_cos_lat(j)*(psi(2,j,k)-psi(m,j,k))/two_delta_landa
		   		te3x(1,j,k)=one_cos_lat(j)*(test(2,j)-test(m,j))/two_delta_landa
		   	DO  i=2,m-1 
		     		vv3(i,j,k)=one_cos_lat(j)*(psi(i+1,j,k)-psi(i-1,j,k))/two_delta_landa 
           			te3x(i,j,k)=one_cos_lat(j)*(test(i+1,j)-test(i-1,j))/two_delta_landa 
		  	END DO
           			vv3(m,j,k)=one_cos_lat(j)*(psi(1,j,k)-psi(m-1,j,k))/two_delta_landa           
           			te3x(m,j,k)=one_cos_lat(j)*(test(1,j)-test(m-1,j))/two_delta_landa  
			 END DO
     
! computing the zonal derivative of u and v
         		DO  j=n_st,n_end
		   		uu3x(1,j,k)=one_cos_lat(j)*(uu3(2,j,k)-uu3(m,j,k))/two_delta_landa
		   		vv3x(1,j,k)=one_cos_lat(j)*(vv3(2,j,k)-vv3(m,j,k))/two_delta_landa
		   	DO  i=2,m-1 
		     		uu3x(i,j,k)=one_cos_lat(j)*(uu3(i+1,j,k)-uu3(i-1,j,k))/two_delta_landa 
           			vv3x(i,j,k)=one_cos_lat(j)*(vv3(i+1,j,k)-vv3(i-1,j,k))/two_delta_landa 
		   	END DO
           			uu3x(m,j,k)=one_cos_lat(j)*(uu3(1,j,k)-uu3(m-1,j,k))/two_delta_landa           
           			vv3x(m,j,k)=one_cos_lat(j)*(vv3(1,j,k)-vv3(m-1,j,k))/two_delta_landa           
		 	END DO

!   
!! q' calculation (3 terms!!)
! computing perturbation vertical vorticity      
         		DO  j=n_st,n_end
           			dd(1,j)=(psi(m,j,k)-two*psi(1,j,k)+psi(2,j,k))*alanda(j)
           			DO  i=2,m-1
             			dd(i,j)=(psi(i-1,j,k)-two*psi(i,j,k)+psi(i+1,j,k))*alanda(j)
           			END DO
           			dd(m,j)=(psi(m-1,j,k)-two*psi(m,j,k)+psi(1,j,k))*alanda(j)
         		END DO 
!
         		DO  j=n_st,n_end
           			DO  i=1,m
             			dd(i,j)=dd(i,j)+(psi(i,j-1,k)-two*psi(i,j,k)+psi(i,j+1,k))*delta_phi_sqi    
             			dd(i,j)=dd(i,j)+rr*tang(j)*uu3(i,j,k)               
             			dd(i,j)=((1/rr)**2)*dd(i,j)
           			END DO
         		END DO

! computing the stretching term
!		 
         		DO  j=n_st,n_end
           			DO  i=1,m
		          	vertical=(Te3(i,j,k)-Te3(i,j,k-1))*1.D-4
                		vertical2=(teta_r(k-1)-two*teta_r(k)+teta_r(k+1))*1.D-8
 !              	 	vertical2=0.d0
		          	vertical2=vertical2/dtetar_dp(k)**2
		          	vertical=vertical/dtetar_dp(k)-vertical2*test(i,j)
                		dd(i,j)=dd(i,j)+ffm*vertical      !! now dd is q' 
           			END DO
         		END DO

!HORIZONTAL AVEARAGING
         		DO  j=n_st,n_end   
         	  		DO  i=1,m
             			temp1(i,j)=dd(i,j)
		   		ENDDO
		 	ENDDO   

         		CALL Averaging2D(temp1,n_st,n_end,temp2)		   !! temp2 is output of averaging (pql),temp1 is input!   
		   
         		DO  j=n_st,n_end	   !! end at 670	   
		   		pqlza(j,k)=zero         
		   		DO  i=1,m	 
             			pql(i,j,k)=temp2(i,j)
!	         		out(i,j,k)=pql(i,j,k)
             			mpql(i,j,k)=mpql(i,j,k)+pql(i,j,k)**2        ! here mpql is proportional to eddy potential enstrophy
           			END DO
		   			pqlza(j,k)=pqlza(j,k)/dble(m)               !! zonal mean of pql
		   			chsi(1,j,k)=zero                            
		   			chsi_ave=zero
!
           			DO  i=2,m
			 		chsi(i,j,k)=chsi(i-1,j,k)+half*(pql(i,j,k)+pql(i-1,j,k))-pqlza(j,k)        !!??chsi= sum(over i) of derivation of pql(mean over i,i-1)and zonal mean-pql
			 		chsi_ave=chsi_ave+chsi(i,j,k)
           			END DO
           			chsi_ave=chsi_ave/dble(m)
           			DO  i=1,m
              		chsi(i,j,k)=chsi(i,j,k)-chsi_ave      !! derivation of chsi from it zonal mean (chsi=chsi')
           			END DO
!
           			DO  i=2,m
			 		chsi(i,j,k)=chsi(i-1,j,k)+half*(pql(i,j,k)+pql(i-1,j,k))-pqlza(j,k) 
           			END DO
!*************************************************************
! CORRECTED
!           DO  i=2,m
           			DO  i=1,m
			 	     chsi(i,j,k)=coef_chsi(j)*chsi(i,j,k)               !! why? to antegrate??
           			END DO
		   		pqlx(1,j,k)=(pql(2,j,k)-pql(m,j,k))/two_delta_landa    
		   		DO  i=2,m-1
		     		pqlx(i,j,k)=(pql(i+1,j,k)-pql(i-1,j,k))/two_delta_landa
		   		ENDDO 
		   		pqlx(m,j,k)=(pql(1,j,k)-pql(m-1,j,k))/two_delta_landa
!
		   		DO  i=1,m
		     		eddy_act(i,j,k)=quart*(rr*cos_lat(j)*pql(i,j,k)**2-chsi(i,j,k)*pqlx(i,j,k)) 
		     		eddy_act(i,j,k)=eddy_act(i,j,k)/mqly(j,k)
				   !write(*,*) "ea:", eddy_act(i,j,k)
              	write(600,*) "eddy activity=",i,",",j,",",k,"=", eddy_act(i,j,k)   
             			!if (c==4 .and. mod(tt-1,4)==0)	write(660,rec=tt) eddy_act(i,j,k)   
		   		ENDDO
!
           			coef= half*cos_lat(j)     
		   		DO  i=1,m
			 		f1(i,j,k)=half*(psi(i,j,k)*pql(i,j,k)-vv3(i,j,k)*chsi(i,j,k)) + &
                       		vv3(i,j,k)**2-psi(i,j,k)*vv3x(i,j,k)                            
             				f2(i,j,k)=-uu3(i,j,k)*vv3(i,j,k)+psi(i,j,k)*uu3x(i,j,k)                          
			 		f3(i,j,k)=(fc(j)/dTetar_dp(k))*(vv3(i,j,k)*test(i,j)-   &
      			      		psi(i,j,k)*te3x(i,j,k))
     
             				f1(i,j,k)=u_zm(j,k)*eddy_act(i,j,k)+coef*f1(i,j,k)   
!             			f1(i,j,k)=(coef*f1(i,j,k))             
			 		f2(i,j,k)=(coef*f2(i,j,k))      
			 		f3(i,j,k)=(coef*f3(i,j,k))     
		   		ENDDO 
         		END DO   !! do j at line 608

!************ computing scaled 3D flux divergence ***
    			DO  j=n_st,n_end
     				divx(1,j)=(f1(2,j,k)-f1(m,j,k))/two_delta_landa
     				DO  i=2,m-1     
       			divx(i,j)=(f1(i+1,j,k)-f1(i-1,j,k))/two_delta_landa
     				ENDDO
     				divx(m,j)=(f1(1,j,k)-f1(m-1,j,k))/two_delta_landa
   			ENDDO
   			j=n_st
   			DO  i=1,m 
     			divy(i,j)=(f2(i,j+1,k)*cos_lat(j+1)-f2(i,j,k)*cos_lat(j))/delta_phi        !we add coslat(j) because for y component we add 1/coslat later...) 
   			ENDDO
   			DO  j=n_st+1,n_end-1
     				DO  i=1,m
       			divy(i,j)=(f2(i,j+1,k)*cos_lat(j+1)-f2(i,j-1,k)*cos_lat(j-1))/two_delta_phi
     				ENDDO
   			ENDDO
   			j=n_end
   			DO  i=1,m 
     			divy(i,j)=(f2(i,j,k)*cos_lat(j)-f2(i,j-1,k)*cos_lat(j-1))/delta_phi
   			ENDDO
   			DO  j=n_st,n_end
     				DO  i=1,m
       			if (k .eq. 2) then                                              !!so what about k==1, 10??? (k loop start at 2, so k=2 is a boundary)
         				divp(i,j)=(f3(i,j,k+1)-f3(i,j,k))*1.D-4              !! delta p=10^4
       			elseif (k .eq. 9) then
         				divp(i,j)=(f3(i,j,k)-f3(i,j,k-1))*1.D-4  
       			else                                                        !! 2<k<9
         				divp(i,j)=(f3(i,j,k+1)-f3(i,j,k-1))*5.D-5                  ! 1/2 * 10^4
       			endif 
     				ENDDO
   			ENDDO
! 
   			DO  j=n_st,n_end
     				DO  i=1,m
       				div(i,j,k)=one_cos_lat(j)*(divx(i,j)+divy(i,j))
       				divyp(i,j,k)=one_cos_lat(j)*divy(i,j)+divp(i,j)
       				divxyp(i,j,k)=div(i,j,k)+divp(i,j)
     				ENDDO
   			ENDDO   

!*************************************************************************
 		ENDDO                ! END OF VERTICAL LOOP,536
    
!*******************COMPUTING Horizontal and Vertical Mean
   	ave_check=.TRUE.      !This check is applied for using horizontal average of parameters in vertical averaging
   	CALL Averaging(ave_check,eddy_act,n_st,n_end,HMeddy_act,VMeddy_act)
   	CALL Averaging(ave_check,f1,n_st,n_end,HMf1,VMf1)
!   	ave_check=.FALSE.
   	CALL Averaging(ave_check,f2,n_st,n_end,HMf2,VMf2)
!   	ave_check=.TRUE.
   	CALL Averaging(ave_check,f3,n_st,n_end,HMf3,VMf3)
   	CALL Averaging(ave_check,div,n_st,n_end,HMdiv,VMdiv)
   	CALL Averaging(ave_check,divxyp,n_st,n_end,HMdivxyp,VMdivxyp)
   	CALL Averaging(ave_check,divyp,n_st,n_end,HMdivyp,VMdivyp)  !!?? I ADDED!!NOTICE t??************888

!**********calculate 3D input flux to the Mediterranean area from the Atlantic
    !IF ((k.GE.2).and.(k.LE.6)) THEN  
    	i1=280       !left bound of rectangular   (20-40n,10w-60e )
    	i2=11        !right bound of west rectangular and left bound of central rectangular
    	i3=30        !right bound of central rectangular and left bound of east rectangular
    	i4=48       !right bound of east rectangular
    
    	J_lower=89   !lower bound of rectangular
    	J_upper= 105     !upper bound of rectangular

       sa_w=0.0           !!??? for what???
       sa_c=0.0
       sa_e=0.0

       	DO  j=J_lower,J_upper-1
          	sa_w=sa_w+cos_hlat(j)
       	ENDDO 
       sa_w=sa_w*16.0
       sa_c=sa_w
       sa_e=sa_w
       sa_t=sa_w+sa_c+sa_e
    
    		Do k=2,6    
	  	!----------Over Total Rectangular
	  	j=J_lower  !f_input1
	  		DO i=i1,287
	    		f_input1(k)=f_input1(k)-half*cos_lat(j)*(f2(i,j,k)+f2(i+1,j,k))               !! why ----????
	  		ENDDO
	     		f_input1(k)=f_input1(k)-half*cos_lat(j)*(f2(288,j,k)+f2(1,j,k))          !! to calculate the F for i=288!!we nees coslatj to antegrate
	  		DO i=1,i4
	    		f_input1(k)=f_input1(k)-half*cos_lat(j)*(f2(i,j,k)+f2(i+1,j,k))
	  		ENDDO
	    
	  	i=i4  !f_input2
	  		DO j=J_lower,J_upper
				f_input2(k)=f_input2(k)+half*(f1(i,j,k)+f1(i,j+1,k))                             !!---
      			ENDDO

	  	j=J_upper  !f_input3
	  		DO i=i1,287
	    		f_input3(k)=f_input3(k)+half*cos_lat(j)*(f2(i,j,k)+f2(i+1,j,k))               !! it should be -!!!!
	  		ENDDO
	  		f_input3(k)=f_input3(k)+half*cos_lat(j)*(f2(288,j,k)+f2(1,j,k))      
	  		DO i=1,i4
	    		f_input3(k)=f_input3(k)+half*cos_lat(j)*(f2(i,j,k)+f2(i+1,j,k))
	  		ENDDO

	  	i=i1  !f_input4
	  		DO j=J_lower,J_upper
			f_input4(k)=f_input4(k)-half*(f1(i,j,k)+f1(i,j+1,k))          
      			ENDDO

	  	!----------Over West Rectangular
	  	j=J_lower  !f_inputw1    lower boundary 
	  		DO i=i1,287
	    		f_input1w(k)=f_input1w(k)-half*cos_lat(j)*(f2(i,j,k)+f2(i+1,j,k))
	  		ENDDO
	  		f_input1w(k)=f_input1w(k)-half*cos_lat(j)*(f2(288,j,k)+f2(1,j,k))
	  		DO i=1,i2-1
	    		f_input1w(k)=f_input1w(k)-half*cos_lat(j)*(f2(i,j,k)+f2(i+1,j,k))
	  		ENDDO
	    
	  	i=i2  !f_inputw2          right boundary 
	  		DO j=J_lower,J_upper
			f_input2w(k)=f_input2w(k)+half*(f1(i,j,k)+f1(i,j+1,k))   !!should be --
      			ENDDO

	  	j=J_upper  !f_inputw3        upper boundary 
	  		DO i=i1,287
	    		f_input3w(k)=f_input3w(k)+half*cos_lat(j)*(f2(i,j,k)+f2(i+1,j,k))
	  		ENDDO
	  		f_input3w(k)=f_input3w(k)+half*cos_lat(j)*(f2(288,j,k)+f2(1,j,k))
	  		DO i=1,i2-1
	    		f_input3w(k)=f_input3w(k)+half*cos_lat(j)*(f2(i,j,k)+f2(i+1,j,k))
	  		ENDDO

	  	i=i1  !f_inputw4          left boundary
	  		DO j=J_lower,J_upper
			f_input4w(k)=f_input4w(k)-half*(f1(i,j,k)+f1(i,j+1,k))
      			ENDDO
! compute area mean of eddy activity and 3D flux divergence
          		DO  j=J_lower,J_upper-1
            		DO  i=i1,m
              		eamean_w(k)=eamean_w(k)+cos_hlat(j)*eddy_act(i,j,k)           !!! why cos_hlat????
              		fdmean_w(k)=fdmean_w(k)+cos_hlat(j)*divxyp(i,j,k)                
            		ENDDO
            		DO  i=1,i2-1
              		eamean_w(k)=eamean_w(k)+cos_hlat(j)*eddy_act(i,j,k)
              		fdmean_w(k)=fdmean_w(k)+cos_hlat(j)*divxyp(i,j,k)
            		ENDDO
          		ENDDO

	  	!----------Over Central Rectangular
	  	j=J_lower       !lower boundary 
	  		DO i=i2,i3-1
	    		f_input1c(k)=f_input1c(k)-half*cos_lat(j)*(f2(i,j,k)+f2(i+1,j,k))
	  		ENDDO
	    
	  	i=i3  !   right boundary 
	  		DO j=J_lower,J_upper
	    		f_input2c(k)=f_input2c(k)+half*(f1(i,j,k)+f1(i,j+1,k))
      			ENDDO

	  	j=J_upper  !   upper boundary 
	  		DO i=i2,i3-1
	    		f_input3c(k)=f_input3c(k)+half*cos_lat(j)*(f2(i,j,k)+f2(i+1,j,k))
	  		ENDDO

	  	i=i2  !  left boundary
	  		DO j=J_lower,J_upper
			f_input4c(k)=f_input4c(k)-half*(f1(i,j,k)+f1(i,j+1,k))
      			ENDDO
! compute area mean of eddy activity and 3D flux divergence
          		DO  j=J_lower,J_upper-1
            		DO  i=i2,i3-1
              		eamean_c(k)=eamean_c(k)+cos_hlat(j)*eddy_act(i,j,k)
              		fdmean_c(k)=fdmean_c(k)+cos_hlat(j)*divxyp(i,j,k)
            		ENDDO
          		ENDDO

			   !------------over east rectangular
	  	j=J_lower  !f_inputw1      lower boundary 
	  		DO i=i3,i4
	    		f_input1e(k)=f_input1e(k)-half*cos_lat(j)*(f2(i,j,k)+f2(i+1,j,k))
	  		ENDDO
	    
	  	i=i4  !f_inputw2         right boundary 
	  		DO j=J_lower,J_upper
	    		f_input2e(k)=f_input2e(k)+half*(f1(i,j,k)+f1(i,j+1,k))
      			ENDDO

	  	j=J_upper  !f_inputw3      upper boundary 
	  		DO i=i3,i4
	    		f_input3e(k)=f_input3e(k)+half*cos_lat(j)*(f2(i,j,k)+f2(i+1,j,k))
	  		ENDDO

	  	i=i3  !f_inputw4       left boundary
	  		DO j=J_lower,J_upper
			f_input4e(k)=f_input4e(k)-half*(f1(i,j,k)+f1(i,j+1,k))
      			ENDDO
! compute area mean of eddy activity and 3D flux divergence
          		DO  j=J_lower,J_upper-1
            		DO  i=i3,i4
              		eamean_e(k)=eamean_e(k)+cos_hlat(j)*eddy_act(i,j,k)
              		fdmean_e(k)=fdmean_e(k)+cos_hlat(j)*divxyp(i,j,k)
            		ENDDO
          		ENDDO

          	eamean(k)=(eamean_w(k)+eamean_c(k)+eamean_e(k))/sa_t
          	eamean_w(k)=eamean_w(k)/sa_w
          	eamean_c(k)=eamean_c(k)/sa_c
          	eamean_e(k)=eamean_e(k)/sa_e

          	fdmean(k)=(fdmean_w(k)+fdmean_c(k)+fdmean_e(k))/sa_t
          	fdmean_w(k)=fdmean_w(k)/sa_w
          	fdmean_c(k)=fdmean_c(k)/sa_c
          	fdmean_e(k)=fdmean_e(k)/sa_e			   
   ! 
   		f_input1(k)=rr*delta_landa*f_input1(k)         !! we need a coslatj(for dlanda) that we add when we calculate the sum!:)
   		f_input2(k)=rr*delta_phi*f_input2(k)
   		f_input3(k)=rr*delta_landa*f_input3(k)
   		f_input4(k)=rr*delta_phi*f_input4(k)
!
   		f_input1w(k)=rr*delta_landa*f_input1w(k)        !! to antegrate
   		f_input2w(k)=rr*delta_phi*f_input2w(k)
   		f_input3w(k)=rr*delta_landa*f_input3w(k)
   		f_input4w(k)=rr*delta_phi*f_input4w(k)
!
   		f_input1c(k)=rr*delta_landa*f_input1c(k)
   		f_input2c(k)=rr*delta_phi*f_input2c(k)
   		f_input3c(k)=rr*delta_landa*f_input3c(k)
   		f_input4c(k)=rr*delta_phi*f_input4c(k)
!
   		f_input1e(k)=rr*delta_landa*f_input1e(k)
   		f_input2e(k)=rr*delta_phi*f_input2e(k)
   		f_input3e(k)=rr*delta_landa*f_input3e(k)
   		f_input4e(k)=rr*delta_phi*f_input4e(k)

   	END DO !for k
    	write(*,*) "I compelete ave of each part(E,C,W,total) of flux(k)"  

   	f_input1_3D=0.;   f_input2_3D=0.
   	f_input3_3D=0.;   f_input4_3D=0.
   
   	f_input1w_3D=0.;  f_input2w_3D=0.
   	f_input3w_3D=0.;  f_input4w_3D=0.
   !
   	f_input1c_3D=0.;  f_input2c_3D=0.
   	f_input3c_3D=0.;  f_input4c_3D=0.
   !
   	f_input1e_3D=0.;  f_input2e_3D=0.
   	f_input3e_3D=0.;  f_input4e_3D=0.
   !
   	eamean_3D=0.;     eamean_w_3D=0.
   	eamean_c_3D=0.;   eamean_e_3D=0.
   !
   	fdmean_3D=0.;     fdmean_w_3D=0.
   	fdmean_c_3D=0.;   fdmean_e_3D=0.
  
	Do k=2,6 

!   	DO k=2,4          
     	f_input1_3D=f_input1_3D+(f_input1(k)*fth)     
     	f_input2_3D=f_input2_3D+(f_input2(k)*fth)
     	f_input3_3D=f_input3_3D+(f_input3(k)*fth)
     	f_input4_3D=f_input4_3D+(f_input4(k)*fth)
 !
     	f_input1w_3D=f_input1w_3D+(f_input1w(k)*fth)
     	f_input2w_3D=f_input2w_3D+(f_input2w(k)*fth)
     	f_input3w_3D=f_input3w_3D+(f_input3w(k)*fth)
     	f_input4w_3D=f_input4w_3D+(f_input4w(k)*fth)
!
     	f_input1c_3D=f_input1c_3D+(f_input1c(k)*fth)
     	f_input2c_3D=f_input2c_3D+(f_input2c(k)*fth)
     	f_input3c_3D=f_input3c_3D+(f_input3c(k)*fth)
     	f_input4c_3D=f_input4c_3D+(f_input4c(k)*fth)
!
     	f_input1e_3D=f_input1e_3D+(f_input1e(k)*fth)
     	f_input2e_3D=f_input2e_3D+(f_input2e(k)*fth)
     	f_input3e_3D=f_input3e_3D+(f_input3e(k)*fth)
     	f_input4e_3D=f_input4e_3D+(f_input4e(k)*fth)
!
     	eamean_3D=eamean_3D+(eamean(k)*fth) 
      write(*,*) "ea3d",eamean_3D   
     	eamean_w_3D=eamean_w_3D+(eamean_w(k)*fth)  
      write(*,*) "ea3dw",eamean_3D     
     	eamean_c_3D=eamean_c_3D+(eamean_c(k)*fth) 
      write(*,*) "ea3dc",eamean_3D     
     	eamean_e_3D=eamean_e_3D+(eamean_e(k)*fth)
      write(*,*) "ea3de",eamean_3D     
!
     	fdmean_3D=fdmean_3D+(fdmean(k)*fth)    
     	fdmean_w_3D=fdmean_w_3D+(fdmean_w(k)*fth)    
     	fdmean_c_3D=fdmean_c_3D+(fdmean_c(k)*fth) 
     	fdmean_e_3D=fdmean_e_3D+(fdmean_e(k)*fth)
   	ENDDO

   	f_input_total= f_input1_3D  + f_input2_3D  + f_input3_3D  + f_input4_3D
   	f_input_totalw=f_input1w_3D + f_input2w_3D + f_input3w_3D + f_input4w_3D
   	f_input_totalc=f_input1c_3D + f_input2c_3D + f_input3c_3D + f_input4c_3D
   	f_input_totale=f_input1e_3D + f_input2e_3D + f_input3e_3D + f_input4e_3D
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

     	Ave50=Ave50+f_input_total
	!write (*,*) "sfinputt=" , ave50
     	Ave51=Ave51+f_input1_3D
     	Ave52=Ave52+f_input2_3D
     	Ave53=Ave53+f_input3_3D
     	Ave54=Ave54+f_input4_3D
    !
     	Ave55=Ave55+f_input_totalw
	!write (*,*) "sfinputtw=" , ave55
     	Ave56=Ave56+f_input1w_3D
     	Ave57=Ave57+f_input2w_3D
     	Ave58=Ave58+f_input3w_3D
     	Ave59=Ave59+f_input4w_3D
     !
     	Ave60=Ave60+f_input_totalc
	!write (*,*) "sfinputtc=" , ave60
     	Ave61=Ave61+f_input1c_3D
     	Ave62=Ave62+f_input2c_3D
     	Ave63=Ave63+f_input3c_3D
     	Ave64=Ave64+f_input4c_3D
     !
     	Ave65=Ave65+f_input_totale
	!write (*,*) "sfinputte=" , ave65            
     	Ave66=Ave66+f_input1e_3D
     	Ave67=Ave67+f_input2e_3D
     	Ave68=Ave68+f_input3e_3D
     	Ave69=Ave69+f_input4e_3D
     !
     	Ave70=Ave70+eamean_3D
	!write (*,*) "sea3d=" , ave70            
     	Ave71=Ave71+eamean_w_3D
     	Ave72=Ave72+eamean_c_3D
     	Ave73=Ave73+eamean_e_3D
     !
     	Ave74=Ave74+fdmean_3D
	!write (*,*) "sfd3d=" , ave74            
     	Ave75=Ave75+fdmean_w_3D
     	Ave76=Ave76+fdmean_c_3D
     	Ave77=Ave77+fdmean_e_3D

 !**********************************************************************************************
 
  	DO  k=2,9
   		DO  j=n_st,n_end
     		DO  i=1,m
       	Ave1(i,j,k)=Ave1(i,j,k)+HMeddy_act(i,j,k) 
       		write (*,*) "ave1,i,j,k= " , ave1(i,j,k)   
       	Ave3(i,j,k)=Ave3(i,j,k)+HMf1(i,j,k) 
       	Ave5(i,j,k)=Ave5(i,j,k)+HMf2(i,j,k) 
       	Ave7(i,j,k)=Ave7(i,j,k)+HMf3(i,j,k) 
       	Ave9(i,j,k)=Ave9(i,j,k)+HMdiv(i,j,k) 
       	Ave11(i,j,k)=Ave11(i,j,k)+HMdivyp(i,j,k) 
       	Ave13(i,j,k)=Ave13(i,j,k)+HMdivxyp(i,j,k) 
      		END DO
   		END DO
 	END DO
  
  		DO  j=n_st,n_end
     		DO  i=1,m
       	Ave2(i,j)=Ave2(i,j)+VMeddy_act(i,j) 
       	Ave4(i,j)=Ave4(i,j)+VMf1(i,j) 
       	Ave6(i,j)=Ave6(i,j)+VMf2(i,j) 
       	Ave8(i,j)=Ave8(i,j)+VMf3(i,j) 
       	Ave10(i,j)=Ave10(i,j)+VMdiv(i,j) 
       	Ave12(i,j)=Ave12(i,j)+VMdivyp(i,j) 
      	 	Ave14(i,j)=Ave14(i,j)+VMdivxyp(i,j) 
      		END DO
   		END DO
    
END DO                !END OF TIME LOOP,tt
 
	write (*,*) "tmax ,c =" , tmax(c) ,c 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     	f_input_total=Ave50/(REAL(tmax(c))) 
     	f_input1_3D=Ave51/(REAL(tmax(c)))  
     	f_input2_3D=Ave52/(REAL(tmax(c)))
     	f_input3_3D=Ave53/(REAL(tmax(c)))
     	f_input4_3D=Ave54/(REAL(tmax(c)))        
     !
     	f_input_totalw=Ave55/(REAL(tmax(c)))
	!write (*,*) "mfinputtw=" , ave55        
     	f_input1w_3D=Ave56/(REAL(tmax(c)))     !! this is time mean for each duration
     	f_input2w_3D=Ave57/(REAL(tmax(c)))
     	f_input3w_3D=Ave58/(REAL(tmax(c)))
     	f_input4w_3D=Ave59/(REAL(tmax(c)))   
     !
     	f_input_totalc=Ave60/(REAL(tmax(c)))
	!write (*,*) "mfinputtc=" , ave60       
     	f_input1c_3D=Ave61/(REAL(tmax(c)))     
     	f_input2c_3D=Ave62/(REAL(tmax(c)))
     	f_input3c_3D=Ave63/(REAL(tmax(c)))
     	f_input4c_3D=Ave64/(REAL(tmax(c)))   
     !
     	f_input_totale=Ave65/(REAL(tmax(c)))
	!write (*,*) "mfinputte=" , ave65            
     	f_input1e_3D=Ave66/(REAL(tmax(c)))    
     	f_input2e_3D=Ave67/(REAL(tmax(c)))
     	f_input3e_3D=Ave68/(REAL(tmax(c)))
     	f_input4e_3D=Ave69/(REAL(tmax(c)))   
     !
     	eamean_3D=Ave70/(REAL(tmax(c)))    
     	eamean_w_3D=Ave71/(REAL(tmax(c)))     
     	eamean_c_3D=Ave72/(REAL(tmax(c)))
     	eamean_e_3D=Ave73/(REAL(tmax(c)))
     !
     	fdmean_3D=Ave74/(REAL(tmax(c)))   
     	fdmean_w_3D=Ave75/(REAL(tmax(c)))     
     	fdmean_c_3D=Ave76/(REAL(tmax(c)))    
     	fdmean_e_3D=Ave77/(REAL(tmax(c)))
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	write (*,*) "sHMEddyact= " , ave1(41,30,2)
       HMeddy_act(:,n_st:n_end,2:9)= Ave1(:,n_st:n_end,2:9)/(REAL(tmax(c)))     !! this is time mean for each duration
	write (*,*) "meanHMEddyact= " , HMeddy_act(41,30,2)
	HMf1(:,n_st:n_end,2:9) = Ave3(:,n_st:n_end,2:9)/(REAL(tmax(c)))
	HMf2(:,n_st:n_end,2:9) = Ave5(:,n_st:n_end,2:9)/(REAL(tmax(c)))
	HMf3(:,n_st:n_end,2:9) = Ave7(:,n_st:n_end,2:9)/(REAL(tmax(c)))
	HMdiv(:,n_st:n_end,2:9) = Ave9(:,n_st:n_end,2:9)/(REAL(tmax(c)))
      	HMdivyp(:,n_st:n_end,2:9)=Ave11(:,n_st:n_end,2:9)/(REAL(tmax(c)))
      	HMdivxyp(:,n_st:n_end,2:9)=Ave13(:,n_st:n_end,2:9)/(REAL(tmax(c)))
 
   		DO  j=n_st,n_end
     		DO  i=1,m
       	VMeddy_act(i,j)= Ave2(i,j)/(REAL(tmax(c))) 

       	VMf1(i,j)= Ave4(i,j)/(REAL(tmax(c))) 
       	VMf2(i,j)= Ave6(i,j)/(REAL(tmax(c))) 
       	VMf3(i,j)= Ave8(i,j)/(REAL(tmax(c)))
       	VMdiv(i,j)= Ave10(i,j)/(REAL(tmax(c)))  
       	VMdivyp(i,j)= Ave12(i,j)/(REAL(tmax(c))) 
       	VMdivxyp(i,j)= Ave14(i,j)/(REAL(tmax(c))) 
      		END DO
   		END DO
	write (*,*) "meanVMEddyact= " , VMeddy_act(41,30)
    	write(*,*) "I compelete flux3D-timeave"  
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	800 format (A8,6x,A12,1X,A12,1X,A12,1X,A12,2X,A14)
	!800 format (A4,2X,3(A13),2X,A13,3X,A13)
	open(100,file=location2//'analysis-BigRec3D.txt',iostat=ioerr)
					if (ioerr /= 0) then
					write(*,*) "ioerr100=", ioerr
				   	!stop
					end if
!
	open(101,file=location2//'analysis-WestRec3D.txt')
	open(102,file=location2//'analysis-CentralRec3D.txt')
	open(103,file=location2//'analysis-EastRec3D.txt')

	open(110,file=location2//'mean-ea_all.txt')  
	open(111,file=location2//'mean-fd_all.txt') 

	OPEN(104,FILE=location2//'Analysis3D_total.dat',FORM='unformatted',STATUS='unknown',iostat=ioerr,ACCESS='direct',RECL=4)
					if (ioerr /= 0) then
					write(*,*) "ioerr104=", ioerr
				   	!stop
					end if
	OPEN(105,FILE=location2//'Analysis3D_WRect.dat',FORM='unformatted',STATUS='unknown',ACCESS='direct',RECL=4)
	OPEN(106,FILE=location2//'Analysis3D_CRect.dat',FORM='unformatted',STATUS='unknown',ACCESS='direct',RECL=4)
	OPEN(107,FILE=location2//'Analysis3D_ERect.dat',FORM='unformatted',STATUS='unknown',ACCESS='direct',RECL=4)

	OPEN(114,FILE=location2//'mean3D_total.dat',FORM='unformatted',STATUS='unknown',ACCESS='direct',RECL=4)
	OPEN(115,FILE=location2//'mean3D_WRect.dat',FORM='unformatted',STATUS='unknown',ACCESS='direct',RECL=4)
	OPEN(116,FILE=location2//'mean3D_CRect.dat',FORM='unformatted',STATUS='unknown',ACCESS='direct',RECL=4)
	OPEN(117,FILE=location2//'mean3D_ERect.dat',FORM='unformatted',STATUS='unknown',ACCESS='direct',RECL=4)

	OPEN(124,FILE=location2//'mean3Dfd_total.dat',FORM='unformatted',STATUS='unknown',ACCESS='direct',RECL=4)
	OPEN(125,FILE=location2//'mean3Dfd_WRect.dat',FORM='unformatted',STATUS='unknown',ACCESS='direct',RECL=4)
	OPEN(126,FILE=location2//'mean3Dfd_CRect.dat',FORM='unformatted',STATUS='unknown',ACCESS='direct',RECL=4)
	OPEN(127,FILE=location2//'mean3Dfd_ERect.dat',FORM='unformatted',STATUS='unknown',ACCESS='direct',RECL=4)

	write(100,800) 'duration','f_input1_3D','f_input2_3D','f_input3_3D','f_input4_3D','f_input_total'
	write(101,800) 'duration','f_inputW1_3D','f_inputW2_3D','f_inputW3_3D','f_inputW4_3D','f_input_totalW'
	write(102,800) 'duration','f_inputC1_3D','f_inputC2_3D','f_inputC3_3D','f_inputC4_3D','f_input_totalC'
	write(103,800) 'duration','f_inputE1_3D','f_inputE2_3D','f_inputE3_3D','f_inputE4_3D','f_input_totalE'

	write(110,800) 'duration','eamean_3D','eamean_w_3D','eamean_c_3D','eamean_e_3D'
	write(111,800) 'duration','fdmean_3D','fdmean_w_3D','fdmean_c_3D','fdmean_e_3D'

!------------------------Prepare output file
      	undef=1e+15
!      

		write(out_H,1001) "Hd_", duration(c), ".dat"
		write(out_V,1001) "Vd_", duration(c), ".dat"
		!write(*,*) out_H, out_V

      		OPEN(1,FILE = location2//out_H,FORM='unformatted',&
	      	STATUS='unknown',ACCESS='direct',iostat=ioerr,RECL=4*m*n_tot)
					if (ioerr /= 0) then
					write(*,*) "ioerr1=", ioerr
				   	!stop
					end if

      		OPEN(2,FILE = location2//out_V,FORM='unformatted',&
	      	STATUS='unknown',ACCESS='direct',iostat=ioerr,RECL=4*m*n_tot)
					if (ioerr /= 0) then
					write(*,*) "ioerr2=", ioerr
				   	!stop
					end if
!----------------------------------------------------------
		!write(*,*) "PO=", PO   !! what is PO?????
!!		write(101,'(i5,5e14.6)') PO,f_input1w_3D,f_input2w_3D,f_input3w_3D,f_input4w_3D,f_input_totalw (I change it by myself)
		write(100,'(A2,5e14.6)') duration(c),f_input1_3D,f_input2_3D,f_input3_3D,f_input4_3D,f_input_total         
   		write(101,'(A2,5e14.6)') duration(c),f_input1w_3D,f_input2w_3D,f_input3w_3D,f_input4w_3D,f_input_totalw
   		write(102,'(A2,5e14.6)') duration(c),f_input1c_3D,f_input2c_3D,f_input3c_3D,f_input4c_3D,f_input_totalc
   		write(103,'(A2,5e14.6)') duration(c),f_input1e_3D,f_input2e_3D,f_input3e_3D,f_input4e_3D,f_input_totale
   
		write(110,'(A2,4e14.6)') duration(c),eamean_3D,eamean_w_3D,eamean_c_3D,eamean_e_3D
		write(111,'(A2,4e14.6)') duration(c),fdmean_3D,fdmean_w_3D,fdmean_c_3D,fdmean_e_3D

   		write(104,Rec=(tt-1)*5+1) REAL(f_input_total)  !rec chrjur hesab shode???
   		write(104,Rec=(tt-1)*5+2) REAL(f_input1_3D)
   		write(104,Rec=(tt-1)*5+3) REAL(f_input2_3D)
   		write(104,Rec=(tt-1)*5+4) REAL(f_input3_3D)
   		write(104,Rec=(tt-1)*5+5) REAL(f_input4_3D)
 
   		write(105,Rec=(tt-1)*5+1) REAL(f_input_totalW)
   		write(105,Rec=(tt-1)*5+2) REAL(f_input1W_3D)
   		write(105,Rec=(tt-1)*5+3) REAL(f_input2W_3D)
   		write(105,Rec=(tt-1)*5+4) REAL(f_input3W_3D)
   		write(105,Rec=(tt-1)*5+5) REAL(f_input4W_3D)

   		write(106,Rec=(tt-1)*5+1) REAL(f_input_totalC)
   		write(106,Rec=(tt-1)*5+2) REAL(f_input1C_3D)
   		write(106,Rec=(tt-1)*5+3) REAL(f_input2C_3D)
   		write(106,Rec=(tt-1)*5+4) REAL(f_input3C_3D)
   		write(106,Rec=(tt-1)*5+5) REAL(f_input4C_3D)
   
   
   		write(107,Rec=(tt-1)*5+1) REAL(f_input_totalE)
   		write(107,Rec=(tt-1)*5+2) REAL(f_input1E_3D)
   		write(107,Rec=(tt-1)*5+3) REAL(f_input2E_3D)
   		write(107,Rec=(tt-1)*5+4) REAL(f_input3E_3D)
   		write(107,Rec=(tt-1)*5+5) REAL(f_input4E_3D)

   		write(114,Rec=tt) REAL(eamean_3D)
   		write(115,Rec=tt) REAL(eamean_w_3D)
   		write(116,Rec=tt) REAL(eamean_c_3D)
   		write(117,Rec=tt) REAL(eamean_e_3D)

   		write(124,Rec=tt) REAL(fdmean_3D)
   		write(125,Rec=tt) REAL(fdmean_w_3D)
   		write(126,Rec=tt) REAL(fdmean_c_3D)
   		write(127,Rec=tt) REAL(fdmean_e_3D)

 !####################################################################################
 !_________________________________ WRITE OUTPUT ______________________________________
 !#####################################################################################

 		jj=1 
 
  		DO  k=2,9
    			DO  j=n_st,n_end
     			DO  i=1,m
       		out(i,j)=HMeddy_act(i,j,k)
      			END DO
   			END DO
   		WRITE(1,REC=jj*8+2-k) ((out(i,j),i=1,m),j=n_st,n_end)     !! how is the rec?
   		END DO
  		WRITE(2,REC=jj) ((SNGL(VMeddy_act(i,j)),i=1,m),j=n_st,n_end)
   		jj=jj+1
    
   
   		DO  k=2,9
    			DO  j=n_st,n_end
     			DO  i=1,m
       		out(i,j)= HMf1(i,j,k)
      			END DO
   			END DO
    		WRITE(1,REC=jj*8+2-k) ((out(i,j),i=1,m),j=n_st,n_end)
    		END DO
  		WRITE(2,REC=jj) ((SNGL(VMf1(i,j)),i=1,m),j=n_st,n_end) 
		jj=jj+1

		DO  k=2,9
   			DO  j=n_st,n_end
     			DO  i=1,m
       		out(i,j)= HMf2(i,j,k)
      			END DO
   			END DO
    		WRITE(1,REC=jj*8+2-k) ((out(i,j),i=1,m),j=n_st,n_end)
    		END DO
  		WRITE(2,REC=jj) ((SNGL(VMf2(i,j)),i=1,m),j=n_st,n_end)
		jj=jj+1


		DO  k=2,9
    			DO  j=n_st,n_end
     			DO  i=1,m
       		out(i,j)= HMf3(i,j,k)
      			END DO
   			END DO
    		WRITE(1,REC=jj*8+2-k) ((out(i,j),i=1,m),j=n_st,n_end)
    		END DO
 		WRITE(2,REC=jj) ((SNGL(VMf3(i,j)),i=1,m),j=n_st,n_end)
		jj=jj+1

		DO  k=2,9
    			DO  j=n_st,n_end
     			DO  i=1,m
       		out(i,j)= HMdiv(i,j,k)
      			END DO
   			END DO
    		WRITE(1,REC=jj*8+2-k) ((out(i,j),i=1,m),j=n_st,n_end)
    		END DO
 		WRITE(2,REC=jj) ((SNGL(VMdiv(i,j)),i=1,m),j=n_st,n_end)
		jj=jj+1

		DO  k=2,9
    			DO  j=n_st,n_end
     			DO  i=1,m
       		out(i,j)= HMdivyp(i,j,k)
      			END DO
   			END DO
    		WRITE(1,REC=jj*8+2-k) ((out(i,j),i=1,m),j=n_st,n_end)
    		END DO
 		WRITE(2,REC=jj) ((SNGL(VMdivyp(i,j)),i=1,m),j=n_st,n_end)
		jj=jj+1


		DO  k=2,9
    			DO  j=n_st,n_end
    	 		DO  i=1,m
     	  		out(i,j)= HMdivxyp(i,j,k)
         		END DO
  	    		END DO
   	    	WRITE(1,REC=jj*8+2-k) ((out(i,j),i=1,m),j=n_st,n_end)
    		END DO
    		WRITE(2,REC=jj) ((SNGL(VMdivxyp(i,j)),i=1,m),j=n_st,n_end)
		jj=jj+1

	c=c+1 
	end do !! this do is for duration c  

END program

!==========================================================================
!SUBROUTINE FOR COMPUTING HORIZONTAL AVERAGING
SUBROUTINE Averaging(ave_check,inparam,n_st,n_end,HMparam,VMparam)

PARAMETER (m=288,n=145)

INTEGER Dphi,Dlanda

LOGICAL ave_check

integer, PARAMETER :: Dlat=6,DLon=12    !!! what is this??

REAL*8 inparam(m,n,10),HMparam(m,n,10),VMparam(m,n)
REAL*8 inparamE(1-Dlon:m+Dlon,n)
REAL*8 phi_a(n),Area(n)

REAL*8 delta_phi,delta_landa,tem

COMMON /Average_S/ phi_a,delta_phi,delta_landa

!COMPUTING SPACE AVERAGING
DO j=1,n
  Area(j)=0.D0
ENDDO

DO j=1+Dlat,n-Dlat

  DO jj=j-Dlat,j+Dlat-1
    DO i=i-Dlon,i+Dlon-1
      Area(j)=Area(j)+DCOS((phi_a(jj)+phi_a(jj+1))/2.D0)*delta_phi*delta_landa
    ENDDO  
  ENDDO

ENDDO    

DO k=2,9
  
  DO j=n_st,n_end
    DO i=1,Dlon
      inparamE(1-i,j)=inparam(m-i+1,j,k)
      inparamE(m+i,j)=inparam(i,j,k)
    ENDDO
    DO i=1,m
      inparamE(i,j)=inparam(i,j,k)
    ENDDO
  ENDDO

  DO j=1,n
    DO i=1,m
      HMparam(i,j,k)=0.D0
    ENDDO
  ENDDO    
  
  DO j=n_st+Dlat,n_end-Dlat
    DO i=1,m

      DO jj=j-Dlat,j+Dlat-1
        DO ii=i-Dlon,i+Dlon-1
          tem=(inparamE(ii,jj)+inparamE(ii+1,jj)+inparamE(ii,jj+1)+inparamE(ii+1,jj+1))/4.D0
          HMparam(i,j,k)=HMparam(i,j,k)+tem*DCOS((phi_a(jj)+phi_a(jj+1))/2.D0)*delta_phi*delta_landa
        ENDDO
      ENDDO    
      HMparam(i,j,k)=HMparam(i,j,k)/Area(j)

    ENDDO
  ENDDO

ENDDO  

!Computing vertical mean
   DO  j=n_st,n_end
     DO  i=1,m
         VMparam(i,j)=0.D0
	 END DO
   ENDDO

  IF (ave_check) THEN
     DO  j=n_st,n_end
      DO  i=1,m
	    DO  k=2,4
	      VMparam(i,j)=VMparam(i,j)+HMparam(i,j,k) 
	    END DO
	  END DO
     ENDDO
   ELSE 
     DO  j=n_st,n_end
      DO  i=1,m
	    DO  k=2,4
	      VMparam(i,j)=VMparam(i,j)+inparam(i,j,k) 
	    END DO
	  END DO
     ENDDO
   ENDIF

  DO  j=n_st,n_end
    DO  i=1,m
      VMparam(i,j)=VMparam(i,j)/3.D0
    END DO
  ENDDO

END SUBROUTINE

!==========================================================================
!SUBROUTINE FOR COMPUTING HORIZONTAL AVERAGING TWO DIMENSION
SUBROUTINE Averaging2D(inparam,n_st,n_end,HMparam)

PARAMETER (m=288,n=145)

INTEGER Dphi,Dlanda

integer, PARAMETER :: Dlat=10,DLon=20

REAL*8 inparam(m,n),HMparam(m,n)
REAL*8 inparamE(1-Dlon:m+Dlon,n)
REAL*8 phi_a(n),Area(n)
REAL*8 delta_phi,delta_landa,tem
COMMON /Average_S/ phi_a,delta_phi,delta_landa

!COMPUTING SPACE AVERAGING
DO j=1,n
  Area(j)=0.D0
ENDDO

DO j=1+Dlat,n-Dlat

  DO jj=j-Dlat,j+Dlat-1
    DO i=i-Dlon,i+Dlon-1
      Area(j)=Area(j)+DCOS((phi_a(jj)+phi_a(jj+1))/2.D0)*delta_phi*delta_landa
    ENDDO  
  ENDDO

ENDDO    

!PREPARE EXTENDED PARAMETER IN LONGITUDE
DO j=n_st,n_end
  DO i=1,Dlon
    inparamE(1-i,j)=inparam(m-i+1,j)
    inparamE(m+i,j)=inparam(i,j)
  ENDDO
  DO i=1,m
    inparamE(i,j)=inparam(i,j)
  ENDDO
ENDDO

DO j=1,n
  DO i=1,m
    HMparam(i,j)=0.D0
  ENDDO
ENDDO    
  
DO j=n_st+Dlat,n_end-Dlat
  DO i=1,m

      DO jj=j-Dlat,j+Dlat-1
        DO ii=i-Dlon,i+Dlon-1
          tem=(inparamE(ii,jj)+inparamE(ii+1,jj)+inparamE(ii,jj+1)+inparamE(ii+1,jj+1))/4.D0
          HMparam(i,j)=HMparam(i,j)+tem*DCOS((phi_a(jj)+phi_a(jj+1))/2.D0)*delta_phi*delta_landa
        ENDDO
      ENDDO    
      HMparam(i,j)=HMparam(i,j)/Area(j)
  ENDDO
ENDDO

RETURN
END SUBROUTINE
