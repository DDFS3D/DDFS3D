!*********************************************************************************************
!Part 1: Solve for the displacement discontinuities at the centriod of each element.
!*********************************************************************************************	

program DDM3D_qua_part1
 
use prog_bar
use subs_dd3d_qua_num
use solver
use omp_lib

implicit none 

integer:: num_threads = 8
real*8 ::E,v,KN,KS,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4
real*8 ::gpx,gpy,gpz,COEFF1(3,3),COEFF2(3,3),IEvect(3,3),JEvect(3,3)
integer::M,NELE,NNOD,AD,AD1,TY,i,j,k,ii,jj,RES(8),iter
real*8::t1,t2,D1,D2,ratio,L(8),GPXX,GPYY,GPZZ,pxx,pyy,pzz,pxy,pxz,pyz,EI(3,3),Evect(3,3),trac(3) 
real*8,allocatable::NOD(:,:),VOB(:,:),COE(:,:),X(:),B(:),nor(:,:),gx(:),gy(:),gz(:)
integer,allocatable::ELE(:,:),E_type(:),SOB(:,:)
type (type_prog):: progress
Character*100::Title

!*****************************************************************************
!$ call omp_set_num_threads(num_threads)

write (*,*)'*************************************************************************'
write (*,*)'************		    DDM3D_qua_part1		        **********'
write (*,*)'************	    developed by Wong and Cui (2021)            **********'
write (*,*)'************	       DDFS3D (version 1.0), 31-5-2021		**********'
write (*,*)'*************************************************************************'

call cpu_time(t1)

open (21,file='input_part1.txt')

read (21,'(A100)') Title
read (21,*) E,v,KS,KN
read (21,*) pxx,pyy,pzz,pxy,pxz,pyz                
read (21,*) NELE,NNOD

allocate (NOD(NNOD,3),VOB(NELE,3),COE(3*NELE,3*NELE),X(3*NELE),B(3*NELE),ELE(NELE,4),E_type(NELE),SOB(NELE,3),nor(NELE,3),gx(NELE),gy(NELE),gz(NELE))

write (*,*)'Reading element information ...'
do i=1,NELE
    read (21,*) AD,ELE(i,1:4)
    call progress%input(NELE,20)
    call progress%output(i)
end do 

write (*,*)'Reading node information ...'
do i=1,NNOD
    read (21,*) AD,NOD(i,1:3)
    call progress%input(NNOD,20)
    call progress%output(i)
end do

write (*,*)'Reading boundary condition ...'
do i=1,NELE
    read (21,*) AD,E_type(i),SOB(i,1:3),VOB(i,1:3)
	call progress%input(NELE,20)
    call progress%output(i)
end do 

close(21)

!$omp parallel do default(shared)
do i=1,3*NELE
       X(i)=0.
       B(i)=0.
    do j=1,3*NELE
        COE(i,j)=0.
    end do 
end do 
!$omp end parallel do 
ITER=0
AD1=0
!*********************************************************
write (*,*)'Calculating influence coefficients ...'
!$omp parallel do default(none) shared(COE,nor,AD1,progress,gx,gy,gz) shared(NELE,NOD,ELE,SOB,E,v,Pxx,Pyy,Pzz,Pxy,Pxz,Pyz) &
!$omp private (x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,gpx,gpy,gpz,xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4, &
!$omp		  gpxx,gpyy,gpzz,l,AD,res,D1,D2,ratio,M,IEvect,JEvect,COEFF1,COEFF2,EI,Trac,Evect)

do 100 i=1,NELE
	AD1=AD1+1
	!$omp critical 
	call progress%input(NELE,20)
	call progress%output(AD1)
	!$omp end critical    
	
	x1=NOD(ELE(I,1),1)
	y1=NOD(ELE(I,1),2)
	z1=NOD(ELE(I,1),3)
       
	x2=NOD(ELE(I,2),1)
	y2=NOD(ELE(I,2),2)
	z2=NOD(ELE(I,2),3)
       
	x3=NOD(ELE(I,3),1)
	y3=NOD(ELE(I,3),2)
	z3=NOD(ELE(I,3),3)
       
	x4=NOD(ELE(I,4),1)
	y4=NOD(ELE(I,4),2)
	z4=NOD(ELE(I,4),3)
       
	gpx=(x1+x2+x3+x4)/4.
	gpy=(y1+y2+y3+y4)/4.
	gpz=(z1+z2+z3+z4)/4.
       
	gx(i)=gpx
	gy(i)=gpy
	gz(i)=gpz
       
    do 200 j=1,NELE

        xx1=NOD(ELE(J,1),1)
        yy1=NOD(ELE(J,1),2)
        zz1=NOD(ELE(J,1),3)
       
        xx2=NOD(ELE(J,2),1)
        yy2=NOD(ELE(J,2),2)
        zz2=NOD(ELE(J,2),3)
       
        xx3=NOD(ELE(J,3),1)
        yy3=NOD(ELE(J,3),2)
        zz3=NOD(ELE(J,3),3)
       
        xx4=NOD(ELE(J,4),1)
        yy4=NOD(ELE(J,4),2)
        zz4=NOD(ELE(J,4),3)
       
        gpxx=(xx1+xx2+xx3+xx4)/4.
        gpyy=(yy1+yy2+yy3+yy4)/4.
        gpzz=(zz1+zz2+zz3+zz4)/4.
       
!*****************************************************
! Determine the number of Gaussian points
!*****************************************************
        l(1)=sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2) 
		l(2)=sqrt((x3-x2)**2+(y3-y2)**2+(z3-z2)**2) 
		l(3)=sqrt((x4-x3)**2+(y4-y3)**2+(z4-z3)**2) 
		l(4)=sqrt((x1-x4)**2+(y1-y4)**2+(z1-z4)**2) 
		l(5)=sqrt((xx2-xx1)**2+(yy2-yy1)**2+(zz2-zz1)**2)
		l(6)=sqrt((xx3-xx2)**2+(yy3-yy2)**2+(zz3-zz2)**2) 
		l(7)=sqrt((xx4-xx3)**2+(yy4-yy3)**2+(zz4-zz3)**2) 
		l(8)=sqrt((xx1-xx4)**2+(yy1-yy4)**2+(zz1-zz4)**2)  
 
		do ii=1,8
			AD=0
			do jj=1,8
				if ((l(ii)-l(jj))>=-1.e-8)then
					AD=AD+1
			end if 
		end do  
			res(ii)=AD
		end do 

		do ii=1,8
			if (res(ii)==8)then
				D1=l(ii)
			end if 
		end do 

		D2=sqrt((gpx-gpxx)**2+(gpy-gpyy)**2+(gpz-gpzz)**2)

		ratio=D2/D1

		if (ratio<=0.5)then
			M=24
		elseif (0.5<ratio.and.ratio<3)then
			M=ceiling(27-7.6*ratio)
		else 
			M=5
		end if 
		
		if (i==j)then
			M=18
		end if 
!************************************************************************
! Calculate coefficient matirx
!************************************************************************        
        call convert (x1,y1,z1,x2,y2,z2,x3,y3,z3,IEvect)
        call convert (xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,JEvect)
      
        if (SOB(i,1)==0.and.SOB(i,2)==0.and.SOB(i,3)==0) then 
            if (i/=j)then
                call stress_regular_dd3d_qua_num_part1 (xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4,gpx,gpy,gpz,E,v,M,IEvect,JEvect,COEFF1)
                COE(3*(i-1)+1,3*(j-1)+1)=COEFF1(1,1)
                COE(3*(i-1)+1,3*(j-1)+2)=COEFF1(1,2)
                COE(3*(i-1)+1,3*(j-1)+3)=COEFF1(1,3)           
                COE(3*(i-1)+2,3*(j-1)+1)=COEFF1(2,1)
                COE(3*(i-1)+2,3*(j-1)+2)=COEFF1(2,2)
                COE(3*(i-1)+2,3*(j-1)+3)=COEFF1(2,3) 
                COE(3*(i-1)+3,3*(j-1)+1)=COEFF1(3,1)
                COE(3*(i-1)+3,3*(j-1)+2)=COEFF1(3,2)
                COE(3*(i-1)+3,3*(j-1)+3)=COEFF1(3,3) 
		    else  
                call stress_singular_dd3d_qua_num (xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4,E,v,M,IEvect,COEFF1)
                COE(3*(i-1)+1,3*(j-1)+1)=COEFF1(1,1)
                COE(3*(i-1)+1,3*(j-1)+2)=COEFF1(1,2)
                COE(3*(i-1)+1,3*(j-1)+3)=COEFF1(1,3)           
                COE(3*(i-1)+2,3*(j-1)+1)=COEFF1(2,1)
                COE(3*(i-1)+2,3*(j-1)+2)=COEFF1(2,2)
                COE(3*(i-1)+2,3*(j-1)+3)=COEFF1(2,3) 
                COE(3*(i-1)+3,3*(j-1)+1)=COEFF1(3,1)
                COE(3*(i-1)+3,3*(j-1)+2)=COEFF1(3,2)
                COE(3*(i-1)+3,3*(j-1)+3)=COEFF1(3,3)  
            end if 
        else if (SOB(i,1)==1.and.SOB(i,2)==1.and.SOB(i,3)==1)then
            if (i/=j)then  
                call dis_regular_dd3d_qua_num (xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4,gpx,gpy,gpz,E,v,M,IEvect,JEvect,COEFF1) 
                COE(3*(i-1)+1,3*(j-1)+1)=COEFF1(1,1)
                COE(3*(i-1)+1,3*(j-1)+2)=COEFF1(1,2)
                COE(3*(i-1)+1,3*(j-1)+3)=COEFF1(1,3)           
                COE(3*(i-1)+2,3*(j-1)+1)=COEFF1(2,1)
                COE(3*(i-1)+2,3*(j-1)+2)=COEFF1(2,2)
                COE(3*(i-1)+2,3*(j-1)+3)=COEFF1(2,3) 
                COE(3*(i-1)+3,3*(j-1)+1)=COEFF1(3,1)
                COE(3*(i-1)+3,3*(j-1)+2)=COEFF1(3,2)
                COE(3*(i-1)+3,3*(j-1)+3)=COEFF1(3,3) 
            else 
                call dis_singular_dd3d_qua_num (xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4,E,v,M,IEvect,COEFF1)
                COE(3*(i-1)+1,3*(j-1)+1)=COEFF1(1,1)
                COE(3*(i-1)+1,3*(j-1)+2)=COEFF1(1,2)
                COE(3*(i-1)+1,3*(j-1)+3)=COEFF1(1,3)           
                COE(3*(i-1)+2,3*(j-1)+1)=COEFF1(2,1)
                COE(3*(i-1)+2,3*(j-1)+2)=COEFF1(2,2)
                COE(3*(i-1)+2,3*(j-1)+3)=COEFF1(2,3) 
                COE(3*(i-1)+3,3*(j-1)+1)=COEFF1(3,1)
                COE(3*(i-1)+3,3*(j-1)+2)=COEFF1(3,2)
                COE(3*(i-1)+3,3*(j-1)+3)=COEFF1(3,3) 
            end if 
        else
            if (i/=j)then
                call stress_regular_dd3d_qua_num_part1 (xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4,gpx,gpy,gpz,E,v,M,IEvect,JEvect,COEFF1)  
                call          dis_regular_dd3d_qua_num (xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4,gpx,gpy,gpz,E,v,M,IEvect,JEvect,COEFF2) 
          
                if (SOB(i,1)==0)then 
                    COE(3*(i-1)+1,3*(j-1)+1)=COEFF1(1,1)
                    COE(3*(i-1)+1,3*(j-1)+2)=COEFF1(1,2)
                    COE(3*(i-1)+1,3*(j-1)+3)=COEFF1(1,3) 
                else 
                    COE(3*(i-1)+1,3*(j-1)+1)=COEFF2(1,1)
                    COE(3*(i-1)+1,3*(j-1)+2)=COEFF2(1,2)
                    COE(3*(i-1)+1,3*(j-1)+3)=COEFF2(1,3)   
                end if 
                if (SOB(i,2)==0)then 
                    COE(3*(i-1)+2,3*(j-1)+1)=COEFF1(2,1)
                    COE(3*(i-1)+2,3*(j-1)+2)=COEFF1(2,2)
                    COE(3*(i-1)+2,3*(j-1)+3)=COEFF1(2,3) 
                else 
                    COE(3*(i-1)+2,3*(j-1)+1)=COEFF2(2,1)
                    COE(3*(i-1)+2,3*(j-1)+2)=COEFF2(2,2)
                    COE(3*(i-1)+2,3*(j-1)+3)=COEFF2(2,3)   
                end if 
                if (SOB(i,3)==0)then 
                    COE(3*(i-1)+3,3*(j-1)+1)=COEFF1(3,1)
                    COE(3*(i-1)+3,3*(j-1)+2)=COEFF1(3,2)
                    COE(3*(i-1)+3,3*(j-1)+3)=COEFF1(3,3) 
                else 
                    COE(3*(i-1)+3,3*(j-1)+1)=COEFF2(3,1)
                    COE(3*(i-1)+3,3*(j-1)+2)=COEFF2(3,2)
                    COE(3*(i-1)+3,3*(j-1)+3)=COEFF2(3,3)   
                end if 
            else 
                call stress_singular_dd3d_qua_num (xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4,E,v,M,IEvect,COEFF1)
                call    dis_singular_dd3d_qua_num (xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4,E,v,M,IEvect,COEFF2)
               
                if (SOB(i,1)==0)then 
                    COE(3*(i-1)+1,3*(j-1)+1)=COEFF1(1,1)
                    COE(3*(i-1)+1,3*(j-1)+2)=COEFF1(1,2)
                    COE(3*(i-1)+1,3*(j-1)+3)=COEFF1(1,3) 
                else 
                    COE(3*(i-1)+1,3*(j-1)+1)=COEFF2(1,1)
                    COE(3*(i-1)+1,3*(j-1)+2)=COEFF2(1,2)
                    COE(3*(i-1)+1,3*(j-1)+3)=COEFF2(1,3)   
                end if  
                if (SOB(i,2)==0)then 
                    COE(3*(i-1)+2,3*(j-1)+1)=COEFF1(2,1)
                    COE(3*(i-1)+2,3*(j-1)+2)=COEFF1(2,2)
                    COE(3*(i-1)+2,3*(j-1)+3)=COEFF1(2,3) 
                else 
                    COE(3*(i-1)+2,3*(j-1)+1)=COEFF2(2,1)
                    COE(3*(i-1)+2,3*(j-1)+2)=COEFF2(2,2)
                    COE(3*(i-1)+2,3*(j-1)+3)=COEFF2(2,3)   
                end if  
                if (SOB(i,3)==0)then 
                    COE(3*(i-1)+3,3*(j-1)+1)=COEFF1(3,1)
                    COE(3*(i-1)+3,3*(j-1)+2)=COEFF1(3,2)
                    COE(3*(i-1)+3,3*(j-1)+3)=COEFF1(3,3) 
                else 
                    COE(3*(i-1)+3,3*(j-1)+1)=COEFF2(3,1)
                    COE(3*(i-1)+3,3*(j-1)+2)=COEFF2(3,2)
                    COE(3*(i-1)+3,3*(j-1)+3)=COEFF2(3,3)   
			    end if
		    end if 
	    end if
!*************************************************************************************
! Transform the initial stress field to the local coordinate system of each element.
!*************************************************************************************    
		do ii=1,3
			do jj=1,3
				Evect(ii,jj)=0.
				if (ii==jj)then
					Evect(ii,jj)=1.
				end if 
			end do
		end do 

		EI(1,1)=Evect(1,1)*IEvect(1,1)+Evect(1,2)*IEvect(1,2)+Evect(1,3)*IEvect(1,3)
		EI(2,1)=Evect(2,1)*IEvect(1,1)+Evect(2,2)*IEvect(1,2)+Evect(2,3)*IEvect(1,3)
		EI(3,1)=Evect(3,1)*IEvect(1,1)+Evect(3,2)*IEvect(1,2)+Evect(3,3)*IEvect(1,3)
      
		EI(1,2)=Evect(1,1)*IEvect(2,1)+Evect(1,2)*IEvect(2,2)+Evect(1,3)*IEvect(2,3)
		EI(2,2)=Evect(2,1)*IEvect(2,1)+Evect(2,2)*IEvect(2,2)+Evect(2,3)*IEvect(2,3)
		EI(3,2)=Evect(3,1)*IEvect(2,1)+Evect(3,2)*IEvect(2,2)+Evect(3,3)*IEvect(2,3) 

		EI(1,3)=Evect(1,1)*IEvect(3,1)+Evect(1,2)*IEvect(3,2)+Evect(1,3)*IEvect(3,3)
		EI(2,3)=Evect(2,1)*IEvect(3,1)+Evect(2,2)*IEvect(3,2)+Evect(2,3)*IEvect(3,3)
		EI(3,3)=Evect(3,1)*IEvect(3,1)+Evect(3,2)*IEvect(3,2)+Evect(3,3)*IEvect(3,3)
        
		Trac(1)=Pxx*EI(1,3)+Pxy*EI(2,3)+Pxz*EI(3,3)
		Trac(2)=Pxy*EI(1,3)+Pyy*EI(2,3)+Pyz*EI(3,3)
		Trac(3)=Pxz*EI(1,3)+Pyz*EI(2,3)+Pzz*EI(3,3)
     
		nor(i,1)=Trac(1)*EI(1,1)+Trac(2)*EI(2,1)+Trac(3)*EI(3,1)
		nor(i,2)=Trac(1)*EI(1,2)+Trac(2)*EI(2,2)+Trac(3)*EI(3,2)
		nor(i,3)=Trac(1)*EI(1,3)+Trac(2)*EI(2,3)+Trac(3)*EI(3,3)     
!*************************************************************************************        
200 continue 
100 continue 
!$omp end parallel do 
!*****************************************************************************    
! Add the normal and tangential stiffness coefficients 
!*****************************************************************************
!$omp parallel do default(none) shared(COE,NELE) shared(KS,KN,E_Type)
do i=1,NELE 
    if (E_type(i)==2)then
      COE(3*i-2,3*i-2)=COE(3*i-2,3*i-2)+KS
      COE(3*i-1,3*i-1)=COE(3*i-1,3*i-1)+KS
      COE(3*i  ,3*i  )=COE(3*i  ,3*i  )+KN
    end if 
end do
!$omp end parallel do
!**************************************************************
! Evaluate the righ-hand side vector
!**************************************************************
!$omp parallel do default(none) shared(B,NELE,SOB) shared(VOB,nor)
do i=1,NELE
    if (SOB(i,1)==0)then
        B(3*i-2)=VOB(i,1)-nor(i,1)
    else
        B(3*i-2)=VOB(i,1)
    end if 
    if (SOB(i,2)==0)then
        B(3*i-1)=VOB(i,2)-nor(i,2)
    else
        B(3*i-1)=VOB(i,2)
    end if 
    if (SOB(i,3)==0)then
        B(3*i)=VOB(i,3)-nor(i,3)
    else
        B(3*i)=VOB(i,3)
    end if    
end do 
!$omp end parallel do 
!**************************************************************
! Solve the system equations
!**************************************************************
call GS(COE,B,X,ITER) 

call cpu_time(t2)

write (*,*)
if (iter<500)then
	AD=0
	!$omp parallel default(none) shared(AD)
		!$omp master
			!$ if (omp_in_parallel())then
				!$ AD=1
			!$ end if 
		!$omp end master	
	!$omp end parallel
	if (AD==0)then
		write (*,'(a18,f7.2,a2)')' Computing time = ',t2-t1,' s'
	else 
		!$ write (*,'(a18,f7.2,a2)')' Computing time = ',(t2-t1)/num_threads,' s'
	end if 
else 
    write (*,*)'Times of iteration exceed the maximum value (500),iteration fails'
	stop
end if 
!**************************************************************
! Write output files
!**************************************************************
open (41,file='output_part1.txt')
write (41,*)'*************************************************************************'
write (41,*)'************		    DDM3D_qua_part1		        **********'
write (41,*)'************	    developed by Wong and Cui (2021)            **********'
write (41,*)'************	       DDFS3D (version 1.0), 31-5-2021		**********'
write (41,*)'*************************************************************************'
write (41,'(2A)')' Title: ', adjustl(Title)
write (41,'(A9,es12.4,A16,f10.4,A5,es12.4,A5,es12.4)')'Modulus=',E,'Poisson ratio=',v,'  KS=',KS,'  KN=',KN
write (41,*)'Initial stress field:'
write (41,'(A5,es12.4,A5,es12.4,A5,es12.4)')'Pxx=',pxx,'Pyy=',pyy,'Pzz=',pzz
write (41,'(A5,es12.4,A5,es12.4,A5,es12.4)')'Pxy=',pxy,'Pxz=',pxz,'Pyz=',pyz
write (41,'(A20,I10)')'Number of elements:',NELE
write (41,'(A20,I10)')'Number of nodes   :',NNOD
write (41,'(A20,I10)')'Times of iteration:',ITER
if (AD==0)then
	write (41,'(A20,f10.2,A2)')'Computing time    :',t2-t1,'s'
else 
	!$ write (41,'(A20,f10.2,A2)')'Computing time    :',(t2-t1)/num_threads,' s'
end if 
write (41,'(A83)')'Coordinates of center of gravity and displacement discontinuities of each element:'
write (41,'(A83)')' No.	     Gx            Gy            Gz            Dx            Dy           Dz '
do i=1,NELE
  write (41,'(I5,6es14.4)')i,gx(i),gy(i),gz(i),X(3*i-2),X(3*i-1),X(3*i)
end do 
close(41)

open (42,file='D.txt')
do i=1,NELE
   write (42,'(I8,3ES16.6)')i,X(3*i-2),X(3*i-1),X(3*i)
end do 
close (42)

end program DDM3D_qua_part1
