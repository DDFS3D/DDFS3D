!***************************************************************************************************************
!Part 3: Calculate the stresses and displacements at arbitrary points in the prescribed coordinate system(s).
!***************************************************************************************************************

program FSM3D_qua_part3
 
use prog_bar
use subs_fs3d_qua_num
    
implicit none 

real*8 ::E,v,KN,KS,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,ze,xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4
real*8 ::gpx,gpy,gpz,COEFF1(6,3),COEFF2(3,3),IEvect(3,3),JEvect(3,3)
integer::M,NELE,NNOD,AD,TY,i,j,k,ii,jj,RES(6),NOP
real*8::t1,t2,D1,D2,ratio,L(6),GPXX,GPYY,GPZZ,pxx,pyy,pzz,pxy,pxz,pyz,EI(3,3),Evect(3,3),trac(3) 
real*8,allocatable::NOD(:,:),VOB(:,:),nor(:,:),X(:),CFP(:,:),LCS(:,:)
real*8,allocatable::PPxx(:),PPxy(:),PPxz(:),PPyy(:),PPyz(:),PPzz(:),Ux(:),Uy(:),Uz(:)
integer,allocatable::ELE(:,:),SOB(:,:),E_type(:)
type (type_prog)::progress
character*100::Title

!*****************************************************************************
write (*,*)'*************************************************************************'
write (*,*)'************		    FSM3D_qua_part3		        **********'
write (*,*)'************	    developed by Wong and Cui (2021)            **********'
write (*,*)'************	       DDFS3D (version 1.0), 31-5-2021		**********'
write (*,*)'*************************************************************************'

call cpu_time(t1)

open (21,file='input_part1.txt')
read (21,'(A100)')Title
read (21,*) E,V,KS,KN
read (21,*) pxx,pyy,pzz,pxy,pxz,pyz                
read (21,*) NELE,NNOD

allocate(NOD(NNOD,3),VOB(NELE,3),X(3*NELE),ELE(NELE,4),SOB(NELE,3),E_type(NELE))

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

write (*,*)'Reading fictitious stresses ...'
open (22,file='P.txt')
do i=1,NELE
    read (22,*)AD,X(3*i-2:3*i)
	call progress%input(NELE,20)
    call progress%output(i)
end do 
close (22)

write (*,*)'Reading points and coordinate system(s) where S and U are to be calculated ...'
open (23,file='input_part3.txt')
read (23,*)NOP

allocate(CFP(NOP,3),LCS(NOP,9),PPxx(NOP),PPxy(NOP),PPxz(NOP),PPyy(NOP),PPyz(NOP),PPzz(NOP),Ux(NOP),Uy(NOP),Uz(NOP),nor(NOP,6))
do i=1,NOP
  read(23,*)AD,CFP(i,1:3),LCS(i,1:9)
end do 
close (23)

do i=1,NOP
    PPxx(i)=0.
    PPxy(i)=0.
    PPxz(i)=0.
    PPyy(i)=0.
    PPyz(i)=0.
    PPzz(i)=0.
    Ux(i)=0.
    Uy(i)=0.
    Uz(i)=0.
end do 

write (*,*)'Calculating S and U ...'
do 100 i=1,NOP     
	call progress%input(NOP,20)
    call progress%output(i)
	
	gpx=CFP(i,1)
	gpy=CFP(i,2)
	gpz=CFP(i,3)
	
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
		l(1)=sqrt((xx2-xx1)**2+(yy2-yy1)**2+(zz2-zz1)**2)
		l(2)=sqrt((xx3-xx2)**2+(yy3-yy2)**2+(zz3-zz2)**2) 
		l(3)=sqrt((xx4-xx3)**2+(yy4-yy3)**2+(zz4-zz3)**2) 
		l(4)=sqrt((xx1-xx4)**2+(yy1-yy4)**2+(zz1-zz4)**2) 

		do ii=1,4
			AD=0
			do jj=1,4
				if ((l(ii)-l(jj))>=-1.e-8)then
					AD=AD+1
			end if 
		end do  
			res(ii)=AD
		end do 
	  
		do ii=1,4
			if (res(ii)==4)then
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
!************************************************************************
! Calculate influence coefficients
!************************************************************************        
     IEvect(1,1)=LCS(i,1)
     IEvect(1,2)=LCS(i,2)
	 IEvect(1,3)=LCS(i,3)
            
     IEvect(2,1)=LCS(i,4)
     IEvect(2,2)=LCS(i,5)
     IEvect(2,3)=LCS(i,6)           
            
     IEvect(3,1)=LCS(i,7)
     IEvect(3,2)=LCS(i,8)
     IEvect(3,3)=LCS(i,9)

     call convert (xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,JEvect)   
     call stress_regular_fs3d_qua_num_part3 (xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4,gpx,gpy,gpz,E,v,M,IEvect,JEvect,COEFF1)  
     call          dis_regular_fs3d_qua_num (xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4,gpx,gpy,gpz,E,v,M,IEvect,JEvect,COEFF2)  
        
     PPxx(i)= PPxx(i)+COEFF1(1,1)*X(3*J-2)+COEFF1(1,2)*X(3*J-1)+COEFF1(1,3)*X(3*J)
	 PPxy(i)= PPxy(i)+COEFF1(2,1)*X(3*J-2)+COEFF1(2,2)*X(3*J-1)+COEFF1(2,3)*X(3*J)
	 PPxz(i)= PPxz(i)+COEFF1(3,1)*X(3*J-2)+COEFF1(3,2)*X(3*J-1)+COEFF1(3,3)*X(3*J)                     
	 PPyy(i)= PPyy(i)+COEFF1(4,1)*X(3*J-2)+COEFF1(4,2)*X(3*J-1)+COEFF1(4,3)*X(3*J)
	 PPyz(i)= PPyz(i)+COEFF1(5,1)*X(3*J-2)+COEFF1(5,2)*X(3*J-1)+COEFF1(5,3)*X(3*J)
	 PPzz(i)= PPzz(i)+COEFF1(6,1)*X(3*J-2)+COEFF1(6,2)*X(3*J-1)+COEFF1(6,3)*X(3*J)        
                       
     Ux(i)= Ux(i)+COEFF2(1,1)*X(3*J-2)+COEFF2(1,2)*X(3*J-1)+COEFF2(1,3)*X(3*J)
     Uy(i)= Uy(i)+COEFF2(2,1)*X(3*J-2)+COEFF2(2,2)*X(3*J-1)+COEFF2(2,3)*X(3*J)
     Uz(i)= Uz(i)+COEFF2(3,1)*X(3*J-2)+COEFF2(3,2)*X(3*J-1)+COEFF2(3,3)*X(3*J)                                 
200 continue 
    
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
 
    Trac(1)=Pxx*EI(1,1)+Pxy*EI(2,1)+Pxz*EI(3,1)
    Trac(2)=Pxy*EI(1,1)+Pyy*EI(2,1)+Pyz*EI(3,1)
    Trac(3)=Pxz*EI(1,1)+Pyz*EI(2,1)+Pzz*EI(3,1)
      
    nor(i,1)=Trac(1)*EI(1,1)+Trac(2)*EI(2,1)+Trac(3)*EI(3,1)
    nor(i,2)=Trac(1)*EI(1,2)+Trac(2)*EI(2,2)+Trac(3)*EI(3,2)
    nor(i,3)=Trac(1)*EI(1,3)+Trac(2)*EI(2,3)+Trac(3)*EI(3,3)  
     
    Trac(1)=Pxx*EI(1,2)+Pxy*EI(2,2)+Pxz*EI(3,2)
    Trac(2)=Pxy*EI(1,2)+Pyy*EI(2,2)+Pyz*EI(3,2)
    Trac(3)=Pxz*EI(1,2)+Pyz*EI(2,2)+Pzz*EI(3,2)
     
    nor(i,4)=Trac(1)*EI(1,2)+Trac(2)*EI(2,2)+Trac(3)*EI(3,2)
    nor(i,5)=Trac(1)*EI(1,3)+Trac(2)*EI(2,3)+Trac(3)*EI(3,3)  
     
    Trac(1)=Pxx*EI(1,3)+Pxy*EI(2,3)+Pxz*EI(3,3)
    Trac(2)=Pxy*EI(1,3)+Pyy*EI(2,3)+Pyz*EI(3,3)
    Trac(3)=Pxz*EI(1,3)+Pyz*EI(2,3)+Pzz*EI(3,3)
     
    nor(i,6)=Trac(1)*EI(1,3)+Trac(2)*EI(2,3)+Trac(3)*EI(3,3)   
     
100 continue 
!*******************************************************************
! Superimpose the initial stress field    
!*******************************************************************
do i=1,NOP
    PPxx(i)=PPxx(i)+nor(i,1)
    PPxy(i)=PPxy(i)+nor(i,2)
    PPxz(i)=PPxz(i)+nor(i,3)
    PPyy(i)=PPyy(i)+nor(i,4)
    PPyz(i)=PPyz(i)+nor(i,5)
    PPzz(i)=PPzz(i)+nor(i,6)
end do   

call cpu_time(t2)

write (*,*)'Successful !!!'
write(*,'(A18,f7.2,a2)')' Computing time = ',t2-t1,' s'
!*******************************************************************
!  Write the output file
!*******************************************************************
!open (40, file='dis.txt')
!do i=1,NOP
!	write (40,'(I5,6es12.4)')i,CFP(i,1:3),Ux(i),Uy(i),Uz(i)
!end do 
!close (40)
!
!open (40, file='stress.txt')
!do i=1,NOP
!	write (40,'(I5,9es12.4)')i,CFP(i,1:3),PPxx(i),PPyy(i),PPzz(i),PPxy(i),PPxz(i),PPyz(i)
!end do 
!close (40)

open (41,file='output_part3.txt')
write (41,*)'*************************************************************************'
write (41,*)'************		    FSM3D_qua_part3		        **********'
write (41,*)'************	    developed by Wong and Cui (2021)            **********'
write (41,*)'************	       DDFS3D (version 1.0), 31-5-2021		**********'
write (41,*)'*************************************************************************'
write (41,'(2A)')' Title: ',adjustl(Title)
write (41,'(A9,es12.4,A15,f10.4,A5,es12.4,A5,es12.4)')'Modulus=',E,'Poison ratio=',v,'  KS=',KS,'  KN=',KN
write (41,*)'Initial stress field:'
write (41,'(A5,es12.4,A5,es12.4,A5,es12.4)')'Pxx=',pxx,'Pyy=',pyy,'Pzz=',pzz
write (41,'(A5,es12.4,A5,es12.4,A5,es12.4)')'Pxy=',pxy,'Pxz=',pxz,'Pyz=',pyz
write (41,'(A20,I10)')'Number of elements:',NELE
write (41,'(A20,I10)')'Number of nodes   :',NNOD
write (41,'(A20,f10.2,A2)')'Computing time    :',t2-t1,'s'
write (41,'(A86)')'Stresses and displacements at specific points in the prescribed coordinate system(s):' 
do i=1,NOP
    write (41,'(I5,A6,f15.6,A6,f15.6,A6,f15.6)')i,'x=',CFP(i,1),'y=',CFP(i,2),'z=',CFP(i,3)
	write (41,'(A12,ES15.6,A8,ES15.6,A8,ES15.6)')'         l1=',LCS(i,1),' m1=',LCS(i,2),' n1=',LCS(i,3)
	write (41,'(A12,ES15.6,A8,ES15.6,A8,ES15.6)')'         l2=',LCS(i,4),' m2=',LCS(i,5),' n2=',LCS(i,6)
	write (41,'(A12,ES15.6,A8,ES15.6,A8,ES15.6)')'         l3=',LCS(i,7),' m3=',LCS(i,8),' n3=',LCS(i,9)
    write (41,'(A12,ES15.6,A8,ES15.6,A8,ES15.6)')'        Pxx=',PPxx(i),'Pxy=',PPxy(i),'Pxz=',PPxz(i)
    write (41,'(A12,ES15.6,A8,ES15.6,A8,ES15.6)')'        Pyy=',PPyy(i),'Pyz=',PPyz(i),'Pzz=',PPzz(i)
    write (41,'(A12,ES15.6,A8,ES15.6,A8,ES15.6)')'        Ux=',Ux(i),'Uy=',Uy(i),'Uz=',Uz(i)
end do
close(41)

end program FSM3D_qua_part3 