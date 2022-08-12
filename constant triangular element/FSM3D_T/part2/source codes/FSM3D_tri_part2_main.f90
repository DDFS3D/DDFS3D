!**************************************************************************************
!Part 2: Calculate the stresses and displacements at the centriod of each element.
!**************************************************************************************

program FSM3D_tri_Part2
 
use prog_bar
use subs_fs3d_tri_num
use subs_fs3d_tri_ana
    
implicit none 

real*8::E,v,KN,KS,x1,y1,z1,x2,y2,z2,x3,y3,z3,xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3
real*8::gpx,gpy,gpz,COEFF1(3,3),COEFF2(3,3),IEvect(3,3),JEvect(3,3),pp
integer::Self_T,M,NELE,NNOD,AD,TY,i,j,k,ii,jj,RES(6)
real*8::t1,t2,D1,D2,ratio,L(6),gpxx,gpyy,gpzz,pxx,pyy,pzz,pxy,pxz,pyz,EI(3,3),Evect(3,3),trac(3) 
real*8,allocatable::NOD(:,:),VOB(:,:),nor(:,:),str(:,:),dis(:,:),X(:)
integer,allocatable::ELE(:,:),SOB(:,:),E_type(:)
type (type_prog)::progress
character*100::Title

!*****************************************************************************
write(*,*)'*************************************************************************'
write(*,*)'************		    FSM3D_tri_part2		        **********'
write(*,*)'************	      Algorithm developed by                    **********'
write(*,*)'************	  Kuriyama et al.(1995) and Wong & Cui (2021)   **********'
write(*,*)'************	   Codes developed by Wong & Cui (2021)	        **********'
write(*,*)'************	       DDFS3D (version 1.0), 31-5-2021		**********'
write(*,*)'*************************************************************************'

call cpu_time(t1)

open (21,file='input_part1.txt')
read (21,'(A100)') Title 
read (21,*) E,V,KS,KN,Self_T
read (21,*) pxx,pyy,pzz,pxy,pxz,pyz                
read (21,*) NELE,NNOD

allocate(NOD(NNOD,3),VOB(NELE,3),X(3*NELE),ELE(NELE,3),SOB(NELE,3),E_type(NELE),nor(NELE,3),str(NELE,3),dis(NELE,3))

write (*,*)'Reading element information ...'
do i=1,NELE
   read (21,*) AD,ELE(i,1:3)
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

do i=1,NELE
    str(i,1)=0.
    str(i,2)=0.
    str(i,3)=0.
    dis(i,1)=0.
    dis(i,2)=0.
    dis(i,3)=0.
end do 
!*********************************************************
write (*,*)'Calculating stresses and displacements at the centriod of each element ...'
do 100 i=1,NELE     
	call progress%input(NELE,20)
    call progress%output(i)
	
	x1=NOD(ELE(I,1),1)
	y1=NOD(ELE(I,1),2)
	z1=NOD(ELE(I,1),3)
       
	x2=NOD(ELE(I,2),1)
	y2=NOD(ELE(I,2),2)
	z2=NOD(ELE(I,2),3)
       
	x3=NOD(ELE(I,3),1)
	y3=NOD(ELE(I,3),2)
	z3=NOD(ELE(I,3),3)
	
	gpx=(x1+x2+x3)/3.
	gpy=(y1+y2+y3)/3.
	gpz=(z1+z2+z3)/3.
	
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
           
        gpxx=(xx1+xx2+xx3)/3.
        gpyy=(yy1+yy2+yy3)/3.
        gpzz=(zz1+zz2+zz3)/3.
       
!*****************************************************
! Determine the number of Gaussian points
!*****************************************************
        l(1)=sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2) 
        l(2)=sqrt((x3-x2)**2+(y3-y2)**2+(z3-z2)**2) 
        l(3)=sqrt((x3-x1)**2+(y3-y1)**2+(z3-z1)**2) 
        l(4)=sqrt((xx2-xx1)**2+(yy2-yy1)**2+(zz2-zz1)**2)
        l(5)=sqrt((xx3-xx2)**2+(yy3-yy2)**2+(zz3-zz2)**2) 
        l(6)=sqrt((xx3-xx1)**2+(yy3-yy1)**2+(zz3-zz1)**2) 
 
        pp=(l(1)+l(2)+l(3))/2.

        do ii=1,6
           AD=0
           do jj=1,6
               if ((l(ii)-l(jj))>=-1.e-8)then
                   AD=AD+1
			   end if 
		   end do  
		   res(ii)=AD
		end do 

        do ii=1,6
            if (res(ii)==6)then
                D1=l(ii)
	        end if 
	    end do 

        D2=sqrt((gpx-gpxx)**2+(gpy-gpyy)**2+(gpz-gpzz)**2)
        ratio=D2/D1

		if (ratio<=0.5) then
			M=21
		elseif (0.5<ratio.and.ratio<3)then
			M=ceiling(8+0.1923*ratio)
		else 
			M=8
		end if 

		if (i==j) then
			M=17
		end if 
!************************************************************************
! Calculate the coefficient matirx
!************************************************************************        
        call convert (x1,y1,z1,x2,y2,z2,x3,y3,z3,IEvect)
        call convert (xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,JEvect)
        if (i/=j)then
            call stress_regular_fs3d_tri_num_part1 (xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,gpx,gpy,gpz,E,v,M,IEvect,JEvect,COEFF1)  
            call          dis_regular_fs3d_tri_num (xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,gpx,gpy,gpz,E,v,M,IEvect,JEvect,COEFF2)  
        else 
            if (self_T==0)then
                call diag_stress_fs3d_tri_ana (xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,E,v,IEvect,COEFF1)
                call    diag_dis_fs3d_tri_ana (xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,E,v,IEvect,COEFF2)
			else 
				call stress_singular_fs3d_tri_num (xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,E,v,M,IEvect,COEFF1)
                call    dis_singular_fs3d_tri_num (xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,E,v,M,IEvect,COEFF2)
			end if 
        end if 

		Str(i,1)= Str(i,1)+COEFF1(1,1)*X(3*J-2)+COEFF1(1,2)*X(3*J-1)+COEFF1(1,3)*X(3*J)
        Str(i,2)= Str(i,2)+COEFF1(2,1)*X(3*J-2)+COEFF1(2,2)*X(3*J-1)+COEFF1(2,3)*X(3*J)
        Str(i,3)= Str(i,3)+COEFF1(3,1)*X(3*J-2)+COEFF1(3,2)*X(3*J-1)+COEFF1(3,3)*X(3*J)
                          
        Dis(i,1)=Dis(i,1)+COEFF2(1,1)*X(3*J-2)+COEFF2(1,2)*X(3*J-1)+COEFF2(1,3)*X(3*J)
        Dis(i,2)=Dis(i,2)+COEFF2(2,1)*X(3*J-2)+COEFF2(2,2)*X(3*J-1)+COEFF2(2,3)*X(3*J)
        Dis(i,3)=Dis(i,3)+COEFF2(3,1)*X(3*J-2)+COEFF2(3,2)*X(3*J-1)+COEFF2(3,3)*X(3*J)
     
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
        
	Trac(1)=Pxx*EI(1,3)+Pxy*EI(2,3)+Pxz*EI(3,3)
	Trac(2)=Pxy*EI(1,3)+Pyy*EI(2,3)+Pyz*EI(3,3)
	Trac(3)=Pxz*EI(1,3)+Pyz*EI(2,3)+Pzz*EI(3,3)
     
	nor(i,1)=Trac(1)*EI(1,1)+Trac(2)*EI(2,1)+Trac(3)*EI(3,1)
	nor(i,2)=Trac(1)*EI(1,2)+Trac(2)*EI(2,2)+Trac(3)*EI(3,2)
	nor(i,3)=Trac(1)*EI(1,3)+Trac(2)*EI(2,3)+Trac(3)*EI(3,3)   
     
100 continue 

!*******************************************************************
! Superimpose the initial stress field     
!*******************************************************************
do i=1,NELE
    str(i,1)=str(i,1)+nor(i,1)
    str(i,2)=str(i,2)+nor(i,2)
    str(i,3)=str(i,3)+nor(i,3)
end do   

call cpu_time(t2)

write (*,*)'Successful !!!'
write (*,'(A18,f7.2,a2)')' Computing time = ',t2-t1,' s'
!*******************************************************************
!  Write the output file 
!*******************************************************************
open (41,file='output_part2.txt')
write (41,*)'*************************************************************************'
write (41,*)'************		    FSM3D_tri_part2		        **********'
write (41,*)'************	      Algorithm developed by                    **********'
write (41,*)'************	  Kuriyama et al.(1995) and Wong & Cui (2021)   **********'
write (41,*)'************	   Codes developed by Wong & Cui (2021)	        **********'
write (41,*)'************	       DDFS3D (version 1.0), 31-5-2021		**********'
write (41,*)'*************************************************************************'
write (41,'(2A)')' Title: ', adjustl(Title)
write (41,'(A9,es12.4,A15,f10.4,A5,es12.4,A5,es12.4)')'Modulus=',E,'Poison ratio=',v,'  KS=',KS,'  KN=',KN
write (41,*)'Initial stress field:'
write (41,'(A5,es12.4,A5,es12.4,A5,es12.4)')'Pxx=',pxx,'Pyy=',pyy,'Pzz=',pzz
write (41,'(A5,es12.4,A5,es12.4,A5,es12.4)')'Pxy=',pxy,'Pxz=',pxz,'Pyz=',pyz
write (41,'(A20,I10)')'Number of elements:',NELE
write (41,'(A20,I10)')'Number of nodes   :',NNOD
write (41,'(A20,f10.2,A2)')'Computing time    :',t2-t1,'s'
write (41,'(A83)')'Stresses and displacments at the controid of each element in the respective local coordinate system:' 
write (41,'(A83)')' No.	    Szx            Szy          Szz            Ux            Uy           Uz '
do i=1,NELE
    write(41,'(I5,6es14.4)')i,str(i,1:3),dis(i,1:3)
end do
close(41)

end program FSM3D_tri_Part2 