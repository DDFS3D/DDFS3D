module subs_fs3d_qua_num
    
contains 

!**************************************************************

Subroutine convert (x1,y1,z1,x2,y2,z2,x3,y3,z3,Evect)
implicit none 

real*8::x1,y1,z1,x2,y2,z2,x3,y3,z3,Evect(3,3)
real*8::VecA_x,VecA_y,VecA_z,VecB_x,VecB_y,VecB_z,VecC_x,VecC_y,VecC_z,PX,PY,PZ,cx,cy,cz
real*8::lengthA,lengthC,AZ,AY,AX

    VecA_x=x2-x1
    VecA_y=y2-y1
    VecA_z=z2-z1
    
    VecB_x=x3-x1
    VecB_y=y3-y1
    VecB_z=z3-z1
    
    VecC_x=VecA_y*VecB_z-VecA_z*VecB_y
    VecC_y=VecA_z*VecB_x-VecA_x*VecB_z
    VecC_z=VecA_x*VecB_y-VecA_y*VecB_x
    
    lengthA=sqrt(VecA_x**2+VecA_y**2+VecA_z**2)
    lengthC=sqrt(VecC_x**2+VecC_y**2+VecC_z**2)
    
    Ax=VecA_x/lengthA
    Ay=VecA_y/lengthA
    Az=VecA_z/lengthA
    
    Cx=VecC_x/lengthC
    Cy=VecC_y/lengthC
    Cz=VecC_z/lengthC
    
    Px=Cy*Az-Cz*Ay
    Py=Cz*Ax-Cx*Az
    Pz=Cx*Ay-Cy*Ax
    
    Evect(1,1)=Ax
    Evect(1,2)=Ay
    Evect(1,3)=Az
    
    Evect(2,1)=Px
    Evect(2,2)=Py
    Evect(2,3)=Pz
    
    Evect(3,1)=Cx
    Evect(3,2)=Cy
    Evect(3,3)=Cz
    
end subroutine convert


subroutine dis_regular_fs3d_qua_num (xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4,gpx,gpy,gpz,E,v,M,IEvect,JEvect,COEFF)
implicit none 

real*8::E,v,cond,l1,l2,l3,p,area,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4,gpx,gpy,gpz
real*8::VecIJ_x,VecIJ_y,VecIJ_z,localx,localy,localz,IEvect(3,3),JEvect(3,3),COEFF(3,3)
real*8::f1,f2,f3,f4,f5,f6,f7,G,kesi,yita,A,B,U(3,3),JI(3,3)
real*8::N1,N2,N3,N4,jaco_11,jaco_12,jaco_21,jaco_22,jaco,gpxx,gpyy,gpzz,vax,vay,vaz,vbx,vby,vbz,vcx,vcy,vcz,vdx,vdy,vdz
integer::i,j,k,M
real*8,allocatable::pw(:,:)

    G=E/2/(1+v)
    cond=1./16./3.141592653/(1-v)/G
    
    gpxx=(xx1+xx2+xx3+xx4)/4.
    gpyy=(yy1+yy2+yy3+yy4)/4.
    gpzz=(zz1+zz2+zz3+zz4)/4.
    
    VAX=xx1-gpxx
    VAY=yy1-gpyy
    VAZ=zz1-gpzz

    x1=VAX*JEvect(1,1)+VAY*JEvect(1,2)+VAZ*JEvect(1,3)
    y1=VAX*JEvect(2,1)+VAY*JEvect(2,2)+VAZ*JEvect(2,3)
   
    VBX=xx2-gpxx
    VBY=yy2-gpyy
    VBZ=zz2-gpzz

    x2=VBX*JEvect(1,1)+VBY*JEvect(1,2)+VBZ*JEvect(1,3)
    y2=VBX*JEvect(2,1)+VBY*JEvect(2,2)+VBZ*JEvect(2,3) 
    
    VCX=xx3-gpxx
    VCY=yy3-gpyy
    VCZ=zz3-gpzz

    x3=VCX*JEvect(1,1)+VCY*JEvect(1,2)+VCZ*JEvect(1,3)
    y3=VCX*JEvect(2,1)+VCY*JEvect(2,2)+VCZ*JEvect(2,3)
    
    VDX=xx4-gpxx
    VDY=yy4-gpyy
    VDZ=zz4-gpzz

    x4=VDX*JEvect(1,1)+VDY*JEvect(1,2)+VDZ*JEvect(1,3)
    y4=VDX*JEvect(2,1)+VDY*JEvect(2,2)+VDZ*JEvect(2,3)
    
    VecIJ_x=gpx-gpxx
    VecIJ_y=gpy-gpyy
    VecIJ_z=gpz-gpzz
    
    Localx=VecIJ_x*JEvect(1,1)+VecIJ_y*JEvect(1,2)+VecIJ_z*JEvect(1,3)
    Localy=VecIJ_x*JEvect(2,1)+VecIJ_y*JEvect(2,2)+VecIJ_z*JEvect(2,3)
    Localz=VecIJ_x*JEvect(3,1)+VecIJ_y*JEvect(3,2)+VecIJ_z*JEvect(3,3)
    
!**************************************************************  
allocate (pw(M,2))
call Gauss_ponits_weights(M,pw)
!*****************************************************************************
!Calculate f1   
    f1=0.
    do i=1,M
        do j=1,M 
             N1=(1+pw(i,1))*(1+pw(j,1))/4.
             N2=(1-pw(i,1))*(1+pw(j,1))/4.
             N3=(1-pw(i,1))*(1-pw(j,1))/4.
             N4=(1+pw(i,1))*(1-pw(j,1))/4.
             JACO_11=(1+pw(j,1))/4.*(x1-x2)+(1-pw(j,1))/4.*(x4-x3)
             JACO_12=(1+pw(i,1))/4.*(x1-x4)+(1-pw(i,1))/4.*(x2-x3)
             JACO_21=(1+pw(j,1))/4.*(y1-y2)+(1-pw(j,1))/4.*(y4-y3)
             JACO_22=(1+pw(i,1))/4.*(y1-y4)+(1-pw(i,1))/4.*(y2-y3)
             JACO=JACO_11*JACO_22-JACO_12*JACO_21
             
             kesi=x1*N1+x2*N2+x3*N3+x4*N4
             yita=y1*N1+y2*N2+y3*N3+y4*N4
            
             A=(localx-kesi)**2+(localy-yita)**2+localz**2
             B=pw(i,2)*pw(j,2)*JACO
             
             f1=f1+B*A**(-0.5)               
        end do 
    end do   

!Calculate f2
    f2=0.
    do i=1,M
        do j=1,M
            
             N1=(1+pw(i,1))*(1+pw(j,1))/4.
             N2=(1-pw(i,1))*(1+pw(j,1))/4.
             N3=(1-pw(i,1))*(1-pw(j,1))/4.
             N4=(1+pw(i,1))*(1-pw(j,1))/4.
             JACO_11=(1+pw(j,1))/4.*(x1-x2)+(1-pw(j,1))/4.*(x4-x3)
             JACO_12=(1+pw(i,1))/4.*(x1-x4)+(1-pw(i,1))/4.*(x2-x3)
             JACO_21=(1+pw(j,1))/4.*(y1-y2)+(1-pw(j,1))/4.*(y4-y3)
             JACO_22=(1+pw(i,1))/4.*(y1-y4)+(1-pw(i,1))/4.*(y2-y3)
             JACO=JACO_11*JACO_22-JACO_12*JACO_21
             
             kesi=x1*N1+x2*N2+x3*N3+x4*N4
             yita=y1*N1+y2*N2+y3*N3+y4*N4
            
             A=(localx-kesi)**2+(localy-yita)**2+localz**2
             B=pw(i,2)*pw(j,2)*JACO
            
             f2=f2+B*A**(-1.5)*localz
        end do 
    end do

!Calculate f3 
    f3=0.
    do i=1,M
        do j=1,M
             N1=(1+pw(i,1))*(1+pw(j,1))/4.
             N2=(1-pw(i,1))*(1+pw(j,1))/4.
             N3=(1-pw(i,1))*(1-pw(j,1))/4.
             N4=(1+pw(i,1))*(1-pw(j,1))/4.
             JACO_11=(1+pw(j,1))/4.*(x1-x2)+(1-pw(j,1))/4.*(x4-x3)
             JACO_12=(1+pw(i,1))/4.*(x1-x4)+(1-pw(i,1))/4.*(x2-x3)
             JACO_21=(1+pw(j,1))/4.*(y1-y2)+(1-pw(j,1))/4.*(y4-y3)
             JACO_22=(1+pw(i,1))/4.*(y1-y4)+(1-pw(i,1))/4.*(y2-y3)
             JACO=JACO_11*JACO_22-JACO_12*JACO_21
             
             kesi=x1*N1+x2*N2+x3*N3+x4*N4
             yita=y1*N1+y2*N2+y3*N3+y4*N4
            
             A=(localx-kesi)**2+(localy-yita)**2+localz**2
             B=pw(i,2)*pw(j,2)*JACO
             
             f3=f3+B*A**(-1.5)*(localx-kesi)
        end do 
    end do
  
!Calculate f4 
    f4=0. 
    do i=1,M
        do j=1,M
             N1=(1+pw(i,1))*(1+pw(j,1))/4.
             N2=(1-pw(i,1))*(1+pw(j,1))/4.
             N3=(1-pw(i,1))*(1-pw(j,1))/4.
             N4=(1+pw(i,1))*(1-pw(j,1))/4.
             JACO_11=(1+pw(j,1))/4.*(x1-x2)+(1-pw(j,1))/4.*(x4-x3)
             JACO_12=(1+pw(i,1))/4.*(x1-x4)+(1-pw(i,1))/4.*(x2-x3)
             JACO_21=(1+pw(j,1))/4.*(y1-y2)+(1-pw(j,1))/4.*(y4-y3)
             JACO_22=(1+pw(i,1))/4.*(y1-y4)+(1-pw(i,1))/4.*(y2-y3)
             JACO=JACO_11*JACO_22-JACO_12*JACO_21
             
             kesi=x1*N1+x2*N2+x3*N3+x4*N4
             yita=y1*N1+y2*N2+y3*N3+y4*N4
            
             A=(localx-kesi)**2+(localy-yita)**2+localz**2
             B=pw(i,2)*pw(j,2)*JACO
             
             f4=f4+B*A**(-1.5)*(localy-yita)
        end do 
    end do      

!Calculate f5 
    f5=0. 
    do i=1,M
        do j=1,M
             N1=(1+pw(i,1))*(1+pw(j,1))/4.
             N2=(1-pw(i,1))*(1+pw(j,1))/4.
             N3=(1-pw(i,1))*(1-pw(j,1))/4.
             N4=(1+pw(i,1))*(1-pw(j,1))/4.
             JACO_11=(1+pw(j,1))/4.*(x1-x2)+(1-pw(j,1))/4.*(x4-x3)
             JACO_12=(1+pw(i,1))/4.*(x1-x4)+(1-pw(i,1))/4.*(x2-x3)
             JACO_21=(1+pw(j,1))/4.*(y1-y2)+(1-pw(j,1))/4.*(y4-y3)
             JACO_22=(1+pw(i,1))/4.*(y1-y4)+(1-pw(i,1))/4.*(y2-y3)
             JACO=JACO_11*JACO_22-JACO_12*JACO_21
             
             kesi=x1*N1+x2*N2+x3*N3+x4*N4
             yita=y1*N1+y2*N2+y3*N3+y4*N4
            
             A=(localx-kesi)**2+(localy-yita)**2+localz**2
             B=pw(i,2)*pw(j,2)*JACO
             
             f5=f5+B*A**(-1.5)*(localx-kesi)**2
        end do 
    end do
  
!Calculate f6 
    f6=0. 
    do i=1,M
        do j=1,M
             N1=(1+pw(i,1))*(1+pw(j,1))/4.
             N2=(1-pw(i,1))*(1+pw(j,1))/4.
             N3=(1-pw(i,1))*(1-pw(j,1))/4.
             N4=(1+pw(i,1))*(1-pw(j,1))/4.
             JACO_11=(1+pw(j,1))/4.*(x1-x2)+(1-pw(j,1))/4.*(x4-x3)
             JACO_12=(1+pw(i,1))/4.*(x1-x4)+(1-pw(i,1))/4.*(x2-x3)
             JACO_21=(1+pw(j,1))/4.*(y1-y2)+(1-pw(j,1))/4.*(y4-y3)
             JACO_22=(1+pw(i,1))/4.*(y1-y4)+(1-pw(i,1))/4.*(y2-y3)
             JACO=JACO_11*JACO_22-JACO_12*JACO_21
             
             kesi=x1*N1+x2*N2+x3*N3+x4*N4
             yita=y1*N1+y2*N2+y3*N3+y4*N4
            
             A=(localx-kesi)**2+(localy-yita)**2+localz**2
             B=pw(i,2)*pw(j,2)*JACO
             
             f6=f6+B*A**(-1.5)*(localy-yita)**2
        end do 
    end do  

!Calculate f7  
    f7=0. 
    do i=1,M
        do j=1,M
             N1=(1+pw(i,1))*(1+pw(j,1))/4.
             N2=(1-pw(i,1))*(1+pw(j,1))/4.
             N3=(1-pw(i,1))*(1-pw(j,1))/4.
             N4=(1+pw(i,1))*(1-pw(j,1))/4.
             JACO_11=(1+pw(j,1))/4.*(x1-x2)+(1-pw(j,1))/4.*(x4-x3)
             JACO_12=(1+pw(i,1))/4.*(x1-x4)+(1-pw(i,1))/4.*(x2-x3)
             JACO_21=(1+pw(j,1))/4.*(y1-y2)+(1-pw(j,1))/4.*(y4-y3)
             JACO_22=(1+pw(i,1))/4.*(y1-y4)+(1-pw(i,1))/4.*(y2-y3)
             JACO=JACO_11*JACO_22-JACO_12*JACO_21
             
             kesi=x1*N1+x2*N2+x3*N3+x4*N4
             yita=y1*N1+y2*N2+y3*N3+y4*N4
            
             A=(localx-kesi)**2+(localy-yita)**2+localz**2
             B=pw(i,2)*pw(j,2)*JACO
             
             f7=f7+B*A**(-1.5)*(localy-yita)*(localx-kesi)
        end do 
    end do    
        
!!Calculate f8 
!    f8=0. 
!    do i=1,M
!        JJ=M+2-i
!        allocate (pw2(JJ,2))
!        call Gauss_ponits_weights(JJ,pw2)
!        do j=1,JJ
!             kesi=x1+0.5*(x2-x1)*(1+pw1(i,1))+0.25*(x3-x1)*(1-pw1(i,1))*(1+pw2(j,1))
!             yita=y1+0.5*(y2-y1)*(1+pw1(i,1))+0.25*(y3-y1)*(1-pw1(i,1))*(1+pw2(j,1))
!             A=(localx-kesi)**2+(localy-yita)**2+localz**2
!             B=2*Area*pw1(i,2)*pw2(j,2)*(1-pw1(i,1))/8.
!             
!             f8=f8+B*A**(-2.5)*localz**3
!             
!        end do 
!             deallocate(pw2)
!    end do   
!   
!!Calculate f9    
!    f9=0. 
!    do i=1,M
!        JJ=M+2-i
!        allocate (pw2(JJ,2))
!        call Gauss_ponits_weights(JJ,pw2)
!        do j=1,JJ
!             kesi=x1+0.5*(x2-x1)*(1+pw1(i,1))+0.25*(x3-x1)*(1-pw1(i,1))*(1+pw2(j,1))
!             yita=y1+0.5*(y2-y1)*(1+pw1(i,1))+0.25*(y3-y1)*(1-pw1(i,1))*(1+pw2(j,1))
!             A=(localx-kesi)**2+(localy-yita)**2+localz**2
!             B=2*Area*pw1(i,2)*pw2(j,2)*(1-pw1(i,1))/8.
!             
!             f9=f9+B*A**(-2.5)*(localx-kesi)*localz**2
!      
!        end do    
!            deallocate(pw2)
!    end do 
!*****************************************************************
! Calculate influence coefficients   
!*****************************************************************   
      U(1,1)= F5          + ( 3.0D0 - 4.0D0 * V ) * F1
      U(1,2)= F7
      U(1,3)= F3 * localZ
      U(2,1)= F7
      U(2,2)= F6          + ( 3.0D0 - 4.0D0 * V ) * F1
      U(2,3)= F4 * localZ
      U(3,1)= F3 * localZ
      U(3,2)= F4 * localZ
      U(3,3)= F2 * localZ + ( 3.0D0 - 4.0D0 * V ) * F1
        
      JI(1,1)=JEvect(1,1)*IEvect(1,1)+JEvect(1,2)*IEvect(1,2)+JEvect(1,3)*IEvect(1,3)
      JI(2,1)=JEvect(2,1)*IEvect(1,1)+JEvect(2,2)*IEvect(1,2)+JEvect(2,3)*IEvect(1,3)
      JI(3,1)=JEvect(3,1)*IEvect(1,1)+JEvect(3,2)*IEvect(1,2)+JEvect(3,3)*IEvect(1,3)
      
      JI(1,2)=JEvect(1,1)*IEvect(2,1)+JEvect(1,2)*IEvect(2,2)+JEvect(1,3)*IEvect(2,3)
      JI(2,2)=JEvect(2,1)*IEvect(2,1)+JEvect(2,2)*IEvect(2,2)+JEvect(2,3)*IEvect(2,3)
      JI(3,2)=JEvect(3,1)*IEvect(2,1)+JEvect(3,2)*IEvect(2,2)+JEvect(3,3)*IEvect(2,3) 
      
      JI(1,3)=JEvect(1,1)*IEvect(3,1)+JEvect(1,2)*IEvect(3,2)+JEvect(1,3)*IEvect(3,3)
      JI(2,3)=JEvect(2,1)*IEvect(3,1)+JEvect(2,2)*IEvect(3,2)+JEvect(2,3)*IEvect(3,3)
      JI(3,3)=JEvect(3,1)*IEvect(3,1)+JEvect(3,2)*IEvect(3,2)+JEvect(3,3)*IEvect(3,3)  
      
      COEFF(1,1)=cond*(U(1,1)*JI(1,1)+U(2,1)*JI(2,1)+U(3,1)*JI(3,1))
      COEFF(2,1)=cond*(U(1,1)*JI(1,2)+U(2,1)*JI(2,2)+U(3,1)*JI(3,2))
      COEFF(3,1)=cond*(U(1,1)*JI(1,3)+U(2,1)*JI(2,3)+U(3,1)*JI(3,3))
     
      COEFF(1,2)=cond*(U(1,2)*JI(1,1)+U(2,2)*JI(2,1)+U(3,2)*JI(3,1))
      COEFF(2,2)=cond*(U(1,2)*JI(1,2)+U(2,2)*JI(2,2)+U(3,2)*JI(3,2))
      COEFF(3,2)=cond*(U(1,2)*JI(1,3)+U(2,2)*JI(2,3)+U(3,2)*JI(3,3))
     
      COEFF(1,3)=cond*(U(1,3)*JI(1,1)+U(2,3)*JI(2,1)+U(3,3)*JI(3,1))
      COEFF(2,3)=cond*(U(1,3)*JI(1,2)+U(2,3)*JI(2,2)+U(3,3)*JI(3,2))
      COEFF(3,3)=cond*(U(1,3)*JI(1,3)+U(2,3)*JI(2,3)+U(3,3)*JI(3,3))
	      
end subroutine dis_regular_fs3d_qua_num
     

subroutine dis_singular_fs3d_qua_num (xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4,EE,v,M,IEvect,COEFF)
implicit none

real*8::EE,v,GG,cons,l1,l2,l3,p,area,x(5),y(5),z(5),xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4,xx(4),yy(4),zz(4)
real*8::IEvect(3,3),JEvect(3,3),COEFF(3,3),A,B,C,D,E,F,yip,xm,ym,r,localz,gpx,gpy,gpz
real*8::f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,cond,U(3,3),JI(3,3)
integer::M,I,J,ii,jj
real*8,allocatable::pw1(:,:)

GG=EE/2/(1+v)
cond=1./16./3.141592653/(1-v)/GG
localz=0.
    
allocate (pw1(M,2))

do i=1,3
	do j=1,3
		JEvect(i,j)=IEvect(i,j)
	end do 
end do 
    
gpx=(xx1+xx2+xx3+xx4)/4.
gpy=(yy1+yy2+yy3+yy4)/4.
gpz=(zz1+zz2+zz3+zz4)/4.

xx(1)=xx1-gpx
yy(1)=yy1-gpy
zz(1)=zz1-gpz

xx(2)=xx2-gpx
yy(2)=yy2-gpy
zz(2)=zz2-gpz

xx(3)=xx3-gpx
yy(3)=yy3-gpy
zz(3)=zz3-gpz

xx(4)=xx4-gpx
yy(4)=yy4-gpy
zz(4)=zz4-gpz

x(1)=xx(1)*IEvect(1,1)+yy(1)*IEvect(1,2)+zz(1)*IEvect(1,3)
y(1)=xx(1)*IEvect(2,1)+yy(1)*IEvect(2,2)+zz(1)*IEvect(2,3)
z(1)=xx(1)*IEvect(3,1)+yy(1)*IEvect(3,2)+zz(1)*IEvect(3,3)

x(2)=xx(2)*IEvect(1,1)+yy(2)*IEvect(1,2)+zz(2)*IEvect(1,3)
y(2)=xx(2)*IEvect(2,1)+yy(2)*IEvect(2,2)+zz(2)*IEvect(2,3)
z(2)=xx(2)*IEvect(3,1)+yy(2)*IEvect(3,2)+zz(2)*IEvect(3,3)

x(3)=xx(3)*IEvect(1,1)+yy(3)*IEvect(1,2)+zz(3)*IEvect(1,3)
y(3)=xx(3)*IEvect(2,1)+yy(3)*IEvect(2,2)+zz(3)*IEvect(2,3)
z(3)=xx(3)*IEvect(3,1)+yy(3)*IEvect(3,2)+zz(3)*IEvect(3,3)

x(4)=xx(4)*IEvect(1,1)+yy(4)*IEvect(1,2)+zz(4)*IEvect(1,3)
y(4)=xx(4)*IEvect(2,1)+yy(4)*IEvect(2,2)+zz(4)*IEvect(2,3)
z(4)=xx(4)*IEvect(3,1)+yy(4)*IEvect(3,2)+zz(4)*IEvect(3,3)

x(5)=x(1)
y(5)=y(1)
z(5)=z(1)

call Gauss_ponits_weights(m,pw1)

f1=0.
f2=0.
f3=0.
f4=0.
f5=0.
f6=0.
f7=0.
f8=0.
f9=0.
f10=0.
f11=0.
f12=0.
f13=0.
f14=0.
f15=0.
f16=0.
f17=0.

do 100 ii=1,4
    
    l1=sqrt(x(ii)**2+y(ii)**2+z(ii)**2)
    l2=sqrt((x(ii+1)-x(ii))**2+(y(ii+1)-y(ii))**2+(z(ii+1)-z(ii))**2)
    l3=sqrt(x(ii+1)**2+y(ii+1)**2+z(ii+1)**2)
    p=(l1+l2+l3)/2.
    area=sqrt(p*(p-l1)*(p-l2)*(p-l3))
       
    A=0.
    B=0.
    C=0.
    D=0.
    E=0.
    F=0.
    !G=0.
    !H=0.
    !I=0.
    !J=0.

    do jj=1,m
        
        yip=0.5*(pw1(jj,1)+1) 
        xm=(1-yip)*x(ii)+yip*x(ii+1)
        ym=(1-yip)*y(ii)+yip*y(ii+1)
        r=sqrt(xm**2+ym**2) 
 !**************************************************************************************        
 ! For self-influence coefficients, only sencond-order partial derivatives are needed.  
 !**************************************************************************************
        A=A+pw1(jj,2)*area/r
        B=B-pw1(jj,2)*area*xm*log(r)/r**3
        C=C-pw1(jj,2)*area*ym*log(r)/r**3
        D=D+pw1(jj,2)*area*xm**2/r**3
        E=E+pw1(jj,2)*area*ym**2/r**3
        F=F+pw1(jj,2)*area*xm*ym/r**3
        !G=G-pw1(jj,2)*area*xm**3*log(r)/r**5
        !H=H-pw1(jj,2)*area*ym**3*log(r)/r**5
        !I=I-pw1(jj,2)*area*ym**2*xm*log(r)/r**5
        !J=J-pw1(jj,2)*area*xm**2*ym*log(r)/r**5
        !K=K+pw1(jj,2)*area*xm*ym/r**5
        !write(*,*)jj,'A=',A,'B=',B,'C=',C,'G=',G
    end do
    
       f1=f1+A
       f2=0
       f3=f3+B
       f4=f4+C
       f5=f5+D
       f6=f6+E
       f7=f7+F
       !f8=0.
       !f9=0.
       !f10=0.
       !f11=0.
       !f12=0.
       !f13=0.
       !f14=f14+G
       !f15=f15+H
       !f16=f16+I
       !f17=f17+J
 
100 continue  
     
!*****************************************************************
! Calculate influence coefficients   
!*****************************************************************   
      U(1,1)= F5          + ( 3.0D0 - 4.0D0 * V ) * F1
      U(1,2)= F7
      U(1,3)= F3 * localZ
      U(2,1)= F7
      U(2,2)= F6          + ( 3.0D0 - 4.0D0 * V ) * F1
      U(2,3)= F4 * localZ
      U(3,1)= F3 * localZ
      U(3,2)= F4 * localZ
      U(3,3)= F2 * localZ + ( 3.0D0 - 4.0D0 * V ) * F1
       
      JI(1,1)=JEvect(1,1)*IEvect(1,1)+JEvect(1,2)*IEvect(1,2)+JEvect(1,3)*IEvect(1,3)
      JI(2,1)=JEvect(2,1)*IEvect(1,1)+JEvect(2,2)*IEvect(1,2)+JEvect(2,3)*IEvect(1,3)
      JI(3,1)=JEvect(3,1)*IEvect(1,1)+JEvect(3,2)*IEvect(1,2)+JEvect(3,3)*IEvect(1,3)
      
      JI(1,2)=JEvect(1,1)*IEvect(2,1)+JEvect(1,2)*IEvect(2,2)+JEvect(1,3)*IEvect(2,3)
      JI(2,2)=JEvect(2,1)*IEvect(2,1)+JEvect(2,2)*IEvect(2,2)+JEvect(2,3)*IEvect(2,3)
      JI(3,2)=JEvect(3,1)*IEvect(2,1)+JEvect(3,2)*IEvect(2,2)+JEvect(3,3)*IEvect(2,3) 
      
      JI(1,3)=JEvect(1,1)*IEvect(3,1)+JEvect(1,2)*IEvect(3,2)+JEvect(1,3)*IEvect(3,3)
      JI(2,3)=JEvect(2,1)*IEvect(3,1)+JEvect(2,2)*IEvect(3,2)+JEvect(2,3)*IEvect(3,3)
      JI(3,3)=JEvect(3,1)*IEvect(3,1)+JEvect(3,2)*IEvect(3,2)+JEvect(3,3)*IEvect(3,3)   
     
      COEFF(1,1)=cond*(U(1,1)*JI(1,1)+U(2,1)*JI(2,1)+U(3,1)*JI(3,1))
      COEFF(2,1)=cond*(U(1,1)*JI(1,2)+U(2,1)*JI(2,2)+U(3,1)*JI(3,2))
      COEFF(3,1)=cond*(U(1,1)*JI(1,3)+U(2,1)*JI(2,3)+U(3,1)*JI(3,3))
     
      COEFF(1,2)=cond*(U(1,2)*JI(1,1)+U(2,2)*JI(2,1)+U(3,2)*JI(3,1))
      COEFF(2,2)=cond*(U(1,2)*JI(1,2)+U(2,2)*JI(2,2)+U(3,2)*JI(3,2))
      COEFF(3,2)=cond*(U(1,2)*JI(1,3)+U(2,2)*JI(2,3)+U(3,2)*JI(3,3))
     
      COEFF(1,3)=cond*(U(1,3)*JI(1,1)+U(2,3)*JI(2,1)+U(3,3)*JI(3,1))
      COEFF(2,3)=cond*(U(1,3)*JI(1,2)+U(2,3)*JI(2,2)+U(3,3)*JI(3,2))
      COEFF(3,3)=cond*(U(1,3)*JI(1,3)+U(2,3)*JI(2,3)+U(3,3)*JI(3,3))
            
end subroutine dis_singular_fs3d_qua_num


subroutine stress_regular_fs3d_qua_num_part1 (xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4,gpx,gpy,gpz,E,v,M,IEvect,JEvect,COEFF)
implicit none 

real*8::E,v,G,cons,l1,l2,l3,p,area,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4,gpx,gpy,gpz
real*8::VecIJ_x,VecIJ_y,VecIJ_z,localx,localy,localz,IEvect(3,3),JEvect(3,3),COEFF(3,3)
real*8::N1,N2,N3,N4,jaco_11,jaco_12,jaco_21,jaco_22,jaco,gpxx,gpyy,gpzz,vax,vay,vaz,vbx,vby,vbz,vcx,vcy,vcz,vdx,vdy,vdz
real*8::f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,kesi,yita,A,B,S(6,3),JI(3,3),TRAC(3)
integer::i,j,k,M

real*8,allocatable::pw(:,:)

    G=E/2./(1+v)
    cons=1./8./3.141592653/(1-v)
    
    gpxx=(xx1+xx2+xx3+xx4)/4.
    gpyy=(yy1+yy2+yy3+yy4)/4.
    gpzz=(zz1+zz2+zz3+zz4)/4.
    
    VAX=xx1-gpxx
    VAY=yy1-gpyy
    VAZ=zz1-gpzz

    x1=VAX*JEvect(1,1)+VAY*JEvect(1,2)+VAZ*JEvect(1,3)
    y1=VAX*JEvect(2,1)+VAY*JEvect(2,2)+VAZ*JEvect(2,3)
   
    VBX=xx2-gpxx
    VBY=yy2-gpyy
    VBZ=zz2-gpzz

    x2=VBX*JEvect(1,1)+VBY*JEvect(1,2)+VBZ*JEvect(1,3)
    y2=VBX*JEvect(2,1)+VBY*JEvect(2,2)+VBZ*JEvect(2,3) 
    
    VCX=xx3-gpxx
    VCY=yy3-gpyy
    VCZ=zz3-gpzz

    x3=VCX*JEvect(1,1)+VCY*JEvect(1,2)+VCZ*JEvect(1,3)
    y3=VCX*JEvect(2,1)+VCY*JEvect(2,2)+VCZ*JEvect(2,3)
    
    VDX=xx4-gpxx
    VDY=yy4-gpyy
    VDZ=zz4-gpzz

    x4=VDX*JEvect(1,1)+VDY*JEvect(1,2)+VDZ*JEvect(1,3)
    y4=VDX*JEvect(2,1)+VDY*JEvect(2,2)+VDZ*JEvect(2,3)
    
    VecIJ_x=gpx-gpxx
    VecIJ_y=gpy-gpyy
    VecIJ_z=gpz-gpzz
    
    Localx=VecIJ_x*JEvect(1,1)+VecIJ_y*JEvect(1,2)+VecIJ_z*JEvect(1,3)
    Localy=VecIJ_x*JEvect(2,1)+VecIJ_y*JEvect(2,2)+VecIJ_z*JEvect(2,3)
    Localz=VecIJ_x*JEvect(3,1)+VecIJ_y*JEvect(3,2)+VecIJ_z*JEvect(3,3)

!**************************************************************
allocate (pw(M,2))
call Gauss_ponits_weights(M,pw)
!*****************************************************************************
!Calculate f1   
    f1=0.
    do i=1,M
        do j=1,M
             N1=(1+pw(i,1))*(1+pw(j,1))/4.
             N2=(1-pw(i,1))*(1+pw(j,1))/4.
             N3=(1-pw(i,1))*(1-pw(j,1))/4.
             N4=(1+pw(i,1))*(1-pw(j,1))/4.
             JACO_11=(1+pw(j,1))/4.*(x1-x2)+(1-pw(j,1))/4.*(x4-x3)
             JACO_12=(1+pw(i,1))/4.*(x1-x4)+(1-pw(i,1))/4.*(x2-x3)
             JACO_21=(1+pw(j,1))/4.*(y1-y2)+(1-pw(j,1))/4.*(y4-y3)
             JACO_22=(1+pw(i,1))/4.*(y1-y4)+(1-pw(i,1))/4.*(y2-y3)
             JACO=JACO_11*JACO_22-JACO_12*JACO_21
             
             kesi=x1*N1+x2*N2+x3*N3+x4*N4
             yita=y1*N1+y2*N2+y3*N3+y4*N4
            
             A=(localx-kesi)**2+(localy-yita)**2+localz**2
             B=pw(i,2)*pw(j,2)*JACO
             
             f1=f1+B*A**(-0.5)          
        end do 
    end do   

!Calculate f2
    f2=0.
    do i=1,M
        do j=1,M 
             N1=(1+pw(i,1))*(1+pw(j,1))/4.
             N2=(1-pw(i,1))*(1+pw(j,1))/4.
             N3=(1-pw(i,1))*(1-pw(j,1))/4.
             N4=(1+pw(i,1))*(1-pw(j,1))/4.
             JACO_11=(1+pw(j,1))/4.*(x1-x2)+(1-pw(j,1))/4.*(x4-x3)
             JACO_12=(1+pw(i,1))/4.*(x1-x4)+(1-pw(i,1))/4.*(x2-x3)
             JACO_21=(1+pw(j,1))/4.*(y1-y2)+(1-pw(j,1))/4.*(y4-y3)
             JACO_22=(1+pw(i,1))/4.*(y1-y4)+(1-pw(i,1))/4.*(y2-y3)
             JACO=JACO_11*JACO_22-JACO_12*JACO_21
             
             kesi=x1*N1+x2*N2+x3*N3+x4*N4
             yita=y1*N1+y2*N2+y3*N3+y4*N4
            
             A=(localx-kesi)**2+(localy-yita)**2+localz**2
             B=pw(i,2)*pw(j,2)*JACO
             
             f2=f2+B*A**(-1.5)*localz
        end do 
    end do


!Calculate f3 
    f3=0.
    do i=1,M
        do j=1,M
             N1=(1+pw(i,1))*(1+pw(j,1))/4.
             N2=(1-pw(i,1))*(1+pw(j,1))/4.
             N3=(1-pw(i,1))*(1-pw(j,1))/4.
             N4=(1+pw(i,1))*(1-pw(j,1))/4.
             JACO_11=(1+pw(j,1))/4.*(x1-x2)+(1-pw(j,1))/4.*(x4-x3)
             JACO_12=(1+pw(i,1))/4.*(x1-x4)+(1-pw(i,1))/4.*(x2-x3)
             JACO_21=(1+pw(j,1))/4.*(y1-y2)+(1-pw(j,1))/4.*(y4-y3)
             JACO_22=(1+pw(i,1))/4.*(y1-y4)+(1-pw(i,1))/4.*(y2-y3)
             JACO=JACO_11*JACO_22-JACO_12*JACO_21
             
             kesi=x1*N1+x2*N2+x3*N3+x4*N4
             yita=y1*N1+y2*N2+y3*N3+y4*N4
            
             A=(localx-kesi)**2+(localy-yita)**2+localz**2
             B=pw(i,2)*pw(j,2)*JACO
             
             f3=f3+B*A**(-1.5)*(localx-kesi)
        end do 
    end do
  
!Calculate f4 
    f4=0. 
    do i=1,M
        do j=1,M
             N1=(1+pw(i,1))*(1+pw(j,1))/4.
             N2=(1-pw(i,1))*(1+pw(j,1))/4.
             N3=(1-pw(i,1))*(1-pw(j,1))/4.
             N4=(1+pw(i,1))*(1-pw(j,1))/4.
             JACO_11=(1+pw(j,1))/4.*(x1-x2)+(1-pw(j,1))/4.*(x4-x3)
             JACO_12=(1+pw(i,1))/4.*(x1-x4)+(1-pw(i,1))/4.*(x2-x3)
             JACO_21=(1+pw(j,1))/4.*(y1-y2)+(1-pw(j,1))/4.*(y4-y3)
             JACO_22=(1+pw(i,1))/4.*(y1-y4)+(1-pw(i,1))/4.*(y2-y3)
             JACO=JACO_11*JACO_22-JACO_12*JACO_21
             
             kesi=x1*N1+x2*N2+x3*N3+x4*N4
             yita=y1*N1+y2*N2+y3*N3+y4*N4
            
             A=(localx-kesi)**2+(localy-yita)**2+localz**2
             B=pw(i,2)*pw(j,2)*JACO
             
             f4=f4+B*A**(-1.5)*(localy-yita)
        end do 
    end do      

!Calculate f5 
    f5=0. 
    do i=1,M
        do j=1,M
             N1=(1+pw(i,1))*(1+pw(j,1))/4.
             N2=(1-pw(i,1))*(1+pw(j,1))/4.
             N3=(1-pw(i,1))*(1-pw(j,1))/4.
             N4=(1+pw(i,1))*(1-pw(j,1))/4.
             JACO_11=(1+pw(j,1))/4.*(x1-x2)+(1-pw(j,1))/4.*(x4-x3)
             JACO_12=(1+pw(i,1))/4.*(x1-x4)+(1-pw(i,1))/4.*(x2-x3)
             JACO_21=(1+pw(j,1))/4.*(y1-y2)+(1-pw(j,1))/4.*(y4-y3)
             JACO_22=(1+pw(i,1))/4.*(y1-y4)+(1-pw(i,1))/4.*(y2-y3)
             JACO=JACO_11*JACO_22-JACO_12*JACO_21
             
             kesi=x1*N1+x2*N2+x3*N3+x4*N4
             yita=y1*N1+y2*N2+y3*N3+y4*N4
            
             A=(localx-kesi)**2+(localy-yita)**2+localz**2
             B=pw(i,2)*pw(j,2)*JACO
             
             f5=f5+B*A**(-1.5)*(localx-kesi)**2
        end do 
    end do
  
!Calculate f6 
    f6=0. 
    do i=1,M
        do j=1,M
             N1=(1+pw(i,1))*(1+pw(j,1))/4.
             N2=(1-pw(i,1))*(1+pw(j,1))/4.
             N3=(1-pw(i,1))*(1-pw(j,1))/4.
             N4=(1+pw(i,1))*(1-pw(j,1))/4.
             JACO_11=(1+pw(j,1))/4.*(x1-x2)+(1-pw(j,1))/4.*(x4-x3)
             JACO_12=(1+pw(i,1))/4.*(x1-x4)+(1-pw(i,1))/4.*(x2-x3)
             JACO_21=(1+pw(j,1))/4.*(y1-y2)+(1-pw(j,1))/4.*(y4-y3)
             JACO_22=(1+pw(i,1))/4.*(y1-y4)+(1-pw(i,1))/4.*(y2-y3)
             JACO=JACO_11*JACO_22-JACO_12*JACO_21
             
             kesi=x1*N1+x2*N2+x3*N3+x4*N4
             yita=y1*N1+y2*N2+y3*N3+y4*N4
            
             A=(localx-kesi)**2+(localy-yita)**2+localz**2
             B=pw(i,2)*pw(j,2)*JACO
             
             f6=f6+B*A**(-1.5)*(localy-yita)**2
        end do 
    end do  

!Calculate f7  
    f7=0. 
    do i=1,M
        do j=1,M
             N1=(1+pw(i,1))*(1+pw(j,1))/4.
             N2=(1-pw(i,1))*(1+pw(j,1))/4.
             N3=(1-pw(i,1))*(1-pw(j,1))/4.
             N4=(1+pw(i,1))*(1-pw(j,1))/4.
             JACO_11=(1+pw(j,1))/4.*(x1-x2)+(1-pw(j,1))/4.*(x4-x3)
             JACO_12=(1+pw(i,1))/4.*(x1-x4)+(1-pw(i,1))/4.*(x2-x3)
             JACO_21=(1+pw(j,1))/4.*(y1-y2)+(1-pw(j,1))/4.*(y4-y3)
             JACO_22=(1+pw(i,1))/4.*(y1-y4)+(1-pw(i,1))/4.*(y2-y3)
             JACO=JACO_11*JACO_22-JACO_12*JACO_21
             
             kesi=x1*N1+x2*N2+x3*N3+x4*N4
             yita=y1*N1+y2*N2+y3*N3+y4*N4
            
             A=(localx-kesi)**2+(localy-yita)**2+localz**2
             B=pw(i,2)*pw(j,2)*JACO
             
             f7=f7+B*A**(-1.5)*(localy-yita)*(localx-kesi)
        end do 
    end do    
        
!Calculate f8 
    f8=0. 
    do i=1,M
        do j=1,M
             N1=(1+pw(i,1))*(1+pw(j,1))/4.
             N2=(1-pw(i,1))*(1+pw(j,1))/4.
             N3=(1-pw(i,1))*(1-pw(j,1))/4.
             N4=(1+pw(i,1))*(1-pw(j,1))/4.
             JACO_11=(1+pw(j,1))/4.*(x1-x2)+(1-pw(j,1))/4.*(x4-x3)
             JACO_12=(1+pw(i,1))/4.*(x1-x4)+(1-pw(i,1))/4.*(x2-x3)
             JACO_21=(1+pw(j,1))/4.*(y1-y2)+(1-pw(j,1))/4.*(y4-y3)
             JACO_22=(1+pw(i,1))/4.*(y1-y4)+(1-pw(i,1))/4.*(y2-y3)
             JACO=JACO_11*JACO_22-JACO_12*JACO_21
             
             kesi=x1*N1+x2*N2+x3*N3+x4*N4
             yita=y1*N1+y2*N2+y3*N3+y4*N4
            
             A=(localx-kesi)**2+(localy-yita)**2+localz**2
             B=pw(i,2)*pw(j,2)*JACO
             
             f8=f8+B*A**(-2.5)*localz**3
        end do 
    end do   
   
!Calculate f9    
    f9=0. 
    do i=1,M
        do j=1,M
             N1=(1+pw(i,1))*(1+pw(j,1))/4.
             N2=(1-pw(i,1))*(1+pw(j,1))/4.
             N3=(1-pw(i,1))*(1-pw(j,1))/4.
             N4=(1+pw(i,1))*(1-pw(j,1))/4.
             JACO_11=(1+pw(j,1))/4.*(x1-x2)+(1-pw(j,1))/4.*(x4-x3)
             JACO_12=(1+pw(i,1))/4.*(x1-x4)+(1-pw(i,1))/4.*(x2-x3)
             JACO_21=(1+pw(j,1))/4.*(y1-y2)+(1-pw(j,1))/4.*(y4-y3)
             JACO_22=(1+pw(i,1))/4.*(y1-y4)+(1-pw(i,1))/4.*(y2-y3)
             JACO=JACO_11*JACO_22-JACO_12*JACO_21
             
             kesi=x1*N1+x2*N2+x3*N3+x4*N4
             yita=y1*N1+y2*N2+y3*N3+y4*N4
            
             A=(localx-kesi)**2+(localy-yita)**2+localz**2
             B=pw(i,2)*pw(j,2)*JACO
             
             f9=f9+B*A**(-2.5)*(localx-kesi)*localz**2
        end do    
	end do 
	
!Calculate f10
    f10=0. 
    do i=1,M
        do j=1,M
             N1=(1+pw(i,1))*(1+pw(j,1))/4.
             N2=(1-pw(i,1))*(1+pw(j,1))/4.
             N3=(1-pw(i,1))*(1-pw(j,1))/4.
             N4=(1+pw(i,1))*(1-pw(j,1))/4.
             JACO_11=(1+pw(j,1))/4.*(x1-x2)+(1-pw(j,1))/4.*(x4-x3)
             JACO_12=(1+pw(i,1))/4.*(x1-x4)+(1-pw(i,1))/4.*(x2-x3)
             JACO_21=(1+pw(j,1))/4.*(y1-y2)+(1-pw(j,1))/4.*(y4-y3)
             JACO_22=(1+pw(i,1))/4.*(y1-y4)+(1-pw(i,1))/4.*(y2-y3)
             JACO=JACO_11*JACO_22-JACO_12*JACO_21
             
             kesi=x1*N1+x2*N2+x3*N3+x4*N4
             yita=y1*N1+y2*N2+y3*N3+y4*N4
            
             A=(localx-kesi)**2+(localy-yita)**2+localz**2
             B=pw(i,2)*pw(j,2)*JACO
             
             f10=f10+B*A**(-2.5)*(localy-yita)*localz**2
        end do    
	end do 
	
!Calculate f11    
	f11=0. 
    do i=1,M
        do j=1,M 
             N1=(1+pw(i,1))*(1+pw(j,1))/4.
             N2=(1-pw(i,1))*(1+pw(j,1))/4.
             N3=(1-pw(i,1))*(1-pw(j,1))/4.
             N4=(1+pw(i,1))*(1-pw(j,1))/4.
             JACO_11=(1+pw(j,1))/4.*(x1-x2)+(1-pw(j,1))/4.*(x4-x3)
             JACO_12=(1+pw(i,1))/4.*(x1-x4)+(1-pw(i,1))/4.*(x2-x3)
             JACO_21=(1+pw(j,1))/4.*(y1-y2)+(1-pw(j,1))/4.*(y4-y3)
             JACO_22=(1+pw(i,1))/4.*(y1-y4)+(1-pw(i,1))/4.*(y2-y3)
             JACO=JACO_11*JACO_22-JACO_12*JACO_21
             
             kesi=x1*N1+x2*N2+x3*N3+x4*N4
             yita=y1*N1+y2*N2+y3*N3+y4*N4
            
             A=(localx-kesi)**2+(localy-yita)**2+localz**2
             B=pw(i,2)*pw(j,2)*JACO
             
             f11=f11+B*A**(-2.5)*(localx-kesi)**2*localz
        end do    
	end do   
	
!Calculate f12   
    f12=0. 
    do i=1,M
        do j=1,M
             N1=(1+pw(i,1))*(1+pw(j,1))/4.
             N2=(1-pw(i,1))*(1+pw(j,1))/4.
             N3=(1-pw(i,1))*(1-pw(j,1))/4.
             N4=(1+pw(i,1))*(1-pw(j,1))/4.
             JACO_11=(1+pw(j,1))/4.*(x1-x2)+(1-pw(j,1))/4.*(x4-x3)
             JACO_12=(1+pw(i,1))/4.*(x1-x4)+(1-pw(i,1))/4.*(x2-x3)
             JACO_21=(1+pw(j,1))/4.*(y1-y2)+(1-pw(j,1))/4.*(y4-y3)
             JACO_22=(1+pw(i,1))/4.*(y1-y4)+(1-pw(i,1))/4.*(y2-y3)
             JACO=JACO_11*JACO_22-JACO_12*JACO_21
             
             kesi=x1*N1+x2*N2+x3*N3+x4*N4
             yita=y1*N1+y2*N2+y3*N3+y4*N4
            
             A=(localx-kesi)**2+(localy-yita)**2+localz**2
             B=pw(i,2)*pw(j,2)*JACO
             
             f12=f12+B*A**(-2.5)*(localy-yita)**2*localz
        end do    
	end do 
	
!Calculate f13
	f13=0. 
    do i=1,M
        do j=1,M
             N1=(1+pw(i,1))*(1+pw(j,1))/4.
             N2=(1-pw(i,1))*(1+pw(j,1))/4.
             N3=(1-pw(i,1))*(1-pw(j,1))/4.
             N4=(1+pw(i,1))*(1-pw(j,1))/4.
             JACO_11=(1+pw(j,1))/4.*(x1-x2)+(1-pw(j,1))/4.*(x4-x3)
             JACO_12=(1+pw(i,1))/4.*(x1-x4)+(1-pw(i,1))/4.*(x2-x3)
             JACO_21=(1+pw(j,1))/4.*(y1-y2)+(1-pw(j,1))/4.*(y4-y3)
             JACO_22=(1+pw(i,1))/4.*(y1-y4)+(1-pw(i,1))/4.*(y2-y3)
             JACO=JACO_11*JACO_22-JACO_12*JACO_21
             
             kesi=x1*N1+x2*N2+x3*N3+x4*N4
             yita=y1*N1+y2*N2+y3*N3+y4*N4
            
             A=(localx-kesi)**2+(localy-yita)**2+localz**2
             B=pw(i,2)*pw(j,2)*JACO
             
             f13=f13+B*A**(-2.5)*(localx-kesi)*(localy-yita)
        end do    
    end do    

!Calculate f14
    f14=0. 
    do i=1,M
        do j=1,M
             N1=(1+pw(i,1))*(1+pw(j,1))/4.
             N2=(1-pw(i,1))*(1+pw(j,1))/4.
             N3=(1-pw(i,1))*(1-pw(j,1))/4.
             N4=(1+pw(i,1))*(1-pw(j,1))/4.
             JACO_11=(1+pw(j,1))/4.*(x1-x2)+(1-pw(j,1))/4.*(x4-x3)
             JACO_12=(1+pw(i,1))/4.*(x1-x4)+(1-pw(i,1))/4.*(x2-x3)
             JACO_21=(1+pw(j,1))/4.*(y1-y2)+(1-pw(j,1))/4.*(y4-y3)
             JACO_22=(1+pw(i,1))/4.*(y1-y4)+(1-pw(i,1))/4.*(y2-y3)
             JACO=JACO_11*JACO_22-JACO_12*JACO_21
             
             kesi=x1*N1+x2*N2+x3*N3+x4*N4
             yita=y1*N1+y2*N2+y3*N3+y4*N4
            
             A=(localx-kesi)**2+(localy-yita)**2+localz**2
             B=pw(i,2)*pw(j,2)*JACO
             
             f14=f14+B*A**(-2.5)*(localx-kesi)**3
        end do    
    end do  
    
!Calculate f15
    f15=0. 
    do i=1,M
        do j=1,M
             N1=(1+pw(i,1))*(1+pw(j,1))/4.
             N2=(1-pw(i,1))*(1+pw(j,1))/4.
             N3=(1-pw(i,1))*(1-pw(j,1))/4.
             N4=(1+pw(i,1))*(1-pw(j,1))/4.
             JACO_11=(1+pw(j,1))/4.*(x1-x2)+(1-pw(j,1))/4.*(x4-x3)
             JACO_12=(1+pw(i,1))/4.*(x1-x4)+(1-pw(i,1))/4.*(x2-x3)
             JACO_21=(1+pw(j,1))/4.*(y1-y2)+(1-pw(j,1))/4.*(y4-y3)
             JACO_22=(1+pw(i,1))/4.*(y1-y4)+(1-pw(i,1))/4.*(y2-y3)
             JACO=JACO_11*JACO_22-JACO_12*JACO_21
             
             kesi=x1*N1+x2*N2+x3*N3+x4*N4
             yita=y1*N1+y2*N2+y3*N3+y4*N4
            
             A=(localx-kesi)**2+(localy-yita)**2+localz**2
             B=pw(i,2)*pw(j,2)*JACO
             
             f15=f15+B*A**(-2.5)*(localy-yita)**3
        end do    
    end do    

!Calculate f16
    f16=0. 
    do i=1,M
        do j=1,M
             N1=(1+pw(i,1))*(1+pw(j,1))/4.
             N2=(1-pw(i,1))*(1+pw(j,1))/4.
             N3=(1-pw(i,1))*(1-pw(j,1))/4.
             N4=(1+pw(i,1))*(1-pw(j,1))/4.
             JACO_11=(1+pw(j,1))/4.*(x1-x2)+(1-pw(j,1))/4.*(x4-x3)
             JACO_12=(1+pw(i,1))/4.*(x1-x4)+(1-pw(i,1))/4.*(x2-x3)
             JACO_21=(1+pw(j,1))/4.*(y1-y2)+(1-pw(j,1))/4.*(y4-y3)
             JACO_22=(1+pw(i,1))/4.*(y1-y4)+(1-pw(i,1))/4.*(y2-y3)
             JACO=JACO_11*JACO_22-JACO_12*JACO_21
             
             kesi=x1*N1+x2*N2+x3*N3+x4*N4
             yita=y1*N1+y2*N2+y3*N3+y4*N4
            
             A=(localx-kesi)**2+(localy-yita)**2+localz**2
             B=pw(i,2)*pw(j,2)*JACO
             
             f16=f16+B*A**(-2.5)*(localx-kesi)*(localy-yita)**2
        end do    
    end do    

!Calculate f17
    f17=0. 
    do i=1,M
        do j=1,M  
             N1=(1+pw(i,1))*(1+pw(j,1))/4.
             N2=(1-pw(i,1))*(1+pw(j,1))/4.
             N3=(1-pw(i,1))*(1-pw(j,1))/4.
             N4=(1+pw(i,1))*(1-pw(j,1))/4.
             JACO_11=(1+pw(j,1))/4.*(x1-x2)+(1-pw(j,1))/4.*(x4-x3)
             JACO_12=(1+pw(i,1))/4.*(x1-x4)+(1-pw(i,1))/4.*(x2-x3)
             JACO_21=(1+pw(j,1))/4.*(y1-y2)+(1-pw(j,1))/4.*(y4-y3)
             JACO_22=(1+pw(i,1))/4.*(y1-y4)+(1-pw(i,1))/4.*(y2-y3)
             JACO=JACO_11*JACO_22-JACO_12*JACO_21
             
             kesi=x1*N1+x2*N2+x3*N3+x4*N4
             yita=y1*N1+y2*N2+y3*N3+y4*N4
            
             A=(localx-kesi)**2+(localy-yita)**2+localz**2
             B=pw(i,2)*pw(j,2)*JACO
             
             f17=f17+B*A**(-2.5)*(localx-kesi)**2*(localy-yita)
        end do    
    end do    
!************************************************************
! calculate influence coefficients 
!************************************************************
      S(1,1) =-3.0D0 * F14 - ( 1.0D0 - 2.0D0 * V ) * F3
      S(1,2) =-3.0D0 * F17 + ( 1.0D0 - 2.0D0 * V ) * F4
      S(1,3) =-3.0D0 * F11 + ( 1.0D0 - 2.0D0 * V ) * F2
      S(2,1) =-3.0D0 * F16 + ( 1.0D0 - 2.0D0 * V ) * F3
      S(2,2) =-3.0D0 * F15 - ( 1.0D0 - 2.0D0 * V ) * F4
      S(2,3) =-3.0D0 * F12 + ( 1.0D0 - 2.0D0 * V ) * F2
      S(3,1) =-3.0D0 * F9  + ( 1.0D0 - 2.0D0 * V ) * F3
      S(3,2) =-3.0D0 * F10 + ( 1.0D0 - 2.0D0 * V ) * F4
      S(3,3) =-3.0D0 * F8  - ( 1.0D0 - 2.0D0 * V ) * F2
      S(4,1) =-3.0D0 * F17 - ( 1.0D0 - 2.0D0 * V ) * F4
      S(4,2) =-3.0D0 * F16 - ( 1.0D0 - 2.0D0 * V ) * F3
      S(4,3) =-3.0D0 * F13 * localZ
      S(5,1) =-3.0D0 * F13 * localZ
      S(5,2) =-3.0D0 * F12 - ( 1.0D0 - 2.0D0 * V ) * F2
      S(5,3) =-3.0D0 * F10 - ( 1.0D0 - 2.0D0 * V ) * F4
      S(6,1) =-3.0D0 * F11 - ( 1.0D0 - 2.0D0 * V ) * F2
      S(6,2) =-3.0D0 * F13 * localZ
      S(6,3) =-3.0D0 * F9  - ( 1.0D0 - 2.0D0 * V ) * F3
   
      JI(1,1)=JEvect(1,1)*IEvect(1,1)+JEvect(1,2)*IEvect(1,2)+JEvect(1,3)*IEvect(1,3)
	  JI(2,1)=JEvect(2,1)*IEvect(1,1)+JEvect(2,2)*IEvect(1,2)+JEvect(2,3)*IEvect(1,3)
      JI(3,1)=JEvect(3,1)*IEvect(1,1)+JEvect(3,2)*IEvect(1,2)+JEvect(3,3)*IEvect(1,3)
      
      JI(1,2)=JEvect(1,1)*IEvect(2,1)+JEvect(1,2)*IEvect(2,2)+JEvect(1,3)*IEvect(2,3)
      JI(2,2)=JEvect(2,1)*IEvect(2,1)+JEvect(2,2)*IEvect(2,2)+JEvect(2,3)*IEvect(2,3)
      JI(3,2)=JEvect(3,1)*IEvect(2,1)+JEvect(3,2)*IEvect(2,2)+JEvect(3,3)*IEvect(2,3) 
      
      JI(1,3)=JEvect(1,1)*IEvect(3,1)+JEvect(1,2)*IEvect(3,2)+JEvect(1,3)*IEvect(3,3)
      JI(2,3)=JEvect(2,1)*IEvect(3,1)+JEvect(2,2)*IEvect(3,2)+JEvect(2,3)*IEvect(3,3)
      JI(3,3)=JEvect(3,1)*IEvect(3,1)+JEvect(3,2)*IEvect(3,2)+JEvect(3,3)*IEvect(3,3)  
  
      Trac(1)=S(1,1)*JI(1,3)+S(4,1)*JI(2,3)+S(6,1)*JI(3,3)
      Trac(2)=S(4,1)*JI(1,3)+S(2,1)*JI(2,3)+S(5,1)*JI(3,3)
      Trac(3)=S(6,1)*JI(1,3)+S(5,1)*JI(2,3)+S(3,1)*JI(3,3)
     
      COEFF(1,1)=cons*(Trac(1)*JI(1,1)+Trac(2)*JI(2,1)+Trac(3)*JI(3,1))
      COEFF(2,1)=cons*(Trac(1)*JI(1,2)+Trac(2)*JI(2,2)+Trac(3)*JI(3,2))
      COEFF(3,1)=cons*(Trac(1)*JI(1,3)+Trac(2)*JI(2,3)+Trac(3)*JI(3,3))
     
      Trac(1)=S(1,2)*JI(1,3)+S(4,2)*JI(2,3)+S(6,2)*JI(3,3)
      Trac(2)=S(4,2)*JI(1,3)+S(2,2)*JI(2,3)+S(5,2)*JI(3,3)
      Trac(3)=S(6,2)*JI(1,3)+S(5,2)*JI(2,3)+S(3,2)*JI(3,3)
     
      COEFF(1,2)=cons*(Trac(1)*JI(1,1)+Trac(2)*JI(2,1)+Trac(3)*JI(3,1))
      COEFF(2,2)=cons*(Trac(1)*JI(1,2)+Trac(2)*JI(2,2)+Trac(3)*JI(3,2))
      COEFF(3,2)=cons*(Trac(1)*JI(1,3)+Trac(2)*JI(2,3)+Trac(3)*JI(3,3))
     
      Trac(1)=S(1,3)*JI(1,3)+S(4,3)*JI(2,3)+S(6,3)*JI(3,3)
      Trac(2)=S(4,3)*JI(1,3)+S(2,3)*JI(2,3)+S(5,3)*JI(3,3)
      Trac(3)=S(6,3)*JI(1,3)+S(5,3)*JI(2,3)+S(3,3)*JI(3,3)
     
      COEFF(1,3)=cons*(Trac(1)*JI(1,1)+Trac(2)*JI(2,1)+Trac(3)*JI(3,1))
      COEFF(2,3)=cons*(Trac(1)*JI(1,2)+Trac(2)*JI(2,2)+Trac(3)*JI(3,2))
      COEFF(3,3)=cons*(Trac(1)*JI(1,3)+Trac(2)*JI(2,3)+Trac(3)*JI(3,3))

end subroutine stress_regular_fs3d_qua_num_part1


subroutine stress_regular_fs3d_qua_num_part3 (xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4,gpx,gpy,gpz,E,v,M,IEvect,JEvect,COEFF)
implicit none 

real*8::E,v,G,cons,l1,l2,l3,p,area,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4,gpx,gpy,gpz
real*8::VecIJ_x,VecIJ_y,VecIJ_z,localx,localy,localz,IEvect(3,3),JEvect(3,3),COEFF(6,3)
real*8::N1,N2,N3,N4,jaco_11,jaco_12,jaco_21,jaco_22,jaco,gpxx,gpyy,gpzz,vax,vay,vaz,vbx,vby,vbz,vcx,vcy,vcz,vdx,vdy,vdz
real*8::f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,kesi,yita,A,B,S(6,3),JI(3,3),TRAC(3)
integer::i,j,k,M

real*8,allocatable::pw(:,:)

    G=E/2./(1+v)
    cons=1./8./3.141592653/(1-v)
    
    gpxx=(xx1+xx2+xx3+xx4)/4.
    gpyy=(yy1+yy2+yy3+yy4)/4.
    gpzz=(zz1+zz2+zz3+zz4)/4.
    
    VAX=xx1-gpxx
    VAY=yy1-gpyy
    VAZ=zz1-gpzz

    x1=VAX*JEvect(1,1)+VAY*JEvect(1,2)+VAZ*JEvect(1,3)
    y1=VAX*JEvect(2,1)+VAY*JEvect(2,2)+VAZ*JEvect(2,3)
   
    VBX=xx2-gpxx
    VBY=yy2-gpyy
    VBZ=zz2-gpzz

    x2=VBX*JEvect(1,1)+VBY*JEvect(1,2)+VBZ*JEvect(1,3)
    y2=VBX*JEvect(2,1)+VBY*JEvect(2,2)+VBZ*JEvect(2,3) 
    
    VCX=xx3-gpxx
    VCY=yy3-gpyy
    VCZ=zz3-gpzz

    x3=VCX*JEvect(1,1)+VCY*JEvect(1,2)+VCZ*JEvect(1,3)
    y3=VCX*JEvect(2,1)+VCY*JEvect(2,2)+VCZ*JEvect(2,3)
    
    VDX=xx4-gpxx
    VDY=yy4-gpyy
    VDZ=zz4-gpzz

    x4=VDX*JEvect(1,1)+VDY*JEvect(1,2)+VDZ*JEvect(1,3)
    y4=VDX*JEvect(2,1)+VDY*JEvect(2,2)+VDZ*JEvect(2,3)
    
    VecIJ_x=gpx-gpxx
    VecIJ_y=gpy-gpyy
    VecIJ_z=gpz-gpzz
    
    Localx=VecIJ_x*JEvect(1,1)+VecIJ_y*JEvect(1,2)+VecIJ_z*JEvect(1,3)
    Localy=VecIJ_x*JEvect(2,1)+VecIJ_y*JEvect(2,2)+VecIJ_z*JEvect(2,3)
    Localz=VecIJ_x*JEvect(3,1)+VecIJ_y*JEvect(3,2)+VecIJ_z*JEvect(3,3)
    
!************************************************************** 
allocate (pw(M,2))
call Gauss_ponits_weights(M,pw)
!*****************************************************************************
!Calculate f1   
    f1=0.
    do i=1,M
        do j=1,M
             N1=(1+pw(i,1))*(1+pw(j,1))/4.
             N2=(1-pw(i,1))*(1+pw(j,1))/4.
             N3=(1-pw(i,1))*(1-pw(j,1))/4.
             N4=(1+pw(i,1))*(1-pw(j,1))/4.
             JACO_11=(1+pw(j,1))/4.*(x1-x2)+(1-pw(j,1))/4.*(x4-x3)
             JACO_12=(1+pw(i,1))/4.*(x1-x4)+(1-pw(i,1))/4.*(x2-x3)
             JACO_21=(1+pw(j,1))/4.*(y1-y2)+(1-pw(j,1))/4.*(y4-y3)
             JACO_22=(1+pw(i,1))/4.*(y1-y4)+(1-pw(i,1))/4.*(y2-y3)
             JACO=JACO_11*JACO_22-JACO_12*JACO_21
             
             kesi=x1*N1+x2*N2+x3*N3+x4*N4
             yita=y1*N1+y2*N2+y3*N3+y4*N4
            
             A=(localx-kesi)**2+(localy-yita)**2+localz**2
             B=pw(i,2)*pw(j,2)*JACO
             
             f1=f1+B*A**(-0.5)            
        end do 
    end do   

!Calculate f2
    f2=0.
    do i=1,M
        do j=1,M
             N1=(1+pw(i,1))*(1+pw(j,1))/4.
             N2=(1-pw(i,1))*(1+pw(j,1))/4.
             N3=(1-pw(i,1))*(1-pw(j,1))/4.
             N4=(1+pw(i,1))*(1-pw(j,1))/4.
             JACO_11=(1+pw(j,1))/4.*(x1-x2)+(1-pw(j,1))/4.*(x4-x3)
             JACO_12=(1+pw(i,1))/4.*(x1-x4)+(1-pw(i,1))/4.*(x2-x3)
             JACO_21=(1+pw(j,1))/4.*(y1-y2)+(1-pw(j,1))/4.*(y4-y3)
             JACO_22=(1+pw(i,1))/4.*(y1-y4)+(1-pw(i,1))/4.*(y2-y3)
             JACO=JACO_11*JACO_22-JACO_12*JACO_21
             
             kesi=x1*N1+x2*N2+x3*N3+x4*N4
             yita=y1*N1+y2*N2+y3*N3+y4*N4
            
             A=(localx-kesi)**2+(localy-yita)**2+localz**2
             B=pw(i,2)*pw(j,2)*JACO
             
             f2=f2+B*A**(-1.5)*localz
        end do 
    end do

!Calculate f3 
    f3=0.
    do i=1,M
        do j=1,M
             N1=(1+pw(i,1))*(1+pw(j,1))/4.
             N2=(1-pw(i,1))*(1+pw(j,1))/4.
             N3=(1-pw(i,1))*(1-pw(j,1))/4.
             N4=(1+pw(i,1))*(1-pw(j,1))/4.
             JACO_11=(1+pw(j,1))/4.*(x1-x2)+(1-pw(j,1))/4.*(x4-x3)
             JACO_12=(1+pw(i,1))/4.*(x1-x4)+(1-pw(i,1))/4.*(x2-x3)
             JACO_21=(1+pw(j,1))/4.*(y1-y2)+(1-pw(j,1))/4.*(y4-y3)
             JACO_22=(1+pw(i,1))/4.*(y1-y4)+(1-pw(i,1))/4.*(y2-y3)
             JACO=JACO_11*JACO_22-JACO_12*JACO_21
             
             kesi=x1*N1+x2*N2+x3*N3+x4*N4
             yita=y1*N1+y2*N2+y3*N3+y4*N4
            
             A=(localx-kesi)**2+(localy-yita)**2+localz**2
             B=pw(i,2)*pw(j,2)*JACO
             
             f3=f3+B*A**(-1.5)*(localx-kesi)
        end do 
    end do
  
!Calculate f4 
    f4=0. 
    do i=1,M
        do j=1,M
             N1=(1+pw(i,1))*(1+pw(j,1))/4.
             N2=(1-pw(i,1))*(1+pw(j,1))/4.
             N3=(1-pw(i,1))*(1-pw(j,1))/4.
             N4=(1+pw(i,1))*(1-pw(j,1))/4.
             JACO_11=(1+pw(j,1))/4.*(x1-x2)+(1-pw(j,1))/4.*(x4-x3)
             JACO_12=(1+pw(i,1))/4.*(x1-x4)+(1-pw(i,1))/4.*(x2-x3)
             JACO_21=(1+pw(j,1))/4.*(y1-y2)+(1-pw(j,1))/4.*(y4-y3)
             JACO_22=(1+pw(i,1))/4.*(y1-y4)+(1-pw(i,1))/4.*(y2-y3)
             JACO=JACO_11*JACO_22-JACO_12*JACO_21
             
             kesi=x1*N1+x2*N2+x3*N3+x4*N4
             yita=y1*N1+y2*N2+y3*N3+y4*N4
            
             A=(localx-kesi)**2+(localy-yita)**2+localz**2
             B=pw(i,2)*pw(j,2)*JACO
             
             f4=f4+B*A**(-1.5)*(localy-yita)
        end do 
    end do      

!Calculate f5 
    f5=0. 
    do i=1,M
        do j=1,M
             N1=(1+pw(i,1))*(1+pw(j,1))/4.
             N2=(1-pw(i,1))*(1+pw(j,1))/4.
             N3=(1-pw(i,1))*(1-pw(j,1))/4.
             N4=(1+pw(i,1))*(1-pw(j,1))/4.
             JACO_11=(1+pw(j,1))/4.*(x1-x2)+(1-pw(j,1))/4.*(x4-x3)
             JACO_12=(1+pw(i,1))/4.*(x1-x4)+(1-pw(i,1))/4.*(x2-x3)
             JACO_21=(1+pw(j,1))/4.*(y1-y2)+(1-pw(j,1))/4.*(y4-y3)
             JACO_22=(1+pw(i,1))/4.*(y1-y4)+(1-pw(i,1))/4.*(y2-y3)
             JACO=JACO_11*JACO_22-JACO_12*JACO_21
             
             kesi=x1*N1+x2*N2+x3*N3+x4*N4
             yita=y1*N1+y2*N2+y3*N3+y4*N4
            
             A=(localx-kesi)**2+(localy-yita)**2+localz**2
             B=pw(i,2)*pw(j,2)*JACO
             
             f5=f5+B*A**(-1.5)*(localx-kesi)**2
        end do 
    end do
  
!Calculate f6 
    f6=0. 
    do i=1,M
        do j=1,M
             N1=(1+pw(i,1))*(1+pw(j,1))/4.
             N2=(1-pw(i,1))*(1+pw(j,1))/4.
             N3=(1-pw(i,1))*(1-pw(j,1))/4.
             N4=(1+pw(i,1))*(1-pw(j,1))/4.
             JACO_11=(1+pw(j,1))/4.*(x1-x2)+(1-pw(j,1))/4.*(x4-x3)
             JACO_12=(1+pw(i,1))/4.*(x1-x4)+(1-pw(i,1))/4.*(x2-x3)
             JACO_21=(1+pw(j,1))/4.*(y1-y2)+(1-pw(j,1))/4.*(y4-y3)
             JACO_22=(1+pw(i,1))/4.*(y1-y4)+(1-pw(i,1))/4.*(y2-y3)
             JACO=JACO_11*JACO_22-JACO_12*JACO_21
             
             kesi=x1*N1+x2*N2+x3*N3+x4*N4
             yita=y1*N1+y2*N2+y3*N3+y4*N4
            
             A=(localx-kesi)**2+(localy-yita)**2+localz**2
             B=pw(i,2)*pw(j,2)*JACO
             
             f6=f6+B*A**(-1.5)*(localy-yita)**2
        end do 
    end do  

!Calculate f7  
    f7=0. 
    do i=1,M
        do j=1,M
             N1=(1+pw(i,1))*(1+pw(j,1))/4.
             N2=(1-pw(i,1))*(1+pw(j,1))/4.
             N3=(1-pw(i,1))*(1-pw(j,1))/4.
             N4=(1+pw(i,1))*(1-pw(j,1))/4.
             JACO_11=(1+pw(j,1))/4.*(x1-x2)+(1-pw(j,1))/4.*(x4-x3)
             JACO_12=(1+pw(i,1))/4.*(x1-x4)+(1-pw(i,1))/4.*(x2-x3)
             JACO_21=(1+pw(j,1))/4.*(y1-y2)+(1-pw(j,1))/4.*(y4-y3)
             JACO_22=(1+pw(i,1))/4.*(y1-y4)+(1-pw(i,1))/4.*(y2-y3)
             JACO=JACO_11*JACO_22-JACO_12*JACO_21
             
             kesi=x1*N1+x2*N2+x3*N3+x4*N4
             yita=y1*N1+y2*N2+y3*N3+y4*N4
            
             A=(localx-kesi)**2+(localy-yita)**2+localz**2
             B=pw(i,2)*pw(j,2)*JACO
             
             f7=f7+B*A**(-1.5)*(localy-yita)*(localx-kesi)
        end do 
    end do    
        
!Calculate f8 
    f8=0. 
    do i=1,M
        do j=1,M
             N1=(1+pw(i,1))*(1+pw(j,1))/4.
             N2=(1-pw(i,1))*(1+pw(j,1))/4.
             N3=(1-pw(i,1))*(1-pw(j,1))/4.
             N4=(1+pw(i,1))*(1-pw(j,1))/4.
             JACO_11=(1+pw(j,1))/4.*(x1-x2)+(1-pw(j,1))/4.*(x4-x3)
             JACO_12=(1+pw(i,1))/4.*(x1-x4)+(1-pw(i,1))/4.*(x2-x3)
             JACO_21=(1+pw(j,1))/4.*(y1-y2)+(1-pw(j,1))/4.*(y4-y3)
             JACO_22=(1+pw(i,1))/4.*(y1-y4)+(1-pw(i,1))/4.*(y2-y3)
             JACO=JACO_11*JACO_22-JACO_12*JACO_21
             
             kesi=x1*N1+x2*N2+x3*N3+x4*N4
             yita=y1*N1+y2*N2+y3*N3+y4*N4
            
             A=(localx-kesi)**2+(localy-yita)**2+localz**2
             B=pw(i,2)*pw(j,2)*JACO
             
             f8=f8+B*A**(-2.5)*localz**3
        end do 
    end do   
   
!Calculate f9    
    f9=0. 
    do i=1,M
        do j=1,M
             N1=(1+pw(i,1))*(1+pw(j,1))/4.
             N2=(1-pw(i,1))*(1+pw(j,1))/4.
             N3=(1-pw(i,1))*(1-pw(j,1))/4.
             N4=(1+pw(i,1))*(1-pw(j,1))/4.
             JACO_11=(1+pw(j,1))/4.*(x1-x2)+(1-pw(j,1))/4.*(x4-x3)
             JACO_12=(1+pw(i,1))/4.*(x1-x4)+(1-pw(i,1))/4.*(x2-x3)
             JACO_21=(1+pw(j,1))/4.*(y1-y2)+(1-pw(j,1))/4.*(y4-y3)
             JACO_22=(1+pw(i,1))/4.*(y1-y4)+(1-pw(i,1))/4.*(y2-y3)
             JACO=JACO_11*JACO_22-JACO_12*JACO_21
             
             kesi=x1*N1+x2*N2+x3*N3+x4*N4
             yita=y1*N1+y2*N2+y3*N3+y4*N4
            
             A=(localx-kesi)**2+(localy-yita)**2+localz**2
             B=pw(i,2)*pw(j,2)*JACO
             
             f9=f9+B*A**(-2.5)*(localx-kesi)*localz**2
        end do    
	end do 
	
!Calculate f10
    f10=0. 
    do i=1,M
        do j=1,M
             N1=(1+pw(i,1))*(1+pw(j,1))/4.
             N2=(1-pw(i,1))*(1+pw(j,1))/4.
             N3=(1-pw(i,1))*(1-pw(j,1))/4.
             N4=(1+pw(i,1))*(1-pw(j,1))/4.
             JACO_11=(1+pw(j,1))/4.*(x1-x2)+(1-pw(j,1))/4.*(x4-x3)
             JACO_12=(1+pw(i,1))/4.*(x1-x4)+(1-pw(i,1))/4.*(x2-x3)
             JACO_21=(1+pw(j,1))/4.*(y1-y2)+(1-pw(j,1))/4.*(y4-y3)
             JACO_22=(1+pw(i,1))/4.*(y1-y4)+(1-pw(i,1))/4.*(y2-y3)
             JACO=JACO_11*JACO_22-JACO_12*JACO_21
             
             kesi=x1*N1+x2*N2+x3*N3+x4*N4
             yita=y1*N1+y2*N2+y3*N3+y4*N4
            
             A=(localx-kesi)**2+(localy-yita)**2+localz**2
             B=pw(i,2)*pw(j,2)*JACO
             
             f10=f10+B*A**(-2.5)*(localy-yita)*localz**2
        end do    
	end do 
	
!Calculate f11    
    f11=0. 
    do i=1,M
        do j=1,M
             N1=(1+pw(i,1))*(1+pw(j,1))/4.
             N2=(1-pw(i,1))*(1+pw(j,1))/4.
             N3=(1-pw(i,1))*(1-pw(j,1))/4.
             N4=(1+pw(i,1))*(1-pw(j,1))/4.
             JACO_11=(1+pw(j,1))/4.*(x1-x2)+(1-pw(j,1))/4.*(x4-x3)
             JACO_12=(1+pw(i,1))/4.*(x1-x4)+(1-pw(i,1))/4.*(x2-x3)
             JACO_21=(1+pw(j,1))/4.*(y1-y2)+(1-pw(j,1))/4.*(y4-y3)
             JACO_22=(1+pw(i,1))/4.*(y1-y4)+(1-pw(i,1))/4.*(y2-y3)
             JACO=JACO_11*JACO_22-JACO_12*JACO_21
             
             kesi=x1*N1+x2*N2+x3*N3+x4*N4
             yita=y1*N1+y2*N2+y3*N3+y4*N4
            
             A=(localx-kesi)**2+(localy-yita)**2+localz**2
             B=pw(i,2)*pw(j,2)*JACO
             
             f11=f11+B*A**(-2.5)*(localx-kesi)**2*localz
        end do    
	end do  
	
!Calculate f12   
    f12=0. 
    do i=1,M
        do j=1,M
             N1=(1+pw(i,1))*(1+pw(j,1))/4.
             N2=(1-pw(i,1))*(1+pw(j,1))/4.
             N3=(1-pw(i,1))*(1-pw(j,1))/4.
             N4=(1+pw(i,1))*(1-pw(j,1))/4.
             JACO_11=(1+pw(j,1))/4.*(x1-x2)+(1-pw(j,1))/4.*(x4-x3)
             JACO_12=(1+pw(i,1))/4.*(x1-x4)+(1-pw(i,1))/4.*(x2-x3)
             JACO_21=(1+pw(j,1))/4.*(y1-y2)+(1-pw(j,1))/4.*(y4-y3)
             JACO_22=(1+pw(i,1))/4.*(y1-y4)+(1-pw(i,1))/4.*(y2-y3)
             JACO=JACO_11*JACO_22-JACO_12*JACO_21
             
             kesi=x1*N1+x2*N2+x3*N3+x4*N4
             yita=y1*N1+y2*N2+y3*N3+y4*N4
            
             A=(localx-kesi)**2+(localy-yita)**2+localz**2
             B=pw(i,2)*pw(j,2)*JACO
             
             f12=f12+B*A**(-2.5)*(localy-yita)**2*localz
        end do    
	end do  
	
!Calculate f13
    f13=0. 
    do i=1,M
        do j=1,M 
             N1=(1+pw(i,1))*(1+pw(j,1))/4.
             N2=(1-pw(i,1))*(1+pw(j,1))/4.
             N3=(1-pw(i,1))*(1-pw(j,1))/4.
             N4=(1+pw(i,1))*(1-pw(j,1))/4.
             JACO_11=(1+pw(j,1))/4.*(x1-x2)+(1-pw(j,1))/4.*(x4-x3)
             JACO_12=(1+pw(i,1))/4.*(x1-x4)+(1-pw(i,1))/4.*(x2-x3)
             JACO_21=(1+pw(j,1))/4.*(y1-y2)+(1-pw(j,1))/4.*(y4-y3)
             JACO_22=(1+pw(i,1))/4.*(y1-y4)+(1-pw(i,1))/4.*(y2-y3)
             JACO=JACO_11*JACO_22-JACO_12*JACO_21
             
             kesi=x1*N1+x2*N2+x3*N3+x4*N4
             yita=y1*N1+y2*N2+y3*N3+y4*N4
            
             A=(localx-kesi)**2+(localy-yita)**2+localz**2
             B=pw(i,2)*pw(j,2)*JACO
             
             f13=f13+B*A**(-2.5)*(localx-kesi)*(localy-yita)
        end do    
    end do    

!Calculate f14
    f14=0. 
    do i=1,M
        do j=1,M
             N1=(1+pw(i,1))*(1+pw(j,1))/4.
             N2=(1-pw(i,1))*(1+pw(j,1))/4.
             N3=(1-pw(i,1))*(1-pw(j,1))/4.
             N4=(1+pw(i,1))*(1-pw(j,1))/4.
             JACO_11=(1+pw(j,1))/4.*(x1-x2)+(1-pw(j,1))/4.*(x4-x3)
             JACO_12=(1+pw(i,1))/4.*(x1-x4)+(1-pw(i,1))/4.*(x2-x3)
             JACO_21=(1+pw(j,1))/4.*(y1-y2)+(1-pw(j,1))/4.*(y4-y3)
             JACO_22=(1+pw(i,1))/4.*(y1-y4)+(1-pw(i,1))/4.*(y2-y3)
             JACO=JACO_11*JACO_22-JACO_12*JACO_21
             
             kesi=x1*N1+x2*N2+x3*N3+x4*N4
             yita=y1*N1+y2*N2+y3*N3+y4*N4
            
             A=(localx-kesi)**2+(localy-yita)**2+localz**2
             B=pw(i,2)*pw(j,2)*JACO
             
             f14=f14+B*A**(-2.5)*(localx-kesi)**3
        end do    
    end do  
    
!Calculate f15
    f15=0. 
    do i=1,M
        do j=1,M
             N1=(1+pw(i,1))*(1+pw(j,1))/4.
             N2=(1-pw(i,1))*(1+pw(j,1))/4.
             N3=(1-pw(i,1))*(1-pw(j,1))/4.
             N4=(1+pw(i,1))*(1-pw(j,1))/4.
             JACO_11=(1+pw(j,1))/4.*(x1-x2)+(1-pw(j,1))/4.*(x4-x3)
             JACO_12=(1+pw(i,1))/4.*(x1-x4)+(1-pw(i,1))/4.*(x2-x3)
             JACO_21=(1+pw(j,1))/4.*(y1-y2)+(1-pw(j,1))/4.*(y4-y3)
             JACO_22=(1+pw(i,1))/4.*(y1-y4)+(1-pw(i,1))/4.*(y2-y3)
             JACO=JACO_11*JACO_22-JACO_12*JACO_21
             
             kesi=x1*N1+x2*N2+x3*N3+x4*N4
             yita=y1*N1+y2*N2+y3*N3+y4*N4
            
             A=(localx-kesi)**2+(localy-yita)**2+localz**2
             B=pw(i,2)*pw(j,2)*JACO
             
             f15=f15+B*A**(-2.5)*(localy-yita)**3
        end do    
    end do    

!Calculate f16
    f16=0. 
    do i=1,M
        do j=1,M
             N1=(1+pw(i,1))*(1+pw(j,1))/4.
             N2=(1-pw(i,1))*(1+pw(j,1))/4.
             N3=(1-pw(i,1))*(1-pw(j,1))/4.
             N4=(1+pw(i,1))*(1-pw(j,1))/4.
             JACO_11=(1+pw(j,1))/4.*(x1-x2)+(1-pw(j,1))/4.*(x4-x3)
             JACO_12=(1+pw(i,1))/4.*(x1-x4)+(1-pw(i,1))/4.*(x2-x3)
             JACO_21=(1+pw(j,1))/4.*(y1-y2)+(1-pw(j,1))/4.*(y4-y3)
             JACO_22=(1+pw(i,1))/4.*(y1-y4)+(1-pw(i,1))/4.*(y2-y3)
             JACO=JACO_11*JACO_22-JACO_12*JACO_21
             
             kesi=x1*N1+x2*N2+x3*N3+x4*N4
             yita=y1*N1+y2*N2+y3*N3+y4*N4
            
             A=(localx-kesi)**2+(localy-yita)**2+localz**2
             B=pw(i,2)*pw(j,2)*JACO
             
             f16=f16+B*A**(-2.5)*(localx-kesi)*(localy-yita)**2
        end do    
    end do    

!Calculate f17
    f17=0. 
    do i=1,M
        do j=1,M
             N1=(1+pw(i,1))*(1+pw(j,1))/4.
             N2=(1-pw(i,1))*(1+pw(j,1))/4.
             N3=(1-pw(i,1))*(1-pw(j,1))/4.
             N4=(1+pw(i,1))*(1-pw(j,1))/4.
             JACO_11=(1+pw(j,1))/4.*(x1-x2)+(1-pw(j,1))/4.*(x4-x3)
             JACO_12=(1+pw(i,1))/4.*(x1-x4)+(1-pw(i,1))/4.*(x2-x3)
             JACO_21=(1+pw(j,1))/4.*(y1-y2)+(1-pw(j,1))/4.*(y4-y3)
             JACO_22=(1+pw(i,1))/4.*(y1-y4)+(1-pw(i,1))/4.*(y2-y3)
             JACO=JACO_11*JACO_22-JACO_12*JACO_21
             
             kesi=x1*N1+x2*N2+x3*N3+x4*N4
             yita=y1*N1+y2*N2+y3*N3+y4*N4
            
             A=(localx-kesi)**2+(localy-yita)**2+localz**2
             B=pw(i,2)*pw(j,2)*JACO
             
             f17=f17+B*A**(-2.5)*(localx-kesi)**2*(localy-yita)
        end do    
    end do    

!************************************************************
! Calculate influence coefficients 
!************************************************************
      S(1,1) =-3.0D0 * F14 - ( 1.0D0 - 2.0D0 * V ) * F3
      S(1,2) =-3.0D0 * F17 + ( 1.0D0 - 2.0D0 * V ) * F4
      S(1,3) =-3.0D0 * F11 + ( 1.0D0 - 2.0D0 * V ) * F2
      S(2,1) =-3.0D0 * F16 + ( 1.0D0 - 2.0D0 * V ) * F3
      S(2,2) =-3.0D0 * F15 - ( 1.0D0 - 2.0D0 * V ) * F4
      S(2,3) =-3.0D0 * F12 + ( 1.0D0 - 2.0D0 * V ) * F2
      S(3,1) =-3.0D0 * F9  + ( 1.0D0 - 2.0D0 * V ) * F3
      S(3,2) =-3.0D0 * F10 + ( 1.0D0 - 2.0D0 * V ) * F4
      S(3,3) =-3.0D0 * F8  - ( 1.0D0 - 2.0D0 * V ) * F2
      S(4,1) =-3.0D0 * F17 - ( 1.0D0 - 2.0D0 * V ) * F4
      S(4,2) =-3.0D0 * F16 - ( 1.0D0 - 2.0D0 * V ) * F3
      S(4,3) =-3.0D0 * F13 * localZ
      S(5,1) =-3.0D0 * F13 * localZ
      S(5,2) =-3.0D0 * F12 - ( 1.0D0 - 2.0D0 * V ) * F2
      S(5,3) =-3.0D0 * F10 - ( 1.0D0 - 2.0D0 * V ) * F4
      S(6,1) =-3.0D0 * F11 - ( 1.0D0 - 2.0D0 * V ) * F2
      S(6,2) =-3.0D0 * F13 * localZ
      S(6,3) =-3.0D0 * F9  - ( 1.0D0 - 2.0D0 * V ) * F3
   
      JI(1,1)=JEvect(1,1)*IEvect(1,1)+JEvect(1,2)*IEvect(1,2)+JEvect(1,3)*IEvect(1,3)
      JI(2,1)=JEvect(2,1)*IEvect(1,1)+JEvect(2,2)*IEvect(1,2)+JEvect(2,3)*IEvect(1,3)
      JI(3,1)=JEvect(3,1)*IEvect(1,1)+JEvect(3,2)*IEvect(1,2)+JEvect(3,3)*IEvect(1,3)
      
      JI(1,2)=JEvect(1,1)*IEvect(2,1)+JEvect(1,2)*IEvect(2,2)+JEvect(1,3)*IEvect(2,3)
      JI(2,2)=JEvect(2,1)*IEvect(2,1)+JEvect(2,2)*IEvect(2,2)+JEvect(2,3)*IEvect(2,3)
      JI(3,2)=JEvect(3,1)*IEvect(2,1)+JEvect(3,2)*IEvect(2,2)+JEvect(3,3)*IEvect(2,3) 
      
      JI(1,3)=JEvect(1,1)*IEvect(3,1)+JEvect(1,2)*IEvect(3,2)+JEvect(1,3)*IEvect(3,3)
      JI(2,3)=JEvect(2,1)*IEvect(3,1)+JEvect(2,2)*IEvect(3,2)+JEvect(2,3)*IEvect(3,3)
      JI(3,3)=JEvect(3,1)*IEvect(3,1)+JEvect(3,2)*IEvect(3,2)+JEvect(3,3)*IEvect(3,3)  
!*******************************************************************************************
	  Trac(1)=S(1,1)*JI(1,1)+S(4,1)*JI(2,1)+S(6,1)*JI(3,1)
      Trac(2)=S(4,1)*JI(1,1)+S(2,1)*JI(2,1)+S(5,1)*JI(3,1)
      Trac(3)=S(6,1)*JI(1,1)+S(5,1)*JI(2,1)+S(3,1)*JI(3,1)
      
      COEFF(1,1)=cons*(Trac(1)*JI(1,1)+Trac(2)*JI(2,1)+Trac(3)*JI(3,1))
      COEFF(2,1)=cons*(Trac(1)*JI(1,2)+Trac(2)*JI(2,2)+Trac(3)*JI(3,2))
      COEFF(3,1)=cons*(Trac(1)*JI(1,3)+Trac(2)*JI(2,3)+Trac(3)*JI(3,3))
     
      Trac(1)=S(1,1)*JI(1,2)+S(4,1)*JI(2,2)+S(6,1)*JI(3,2)
      Trac(2)=S(4,1)*JI(1,2)+S(2,1)*JI(2,2)+S(5,1)*JI(3,2)
      Trac(3)=S(6,1)*JI(1,2)+S(5,1)*JI(2,2)+S(3,1)*JI(3,2)
     
      COEFF(4,1)=cons*(Trac(1)*JI(1,2)+Trac(2)*JI(2,2)+Trac(3)*JI(3,2))
      COEFF(5,1)=cons*(Trac(1)*JI(1,3)+Trac(2)*JI(2,3)+Trac(3)*JI(3,3))     
     
      Trac(1)=S(1,1)*JI(1,3)+S(4,1)*JI(2,3)+S(6,1)*JI(3,3)
      Trac(2)=S(4,1)*JI(1,3)+S(2,1)*JI(2,3)+S(5,1)*JI(3,3)
      Trac(3)=S(6,1)*JI(1,3)+S(5,1)*JI(2,3)+S(3,1)*JI(3,3)
     
      COEFF(6,1)=cons*(Trac(1)*JI(1,3)+Trac(2)*JI(2,3)+Trac(3)*JI(3,3))  
     
!************************************************************************************    
      Trac(1)=S(1,2)*JI(1,1)+S(4,2)*JI(2,1)+S(6,2)*JI(3,1)
      Trac(2)=S(4,2)*JI(1,1)+S(2,2)*JI(2,1)+S(5,2)*JI(3,1)
      Trac(3)=S(6,2)*JI(1,1)+S(5,2)*JI(2,1)+S(3,2)*JI(3,1)
     
      COEFF(1,2)=cons*(Trac(1)*JI(1,1)+Trac(2)*JI(2,1)+Trac(3)*JI(3,1))
      COEFF(2,2)=cons*(Trac(1)*JI(1,2)+Trac(2)*JI(2,2)+Trac(3)*JI(3,2))
      COEFF(3,2)=cons*(Trac(1)*JI(1,3)+Trac(2)*JI(2,3)+Trac(3)*JI(3,3))
     
      Trac(1)=S(1,2)*JI(1,2)+S(4,2)*JI(2,2)+S(6,2)*JI(3,2)
      Trac(2)=S(4,2)*JI(1,2)+S(2,2)*JI(2,2)+S(5,2)*JI(3,2)
      Trac(3)=S(6,2)*JI(1,2)+S(5,2)*JI(2,2)+S(3,2)*JI(3,2)
     
      COEFF(4,2)=cons*(Trac(1)*JI(1,2)+Trac(2)*JI(2,2)+Trac(3)*JI(3,2))
      COEFF(5,2)=cons*(Trac(1)*JI(1,3)+Trac(2)*JI(2,3)+Trac(3)*JI(3,3))
          
      Trac(1)=S(1,2)*JI(1,3)+S(4,2)*JI(2,3)+S(6,2)*JI(3,3)
      Trac(2)=S(4,2)*JI(1,3)+S(2,2)*JI(2,3)+S(5,2)*JI(3,3)
      Trac(3)=S(6,2)*JI(1,3)+S(5,2)*JI(2,3)+S(3,2)*JI(3,3)
     
      COEFF(6,2)=cons*(Trac(1)*JI(1,3)+Trac(2)*JI(2,3)+Trac(3)*JI(3,3))
!*************************************************************************************
      Trac(1)=S(1,3)*JI(1,1)+S(4,3)*JI(2,1)+S(6,3)*JI(3,1)
      Trac(2)=S(4,3)*JI(1,1)+S(2,3)*JI(2,1)+S(5,3)*JI(3,1)
      Trac(3)=S(6,3)*JI(1,1)+S(5,3)*JI(2,1)+S(3,3)*JI(3,1)
     
      COEFF(1,3)=cons*(Trac(1)*JI(1,1)+Trac(2)*JI(2,1)+Trac(3)*JI(3,1))
	  COEFF(2,3)=cons*(Trac(1)*JI(1,2)+Trac(2)*JI(2,2)+Trac(3)*JI(3,2))
      COEFF(3,3)=cons*(Trac(1)*JI(1,3)+Trac(2)*JI(2,3)+Trac(3)*JI(3,3))
     
      Trac(1)=S(1,3)*JI(1,2)+S(4,3)*JI(2,2)+S(6,3)*JI(3,2)
      Trac(2)=S(4,3)*JI(1,2)+S(2,3)*JI(2,2)+S(5,3)*JI(3,2)
      Trac(3)=S(6,3)*JI(1,2)+S(5,3)*JI(2,2)+S(3,3)*JI(3,2)
     
      COEFF(4,3)=cons*(Trac(1)*JI(1,2)+Trac(2)*JI(2,2)+Trac(3)*JI(3,2))
      COEFF(5,3)=cons*(Trac(1)*JI(1,3)+Trac(2)*JI(2,3)+Trac(3)*JI(3,3))
     
      Trac(1)=S(1,3)*JI(1,3)+S(4,3)*JI(2,3)+S(6,3)*JI(3,3)
      Trac(2)=S(4,3)*JI(1,3)+S(2,3)*JI(2,3)+S(5,3)*JI(3,3)
      Trac(3)=S(6,3)*JI(1,3)+S(5,3)*JI(2,3)+S(3,3)*JI(3,3)
     
      COEFF(6,3)=cons*(Trac(1)*JI(1,3)+Trac(2)*JI(2,3)+Trac(3)*JI(3,3))

end subroutine stress_regular_fs3d_qua_num_part3


subroutine stress_singular_fs3d_qua_num (xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4,EE,v,M,IEvect,COEFF)
implicit none 

real*8::EE,v,GG,cons,l1,l2,l3,p,area,x(5),y(5),z(5),xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4,xx(4),yy(4),zz(4)
real*8::JEvect(3,3),IEvect(3,3),COEFF(3,3),A,B,C,D,E,F,G,H,I,J,k,yip,xm,ym,r,localz,gpx,gpy,gpz
real*8::f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,S(6,3),JI(3,3),TRAC(3)
integer::ii,jj,M
real*8,allocatable::pw1(:,:)

GG=EE/2./(1+v)
cons=1./8./3.141592653/(1-v)
localz=0.
    
allocate (pw1(M,2))

do i=1,3
	do j=1,3
		JEvect(i,j)=IEvect(i,j)
	end do 
end do 
    
gpx=(xx1+xx2+xx3+xx4)/4.
gpy=(yy1+yy2+yy3+yy4)/4.
gpz=(zz1+zz2+zz3+zz4)/4.

xx(1)=xx1-gpx
yy(1)=yy1-gpy
zz(1)=zz1-gpz

xx(2)=xx2-gpx
yy(2)=yy2-gpy
zz(2)=zz2-gpz

xx(3)=xx3-gpx
yy(3)=yy3-gpy
zz(3)=zz3-gpz

xx(4)=xx4-gpx
yy(4)=yy4-gpy
zz(4)=zz4-gpz

x(1)=xx(1)*IEvect(1,1)+yy(1)*IEvect(1,2)+zz(1)*IEvect(1,3)
y(1)=xx(1)*IEvect(2,1)+yy(1)*IEvect(2,2)+zz(1)*IEvect(2,3)
z(1)=xx(1)*IEvect(3,1)+yy(1)*IEvect(3,2)+zz(1)*IEvect(3,3)

x(2)=xx(2)*IEvect(1,1)+yy(2)*IEvect(1,2)+zz(2)*IEvect(1,3)
y(2)=xx(2)*IEvect(2,1)+yy(2)*IEvect(2,2)+zz(2)*IEvect(2,3)
z(2)=xx(2)*IEvect(3,1)+yy(2)*IEvect(3,2)+zz(2)*IEvect(3,3)

x(3)=xx(3)*IEvect(1,1)+yy(3)*IEvect(1,2)+zz(3)*IEvect(1,3)
y(3)=xx(3)*IEvect(2,1)+yy(3)*IEvect(2,2)+zz(3)*IEvect(2,3)
z(3)=xx(3)*IEvect(3,1)+yy(3)*IEvect(3,2)+zz(3)*IEvect(3,3)

x(4)=xx(4)*IEvect(1,1)+yy(4)*IEvect(1,2)+zz(4)*IEvect(1,3)
y(4)=xx(4)*IEvect(2,1)+yy(4)*IEvect(2,2)+zz(4)*IEvect(2,3)
z(4)=xx(4)*IEvect(3,1)+yy(4)*IEvect(3,2)+zz(4)*IEvect(3,3)

x(5)=x(1)
y(5)=y(1)
z(5)=z(1)

call Gauss_ponits_weights(m,pw1)

f1=0.
f2=0.
f3=0.
f4=0.
f5=0.
f6=0.
f7=0.
f8=0.
f9=0.
f10=0.
f11=0.
f12=0.
f13=0.
f14=0.
f15=0.
f16=0.
f17=0.

do 100 ii=1,4
    
    l1=sqrt(x(ii)**2+y(ii)**2+z(ii)**2)
    l2=sqrt((x(ii+1)-x(ii))**2+(y(ii+1)-y(ii))**2+(z(ii+1)-z(ii))**2)
    l3=sqrt(x(ii+1)**2+y(ii+1)**2+z(ii+1)**2)
    p=(l1+l2+l3)/2.
    area=sqrt(p*(p-l1)*(p-l2)*(p-l3))
       
    A=0.
    B=0.
    C=0.
    D=0.
    E=0.
    F=0.
    G=0.
    H=0.
    I=0.
    J=0.
    
    do jj=1,m
        
        yip=0.5*(pw1(jj,1)+1) 
        xm=(1-yip)*x(ii)+yip*x(ii+1)
        ym=(1-yip)*y(ii)+yip*y(ii+1)
        r=sqrt(xm**2+ym**2) 
        
 !**************************************************************************************        
 ! For self-influence coefficients, only sencond-order partial derivatives are needed.  
 !**************************************************************************************
        A=A+pw1(jj,2)*area/r
        B=B-pw1(jj,2)*area*xm*log(r)/r**3
        C=C-pw1(jj,2)*area*ym*log(r)/r**3
        D=D+pw1(jj,2)*area*xm**2/r**3
        E=E+pw1(jj,2)*area*ym**2/r**3
        F=F+pw1(jj,2)*area*xm*ym/r**3
        G=G-pw1(jj,2)*area*xm**3*log(r)/r**5
        H=H-pw1(jj,2)*area*ym**3*log(r)/r**5
        I=I-pw1(jj,2)*area*ym**2*xm*log(r)/r**5
        J=J-pw1(jj,2)*area*xm**2*ym*log(r)/r**5
        !K=K+pw1(jj,2)*area*xm*ym/r**5
        !write(*,*)jj,'A=',A,'B=',B,'C=',C,'G=',G
    end do
    
       f1=f1+A
       f2=0
       f3=f3+B
       f4=f4+C
       f5=f5+D
       f6=f6+E
       f7=f7+F
       f8=0.
       f9=0.
       f10=0.
       f11=0.
       f12=0.
       f13=0.
       f14=f14+G
       f15=f15+H
       f16=f16+I
       f17=f17+J
 
100 continue  
 
!*******************************************************************       
      S(1,1) =-3.0D0 * F14 - ( 1.0D0 - 2.0D0 * V ) * F3
      S(1,2) =-3.0D0 * F17 + ( 1.0D0 - 2.0D0 * V ) * F4
      S(1,3) =-3.0D0 * F11 + ( 1.0D0 - 2.0D0 * V ) * F2
      S(2,1) =-3.0D0 * F16 + ( 1.0D0 - 2.0D0 * V ) * F3
      S(2,2) =-3.0D0 * F15 - ( 1.0D0 - 2.0D0 * V ) * F4
      S(2,3) =-3.0D0 * F12 + ( 1.0D0 - 2.0D0 * V ) * F2
      S(3,1) =-3.0D0 * F9  + ( 1.0D0 - 2.0D0 * V ) * F3
      S(3,2) =-3.0D0 * F10 + ( 1.0D0 - 2.0D0 * V ) * F4
      S(3,3) =-3.0D0 * F8  - ( 1.0D0 - 2.0D0 * V ) * F2
      S(4,1) =-3.0D0 * F17 - ( 1.0D0 - 2.0D0 * V ) * F4
      S(4,2) =-3.0D0 * F16 - ( 1.0D0 - 2.0D0 * V ) * F3
      S(4,3) =-3.0D0 * F13 * localZ
      S(5,1) =-3.0D0 * F13 * localZ
      S(5,2) =-3.0D0 * F12 - ( 1.0D0 - 2.0D0 * V ) * F2
      S(5,3) =-3.0D0 * F10 - ( 1.0D0 - 2.0D0 * V ) * F4
      S(6,1) =-3.0D0 * F11 - ( 1.0D0 - 2.0D0 * V ) * F2
      S(6,2) =-3.0D0 * F13 * localZ
      S(6,3) =-3.0D0 * F9  - ( 1.0D0 - 2.0D0 * V ) * F3  
     
      JI(1,1)=JEvect(1,1)*IEvect(1,1)+JEvect(1,2)*IEvect(1,2)+JEvect(1,3)*IEvect(1,3)
      JI(2,1)=JEvect(2,1)*IEvect(1,1)+JEvect(2,2)*IEvect(1,2)+JEvect(2,3)*IEvect(1,3)
      JI(3,1)=JEvect(3,1)*IEvect(1,1)+JEvect(3,2)*IEvect(1,2)+JEvect(3,3)*IEvect(1,3)
      
      JI(1,2)=JEvect(1,1)*IEvect(2,1)+JEvect(1,2)*IEvect(2,2)+JEvect(1,3)*IEvect(2,3)
      JI(2,2)=JEvect(2,1)*IEvect(2,1)+JEvect(2,2)*IEvect(2,2)+JEvect(2,3)*IEvect(2,3)
      JI(3,2)=JEvect(3,1)*IEvect(2,1)+JEvect(3,2)*IEvect(2,2)+JEvect(3,3)*IEvect(2,3) 
      
      JI(1,3)=JEvect(1,1)*IEvect(3,1)+JEvect(1,2)*IEvect(3,2)+JEvect(1,3)*IEvect(3,3)
      JI(2,3)=JEvect(2,1)*IEvect(3,1)+JEvect(2,2)*IEvect(3,2)+JEvect(2,3)*IEvect(3,3)
      JI(3,3)=JEvect(3,1)*IEvect(3,1)+JEvect(3,2)*IEvect(3,2)+JEvect(3,3)*IEvect(3,3)  

      Trac(1)=S(1,1)*JI(1,3)+S(4,1)*JI(2,3)+S(6,1)*JI(3,3)
      Trac(2)=S(4,1)*JI(1,3)+S(2,1)*JI(2,3)+S(5,1)*JI(3,3)
      Trac(3)=S(6,1)*JI(1,3)+S(5,1)*JI(2,3)+S(3,1)*JI(3,3)
     
      COEFF(1,1)=cons*(Trac(1)*JI(1,1)+Trac(2)*JI(2,1)+Trac(3)*JI(3,1))
      COEFF(2,1)=cons*(Trac(1)*JI(1,2)+Trac(2)*JI(2,2)+Trac(3)*JI(3,2))
      COEFF(3,1)=cons*(Trac(1)*JI(1,3)+Trac(2)*JI(2,3)+Trac(3)*JI(3,3))
     
      Trac(1)=S(1,2)*JI(1,3)+S(4,2)*JI(2,3)+S(6,2)*JI(3,3)
      Trac(2)=S(4,2)*JI(1,3)+S(2,2)*JI(2,3)+S(5,2)*JI(3,3)
      Trac(3)=S(6,2)*JI(1,3)+S(5,2)*JI(2,3)+S(3,2)*JI(3,3)
     
      COEFF(1,2)=cons*(Trac(1)*JI(1,1)+Trac(2)*JI(2,1)+Trac(3)*JI(3,1))
	  COEFF(2,2)=cons*(Trac(1)*JI(1,2)+Trac(2)*JI(2,2)+Trac(3)*JI(3,2))
      COEFF(3,2)=cons*(Trac(1)*JI(1,3)+Trac(2)*JI(2,3)+Trac(3)*JI(3,3))
     
      Trac(1)=S(1,3)*JI(1,3)+S(4,3)*JI(2,3)+S(6,3)*JI(3,3)
      Trac(2)=S(4,3)*JI(1,3)+S(2,3)*JI(2,3)+S(5,3)*JI(3,3)
      Trac(3)=S(6,3)*JI(1,3)+S(5,3)*JI(2,3)+S(3,3)*JI(3,3)
     
      COEFF(1,3)=cons*(Trac(1)*JI(1,1)+Trac(2)*JI(2,1)+Trac(3)*JI(3,1))
      COEFF(2,3)=cons*(Trac(1)*JI(1,2)+Trac(2)*JI(2,2)+Trac(3)*JI(3,2))
      COEFF(3,3)=cons*(Trac(1)*JI(1,3)+Trac(2)*JI(2,3)+Trac(3)*JI(3,3))
     
      COEFF(1,1)=0.5
      COEFF(2,2)=0.5
      COEFF(3,3)=0.5
     
end subroutine stress_singular_fs3d_qua_num

  
subroutine Gauss_ponits_weights (M,pw)
implicit none 
integer,intent(in)::M
real*8,intent(out)::pw(M,2)

select case (M)
case (:1)
    WRITE (*,*)'Error, The muber of Gauss points must be >=2 and <=30.'      
case (2)
    pw(1,1)=-0.577350269189625
    Pw(1,2)=1.0
    pw(2,1)=0.577350269189625
    Pw(2,2)=1.0 
case (3)
    pw(1,1)=-0.774596669241483
    pw(1,2)=0.555555555555555
    pw(2,1)=0.000000000000000	
    Pw(2,2)=0.888888888888888
    pw(3,1)=0.774596669241483
    pw(3,2)=0.555555555555555
case (4)
    pw(1,1)=-0.861136311594052	
    pw(1,2)= 0.347854845137453
    pw(2,1)=-0.339981043584856
	pw(2,2)= 0.652145154862546
    pw(3,1)= 0.339981043584856	
    pw(3,2)= 0.652145154862546
    pw(4,1)= 0.861136311594052	
    pw(4,2)= 0.347854845137453
case (5)
    pw(1,1)=-0.906179845938664	
    pw(1,2)=0.236926885056189
    pw(2,1)=-0.538469310105683	
    pw(2,2)=0.478628670499366
    pw(3,1)=0.000000000000000	
    pw(3,2)=0.568888888888888
    pw(4,1)=0.538469310105683	
    pw(4,2)=0.478628670499366
    pw(5,1)=0.906179845938664	
    pw(5,2)=0.236926885056189		
case (6)
    pw(1,1)=-0.932469514203152
	pw(1,2)= 0.171324492379170
    pw(2,1)=-0.661209386466264
	pw(2,2)= 0.360761573048138
    pw(3,1)=-0.238619186083196
	pw(3,2)= 0.467913934572691
    pw(4,1)= 0.238619186083196
	pw(4,2)= 0.467913934572691
    pw(5,1)= 0.661209386466264
	pw(5,2)= 0.360761573048138
    pw(6,1)= 0.932469514203152
	pw(6,2)= 0.171324492379170
case(7)
    pw(1,1)=-0.949107912342758
	pw(1,2)= 0.129484966168869
    pw(2,1)=-0.741531185599394	
    pw(2,2)= 0.279705391489276
    pw(3,1)=-0.405845151377397
	pw(3,2)= 0.381830050505118
    pw(4,1)= 0.000000000000000	
    pw(4,2)= 0.417959183673469
    pw(5,1)= 0.405845151377397
	pw(5,2)= 0.381830050505118
    pw(6,1)= 0.741531185599394	
    pw(6,2)= 0.279705391489276
    pw(7,1)= 0.949107912342758
	pw(7,2)= 0.129484966168869
case(8)
    pw(1,1)=-0.960289856497536	
    pw(1,2)= 0.101228536290376
    pw(2,1)=-0.796666477413626	
    pw(2,2)= 0.222381034453374
    pw(3,1)=-0.525532409916329	
    pw(3,2)= 0.313706645877887
    pw(4,1)=-0.183434642495649	
    pw(4,2)= 0.362683783378362
    pw(5,1)= 0.183434642495649	
    pw(5,2)= 0.362683783378362
    pw(6,1)= 0.525532409916329	
    pw(6,2)= 0.313706645877887
    pw(7,1)= 0.796666477413626	
    pw(7,2)= 0.222381034453374
    pw(8,1)= 0.960289856497536	
    pw(8,2)= 0.101228536290376
case(9)
    pw(1,1)=-0.968160239507626	
    pw(1,2)= 0.081274388361574
    pw(2,1)=-0.836031107326635	
    pw(2,2)= 0.180648160694857
    pw(3,1)=-0.613371432700590	
    pw(3,2)= 0.260610696402935
    pw(4,1)=-0.324253423403808
	pw(4,2)= 0.312347077040002
    pw(5,1)= 0.000000000000000	
    pw(5,2)= 0.330239355001259
    pw(6,1)= 0.324253423403808
	pw(6,2)= 0.312347077040002
    pw(7,1)= 0.613371432700590	
    pw(7,2)= 0.260610696402935
    pw(8,1)= 0.836031107326635	
    pw(8,2)= 0.180648160694857
    pw(9,1)= 0.968160239507626	
    pw(9,2)= 0.081274388361574
case(10)
    pw(1,1)=-0.973906528517171
	pw(1,2)= 0.066671344308688
    pw(2,1)=-0.865063366688984
	pw(2,2)= 0.149451349150580
    pw(3,1)=-0.679409568299024
	pw(3,2)= 0.219086362515982
    pw(4,1)=-0.433395394129247
	pw(4,2)= 0.269266719309996
    pw(5,1)=-0.148874338981631
	pw(5,2)= 0.295524224714752
    pw(6,1)= 0.148874338981631
	pw(6,2)= 0.295524224714752
    pw(7,1)= 0.433395394129247
	pw(7,2)= 0.269266719309996
    pw(8,1)= 0.679409568299024
	pw(8,2)= 0.219086362515982
    pw(9,1)= 0.865063366688984
	pw(9,2)= 0.149451349150580
    pw(10,1)= 0.973906528517171
	pw(10,2)= 0.066671344308688
case(11)
    pw(1,1)=-0.978228658146057
	pw(1,2)= 0.055668567116174
    pw(2,1)=-0.887062599768095
	pw(2,2)= 0.125580369464904
    pw(3,1)=-0.730152005574049	
    pw(3,2)= 0.186290210927734
    pw(4,1)=-0.519096129206811
	pw(4,2)= 0.233193764591990
    pw(5,1)=-0.269543155952345
	pw(5,2)= 0.262804544510246
    pw(6,1)= 0.000000000000000
	pw(6,2)= 0.272925086777900
    pw(7,1)= 0.269543155952345
	pw(7,2)= 0.262804544510246
    pw(8,1)= 0.519096129206811
	pw(8,2)= 0.233193764591990
    pw(9,1)= 0.730152005574049	
    pw(9,2)= 0.186290210927734
    pw(10,1)= 0.887062599768095
	pw(10,2)= 0.125580369464904
    pw(11,1)= 0.978228658146057
	pw(11,2)= 0.055668567116174
case(12)
    pw(1,1)=-0.981560634246719
	pw(1,2)= 0.047175336386512
    pw(2,1)=-0.904117256370474	
    pw(2,2)= 0.106939325995318
    pw(3,1)=-0.769902674194304
	pw(3,2)= 0.160078328543346
    pw(4,1)=-0.587317954286617
	pw(4,2)= 0.203167426723065
    pw(5,1)=-0.367831498998180	
    pw(5,2)= 0.233492536538354
    pw(6,1)=-0.125233408511468
	pw(6,2)= 0.249147045813402
    pw(7,1)= 0.125233408511468
	pw(7,2)= 0.249147045813402
    pw(8,1)= 0.367831498998180	
    pw(8,2)= 0.233492536538354
    pw(9,1)= 0.587317954286617
	pw(9,2)= 0.203167426723065
    pw(10,1)= 0.769902674194304
	pw(10,2)= 0.160078328543346
    pw(11,1)= 0.904117256370474	
    pw(11,2)= 0.106939325995318
    pw(12,1)= 0.981560634246719
	pw(12,2)= 0.047175336386512
case(13)
    pw(1,1)=-0.984183054718588	
    pw(1,2)= 0.040484004765316
    pw(2,1)=-0.917598399222977	
    pw(2,2)= 0.092121499837729
    pw(3,1)=-0.801578090733309
	pw(3,2)= 0.138873510219787
    pw(4,1)=-0.642349339440340	
    pw(4,2)= 0.178145980761945    
    pw(5,1)=-0.448492751036446	
    pw(5,2)= 0.207816047536888    
    pw(6,1)=-0.230458315955134
	pw(6,2)= 0.226283180262897    
    pw(7,1)= 0.000000000000000
	pw(7,2)= 0.232551553230873
    pw(8,1)= 0.230458315955134
	pw(8,2)= 0.226283180262897
    pw(9,1)= 0.448492751036446	
    pw(9,2)= 0.207816047536888
    pw(10,1)= 0.642349339440340	
    pw(10,2)= 0.178145980761945
    pw(11,1)= 0.801578090733309
	pw(11,2)= 0.138873510219787
    pw(12,1)= 0.917598399222977	
    pw(12,2)= 0.092121499837729
    pw(13,1)= 0.984183054718588	
    pw(13,2)= 0.040484004765316
case(14)
    pw(1,1)=-0.986283808696812	
    pw(1,2)= 0.035119460331752
    pw(2,1)=-0.928434883663573
	pw(2,2)= 0.080158087159760
    pw(3,1)=-0.827201315069765	
    pw(3,2)= 0.121518570687903    
    pw(4,1)=-0.687292904811685
	pw(4,2)= 0.157203167158193
    pw(5,1)=-0.515248636358154
	pw(5,2)= 0.185538397477937
    pw(6,1)=-0.319112368927889
	pw(6,2)= 0.205198463721295     
    pw(7,1)=-0.108054948707343	
    pw(7,2)= 0.215263853463157
    pw(8,1)= 0.108054948707343	
    pw(8,2)= 0.215263853463157
    pw(9,1)= 0.319112368927889
	pw(9,2)= 0.205198463721295
    pw(10,1)= 0.515248636358154
	pw(10,2)= 0.185538397477937
    pw(11,1)= 0.687292904811685
	pw(11,2)= 0.157203167158193
    pw(12,1)= 0.827201315069765	
    pw(12,2)= 0.121518570687903
    pw(13,1)= 0.928434883663573
	pw(13,2)= 0.080158087159760
    pw(14,1)= 0.986283808696812	
    pw(14,2)= 0.035119460331752
case(15)
    pw(1,1)=-0.987992518020485
    pw(1,2)=0.030753241996117
    pw(2,1)=-0.937273392400706	
    pw(2,2)=0.070366047488108
    pw(3,1)=-0.848206583410427	
    pw(3,2)=0.107159220467171
    pw(4,1)=-0.724417731360170	
    pw(4,2)=0.139570677926154
    pw(5,1)=-0.570972172608538	
    pw(5,2)=0.166269205816993
    pw(6,1)=-0.394151347077563	
    pw(6,2)=0.186161000015562
    pw(7,1)=-0.201194093997434	
    pw(7,2)=0.198431485327111
    pw(8,1)=0.000000000000000	
    pw(8,2)=0.202578241925561
    pw(9,1)=0.201194093997434	
    pw(9,2)=0.198431485327111
    pw(10,1)=0.394151347077563	
    pw(10,2)=0.186161000015562
    pw(11,1)=0.570972172608538	
    pw(11,2)=0.166269205816993
    pw(12,1)=0.724417731360170	
    pw(12,2)=0.139570677926154
    pw(13,1)=0.848206583410427	
    pw(13,2)=0.107159220467171
    pw(14,1)=0.937273392400706	
    pw(14,2)=0.070366047488108
    pw(15,1)=0.987992518020485	
    pw(15,2)=0.030753241996117
case (16)
    pw(1,1)=-0.989400934991649
	pw(1,2)= 0.027152459411754
    pw(2,1)=-0.944575023073232
	pw(2,2)= 0.062253523938648
    pw(3,1)=-0.865631202387831	
    pw(3,2)= 0.095158511682493
    pw(4,1)=-0.755404408355003
	pw(4,2)= 0.124628971255533
    pw(5,1)=-0.617876244402643
	pw(5,2)= 0.149595988816576
    pw(6,1)=-0.458016777657227
    pw(6,2)= 0.169156519395002
    pw(7,1)=-0.281603550779258	
    pw(7,2)= 0.182603415044923
    pw(8,1)=-0.095012509837637	
    pw(8,2)=0.189450610455068
    pw(9,1)=0.095012509837637	
    pw(9,2)=0.189450610455068
    pw(10,1)=0.281603550779258	
    pw(10,2)=0.182603415044923
    pw(11,1)=0.458016777657227	
    pw(11,2)=0.169156519395002
    pw(12,1)=0.617876244402643	
    pw(12,2)=0.149595988816576
    pw(13,1)=0.755404408355003	
    pw(13,2)=0.124628971255533
    pw(14,1)=0.865631202387831	
    pw(14,2)=0.095158511682493
    pw(15,1)=0.944575023073232	
    pw(15,2)=0.062253523938648
    pw(16,1)=0.989400934991649	
    pw(16,2)=0.027152459411754
case (17)
    pw(1,1)=-0.990575475314417	
    pw(1,2)=0.024148302868548
    pw(2,1)=-0.950675521768767	
    pw(2,2)=0.055459529373987
    pw(3,1)=-0.880239153726985	
    pw(3,2)=0.085036148317179
    pw(4,1)=-0.781514003896801	
    pw(4,2)=0.111883847193404
    pw(5,1)=-0.657671159216690	
    pw(5,2)=0.135136368468525
    pw(6,1)=-0.512690537086476	
    pw(6,2)=0.154045761076810
    pw(7,1)=-0.351231763453876	
    pw(7,2)=0.168004102156450
    pw(8,1)=-0.178484181495847
    pw(8,2)=0.176562705366992
    pw(9,1)=0.000000000000000	
    pw(9,2)=0.179446470356206
    pw(10,1)=0.178484181495847
    pw(10,2)=0.176562705366992
    pw(11,1)=0.351231763453876	
    pw(11,2)=0.168004102156450
    pw(12,1)=0.512690537086476	
    pw(12,2)=0.154045761076810
    pw(13,1)=0.657671159216690	
    pw(13,2)=0.135136368468525
    pw(14,1)=0.781514003896801	
    pw(14,2)=0.111883847193404
    pw(15,1)=0.880239153726985	
    pw(15,2)=0.085036148317179
    pw(16,1)=0.950675521768767
    pw(16,2)=0.055459529373987
    pw(17,1)=0.990575475314417	
    pw(17,2)=0.024148302868548
case (18)
    pw(1,1)=-0.991565168420930
    pw(1,2)= 0.021616013526483
    pw(2,1)=-0.955823949571397
    pw(2,2)= 0.049714548894970
    pw(3,1)=-0.892602466497555
    pw(3,2)= 0.076425730254889
    pw(4,1)=-0.803704958972523	
    pw(4,2)= 0.100942044106287
    pw(5,1)=-0.691687043060353	
    pw(5,2)= 0.122555206711478
    pw(6,1)=-0.559770831073947	
    pw(6,2)= 0.140642914670650
    pw(7,1)=-0.411751161462842
    pw(7,2)= 0.154684675126265
    pw(8,1)=-0.251886225691505
    pw(8,2)=0.164276483745832
    pw(9,1)=-0.084775013041735
    pw(9,2)=0.169142382963143
    pw(10,1)=0.084775013041735
    pw(10,2)=0.169142382963143
    pw(11,1)=0.251886225691505	
    pw(11,2)=0.164276483745832
    pw(12,1)=0.411751161462842	
    pw(12,2)=0.154684675126265
    pw(13,1)=0.559770831073947	
    pw(13,2)=0.140642914670650
    pw(14,1)=0.691687043060353	
    pw(14,2)=0.122555206711478
    pw(15,1)=0.803704958972523	
    pw(15,2)=0.100942044106287
    pw(16,1)=0.892602466497555	
    pw(16,2)=0.076425730254889
    pw(17,1)=0.955823949571397	
    pw(17,2)=0.049714548894970
    pw(18,1)=0.991565168420930	
    pw(18,2)=0.021616013526483
case(19)
    pw(1,1)=-0.992406843843584
    pw(1,2)= 0.019461788229727
    pw(2,1)=-0.960208152134830
    pw(2,2)= 0.044814226765700
    pw(3,1)=-0.903155903614817
    pw(3,2)= 0.069044542737641
    pw(4,1)=-0.822714656537142
    pw(4,2)= 0.091490021622450
    pw(5,1)=-0.720966177335229	
    pw(5,2)= 0.111566645547334
    pw(6,1)=-0.600545304661681	
    pw(6,2)= 0.128753962539336
    pw(7,1)=-0.464570741375960	
    pw(7,2)= 0.142606702173606
    pw(8,1)=-0.316564099963629	
    pw(8,2)= 0.152766042065859
    pw(9,1)=-0.160358645640225	
    pw(9,2)= 0.158968843393954
    pw(10,1)=0.000000000000000	
    pw(10,2)=0.161054449848783
    pw(11,1)=0.160358645640225
    pw(11,2)=0.158968843393954
    pw(12,1)=0.316564099963629	
    pw(12,2)=0.152766042065859
    pw(13,1)=0.464570741375960	
    pw(13,2)=0.142606702173606
    pw(14,1)=0.600545304661681	
    pw(14,2)=0.128753962539336
    pw(15,1)=0.720966177335229
    pw(15,2)=0.111566645547334
    pw(16,1)=0.822714656537142	
    pw(16,2)=0.091490021622450
    pw(17,1)=0.903155903614817
    pw(17,2)=0.069044542737641
    pw(18,1)=0.960208152134830	
    pw(18,2)=0.044814226765700
    pw(19,1)=0.992406843843584	
    pw(19,2)=0.019461788229727
case(20)
    pw(1,1)=-0.993128599185094	
    pw(1,2)= 0.017614007139152
    pw(2,1)=-0.963971927277913	
    pw(2,2)= 0.040601429800387
    pw(3,1)=-0.912234428251325	
    pw(3,2)= 0.062672048334109
    pw(4,1)=-0.839116971822218	
    pw(4,2)= 0.083276741576705
    pw(5,1)=-0.746331906460150	
    pw(5,2)= 0.101930119817240
    pw(6,1)=-0.636053680726515	
    pw(6,2)= 0.118194531961518
    pw(7,1)=-0.510867001950827	
    pw(7,2)= 0.131688638449176
    pw(8,1)=-0.373706088715419	
    pw(8,2)= 0.142096109318382
    pw(9,1)=-0.227785851141645	
    pw(9,2)= 0.149172986472603
    pw(10,1)=-0.076526521133497	
    pw(10,2)=0.152753387130725
    pw(11,1)=0.076526521133497	
    pw(11,2)=0.152753387130725
    pw(12,1)=0.227785851141645	
    pw(12,2)=0.149172986472603
    pw(13,1)=0.373706088715419	
    pw(13,2)=0.142096109318382
    pw(14,1)=0.510867001950827	
    pw(14,2)=0.131688638449176
    pw(15,1)=0.636053680726515	
    pw(15,2)=0.118194531961518
    pw(16,1)=0.746331906460150	
    pw(16,2)=0.101930119817240
    pw(17,1)=0.839116971822218	
    pw(17,2)=0.083276741576705
    pw(18,1)=0.912234428251325	
    pw(18,2)=0.062672048334109
    pw(19,1)=0.963971927277913	
    pw(19,2)=0.040601429800387
    pw(20,1)=0.993128599185094	
    pw(20,2)=0.017614007139152
case(21)
    pw(1,1)=-0.993752170620389
    pw(1,2)=0.016017228257774
    pw(2,1)=-0.967226838566306
    pw(2,2)=0.036953789770853
    pw(3,1)=-0.920099334150400
    pw(3,2)=0.057134425426857
    pw(4,1)=-0.853363364583317
    pw(4,2)=0.076100113628379
    pw(5,1)=-0.768439963475677
    pw(5,2)=0.093444423456034
    pw(6,1)=-0.667138804197412
    pw(6,2)=0.108797299167148
    pw(7,1)=-0.551618835887219
    pw(7,2)=0.121831416053728
    pw(8,1)=-0.424342120207438
    pw(8,2)=0.132268938633337
    pw(9,1)=-0.288021316802401
    pw(9,2)=0.139887394791073
    pw(10,1)=-0.145561854160895	
    pw(10,2)=0.144524403989970
    pw(11,1)=0.000000000000000
    pw(11,2)=0.146081133649690
    pw(12,1)=0.145561854160895
    pw(12,2)=0.144524403989970
    pw(13,1)=0.288021316802401
    pw(13,2)=0.139887394791073
    pw(14,1)=0.424342120207438
    pw(14,2)=0.132268938633337
    pw(15,1)=0.551618835887219	
    pw(15,2)=0.121831416053728
    pw(16,1)=0.667138804197412
    pw(16,2)=0.108797299167148
    pw(17,1)=0.768439963475677	
    pw(17,2)=0.093444423456034
    pw(18,1)=0.853363364583317	
    pw(18,2)=0.076100113628379
    pw(19,1)=0.920099334150400	
    pw(19,2)=0.057134425426857
    pw(20,1)=0.967226838566306	
    pw(20,2)=0.036953789770853
    pw(21,1)=0.993752170620389	
    pw(21,2)=0.016017228257774
case(22)
    pw(1,1)=-0.994294585482399	
    pw(1,2)=0.014627995298272
    pw(2,1)=-0.970060497835428	
    pw(2,2)=0.033774901584814
    pw(3,1)=-0.926956772187174	
    pw(3,2)=0.052293335152683
    pw(4,1)=-0.865812577720300	
    pw(4,2)=0.069796468424521
    pw(5,1)=-0.787816805979208	
    pw(5,2)=0.085941606217068
    pw(6,1)=-0.694487263186682	
    pw(6,2)=0.100414144442881
    pw(7,1)=-0.587640403506911	
    pw(7,2)=0.112932296080539
    pw(8,1)=-0.469355837986757	
    pw(8,2)=0.123252376810512
    pw(9,1)=-0.341935820892084	
    pw(9,2)=0.131173504787062
    pw(10,1)=-0.207860426688221	
    pw(10,2)=0.136541498346015
    pw(11,1)=-0.069739273319722	
    pw(11,2)=0.139251872855632
    pw(12,1)=0.069739273319722	
    pw(12,2)=0.139251872855632
    pw(13,1)=0.207860426688221	
    pw(13,2)=0.136541498346015
    pw(14,1)=0.341935820892084	
    pw(14,2)=0.131173504787062
    pw(15,1)=0.469355837986757	
    pw(15,2)=0.123252376810512
    pw(16,1)=0.587640403506911	
    pw(16,2)=0.112932296080539
    pw(17,1)=0.694487263186682
    pw(17,2)=0.100414144442881
    pw(18,1)=0.787816805979208	
    pw(18,2)=0.085941606217068
    pw(19,1)=0.865812577720300	
    pw(19,2)=0.069796468424521
    pw(20,1)=0.926956772187174	
    pw(20,2)=0.052293335152683
    pw(21,1)=0.970060497835428	
    pw(21,2)=0.033774901584814
    pw(22,1)=0.994294585482399	
    pw(22,2)=0.014627995298272
case(23)
    pw(1,1)=-0.994769334997552	
    pw(1,2)=0.013411859487142
    pw(2,1)=-0.972542471218115	
    pw(2,2)=0.030988005856979
    pw(3,1)=-0.932971086826016	
    pw(3,2)=0.048037671731085
    pw(4,1)=-0.876752358270441
    pw(4,2)=0.064232421408526
    pw(5,1)=-0.804888401618839	
    pw(5,2)=0.079281411776719
    pw(6,1)=-0.718661363131950	
    pw(6,2)=0.092915766060035
    pw(7,1)=-0.619609875763646	
    pw(7,2)=0.104892091464541
    pw(8,1)=-0.509501477846007	
    pw(8,2)=0.114996640222411
    pw(9,1)=-0.390301038030290	
    pw(9,2)=0.123049084306729
    pw(10,1)=-0.264135680970345	
    pw(10,2)=0.128905722188082
    pw(11,1)=-0.133256824298466
    pw(11,2)=0.132462039404696
    pw(12,1)=0.000000000000000
    pw(12,2)=0.133654572186106
    pw(13,1)=0.133256824298466	
    pw(13,2)=0.132462039404696
    pw(14,1)=0.264135680970345	
    pw(14,2)=0.128905722188082
    pw(15,1)=0.390301038030290	
    pw(15,2)=0.123049084306729
    pw(16,1)=0.509501477846007	
    pw(16,2)=0.114996640222411
    pw(17,1)=0.619609875763646	
    pw(17,2)=0.104892091464541
    pw(18,1)=0.718661363131950	
    pw(18,2)=0.092915766060035
    pw(19,1)=0.804888401618839	
    pw(19,2)=0.079281411776719
    pw(20,1)=0.876752358270441	
    pw(20,2)=0.064232421408526
    pw(21,1)=0.932971086826016	
    pw(21,2)=0.048037671731085
    pw(22,1)=0.972542471218115	
    pw(22,2)=0.030988005856979
    pw(23,1)=0.994769334997552		
    pw(23,2)=0.013411859487142
case(24)
    pw(1,1)=-0.995187219997021	
    pw(1,2)= 0.012341229799987
    pw(2,1)=-0.974728555971309	
    pw(2,2)= 0.028531388628934
    pw(3,1)=-0.938274552002732	
    pw(3,2)= 0.044277438817420
    pw(4,1)=-0.886415527004401	
    pw(4,2)= 0.059298584915437
    pw(5,1)=-0.820001985973902	
    pw(5,2)= 0.073346481411080
    pw(6,1)=-0.740124191578554	
    pw(6,2)= 0.086190161531953
    pw(7,1)=-0.648093651936975	
    pw(7,2)= 0.097618652104114
    pw(8,1)=-0.545421471388839	
    pw(8,2)= 0.107444270115965
    pw(9,1)=-0.433793507626045	
    pw(9,2)= 0.115505668053725
    pw(10,1)=-0.315042679696163	
    pw(10,2)= 0.121670472927803
    pw(11,1)=-0.191118867473616	
    pw(11,2)= 0.125837456346828
    pw(12,1)=-0.064056892862606	
    pw(12,2)= 0.127938195346752
    pw(13,1)= 0.064056892862606	
    pw(13,2)= 0.127938195346752
    pw(14,1)= 0.191118867473616	
    pw(14,2)= 0.125837456346828
    pw(15,1)= 0.315042679696163	
    pw(15,2)= 0.121670472927803
    pw(16,1)= 0.433793507626045	
    pw(16,2)= 0.115505668053725
    pw(17,1)= 0.545421471388839	
    pw(17,2)= 0.107444270115965
    pw(18,1)= 0.648093651936975	
    pw(18,2)= 0.097618652104114
    pw(19,1)= 0.740124191578554	
    pw(19,2)= 0.086190161531953
    pw(20,1)= 0.820001985973902	
    pw(20,2)= 0.073346481411080
    pw(21,1)= 0.886415527004401	
    pw(21,2)= 0.059298584915437
    pw(22,1)= 0.938274552002732	
    pw(22,2)= 0.044277438817420
    pw(23,1)= 0.974728555971309	
    pw(23,2)= 0.028531388628934
    pw(24,1)= 0.995187219997021	
    pw(24,2)= 0.012341229799987
case(25)
    pw(1,1)=-0.995556969790498
    pw(1,2)= 0.011393798501026
    pw(2,1)=-0.976663921459517	
    pw(2,2)= 0.026354986615032
    pw(3,1)=-0.942974571228974	
    pw(3,2)= 0.040939156701306
    pw(4,1)=-0.894991997878275	
    pw(4,2)= 0.054904695975835
    pw(5,1)=-0.833442628760834	
    pw(5,2)= 0.068038333812357
    pw(6,1)=-0.759259263037357	
    pw(6,2)= 0.080140700335001
    pw(7,1)=-0.673566368473468	
    pw(7,2)= 0.091028261982964
    pw(8,1)=-0.577662930241222	
    pw(8,2)= 0.100535949067050
    pw(9,1)=-0.473002731445715	
    pw(9,2)= 0.108519624474263
    pw(10,1)=-0.361172305809387	
    pw(10,2)= 0.114858259145711
    pw(11,1)=-0.243866883720988	
    pw(11,2)= 0.119455763535784
    pw(12,1)=-0.122864692610710	
    pw(12,2)= 0.122242442990310
    pw(13,1)=0.000000000000000	
    pw(13,2)= 0.123176053726715
    pw(14,1)= 0.122864692610710	
    pw(14,2)= 0.122242442990310
    pw(15,1)=0.243866883720988	
    pw(15,2)=0.119455763535784
    pw(16,1)=0.361172305809387	
    pw(16,2)=0.114858259145711
    pw(17,1)=0.473002731445715	
    pw(17,2)=0.108519624474263
    pw(18,1)=0.577662930241222	
    pw(18,2)=0.100535949067050
    pw(19,1)=0.673566368473468	
    pw(19,2)=0.091028261982964
    pw(20,1)=0.759259263037357	
    pw(20,2)=0.080140700335001
    pw(21,1)=0.833442628760834	
    pw(21,2)=0.068038333812357
    pw(22,1)=0.894991997878275	
    pw(22,2)=0.054904695975835
    pw(23,1)=0.942974571228974	
    pw(23,2)=0.040939156701306
    pw(24,1)=0.976663921459517	
    pw(24,2)=0.026354986615032
    pw(25,1)=0.995556969790498	
    pw(25,2)=0.011393798501026
case(26)
    pw(1,1)=-0.995885701145616	
    pw(1,2)=0.010551372617343
    pw(2,1)=-0.978385445956471	
    pw(2,2)=0.024417851092632
    pw(3,1)=-0.947159066661714	
    pw(3,2)=0.037962383294363
    pw(4,1)=-0.902637861984307	
    pw(4,2)=0.050975825297148
    pw(5,1)=-0.845445942788498	
    pw(5,2)=0.063274046329575
    pw(6,1)=-0.776385948820678	
    pw(6,2)=0.074684149765660
    pw(7,1)=-0.696427260419957	
    pw(7,2)=0.085045894313485
    pw(8,1)=-0.606692293017618	
    pw(8,2)=0.094213800355914
    pw(9,1)=-0.508440714824505	
    pw(9,2)=0.102059161094425
    pw(10,1)=-0.403051755123486	
    pw(10,2)=0.108471840528576
    pw(11,1)=-0.292004839485956	
    pw(11,2)=0.113361816546319
    pw(12,1)=-0.176858820356890	
    pw(12,2)=0.116660443485296
    pw(13,1)=-0.059230093429313	
    pw(13,2)=0.118321415279262
    pw(14,1)=0.059230093429313	
    pw(14,2)=0.118321415279262
    pw(15,1)=0.176858820356890	
    pw(15,2)=0.116660443485296
    pw(16,1)=0.292004839485956	
    pw(16,2)=0.113361816546319
    pw(17,1)=0.403051755123486	
    pw(17,2)=0.108471840528576
    pw(18,1)=0.508440714824505	
    pw(18,2)=0.102059161094425
    pw(19,1)=0.606692293017618	
    pw(19,2)=0.094213800355914
    pw(20,1)=0.696427260419957	
    pw(20,2)=0.085045894313485
    pw(21,1)=0.776385948820678	
    pw(21,2)=0.074684149765660
    pw(22,1)=0.845445942788498	
    pw(22,2)=0.063274046329575
    pw(23,1)=0.902637861984307	
    pw(23,2)=0.050975825297148
    pw(24,1)=0.947159066661714	
    pw(24,2)=0.037962383294363
    pw(25,1)=0.978385445956471	
    pw(25,2)=0.024417851092632
    pw(26,1)=0.995885701145616	
    pw(26,2)=0.010551372617343
case(27)
    pw(1,1)=-0.996179262888988	
    pw(1,2)=0.009798996051294
    pw(2,1)=-0.979923475961501	
    pw(2,2)=0.022686231596181
    pw(3,1)=-0.950900557814705	
    pw(3,2)=0.035297053757420
    pw(4,1)=-0.909482320677491	
    pw(4,2)=0.047449412520615
    pw(5,1)=-0.856207908018294	
    pw(5,2)=0.058983536859834
    pw(6,1)=-0.791771639070508	
    pw(6,2)=0.069748823766246
    pw(7,1)=-0.717013473739423	
    pw(7,2)=0.079604867773058
    pw(8,1)=-0.632907971946495	
    pw(8,2)=0.088423158543757
    pw(9,1)=-0.540551564579456	
    pw(9,2)=0.096088727370029
    pw(10,1)=-0.441148251750026	
    pw(10,2)=0.102501637817745
    pw(11,1)=-0.335993903638508	
    pw(11,2)=0.107578285788533
    pw(12,1)=-0.226459365439536	
    pw(12,2)=0.111252488356845
    pw(13,1)=-0.113972585609530	
    pw(13,2)=0.113476346108965
    pw(14,1)=0.000000000000000	
    pw(14,2)=0.114220867378957
    pw(15,1)=0.113972585609530	
    pw(15,2)=0.113476346108965
    pw(16,1)=0.226459365439536	
    pw(16,2)=0.111252488356845
    pw(17,1)=0.335993903638508	
    pw(17,2)=0.107578285788533
    pw(18,1)=0.441148251750026	
    pw(18,2)=0.102501637817745
    pw(19,1)=0.540551564579456	
    pw(19,2)=0.096088727370029
    pw(20,1)=0.632907971946495	
    pw(20,2)=0.088423158543757
    pw(21,1)=0.717013473739423	
    pw(21,2)=0.079604867773058
    pw(22,1)=0.791771639070508	
    pw(22,2)=0.069748823766246
    pw(23,1)=0.856207908018294	
    pw(23,2)=0.058983536859834
    pw(24,1)=0.909482320677491	
    pw(24,2)=0.047449412520615
    pw(25,1)=0.950900557814705	
    pw(25,2)=0.035297053757420
    pw(26,1)=0.979923475961501	
    pw(26,2)=0.022686231596181
	pw(27,1)=0.996179262888988  
    pw(27,2)=0.009798996051294
case(28)
    pw(1,1)=-0.996442497573954	
    pw(1,2)=0.009124282593095
    pw(2,1)=-0.981303165370872	
    pw(2,2)=0.021132112592771
    pw(3,1)=-0.954259280628938	
    pw(3,2)=0.032901427782304
    pw(4,1)=-0.915633026392132	
    pw(4,2)=0.044272934759004
    pw(5,1)=-0.865892522574395	
    pw(5,2)=0.055107345675717
    pw(6,1)=-0.805641370917179	
    pw(6,2)=0.065272923967000
    pw(7,1)=-0.735610878013631	
    pw(7,2)=0.074646214234569
    pw(8,1)=-0.656651094038865	
    pw(8,2)=0.083113417228901
    pw(9,1)=-0.569720471811401	
    pw(9,2)=0.090571744393033
    pw(10,1)=-0.475874224955118	
    pw(10,2)=0.096930657997930
    pw(11,1)=-0.376251516089078	
    pw(11,2)=0.102112967578060
    pw(12,1)=-0.272061627635178	
    pw(12,2)=0.106055765922846
    pw(13,1)=-0.164569282133380	
    pw(13,2)=0.108711192258294
    pw(14,1)=-0.055079289884034	
    pw(14,2)=0.110047013016475
    pw(15,1)=0.055079289884034	
    pw(15,2)=0.110047013016475
    pw(16,1)=0.164569282133380	
    pw(16,2)=0.108711192258294
    pw(17,1)=0.272061627635178	
    pw(17,2)=0.106055765922846
    pw(18,1)=0.376251516089078	
    pw(18,2)=0.102112967578060
    pw(19,1)=0.475874224955118	
    pw(19,2)=0.096930657997930
    pw(20,1)=0.569720471811401	
    pw(20,2)=0.090571744393033
    pw(21,1)=0.656651094038865	
    pw(21,2)=0.083113417228901
    pw(22,1)=0.735610878013631	
    pw(22,2)=0.074646214234569
    pw(23,1)=0.805641370917179	
    pw(23,2)=0.065272923967000
    pw(24,1)=0.865892522574395	
    pw(24,2)=0.055107345675717
    pw(25,1)=0.915633026392132	
    pw(25,2)=0.044272934759004
    pw(26,1)=0.954259280628938	
    pw(26,2)=0.032901427782304
    pw(27,1)=0.981303165370872	
    pw(27,2)=0.021132112592771
    pw(28,1)=0.996442497573954	
    pw(28,2)=0.009124282593095
case(29)
    pw(1,1)=-0.996679442260596
    pw(1,2)=0.008516903878746
    pw(2,1)=-0.982545505261413	
    pw(2,2)=0.019732085056123
    pw(3,1)=-0.957285595778087	
    pw(3,2)=0.030740492202094
    pw(4,1)=-0.921180232953058	
    pw(4,2)=0.041402062518683
    pw(5,1)=-0.874637804920102	
    pw(5,2)=0.051594826902498
    pw(6,1)=-0.818185487615252	
    pw(6,2)=0.061203090657079
    pw(7,1)=-0.752462851734477	
    pw(7,2)=0.070117933255051
    pw(8,1)=-0.678214537602686	
    pw(8,2)=0.078238327135764
    pw(9,1)=-0.596281797138227	
    pw(9,2)=0.085472257366173
    pw(10,1)=-0.507592955124227	
    pw(10,2)=0.091737757139259
    pw(11,1)=-0.413152888174008	
    pw(11,2)=0.096963834094409
    pw(12,1)=-0.314031637867639	
    pw(12,2)=0.101091273759915
    pw(13,1)=-0.211352286166001	
    pw(13,2)=0.104073310077729
    pw(14,1)=-0.106278230132679	
    pw(14,2)=0.105876155097320
    pw(15,1)=0.000000000000000	
    pw(15,2)=0.106479381718314
    pw(16,1)=0.106278230132679	
    pw(16,2)=0.105876155097320
    pw(17,1)=0.211352286166001	
    pw(17,2)=0.104073310077729
    pw(18,1)=0.314031637867639	
    pw(18,2)=0.101091273759915
    pw(19,1)=0.413152888174008	
    pw(19,2)=0.096963834094409
    pw(20,1)=0.507592955124227	
    pw(20,2)=0.091737757139259
    pw(21,1)=0.596281797138227	
    pw(21,2)=0.085472257366173
    pw(22,1)=0.678214537602686	
    pw(22,2)=0.078238327135764
    pw(23,1)=0.752462851734477	
    pw(23,2)=0.070117933255051
    pw(24,1)=0.818185487615252	
    pw(24,2)=0.061203090657079
    pw(25,1)=0.874637804920102	
    pw(25,2)=0.051594826902498
    pw(26,1)=0.921180232953058	
    pw(26,2)=0.041402062518683
    pw(27,1)=0.957285595778087	
    pw(27,2)=0.030740492202094
    pw(28,1)=0.982545505261413	
    pw(28,2)=0.019732085056123
    pw(29,1)=0.996679442260596	
    pw(29,2)=0.008516903878746
case(30)
    pw(1,1)=-0.996893484074649	
    pw(1,2)=0.007968192496167
    pw(2,1)=-0.983668123279747	
    pw(2,2)=0.018466468311091
    pw(3,1)=-0.960021864968307	
    pw(3,2)=0.028784707883323
    pw(4,1)=-0.926200047429274	
    pw(4,2)=0.038799192569627
    pw(5,1)=-0.882560535792052	
    pw(5,2)=0.048402672830594
    pw(6,1)=-0.829565762382768	
    pw(6,2)=0.057493156217619
    pw(7,1)=-0.767777432104826	
    pw(7,2)=0.065974229882181
    pw(8,1)=-0.697850494793315	
    pw(8,2)=0.073755974737705
    pw(9,1)=-0.620526182989242	
    pw(9,2)=0.080755895229420
    pw(10,1)=-0.536624148142019	
    pw(10,2)=0.086899787201083
    pw(11,1)=-0.447033769538089	
    pw(11,2)=0.092122522237786
    pw(12,1)=-0.352704725530878	
    pw(12,2)=0.096368737174644
    pw(13,1)=-0.254636926167889	
    pw(13,2)=0.099593420586795
    pw(14,1)=-0.153869913608583	
    pw(14,2)=0.101762389748405
    pw(15,1)=-0.051471842555318	
    pw(15,2)=0.102852652893558
    pw(16,1)=0.051471842555318	
    pw(16,2)=0.102852652893558
    pw(17,1)=0.153869913608583	
    pw(17,2)=0.101762389748405
    pw(18,1)=0.254636926167889	
    pw(18,2)=0.099593420586795
    pw(19,1)=0.352704725530878	
    pw(19,2)=0.096368737174644
    pw(20,1)=0.447033769538089	
    pw(20,2)=0.092122522237786
    pw(21,1)=0.536624148142019	
    pw(21,2)=0.086899787201083
    pw(22,1)=0.620526182989242	
    pw(22,2)=0.080755895229420
    pw(23,1)=0.697850494793315	
    pw(23,2)=0.073755974737705
    pw(24,1)=0.767777432104826	
    pw(24,2)=0.065974229882181
    pw(25,1)=0.829565762382768	
    pw(25,2)=0.057493156217619
    pw(26,1)=0.882560535792052	
    pw(26,2)=0.048402672830594
    pw(27,1)=0.926200047429274	
    pw(27,2)=0.038799192569627
    pw(28,1)=0.960021864968307	
    pw(28,2)=0.028784707883323
    pw(29,1)=0.983668123279747	
    pw(29,2)=0.018466468311091
	pw(30,1)=0.996893484074649		
    pw(30,2)=0.007968192496167
case(31:)
    write (*,*)'Error, The number of Gauss points must be >=2 and <=30.'
end select 
  
end subroutine Gauss_ponits_weights


end module subs_fs3d_qua_num













