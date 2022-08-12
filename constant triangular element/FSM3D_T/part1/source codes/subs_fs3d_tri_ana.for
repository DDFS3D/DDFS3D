      module subs_fs3d_tri_ana
      contains 
      
      subroutine diag_stress_fs3d_tri_ana (xx1,yy1,zz1,xx2,yy2,zz2,xx3,
     &yy3,zz3,E,v,IEvect,COEFF)
      
      real*8 x(5),y(5),z(5),a(4),b(4),pp(4,5),lxi,lyi,mi,li
      real*8 f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,
     &rii,rxi,rxii,ryii,d1,d2,d3,d4,d5,d6,d7,d0,ri,ryi
      real*8 s1,t1,s21,t21,s22,t22,g1,h1,g3,h3,g4,h4
      real*8 x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
      
      real*8::E,v,IEvect(3,3),JEvect(3,3),COEFF(3,3),JI(3,3),trac(3)
      real*8::S(6,3),xx(4),yy(4),zz(4),gpx,gpy,gpz,GG,LOCALZ
      real*8::xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,cons
      integer::i,j
*-----(x1,y2,z3)
      
       GG=E/(1+v)/2.
       cons=1./8./3.141592653/(1-v)
       LOCALZ=0.
      
       
      DO I=1,3
        DO J=1,3
            JEvect(i,j)=IEvect(i,j)
        end do 
      end do 
      
       gpx=(xx1+xx2+xx3)/3.
       gpy=(yy1+yy2+yy3)/3.
       gpz=(zz1+zz2+zz3)/3.

       xx(1)=xx1-xx1
       yy(1)=yy1-yy1
       zz(1)=zz1-zz1

       xx(2)=xx2-xx1
       yy(2)=yy2-yy1
       zz(2)=zz2-zz1

       xx(3)=xx3-xx1
       yy(3)=yy3-yy1
       zz(3)=zz3-zz1
       
      
       x(2)=xx(1)*IEvect(1,1)+yy(1)*IEvect(1,2)+zz(1)*IEvect(1,3)
       y(2)=xx(1)*IEvect(2,1)+yy(1)*IEvect(2,2)+zz(1)*IEvect(2,3)
       z(2)=xx(1)*IEvect(3,1)+yy(1)*IEvect(3,2)+zz(1)*IEvect(3,3)

       x(3)=xx(2)*IEvect(1,1)+yy(2)*IEvect(1,2)+zz(2)*IEvect(1,3)
       y(3)=xx(2)*IEvect(2,1)+yy(2)*IEvect(2,2)+zz(2)*IEvect(2,3)
       z(3)=xx(2)*IEvect(3,1)+yy(2)*IEvect(3,2)+zz(2)*IEvect(3,3)

       x(4)=xx(3)*IEvect(1,1)+yy(3)*IEvect(1,2)+zz(3)*IEvect(1,3)
       y(4)=xx(3)*IEvect(2,1)+yy(3)*IEvect(2,2)+zz(3)*IEvect(2,3)
       z(4)=xx(3)*IEvect(3,1)+yy(3)*IEvect(3,2)+zz(3)*IEvect(3,3)   
      
       
            
       x(1)=(x(2)+x(3)+x(4))/3.
       y(1)=(y(2)+y(3)+y(4))/3.
       z(1)=(z(2)+z(3)+z(4))/3.  
       
      x(5)=x(2)
      y(5)=y(2)
      z(5)=z(2)
      
      !do i=1,5
      !    write(*,*)i,'x,y,z=',x(i),y(i),z(i)
      !end do 
      !
      
      
      do 10 i=1,4
         do 15 j=1,5
            pp(i,j)=dsqrt((x(j)-x(i))**2.D0+(y(j)-y(i))**2.D0
     &	          +(z(j)-z(i))**2.D0)
   15    continue
   10 continue
      f1=0.0D0
      f2=0.0D0
      f3=0.0D0
      f4=0.0D0
      f5=0.0D0
      f6=0.0D0
      f7=0.0D0
      f8=0.0D0
      f9=0.0D0
      f10=0.0D0
      f11=0.0D0
      f12=0.0D0
      f13=0.0D0
      f14=0.0D0
      f15=0.0D0
      f16=0.0D0
      f17=0.0D0
      do 20 i=2,4
         a(i)=((x(i)-x(1))*(x(i+1)-x(i))+(y(i)-y(1))*(y(i+1)-y(i)))
     &   /pp(i,i+1)
 
         b(i)=((x(i+1)-x(1))*(x(i+1)-x(i))+(y(i+1)-y(1))*(y(i+1)-y(
     &   i)))/pp(i,i+1)
   20 continue
 
      do 30 i=2,4
         mi=(x(i)-x(1))*(y(i+1)-y(1))-(x(i+1)-x(1))*(y(i)-y(1))
         li=pp(i,i+1)
         ri=pp(1,i)
         rii=pp(1,i+1)
         rxi=x(i)-x(1)
         ryi=y(i)-y(1)
         lxi=x(i+1)-x(i)
         lyi=y(i+1)-y(i)
         rxii=x(i+1)-x(1)
         ryii=y(i+1)-y(1)
*
         d0=dsign(1.0D0,z(1))
         d1=dsign(1.0D0,lxi)
         d2=dsign(1.0D0,lyi)
         d3=dsign(1.0D0,mi)
         d4=dsign(1.0D0,rxi)
         d5=dsign(1.0D0,ryi)
         d6=dsign(1.0D0,rxii)
         d7=dsign(1.0D0,ryii)
*
         if (dabs(z(1)).le.1.0D-3) then 
            if (dabs(lxi).ge.1.0D-3) then 
               s1 =0.0D0
               s21=(-1)*d1*d3*d6*3.14D0/2.0D0
               t21=(-1)*d1*d3*d4*3.14D0/2.0D0
            else
               s1 =0.0D0
               s21=(-1)*d3*d6*3.14D0/2.0D0
               t21=(-1)*d3*d4*3.14D0/2.0D0
            end if
         end if
         if (dabs(z(1)).ge.1.0D-3) then 
            if (dabs(lxi).ge.1.0D-3) then 
               s1=-d1*dabs(z(1))*(datan(li/dabs(lxi)/dabs(z(1))*(ryii
     &          +b(i)+(lyi/li+1)*rii))-datan(li/dabs(lxi)/dabs(z(
     &         1))*(ryi+a(i)+(lyi/li+1)*ri)))
               s21=datan(-(rxii*mi+lyi*z(1)**2.D0)/(dabs(z(1)))/lxi/rii)
               t21=datan(-(rxi *mi+lyi*z(1)**2.D0)/(dabs(z(1)))/lxi/ri )
            else
               s1 =0.0D0
               ds21=dsign(1.0D0,rxii*mi+lyi*z(1)**2.D0)
               if ((rxii*mi+lyi*z(1)**2.D0).le.1.0D-3) then ds21=0.0D0
               dt21=dsign(1.0D0,rxi *mi+lyi*z(1)**2.D0)
               if ((rxi *mi+lyi*z(1)**2).le.1.0D-3) then dt21=0.0D0
               s21=(-1.D0)*3.14D0/2.0D0*ds21
               t21=(-1.D0)*3.14D0/2.0D0*dt21
            end if
         end if
 
         if (dabs(z(1)).le.1.0D-3) then 
            if (dabs(lyi).ge.1.0D-3) then 
               t1 =0.0D0
               s22=(-1.D0)*d2*d3*d7*3.14D0/2.0D0
               t22=(-1.D0)*d2*d3*d5*3.14D0/2.0D0
            else
               t1 =0.0D0
               s22=(-1.D0)*d3*d7*3.14D0/2.0D0
               t22=(-1.D0)*d3*d5*3.14D0/2.0D0
            end if
         end if
         if (dabs(z(1)).ge.1.0D-3) then 
            if (dabs(lyi).ge.1.0D-3) then 
               t1= d2*dabs(z(1))*(datan(li/dabs(lyi)/dabs(z(1))*(rxii+
     &         b(i)+(lxi/li+1)*rii))-datan(li/dabs(lyi)/dabs(z(1)
     &         )*(rxi+a(i)+(lxi/li+1)*ri)))
               s22=datan(-(ryii*mi-lxi*z(1)**2.D0)/(dabs(z(1)))/lyi/rii)
               t22=datan(-(ryi *mi-lxi*z(1)**2.D0)/(dabs(z(1)))/lyi/ri )
            else
               t1 =0.0D0
               ds22=dsign(1.0D0,ryii*mi-lxi*z(1)**2)
               if ((ryii*mi-lxi*z(1)**2.D0).le.1.0D-3) then ds22=0.0D0
               dt22=dsign(1.0D0,ryi *mi-lxi*z(1)**2)
               if ((ryi *mi-lxi*z(1)**2.D0).le.1.0D-3) then dt22=0.0D0
               s22=(-1.D0)*3.14/2.0D0*ds22
               t22=(-1.D0)*3.14/2.0D0*dt22
            end if
         end if
         g1=(b(i)/rii-a(i)/ri)/(mi**2.D0+z(1)**2.D0*li**2.D0)
         h1=dlog(dabs((b(i)+rii)/(a(i)+ri)))
 
         g3=-d7*dlog((rii-dabs(ryii))/(rii+dabs(ryii)))
         h3=-d5*dlog((ri -dabs(ryi ))/(ri +dabs(ryi )))
         g4= d6*dlog((rii-dabs(rxii))/(rii+dabs(rxii)))
         h4= d4*dlog((ri -dabs(rxi ))/(ri +dabs(rxi )))
         
       
 
         f1=f1+s1+t1+mi/li*log(dabs((b(i)+rii)/(a(i)+ri)))
 
         f2=f2-d0*((s21-t21)+(s22-t22))/2.D0
 
         f3=f3+lyi/li*dlog(dabs((b(i)+rii)/(a(i)+ri)))+(g3-h3)/4.D0
 
         f4=f4-lxi/li*dlog(dabs((b(i)+rii)/(a(i)+ri)))+(g4-h4)/4.D0
 
         
         
         f5=f5-lxi*lyi/li**2.D0*(rii-ri)-
     &                    lyi**2.D0*mi/li**3.D0*h1+mi/li*h1+2.D0*s1
 
         f6=f6+lxi*lyi/li**2.D0*(rii-ri)-
     &                    lxi**2.D0*mi/li**3.D0*h1+mi/li*h1+2.D0*t1
 
         f7=f7-(lyi**2.D0-lxi**2.D0)/li**2.D0*(rii-ri)/2.D0+
     &                                        lxi*lyi*mi/li**3.D0*h1
 
         f8=f8+d0*(-s21+t21-s22+t22)/6.D0+li*mi*g1/3.D0*z(1)
 
         f9 =f9 -mi*lxi/li**2.D0*(1.D0/rii-1.D0/ri)/3.D0+
     &                  lyi*g1*(2.D0*li*z(1)**2.D0+mi**2.D0/li)/3.D0
 
         f10=f10-lyi*mi/li**2.D0*(1.D0/rii-1.D0/ri)/3.D0-
     &                  lxi*g1*(2.D0*li*z(1)**2.D0+mi**2.D0/li)/3.D0
 
         f11=f11+d0*(-s21+t21+(s22-t22)*2.D0/3.D0)+z(1)*(lyi*lxi/
     &             li**2.D0*(1.D0/rii-1.D0/ri)-lyi**2.D0*mi/li*g1)/3.D0

         f12=f12+d0*((s21-t21)*2/3-s22+t22)-z(1)*(lyi*lxi/
     &             li**2.D0*(1.D0/rii-1.D0/ri)+lxi**2.D0*mi/li*g1)/3.D0
 
         f13=f13-(lxi**2.D0-lyi**2.D0)/li**2.D0*(1.D0/rii-1.D0/ri)/6.D0
     &             +2.D0*lxi*lyi*mi/li*g1/6.D0
 
         f16=f16+(lyi/li*h1+(g3-h3)/2.D0
     &           +lyi*(lyi**2.D0-lxi**2.D0)/li**3.D0*h1)/6.D0
     &           +lxi*mi/li**4.D0*(3.D0*lyi**2.D0-lxi**2.D0)
     &           *(1.D0/rii-1.D0/ri)/6.D0+lyi*mi**2.D0*(3.D0*lxi**2.D0
     &           -lyi**2.D0)/li**3.D0*g1/6+z(1)**2.D0*lyi*(lxi**2.D0-
     &           lyi**2.D0)/li*g1/6.D0
 
         f17=f17-(lxi/li*h1-(g4-h4)/2.D0+lxi*(lxi**2-lyi**2)/li**3*h1)
     &       /6.D0+lyi*mi/li**4*(3.D0*lxi**2-lyi**2)*(1.D0/rii-1.D0/ri)
     &       /6.D0-lxi*mi**2*(3.D0*lyi**2-lxi**2)/li**3*g1/6.D0
     &       -z(1)**2*lxi*(lyi**2-lxi**2)/li*g1/6.D0 
 
            !if (isnan(f1))then
            !    write(*,*)i,j,'ai=',a(i),'ri=',ri
            !end if 
                   
   30 continue
 
      f14=f3-f16-f9
      f15=f4-f17-f10
 
      f2 = -dabs(f2)
      f8 = -dabs(f8)
      f11= -dabs(f11)
      f12= -dabs(f12)
      
      !write(*,*)'f1=',f1
      !write(*,*)'f2=',f2
      !write(*,*)'f3=',f3
      !write(*,*)'f4=',f4
      !write(*,*)'f5=',f5
      !write(*,*)'f6=',f6
      !write(*,*)'f7=',f7
      !write(*,*)'f8=',f8
      !write(*,*)'f9=',f9
      !write(*,*)'f10=',f10
      !write(*,*)'f11=',f11
      !write(*,*)'f12=',f12
      !write(*,*)'f13=',f13
      !write(*,*)'f14=',f14
      !write(*,*)'f15=',f15
      !write(*,*)'f16=',f16
      !write(*,*)'f17=',f17
      
         
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
      !
      !write(*,*) 'S(1,1)', S(1,1)
      !write(*,*) 'S(1,2)', S(1,2)
      !write(*,*) 'S(1,3)', S(1,3)
      !write(*,*) 'S(2,1)', S(2,1)
      !write(*,*) 'S(2,2)', S(2,2)
      !write(*,*) 'S(2,3)', S(2,3)
      !write(*,*) 'S(3,1)', S(3,1)
      !write(*,*) 'S(3,2)', S(3,2)
      !write(*,*) 'S(3,3)', S(3,3)
      !write(*,*) 'S(4,1)', S(4,1)
      !write(*,*) 'S(4,2)', S(4,2)
      !write(*,*) 'S(4,3)', S(4,3)
      !write(*,*) 'S(5,1)', S(5,1)
      !write(*,*) 'S(5,2)', S(5,2)
      !write(*,*) 'S(5,3)', S(5,3)
      !write(*,*) 'S(6,1)', S(6,1)
      !write(*,*) 'S(6,2)', S(6,2)
      !write(*,*) 'S(6,3)', S(6,3)
      !!! 
      
      JI(1,1)=JEvect(1,1)*IEvect(1,1)+JEvect(1,2)*IEvect(1,2)+JEvect(1,3
     &)*IEvect(1,3)
      JI(2,1)=JEvect(2,1)*IEvect(1,1)+JEvect(2,2)*IEvect(1,2)+JEvect(2,3
     &)*IEvect(1,3)
      JI(3,1)=JEvect(3,1)*IEvect(1,1)+JEvect(3,2)*IEvect(1,2)+JEvect(3,3
     &)*IEvect(1,3)
      
      JI(1,2)=JEvect(1,1)*IEvect(2,1)+JEvect(1,2)*IEvect(2,2)+JEvect(1,3
     &)*IEvect(2,3)
      JI(2,2)=JEvect(2,1)*IEvect(2,1)+JEvect(2,2)*IEvect(2,2)+JEvect(2,3
     &)*IEvect(2,3)
      JI(3,2)=JEvect(3,1)*IEvect(2,1)+JEvect(3,2)*IEvect(2,2)+JEvect(3,3
     &)*IEvect(2,3) 
      
      JI(1,3)=JEvect(1,1)*IEvect(3,1)+JEvect(1,2)*IEvect(3,2)+JEvect(1,3
     &)*IEvect(3,3)
      JI(2,3)=JEvect(2,1)*IEvect(3,1)+JEvect(2,2)*IEvect(3,2)+JEvect(2,3
     &)*IEvect(3,3)
      JI(3,3)=JEvect(3,1)*IEvect(3,1)+JEvect(3,2)*IEvect(3,2)+JEvect(3,3
     &)*IEvect(3,3)  

      !write(*,*)'JI(1,1)',JI(1,1)
      !write(*,*)'JI(2,1)',JI(2,1)
      !write(*,*)'JI(3,1)',JI(3,1)
      !write(*,*)'JI(1,2)',JI(1,2)
      !write(*,*)'JI(2,2)',JI(2,2)
      !write(*,*)'JI(3,2)',JI(3,2) 
      !write(*,*)'JI(1,3)',JI(1,3)
      !write(*,*)'JI(2,3)',JI(2,3)
      !write(*,*)'JI(3,3)',JI(3,3)
     
     
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
     
      !write(*,*)'coeff(1,1)',coeff(1,1)
      !write(*,*)'coeff(2,1)',coeff(2,1)
      !write(*,*)'coeff(3,1)',coeff(3,1)
      !write(*,*)'coeff(1,2)',coeff(1,2)
      !write(*,*)'coeff(2,2)',coeff(2,2)
      !write(*,*)'coeff(3,2)',coeff(3,2)
      !write(*,*)'coeff(1,3)',coeff(1,3)
      !write(*,*)'coeff(2,3)',coeff(2,3)
      !write(*,*)'coeff(3,3)',coeff(3,3) 

      end subroutine diag_stress_fs3d_tri_ana
      
      
      subroutine diag_dis_fs3d_tri_ana (xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,
     &zz3,E,v,IEvect,COEFF)

       real*8 x(5),y(5),z(5),a(4),b(4),pp(4,5),lxi,lyi,mi,li
      real*8 f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,
     &rii,rxi,rxii,ryii,d1,d2,d3,d4,d5,d6,d7,d0,ri,ryi
      real*8 s1,t1,s21,t21,s22,t22,g1,h1,g3,h3,g4,h4
      real*8 x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
      
      real*8::E,v,IEvect(3,3),JEvect(3,3),COEFF(3,3),JI(3,3),trac(3)
      real*8::u(3,3),xx(4),yy(4),zz(4),gpx,gpy,gpz,GG,LOCALZ
      real*8::xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,cond
      integer::i,j
*-----(x1,y2,z3)
!      
       GG=E/2/(1+v)
       cond=1./16./3.141592653/(1-v)/GG
       LOCALZ=0.
!      
      DO I=1,3
        DO J=1,3
            JEvect(i,j)=IEvect(i,j)
        end do 
      end do 
!      
       gpx=(xx1+xx2+xx3)/3.
       gpy=(yy1+yy2+yy3)/3.
       gpz=(zz1+zz2+zz3)/3.

       xx(1)=xx1-xx1
       yy(1)=yy1-yy1
       zz(1)=zz1-zz1

       xx(2)=xx2-xx1
       yy(2)=yy2-yy1
       zz(2)=zz2-zz1

       xx(3)=xx3-xx1
       yy(3)=yy3-yy1
       zz(3)=zz3-zz1
       
      
       x(2)=xx(1)*IEvect(1,1)+yy(1)*IEvect(1,2)+zz(1)*IEvect(1,3)
       y(2)=xx(1)*IEvect(2,1)+yy(1)*IEvect(2,2)+zz(1)*IEvect(2,3)
       z(2)=xx(1)*IEvect(3,1)+yy(1)*IEvect(3,2)+zz(1)*IEvect(3,3)

       x(3)=xx(2)*IEvect(1,1)+yy(2)*IEvect(1,2)+zz(2)*IEvect(1,3)
       y(3)=xx(2)*IEvect(2,1)+yy(2)*IEvect(2,2)+zz(2)*IEvect(2,3)
       z(3)=xx(2)*IEvect(3,1)+yy(2)*IEvect(3,2)+zz(2)*IEvect(3,3)

       x(4)=xx(3)*IEvect(1,1)+yy(3)*IEvect(1,2)+zz(3)*IEvect(1,3)
       y(4)=xx(3)*IEvect(2,1)+yy(3)*IEvect(2,2)+zz(3)*IEvect(2,3)
       z(4)=xx(3)*IEvect(3,1)+yy(3)*IEvect(3,2)+zz(3)*IEvect(3,3)   
      
       
            
       x(1)=(x(2)+x(3)+x(4))/3.
       y(1)=(y(2)+y(3)+y(4))/3.
       z(1)=(z(2)+z(3)+z(4))/3.  
       
      x(5)=x(2)
      y(5)=y(2)
      z(5)=z(2)
      
      !do i=1,5
      !    write(*,*)i,'x,y,z=',x(i),y(i),z(i)
      !end do 
      !
      
!      
      do 10 i=1,4
         do 15 j=1,5
            pp(i,j)=dsqrt((x(j)-x(i))**2.D0+(y(j)-y(i))**2.D0
     &	          +(z(j)-z(i))**2.D0)
   15    continue
   10 continue
      f1=0.0D0
      f2=0.0D0
      f3=0.0D0
      f4=0.0D0
      f5=0.0D0
      f6=0.0D0
      f7=0.0D0
      f8=0.0D0
      f9=0.0D0
      f10=0.0D0
      f11=0.0D0
      f12=0.0D0
      f13=0.0D0
      f14=0.0D0
      f15=0.0D0
      f16=0.0D0
      f17=0.0D0
      do 20 i=2,4
         a(i)=((x(i)-x(1))*(x(i+1)-x(i))+(y(i)-y(1))*(y(i+1)-y(i)))
     &   /pp(i,i+1)
 
         b(i)=((x(i+1)-x(1))*(x(i+1)-x(i))+(y(i+1)-y(1))*(y(i+1)-y(
     &   i)))/pp(i,i+1)
   20 continue
! 
      do 30 i=2,4
         mi=(x(i)-x(1))*(y(i+1)-y(1))-(x(i+1)-x(1))*(y(i)-y(1))
         li=pp(i,i+1)
         ri=pp(1,i)
         rii=pp(1,i+1)
         rxi=x(i)-x(1)
         ryi=y(i)-y(1)
         lxi=x(i+1)-x(i)
         lyi=y(i+1)-y(i)
         rxii=x(i+1)-x(1)
         ryii=y(i+1)-y(1)
*
         d0=dsign(1.0D0,z(1))
         d1=dsign(1.0D0,lxi)
         d2=dsign(1.0D0,lyi)
         d3=dsign(1.0D0,mi)
         d4=dsign(1.0D0,rxi)
         d5=dsign(1.0D0,ryi)
         d6=dsign(1.0D0,rxii)
         d7=dsign(1.0D0,ryii)
*
         if (dabs(z(1)).le.1.0D-3) then 
            if (dabs(lxi).ge.1.0D-3) then 
               s1 =0.0D0
               s21=(-1)*d1*d3*d6*3.14D0/2.0D0
               t21=(-1)*d1*d3*d4*3.14D0/2.0D0
            else
               s1 =0.0D0
               s21=(-1)*d3*d6*3.14D0/2.0D0
               t21=(-1)*d3*d4*3.14D0/2.0D0
            end if
         end if
         if (dabs(z(1)).ge.1.0D-3) then 
            if (dabs(lxi).ge.1.0D-3) then 
               s1=-d1*dabs(z(1))*(datan(li/dabs(lxi)/dabs(z(1))*(ryii
     &          +b(i)+(lyi/li+1)*rii))-datan(li/dabs(lxi)/dabs(z(
     &         1))*(ryi+a(i)+(lyi/li+1)*ri)))
               s21=datan(-(rxii*mi+lyi*z(1)**2.D0)/(dabs(z(1)))/lxi/rii)
               t21=datan(-(rxi *mi+lyi*z(1)**2.D0)/(dabs(z(1)))/lxi/ri )
            else
               s1 =0.0D0
               ds21=dsign(1.0D0,rxii*mi+lyi*z(1)**2.D0)
               if ((rxii*mi+lyi*z(1)**2.D0).le.1.0D-3) then ds21=0.0D0
               dt21=dsign(1.0D0,rxi *mi+lyi*z(1)**2.D0)
               if ((rxi *mi+lyi*z(1)**2).le.1.0D-3) then dt21=0.0D0
               s21=(-1.D0)*3.14D0/2.0D0*ds21
               t21=(-1.D0)*3.14D0/2.0D0*dt21
            end if
         end if
 
         if (dabs(z(1)).le.1.0D-3) then 
            if (dabs(lyi).ge.1.0D-3) then 
               t1 =0.0D0
               s22=(-1.D0)*d2*d3*d7*3.14D0/2.0D0
               t22=(-1.D0)*d2*d3*d5*3.14D0/2.0D0
            else
               t1 =0.0D0
               s22=(-1.D0)*d3*d7*3.14D0/2.0D0
               t22=(-1.D0)*d3*d5*3.14D0/2.0D0
            end if
         end if
         if (dabs(z(1)).ge.1.0D-3) then 
            if (dabs(lyi).ge.1.0D-3) then 
               t1= d2*dabs(z(1))*(datan(li/dabs(lyi)/dabs(z(1))*(rxii+
     &         b(i)+(lxi/li+1)*rii))-datan(li/dabs(lyi)/dabs(z(1)
     &         )*(rxi+a(i)+(lxi/li+1)*ri)))
               s22=datan(-(ryii*mi-lxi*z(1)**2.D0)/(dabs(z(1)))/lyi/rii)
               t22=datan(-(ryi *mi-lxi*z(1)**2.D0)/(dabs(z(1)))/lyi/ri )
            else
               t1 =0.0D0
               ds22=dsign(1.0D0,ryii*mi-lxi*z(1)**2)
               if ((ryii*mi-lxi*z(1)**2.D0).le.1.0D-3) then ds22=0.0D0
               dt22=dsign(1.0D0,ryi *mi-lxi*z(1)**2)
               if ((ryi *mi-lxi*z(1)**2.D0).le.1.0D-3) then dt22=0.0D0
               s22=(-1.D0)*3.14/2.0D0*ds22
               t22=(-1.D0)*3.14/2.0D0*dt22
            end if
         end if
         g1=(b(i)/rii-a(i)/ri)/(mi**2.D0+z(1)**2.D0*li**2.D0)
         h1=dlog(dabs((b(i)+rii)/(a(i)+ri)))
 
         g3=-d7*dlog((rii-dabs(ryii))/(rii+dabs(ryii)))
         h3=-d5*dlog((ri -dabs(ryi ))/(ri +dabs(ryi )))
         g4= d6*dlog((rii-dabs(rxii))/(rii+dabs(rxii)))
         h4= d4*dlog((ri -dabs(rxi ))/(ri +dabs(rxi )))
         
       
 
         f1=f1+s1+t1+mi/li*log(dabs((b(i)+rii)/(a(i)+ri)))
 
         f2=f2-d0*((s21-t21)+(s22-t22))/2.D0
 
         f3=f3+lyi/li*dlog(dabs((b(i)+rii)/(a(i)+ri)))+(g3-h3)/4.D0
 
         f4=f4-lxi/li*dlog(dabs((b(i)+rii)/(a(i)+ri)))+(g4-h4)/4.D0
 
         
         
         f5=f5-lxi*lyi/li**2.D0*(rii-ri)-
     &                    lyi**2.D0*mi/li**3.D0*h1+mi/li*h1+2.D0*s1
 
         f6=f6+lxi*lyi/li**2.D0*(rii-ri)-
     &                    lxi**2.D0*mi/li**3.D0*h1+mi/li*h1+2.D0*t1
 
         f7=f7-(lyi**2.D0-lxi**2.D0)/li**2.D0*(rii-ri)/2.D0+
     &                                        lxi*lyi*mi/li**3.D0*h1
 
         f8=f8+d0*(-s21+t21-s22+t22)/6.D0+li*mi*g1/3.D0*z(1)
 
         f9 =f9 -mi*lxi/li**2.D0*(1.D0/rii-1.D0/ri)/3.D0+
     &                  lyi*g1*(2.D0*li*z(1)**2.D0+mi**2.D0/li)/3.D0
 
         f10=f10-lyi*mi/li**2.D0*(1.D0/rii-1.D0/ri)/3.D0-
     &                  lxi*g1*(2.D0*li*z(1)**2.D0+mi**2.D0/li)/3.D0
 
         f11=f11+d0*(-s21+t21+(s22-t22)*2.D0/3.D0)+z(1)*(lyi*lxi/
     &             li**2.D0*(1.D0/rii-1.D0/ri)-lyi**2.D0*mi/li*g1)/3.D0

         f12=f12+d0*((s21-t21)*2/3-s22+t22)-z(1)*(lyi*lxi/
     &             li**2.D0*(1.D0/rii-1.D0/ri)+lxi**2.D0*mi/li*g1)/3.D0
 
         f13=f13-(lxi**2.D0-lyi**2.D0)/li**2.D0*(1.D0/rii-1.D0/ri)/6.D0
     &             +2.D0*lxi*lyi*mi/li*g1/6.D0
 
         f16=f16+(lyi/li*h1+(g3-h3)/2.D0
     &           +lyi*(lyi**2.D0-lxi**2.D0)/li**3.D0*h1)/6.D0
     &           +lxi*mi/li**4.D0*(3.D0*lyi**2.D0-lxi**2.D0)
     &           *(1.D0/rii-1.D0/ri)/6.D0+lyi*mi**2.D0*(3.D0*lxi**2.D0
     &           -lyi**2.D0)/li**3.D0*g1/6+z(1)**2.D0*lyi*(lxi**2.D0-
     &           lyi**2.D0)/li*g1/6.D0
 
         f17=f17-(lxi/li*h1-(g4-h4)/2.D0+lxi*(lxi**2-lyi**2)/li**3*h1)
     &       /6.D0+lyi*mi/li**4*(3.D0*lxi**2-lyi**2)*(1.D0/rii-1.D0/ri)
     &       /6.D0-lxi*mi**2*(3.D0*lyi**2-lxi**2)/li**3*g1/6.D0
     &       -z(1)**2*lxi*(lyi**2-lxi**2)/li*g1/6.D0 
 
            !if (isnan(f1))then
            !    write(*,*)i,j,'ai=',a(i),'ri=',ri
            !end if 
                   
   30 continue
! 
      f14=f3-f16-f9
      f15=f4-f17-f10
 
      f2 = -dabs(f2)
      f8 = -dabs(f8)
      f11= -dabs(f11)
      f12= -dabs(f12)
!      
      !write(*,*)'f1=',f1
      !write(*,*)'f2=',f2
      !write(*,*)'f3=',f3
      !write(*,*)'f4=',f4
      !write(*,*)'f5=',f5
      !write(*,*)'f6=',f6
      !write(*,*)'f7=',f7
      !write(*,*)'f8=',f8
      !write(*,*)'f9=',f9
      !write(*,*)'f10=',f10
      !write(*,*)'f11=',f11
      !write(*,*)'f12=',f12
      !write(*,*)'f13=',f13
      !write(*,*)'f14=',f14
      !write(*,*)'f15=',f15
      !write(*,*)'f16=',f16
      !write(*,*)'f17=',f17
!      
!      
      U(1,1)= F5          + ( 3.0D0 - 4.0D0 * V ) * F1
      U(1,2)= F7
      U(1,3)= F3 * localZ
      U(2,1)= F7
      U(2,2)= F6          + ( 3.0D0 - 4.0D0 * V ) * F1
      U(2,3)= F4 * localZ
      U(3,1)= F3 * localZ
      U(3,2)= F4 * localZ
      U(3,3)= F2 * localZ + ( 3.0D0 - 4.0D0 * V ) * F1
!      !
!      
      !write(*,*) 'U(1,1)', U(1,1)
      !write(*,*) 'U(1,2)', U(1,2)
      !write(*,*) 'U(1,3)', U(1,3)
      !write(*,*) 'U(2,1)', U(2,1)
      !write(*,*) 'U(2,2)', U(2,2)
      !write(*,*) 'U(2,3)', U(2,3)
      !write(*,*) 'U(3,1)', U(3,1)
      !write(*,*) 'U(3,2)', U(3,2)
      !write(*,*) 'U(3,3)', U(3,3)
!      
!      
!      
      JI(1,1)=JEvect(1,1)*IEvect(1,1)+JEvect(1,2)*IEvect(1,2)+JEvect(1,3
     &)*IEvect(1,3)
      JI(2,1)=JEvect(2,1)*IEvect(1,1)+JEvect(2,2)*IEvect(1,2)+JEvect(2,3
     &)*IEvect(1,3)
      JI(3,1)=JEvect(3,1)*IEvect(1,1)+JEvect(3,2)*IEvect(1,2)+JEvect(3,3
     &)*IEvect(1,3)
      
      JI(1,2)=JEvect(1,1)*IEvect(2,1)+JEvect(1,2)*IEvect(2,2)+JEvect(1,3
     &)*IEvect(2,3)
      JI(2,2)=JEvect(2,1)*IEvect(2,1)+JEvect(2,2)*IEvect(2,2)+JEvect(2,3
     &)*IEvect(2,3)
      JI(3,2)=JEvect(3,1)*IEvect(2,1)+JEvect(3,2)*IEvect(2,2)+JEvect(3,3
     &)*IEvect(2,3) 
      
      JI(1,3)=JEvect(1,1)*IEvect(3,1)+JEvect(1,2)*IEvect(3,2)+JEvect(1,3
     &)*IEvect(3,3)
      JI(2,3)=JEvect(2,1)*IEvect(3,1)+JEvect(2,2)*IEvect(3,2)+JEvect(2,3
     &)*IEvect(3,3)
      JI(3,3)=JEvect(3,1)*IEvect(3,1)+JEvect(3,2)*IEvect(3,2)+JEvect(3,3
     &)*IEvect(3,3)  
!      
!     
      COEFF(1,1)=cond*(U(1,1)*JI(1,1)+U(2,1)*JI(2,1)+U(3,1)*JI(3,1))
      COEFF(2,1)=cond*(U(1,1)*JI(1,2)+U(2,1)*JI(2,2)+U(3,1)*JI(3,2))
      COEFF(3,1)=cond*(U(1,1)*JI(1,3)+U(2,1)*JI(2,3)+U(3,1)*JI(3,3))
     
      COEFF(1,2)=cond*(U(1,2)*JI(1,1)+U(2,2)*JI(2,1)+U(3,2)*JI(3,1))
      COEFF(2,2)=cond*(U(1,2)*JI(1,2)+U(2,2)*JI(2,2)+U(3,2)*JI(3,2))
      COEFF(3,2)=cond*(U(1,2)*JI(1,3)+U(2,2)*JI(2,3)+U(3,2)*JI(3,3))
     
      COEFF(1,3)=cond*(U(1,3)*JI(1,1)+U(2,3)*JI(2,1)+U(3,3)*JI(3,1))
      COEFF(2,3)=cond*(U(1,3)*JI(1,2)+U(2,3)*JI(2,2)+U(3,3)*JI(3,2))
      COEFF(3,3)=cond*(U(1,3)*JI(1,3)+U(2,3)*JI(2,3)+U(3,3)*JI(3,3))
!     
      !write(*,*)'coeff(1,1)',coeff(1,1)
      !write(*,*)'coeff(2,1)',coeff(2,1)
      !write(*,*)'coeff(3,1)',coeff(3,1)
      !write(*,*)'coeff(1,2)',coeff(1,2)
      !write(*,*)'coeff(2,2)',coeff(2,2)
      !write(*,*)'coeff(3,2)',coeff(3,2)
      !write(*,*)'coeff(1,3)',coeff(1,3)
      !write(*,*)'coeff(2,3)',coeff(2,3)
      !write(*,*)'coeff(3,3)',coeff(3,3)      
   
		  end subroutine diag_dis_fs3d_tri_ana
      
		  
      end module subs_fs3d_tri_ana