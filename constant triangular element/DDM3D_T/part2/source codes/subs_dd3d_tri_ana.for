      module subs_dd3d_tri_ana
      contains 
      
      subroutine diag_stress_dd3d_tri_ana (xx1,yy1,zz1,xx2,yy2,zz2,xx3,
     &yy3,zz3,E,v,IEvect,COEFF)

      real*8 x(5),y(5),z(5),a(5),b(5),pp(4,5),sb(5),sa(5),sbx(5),sax(5),
     &m(5),g(5),h(5),gx(5),hx(5),gy(5),hy(5),sby(5),say(5),gz(5),hz(5),
     &sbz(5),saz(5),gxx(5),hxx(5),tb(5),ta(5),tbx(5),tax(5),sbxx(5),
     &saxx(5),gxy(5),hxy(5),sbxy(5),saxy(5),tby(5),tay(5),sbyy(5),
     &sayy(5),gyy(5),hyy(5),gyz(5),hyz(5),sbyz(5),sayz(5),tbz(5),taz(5),
     &gxz(5),hxz(5),sbxz(5),saxz(5),gzz(5),hzz(5),sbzz(5),sazz(5),
     &gxxx(5),hxxx(5),sbxxx(5),saxxx(5),tbxx(5),taxx(5),gyyy(5),hyyy(5),
     &sbyyy(5),sayyy(5),tbyy(5),tayy(5),gzzz(5),hzzz(5),sbzzz(5),
     &sazzz(5),tbzz(5),tazz(5),gxxy(5),hxxy(5),sbxxy(5),saxxy(5),tbxy(5)
     &,taxy(5),gxxz(5),hxxz(5),sbxxz(5),saxxz(5),tbxz(5),taxz(5),gyyz(5)
     &,hyyz(5),sbyyz(5),sayyz(5),gxyy(5),hxyy(5),sbxyy(5),saxyy(5),
     &gxyz(5),hxyz(5),sbxyz(5),saxyz(5),tbyz(5),tayz(5),gxzz(5),hxzz(5),
     &sbxzz(5),saxzz(5),gyzz(5),hyzz(5),sbyzz(5),sayzz(5),rsb(5),rsa(5),
     &rtb(5),rta(5),rsbx(5),rsax(5),rsaz(5),rsbz(5),rtbx(5),rtax(5),
     &rtbz(5),rtaz(5),rtby(5),rtay(5),rsby(5),rsay(5) 
      real*8 lx,ly,f,fx,fy,fz,fxx,fxy,fyy,fyz,fxz,fzz,fxxx,fyyy,fzzz,
     &fxxy,fxxz,fyyz,fxyy,fxyz,fxzz,fyzz,EE,x1,x2,x3,y1,y2,y3,z1,z2,z3,
     &x4,y4,z4
      real*8::E,v,IEvect(3,3),JEvect(3,3),COEFF(3,3),JI(3,3),trac(3)
      real*8::S(6,3),xx(4),yy(4),zz(4),gpx,gpy,gpz,GG,LOCALZ
      real*8::xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,cons
      integer::i,j
      
       GG=E/(1+v)/2.
       cons=1./4./3.141592653/(1-v)*GG
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
      
      do 10 i=1,4
         do 15 j=1,5
            pp(i,j)=dsqrt((x(j)-x(i))**2+(y(j)-y(i))**2+(z(j)-z(i))**2)
   15    continue
   10 continue
      EE=DBLE(0.0D0)
      f =0.0D0
      fx =0.0D0
      fy =0.0d0
      fz =0.0D0
      fxx =0.0D0
      fxy =0.0D0
      fyy =0.0D0
      fyz =0.0D0
      fxz =0.0D0
      fzz =0.0D0
      fxxx=0.0D0
      fyyy=0.0D0
      fzzz=0.0D0
      fxxy=0.0D0
      fxxz=0.0D0
      fyyz=0.0D0
      fxyy=0.0D0
      fxyz=0.0D0
      fxzz=0.0D0
      fyzz=0.0D0
      do 20 i=2,4
         a(i)=((x(i)-x(1))*(x(i+1)-x(i))+(y(i)-y(1))*(y(i+1)-y(i)))/pp(
     &   i,i+1)
 
         b(i)=((x(i+1)-x(1))*(x(i+1)-x(i))+(y(i+1)-y(1))*(y(i+1)-y(i)))
     &   /pp(i,i+1)
 
         lx=x(i)-x(i+1) 
         ly=y(i)-y(i+1)   
         m(i)=(x(i)-x(1))*(y(i+1)-y(1))-(x(i+1)-x(1))*(y(i)-y(1))
         
         
c         write (210,*)'bi=',b(i),'pp=',pp(1,i+1)
         
         g(i)=dlog(pp(1,i+1)+b(i))
         h(i)=dlog(pp(1,i)+a(i))
         
         gx(i)=((x(1)-x(i+1))/pp(1,i+1)+lx/pp(i,i+1))/(pp(1,i+1)+b(i))
         hx(i)=((x(1)-x(i))/pp(1,i)+lx/pp(i,i+1))/(pp(1,i)+a(i))
 
         gy(i)=((y(1)-y(i+1))/pp(1,i+1)+ly/pp(i,i+1))/(pp(1,i+1)+b(i))
         hy(i)=((y(1)-y(i))/pp(1,i)+ly/pp(i,i+1))/(pp(1,i)+a(i))
 
         gz(i)=(z(1)/pp(1,i+1))/(pp(1,i+1)+b(i))
         hz(i)=(z(1)/pp(1,i))/(pp(1,i)+a(i))
 
         gxx(i)=-gx(i)*gx(i)+(((y(1)-y(i+1))**2+z(1)**2)/pp(1,i+1)**3
     &    )/(pp(1,i+1)+b(i))
         hxx(i)=-hx(i)*hx(i)+(((y(1)-y(i))**2+z(1)**2)/pp(1,i)**3
     &    )/(pp(1,i)+a(i))
 
         gxy(i)=-(gx(i)*((y(1)-y(i+1))/pp(1,i+1)+ly/pp(i,i+1))+
     &   (x(1)-x(i+1))*(y(1)-y(i+1))/pp(1,i+1)**3)/(pp(1,i+1)+b(i))
         hxy(i)=-(hx(i)*((y(1)-y(i))/pp(1,i)+ly/pp(i,i+1))+
     &   (x(1)-x(i))*(y(1)-y(i))/pp(1,i)**3)/(pp(1,i)+a(i))
 
         gyy(i)=-gy(i)*gy(i)+(((x(1)-x(i+1))**2+z(1)**2)/pp(1,i+1)**3
     &    )/(pp(1,i+1)+b(i))
         hyy(i)=-hy(i)*hy(i)+(((x(1)-x(i))**2+z(1)**2)/pp(1,i)**3)
     &   /(pp(1,i)+a(i))
         gyz(i)=-(gy(i)*z(1)/pp(1,i+1)+(y(1)-y(i+1))*z(1)/pp(1,i+1)**3)
     &   /(pp(1,i+1)+b(i))
         hyz(i)=-(hy(i)*z(1)/pp(1,i)+(y(1)-y(i))*z(1)/pp(1,i)**3)
     &   /(pp(1,i)+a(i))
         gxz(i)=-(gx(i)*z(1)/pp(1,i+1)+(x(1)-x(i+1))*z(1)/pp(1,i+1)**3)
     &   /(pp(1,i+1)+b(i))
         hxz(i)=-(hx(i)*z(1)/pp(1,i)+(x(1)-x(i))*z(1)/pp(1,i)**3)
     &   /(pp(1,i)+a(i))
         gzz(i)=-gz(i)*gz(i)+(((x(1)-x(i+1))**2+(y(1)-y(i+1))**2)
     &   /pp(1,i+1)**3)/(pp(1,i+1)+b(i))
         hzz(i)=-hz(i)*hz(i)+(((x(1)-x(i))**2+(y(1)-y(i))**2)
     &   /pp(1,i)**3)/(pp(1,i)+a(i))
         gxxx(i)=2.*gx(i)**3-3.*gx(i)/(pp(1,i+1)+b(i))*((y(1)-y(i+1))**2
     &    +z(1)**2)/pp(1,i+1)**3-3.*((x(1)-x(i+1))*((y(1)-y(i+1))**2
     &    +z(1)**2)/pp(1,i+1)**5)/(pp(1,i+1)+b(i))
         hxxx(i)=2.*hx(i)**3-3.*hx(i)/(pp(1,i)+a(i))*((y(1)-y(i))**2
     &    +z(1)**2)/pp(1,i)**3-3.*((x(1)-x(i))*((y(1)-y(i))**2+z(1)**2)
     &   /pp(1,i)**5)/(pp(1,i)+a(i))
         gyyy(i)=2.*gy(i)**3-3.*gy(i)/(pp(1,i+1)+b(i))*((x(1)-x(i+1))**2
     &    +z(1)**2)/pp(1,i+1)**3-3.*((y(1)-y(i+1))*((x(1)-x(i+1))**2
     &    +z(1)**2)/pp(1,i+1)**5)/(pp(1,i+1)+b(i))
         hyyy(i)=2.*hy(i)**3-3.*hy(i)/(pp(1,i)+a(i))*((x(1)-x(i))**2
     &    +z(1)**2)/pp(1,i)**3-3.*((y(1)-y(i))*((x(1)-x(i))**2
     &    +z(1)**2)/pp(1,i)**5)/(pp(1,i)+a(i))
         gzzz(i)=2.*gz(i)**3-3.*gz(i)/(pp(1,i+1)+b(i))*((x(1)-x(i+1))**2
     &    +(y(1)-y(i+1))**2)/pp(1,i+1)**3-3.*(z(1)*((x(1)-x(i+1))**2
     &    +(y(1)-y(i+1))**2)/pp(1,i+1)**5)/(pp(1,i+1)+b(i))
         hzzz(i)=2.*hz(i)**3-3.*hz(i)/(pp(1,i)+a(i))*((x(1)-x(i))**2
     &   +(y(1)-y(i))**2)/pp(1,i)**3-3.*(z(1)*((x(1)-x(i))**2+(y(1)-y(i)
     &   )**2)/pp(1,i)**5)/(pp(1,i)+a(i))
         gxxy(i)=2.*gx(i)**2*gy(i)+2.*gx(i)*(x(1)-x(i+1))*(y(1)-y(i+1))
     &   /pp(1,i+1)**3/(pp(1,i+1)+b(i))-gy(i)*((y(1)-y(i+1))**2+z(1)**2
     &    )/pp(1,i+1)**3/(pp(1,i+1)+b(i))+(y(1)-y(i+1))*(2.*pp(1,i+1)**2
     &    -3.*((y(1)-y(i+1))**2+z(1)**2))/pp(1,i+1)**5/(pp(1,i+1)+b(i))
         hxxy(i)=2.*hx(i)**2*hy(i)+2.*hx(i)*(x(1)-x(i))*(y(1)-y(i))
     &   /pp(1,i)**3/(pp(1,i)+a(i))-hy(i)*((y(1)-y(i))**2+z(1)**2)
     &   /pp(1,i)**3/(pp(1,i)+a(i))+(y(1)-y(i))*(2.*pp(1,i)**2-3.*(
     &   (y(1)-y(i))**2+z(1)**2))/pp(1,i)**5/(pp(1,i)+a(i))
         gxxz(i)=2.*gx(i)**2*gz(i)+2.*gx(i)*z(1)*(x(1)-x(i+1))/
     &   pp(1,i+1)**3/(pp(1,i+1)+b(i))-gz(i)*((y(1)-y(i+1))**2+z(1)**2)
     &   /pp(1,i+1)**3/(pp(1,i+1)+b(i))+z(1)*(2.*pp(1,i+1)**2-3.*
     &   ((y(1)-y(i+1))**2+z(1)**2))/pp(1,i+1)**5/(pp(1,i+1)+b(i))
         hxxz(i)=2.*hx(i)**2*hz(i)+2.*hx(i)*z(1)*(x(1)-x(i))/pp(1,i)**3
     &    /(pp(1,i)+a(i))-hz(i)*((y(1)-y(i))**2+z(1)**2)/pp(1,i)**3/
     &   (pp(1,i)+a(i))+z(1)*(2.*pp(1,i)**2-3.*((y(1)-y(i))**2+z(1)**2))
     &   /pp(1,i)**5/(pp(1,i)+a(i))
         gyyz(i)=2.*gy(i)**2*gz(i)+2.*gy(i)*z(1)*(y(1)-y(i+1))/
     &   pp(1,i+1)**3/(pp(1,i+1)+b(i))-gz(i)*((x(1)-x(i+1))**2+z(1)**2)
     &   /pp(1,i+1)**3/(pp(1,i+1)+b(i))+z(1)*(2*pp(1,i+1)**2-3.*
     &   ((x(1)-x(i+1))**2+z(1)**2))/pp(1,i+1)**5/(pp(1,i+1)+b(i))
         hyyz(i)=2.*hy(i)**2*hz(i)+2.*hy(i)*z(1)*(y(1)-y(i))/pp(1,i)**3
     &    /(pp(1,i)+a(i))-hz(i)*((x(1)-x(i))**2+z(1)**2)/pp(1,i)**3/
     &   (pp(1,i)+a(i))+z(1)*(2.*pp(1,i)**2-3.*((x(1)-x(i))**2+z(1)**2))
     &   /pp(1,i)**5/(pp(1,i)+a(i))
         gxyy(i)=2.*gx(i)*gy(i)**2+2.*gy(i)*(x(1)-x(i+1))*(y(1)-y(i+1))/
     &   pp(1,i+1)**3/(pp(1,i+1)+b(i))-gx(i)*((x(1)-x(i+1))**2+z(1)**2
     &    )/pp(1,i+1)**3/(pp(1,i+1)+b(i))+(x(1)-x(i+1))*(2.*pp(1,i+1)**2
     &    -3.*((x(1)-x(i+1))**2+z(1)**2))/pp(1,i+1)**5/(pp(1,i+1)+b(i))
         hxyy(i)=2.*hx(i)*hy(i)**2+2.*hy(i)*(x(1)-x(i))*(y(1)-y(i))/
     &   pp(1,i)**3/(pp(1,i)+a(i))-hx(i)*((x(1)-x(i))**2+z(1)**2)
     &   /pp(1,i)**3/(pp(1,i)+a(i))+(x(1)-x(i))*(2.*pp(1,i)**2
     &    -3.*((x(1)-x(i))**2+z(1)**2))/pp(1,i)**5/(pp(1,i)+a(i))
         gxyz(i)=2.*gx(i)*gy(i)*gz(i)+z(1)*(((x(1)-x(i+1))*(y(1)-y(i+1))
     &   /pp(1,i+1)+(y(1)-y(i+1))*((x(1)-x(i+1))/pp(1,i+1)+lx/
     &   pp(i,i+1))+(x(1)-x(i+1))*((y(1)-y(i+1))/pp(1,i+1)+ly/pp(i,i+1)
     &   ))/(pp(1,i+1)+b(i))**2/pp(1,i+1)**3+3.*(x(1)-x(i+1))*(y(1)-
     &   y(i+1))/pp(1,i+1)**5/(pp(1,i+1)+b(i)))
         hxyz(i)=2.*hx(i)*hy(i)*hz(i)+z(1)*(((x(1)-x(i))*(y(1)-y(i))
     &   /pp(1,i)+(y(1)-y(i))*((x(1)-x(i))/pp(1,i)+lx/pp(i,i+1))
     &   +(x(1)-x(i))*((y(1)-y(i))/pp(1,i)+ly/pp(i,i+1)))
     &   /(pp(1,i)+a(i))**2/pp(1,i)**3+3.*(x(1)-x(i))*(y(1)-y(i))
     &   /pp(1,i)**5/(pp(1,i)+a(i)))
         gxzz(i)=2.*gx(i)*gz(i)**2+2.*gz(i)*(x(1)-x(i+1))*z(1)/pp(1,i+1)
     &   **3
     &    /(pp(1,i+1)+b(i))-gx(i)*((x(1)-x(i+1))**2+(y(1)-y(i+1))**2)/
     &   pp(1,i+1)**3/(pp(1,i+1)+b(i))+(x(1)-x(i+1))*(2.*pp(1,i+1)**2
     &    -3.*((x(1)-x(i+1))**2+(y(1)-y(i+1))**2))/pp(1,i+1)**5/
     &   (pp(1,i+1)+b(i)) 
         hxzz(i)=2.*hx(i)*hz(i)**2+2.*hz(i)*(x(1)-x(i))*z(1)/pp(1,i)**3
     &    /(pp(1,i)+a(i))-hx(i)*((x(1)-x(i))**2+(y(1)-y(i))**2)/
     &   pp(1,i)**3/(pp(1,i)+a(i))+(x(1)-x(i))*(2.*pp(1,i)**2
     &    -3.*((x(1)-x(i))**2+(y(1)-y(i))**2))/pp(1,i)**5/(pp(1,i)+a(i))
 
         gyzz(i)=2.*gy(i)*gz(i)**2+2.*gz(i)*(y(1)-y(i+1))*z(1)/pp(1,i+1)
     &   **3
     &    /(pp(1,i+1)+b(i))-gy(i)*((x(1)-x(i+1))**2+(y(1)-y(i+1))**2)/
     &   pp(1,i+1)**3/(pp(1,i+1)+b(i))+(y(1)-y(i+1))*(2.*pp(1,i+1)**2
     &    -3.*((x(1)-x(i+1))**2+(y(1)-y(i+1))**2))/pp(1,i+1)**5/
     &   (pp(1,i+1)+b(i))
         hyzz(i)=2.*hy(i)*hz(i)**2+2.*hz(i)*(y(1)-y(i))*z(1)/pp(1,i)**3
     &    /(pp(1,i)+a(i))-hy(i)*((x(1)-x(i))**2+(y(1)-y(i))**2)/
     &   pp(1,i)**3/(pp(1,i)+a(i))+(y(1)-y(i))*(2.*pp(1,i)**2
     &    -3.*((x(1)-x(i))**2+(y(1)-y(i))**2))/pp(1,i)**5/(pp(1,i)+a(i))
 
         if (((z(1).ge.-0.1D-10).and.(z(1).le.0.1D-10)).and.
     &   ((lx.le.-0.1D-10).or.(lx.ge.0.1D-10))) then
            tt=1
            rsb(i)=(pp(i,i+1)/(-lx))*(y(i+1)-y(1)+b(i)+(1-ly
     &       /pp(i,i+1))*pp(1,i+1))
            rsa(i)=(pp(i,i+1)/(-lx))*(y(i)-y(1)+a(i)+(1-ly/pp(i,i+1))
     &      *pp(1,i))
 
            rsbx(i)=pp(i,i+1)/(-lx)*(lx/pp(i,i+1)+(1-ly/pp(i,i+1))
     &      *(-x(i+1)+x(1))/pp(1,i+1))
            rsax(i)=pp(i,i+1)/(-lx)*(lx/pp(i,i+1)+(1-ly/pp(i,i+1))
     &      *(-x(i)+x(1))/pp(1,i))
            rsbz(i)=pp(i,i+1)/(-lx)*(1-ly/pp(i,i+1))*z(1)/pp(1,i+1)
            rsaz(i)=pp(i,i+1)/(-lx)*(1-ly/pp(i,i+1))*z(1)/pp(1,i)
            rsby(i)=pp(i,i+1)/(-lx)*(-1+ly/pp(i,i+1)+(1-ly/pp(i,i+1))
     &      *(-y(i+1)+y(1))/pp(1,i+1))
            rsay(i)=pp(i,i+1)/(-lx)*(-1+ly/pp(i,i+1)+(1-ly/pp(i,i+1))
     &      *(-y(i)+y(1))/pp(1,i))
 
 
            rtb(i)=rsb(i)**2+z(1)**2
            rta(i)=rsa(i)**2+z(1)**2
            rtbx(i)=2*rsb(i)*rsbx(i)
            rtax(i)=2*rsa(i)*rsax(i)
            rtbz(i)=2*(rsb(i)*rsbz(i)+z(1))
            rtaz(i)=2*(rsa(i)*rsaz(i)+z(1))
            rtby(i)=2*rsb(i)*rsby(i)
            rtay(i)=2*rsa(i)*rsay(i)
 
            fzz=fzz+(m(i)*(gzz(i)-hzz(i)))/pp(i,i+1)+4*tt*(rsb(i)/rtb(i)
     &      -rsa(i)/rta(i))
            fxzz=fxzz+(ly*(gzz(i)-hzz(i))+m(i)*(gxzz(i)-hxzz(i)))
     &      /pp(i,i+1)-2*tt*(2*((rsb(i)*rtbx(i)-rsbx(i)*rtb(i))/rtb(i)
     &      **2
     &       -(rsa(i)*rtax(i)-rsax(i)*rta(i))/rta(i)**2)-((rsbz(i)-
     &      rsb(i))
     &      *rtbz(i)/rtb(i)**2-(rsaz(i)-rsa(i))*rtaz(i)/rta(i)**2))
            fyzz=fyzz+(-lx*(gzz(i)-hzz(i))+m(i)*(gyzz(i)-hyzz(i)))
     &      /pp(i,i+1)-2*tt*(2*((rsb(i)*rtby(i)-rsby(i)*rtb(i))/rtb(i)
     &      **2
     &       -(rsa(i)*rtay(i)-rsay(i)*rta(i))/rta(i)**2)-((rsbz(i)-
     &      rsb(i))
     &      *rtbz(i)/rtb(i)**2-(rsaz(i)-rsa(i))*rtaz(i)/rta(i)**2))
         else
            fzz=fzz+(m(i)*(gzz(i)-hzz(i)))/pp(i,i+1)
            fxzz=fxzz+(ly*(gzz(i)-hzz(i))+m(i)*(gxzz(i)-hxzz(i)))/pp(i,
     &      i+1)
            fyzz=fyzz+(-lx*(gzz(i)-hzz(i))+m(i)*(gyzz(i)-hyzz(i)))
     &      /pp(i,i+1)
         end if
 
         f=f+m(i)/pp(i,i+1)*(g(i)-h(i))
         fx=fx+(ly*(g(i)-h(i))+m(i)*(gx(i)-hx(i)))/pp(i,i+1)
              
         fy=fy+(-lx*(g(i)-h(i))+m(i)*(gy(i)-hy(i)))/pp(i,i+1)
         fz=fz+(m(i)*(gz(i)-hz(i)))/pp(i,i+1)
         fxx=fxx+(2*ly*(gx(i)-hx(i))+m(i)*(gxx(i)-hxx(i)))/pp(i,i+1)
         fxy=fxy+(-lx*(gx(i)-hx(i))+ly*(gy(i)-hy(i))+m(i)*(gxy(i)-
     &   hxy(i)))/pp(i,i+1) 
         fyy=fyy+(-2*lx*(gy(i)-hy(i))+m(i)*(gyy(i)-hyy(i)))/pp(i,i+1)
         fyz=fyz+(-lx*(gz(i)-hz(i))+m(i)*(gyz(i)-hyz(i)))/pp(i,i+1)
         fxz=fxz+(ly*(gz(i)-hz(i))+m(i)*(gxz(i)-hxz(i)))/pp(i,i+1)
         fxxx=fxxx+(3*ly*(gxx(i)-hxx(i))+m(i)*(gxxx(i)-hxxx(i)))
     &   /pp(i,i+1)
         fyyy=fyyy+(-3*lx*(gyy(i)-hyy(i))+m(i)*(gyyy(i)-hyyy(i)))
     &   /pp(i,i+1)
         fzzz=fzzz+m(i)*(gzzz(i)-hzzz(i))/pp(i,i+1)
         fxxy=fxxy+(2*ly*(gxy(i)-hxy(i))-lx*(gxx(i)-hxx(i))+m(i)*
     &   (gxxy(i)-hxxy(i)))/pp(i,i+1)
         fxxz=fxxz+(2*ly*(gxz(i)-hxz(i))+m(i)*(gxxz(i)-hxxz(i)))
     &   /pp(i,i+1)
         fyyz=fyyz+(-2*lx*(gyz(i)-hyz(i))+m(i)*(gyyz(i)-hyyz(i)))
     &   /pp(i,i+1)
         fxyy=fxyy+(-2*lx*(gxy(i)-hxy(i))+ly*(gyy(i)-hyy(i))+m(i)*
     &   (gxyy(i)-hxyy(i)))/pp(i,i+1)
         fxyz=fxyz+(-lx*(gxz(i)-hxz(i))+ly*(gyz(i)-hyz(i))+m(i)*
     &   (gxyz(i)-hxyz(i)))/pp(i,i+1)
 
         if ((z(1).le.-0.1D-10).or.(z(1).ge.0.1D-10)) then
            if ((lx.le.-0.1D-10).or.(lx.ge.0.1D-10)) then
               tt=1
               sb(i)=(pp(i,i+1)/((-lx)*z(1)))*(y(i+1)-y(1)+b(i)+(1-ly
     &          /pp(i,i+1))*pp(1,i+1))
               sa(i)=(pp(i,i+1)/((-lx)*z(1)))*(y(i)-y(1)+a(i)+(1-
     &         ly/pp(i,i+1))
     &         *pp(1,i))
               sbx(i)=pp(i,i+1)/((-lx)*z(1))*(lx/pp(i,i+1)+(1-ly/pp(i,i+
     &         1))
     &         *(-x(i+1)+x(1))/pp(1,i+1))
               sax(i)=pp(i,i+1)/((-lx)*z(1))*(lx/pp(i,i+1)+(1-ly/pp(i,i+
     &         1))
     &         *(-x(i)+x(1))/pp(1,i))
 
               sby(i)=pp(i,i+1)/((-lx)*z(1))*(-1+ly/pp(i,i+1)+(1-
     &         ly/pp(i,i+1))
     &         *(-y(i+1)+y(1))/pp(1,i+1))
               say(i)=pp(i,i+1)/((-lx)*z(1))*(-1+ly/pp(i,i+1)+(1-
     &         ly/pp(i,i+1))
     &         *(-y(i)+y(1))/pp(1,i))
 
               sbz(i)=sb(i)*(-1)/z(1)+pp(i,i+1)/((-lx)*z(1))*((1-
     &         ly/pp(i,i+1))
     &         *z(1)/pp(1,i+1))
               saz(i)=sa(i)*(-1)/z(1)+pp(i,i+1)/((-lx)*z(1))*((1-
     &         ly/pp(i,i+1))
     &         *z(1)/pp(1,i))
               sbxx(i)=pp(i,i+1)/((-lx)*z(1))*((1-ly/pp(i,i+1))*(pp(1,i+
     &         1)
     &         **2-(x(1)-x(i+1))**2)/pp(1,i+1)**3)
               saxx(i)=pp(i,i+1)/((-lx)*z(1))*((1-ly/pp(i,i+1))*(pp(1,i)
     &         **2-(x(1)-x(i))**2)/pp(1,i)**3)
               sbxy(i)=pp(i,i+1)/((-lx)*z(1))*((1-ly/pp(i,i+1))*(x(i+1)-
     &         x(1))
     &         *(y(1)-y(i+1))/pp(1,i+1)**3)
               saxy(i)=pp(i,i+1)/((-lx)*z(1))*((1-ly/pp(i,i+1))*(x(i)-
     &         x(1))
     &         *(y(1)-y(i))/pp(1,i)**3)
               sbyy(i)=pp(i,i+1)/((-lx)*z(1))*((1-ly/pp(i,i+1))*(pp(1,i+
     &         1)
     &         **2-(y(1)-y(i+1))**2)/pp(1,i+1)**3)
               sayy(i)=pp(i,i+1)/((-lx)*z(1))*((1-ly/pp(i,i+1))*(pp(1,i)
     &         **2-(y(1)-Y(i))**2)/pp(1,i)**3)
               sbyz(i)=pp(i,i+1)/(-lx)*((-1)/z(1)**2*(-1+ly/pp(i,i+1)
     &         +(1-ly/pp(i,i+1))*(y(1)-y(i+1))/pp(1,i+1))-(1-ly/pp(i,i+
     &         1))
     &         *(y(1)-y(i+1))/pp(1,i+1)**3)
               sayz(i)=pp(i,i+1)/(-lx)*((-1)/z(1)**2*(-1+ly/pp(i,i+1)
     &         +(1-ly/pp(i,i+1))*(y(1)-y(i))/pp(1,i))-(1-ly/pp(i,i+1))
     &         *(y(1)-y(i))/pp(1,i)**3)
               sbxz(i)=pp(i,i+1)/(-lx)*((-1)/z(1)**2*(lx/pp(i,i+1)+(1-
     &         ly/
     &         pp(i,i+1))*(x(1)-x(i+1))/pp(1,i+1))-(1-ly/pp(i,i+1))*
     &         (x(1)-x(i+1))/pp(1,i+1)**3)
               saxz(i)=pp(i,i+1)/(-lx)*((-1)/z(1)**2*(lx/pp(i,i+1)+(1-
     &         ly/
     &         pp(i,i+1))*(x(1)-x(i))/pp(1,i))-(1-ly/pp(i,i+1))*
     &         (x(1)-x(i))/pp(1,i)**3)
 
               sbzz(i)=pp(i,i+1)/(-lx)*(2/z(1)**3*(y(i+1)-y(1)+b(i)+(1-
     &         ly/
     &         pp(i,i+1))*pp(1,i+1))-1/z(1)*((1-ly/pp(i,i+1))/pp(1,i+1))
     &         -z(1)*((1-ly/pp(i,i+1))/pp(1,i+1)**3))
               sazz(i)=pp(i,i+1)/(-lx)*(2/z(1)**3*(y(i)-y(1)+a(i)+(1-ly/
     &         pp(i,i+1))*pp(1,i))-1/z(1)*((1-ly/pp(i,i+1))/pp(1,i))
     &         -z(1)*((1-ly/pp(i,i+1))/pp(1,i)**3))
               sbxxx(i)=pp(i,i+1)/(-lx)/z(1)*((1-ly/pp(i,i+1))*(-3*
     &         (x(1)-X(i+1))*(pp(1,i+1)**2-(x(1)-x(i+1))**2)/pp(1,i+1)
     &         **5))
               saxxx(i)=pp(i,i+1)/(-lx)/z(1)*((1-ly/pp(i,i+1))*(-3*
     &         (x(1)-x(i))*(pp(1,i)**2-(x(1)-x(i))**2)/pp(1,i)**5))
               sbyyy(i)=pp(i,i+1)/(-lx)/z(1)*((1-ly/pp(i,i+1))*(-3*
     &         (y(1)-y(i+1))*(pp(1,i+1)**2-(y(1)-y(i+1))**2)/pp(1,i+1)
     &         **5))
               sayyy(i)=pp(i,i+1)/(-lx)/z(1)*((1-ly/pp(i,i+1))*(-3*
     &         (y(1)-y(i))*(pp(1,i)**2-(y(1)-y(i))**2)/pp(1,i)**5))
               sbzzz(i)=pp(i,i+1)/(-lx)*((-6)/z(1)**4*(-y(1)+y(i+1)+b(i)
     &         +(1-ly/pp(i,i+1))*pp(1,i+1))+3/z(1)**2*((1-ly/pp(i,i+1))
     &         /pp(1,i+1))+3*z(1)**2*((1-ly/pp(i,i+1))/pp(1,i+1)**5))
               sazzz(i)=pp(i,i+1)/(-lx)*((-6)/z(1)**4*(-y(1)+y(i)+a(i)
     &         +(1-ly/pp(i,i+1))*pp(1,i))+3/z(1)**2*((1-ly/pp(i,i+1))
     &         /pp(1,i)
     &         )+3*z(1)**2*((1-ly/pp(i,i+1))/pp(1,i)**5))
               sbxxy(i)=pp(i,i+1)/(-lx)/z(1)*((1-ly/pp(i,i+1))*((y(1)-
     &         y(i+1)
     &         )*(-pp(1,i+1)**2+3*(x(1)-x(i+1))**2)/pp(1,i+1)**5))
               saxxy(i)=pp(i,i+1)/(-lx)/z(1)*((1-ly/pp(i,i+1))*((y(1)-
     &         y(i))
     &         *(-pp(1,i)**2+3*(x(1)-x(i))**2)/pp(1,i)**5))
               sbxxz(i)=pp(i,i+1)/(-lx)*((-1)/z(1)**2*(1-ly/pp(i,i+1))*
     &         (pp(1,i+1)**2-(x(1)-x(i+1))**2)/pp(1,i+1)**3+(1-ly/pp(i,
     &         i+1))
     &         *(-pp(1,i+1)**2+3*(x(1)-x(i+1))**2)/pp(1,i+1)**5)
               saxxz(i)=pp(i,i+1)/(-lx)*((-1)/z(1)**2*(1-ly/pp(i,i+1))*
     &         (pp(1,i)**2-(x(1)-x(i))**2)/pp(1,i)**3+(1-ly/pp(i,i+1))*
     &         (-pp(1,i)**2+3*(x(1)-x(i))**2)/pp(1,i)**5)
               sbyyz(i)=pp(i,i+1)/(-lx)*((-1)/z(1)**2*(1-ly/pp(i,i+1))*
     &         (pp(1,i+1)**2-(y(1)-y(i+1))**2)/pp(1,i+1)**3+(1-ly/pp(i,
     &         i+1))
     &         *(-pp(1,i+1)**2+3*(y(1)-y(i+1))**2)/pp(1,i+1)**5)
               sayyz(i)=pp(i,i+1)/(-lx)*((-1)/z(1)**2*(1-ly/pp(i,i+1))*
     &         (pp(1,i)**2-(y(1)-y(i))**2)/pp(1,i)**3+(1-ly/pp(i,i+1))*
     &         (-pp(1,i)**2+3*(y(1)-y(i))**2)/pp(1,i)**5)
 
               sbxyy(i)=pp(i,i+1)/(-lx)/z(1)*(1-ly/pp(i,i+1))*(-(x(1)-
     &         x(i+1))
     &         *(pp(1,i+1)**2-3*(y(1)-y(i+1))**2))/pp(1,i+1)**5
               saxyy(i)=pp(i,i+1)/(-lx)/z(1)*(1-ly/pp(i,i+1))*(-(x(1)-
     &         x(i))
     &         *(pp(1,i)**2-3*(y(1)-y(i))**2))/pp(1,i)**5
               sbxyz(i)=pp(i,i+1)/(-lx)*((-1)/z(1)**2*(1-ly/pp(i,i+1))*
     &         (-(x(1)-x(i+1))*(y(1)-y(i+1))/pp(1,i+1)**3)+(1-ly/pp(i,i+
     &         1))
     &         *3*(x(1)-x(i+1))*(y(1)-y(i+1))/pp(1,i+1)**5)
               saxyz(i)=pp(i,i+1)/(-lx)*((-1)/z(1)**2*(1-ly/pp(i,i+1))*
     &         (-(x(1)-x(i))*(y(1)-y(i))/pp(1,i)**3)+(1-ly/pp(i,i+1))
     &         *3*(x(1)-x(i))*(y(1)-y(i))/pp(1,i)**5)
 
               sbxzz(i)=pp(i,i+1)/(-lx)*(2/z(1)**3*(lx/pp(i,i+1)+(1-ly/
     &         pp(i,i+1))*(x(1)-x(i+1))/pp(1,i+1))+1/z(1)*(1-ly/pp(i,i+
     &         1))
     &         *(x(1)-x(i+1))/pp(1,i+1)**3+z(1)*(1-ly/pp(i,i+1))*3*
     &         (x(1)-x(i+1))/pp(1,i+1)**5)
               saxzz(i)=pp(i,i+1)/(-lx)*(2/z(1)**3*(lx/pp(i,i+1)+(1-ly/
     &         pp(i,i+1))*(x(1)-x(i))/pp(1,i))+1/z(1)*(1-ly/pp(i,i+1))
     &         *(x(1)-x(i))/pp(1,i)**3+z(1)*(1-ly/pp(i,i+1))*3*
     &         (x(1)-x(i))/pp(1,i)**5)
 
               sbyzz(i)=pp(i,i+1)/(-lx)*(2/z(1)**3*(-1+ly/pp(i,i+1)+(1-
     &         ly/
     &         pp(i,i+1))*(y(1)-y(i+1))/pp(1,i+1))+1/z(1)*(1-ly/pp(i,i+
     &         1))
     &         *(y(1)-y(i+1))/pp(1,i+1)**3+z(1)*(1-ly/pp(i,i+1))*3*
     &         (y(1)-y(i+1))/pp(1,i+1)**5)
               sayzz(i)=pp(i,i+1)/(-lx)*(2/z(1)**3*(-1+ly/pp(i,i+1)+(1-
     &         ly/
     &         pp(i,i+1))*(y(1)-y(i))/pp(1,i))+1/z(1)*(1-ly/pp(i,i+1))
     &         *(y(1)-y(i))/pp(1,i)**3+z(1)*(1-ly/pp(i,i+1))*3*
     &         (y(1)-y(i))/pp(1,i)**5)
               tb(i)=sb(i)**2+1
               ta(i)=sa(i)**2+1
               tbx(i)=2*sb(i)*sbx(i)
               tax(i)=2*sa(i)*sax(i)
               tby(i)=2*sb(i)*sby(i)
               tay(i)=2*sa(i)*say(i)
               tbz(i)=2*sb(i)*sbz(i)
               taz(i)=2*sa(i)*saz(i)
               tbxx(i)=2*(sbx(i)**2+sb(i)*sbxx(i))
               taxx(i)=2*(sax(i)**2+sa(i)*saxx(i))
               tbyy(i)=2*(sby(i)**2+sb(i)*sbyy(i))
               tayy(i)=2*(say(i)**2+sa(i)*sayy(i))
               tbzz(i)=2*(sbz(i)**2+sb(i)*sbzz(i))
               tazz(i)=2*(saz(i)**2+sa(i)*sazz(i))
               tbxy(i)=2*(sby(i)*sbx(i)+sb(i)*sbxy(i))
               taxy(i)=2*(say(i)*sax(i)+sa(i)*saxy(i))
               tbxz(i)=2*(sbz(i)*sbx(i)+sb(i)*sbxz(i))
               taxz(i)=2*(saz(i)*sax(i)+sa(i)*saxz(i))
               tbyz(i)=2*(sbz(i)*sby(i)+sb(i)*sbyz(i))
               tayz(i)=2*(saz(i)*say(i)+sa(i)*sayz(i))
               f=f-2*tt*z(1)*(datan(sb(i))-datan(sa(i)))
               fx=fx-2*tt*z(1)*(sbx(i)/tb(i)-sax(i)/ta(i))
               fy=fy-2*tt*z(1)*(sby(i)/tb(i)-say(i)/ta(i))
               fz=fz-2*tt*(datan(sb(i))-datan(sa(i))+z(1)*(sbz(i)/tb(i)
     &         -saz(i)/ta(i)))
               fxx=fxx-2*tt*z(1)*((sbxx(i)*tb(i)-sbx(i)*tbx(i))/tb(i)**2
     &          -(saxx(i)*ta(i)-sax(i)*tax(i))/ta(i)**2)
               fxy=fxy-2*tt*z(1)*((sbxy(i)*tb(i)-sbx(i)*tby(i))/tb(i)**2
     &          -(saxy(i)*ta(i)-sax(i)*tay(i))/ta(i)**2)
               fyy=fyy-2*tt*z(1)*((sbyy(i)*tb(i)-sby(i)*tby(i))/tb(i)**2
     &          -(sayy(i)*ta(i)-say(i)*tay(i))/ta(i)**2)
               fyz=fyz-2*tt*(sby(i)/tb(i)-say(i)/ta(i)+z(1)*((sbyz(i)
     &         *tb(i)
     &         -sby(i)*tbz(i))/tb(i)**2-(sayz(i)*ta(i)-say(i)*taz(i))
     &         /ta(i)
     &         **2))
               fxz=fxz-2*tt*(sbx(i)/tb(i)-sax(i)/ta(i)+z(1)*((sbxz(i)
     &         *tb(i)
     &         -sbx(i)*tbz(i))/tb(i)**2-(saxz(i)*ta(i)-sax(i)*taz(i))
     &         /ta(i)
     &         **2))
               fzz=fzz-2*tt*(2*(sbz(i)/tb(i)-saz(i)/ta(i))+z(1)
     &         *((sbzz(i)*
     &         tb(i)-sbz(i)*tbz(i))/tb(i)**2-(sazz(i)*ta(i)-saz(i)
     &         *taz(i))
     &         /ta(i)**2))
               fxxx=fxxx-2*tt*z(1)*(((sbxxx(i)*tb(i)-2*sbxx(i)*tbx(i)-
     &         sbx(i)
     &         *tbxx(i))*tb(i)+2*sbx(i)*tbx(i)**2)/tb(i)**3-((saxxx(i)
     &         *ta(i)
     &         -2*saxx(i)*tax(i)-sax(i)*taxx(i))*ta(i)+2*sax(i)*tax(i)
     &         **2)
     &         /ta(i)**3)
               fyyy=fyyy-2*tt*z(1)*(((sbyyy(i)*tb(i)-2*sbyy(i)*tby(i)-
     &         sby(i)
     &         *tbyy(i))*tb(i)+2*sby(i)*tby(i)**2)/tb(i)**3-((sayyy(i)
     &         *ta(i)
     &         -2*sayy(i)*tay(i)-say(i)*tayy(i))*ta(i)+2*say(i)*tay(i)
     &         **2)
     &         /ta(i)**3)
               fzzz=fzzz-2*tt*(3*((sbzz(i)*tb(i)-sbz(i)*tbz(i))/tb(i)**2
     &          -(sazz(i)*ta(i)-saz(i)*taz(i))/ta(i)**2)
     &         +z(1)*(((sbzzz(i)*tb(i)-2*sbzz(i)*tbz(i)-sbz(i)*tbzz(i))*
     &         tb(i)+2*sbz(i)*tbz(i)**2)/tb(i)**3-((sazzz(i)*ta(i)-
     &         2*sazz(i)
     &         *taz(i)-saz(i)*tazz(i))*ta(i)+2*saz(i)*taz(i)**2)/ta(i)
     &         **3))
               fxxy=fxxy-2*tt*z(1)*(((sbxxy(i)*tb(i)-sbxx(i)*tby(i)-
     &         (sbxy(i)
     &         *tbx(i)+sbx(i)*tbxy(i)))*tb(i)+2*sbx(i)*tbx(i)*tby(i))
     &         /tb(i)**3-((saxxy(i)*ta(i)-saxx(i)*tay(i)-(saxy(i)*tax(i)
     &         +sax(i)*taxy(i)))*ta(i)+2*sax(i)*tax(i)*tay(i))/ta(i)**3)
               fxxz=fxxz-2*tt*((sbxx(i)*tb(i)-sbx(i)*tbx(i))/tb(i)**2-
     &         (saxx(i)*ta(i)-sax(i)*tax(i))/ta(i)**2+z(1)*(((sbxxz(i)
     &         *tb(i)
     &         -sbxx(i)*tbz(i)-sbxz(i)*tbx(i)-sbx(i)*tbxz(i))*tb(i)+
     &         2*sbx(i)
     &         *tbx(i)*tbz(i))/tb(i)**3-((saxxz(i)*ta(i)-saxx(i)*taz(i)-
     &         saxz(i)*tax(i)-sax(i)*taxz(i))*ta(i)+2*sax(i)*tax(i)
     &         *taz(i))
     &         /ta(i)**3))
               fyyz=fyyz-2*tt*((sbyy(i)*tb(i)-sby(i)*tby(i))/tb(i)**2-
     &         (sayy(i)*ta(i)-say(i)*tay(i))/ta(i)**2+z(1)*(((sbyyz(i)
     &         *tb(i)
     &         -sbyy(i)*tbz(i)-sbyz(i)*tby(i)-sby(i)*tbyz(i))*tb(i)+
     &         2*sby(i)
     &         *tby(i)*tbz(i))/tb(i)**3-((sayyz(i)*ta(i)-sayy(i)*taz(i)-
     &         sayz(i)*tay(i)-say(i)*tayz(i))*ta(i)+2*say(i)*tay(i)
     &         *taz(i))
     &         /ta(i)**3))
               fxyy=fxyy-2*tt*z(1)*(((sbxyy(i)*tb(i)-2*sbxy(i)*tby(i)-
     &         sbx(i)*
     &         tbyy(i))*tb(i)+2*sbx(i)*tby(i)**2)/tb(i)**3-((saxyy(i)
     &         *ta(i)
     &         -2*saxy(i)*tay(i)-sax(i)*tayy(i))*ta(i)+2*sax(i)*tay(i)
     &         **2)/
     &         ta(i)**3)
               fxyz=fxyz-2*tt*((sbxy(i)*tb(i)-sbx(i)*tby(i))/tb(i)**2-
     &         (saxy(i)
     &         *ta(i)-sax(i)*tay(i))/ta(i)**2+z(1)*(((sbxyz(i)*tb(i)-
     &         sbxy(i)
     &         *tbz(i)-sbxz(i)*tby(i)-sbx(i)*tbyz(i))*tb(i)+2*sbx(i)
     &         *tby(i)
     &         *tbz(i))/tb(i)**3-((saxyz(i)*ta(i)-saxy(i)*taz(i)-saxz(i)
     &         *
     &         tay(i)-sax(i)*tayz(i))*ta(i)+2*sax(i)*tay(i)*taz(i))
     &         /ta(i)**3))
               fxzz=fxzz-2*tt*(2*((sbxz(i)*tb(i)-sbx(i)*tbz(i))/tb(i)**2
     &          -(saxz(i)*ta(i)-sax(i)*taz(i))/ta(i)**2)+z(1)
     &         *(((sbxzz(i)*tb(i)
     &         -2*sbxz(i)*tbz(i)-sbx(i)*tbzz(i))*tb(i)+2*sbx(i)*tbz(i)
     &         **2)/
     &         tb(i)**3-((saxzz(i)*ta(i)-2*saxz(i)*taz(i)-sax(i)*tazz(i)
     &         )
     &         *ta(i)+2*sax(i)*taz(i)**2)/ta(i)**3))
               fyzz=fyzz-2*tt*(2*((sbyz(i)*tb(i)-sby(i)*tbz(i))/tb(i)**2
     &          -(sayz(i)*ta(i)-say(i)*taz(i))/ta(i)**2)+z(1)
     &         *(((sbyzz(i)*tb(i)
     &         -2*sbyz(i)*tbz(i)-sby(i)*tbzz(i))*tb(i)+2*sby(i)*tbz(i)
     &         **2)/
     &         tb(i)**3-((sayzz(i)*ta(i)-2*sayz(i)*taz(i)-say(i)*tazz(i)
     &         )
     &         *ta(i)+2*say(i)*taz(i)**2)/ta(i)**3))
            end if
         end if
   20 continue
      
      S(1,1) =     + ( 2.0D0           ) * FXZ - LOCALZ * FXXX
      S(1,2) =     + ( 2.0D0 * V       ) * FYZ - LOCALZ * FXXY
      S(1,3) = FZZ + (1.0D0 - 2.0D0 * V) * FYY - LOCALZ * FXXZ
      S(2,1) =     + ( 2.0D0 * V       ) * FXZ - LOCALZ * FXYY
      S(2,2) =     + ( 2.0D0           ) * FYZ - LOCALZ * FYYY
      S(2,3) = FZZ + (1.0D0 - 2.0D0 * V) * FXX - LOCALZ * FYYZ
      S(3,1) =     + (0.0D0            )       - LOCALZ * FXZZ
      S(3,2) =     + (0.0D0            )       - LOCALZ * FYZZ
      S(3,3) = FZZ + (0.0D0            )       - LOCALZ * FZZZ
      S(4,1) =     + (1.0D0 - V        ) * FYZ - LOCALZ * FXXY
      S(4,2) =     + (1.0D0 - V        ) * FXZ - LOCALZ * FXYY
      S(4,3) =     - (1.0D0 - 2.0D0 * V) * FXY - LOCALZ * FXYZ
      S(5,1) =     + ( - V             ) * FXY - LOCALZ * FXYZ
      S(5,2) = FZZ + ( + V             ) * FXX - LOCALZ * FYYZ
      S(5,3) =     + (0.0D0            )       - LOCALZ * FYZZ
      S(6,1) = FZZ + ( + V             ) * FYY - LOCALZ * FXXZ
      S(6,2) =     + ( - V             ) * FXY - LOCALZ * FXYZ
      S(6,3) =     + (0.0D0            )       - LOCALZ * FXZZ
      
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
       
      return
      end subroutine diag_stress_dd3d_tri_ana

      
      subroutine diag_dis_dd3d_tri_ana (xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,
     &zz3,E,v,IEvect,COEFF)

      real*8 x(5),y(5),z(5),a(5),b(5),pp(4,5),sb(5),sa(5),sbx(5),sax(5),
     &m(5),g(5),h(5),gx(5),hx(5),gy(5),hy(5),sby(5),say(5),gz(5),hz(5),
     &sbz(5),saz(5),gxx(5),hxx(5),tb(5),ta(5),tbx(5),tax(5),sbxx(5),
     &saxx(5),gxy(5),hxy(5),sbxy(5),saxy(5),tby(5),tay(5),sbyy(5),
     &sayy(5),gyy(5),hyy(5),gyz(5),hyz(5),sbyz(5),sayz(5),tbz(5),taz(5),
     &gxz(5),hxz(5),sbxz(5),saxz(5),gzz(5),hzz(5),sbzz(5),sazz(5),
     &gxxx(5),hxxx(5),sbxxx(5),saxxx(5),tbxx(5),taxx(5),gyyy(5),hyyy(5),
     &sbyyy(5),sayyy(5),tbyy(5),tayy(5),gzzz(5),hzzz(5),sbzzz(5),
     &sazzz(5),tbzz(5),tazz(5),gxxy(5),hxxy(5),sbxxy(5),saxxy(5),tbxy(5)
     &,taxy(5),gxxz(5),hxxz(5),sbxxz(5),saxxz(5),tbxz(5),taxz(5),gyyz(5)
     &,hyyz(5),sbyyz(5),sayyz(5),gxyy(5),hxyy(5),sbxyy(5),saxyy(5),
     &gxyz(5),hxyz(5),sbxyz(5),saxyz(5),tbyz(5),tayz(5),gxzz(5),hxzz(5),
     &sbxzz(5),saxzz(5),gyzz(5),hyzz(5),sbyzz(5),sayzz(5),rsb(5),rsa(5),
     &rtb(5),rta(5),rsbx(5),rsax(5),rsaz(5),rsbz(5),rtbx(5),rtax(5),
     &rtbz(5),rtaz(5),rtby(5),rtay(5),rsby(5),rsay(5) 
      real*8 lx,ly,f,fx,fy,fz,fxx,fxy,fyy,fyz,fxz,fzz,fxxx,fyyy,fzzz,
     &fxxy,fxxz,fyyz,fxyy,fxyz,fxzz,fyzz,EE,x1,x2,x3,y1,y2,y3,z1,z2,z3,
     &x4,y4,z4
      real*8::E,v,IEvect(3,3),JEvect(3,3),COEFF(3,3),JI(3,3),trac(3)
      real*8::u(3,3),xx(4),yy(4),zz(4),gpx,gpy,gpz,GG,LOCALZ
      real*8::xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,cond
      integer::i,j

	
       cond=1./8./3.141592653/(1-v)
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
      
      do 10 i=1,4
         do 15 j=1,5
            pp(i,j)=dsqrt((x(j)-x(i))**2+(y(j)-y(i))**2+(z(j)-z(i))**2)
   15    continue
   10 continue
      EE=DBLE(0.0D0)
      f =0.0D0
      fx =0.0D0
      fy =0.0d0
      fz =0.0D0
      fxx =0.0D0
      fxy =0.0D0
      fyy =0.0D0
      fyz =0.0D0
      fxz =0.0D0
      fzz =0.0D0
      fxxx=0.0D0
      fyyy=0.0D0
      fzzz=0.0D0
      fxxy=0.0D0
      fxxz=0.0D0
      fyyz=0.0D0
      fxyy=0.0D0
      fxyz=0.0D0
      fxzz=0.0D0
      fyzz=0.0D0
      do 20 i=2,4
         a(i)=((x(i)-x(1))*(x(i+1)-x(i))+(y(i)-y(1))*(y(i+1)-y(i)))/pp(
     &   i,i+1)
 
         b(i)=((x(i+1)-x(1))*(x(i+1)-x(i))+(y(i+1)-y(1))*(y(i+1)-y(i)))
     &   /pp(i,i+1)
 
         lx=x(i)-x(i+1) 
         ly=y(i)-y(i+1)   
         m(i)=(x(i)-x(1))*(y(i+1)-y(1))-(x(i+1)-x(1))*(y(i)-y(1))
         
         g(i)=dlog(pp(1,i+1)+b(i))
         h(i)=dlog(pp(1,i)+a(i))
         
         gx(i)=((x(1)-x(i+1))/pp(1,i+1)+lx/pp(i,i+1))/(pp(1,i+1)+b(i))
         hx(i)=((x(1)-x(i))/pp(1,i)+lx/pp(i,i+1))/(pp(1,i)+a(i))
 
         gy(i)=((y(1)-y(i+1))/pp(1,i+1)+ly/pp(i,i+1))/(pp(1,i+1)+b(i))
         hy(i)=((y(1)-y(i))/pp(1,i)+ly/pp(i,i+1))/(pp(1,i)+a(i))
 
         gz(i)=(z(1)/pp(1,i+1))/(pp(1,i+1)+b(i))
         hz(i)=(z(1)/pp(1,i))/(pp(1,i)+a(i))
 
         gxx(i)=-gx(i)*gx(i)+(((y(1)-y(i+1))**2+z(1)**2)/pp(1,i+1)**3
     &    )/(pp(1,i+1)+b(i))
         hxx(i)=-hx(i)*hx(i)+(((y(1)-y(i))**2+z(1)**2)/pp(1,i)**3
     &    )/(pp(1,i)+a(i))
 
         gxy(i)=-(gx(i)*((y(1)-y(i+1))/pp(1,i+1)+ly/pp(i,i+1))+
     &   (x(1)-x(i+1))*(y(1)-y(i+1))/pp(1,i+1)**3)/(pp(1,i+1)+b(i))
         hxy(i)=-(hx(i)*((y(1)-y(i))/pp(1,i)+ly/pp(i,i+1))+
     &   (x(1)-x(i))*(y(1)-y(i))/pp(1,i)**3)/(pp(1,i)+a(i))
 
         gyy(i)=-gy(i)*gy(i)+(((x(1)-x(i+1))**2+z(1)**2)/pp(1,i+1)**3
     &    )/(pp(1,i+1)+b(i))
         hyy(i)=-hy(i)*hy(i)+(((x(1)-x(i))**2+z(1)**2)/pp(1,i)**3)
     &   /(pp(1,i)+a(i))
         gyz(i)=-(gy(i)*z(1)/pp(1,i+1)+(y(1)-y(i+1))*z(1)/pp(1,i+1)**3)
     &   /(pp(1,i+1)+b(i))
         hyz(i)=-(hy(i)*z(1)/pp(1,i)+(y(1)-y(i))*z(1)/pp(1,i)**3)
     &   /(pp(1,i)+a(i))
         gxz(i)=-(gx(i)*z(1)/pp(1,i+1)+(x(1)-x(i+1))*z(1)/pp(1,i+1)**3)
     &   /(pp(1,i+1)+b(i))
         hxz(i)=-(hx(i)*z(1)/pp(1,i)+(x(1)-x(i))*z(1)/pp(1,i)**3)
     &   /(pp(1,i)+a(i))
         gzz(i)=-gz(i)*gz(i)+(((x(1)-x(i+1))**2+(y(1)-y(i+1))**2)
     &   /pp(1,i+1)**3)/(pp(1,i+1)+b(i))
         hzz(i)=-hz(i)*hz(i)+(((x(1)-x(i))**2+(y(1)-y(i))**2)
     &   /pp(1,i)**3)/(pp(1,i)+a(i))
         gxxx(i)=2.*gx(i)**3-3.*gx(i)/(pp(1,i+1)+b(i))*((y(1)-y(i+1))**2
     &    +z(1)**2)/pp(1,i+1)**3-3.*((x(1)-x(i+1))*((y(1)-y(i+1))**2
     &    +z(1)**2)/pp(1,i+1)**5)/(pp(1,i+1)+b(i))
         hxxx(i)=2.*hx(i)**3-3.*hx(i)/(pp(1,i)+a(i))*((y(1)-y(i))**2
     &    +z(1)**2)/pp(1,i)**3-3.*((x(1)-x(i))*((y(1)-y(i))**2+z(1)**2)
     &   /pp(1,i)**5)/(pp(1,i)+a(i))
         gyyy(i)=2.*gy(i)**3-3.*gy(i)/(pp(1,i+1)+b(i))*((x(1)-x(i+1))**2
     &    +z(1)**2)/pp(1,i+1)**3-3.*((y(1)-y(i+1))*((x(1)-x(i+1))**2
     &    +z(1)**2)/pp(1,i+1)**5)/(pp(1,i+1)+b(i))
         hyyy(i)=2.*hy(i)**3-3.*hy(i)/(pp(1,i)+a(i))*((x(1)-x(i))**2
     &    +z(1)**2)/pp(1,i)**3-3.*((y(1)-y(i))*((x(1)-x(i))**2
     &    +z(1)**2)/pp(1,i)**5)/(pp(1,i)+a(i))
         gzzz(i)=2.*gz(i)**3-3.*gz(i)/(pp(1,i+1)+b(i))*((x(1)-x(i+1))**2
     &    +(y(1)-y(i+1))**2)/pp(1,i+1)**3-3.*(z(1)*((x(1)-x(i+1))**2
     &    +(y(1)-y(i+1))**2)/pp(1,i+1)**5)/(pp(1,i+1)+b(i))
         hzzz(i)=2.*hz(i)**3-3.*hz(i)/(pp(1,i)+a(i))*((x(1)-x(i))**2
     &   +(y(1)-y(i))**2)/pp(1,i)**3-3.*(z(1)*((x(1)-x(i))**2+(y(1)-y(i)
     &   )**2)/pp(1,i)**5)/(pp(1,i)+a(i))
         gxxy(i)=2.*gx(i)**2*gy(i)+2.*gx(i)*(x(1)-x(i+1))*(y(1)-y(i+1))
     &   /pp(1,i+1)**3/(pp(1,i+1)+b(i))-gy(i)*((y(1)-y(i+1))**2+z(1)**2
     &    )/pp(1,i+1)**3/(pp(1,i+1)+b(i))+(y(1)-y(i+1))*(2.*pp(1,i+1)**2
     &    -3.*((y(1)-y(i+1))**2+z(1)**2))/pp(1,i+1)**5/(pp(1,i+1)+b(i))
         hxxy(i)=2.*hx(i)**2*hy(i)+2.*hx(i)*(x(1)-x(i))*(y(1)-y(i))
     &   /pp(1,i)**3/(pp(1,i)+a(i))-hy(i)*((y(1)-y(i))**2+z(1)**2)
     &   /pp(1,i)**3/(pp(1,i)+a(i))+(y(1)-y(i))*(2.*pp(1,i)**2-3.*(
     &   (y(1)-y(i))**2+z(1)**2))/pp(1,i)**5/(pp(1,i)+a(i))
         gxxz(i)=2.*gx(i)**2*gz(i)+2.*gx(i)*z(1)*(x(1)-x(i+1))/
     &   pp(1,i+1)**3/(pp(1,i+1)+b(i))-gz(i)*((y(1)-y(i+1))**2+z(1)**2)
     &   /pp(1,i+1)**3/(pp(1,i+1)+b(i))+z(1)*(2.*pp(1,i+1)**2-3.*
     &   ((y(1)-y(i+1))**2+z(1)**2))/pp(1,i+1)**5/(pp(1,i+1)+b(i))
         hxxz(i)=2.*hx(i)**2*hz(i)+2.*hx(i)*z(1)*(x(1)-x(i))/pp(1,i)**3
     &    /(pp(1,i)+a(i))-hz(i)*((y(1)-y(i))**2+z(1)**2)/pp(1,i)**3/
     &   (pp(1,i)+a(i))+z(1)*(2.*pp(1,i)**2-3.*((y(1)-y(i))**2+z(1)**2))
     &   /pp(1,i)**5/(pp(1,i)+a(i))
         gyyz(i)=2.*gy(i)**2*gz(i)+2.*gy(i)*z(1)*(y(1)-y(i+1))/
     &   pp(1,i+1)**3/(pp(1,i+1)+b(i))-gz(i)*((x(1)-x(i+1))**2+z(1)**2)
     &   /pp(1,i+1)**3/(pp(1,i+1)+b(i))+z(1)*(2*pp(1,i+1)**2-3.*
     &   ((x(1)-x(i+1))**2+z(1)**2))/pp(1,i+1)**5/(pp(1,i+1)+b(i))
         hyyz(i)=2.*hy(i)**2*hz(i)+2.*hy(i)*z(1)*(y(1)-y(i))/pp(1,i)**3
     &    /(pp(1,i)+a(i))-hz(i)*((x(1)-x(i))**2+z(1)**2)/pp(1,i)**3/
     &   (pp(1,i)+a(i))+z(1)*(2.*pp(1,i)**2-3.*((x(1)-x(i))**2+z(1)**2))
     &   /pp(1,i)**5/(pp(1,i)+a(i))
         gxyy(i)=2.*gx(i)*gy(i)**2+2.*gy(i)*(x(1)-x(i+1))*(y(1)-y(i+1))/
     &   pp(1,i+1)**3/(pp(1,i+1)+b(i))-gx(i)*((x(1)-x(i+1))**2+z(1)**2
     &    )/pp(1,i+1)**3/(pp(1,i+1)+b(i))+(x(1)-x(i+1))*(2.*pp(1,i+1)**2
     &    -3.*((x(1)-x(i+1))**2+z(1)**2))/pp(1,i+1)**5/(pp(1,i+1)+b(i))
         hxyy(i)=2.*hx(i)*hy(i)**2+2.*hy(i)*(x(1)-x(i))*(y(1)-y(i))/
     &   pp(1,i)**3/(pp(1,i)+a(i))-hx(i)*((x(1)-x(i))**2+z(1)**2)
     &   /pp(1,i)**3/(pp(1,i)+a(i))+(x(1)-x(i))*(2.*pp(1,i)**2
     &    -3.*((x(1)-x(i))**2+z(1)**2))/pp(1,i)**5/(pp(1,i)+a(i))
         gxyz(i)=2.*gx(i)*gy(i)*gz(i)+z(1)*(((x(1)-x(i+1))*(y(1)-y(i+1))
     &   /pp(1,i+1)+(y(1)-y(i+1))*((x(1)-x(i+1))/pp(1,i+1)+lx/
     &   pp(i,i+1))+(x(1)-x(i+1))*((y(1)-y(i+1))/pp(1,i+1)+ly/pp(i,i+1)
     &   ))/(pp(1,i+1)+b(i))**2/pp(1,i+1)**3+3.*(x(1)-x(i+1))*(y(1)-
     &   y(i+1))/pp(1,i+1)**5/(pp(1,i+1)+b(i)))
         hxyz(i)=2.*hx(i)*hy(i)*hz(i)+z(1)*(((x(1)-x(i))*(y(1)-y(i))
     &   /pp(1,i)+(y(1)-y(i))*((x(1)-x(i))/pp(1,i)+lx/pp(i,i+1))
     &   +(x(1)-x(i))*((y(1)-y(i))/pp(1,i)+ly/pp(i,i+1)))
     &   /(pp(1,i)+a(i))**2/pp(1,i)**3+3.*(x(1)-x(i))*(y(1)-y(i))
     &   /pp(1,i)**5/(pp(1,i)+a(i)))
         gxzz(i)=2.*gx(i)*gz(i)**2+2.*gz(i)*(x(1)-x(i+1))*z(1)/pp(1,i+1)
     &   **3
     &    /(pp(1,i+1)+b(i))-gx(i)*((x(1)-x(i+1))**2+(y(1)-y(i+1))**2)/
     &   pp(1,i+1)**3/(pp(1,i+1)+b(i))+(x(1)-x(i+1))*(2.*pp(1,i+1)**2
     &    -3.*((x(1)-x(i+1))**2+(y(1)-y(i+1))**2))/pp(1,i+1)**5/
     &   (pp(1,i+1)+b(i)) 
         hxzz(i)=2.*hx(i)*hz(i)**2+2.*hz(i)*(x(1)-x(i))*z(1)/pp(1,i)**3
     &    /(pp(1,i)+a(i))-hx(i)*((x(1)-x(i))**2+(y(1)-y(i))**2)/
     &   pp(1,i)**3/(pp(1,i)+a(i))+(x(1)-x(i))*(2.*pp(1,i)**2
     &    -3.*((x(1)-x(i))**2+(y(1)-y(i))**2))/pp(1,i)**5/(pp(1,i)+a(i))
 
         gyzz(i)=2.*gy(i)*gz(i)**2+2.*gz(i)*(y(1)-y(i+1))*z(1)/pp(1,i+1)
     &   **3
     &    /(pp(1,i+1)+b(i))-gy(i)*((x(1)-x(i+1))**2+(y(1)-y(i+1))**2)/
     &   pp(1,i+1)**3/(pp(1,i+1)+b(i))+(y(1)-y(i+1))*(2.*pp(1,i+1)**2
     &    -3.*((x(1)-x(i+1))**2+(y(1)-y(i+1))**2))/pp(1,i+1)**5/
     &   (pp(1,i+1)+b(i))
         hyzz(i)=2.*hy(i)*hz(i)**2+2.*hz(i)*(y(1)-y(i))*z(1)/pp(1,i)**3
     &    /(pp(1,i)+a(i))-hy(i)*((x(1)-x(i))**2+(y(1)-y(i))**2)/
     &   pp(1,i)**3/(pp(1,i)+a(i))+(y(1)-y(i))*(2.*pp(1,i)**2
     &    -3.*((x(1)-x(i))**2+(y(1)-y(i))**2))/pp(1,i)**5/(pp(1,i)+a(i))
 
         if (((z(1).ge.-0.1D-10).and.(z(1).le.0.1D-10)).and.
     &   ((lx.le.-0.1D-10).or.(lx.ge.0.1D-10))) then
            tt=1
            rsb(i)=(pp(i,i+1)/(-lx))*(y(i+1)-y(1)+b(i)+(1-ly
     &       /pp(i,i+1))*pp(1,i+1))
            rsa(i)=(pp(i,i+1)/(-lx))*(y(i)-y(1)+a(i)+(1-ly/pp(i,i+1))
     &      *pp(1,i))
 
            rsbx(i)=pp(i,i+1)/(-lx)*(lx/pp(i,i+1)+(1-ly/pp(i,i+1))
     &      *(-x(i+1)+x(1))/pp(1,i+1))
            rsax(i)=pp(i,i+1)/(-lx)*(lx/pp(i,i+1)+(1-ly/pp(i,i+1))
     &      *(-x(i)+x(1))/pp(1,i))
            rsbz(i)=pp(i,i+1)/(-lx)*(1-ly/pp(i,i+1))*z(1)/pp(1,i+1)
            rsaz(i)=pp(i,i+1)/(-lx)*(1-ly/pp(i,i+1))*z(1)/pp(1,i)
            rsby(i)=pp(i,i+1)/(-lx)*(-1+ly/pp(i,i+1)+(1-ly/pp(i,i+1))
     &      *(-y(i+1)+y(1))/pp(1,i+1))
            rsay(i)=pp(i,i+1)/(-lx)*(-1+ly/pp(i,i+1)+(1-ly/pp(i,i+1))
     &      *(-y(i)+y(1))/pp(1,i))
 
 
            rtb(i)=rsb(i)**2+z(1)**2
            rta(i)=rsa(i)**2+z(1)**2
            rtbx(i)=2*rsb(i)*rsbx(i)
            rtax(i)=2*rsa(i)*rsax(i)
            rtbz(i)=2*(rsb(i)*rsbz(i)+z(1))
            rtaz(i)=2*(rsa(i)*rsaz(i)+z(1))
            rtby(i)=2*rsb(i)*rsby(i)
            rtay(i)=2*rsa(i)*rsay(i)
 
            fzz=fzz+(m(i)*(gzz(i)-hzz(i)))/pp(i,i+1)+4*tt*(rsb(i)/rtb(i)
     &      -rsa(i)/rta(i))
            fxzz=fxzz+(ly*(gzz(i)-hzz(i))+m(i)*(gxzz(i)-hxzz(i)))
     &      /pp(i,i+1)-2*tt*(2*((rsb(i)*rtbx(i)-rsbx(i)*rtb(i))/rtb(i)
     &      **2
     &       -(rsa(i)*rtax(i)-rsax(i)*rta(i))/rta(i)**2)-((rsbz(i)-
     &      rsb(i))
     &      *rtbz(i)/rtb(i)**2-(rsaz(i)-rsa(i))*rtaz(i)/rta(i)**2))
            fyzz=fyzz+(-lx*(gzz(i)-hzz(i))+m(i)*(gyzz(i)-hyzz(i)))
     &      /pp(i,i+1)-2*tt*(2*((rsb(i)*rtby(i)-rsby(i)*rtb(i))/rtb(i)
     &      **2
     &       -(rsa(i)*rtay(i)-rsay(i)*rta(i))/rta(i)**2)-((rsbz(i)-
     &      rsb(i))
     &      *rtbz(i)/rtb(i)**2-(rsaz(i)-rsa(i))*rtaz(i)/rta(i)**2))
         else
            fzz=fzz+(m(i)*(gzz(i)-hzz(i)))/pp(i,i+1)
            fxzz=fxzz+(ly*(gzz(i)-hzz(i))+m(i)*(gxzz(i)-hxzz(i)))/pp(i,
     &      i+1)
            fyzz=fyzz+(-lx*(gzz(i)-hzz(i))+m(i)*(gyzz(i)-hyzz(i)))
     &      /pp(i,i+1)
         end if
 
         f=f+m(i)/pp(i,i+1)*(g(i)-h(i))
         fx=fx+(ly*(g(i)-h(i))+m(i)*(gx(i)-hx(i)))/pp(i,i+1)
         
         fy=fy+(-lx*(g(i)-h(i))+m(i)*(gy(i)-hy(i)))/pp(i,i+1)
         fz=fz+(m(i)*(gz(i)-hz(i)))/pp(i,i+1)
         fxx=fxx+(2*ly*(gx(i)-hx(i))+m(i)*(gxx(i)-hxx(i)))/pp(i,i+1)
         fxy=fxy+(-lx*(gx(i)-hx(i))+ly*(gy(i)-hy(i))+m(i)*(gxy(i)-
     &   hxy(i)))/pp(i,i+1) 
         fyy=fyy+(-2*lx*(gy(i)-hy(i))+m(i)*(gyy(i)-hyy(i)))/pp(i,i+1)
         fyz=fyz+(-lx*(gz(i)-hz(i))+m(i)*(gyz(i)-hyz(i)))/pp(i,i+1)
         fxz=fxz+(ly*(gz(i)-hz(i))+m(i)*(gxz(i)-hxz(i)))/pp(i,i+1)
         fxxx=fxxx+(3*ly*(gxx(i)-hxx(i))+m(i)*(gxxx(i)-hxxx(i)))
     &   /pp(i,i+1)
         fyyy=fyyy+(-3*lx*(gyy(i)-hyy(i))+m(i)*(gyyy(i)-hyyy(i)))
     &   /pp(i,i+1)
         fzzz=fzzz+m(i)*(gzzz(i)-hzzz(i))/pp(i,i+1)
         fxxy=fxxy+(2*ly*(gxy(i)-hxy(i))-lx*(gxx(i)-hxx(i))+m(i)*
     &   (gxxy(i)-hxxy(i)))/pp(i,i+1)
         fxxz=fxxz+(2*ly*(gxz(i)-hxz(i))+m(i)*(gxxz(i)-hxxz(i)))
     &   /pp(i,i+1)
         fyyz=fyyz+(-2*lx*(gyz(i)-hyz(i))+m(i)*(gyyz(i)-hyyz(i)))
     &   /pp(i,i+1)
         fxyy=fxyy+(-2*lx*(gxy(i)-hxy(i))+ly*(gyy(i)-hyy(i))+m(i)*
     &   (gxyy(i)-hxyy(i)))/pp(i,i+1)
         fxyz=fxyz+(-lx*(gxz(i)-hxz(i))+ly*(gyz(i)-hyz(i))+m(i)*
     &   (gxyz(i)-hxyz(i)))/pp(i,i+1)
 
         if ((z(1).le.-0.1D-10).or.(z(1).ge.0.1D-10)) then
            if ((lx.le.-0.1D-10).or.(lx.ge.0.1D-10)) then
               tt=1
               sb(i)=(pp(i,i+1)/((-lx)*z(1)))*(y(i+1)-y(1)+b(i)+(1-ly
     &          /pp(i,i+1))*pp(1,i+1))
               sa(i)=(pp(i,i+1)/((-lx)*z(1)))*(y(i)-y(1)+a(i)+(1-
     &         ly/pp(i,i+1))
     &         *pp(1,i))
               sbx(i)=pp(i,i+1)/((-lx)*z(1))*(lx/pp(i,i+1)+(1-ly/pp(i,i+
     &         1))
     &         *(-x(i+1)+x(1))/pp(1,i+1))
               sax(i)=pp(i,i+1)/((-lx)*z(1))*(lx/pp(i,i+1)+(1-ly/pp(i,i+
     &         1))
     &         *(-x(i)+x(1))/pp(1,i))
 
               sby(i)=pp(i,i+1)/((-lx)*z(1))*(-1+ly/pp(i,i+1)+(1-
     &         ly/pp(i,i+1))
     &         *(-y(i+1)+y(1))/pp(1,i+1))
               say(i)=pp(i,i+1)/((-lx)*z(1))*(-1+ly/pp(i,i+1)+(1-
     &         ly/pp(i,i+1))
     &         *(-y(i)+y(1))/pp(1,i))
 
               sbz(i)=sb(i)*(-1)/z(1)+pp(i,i+1)/((-lx)*z(1))*((1-
     &         ly/pp(i,i+1))
     &         *z(1)/pp(1,i+1))
               saz(i)=sa(i)*(-1)/z(1)+pp(i,i+1)/((-lx)*z(1))*((1-
     &         ly/pp(i,i+1))
     &         *z(1)/pp(1,i))
               sbxx(i)=pp(i,i+1)/((-lx)*z(1))*((1-ly/pp(i,i+1))*(pp(1,i+
     &         1)
     &         **2-(x(1)-x(i+1))**2)/pp(1,i+1)**3)
               saxx(i)=pp(i,i+1)/((-lx)*z(1))*((1-ly/pp(i,i+1))*(pp(1,i)
     &         **2-(x(1)-x(i))**2)/pp(1,i)**3)
               sbxy(i)=pp(i,i+1)/((-lx)*z(1))*((1-ly/pp(i,i+1))*(x(i+1)-
     &         x(1))
     &         *(y(1)-y(i+1))/pp(1,i+1)**3)
               saxy(i)=pp(i,i+1)/((-lx)*z(1))*((1-ly/pp(i,i+1))*(x(i)-
     &         x(1))
     &         *(y(1)-y(i))/pp(1,i)**3)
               sbyy(i)=pp(i,i+1)/((-lx)*z(1))*((1-ly/pp(i,i+1))*(pp(1,i+
     &         1)
     &         **2-(y(1)-y(i+1))**2)/pp(1,i+1)**3)
               sayy(i)=pp(i,i+1)/((-lx)*z(1))*((1-ly/pp(i,i+1))*(pp(1,i)
     &         **2-(y(1)-Y(i))**2)/pp(1,i)**3)
               sbyz(i)=pp(i,i+1)/(-lx)*((-1)/z(1)**2*(-1+ly/pp(i,i+1)
     &         +(1-ly/pp(i,i+1))*(y(1)-y(i+1))/pp(1,i+1))-(1-ly/pp(i,i+
     &         1))
     &         *(y(1)-y(i+1))/pp(1,i+1)**3)
               sayz(i)=pp(i,i+1)/(-lx)*((-1)/z(1)**2*(-1+ly/pp(i,i+1)
     &         +(1-ly/pp(i,i+1))*(y(1)-y(i))/pp(1,i))-(1-ly/pp(i,i+1))
     &         *(y(1)-y(i))/pp(1,i)**3)
               sbxz(i)=pp(i,i+1)/(-lx)*((-1)/z(1)**2*(lx/pp(i,i+1)+(1-
     &         ly/
     &         pp(i,i+1))*(x(1)-x(i+1))/pp(1,i+1))-(1-ly/pp(i,i+1))*
     &         (x(1)-x(i+1))/pp(1,i+1)**3)
               saxz(i)=pp(i,i+1)/(-lx)*((-1)/z(1)**2*(lx/pp(i,i+1)+(1-
     &         ly/
     &         pp(i,i+1))*(x(1)-x(i))/pp(1,i))-(1-ly/pp(i,i+1))*
     &         (x(1)-x(i))/pp(1,i)**3)
 
               sbzz(i)=pp(i,i+1)/(-lx)*(2/z(1)**3*(y(i+1)-y(1)+b(i)+(1-
     &         ly/
     &         pp(i,i+1))*pp(1,i+1))-1/z(1)*((1-ly/pp(i,i+1))/pp(1,i+1))
     &         -z(1)*((1-ly/pp(i,i+1))/pp(1,i+1)**3))
               sazz(i)=pp(i,i+1)/(-lx)*(2/z(1)**3*(y(i)-y(1)+a(i)+(1-ly/
     &         pp(i,i+1))*pp(1,i))-1/z(1)*((1-ly/pp(i,i+1))/pp(1,i))
     &         -z(1)*((1-ly/pp(i,i+1))/pp(1,i)**3))
               sbxxx(i)=pp(i,i+1)/(-lx)/z(1)*((1-ly/pp(i,i+1))*(-3*
     &         (x(1)-X(i+1))*(pp(1,i+1)**2-(x(1)-x(i+1))**2)/pp(1,i+1)
     &         **5))
               saxxx(i)=pp(i,i+1)/(-lx)/z(1)*((1-ly/pp(i,i+1))*(-3*
     &         (x(1)-x(i))*(pp(1,i)**2-(x(1)-x(i))**2)/pp(1,i)**5))
               sbyyy(i)=pp(i,i+1)/(-lx)/z(1)*((1-ly/pp(i,i+1))*(-3*
     &         (y(1)-y(i+1))*(pp(1,i+1)**2-(y(1)-y(i+1))**2)/pp(1,i+1)
     &         **5))
               sayyy(i)=pp(i,i+1)/(-lx)/z(1)*((1-ly/pp(i,i+1))*(-3*
     &         (y(1)-y(i))*(pp(1,i)**2-(y(1)-y(i))**2)/pp(1,i)**5))
               sbzzz(i)=pp(i,i+1)/(-lx)*((-6)/z(1)**4*(-y(1)+y(i+1)+b(i)
     &         +(1-ly/pp(i,i+1))*pp(1,i+1))+3/z(1)**2*((1-ly/pp(i,i+1))
     &         /pp(1,i+1))+3*z(1)**2*((1-ly/pp(i,i+1))/pp(1,i+1)**5))
               sazzz(i)=pp(i,i+1)/(-lx)*((-6)/z(1)**4*(-y(1)+y(i)+a(i)
     &         +(1-ly/pp(i,i+1))*pp(1,i))+3/z(1)**2*((1-ly/pp(i,i+1))
     &         /pp(1,i)
     &         )+3*z(1)**2*((1-ly/pp(i,i+1))/pp(1,i)**5))
               sbxxy(i)=pp(i,i+1)/(-lx)/z(1)*((1-ly/pp(i,i+1))*((y(1)-
     &         y(i+1)
     &         )*(-pp(1,i+1)**2+3*(x(1)-x(i+1))**2)/pp(1,i+1)**5))
               saxxy(i)=pp(i,i+1)/(-lx)/z(1)*((1-ly/pp(i,i+1))*((y(1)-
     &         y(i))
     &         *(-pp(1,i)**2+3*(x(1)-x(i))**2)/pp(1,i)**5))
               sbxxz(i)=pp(i,i+1)/(-lx)*((-1)/z(1)**2*(1-ly/pp(i,i+1))*
     &         (pp(1,i+1)**2-(x(1)-x(i+1))**2)/pp(1,i+1)**3+(1-ly/pp(i,
     &         i+1))
     &         *(-pp(1,i+1)**2+3*(x(1)-x(i+1))**2)/pp(1,i+1)**5)
               saxxz(i)=pp(i,i+1)/(-lx)*((-1)/z(1)**2*(1-ly/pp(i,i+1))*
     &         (pp(1,i)**2-(x(1)-x(i))**2)/pp(1,i)**3+(1-ly/pp(i,i+1))*
     &         (-pp(1,i)**2+3*(x(1)-x(i))**2)/pp(1,i)**5)
               sbyyz(i)=pp(i,i+1)/(-lx)*((-1)/z(1)**2*(1-ly/pp(i,i+1))*
     &         (pp(1,i+1)**2-(y(1)-y(i+1))**2)/pp(1,i+1)**3+(1-ly/pp(i,
     &         i+1))
     &         *(-pp(1,i+1)**2+3*(y(1)-y(i+1))**2)/pp(1,i+1)**5)
               sayyz(i)=pp(i,i+1)/(-lx)*((-1)/z(1)**2*(1-ly/pp(i,i+1))*
     &         (pp(1,i)**2-(y(1)-y(i))**2)/pp(1,i)**3+(1-ly/pp(i,i+1))*
     &         (-pp(1,i)**2+3*(y(1)-y(i))**2)/pp(1,i)**5)
 
               sbxyy(i)=pp(i,i+1)/(-lx)/z(1)*(1-ly/pp(i,i+1))*(-(x(1)-
     &         x(i+1))
     &         *(pp(1,i+1)**2-3*(y(1)-y(i+1))**2))/pp(1,i+1)**5
               saxyy(i)=pp(i,i+1)/(-lx)/z(1)*(1-ly/pp(i,i+1))*(-(x(1)-
     &         x(i))
     &         *(pp(1,i)**2-3*(y(1)-y(i))**2))/pp(1,i)**5
               sbxyz(i)=pp(i,i+1)/(-lx)*((-1)/z(1)**2*(1-ly/pp(i,i+1))*
     &         (-(x(1)-x(i+1))*(y(1)-y(i+1))/pp(1,i+1)**3)+(1-ly/pp(i,i+
     &         1))
     &         *3*(x(1)-x(i+1))*(y(1)-y(i+1))/pp(1,i+1)**5)
               saxyz(i)=pp(i,i+1)/(-lx)*((-1)/z(1)**2*(1-ly/pp(i,i+1))*
     &         (-(x(1)-x(i))*(y(1)-y(i))/pp(1,i)**3)+(1-ly/pp(i,i+1))
     &         *3*(x(1)-x(i))*(y(1)-y(i))/pp(1,i)**5)
 
               sbxzz(i)=pp(i,i+1)/(-lx)*(2/z(1)**3*(lx/pp(i,i+1)+(1-ly/
     &         pp(i,i+1))*(x(1)-x(i+1))/pp(1,i+1))+1/z(1)*(1-ly/pp(i,i+
     &         1))
     &         *(x(1)-x(i+1))/pp(1,i+1)**3+z(1)*(1-ly/pp(i,i+1))*3*
     &         (x(1)-x(i+1))/pp(1,i+1)**5)
               saxzz(i)=pp(i,i+1)/(-lx)*(2/z(1)**3*(lx/pp(i,i+1)+(1-ly/
     &         pp(i,i+1))*(x(1)-x(i))/pp(1,i))+1/z(1)*(1-ly/pp(i,i+1))
     &         *(x(1)-x(i))/pp(1,i)**3+z(1)*(1-ly/pp(i,i+1))*3*
     &         (x(1)-x(i))/pp(1,i)**5)
 
               sbyzz(i)=pp(i,i+1)/(-lx)*(2/z(1)**3*(-1+ly/pp(i,i+1)+(1-
     &         ly/
     &         pp(i,i+1))*(y(1)-y(i+1))/pp(1,i+1))+1/z(1)*(1-ly/pp(i,i+
     &         1))
     &         *(y(1)-y(i+1))/pp(1,i+1)**3+z(1)*(1-ly/pp(i,i+1))*3*
     &         (y(1)-y(i+1))/pp(1,i+1)**5)
               sayzz(i)=pp(i,i+1)/(-lx)*(2/z(1)**3*(-1+ly/pp(i,i+1)+(1-
     &         ly/
     &         pp(i,i+1))*(y(1)-y(i))/pp(1,i))+1/z(1)*(1-ly/pp(i,i+1))
     &         *(y(1)-y(i))/pp(1,i)**3+z(1)*(1-ly/pp(i,i+1))*3*
     &         (y(1)-y(i))/pp(1,i)**5)
               tb(i)=sb(i)**2+1
               ta(i)=sa(i)**2+1
               tbx(i)=2*sb(i)*sbx(i)
               tax(i)=2*sa(i)*sax(i)
               tby(i)=2*sb(i)*sby(i)
               tay(i)=2*sa(i)*say(i)
               tbz(i)=2*sb(i)*sbz(i)
               taz(i)=2*sa(i)*saz(i)
               tbxx(i)=2*(sbx(i)**2+sb(i)*sbxx(i))
               taxx(i)=2*(sax(i)**2+sa(i)*saxx(i))
               tbyy(i)=2*(sby(i)**2+sb(i)*sbyy(i))
               tayy(i)=2*(say(i)**2+sa(i)*sayy(i))
               tbzz(i)=2*(sbz(i)**2+sb(i)*sbzz(i))
               tazz(i)=2*(saz(i)**2+sa(i)*sazz(i))
               tbxy(i)=2*(sby(i)*sbx(i)+sb(i)*sbxy(i))
               taxy(i)=2*(say(i)*sax(i)+sa(i)*saxy(i))
               tbxz(i)=2*(sbz(i)*sbx(i)+sb(i)*sbxz(i))
               taxz(i)=2*(saz(i)*sax(i)+sa(i)*saxz(i))
               tbyz(i)=2*(sbz(i)*sby(i)+sb(i)*sbyz(i))
               tayz(i)=2*(saz(i)*say(i)+sa(i)*sayz(i))
               f=f-2*tt*z(1)*(datan(sb(i))-datan(sa(i)))
               fx=fx-2*tt*z(1)*(sbx(i)/tb(i)-sax(i)/ta(i))
               fy=fy-2*tt*z(1)*(sby(i)/tb(i)-say(i)/ta(i))
               fz=fz-2*tt*(datan(sb(i))-datan(sa(i))+z(1)*(sbz(i)/tb(i)
     &         -saz(i)/ta(i)))
               fxx=fxx-2*tt*z(1)*((sbxx(i)*tb(i)-sbx(i)*tbx(i))/tb(i)**2
     &          -(saxx(i)*ta(i)-sax(i)*tax(i))/ta(i)**2)
               fxy=fxy-2*tt*z(1)*((sbxy(i)*tb(i)-sbx(i)*tby(i))/tb(i)**2
     &          -(saxy(i)*ta(i)-sax(i)*tay(i))/ta(i)**2)
               fyy=fyy-2*tt*z(1)*((sbyy(i)*tb(i)-sby(i)*tby(i))/tb(i)**2
     &          -(sayy(i)*ta(i)-say(i)*tay(i))/ta(i)**2)
               fyz=fyz-2*tt*(sby(i)/tb(i)-say(i)/ta(i)+z(1)*((sbyz(i)
     &         *tb(i)
     &         -sby(i)*tbz(i))/tb(i)**2-(sayz(i)*ta(i)-say(i)*taz(i))
     &         /ta(i)
     &         **2))
               fxz=fxz-2*tt*(sbx(i)/tb(i)-sax(i)/ta(i)+z(1)*((sbxz(i)
     &         *tb(i)
     &         -sbx(i)*tbz(i))/tb(i)**2-(saxz(i)*ta(i)-sax(i)*taz(i))
     &         /ta(i)
     &         **2))
               fzz=fzz-2*tt*(2*(sbz(i)/tb(i)-saz(i)/ta(i))+z(1)
     &         *((sbzz(i)*
     &         tb(i)-sbz(i)*tbz(i))/tb(i)**2-(sazz(i)*ta(i)-saz(i)
     &         *taz(i))
     &         /ta(i)**2))
               fxxx=fxxx-2*tt*z(1)*(((sbxxx(i)*tb(i)-2*sbxx(i)*tbx(i)-
     &         sbx(i)
     &         *tbxx(i))*tb(i)+2*sbx(i)*tbx(i)**2)/tb(i)**3-((saxxx(i)
     &         *ta(i)
     &         -2*saxx(i)*tax(i)-sax(i)*taxx(i))*ta(i)+2*sax(i)*tax(i)
     &         **2)
     &         /ta(i)**3)
               fyyy=fyyy-2*tt*z(1)*(((sbyyy(i)*tb(i)-2*sbyy(i)*tby(i)-
     &         sby(i)
     &         *tbyy(i))*tb(i)+2*sby(i)*tby(i)**2)/tb(i)**3-((sayyy(i)
     &         *ta(i)
     &         -2*sayy(i)*tay(i)-say(i)*tayy(i))*ta(i)+2*say(i)*tay(i)
     &         **2)
     &         /ta(i)**3)
               fzzz=fzzz-2*tt*(3*((sbzz(i)*tb(i)-sbz(i)*tbz(i))/tb(i)**2
     &          -(sazz(i)*ta(i)-saz(i)*taz(i))/ta(i)**2)
     &         +z(1)*(((sbzzz(i)*tb(i)-2*sbzz(i)*tbz(i)-sbz(i)*tbzz(i))*
     &         tb(i)+2*sbz(i)*tbz(i)**2)/tb(i)**3-((sazzz(i)*ta(i)-
     &         2*sazz(i)
     &         *taz(i)-saz(i)*tazz(i))*ta(i)+2*saz(i)*taz(i)**2)/ta(i)
     &         **3))
               fxxy=fxxy-2*tt*z(1)*(((sbxxy(i)*tb(i)-sbxx(i)*tby(i)-
     &         (sbxy(i)
     &         *tbx(i)+sbx(i)*tbxy(i)))*tb(i)+2*sbx(i)*tbx(i)*tby(i))
     &         /tb(i)**3-((saxxy(i)*ta(i)-saxx(i)*tay(i)-(saxy(i)*tax(i)
     &         +sax(i)*taxy(i)))*ta(i)+2*sax(i)*tax(i)*tay(i))/ta(i)**3)
               fxxz=fxxz-2*tt*((sbxx(i)*tb(i)-sbx(i)*tbx(i))/tb(i)**2-
     &         (saxx(i)*ta(i)-sax(i)*tax(i))/ta(i)**2+z(1)*(((sbxxz(i)
     &         *tb(i)
     &         -sbxx(i)*tbz(i)-sbxz(i)*tbx(i)-sbx(i)*tbxz(i))*tb(i)+
     &         2*sbx(i)
     &         *tbx(i)*tbz(i))/tb(i)**3-((saxxz(i)*ta(i)-saxx(i)*taz(i)-
     &         saxz(i)*tax(i)-sax(i)*taxz(i))*ta(i)+2*sax(i)*tax(i)
     &         *taz(i))
     &         /ta(i)**3))
               fyyz=fyyz-2*tt*((sbyy(i)*tb(i)-sby(i)*tby(i))/tb(i)**2-
     &         (sayy(i)*ta(i)-say(i)*tay(i))/ta(i)**2+z(1)*(((sbyyz(i)
     &         *tb(i)
     &         -sbyy(i)*tbz(i)-sbyz(i)*tby(i)-sby(i)*tbyz(i))*tb(i)+
     &         2*sby(i)
     &         *tby(i)*tbz(i))/tb(i)**3-((sayyz(i)*ta(i)-sayy(i)*taz(i)-
     &         sayz(i)*tay(i)-say(i)*tayz(i))*ta(i)+2*say(i)*tay(i)
     &         *taz(i))
     &         /ta(i)**3))
               fxyy=fxyy-2*tt*z(1)*(((sbxyy(i)*tb(i)-2*sbxy(i)*tby(i)-
     &         sbx(i)*
     &         tbyy(i))*tb(i)+2*sbx(i)*tby(i)**2)/tb(i)**3-((saxyy(i)
     &         *ta(i)
     &         -2*saxy(i)*tay(i)-sax(i)*tayy(i))*ta(i)+2*sax(i)*tay(i)
     &         **2)/
     &         ta(i)**3)
               fxyz=fxyz-2*tt*((sbxy(i)*tb(i)-sbx(i)*tby(i))/tb(i)**2-
     &         (saxy(i)
     &         *ta(i)-sax(i)*tay(i))/ta(i)**2+z(1)*(((sbxyz(i)*tb(i)-
     &         sbxy(i)
     &         *tbz(i)-sbxz(i)*tby(i)-sbx(i)*tbyz(i))*tb(i)+2*sbx(i)
     &         *tby(i)
     &         *tbz(i))/tb(i)**3-((saxyz(i)*ta(i)-saxy(i)*taz(i)-saxz(i)
     &         *
     &         tay(i)-sax(i)*tayz(i))*ta(i)+2*sax(i)*tay(i)*taz(i))
     &         /ta(i)**3))
               fxzz=fxzz-2*tt*(2*((sbxz(i)*tb(i)-sbx(i)*tbz(i))/tb(i)**2
     &          -(saxz(i)*ta(i)-sax(i)*taz(i))/ta(i)**2)+z(1)
     &         *(((sbxzz(i)*tb(i)
     &         -2*sbxz(i)*tbz(i)-sbx(i)*tbzz(i))*tb(i)+2*sbx(i)*tbz(i)
     &         **2)/
     &         tb(i)**3-((saxzz(i)*ta(i)-2*saxz(i)*taz(i)-sax(i)*tazz(i)
     &         )
     &         *ta(i)+2*sax(i)*taz(i)**2)/ta(i)**3))
               fyzz=fyzz-2*tt*(2*((sbyz(i)*tb(i)-sby(i)*tbz(i))/tb(i)**2
     &          -(sayz(i)*ta(i)-say(i)*taz(i))/ta(i)**2)+z(1)
     &         *(((sbyzz(i)*tb(i)
     &         -2*sbyz(i)*tbz(i)-sby(i)*tbzz(i))*tb(i)+2*sby(i)*tbz(i)
     &         **2)/
     &         tb(i)**3-((sayzz(i)*ta(i)-2*sayz(i)*taz(i)-say(i)*tazz(i)
     &         )
     &         *ta(i)+2*say(i)*taz(i)**2)/ta(i)**3))
            end if
         end if
   20 continue
      
      fz=2*3.141592653
      
      
      U(1,1)=  2.0 * ( 1.0 - V ) * FZ - localZ*FXX
      U(1,2)=                         - localZ*FXY
      U(1,3)= -( 1.0 - 2.0 * V ) * FX - localZ*FXZ
      U(2,1)=                         - localZ*FXY
      U(2,2)=  2.0 * ( 1.0 - V ) * FZ - localZ*FYY
      U(2,3)= -( 1.0 - 2.0 * V ) * FY - localZ*FYZ
      U(3,1)=  ( 1.0 - 2.0 * V ) * FX - localZ*FXZ
      U(3,2)=  ( 1.0 - 2.0 * V ) * FY - localZ*FYZ
      U(3,3)=  2.0 * ( 1.0 - V ) * FZ - localZ*FZZ
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
      
     
      COEFF(1,1)=cond*(U(1,1)*JI(1,1)+U(2,1)*JI(2,1)+U(3,1)*JI(3,1))
      COEFF(2,1)=cond*(U(1,1)*JI(1,2)+U(2,1)*JI(2,2)+U(3,1)*JI(3,2))
      COEFF(3,1)=cond*(U(1,1)*JI(1,3)+U(2,1)*JI(2,3)+U(3,1)*JI(3,3))
     
      COEFF(1,2)=cond*(U(1,2)*JI(1,1)+U(2,2)*JI(2,1)+U(3,2)*JI(3,1))
      COEFF(2,2)=cond*(U(1,2)*JI(1,2)+U(2,2)*JI(2,2)+U(3,2)*JI(3,2))
      COEFF(3,2)=cond*(U(1,2)*JI(1,3)+U(2,2)*JI(2,3)+U(3,2)*JI(3,3))
     
      COEFF(1,3)=cond*(U(1,3)*JI(1,1)+U(2,3)*JI(2,1)+U(3,3)*JI(3,1))
      COEFF(2,3)=cond*(U(1,3)*JI(1,2)+U(2,3)*JI(2,2)+U(3,3)*JI(3,2))
      COEFF(3,3)=cond*(U(1,3)*JI(1,3)+U(2,3)*JI(2,3)+U(3,3)*JI(3,3))
     
      !write(*,*)'coeff(1,1)',coeff(1,1)
      !write(*,*)'coeff(2,1)',coeff(2,1)
      !write(*,*)'coeff(3,1)',coeff(3,1)
      !write(*,*)'coeff(1,2)',coeff(1,2)
      !write(*,*)'coeff(2,2)',coeff(2,2)
      !write(*,*)'coeff(3,2)',coeff(3,2)
      !write(*,*)'coeff(1,3)',coeff(1,3)
      !write(*,*)'coeff(2,3)',coeff(2,3)
      !write(*,*)'coeff(3,3)',coeff(3,3)      
   
      return
      end subroutine diag_dis_dd3d_tri_ana
      
     
      end module subs_dd3d_tri_ana