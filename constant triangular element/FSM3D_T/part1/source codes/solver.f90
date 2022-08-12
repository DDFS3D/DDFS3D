module solver
use prog_bar
	
contains 
    
subroutine GS (A,B,X,ITER)    
implicit none 

integer::N,Iter,npass,i,j 
real*8::sum,Xold,ratio
real*8,allocatable::A(:,:),B(:),X(:)
real*8,parameter::err=1e-4
integer,parameter::iter_max=500
type (type_prog)::progress

    print *,'Iterating with GS solver ...'
	
	iter = 0

	N=size(B)
	
	do i=1,N
		X(i)=0.
	end do 

30	iter=iter+1 
  
	call progress%input(iter_max,20)
	call progress%output(iter)
	
	npass=0 
	do i=1,N
		sum=0.
		Xold=X(i)
		do j = 1, N
			sum=sum-A(i,j)*X(j)
		end do
		sum  = B(i) + sum
		X(i) = X(i) + sum/A(i,i) 
		if (X(i)/=0.) then
			ratio = abs((X(i)-Xold)/X(i))
		else 
			ratio = 0.
		end if
		if (ratio<err) then
			npass=npass+1
		end if
	end do 
	
	if (iter<iter_max) then
	   if (npass.NE.N) then
		  goto 30
	   end if 
	end if 
         
end subroutine GS

end module solver
         