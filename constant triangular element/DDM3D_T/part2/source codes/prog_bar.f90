Module prog_bar
  Implicit None
  private

  Type , public :: type_prog
    Integer , public :: N , lens , i
    Character :: M = '#' , O = '.'
    Character(len=100) :: message
  Contains
    Procedure :: input
	Procedure :: output
  End Type type_prog
  
  contains

  Subroutine input(sub,N,L)
    Class(type_Prog)::sub
    Integer,Intent(in)::N,L
    sub % N    = N
    sub % lens = L
    sub % i = 0
    sub % message = ' Progress: ' 
  End Subroutine input
  
  
  Subroutine output(sub,K)
    Class(type_Prog)::sub
    Integer , Intent(in)::K
    Character(len=1)::br
	character(len=20)::char1,char2
    integer::cn
	
    sub % i = k 
    if ( sub % i > sub % n ) sub % i = sub % n    
    cn = Nint( real( sub%i * sub%lens ) / real( sub%N ) )
    if ( sub%i < sub%n ) then
      br = char(13)
    else
      br = char(10)
	end if
	write (char1,*)sub%i
	write (char2,*)sub%N
	char1=adjustr(char1)
	char2=adjustl(char2)
	
	write( * , '(a4,9a\)')'    ', trim(sub%message) , '[' , &
	repeat(sub%M , cn ) , repeat( sub%O , sub%lens-cn ) , '] ' , char1,'/', char2 , br
  End Subroutine output
  
    
End Module prog_bar
