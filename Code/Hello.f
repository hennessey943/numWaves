c23456789
      program problemSet1Fortran
	 
	 integer i, n
	 real y
	 write(*,*) 'HELLO FROM FORTRAN'
	 do 10 i = 1, 10
          y=i
          write(*,*) i,'', 1/y,'', sin(1/y)
10    continue

      end
