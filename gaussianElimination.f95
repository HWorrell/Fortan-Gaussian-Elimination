!Please Note: In order to run a fortran OpenMP program, the environment variable
!OMP_NUM_THREADS must be set to the number of threads to be used.  I have been
!running the command 'export OMP_NUM_THREADS=10' to set my thread numbers
!through a system call.  Should this not work, please set the environment
!variable before running.
!
!Thank you

program gaussianElimination

implicit none

integer :: x, y, count, n, rowselect

real :: scalar

real(kind=4), dimension(10, 10) :: matrix

real(kind=4), dimension(10) :: rightside, solutions

!this is the number of rows & columns
n = 10

x = 1

y = 1

do while(x <= n)

        y = 1

        do while(y <= n)

!this code generates values from -10 (inclusive) to 10 (inclusive)
!the basis for this code came from the new mexico tech computer program
!website, at http://infohost.nmt.edu/tcc/help/lang/fortran/scaling.html
                matrix(x, y) = int(rand(0)*(10+1+10))-10

                rightside(x) = int(rand(0)*(10+1+10))-10

                count = count + 1

                y = y + 1

        end do

        x = x + 1

end do

x = 1

count = 0

y = 1

print *, "Original Matrix:"

call printmatrix(matrix, n)

print *, "Beginning Elimination"

call elimination(matrix, n)

print *, "End Elimination"

call backsubstitution(matrix, rightside, solutions, n)

call printmatrix(matrix, n)

Print *, "Below is the solutions to each row: "

do count = 1, 10

print *, rightside(count)

end do

print *, "The solution set for this matrix is: "

do count = 1, 10

print *, solutions(count)

end do

end program gaussianElimination






subroutine printmatrix(matrix, n)

implicit none

integer, intent(in) :: n

integer :: x, y

real, intent(in) :: matrix(10, 10)

x = 1

do while(x <= n)

                write(*, 20) matrix(x, :)

        x = x + 1

end do

20      format(10f10.5)

end subroutine



subroutine elimination(matrix, n)

use omp_lib

implicit none

integer :: row, rowadd, column, test, myid

real, intent(inout) :: matrix(10, 10)

real(kind=4), dimension(10) :: swapvector 
 
real :: scalar

integer, intent(in) :: n

logical :: bool

row = 1

column = 1

do while(row <= n)

        rowadd = row + 1

        !if a zero exists, swap rows
        if(matrix(row, column) == 0) then

                swapvector = matrix(row, :)

                matrix(row, :) = matrix(row + 1, :)

                matrix(row + 1, :) = swapvector

        end if

        if(matrix(row, column) /= 1) then

                scalar = (1/matrix(row, column))

                matrix(row, :) = matrix(row, :) * scalar

        end if

!        print *, "Post reduction"

!        call printmatrix(matrix, n)

        !make all numbers below it in the column a zero

!turn this loop into parallel        do while(rowadd <= n)


        call system("export OMP_NUM_THREADS=10")
        !$OMP PARALLEL DEFAULT(NONE) NUM_THREADS(10) SHARED(matrix, n, row, column) PRIVATE(myid)

        myid = OMP_GET_THREAD_NUM() + 1

        if(myid <= n .and. myid > row) then
                
                matrix(myid, :) = matrix(myid, :) - matrix(row, :) * matrix(myid, column)
        
        end if

        !$OMP END PARALLEL

        !move to the next diagonal

        row = row + 1

        column = column + 1

end do

end subroutine



subroutine backsubstitution(matrix, rightside, solutions, n)

use omp_lib

implicit none

real, intent(in) ::matrix(10, 10)

integer, intent(in) :: n

real, intent(inout) :: rightside(10), solutions(10)

integer :: row, col, myid

row = n

do while(row >= 1)

solutions(row) = rightside(row)

do col = row + 1, n

solutions(row) = solutions(row) - (matrix(row, col) * solutions(col))

end do

solutions(row) = solutions(row) / matrix(row, row)

row = row - 1

end do

end subroutine

