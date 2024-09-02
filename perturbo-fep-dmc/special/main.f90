program main
        use special_func
        implicit none 
        real*8 :: x,y
        x = 0.3;y = 0;
        call erf_c(x,y)
        write(*,'(A20,2E15.5)')'x,erf(x) = ',x,y
        call inverf_c(x,y)
        write(*,'(A20,2E15.5)')'x,inverf(x) = ',x,y

end program
