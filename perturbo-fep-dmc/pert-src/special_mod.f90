module special_func
       interface 
                subroutine erf_c(x,y)
                        use iso_c_binding
                        implicit none 
                        real(c_double) :: x,y
                end subroutine 
                subroutine inverf_c(x,y)
                        use iso_c_binding
                        implicit none
                        real(c_double) :: x,y
                end subroutine

        end interface 
 
end module
