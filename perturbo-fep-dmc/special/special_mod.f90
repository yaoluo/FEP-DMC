module special_func
       interface 
                subroutine erf_c(a,b)
                        use iso_c_binding
                        implicit none 
                        real(c_double) :: a,b
                end subroutine 
                subroutine inverf_c(a,b)
                        use iso_c_binding
                        implicit none
                        real(c_double) :: a,b
                end subroutine

        end interface 
 
end module
