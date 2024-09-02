module FFTW_LY
    use, intrinsic :: iso_c_binding
    implicit none 
    INCLUDE 'fftw3.f03'

    contains

    subroutine ZFFT1D(N,t,ft,wn,fw) !N should be an even number 
        use, intrinsic :: iso_c_binding
        implicit none 
        
        integer :: N, i
        real*8 :: t(N), wn(N), dw 
        complex(C_DOUBLE_COMPLEX):: ft(N),fw(N)
        type(C_PTR) :: plan
        
        call dfftw_plan_dft_1d(plan,N,ft,fw,FFTW_FORWARD,FFTW_ESTIMATE)
        call dfftw_execute_dft(plan,ft,fw)
        call dfftw_destroy_plan(plan)

        !normalize 
        fw = fw * abs(t(2) - t(1)) 
    
        !set the frequency domin
        dw = 2.d0*3.141592652359d0  / ( t(N)-t(1) + t(2)-t(1) )
        do i = 1,N/2
            wn(i) = (i-1)*dw
        end do
        do i = N/2+1,N 
            wn(i) = (i-1-N)*dw 
        end do
    end subroutine ZFFT1D

    subroutine ZIFFT1D(N,t,ft,wn,fw) !N should be an even number 
        use, intrinsic :: iso_c_binding
        implicit none 
        
        integer :: N, i
        real*8 :: t(N), wn(N) 
        complex(C_DOUBLE_COMPLEX):: ft(N),fw(N)
        type(C_PTR) :: plan
        
        call dfftw_plan_dft_1d(plan,N,fw,ft,FFTW_BACKWARD,FFTW_ESTIMATE)
        call dfftw_execute_dft(plan,fw,ft)
        call dfftw_destroy_plan(plan)
        
    end subroutine ZIFFT1D

    subroutine FFT_shift(N,wn,fw)
        implicit none 
        integer :: N, iw
        real*8 :: wn(N), wn_t(N)
        complex*16 :: fw(N),fw_t(N)

        wn_t = wn; fw_t = fw 
        wn(1:N/2) = wn_t(N/2+1:N); wn(N/2+1:N) = wn_t(1:N/2)
        fw(1:N/2) = fw_t(N/2+1:N); fw(N/2+1:N) = fw_t(1:N/2)
        
    end subroutine
    
end module FFTW_LY