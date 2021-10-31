module mod_fft
    use kind_m
    implicit none
    private
    public :: init_fft, fft23, fft23han ! subroutine
    integer, public :: indx576(1152), indx192(384)
    complex (kind = kd), public :: omega576(0:1151), omega192(0:383), sqrt3_2
    real (kind = kd), public :: han576(1152), han192(384)
    real (kind = kd), parameter :: pi = 4 * atan(1.0_kd)
contains
    !-----------------------------------------------------------------------------------------------------------------
    subroutine init_fft
        integer :: i, k, n, m
        sqrt3_2 = cmplx(0.0_kd, sqrt(0.75_kd), kind = kd) ! isqrt(3) / 2
        !
        indx576 = 1
        do i = 1, 1152
            do k = 1, 7
                n = 2**(k - 1)
                m = 1152 / 2**k
                indx576(i) = indx576(i) + ( mod(i - 1, 2 * m) / m ) * n
            end do
            do k = 1, 2
                n = 2**6 * 3**(k - 1) 
                m = 1152 / 2**6 / 3**k
                indx576(i) = indx576(i) + ( mod(i - 1, 3 * m) / m ) * n
            end do
            omega576(i - 1) = exp( cmplx(0.0_kd, 2.0_kd * pi / 1152.0_kd * real(i - 1, kind = kd), kind = kd) )
        ! han576(i) = sqrt(8.0d0 / 3.0d0) * 0.5d0 *( 1.0d0 - cos(2.0d0 * pi * real(i - 1, kind = 8) / 1152.0d0 ) )
            han576(i) = 0.5_kd *( 1.0_kd - cos(2.0_kd * pi * real(i - 1, kind = 8) / 1152.0_kd ) )
        end do
        ! 
        indx192 = 1
        do i = 1, 192
            do k = 1, 6
                n = 2**(k - 1)
                m = 384 / 2**k
                indx192(i) = indx192(i) + ( mod(i - 1, 2 * m) / m ) * n
            end do
            do k = 1, 1
                n = 2**6 * 3**(k - 1) 
                m = 384 / 2**6 / 3**k
                indx192(i) = indx192(i) + ( mod(i - 1, 3 * m) / m ) * n
            end do
            omega192(i - 1) = exp( cmplx(0.0_kd, 2.0_kd * pi / 384.0_kd * real(i - 1, kind = kd), kind = kd) )
            ! han192(i) = sqrt(8.0d0 / 3.0d0) * 0.5d0 *( 1.0d0 - cos(2.0d0 * pi * real(i - 1, kind = 8) / 384.0d0 ) )
            han192(i) = 0.5_kd *( 1.0_kd - cos(2.0_kd * pi * real(i - 1, kind = kd) / 384.0_kd ) )
        end do
    end subroutine init_fft
    !-----------------------------------------------------------------------------------------------------------------
    subroutine fft23(np2, np3, indx, omega, fft)
        integer           , intent(in    ) :: np2, np3, indx(:)
        complex (kind = kd), intent(in    ) :: omega(0:)
        complex (kind = kd), intent(in out) :: fft(:)
        complex (kind = kd) :: c1, c2, c3, c4, tmp1, tmp2, tmp3
        integer :: i, j, nn, iphase1, iphase2, m1, m2, m3, k3, kn3, k2, kn2
        nn = 2**np2 * 3**np3
        fft = fft(indx) / real(nn, kind = 8) ! reorder and normalize
        ! 3**np3
        do k3 = 1, np3 ! 3^n (n=2)
            kn3 = 3**(k3 - 1)
            do i = 1, nn, 3 * kn3 
                do j = 1, kn3
                    iphase1 = 2**np2 * 3**(np3 - k3) * (j - 1) 
                    iphase2 = 2 * iphase1
                    c1 = omega( mod(iphase1, nn) )
                    c2 = omega( mod(iphase2, nn) )
                    m1 = i + j - 1
                    m2 = m1 + kn3
                    m3 = m2 + kn3
                    tmp1 =      fft(m1)
                    tmp2 = c1 * fft(m2)
                    tmp3 = c2 * fft(m3)
                    fft(m1) = tmp1 + tmp2 + tmp3
                    c3 = tmp1 - 0.5d0 * ( tmp2 + tmp3 )
                    c4 =      sqrt3_2 * ( tmp2 - tmp3 ) ! sqrt3_2 = i sqrt(3) / 2
                    fft(m2) = c3 + c4
                    fft(m3) = c3 - c4
                end do
            end do
        end do
        ! 2**np2
        do k2 = 1, np2
            kn2 = 2**(k2 - 1) * 3**np3
            do i = 1, nn, 2 * kn2
                do j = 1, kn2  
                    iphase1 = 2**(np2 - k2) * (j - 1) 
                    c1 = omega( mod(iphase1, nn) )
                    m1 = i + j - 1
                    m2 = m1 + kn2
                    tmp2 = c1 * fft(m2)
                    fft(m2) = fft(m1) - tmp2
                    fft(m1) = fft(m1) + tmp2
                end do
            end do
        end do
    end subroutine fft23
    !-----------------------------------------------------------------------------------------------------------------
    subroutine fft23han(np2, np3, indx, omega, fft, han)
        integer            , intent(in    ) :: np2, np3, indx(:)
        complex (kind = kd), intent(in    ) :: omega(0:)
        complex (kind = kd), intent(in out) :: fft(:)
        real    (kind = kd), intent(in    ) :: han(:)
        fft = fft * han
        call fft23(np2, np3, indx, omega, fft)
    end subroutine fft23han
    !-----------------------------------------------------------------------------------------------------------------
end module mod_fft
