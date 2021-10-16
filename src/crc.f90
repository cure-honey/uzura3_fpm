module mod_crc
    implicit none
    private
    public :: crc16, crc16dim
    integer, parameter :: igenerator = 32773 ! z'8005' !b'1000000000000101'  ! x^16+x^15+x^2+1
contains
!--------------------------------------------------------------------------
    subroutine crc16(n, in, icrc)         ! iso 2.4.3.1, table a.9, table b.5 
        integer, intent(in    ) :: n, in
        integer, intent(in out) :: icrc
        integer :: j, ibit1, ibit2
        do j = n - 1, 0, -1
            ibit1 = ibits(in  ,  j, 1)           ! jth bit                   bit [31......0]
            ibit2 = ibits(icrc, 15, 1)           ! sign bit of 16bit crc
            icrc  = ishft(ibits(icrc, 0, 15), 1) ! shift up 1bit 16bit crc
            if (ieor(ibit1, ibit2) == 1) icrc = ieor(icrc, igenerator)
        end do
    end subroutine crc16
!--------------------------------------------------------------------------
    subroutine crc16dim(n, in, icrc)
        integer, intent(in    ) :: n, in(:)
        integer, intent(in out) :: icrc
        integer :: i
        do i = 1, size(in)
            call crc16(n, in(i), icrc)
        end do
    end subroutine crc16dim
!--------------------------------------------------------------------------
end module mod_crc