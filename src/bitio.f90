module bit_io
    use mod_mpg
    use mod_crc
    use mod_layer3
    use mod_huffman
    implicit none
    private
    public :: open_mpg_file, close_mpg_file, clear_bit_buff, put_bits, put_bits_c &
            , put_bits_dim, put_bits_dim2, write_bits_1frame
    integer, save :: iw, ip
    character (len = 10000) :: bit_string
contains
!---------------------------------------------------------------------
    subroutine abend(text)
        character (len = *), intent(in) :: text
        write(*, *) 'abend:: ', text
        stop
    end subroutine abend
!-------------------------------------------------------------------
    subroutine open_mpg_file(iwrite, fname)
        integer, intent(in) :: iwrite
        character (len = *) :: fname
        integer :: io
        iw = iwrite
        open(iw, file = fname, iostat = io, status = 'unknown', access = 'stream') 
        if (io /= 0) then
            write(*, '(a, i3, a, i3, 2a)') ' i/o error ', io, ' occuerred. file =', iw, ' file name ', fname
            call abend('check output file! suggestion: file may be in use.')
        end if
    end subroutine open_mpg_file
!-------------------------------------------------------------------
    subroutine close_mpg_file
        close(iw)
    end subroutine close_mpg_file
!-------------------------------------------------------------------
    subroutine clear_bit_buff()
        ip = 1                      ! initialize buffer ; set position to 1 
        bit_string = repeat(' ', len(bit_string))
    end subroutine clear_bit_buff
!-------------------------------------------------------------------
    subroutine put_bits(n, inp) 
        integer, intent(in) :: n, inp
        integer :: i 
        if (n > 32) call abend('out of range: n must be < 32: subroutine put_bits')
        do i = 1, n
            if (ibits(inp, n - i, 1) == 1) then
                bit_string(ip:ip) = '1'
            else
                bit_string(ip:ip) = '0'
            end if
            ip = ip + 1
        end do
    end subroutine put_bits
!-------------------------------------------------------------------
    subroutine put_bits_dim(n, inp)
        integer, intent(in) :: n, inp(:)
        integer :: i
        do i = 1, size(inp)
            call put_bits(n, inp(i))
        end do
    end subroutine put_bits_dim
!-------------------------------------------------------------------
    subroutine put_bits_dim2(n, inp)
        integer, intent(in) :: n, inp(:, :)
        integer :: i, j
        do i = 1, size(inp, 1)
            do j = 1, size(inp, 2)
                call put_bits(n, inp(i, j))
            end do
        end do
    end subroutine put_bits_dim2
!-------------------------------------------------------------------
    subroutine put_bits_c(str)
        character (len = *) :: str
        integer :: i
        do i = 1, len_trim(str)
            if (str(i:i) /= '0' .and. str(i:i) /= '1') call abend('invalid string: subroutine put_bit_c')
            bit_string(ip:ip) = str(i:i)
            ip = ip + 1
        end do
    end subroutine put_bits_c
!-------------------------------------------------------------------
    subroutine write_bits_1frame(n)
        integer, intent(in) :: n
        integer :: i, j, ipos, m
        character(len = 4) ::cm
        equivalence (m, cm) ! integer*4 assumed for m
        if (mod(n, 8) /= 0) call abend('input error: n must be multiples of 8: subroutine write_bits_1frame')
        ipos = 0
        do i = 1, n, 8
            m = 0
            do j = 1, 8
                ipos = ipos + 1
                if (ipos > len(bit_string)) exit
                if (bit_string(ipos:ipos) == '1') m = m + 2**(8 - j)
            end do
            write(iw) cm(1:1)   ! little endian assumed
        end do
    end subroutine write_bits_1frame
!-------------------------------------------------------------------
end module bit_io