module wav_io
    use mod_mpg
    implicit none
    private
    public :: riff_chunk, fmt_chunk, data_chunk              ! type
    public :: check_riff_chunk, open_wav_file, close_wav_file, read_pcm_1frame, read_pcm0 ! subroutine
    integer, save :: ir
    !
    type :: fmt_chunk
        character (len = 4):: chunk_id  != 'fmt '
        integer :: ichunk_size
        integer :: iformat_type, ichannels
        integer :: isamples_per_sec
        integer :: ibytes_per_sec
        integer :: iblock_size, ibits_per_sample
    end type fmt_chunk
    !
    type :: fact_chunk
        character (len = 4) :: chunk_id  != 'fact'
        integer :: ichunk_size !( bytes )
    !! some unknown data
    end type fact_chunk
    !
    type :: data_chunk
        character (len = 4) :: chunk_id  != 'data'
        integer :: ichunk_size !( bytes )
    !! pcm data follows integer (kind = 2) :: idata(isize / 2)
    end type data_chunk
    !
    type :: riff_chunk
        character (len = 4) :: chunk_id  != 'riff'
        integer :: ichunk_size
        character (len = 4) :: format_type ! = 'wave'
        type ( fmt_chunk) :: fmt
        type (fact_chunk) :: fct
        type (data_chunk) :: dat
    end type riff_chunk
    !
contains
    !--------------------------------------------------------------------
    subroutine init_chunk_names(riff)
        type (riff_chunk), intent(out) :: riff
        riff%chunk_id     = 'RIFF'
        riff%format_type  = 'WAVE'
        riff%fmt%chunk_id = 'fmt '
        riff%fct%chunk_id = 'fact'
        riff%dat%chunk_id = 'data'
    end subroutine init_chunk_names
    !---------------------------------------------------------------------
    function word32() result(res)
        character (len = 4) :: res
        integer :: io
        read(ir, iostat = io) res
    end function word32
    !---------------------------------------------------------------------
    function word16() result(res)
        character (len = 2) :: res
        integer :: io
        read(ir, iostat = io) res
    end function word16
    !---------------------------------------------------------------------
    function word8() result(res)
        character (len = 1) :: res
        integer :: io
        read(ir, iostat = io) res
    end function word8
    !---------------------------------------------------------------------
    function int32() result(ires)
        integer (kind = 4) :: ires
        ires = transfer(word32(), ires) ! little endian assumed
    end function int32
    !---------------------------------------------------------------------
    function int16() result(ires)
        integer (kind = 2) :: ires
        ires = transfer(word16(), ires) ! little endian assumed
    end function int16
    !---------------------------------------------------------------------
    subroutine abort(text)
        character (len = *), intent(in) :: text
        write(*, *) 'abort:: ', text
        stop
    end subroutine abort
    !------------------------------------------------------------------
    subroutine open_wav_file(iread, fname)
        integer            , intent(in    ) :: iread
        character (len = *), intent(in    ) :: fname
        integer :: io
        ir = iread
        open(ir, file = fname, status = 'old', iostat = io, access = 'stream') 
        if (io /= 0) then
            write(*, '(a, i3, a, i3, 2a)' ) ' i/o error ', io, ' occuerred. file =', iread, ' file name ', fname
            call abort('check input file! suggestion: is file name correct?')
        end if
    end subroutine open_wav_file
    !------------------------------------------------------------------
    subroutine close_wav_file
        close(ir)
    end subroutine close_wav_file
    !------------------------------------------------------------------
    subroutine check_riff_chunk(riff)
        type (riff_chunk), intent(in out) :: riff
        call init_chunk_names(riff)
        if ( word32() == riff%chunk_id ) then    ! 'riff'?
            write(*, '(a)', advance = 'no') ' ms riff '
        else
            call abort('this is not ms-riff file!')
        end if
        riff%ichunk_size = int32()
        if ( word32() == riff%format_type ) then ! 'wave'?
            write(*, '(a)', advance = 'no') 'wav audio '
        else
            write(*, *)
            call abort('this is not wav file!')
        end if
        call check_fmt_chunk(riff%fmt)
        call check_dat_chunk(riff%dat, riff%fct)
    end subroutine check_riff_chunk
    !------------------------------------------------------------------
    subroutine check_fmt_chunk(fmt)
        type (fmt_chunk), intent(in out) :: fmt
        if ( word32() /= fmt%chunk_id ) call abort('cannot find format chunk!')
        fmt%ichunk_size      =     int32()
        fmt%iformat_type     = int(int16(), kind = 4)
        fmt%ichannels        = int(int16(), kind = 4)
        fmt%isamples_per_sec =     int32()
        fmt%ibytes_per_sec   =     int32()
        fmt%iblock_size      = int(int16(), kind = 4)
        fmt%ibits_per_sample = int(int16(), kind = 4)
        if ( fmt%iformat_type     /=  1) call abort('unknown wave format!') !linear pcm
        if ( fmt%ibits_per_sample /= 16) call abort('not 16bit data!')
        select case ( fmt%ichannels )
        case (1)
            write(*, '(a, i3, a, i6, a)', advance = 'no') &
             'monoral', fmt%ibits_per_sample, 'bit sampling rate', fmt%isamples_per_sec, 'hz '
        case (2)
            write(*, '(a, i3, a, i6, a)', advance = 'no') &
             'stereo' , fmt%ibits_per_sample, 'bit sampling rate', fmt%isamples_per_sec, 'hz '
        case default
            write(*, '(a, i1)') ' number of wave channels is ', fmt%ichannels
            call abort('wave channel must be 1 or 2!')
        end select
    end subroutine check_fmt_chunk
    !------------------------------------------------------------------
    subroutine check_dat_chunk(dat, fct)
        type (data_chunk), intent(in out) :: dat
        type (fact_chunk), intent(in out) :: fct
        integer :: i
        character (len = 4) :: chnk_id
        character (len = 1) :: dummy
        chnk_id = word32()
        if      ( chnk_id == fct%chunk_id ) then 
            fct%ichunk_size = int32()
            do i = 1, fct%ichunk_size
                dummy = word8()
            end do
            if ( word32() == dat%chunk_id ) then
                dat%ichunk_size = int32()
            end if
        else  if ( chnk_id == dat%chunk_id ) then
            dat%ichunk_size = int32()
        else
            call abort('cannot find fact chunk nor data chunk!')
        end if
    end subroutine check_dat_chunk
!------------------------------------------------------------------
    subroutine wav_read(pcm) ! 16bit pcm assumed
        real (kind = 8), intent(out) :: pcm(:, :)
        real (kind = 8), parameter   :: denom = 32768.0d0 !32768 = 2^15
        integer       , parameter   :: maxbuff = 1152 * 2
        character (len = 2) :: cbuff16(maxbuff)
        integer :: i, nchannel, ndat
        integer  (kind = 2) :: ibuff16(maxbuff)
        equivalence (cbuff16, ibuff16)
        ndat     = size(pcm, 1)
        nchannel = size(pcm, 2)
        if (ndat * nchannel > maxbuff) call abort('check maxbuff: subroutine wav_get')
        ibuff16 = 0
        select case (nchannel)
        case (1) !mono
            call wav_read_sub( cbuff16(1:ndat) )
            do i = 1, ndat
                pcm(i, 1) = real( ibuff16(i), kind = 8) / denom          ! little endian assumed
                pcm(i, 2) = 0.0d0
            end do
        case (2) !stereo
            call wav_read_sub( cbuff16(1:2 * ndat) )
            do i = 1, ndat
                pcm(i, 1) = real( ibuff16(2 * i - 1), kind = 8) / denom  ! little endian assumed
                pcm(i, 2) = real( ibuff16(2 * i    ), kind = 8) / denom  ! little endian assumed
            end do
        case default
            call abort('ichannel must be 1 or 2: subroutine wav_get')
        end select
    end subroutine wav_read
!------------------------------------------------------------------
    subroutine wav_read_sub(cha16)
        character (len = 2), intent(out) :: cha16(:)
        integer :: io
        read(ir, iostat = io) cha16
        select case (io)
        case (0)
            continue
        case (-1)
            continue
        !  write(*, *) 'end of file!'
        case default
            write(*, '(a, i3, a, i4)') ' file ', ir, ' iostat ', io
            call abort('i/o error occurred while reading wav file.')
        end select
    end subroutine wav_read_sub
    !------------------------------------------------------------------
    subroutine read_pcm_1frame(pcm)
        real (kind = 8), intent(out) :: pcm(:, :)
        pcm = eoshift(pcm, 1152, 0.0d0, 1)
        call wav_read(pcm(481:1632, :))
    end subroutine read_pcm_1frame
    !------------------------------------------------------------------
    subroutine read_pcm0(pcm)
        real (kind = 8), intent(out) :: pcm(:, :)
        call wav_read(pcm(1153:1632, :))
    end subroutine read_pcm0
    !------------------------------------------------------------------
end module wav_io
