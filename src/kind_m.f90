module kind_m
!    use, intrinsic :: iso_fortran_env ! f2008
    private 
    public kd 
    integer, parameter :: kd = kind(1.0d0)  ! kd = real64
end module kind_m