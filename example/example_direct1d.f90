program example_direct1d

    use nnps_module, only: nnps_direct1d, wp
    use display_module, only: display
    implicit none

    type(nnps_direct1d) :: nnps
    integer, parameter :: n = 4                            !! 粒子数
    real(wp), dimension(n) :: loc = [0.0_wp, 1.0_wp, 2.0_wp, 1.5_wp]
    integer, pointer :: pairs(:)
    real(wp), pointer :: rdxs(:)

    call nnps%init(loc, n=n)
    call nnps%query(0.6_wp, pairs, rdxs, n=n)

    print *, '*** direct find (1D)'
    call display(pairs, 'pairs index:')
    call display(loc(pairs), 'pairs coordinates:')
    call display(rdxs, 'rdxs:')

end program example_direct1d
!  *** direct find (1D)
! [vector: 4] pairs index:
! 2, 4, 3, 4
! [vector: 4] pairs coordinates:
!  1.000E+00,  1.500E+00,  2.000E+00,  1.500E+00
! [vector: 4] rdxs:
!  5.000E-01, -5.000E-01,  5.000E-01,  5.000E-01
