program example_direct2d

    use nnps_module, only: nnps_direct2d, wp
    use display_module, only: display
    implicit none

    type(nnps_direct2d) :: nnps
    integer, parameter :: n = 4                            !! 粒子数
    real(wp), dimension(2, n) :: loc = reshape([0.0_wp, 1.0_wp, 2.0_wp, 1.5_wp, &
                                                1.0_wp, 1.0_wp, 0.5_wp, 1.0_wp], [2, 4])
    integer, pointer :: pairs(:)
    real(wp), pointer :: rdxs(:)

    call nnps%init(loc, n=n)
    call nnps%query(0.6_wp, pairs, rdxs, n=n)

    print *, '*** direct find (2D)'
    call display(pairs, 'pairs index:', brief=.false.)
    call display(loc(:, pairs), 'pairs coordinates:', brief=.false.)
    call display(rdxs, 'rdxs:', brief=.false.)

end program example_direct2d
!  *** direct find (2D)
! [vector: 4] pairs index:
! 1, 4, 3, 4
! [matrix: 2*4] pairs coordinates:
!  0.000E+00,  5.000E-01,  1.000E+00,  5.000E-01;
!  1.000E+00,  1.000E+00,  1.000E+00,  1.000E+00
! [vector: 6] rdxs:
!  5.000E-01, -5.000E-01,  0.000E+00, ...  0.000E+00
