program example_grid2d

    use nnps_module, only: nnps_grid2d, rk
    use display_module, only: display
    implicit none

    type(nnps_grid2d) :: nnps
    real(rk), dimension(2, 4) :: loc = reshape([0.0_rk, 1.0_rk, 2.0_rk, 1.5_rk, &
                                                1.0_rk, 1.0_rk, 0.5_rk, 1.0_rk], [2, 4])
    integer, pointer :: pairs(:)
    real(rk), pointer :: rdxs(:)

    call nnps%init(loc, n=4)
    call nnps%query(0.6_rk, pairs, rdxs, n=4)

    print *, '*** grid find (2D)'
    call display(pairs, 'pairs index:', brief=.false.)
    call display(loc(:, pairs), 'pairs coordinates:', brief=.false.)
    call display(rdxs, 'rdxs:', brief=.false.)

end program example_grid2d
!  *** grid find (2D)
! [vector: 4] pairs index:
! 3, 4, 4, 1
! [matrix: 2*4] pairs coordinates:
!  1.000E+00,  5.000E-01,  5.000E-01,  0.000E+00;
!  1.000E+00,  1.000E+00,  1.000E+00,  1.000E+00
! [vector: 6] rdxs:
!  5.000E-01,  5.000E-01,  0.000E+00,  5.000E-01,  5.000E-01,  0.000E+00