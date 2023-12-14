program example_grid2d

    use nnps_module, only: nnps_grid2d, wp
    use display_module, only: display
    implicit none

    type(nnps_grid2d) :: nnps
    real(wp), dimension(2, 4) :: loc = reshape([0.0_wp, 1.0_wp, 2.0_wp, 1.5_wp, &
                                                1.0_wp, 1.0_wp, 0.5_wp, 1.0_wp], [2, 4])
    integer, pointer :: pairs(:)
    real(wp), pointer :: rdxs(:)

    call nnps%init(loc, m=[8, 100], n=4)
    call nnps%query(0.6_wp, pairs, rdxs, n=4)

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