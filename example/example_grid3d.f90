program example_grid3d

    use nnps_module, only: nnps_grid3d, rk
    use display_module, only: display
    implicit none

    type(nnps_grid3d) :: nnps
    real(rk), dimension(3, 4) :: loc = reshape([0.0_rk, 1.0_rk, 2.0_rk, &
                                                1.5_rk, 1.0_rk, 1.0_rk, &
                                                0.5_rk, 1.0_rk, 0.0_rk, &
                                                1.0_rk, 1.0_rk, 1.0_rk], [3, 4])
    integer, pointer :: pairs(:)
    real(rk), pointer :: rdxs(:)

    call nnps%init(loc, minval(loc, 2), maxval(loc, 2), 0.6_rk)
    call nnps%build()
    call nnps%query(0.6_rk, pairs, rdxs)

    print *, '*** grid find (3D)'
    call display(pairs, 'pairs index:', brief=.false.)
    call display(loc(:, pairs), 'pairs coordinates:', brief=.false.)
    call display(rdxs, 'rdxs:', brief=.false.)

end program example_grid3d
!  *** grid find (3D)
! [vector: 2] pairs index:
! 4, 2
! [matrix: 3*2] pairs coordinates:
!  1.000E+00,  1.500E+00;
!  1.000E+00,  1.000E+00;
!  1.000E+00,  1.000E+00
! [vector: 4] rdxs:
!  5.000E-01, -5.000E-01,  0.000E+00,  0.000E+00
