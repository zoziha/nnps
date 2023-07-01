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

    call nnps%init(loc, minval(loc, 2), maxval(loc, 2), 0.6_rk)
    call nnps%build()
    call nnps%query(0.6_rk, pairs)

    print *, '*** grid find (3D)'
    call display(real(pairs), 'pairs index:')
    call display(loc(:, pairs), 'pairs coordinates:')

end program example_grid3d
!> [vector: 2] pairs index:
!>  2.000E+00,  4.000E+00
!> [matrix: 3*2] pairs coordinates:
!>  1.500E+00,  1.000E+00;
!>  1.000E+00,  1.000E+00;
!>  1.000E+00,  1.000E+00
