program example_direct3d

    use nnps_module, only: nnps_direct3d, rk
    use display_module, only: display
    implicit none

    type(nnps_direct3d) :: nnps
    real(rk), dimension(3, 4) :: loc = reshape([0.0_rk, 1.0_rk, 2.0_rk, &
                                                1.5_rk, 1.0_rk, 1.0_rk, &
                                                0.5_rk, 1.0_rk, 0.0_rk, &
                                                1.0_rk, 1.0_rk, 1.0_rk], [3, 4])
    integer, pointer :: pairs(:)
    real(rk), pointer :: rdxs(:)

    call nnps%init(loc)
    call nnps%query(0.6_rk, pairs, rdxs)

    print *, '*** direct find (3D)'
    call display(pairs, 'pairs index:')
    call display(loc(:, pairs), 'pairs coordinates:')
    call display(rdxs, 'rdxs:')

end program example_direct3d
!  *** direct find (3D)
! [vector: 2] pairs index:
! 2, 4
! [matrix: 3*2] pairs coordinates:
!  1.500E+00,  1.000E+00;
!  1.000E+00,  1.000E+00;
!  1.000E+00,  1.000E+00
! [vector: 4] rdxs:
!  5.000E-01,  5.000E-01,  0.000E+00,  0.000E+00
