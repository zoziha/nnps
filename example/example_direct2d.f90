program example_direct2d

    use nnps_module, only: nnps_direct2d, rk
    use display_module, only: display
    implicit none

    type(nnps_direct2d) :: nnps
    real(rk), dimension(2, 4) :: loc = reshape([0.0_rk, 1.0_rk, 2.0_rk, 1.5_rk, &
                                                1.0_rk, 1.0_rk, 0.5_rk, 1.0_rk], [2, 4])
    integer, pointer :: pairs(:)
    real(rk), pointer :: rdxs(:)

    call nnps%init(loc)
    call nnps%query(0.6_rk, pairs, rdxs)

    print *, '*** direct find (2D)'
    call display(pairs, 'pairs index:')
    call display(loc(:, pairs), 'pairs coordinates:')
    call display(rdxs, 'rdxs:')

end program example_direct2d
!  *** direct find (2D)
! [vector: 4] pairs index:
! 1, 4, 3, 4
! [matrix: 2*4] pairs coordinates:
!  0.000E+00,  5.000E-01,  1.000E+00,  5.000E-01;
!  1.000E+00,  1.000E+00,  1.000E+00,  1.000E+00
! [vector: 6] rdxs:
!  5.000E-01, -5.000E-01,  0.000E+00, ...  0.000E+00
