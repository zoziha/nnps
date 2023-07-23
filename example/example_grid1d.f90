program example_grid1d

    use nnps_module, only: nnps_grid1d, rk
    use display_module, only: display
    implicit none

    type(nnps_grid1d) :: nnps
    real(rk), dimension(4) :: loc = [0.0_rk, 1.0_rk, 2.0_rk, 1.5_rk]
    integer, pointer :: pairs(:)
    real(rk), pointer :: rdxs(:)

    call nnps%init(loc, minval(loc), maxval(loc), 0.6_rk)
    call nnps%build()
    call nnps%query(0.6_rk, pairs, rdxs)

    print *, '*** grid find (1D)'
    call display(pairs, 'pairs index:')
    call display(loc(pairs), 'pairs coordinates:')
    call display(rdxs, 'rdxs:')

end program example_grid1d
!  *** grid find (1D)
! [vector: 4] pairs index:
! 4, 2, 3, 4
! [vector: 4] pairs coordinates:
!  1.500E+00,  1.000E+00,  2.000E+00,  1.500E+00
! [vector: 4] rdxs:
!  5.000E-01,  5.000E-01,  5.000E-01,  5.000E-01
