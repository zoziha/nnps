program example_tree1d

    use nnps_module, only: nnps_binarytree, rk
    use display_module, only: display
    implicit none

    type(nnps_binarytree) :: nnps
    real(rk), dimension(4) :: loc = [0.0_rk, 1.0_rk, 2.0_rk, 1.5_rk]
    integer, pointer :: pairs(:)

    call nnps%init(loc, minval(loc), maxval(loc))
    call nnps%build()
    call nnps%query(0.6_rk, pairs)

    print *, '*** kd-tree find (1D)'
    call display(real(pairs), 'pairs index:')
    call display(loc(pairs), 'pairs coordinates:')

end program example_tree1d
!> [vector: 4] pairs index:
!>  2.000E+00,  4.000E+00,  3.000E+00,  4.000E+00
!> [vector: 4] pairs coordinates:
!>  1.000E+00,  1.500E+00,  2.000E+00,  1.500E+00
