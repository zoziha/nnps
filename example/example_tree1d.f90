! tree1d has bugs in gfortran using release mode, which is related to gfortran,
! unfortunately, I can't register gcc-bugzilla account very well.
program example_tree1d

    use nnps_module, only: nnps_binarytree, wp
    use display_module, only: display
    implicit none

    type(nnps_binarytree) :: nnps
    real(wp), dimension(4) :: loc = [0.0_wp, 1.0_wp, 2.0_wp, 1.5_wp]
    integer, pointer :: pairs(:)
    real(wp), pointer :: rdxs(:)

    call nnps%init(loc, minval(loc), maxval(loc))
    call nnps%build()
    call nnps%query(0.6_wp, pairs, rdxs)

    print *, '*** kd-tree find (1D)'
    call display(pairs, 'pairs index:')
    call display(loc(pairs), 'pairs coordinates:')
    call display(rdxs, 'rdxs:')

end program example_tree1d
!  *** kd-tree find (1D)
! [vector: 4] pairs index:
! 2, 4, 3, 4
! [vector: 4] pairs coordinates:
!  1.000E+00,  1.500E+00,  2.000E+00,  1.500E+00
! [vector: 4] rdxs:
!  5.000E-01, -5.000E-01,  5.000E-01,  5.000E-01
