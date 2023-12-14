program example_tree3d

    use nnps_module, only: nnps_octree, wp
    use display_module, only: display
    implicit none

    type(nnps_octree) :: nnps
    real(wp), dimension(3, 4) :: loc = reshape([0.0_wp, 1.0_wp, 2.0_wp, &
                                                1.5_wp, 1.0_wp, 1.0_wp, &
                                                0.5_wp, 1.0_wp, 0.0_wp, &
                                                1.0_wp, 1.0_wp, 1.0_wp], [3, 4])
    integer, pointer :: pairs(:)
    real(wp), pointer :: rdxs(:)

    call nnps%init(loc, minval(loc, 2), maxval(loc, 2))
    call nnps%build()
    call nnps%query(0.6_wp, pairs, rdxs)

    print *, '*** kd-tree find (3D)'
    call display(pairs, 'pairs index:', brief=.false.)
    call display(loc(:, pairs), 'pairs coordinates:', brief=.false.)
    call display(rdxs, 'rdxs:', brief=.false.)

end program example_tree3d
!  *** kd-tree find (3D)
! [vector: 2] pairs index:
! 2, 4
! [matrix: 3*2] pairs coordinates:
!  1.500E+00,  1.000E+00;
!  1.000E+00,  1.000E+00;
!  1.000E+00,  1.000E+00
! [vector: 4] rdxs:
!  5.000E-01,  5.000E-01,  0.000E+00,  0.000E+00
