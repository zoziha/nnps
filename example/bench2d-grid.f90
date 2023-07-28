!> 2D nearest neighbor particle search (uniform distribution)
program main

    use random_module, only: randn, randu
    use nnps_module, only: rk, nnps_direct2d, nnps_quadtree, nnps_grid2d, nnps_grid2d_finalizer
    use timer_module, only: timer, sec2hms
    use display_module, only: display
    implicit none
    real(rk), dimension(:, :), allocatable :: loc
    integer :: m, loop, i
    type(timer) :: tmr
    type(nnps_grid2d) :: nnps_grid
    integer, pointer :: pairs(:)
    real(rk), pointer :: rdxs(:)
    real(rk) :: t1, t2

    m = 10000
    loop = 1000
    allocate (loc(2, m))
    call randu(loc, -100.0_rk, 100.0_rk)
    call nnps_grid%init(loc, m)

    call cpu_time(t1)
    call tmr%tic()
    do i = 1, loop
        call nnps_grid%query(1.0_rk, pairs, rdxs, m)
    end do

    print *, "*** grid2d ***"
    call display(sec2hms(tmr%toc()), 'time:', inline=.true.)
    call cpu_time(t2)
    call display(t2 - t1, 'cpu_time:', inline=.true.)
    call display(pairs, 'pairs:')
    call display(real(nnps_grid%storage())/(8*1024*1024), "storage (tbl/pairs/all):")
    call display(nnps_grid%tbl%activated_buckets(), "activated_buckets:")
    call nnps_grid2d_finalizer(nnps_grid)  ! finalize

end program main
!>  grid2d: 00:00:00.401
!> [vector: 76] pairs:
!>  1.610E+02,  8.700E+01,  7.210E+02, ...  2.170E+02
!>  tree2d: 00:00:02.043
!> [vector: 76] pairs:
!>  2.200E+01,  3.390E+02,  2.800E+01, ...  9.830E+02
!>  direct2d: 00:00:02.635
!> [vector: 76] pairs:
!>  2.200E+01,  3.390E+02,  2.800E+01, ...  9.830E+02
