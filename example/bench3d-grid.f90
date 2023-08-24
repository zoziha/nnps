!> 3D nearest neighbor particle search (uniform distribution)
program main

    use random_module, only: randn, randu
    use nnps_module, only: rk, nnps_direct3d, nnps_octree, nnps_grid3d
    use timer_module, only: timer, sec2hms
    use display_module, only: display
    implicit none
    real(rk), dimension(:, :), allocatable :: loc
    integer :: m, loop, i
    type(timer) :: tmr
    type(nnps_grid3d) :: nnps_grid
    integer, pointer :: pairs(:)
    real(rk), pointer :: rdxs(:)
    real(rk) :: t1, t2

    m = 130000
    loop = 100
    allocate (loc(3, m))
    call randu(loc, -20.0_rk, 20.0_rk)
    call nnps_grid%init(loc, m)

    call cpu_time(t1)
    call tmr%tic()
    do i = 1, loop
        call nnps_grid%query(1.0_rk, pairs, rdxs, m)
    end do

    print *, "*** grid3d ***"
    call display(sec2hms(tmr%toc()), 'time:', inline=.true.)
    call cpu_time(t2)
    call display(t2 - t1, 'cpu_time:', inline=.true.)
    call display(pairs, 'pairs:')
    call display(real(nnps_grid%storage())/(8*1024*1024), "storage (tbl/all):")
    call display(nnps_grid%tbl%activated_buckets(), "activated_buckets (activated/recyclable):")

end program main
!>  *** grid3d ***
!> [scalar] time: '00:00:18.672' (186ms/timestep with 130k particles, 1.6GHz[i5-8250u], 1 core, 3D)
!> [scalar] cpu_time:  1.864E+01
!> [vector: 1073602] pairs:
!> 1, 37623, 1, ... 62281
!> [vector: 2] storage (tbl/all):
!>  2.112E+01,  4.567E+01
!> [vector: 2] activated_buckets (activated/recyclable):
!> 50095, 0
