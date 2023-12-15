!> 2D nearest neighbor particle search (uniform distribution)
program main

    use random_module, only: randn, randu
    use nnps_module, only: wp, nnps_direct2d, nnps_quadtree, nnps_grid2d
    use timer_module, only: timer, sec2hms
    use display_module, only: display
    implicit none
    real(wp), dimension(:, :), allocatable :: loc
    integer :: m, loop, i
    type(timer) :: tmr
    type(nnps_grid2d) :: nnps_grid
    integer, pointer :: pairs(:)
    real(wp), pointer :: rdxs(:)
    real(wp) :: t1, t2

    m = 130000
    loop = 100
    allocate (loc(2, m))
    call randu(loc, -100.0_wp, 100.0_wp)
    call nnps_grid%init(loc, n=m)

    call cpu_time(t1)
    call tmr%tic()
    do i = 1, loop
        call nnps_grid%query(1.0_wp, pairs, rdxs, m)
    end do

    print *, "*** grid2d ***"
    call display(sec2hms(tmr%toc()), 'time:', inline=.true.)
    call cpu_time(t2)
    call display(t2 - t1, 'cpu_time:', inline=.true.)
    call display(pairs, 'pairs:')

end program main
!>  *** grid2d *** (i5-8250U, 1 core, 2D)
!> [scalar] time: '00:00: 8.094'
!> [scalar] cpu_time:  8.094E+00
!> [vector: 1320832] pairs:
!> 1, 41589, 1, ... 53742
!>  *** grid2d *** (R5-2500U, 1 core, 2D)
!> [scalar] time: '00:00:17.219'
!> [scalar] cpu_time:  1.716E+01
!> [vector: 1318576] pairs:
!> 1, 96135, 2, ... 75984
