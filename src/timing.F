      module timing

          type t_timer
              sequence
              integer*8 :: num = 0, clock
              logical   :: running = .false.
              real*8    :: total = 0
          end type t_timer

          type (t_timer) timers(6)

          integer*8          :: resolution

          integer, parameter :: ERESRESNONBOND = 1, ERESRESBOND = 2
          integer, parameter :: ERESSVN = 3, ESVNSVN = 4, ERESSOL = 5
          integer, parameter :: EMAIN = 6

          character(len=15)  :: timer_name(6) =
     &    (/"  RESRESNONBOND","     RESRESBOND","         RESSVN",
     &      "         SVNSVN","         RESSOL","           MAIN"/)

          contains 
              subroutine timer_start(timer)
                  integer, intent(in) :: timer

                  if(.not.timers(timer)%running) then
                      call SYSTEM_CLOCK(timers(timer)%clock,resolution)
                  endif

                  timers(timer)%running = .true.
              end subroutine

              subroutine timer_end(timer)
                  integer, intent(in) :: timer
                  integer*8           :: tmp_clock
                  real*8              :: tmp_time

                  if(timers(timer)%running) then
                      call SYSTEM_CLOCK(tmp_clock,resolution)
                      tmp_time = REAL(tmp_clock - timers(timer)%clock) / resolution
                      timers(timer)%total = timers(timer)%total + tmp_time
                      timers(timer)%num = timers(timer)%num + 1
                  endif

                  timers(timer)%running = .false.
                  timers(timer)%clock = 0
              end subroutine

              subroutine print_timing()
                  integer :: i
                  do i=1,6
                      write(*,"(a15,f8.3,i8,f8.3)") timer_name(i),
     &                timers(i)%total,timers(i)%num,1e6*(timers(i)%total/timers(i)%num)
                  end do
              end subroutine

      end module timing

