      module timing
          integer, parameter :: dp = SELECTED_REAL_KIND(15)

          type t_timer
              sequence
              integer*8 :: num = 0
              logical   :: running = .false.
              real*8    :: total = 0, start
          end type t_timer

          type (t_timer) timers(6)

          real*8          :: resolution

          integer, parameter :: ERESRESNONBOND = 1, ERESRESBOND = 2
          integer, parameter :: ERESSVN = 3, ESVNSVN = 4, ERESSOL = 5
          integer, parameter :: EMAIN = 6

          character(len=15)  :: timer_name(6) =
     &    (/"  RESRESNONBOND","     RESRESBOND","         RESSVN",
     &      "         SVNSVN","         RESSOL","           MAIN"/)

          contains 
              subroutine timer_setup()
                  integer*8 :: tmp_clock, tmp_res
                  call SYSTEM_CLOCK(tmp_clock, tmp_res)
                  resolution = 1d0 / tmp_res
              end subroutine


              subroutine timer_start(timer)
                  integer, intent(in) :: timer
                  integer*8 :: tmp_clock

                  if(.not.timers(timer)%running) then
                      call SYSTEM_CLOCK(tmp_clock)
                      timers(timer)%start = tmp_clock * resolution
                  endif

                  timers(timer)%running = .true.
              end subroutine


              subroutine timer_end(timer)
                  integer, intent(in) :: timer
                  integer*8           :: tmp_clock
                  real*8              :: tmp_time

                  if(timers(timer)%running) then
                      call SYSTEM_CLOCK(tmp_clock)
                      tmp_time = (DBLE(tmp_clock) * resolution) - timers(timer)%start
                      timers(timer)%total = timers(timer)%total + tmp_time
                      timers(timer)%num = timers(timer)%num + 1
                  endif

                  timers(timer)%running = .false.
                  timers(timer)%start = 0
              end subroutine


              subroutine print_timing()
                  character(len=300) :: printstring
                  integer            :: i
                  write(printstring,"(EN8.1E2)") resolution
                  call printLine(4,printstring)
                  do i=1,6
                      write(printstring,"(a15,f8.3,i8,f8.3)") timer_name(i),
     &                timers(i)%total,timers(i)%num,1e6*(timers(i)%total/timers(i)%num)
                      call printLine(4,printstring)
                  end do
              end subroutine

      end module timing

