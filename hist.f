      module xhist
      implicit none
      public
      integer,parameter:: samp=200
      type s200
      real*8 ::x(samp),y(samp),z(samp,samp)
      end type s200
      type line200
      real*8 ::x(samp),y(samp)
      end type line200
      
      contains
      subroutine wrt_s(filein,s)
      type(s200) s
      character*200 filein
      integer :: i,j
      if (maxval(s%z) .ne. 0.d0) then
      s%z=s%z/maxval(s%z)
      open(file=filein,unit=35)
      do i=1,samp
      do j=1,samp
      write(35,'(3f20.6)') s%x(i),s%y(j),s%z(i,j)
      enddo
      enddo
      else 
      write(*,*)" 0 data ! "
      endif
      close(35)
      end subroutine

      ! save s in #filein in csv format
      subroutine wrt_sCsv(filein,s)
      type(s200) s
      character*200 filein
      integer :: i,j

      if (maxval(s%z) .ne. 0.d0) then
      s%z=s%z/maxval(s%z) ! normalize

      open(file=filein,unit=35)
      ! first  line
      print*,size(s%x)
      write(35,"(a)",advance='no')',' ! first comma
      do i=1,size(s%x)-1
          write(35,"(a)",advance='no')trim(r2c(s%x(i)))//','
      enddo
      write(35,"(a)",advance='no')trim(r2c(s%x(i)))
      write(35,'(a)')""
      !
      ! second line to last line
      do j=1,size(s%y)
          write(35,"(a)",advance='no')trim(r2c(s%y(j)))
          do i=1,size(s%x)
              write(35,"(a)",advance='no')','//trim(r2c(s%z(i,j)))
          enddo
          write(35,'(a)')""
      enddo
      close(35)
      else 
      write(*,*)" 0 data ! "
      endif
      end subroutine
      ! save line in #filein
      subroutine wrt_l(filein,l)
      type(line200) l
      character*200 filein
      integer :: i
      if (maxval(l%y) .ne. 0.d0) then
          l%y=l%y/maxval(l%y)
          open(file=filein,unit=35)
          do i=1,samp
              write(35,'(2f20.6)') l%x(i),l%y(i)
          enddo
          close(35)
      else
          write(*,*)" 0 data ! "
      endif
      end subroutine

      subroutine init_h2(a1,b1,a2,b2,s)
      type(s200) :: s
      real*8 :: a1,a2,b1,b2,intv1,intv2
      integer:: i,j,k
      s%z=0
      if( a1 .gt.b1 .or. a2.gt.b2 ) then
          write(*,*) "WRONG RANGE!"
          stop
      else
          intv1=(b1-a1)/samp
          s%x(1)=a1+0.5*intv1
          do i=2,samp
          s%x(i)=s%x(i-1)+intv1
          enddo
          intv2=(b2-a2)/samp
          s%y(1)=a2+0.5*intv2
          do i=2,samp
          s%y(i)=s%y(i-1)+intv2
          enddo
      endif
      end subroutine init_h2

      subroutine  checkS(s)
      type(s200) :: s
      integer:: nsize
      nsize=size(s%x)
      print*,"x:",s%x(1:3),s%x(nsize-2:nsize)
      print*,"y:",s%y(1:3),s%y(nsize-2:nsize)
      end subroutine  checkS
      

      subroutine h2(v1,v2,s)
      type(s200) :: s
      real*8 :: a1,a2,b1,b2,v1,v2,intv1,intv2
      integer :: vceil1,vceil2
      intv1 = s%x(2)-s%x(1)
      a1=s%x(1)-0.5*intv1
      b1=s%x(200)+0.5*intv1
      intv2 = s%y(2)-s%y(1)
      a2=s%y(1)-0.5*intv2
      b2=s%y(200)+0.5*intv2
!     write(*,*) a,b,intv
      if ( v1 .gt. b1 .or. v1 .lt. a1 .or. v2.gt.b2 .or. v2.lt.a2 ) then
          write(*,'(a/6f20.8)') "OUT RANGE h2 !",a1,v1,b1,a2,v2,b2
          stop
      else
        vceil1=ceiling((v1-a1)/intv1)
        vceil2=ceiling((v2-a2)/intv2)
        s%z(vceil1,vceil2)=s%z(vceil1,vceil2)+1
!       write(*,*) l%y(vceil)
      endif
      end subroutine h2
      subroutine init_hist(a,b,l)
      type(line200) :: l
      real*8 :: a,b,intv
      integer:: i,j,k
      l%y=0
      if( a .gt.b ) then
          write(*,*) "WRONG RANGE!"
          stop
      else
          ! set ints of line, plus 1.d-6 to avoid small float error 
          intv=(b-a)/samp
          l%x(1)=a+0.5*intv-1.d-6
          do i=2,200
          l%x(i)=l%x(i-1)+intv+1.d-6
          enddo
      endif
      end subroutine init_hist


      ! histogram v with line l
      subroutine hist(v,l)
      type(line200) :: l
      real*8 :: a,b,v,intv
      integer :: vceil
      intv = l%x(2)-l%x(1)
      a=l%x(1)-0.5*intv
      b=l%x(200)+0.5*intv
!     write(*,*) a,b,intv
      ! if v is outrange, throw a exception
      if ( v .gt. b .or. v .lt. a ) then
          write(*,*) "OUT RANGE!", a, v, b
       !  stop
      else
        vceil=ceiling((v-a)/intv)
      ! write(*,*)v , vceil
        l%y(vceil)=l%y(vceil)+1
!       write(*,*) l%y(vceil)
      endif
      end subroutine hist
        function r2c(r)
            implicit none
            real*8 r
            character*20 r2c
            write(r2c,'(f20.6)')r
            r2c=adjustl(r2c)
            return
        end function
      end module xhist
