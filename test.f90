  ! static distance from N to Ti
    subroutine staticR()
        use xhist
        implicit none
        integer i,j,k
        character*200 fin
        type(line200) lR
        fin="dDist.dat"
        ! initialize
        call init_hist(minval(rNxTi),maxval(rNxTi),lR)
        do i=1, nStatic
            call hist(rNxTi(i),lR)
        enddo
        call wrt_l(fin,lR)
    end  subroutine staticR
