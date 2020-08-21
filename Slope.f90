subroutine Slope
use VarModFoamMech
    implicit none

    search: do i=1,iMax-1
       slp=(nf(i+1,k+1)-nf(i,k+1))/dx
       if(slp<0)then
        a(k+1)=x(i)
        iAva=i
        exit  search
       end if
    end do  search


end subroutine Slope
