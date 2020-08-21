subroutine FractionalFoamMech
  use VarModFoamMech
  implicit none
!$OMP PARALLEL
    !$OMP DO
    do i=1,iMax
      Fw(i,1) = lamdaW(i,k+1) / (lamdaW(i,k+1) + lamdaG(i,k+1))
    end do
    !$OMP END DO
!$OMP END PARALLEL
end subroutine FractionalFoamMech
