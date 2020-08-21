subroutine LamellaCreationCoalescenceFoamMech
    use VarModFoamMech
  implicit none

     rg(1,1) = Sw(1,k+1)*( ( PwJ-0.5d0*(P(1,k+1)+P(2,k+1)) ) /dx )**m
     rg(2,1) = Sw(2,k+1)*((P(1,k+1)-P(2,k+1))/dx)**m
     !rg(1,1) = rg(1,1)*relax + (1-relax)* (Sw(1,k+1)*( ( PwJ-0.5d0*(P(1,k+1)+P(2,k+1)) ) /dx )**m )
     !rg(2,1) = rg(2,1)*relax + (1-relax)* ( Sw(2,k+1)*((P(1,k+1)-P(2,k+1))/dx)**m )


    !$OMP PARALLEL
    !$OMP DO
    do i=3,iMax-1
      rg(i,1) = Sw(i,k+1)*((P(i-1,k+1)-P(i+1,k+1))/(2*dx))**m
      !rg(i,1) = rg(i,1)*relax + (1-relax)* (Sw(i,k+1)*( -(3*P(i,k+1)-4*P(i-1,k+1)+P(i-2,k+1))/dx/2 )**m)
    end do

    !$OMP END DO
    rg(iMax,1) = rg(iMax-1,1)

    !$OMP DO
    do i=1,iMax
      if  ( Sw(i,k+1) > SwAster ) then
         rc(i,1) = nf(i,k+1)*( 1.d0/(Sw(i,k+1)-SwAster) )**n
      else
         rc(i,1) = nf(i,k+1)*( 1.d0/eps )**n
      end if
    end do
    !$OMP END DO
!$OMP END PARALLEL
end subroutine LamellaCreationCoalescenceFoamMech
