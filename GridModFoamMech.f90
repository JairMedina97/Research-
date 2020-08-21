subroutine GridModFoamMech
  use VarModFoamMech
  implicit none

     allocate  (x(iMax),t(kMax),P(iMax,kMax),Sw(iMax,kMax), &
               Sg(iMax,kMax),Fw(iMax,0:1), &
               pc(iMax,0:1),lamdaG(iMax,kMax),lamdaW(iMax,kMax), &
               Nug(iMax,kMax),rg(iMax,0:1),rc(iMax,0:1),        &
               nf(iMax,kMax),Uw(iMax,0:1),Ug(iMax,0:1),   &
               Krg0(iMax,0:1),Krw(iMax,0:1),    &
               nfaux(iMax),Swaux(iMax),Paux(iMax),rgaux(iMax), &
               NuGaux(iMax),&
               a(kMax), Kp(iMax),&
               Ugaux(iMax),Uwaux(iMax),rcaux(iMax),poro(iMax), stat=status)


     if(status/=0) then
       write (*,*)'MemoryError'
     end if

!$OMP PARALLEL
  !$OMP DO
     do i=1,iMax
       x(i) = rw+((i-1.d0)/(iMax-1.d0))*(1-rw)
     end do
  !$OMP END DO
     dx=x(2)-x(1)

  !$OMP DO
     do k=1,kMax
       t(k) = tMax*((k-1.d0)/(kMax-1.d0))
     end do
  !$OMP END DO
     dt=t(2)
!$OMP END PARALLEL

end subroutine GridModFoamMech
