subroutine SatModFoamMech
  use VarModFoamMech
  implicit none

    ! Boundary condition at x=0
      Sw(1,k+1) =   Sw(1,k) &
                   + dt*(FwJ-Fw(1,0))/(rw*dx)+&
                     dt*(1-Fw(1,0))*lamdaW(1,k)*(x(2)*pc(2,0)-x(1)*pc(1,0))/(rw**2*dx**2)
    do i=2,iMax-1
       Sw(i,k+1) =   Sw(i,k)+&
                    dt*(Fw(i-1,0)-Fw(i,0))/(x(i)*dx)+&
                    dt*(1-Fw(i,0))*lamdaW(i,k)*(x(i+1)*pc(i+1,0)-x(i  )*pc(i  ,0))/(x(i)*dx**2*rw)-&
                    dt*(1-Fw(i,0))*lamdaW(i,k)*(x(i  )*pc(i  ,0)-x(i-1)*pc(i-1,0))/(x(i)*dx**2*rw)
    end do

    ! Boundary condition at x=L
    UwN = Fw(iMax  ,0)*rw/x(iMax) - &
          lamdaW(iMax,k)*2*pc(iMax,0)/dx+ &
          lamdaW(iMax,k)*Fw(iMax,0)*2*pc(iMax,0)/dx
    if(UwN < 0.d0) UwN = 0.d0 !No inflow from the core end

    Sw(iMax,k+1) =  Sw(iMax,k)-&
                    dt*2*Fw(iMax,0)/(dx*x(iMax))+&
                    dt*Fw(iMax-1,0)/(dx*x(iMax))-&
                    dt*(1-Fw(iMax,0))*lamdaW(iMax,k)*(x(iMax)*pc(iMax,0)-x(iMax-1)*pc(iMax-1,0))/(x(iMax)*dx**2*rw)+&
                    dt*UwN/(dx*rw)

      do i=1,iMax
         Sg(i,k+1)=1-Sw(i,k+1)
      end do

end subroutine SatModFoamMech
