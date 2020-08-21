subroutine SatModFoamMechIter
    use VarModFoamMech
    implicit none

  ! Boundary condition at x=0
      Sw(1,k+1) =    Sw(1,k) &
                   + dt*gama*(FwJ-Fw(1,0))/(rw*dx)+&
                     dt*(1-gama)*(FwJ-Fw(1,1))/(rw*dx)+&
                     dt*gama*(1-Fw(1,0))*lamdaW(1,k)*(x(2)*pc(2,0)-x(1)*pc(1,0))/(rw**2*dx**2)+&
                     dt*(1-gama)*(1-Fw(1,1))*lamdaW(1,k+1)*(x(2)*pc(2,1)-x(1)*pc(1,1))/(rw**2*dx**2)

    do i=2,iMax-1
        Sw(i,k+1) = Sw(i,k) &
                    +dt*gama*(Fw(i-1,0)-Fw(i,0))/(x(i)*dx) &
                    +dt*(1-gama)*(Fw(i-1,1)-Fw(i,1))/(x(i)*dx) &
                    +dt*gama*(1-Fw(i,0))*lamdaW(i,k)*(x(i+1)*pc(i+1,0)-x(i  )*pc(i  ,0))/(x(i)*dx**2*rw) &
                    -dt*gama*(1-Fw(i,0))*lamdaW(i,k)*(x(i  )*pc(i  ,0)-x(i-1)*pc(i-1,0))/(x(i)*dx**2*rw) &
                    +dt*(1-gama)*(1-Fw(i,1))*lamdaW(i,k+1)*(x(i+1)*pc(i+1,1)-x(  i)*pc(  i,1))/(x(i)*dx**2*rw) &
                    -dt*(1-gama)*(1-Fw(i,1))*lamdaW(i,k+1)*(x(i  )*pc(i  ,1)-x(i-1)*pc(i-1,1))/(x(i)*dx**2*rw)

                   if(Sw(i,k+1)/=Sw(i,k+1)) then
                      write(*,*)'Sw',i
                      stop
                   end if
    end do

    ! Boundary condition at x=L
    UwN = Fw(iMax,1)*rw/x(iMax) + &
          lamdaW(iMax,k+1)*(-2*pc(iMax,1))/(dx*rw)- &
          lamdaW(iMax,k+1)*Fw(iMax,1)*(-2*pc(iMax,1))/(dx*rw)
          if(UwN < 0.d0) UwN = 0.d0 !No inflow from the core end

    Sw(iMax,k+1)= Sw(iMax,k)-&
                  dt*    gama*2*Fw(iMax,0)/(dx*x(iMax))-&
                  dt*(1-gama)*2*Fw(iMax,1)/(dx*x(iMax))+&
                  dt*    gama*Fw(iMax-1,0)/(dx*x(iMax))+&
                  dt*(1-gama)*Fw(iMax-1,1)/(dx*x(iMax))-&
                  dt*    gama*(1-Fw(iMax,0))*lamdaW(iMax,k  )*(x(iMax)*pc(iMax,0)-x(iMax-1)*pc(iMax-1,0))/(x(iMax)*dx**2*rw)-&
                  dt*(1-gama)*(1-Fw(iMax,1))*lamdaW(iMax,k+1)*(x(iMax)*pc(iMax,1)-x(iMax-1)*pc(iMax-1,1))/(x(iMax)*dx**2*rw)+&
                  dt*    gama*UwN/(dx*rw)+ &
                  dt*(1-gama)*UwN/(dx*rw)

   Sw(iMax,k+1) = Sw(iMax-1,k+1)

   do i=1,iMax
        if(Sw(i,k+1) < Swc ) Sw(i,k+1) = Swc
         Sg(i,k+1)=1-Sw(i,k+1)
   end do

 end subroutine SatModFoamMechIter
