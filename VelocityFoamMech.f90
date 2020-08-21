subroutine VelocityFoamMech
    use VarModFoamMech
    implicit none

      Uw(1,1) = FwJ*rw/x(1)
    ! The pressure-capilar gradient in x=0 is asumed zero
       Ug(1,1) =Ug(1,1)*relax+(1-relax)*( rw/x(1)-Uw(1,1))
    !$OMP PARALLEL
    !$OMP DO
    do i=2,iMax-1
      if(iter==1)then
          Uw(i,1) = Fw(i,1)*rw/x(i) + &
                     lamdaW(i,k+1)*(pc(i,1)-x(i-1)*pc(i-1,1)/x(i))/dx- &
                     lamdaW(i,k+1)*Fw(i,1)*(pc(i,1)-x(i-1)*pc(i-1,1)/x(i))/dx
      else
          Uw(i,1) = Uw(i,1)*relax+(1-relax)*( Fw(i,1)*rw/x(i) + &
                     lamdaW(i,k+1)*(pc(i,1)-x(i-1)*pc(i-1,1)/x(i))/dx- &
                     lamdaW(i,k+1)*Fw(i,1)*(pc(i,1)-x(i-1)*pc(i-1,1)/x(i))/dx)
          if(Fw(i,1)/=Fw(i,1) .or. lamdaG(i,k+1)/=lamdaG(i,k+1) ) then
              write(*,*)'velocity aqui',i,Fw(i,1),lamdaG(i,k+1),Pc(i,1)
              stop
           end if
      end if
    end do
    !$OMP END DO
    !$OMP DO
    do i=2,iMax-1
       if(iter==1)then
         Ug(i,1) =rw/x(i)-Uw(i,1)
       else
         Ug(i,1) =Ug(i,1)*relax+(1-relax)* (rw/x(i)-Uw(i,1))
       end if
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    if(iter==1)then
       Uw(iMax,1) = Fw(iMax,1)*rw + &
                      lamdaW(iMax,k+1)*(-x(iMax-1)*pc(iMax-1,1))/dx- &
                      lamdaW(iMax,k+1)*Fw(iMax,1)*(-x(iMax-1)*pc(iMax-1,1))/dx
                      !The pressure-capilar in x=L is asumed zero
       if ( Uw(iMax,1) < 0.d0 ) Uw(iMax,1) = 0.d0
       Ug(iMax,1) = rw/x(iMax)-Uw(iMax,1)
    else
       Uw(iMax,1) =Uw(iMax,1)*relax+(1-relax)* (Fw(iMax,1)*rw + &
                     lamdaW(iMax,k+1)*(-x(iMax-1)*pc(iMax-1,1))/dx- &
                     lamdaW(iMax,k+1)*Fw(iMax,1)*(-x(iMax-1)*pc(iMax-1,1))/dx)
                     !The pressure-capilar in x=L is asumed zero
       if ( Uw(iMax,1) < 0.d0 ) Uw(iMax,1) = 0.d0
       Ug(iMax,1) = Ug(iMax,1)*relax+(1-relax)*(rw/x(iMax)-Uw(iMax,1))
    end if

end subroutine VelocityFoamMech
