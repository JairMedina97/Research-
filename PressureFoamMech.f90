subroutine PressureFoamMech
  use VarModFoamMech
  implicit none

     P(1,k+1) = relax*P(1,k+1) + (1-relax)*( PwJ - (dx*UwJ)/(2*lamdaW(1,k+1)) )
     !P(1,k+1) =  PwJ - (dx*UwJ)/(2*lamdaW(1,k+1))
     !P(2,k+1) =  P(1,k+1) - (dx*Uw(1,1))/lamdaW(1,k+1)
     P(2,k+1) = relax*P(2,k+1) + (1-relax)*( P(1,k+1) - (dx*Uw(1,1))/lamdaW(1,k+1) )

     do i=2,iMax-1
        P(i+1,k+1) = relax*P(i+1,k+1) + (1-relax)*( 4*P(i,k+1)/3 - P(i-1,k+1)/3 - (2*dx*Uw(i+1,1))/lamdaW(i+1,k+1)/3 )
        if(Uw(i+1,1)/=Uw(i+1,1) .or. Krw(i+1,1)/=Krw(i+1,1) ) then
           write(*,*)'pressure aqui',i,Uw(i+1,1),Krw(i+1,1)
           stop
        end if
     end do

!     P(1,k+1) =  PwJ - (dx*UwJ)/(2*lamdaW(1,k+1))
!     P(2,k+1) =  P(1,k+1) - (dx*Uw(1,1))/lamdaW(1,k+1)
!     do i=2,iMax-1
!        P(i+1,k+1) = 4*P(i,k+1)/3 - P(i-1,k+1)/3 - (2*dx*Uw(i+1,1))/lamdaW(i+1,k+1)/3
!        if(Uw(i+1,1)/=Uw(i+1,1) .or. Krw(i+1,1)/=Krw(i+1,1) ) then
!           write(*,*)'pressure aqui',i,Uw(i+1,1),Krw(i+1,1)
!           stop
!        end if
!     end do

!     P(1,k+1) = PwJ - (dx*UwJ)/(2*lamdaW(1,k+1))
!     do i=1,iMax-1
!        P(i+1,k+1) = P(i,k+1) - (dx*Uw(i,1))/lamdaW(i,k+1)
!        if(Uw(i,1)/=Uw(i,1) .or. Krw(i,1)/=Krw(i,1) ) then
!           write(*,*)'pressure aqui',i,Uw(i,1),Krw(i,1)
!           stop
!        end if
!     end do


     do i=1,iMax-1
       if (P(i+1,k+1) > P(i,k+1)) P(i+1,k+1)=P(i,k+1)
     end do
     deltap = PwJ - P(iMax,k+1)

end subroutine PressureFoamMech
