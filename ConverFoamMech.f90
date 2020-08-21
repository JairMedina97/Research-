subroutine ConverFoamMech
   use VarModFoamMech
   implicit none

    !$OMP PARALLEL

        SwAverage = 0.d0
        !$OMP DO
        do i=2,iMax-1
           SwAverage = SwAverage + Sw(i,k+1)*x(i)
        end do
        SwAverage = SwAverage + (Sw(1,k+1)*x(1)+Sw(iMax,k+1)*x(iMax))/2.d0
        SwAverage = SwAverage*2*dx/(1.d0-rw**2)
        !$OMP END DO

        nfAverage = 0.d0
        !$OMP DO
        do i=2,iMax-1
            nfAverage = nfAverage + Nf(i,k+1)*x(i)
        end do
        nfAverage = nfAverage + (Nf(1,k+1)*x(1)+Nf(iMax,k+1)*x(iMax))/2.d0
        nfAverage = nfAverage*2*dx/(1.d0-rw**2)
        !$OMP END DO

        NuGAverage = 0.d0
        !$OMP DO
        do i=2,iMax-1
            NuGAverage = NuGAverage + Nug(i,k+1)*x(i)
        end do
        NuGAverage = NuGAverage + (Nug(1,k+1)*x(1)+Nug(iMax,k+1)*x(iMax))/2.d0
        NuGAverage = NuGAverage*2*dx/(1.d0-rw**2)
        !$OMP END DO

    !$OMP END PARALLEL
           nfAverage=sum(nf(1:iMax,k+1))/iMax

           ConverNuG=abs((NuGAverage-NuGAverageOld)/NuGAverageOld)
           NuGAverageOld=NuGAverage

           ConverSw=abs((SwAverage-SwAverageOld)/SwAverageOld)
           SwAverageOld=SwAverage

           if(nfAverageOld /= 0.d0) then
             ConverNf=abs((nfAverage-nfAverageOld)/nfAverageOld)
             nfAverageOld=nfAverage
           else
             ConverNf = ConverMax/2.d0

           end if
             ConverNuG = ConverMax/2.d0
!            ConverNf=abs((nfAverage-nfAverageOld)/nfAverageOld)
!             nfAverageOld=nfAverage

           conver = max(convernug, conversw, convernf)


end subroutine ConverFoamMech
