module VarModFoamMech
    implicit none
    integer(kind=4)                ::  i,kper,k,iMax,kMax,status,iter,iterMax,nfur,nMax,nn,iAva,iAvaOld,mark, &
                                       kcount

    real(kind=8)                   ::  phy,dt,dx,Nug0,Nuw0,PermAbs,Swc,Sgr,Cf,m,n,CgCc, &
                                       SwAster,tMax,L,Sigma,n1Pc,n2Pc,NuGAverage,NuGAverageOld, &
                                       SwAverage,SwAverageOld,nfAverage,nfAverageOld,ConverSw, &
                                       ConverNf,convernug,eps,an,a0,bn,w0,pi,difNf,difSw,AveNforig, &
                                       epsEst,AveSworig,slp,Da,Dam,Ga,Gam,Ca,qt,h,uc,Pcc,rcc,rgc,deltaPc,&
                                       tc,rw,Br,div,nfMax,icont,icontMax,DesStanPoro,MediaPoro,nPer

    real(kind=8)                   ::  DeltaP,NugF,Cg,Cc,Krg00,Krw0,nkrw,nkrg,gx,relax,&
                                       dRho,FwJ,UwN,UwJ,PwJ,Poutlet,Swini,nfini,ConverMax,Conver, &
                                       AveSw,AveNf,AveUg,AveUw,Averg,Averc,x1,x2,x3,nfper,gama,Br2, &
                                       Sgcr


    real(kind=8), allocatable      ::  x(:),t(:),nfaux(:),Swaux(:),Paux(:),rgaux(:),NuGaux(:),&
                                       P(:,:),Sw(:,:),Sg(:,:),nf(:,:),NuG(:,:),&
                                       Fw(:,:),a(:),kp(:), &
                                       pc(:,:),lamdaG(:,:),lamdaW(:,:), &
                                       rg(:,:),Krg0(:,:),Krw(:,:),Uw(:,:),Ug(:,:),&
                                       rc(:,:),Ugaux(:),Uwaux(:),rcaux(:),Poro(:)

end module VarModFoamMech


