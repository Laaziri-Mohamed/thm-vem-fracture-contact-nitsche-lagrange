ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c

      subroutine solveur_iteratif(NbCV,MemBloc,NbIncPrim,
     &       Aij,nligne_Aij,ncol_Aij,ind_ligne_Aij,
     &       Smm,Soll,critere_arret_gmres,MAXL,IPREC,
     &       ITER,IERR,ERR)



	implicit double precision (a-h,o-z)
	implicit integer (i-n)

cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        include 'include_parameter'
c
c       Nb max de vecteurs pour KRYLOV = MAXL  
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     IPREC = 0 si ILU0 seul sinon combinative avec AMG+ILU0
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
        external MATVEC, MSOLVE, D1MACH
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     JACOBIENNE PAR BLOCS 
c
c
c     structure bloc de la jacobienne: on stocke tous les zeros  
c     lies aux decentrage de facon a avoir la structure bloc maille complete 
c
        dimension Aij(MemBloc,NbIncPrim,NbIncPrim)
c
c
        dimension ncol_Aij(MemBloc)
        dimension nligne_Aij(MemBloc)
        dimension ind_ligne_Aij(NbCV+1)
c
c
c     scaling par bloc a gauche 
c
        double precision, dimension(:,:,:), allocatable :: DiagScaling 
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Sm  
c
        double precision, dimension(:), allocatable :: Sm 
        dimension Smm(NbCV,NbIncPrim)
c
c
c     Sol 
c
        double precision, dimension(:), allocatable :: Sol
        dimension Soll(NbCV,NbIncPrim)
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     VARIABLES LOCALES 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Systeme par blocs d'inconnues 
c
c
        double precision, dimension(:), allocatable :: Asyst
        integer, dimension(:), allocatable :: ncol_Asyst
        integer, dimension(:), allocatable :: nligne_Asyst
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        integer, dimension(:), allocatable :: nrow
        integer, dimension(:), allocatable :: ncol
c
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     GMRES 
c
c     tableaux a passer a MSOLVE (prec)
c
        double precision, dimension(:), allocatable :: RWORK
        integer, dimension(:), allocatable :: IWORK
c
c     LIGW = 20
c     LRGW = 1+nsolve*(MAXL+6)+MAXL*(MAXL+3)
c
c     MAXL = nb max de vecteurs pour KRYLOV 
c
c
        double precision, dimension(:), allocatable :: RGWK
        integer IGWK(20)
c
c     scaling for residual (SB) and solution (SX)     
c
c
        double precision, dimension(:), allocatable :: SB
        double precision, dimension(:), allocatable :: SX

        dimension AA(NbIncPrim,NbIncPrim), XX(NbIncPrim)
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc

c        write(*,*)' NbCVI ',NbCV
c
c        allocate(DiagScaling(NbCV,NbIncPrim,NbIncPrim))
c
ccccccccccccccccccc

c
c     scaling diagonal par bloc a gauche 
c
c
c
c        do k=1,NbCV
c           mk = ind_ligne_Aij(k)+1 
c           DiagScaling(k,1,1) = 1.d0
c           DiagScaling(k,1,2) = 0.d0
c           DiagScaling(k,2,1) = 0.d0
c           DiagScaling(k,2,2) = 1.d0
c        enddo


c        do k = 1,NbCV
c           mk = ind_ligne_Aij(k)+1 
c           do m = mk,ind_ligne_Aij(k+1)

c           do i = 1,NbIncPrim
c           do j = 1,NbIncPrim
c              AA(i,j) = Aij(m,i,j) 
c           enddo
c           enddo

c           do i = 1,NbIncPrim
c           do j = 1,NbIncPrim
c              s = 0.d0
c              do l =1,NbIncPrim
c                 s = s + DiagScaling(k,i,l)*AA(l,j)
c              enddo
c              Aij(m,i,j) = s
c           enddo
c           enddo

c           enddo

c           do i =1,NbIncPrim
c              XX(i) = Smm(k,i)
c           enddo

c           do i =1,NbIncPrim
c           s = 0.d0
c           do l =1,NbIncPrim
c              s = s + DiagScaling(k,i,l)*XX(l)
c           enddo
c           Smm(k,i) = s
c           enddo
c        enddo

       

c        deallocate(DiagScaling)


cccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     dimension du vecteur d'inconnues 
c     et nb non zero de la matrice 
c

        
        nsolve = NbCV*NbIncPrim
        mem_syst = MemBloc*NbIncPrim**2 

        IPRECP = 0

c
cccccccccccccccc
c
c
c     espaces memoire pour AMG pour le bloc App
c
        if (IPREC.eq.1) then 
           mem_amg = 20*MemBloc + 5*NbCV 
           mem_IA_amg = 5*NbCV
        else 
           mem_amg = 1 
           mem_IA_amg = 1 
       endif
c
c
c
c     espace memoire mrwork et iwork pour ILU0
c
        mr_ilu0 = 2*mem_syst + 3*nsolve + 1
        mi_ilu0 = 4*mem_syst + 2*nsolve + 8       
c
c
c     espace memoire total pour IWORK et RWORK     
c
c

        miwork = mi_ilu0 + mem_IA_amg*4 + mem_amg + 2*nsolve +1

        mrwork = mr_ilu0 + mem_IA_amg*2 + mem_amg + NbCV 
c
ccccccccccccccccccccc

c        write(*,*)' NbCV ',NbCV
c        write(*,*)' MemBloc ',MemBloc     
c        write(*,*)' mem_syst ',mem_syst
c        write(*,*)' nsolve ',nsolve 
c        write(*,*)' mr_ilu0 ',mr_ilu0
c        write(*,*)' mi_ilu0 ',mi_ilu0
c        write(*,*)' mem_amg ',mem_amg
c        write(*,*)' mem_IA_amg ',mem_IA_amg
c        write(*,*)' miwork ',miwork
c        write(*,*)' mrwork ',mrwork



c
c     allocations 
c
        allocate(RWORK(mrwork))
        allocate(IWORK(miwork))
        allocate(Sol(nsolve))
        allocate(Sm(nsolve))
        allocate(RGWK(1+nsolve*(MAXL+6)+MAXL*(MAXL+3)))

c
c
c
        allocate(Asyst(mem_syst))
        allocate(ncol_Asyst(mem_syst))
        allocate(nligne_Asyst(mem_syst))
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc


        do k=1,NbCV
           do i=1,NbIncPrim
              Sm(k+(i-1)*NbCV) = Smm(k,i)
              Sol(k+(i-1)*NbCV) = Soll(k,i)
           enddo
        enddo


ccccccccccccccccccccccccccccccccccccccccccccccccccccc       
c
c     
c     POINTEURS POUR RWORK ET IWORK
c
c
c
      jdm1 = 1
      IWORK(jdm1) =   IPREC
      jd0 = jdm1 + 1

      IWORK(jd0) =   IPRECP
      jd1 = jd0 + 1

      IWORK(jd1) =   mem_syst
      jd2 = jd1 + 1

      IWORK(jd2)  = NbCV
      jd3 = jd2 + 1        

      IWORK(jd3) = nsolve
      jd4 = jd3 + 1

      IWORK(jd4) = NbIncPrim 
      jd5 = jd4 + 1

c
c
c     ILU0 
c
c      IWORK(jd5): ncol_AsystL(mem_syst)
      jd6 = jd5 + mem_syst

c      IWORK(jd6): nligne_AsystL(mem_syst)
      jd7 = jd6 + mem_syst
    
c      IWORK(jd7): ncol_AsystU(mem_syst)
      jd8 = jd7 + mem_syst

c      IWORK(jd8): nligne_AsystU(mem_syst)
      jd9 = jd8 + mem_syst
c
c
c
c      IWORK(jd9): mem_amg pour App_amg et JA_amg
      IWORK(jd9) =  mem_amg
      jd10 = jd9 + 1

c      IWORK(jd10): mem_IA_amg pour IA_amg
      IWORK(jd10) = mem_IA_amg
      jd11 = jd10 + 1
c
c      IWORK(jd11): IApp_amg(mem_IA_amg)
      jd12 = jd11 + mem_IA_amg
c
c      IWORK(jd12): JApp_amg(mem_amg)
      jd13 = jd12 +  mem_amg
c
c      IWORK(jd13): IG_amg
      jd14 = jd13 + mem_IA_amg*2
c
c      IWORK(jd14): ind_ligne_Asyst(nsolve+1) 
      jd15 = jd14 + nsolve + 1
c
c
cccccccccccccccccc
c
c
c
c
c
        id1 = 1
c
c
c     pour ILU0 COMBINATIVE 
c
c        RWORK(id1): AsystL
        id2 = id1 + mem_syst

c        RWORK(id2): AsystU
        id3 = id2 + mem_syst

c        RWORK(id3): AsystinvDiag
        id4 = id3 + nsolve

c        RWORK(id4): RRR
        id5 = id4 + nsolve

c        RWORK(id5): ZZZ
        id6 = id5 + nsolve
c
c        RWORK(id6): App_amg(mem_amg)
        id7 = id6 + mem_amg

c        RWORK(id7): U_amg(mem_IA_amg)
        id8 = id7 + mem_IA_amg

c        RWORK(id8): F_amg(mem_IA_amg)
        id9 = id8 + mem_IA_amg

c     diagonale de App pour decouplage ATHOS 
c      (supression du scaling)
c
c        RWORK(id9): diagApp
c
        id10 = id9 + NbCV
c
c
c
c
         if (jd15.gt.miwork) then 
            write(*,*)' dim IWORK insuffisante ',jd15,miwork
            stop
         endif
         if (id10.gt.mrwork) then 
            write(*,*)' dim RWORK insuffisante ',id10,mrwork
            stop
         endif           
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Preconditionnement ILU0 du systeme
c
c
c     ILU0 reservoir + jacobi puits 
c     ordonnancement par variables 
c
c

c      write(*,*)' debut ilu '


         allocate(nrow(nsolve))
         allocate(ncol(nsolve))

      call init_prec_ilu0_bloc(
     &       nsolve,NbCV,NbIncPrim,MemBloc,mem_syst,
     &       nligne_Aij,ncol_Aij,ind_ligne_Aij,Aij,
     &       Asyst,nligne_Asyst,ncol_Asyst,
     &    RWORK(id1),IWORK(jd5),IWORK(jd6),RWORK(id3),
     &    RWORK(id2),IWORK(jd7),IWORK(jd8),
     &       nrow,ncol)

      deallocate(nrow)
      deallocate(ncol)

      
c      write(*,*)' fin ilu '
c
c
c     INIT PRECONDITIONNEMENT AMG
c

      if (IPREC.ne.0) then 

      call CPU_TIME(TIME1_initamg) 

      call init_prec_AMG_blocp(NbCV,NbIncPrim,MemBloc,
     &       mem_amg,mem_IA_amg,   
     &       Aij,ind_ligne_Aij,ncol_Aij,
     &       RWORK(id6),IWORK(jd11),IWORK(jd12),
     &       RWORK(id7),RWORK(id8),IWORK(jd13),
     &       RWORK(id9),IPRECP,
     &       newt)

c     &       App_amg,IApp_amg,JApp_amg,U,F,IG,
c     &       diagApp,IPRECP)

      call CPU_TIME(TIME2_initamg) 

      CPU_initamg = TIME2_initamg-TIME1_initamg
     
c      write(*,*)' '
c      write(*,*)' CPU init amg ',CPU_initamg
c      write(*,*)' '

      endif

cccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c     stockage creux du SYSTEME 
c
c     Numerotation par blocs d'inconnues 
c
          nz_syst = mem_syst
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     equations reservoir
c


        l = 0


        do i=1,NbIncPrim
        do k=1,NbCV
        
         ii = k + (i-1)*NbCV
   
         IWORK(jd14+ii-1) = l
      
           do j=1,NbIncPrim
           do m=ind_ligne_Aij(k)+1,ind_ligne_Aij(k+1)

             kv = ncol_Aij(m)

             l = l+1
             nligne_Asyst(l) = ii
             ncol_Asyst(l) = kv + (j-1)*NbCV
             Asyst(l) =  Aij(m,i,j)

           enddo
           enddo


  
         enddo
         enddo


         IWORK( jd14+nsolve+1-1) = l
c
c
c
        if (l.ne.nz_syst) then 
           write(*,*)' pb nz Asyst ',nz_syst,l
           stop
        endif

c     
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     GMRES 
c 
c
c
c
c
c
      ISYM = 0
      ITOL = 0
      TOL = critere_arret_gmres
      
c
c     ne sert pas ITMAX = MAXL*(NRMAX + 1)
c     ie dim Krylov*nb de restart authorise
c
c
      ITMAX = 1000
c
c     numero affichage ecran
c
c     6 si affichage des iterations 0 sinon 
c
c      IUNIT = 6
      IUNIT = 0
c
c     dim du tableau RGWK (workspace)
c
      LRGW = 1+nsolve*(MAXL+6)+MAXL*(MAXL+3)
c
c     dim Krylov space 
c
      IGWK(1) = MAXL
      IGWK(2) = MAXL
c
c     no scaling residu et solution 
c
      IGWK(3) = 0
c
c     IGWK(4) => preconditionnement a droite (>0)
c     a gauche (<0)
c
      IGWK(4) = 1

      if (IPREC.eq.0) then 
c
c     nb max de restart
c
      IGWK(5) = 40

      else 

      IGWK(5) = 5

c
      endif
c
c     dim du tableau IGWK
c
c     
c 
      LIGW = 20
c
c   
c
cccccccccccccc

c      write(*,*)' debut GMRES ' 


c      do i=1,nsolve
c         Sol(i) = 0.d0
c      enddo

c      write(*,*)' nsolve ',nsolve
c      write(*,*)' nz_syst ',nz_syst
c      do i=1,nsolve
c         write(*,*) ' Sm ',Sm(i)
c      enddo
c      write(*,*)
c      do i=1,nsolve
c         write(*,*) ' Sol ',Sol(i)
c      enddo
c      write(*,*)
c      do i=1,nz_syst
c         write(*,*)' i  j Aij ',nligne_Asyst(i),ncol_Asyst(i),Asyst(i)
c      enddo
c      write(*,*)
     

      allocate(SB(nsolve))
      allocate(SX(nsolve))

      CALL DGMRES (nsolve, Sm, Sol, nz_syst, nligne_Asyst, ncol_Asyst, 
     &   Asyst, ISYM, MATVEC, MSOLVE,
     &   ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, SB, SX, RGWK, LRGW,
     &   IGWK, LIGW, RWORK, IWORK) 

      deallocate(SB)
      deallocate(SX)
       
      if (IERR.ne.0) then 
       write(*,*)' ERROR ESTIMATE (residu relatif prec)',ERR
       write(*,*)' FLAG ERROR ',IERR      
       write(*,*)' RHOL ',RGWK(1)     
       write(*,*)' IGWK(6) ',IGWK(6)
      endif

        do k=1,NbCV
           do i=1,NbIncPrim
             Soll(k,i) = Sol(k+(i-1)*NbCV)
           enddo
        enddo

cccccccccccccccccccccccc
c
c     on enleve le scaling des inconnues 
c
c        do k=1,NbCV
c           S1 = Soll(k,1) 
c           S2 = Soll(k,2) 
c           Soll(k,1) = S2
c           Soll(k,2) = S1 
c        enddo
c
ccccccccccccccccccccccccccccccccccccccc
c
        deallocate(RWORK)
        deallocate(IWORK)
        deallocate(Asyst)
        deallocate(ncol_Asyst)
        deallocate(nligne_Asyst)
        deallocate(Sol)
        deallocate(Sm)

ccccccccccccccccccccccccccccccccccccccc

        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     produit matrice vecteur: Y = A*X
c
c     stockage creux de A: A_IA(l),JA(l) = A(l) 
c
c     NELT: nb d'elts non nuls de A 
c
c     N: dimension de la matrice A 
c
c
             subroutine MATVEC(N, X, Y, NELT, IA, JA, A, ISYM)
c
c
c
      DOUBLE PRECISION A(NELT), X(N), Y(N)             
      INTEGER N, NELT, IA(NELT),  JA(NELT), I, ISYM
c
c
      if (ISYM.ne.0) then 
         write(*,*)' on doit avoir ISYM=0 (mat non sym) ',ISYM
         stop
      endif
c
c     ---------------------------------
c
      DO I=1,N
         Y(I) = 0.D0
      ENDDO

      DO I=1,NELT
         Y(IA(I)) = Y(IA(I)) + A(I)*X(JA(I))
      ENDDO


      RETURN
      END

cccccccccccccccccccccccccccccccccccccccccc
c
c
c     Preconditionnement: Z = M^{-1} R 
c
c      pour R donne
c
c
      subroutine MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
c
c
	implicit double precision (a-h,o-z)
	implicit integer (i-n)
c
c
      DOUBLE PRECISION RWORK(*), R(N), Z(N), A(NELT)
      INTEGER IWORK(*),IA(NELT),JA(NELT)
      
      if (ISYM.ne.0) then 
         write(*,*)' on doit avoir ISYM=0 (mat non sym) ',ISYM
         stop
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccc
cccccccccc
c     
c     POINTEURS POUR RWORK ET IWORK
c
c
c

      jdm1 = 1
      IPREC = IWORK(jdm1)  
      jd0 = jdm1 + 1

      IPRECP = IWORK(jd0)  
      jd1 = jd0 + 1
      
      mem_syst = IWORK(jd1)  
      jd2 = jd1 + 1

      ncell = IWORK(jd2) 
      jd3 = jd2 + 1        

      nsolve = IWORK(jd3) 
      jd4 = jd3 + 1

      ninc = IWORK(jd4) 
      jd5 = jd4 + 1

c
c
c     ILU0 
c
c      IWORK(jd5): ncol_AsystL(mem_syst)
      jd6 = jd5 + mem_syst

c      IWORK(jd6): nligne_AsystL(mem_syst)
      jd7 = jd6 + mem_syst
    
c      IWORK(jd7): ncol_AsystU(mem_syst)
      jd8 = jd7 + mem_syst

c      IWORK(jd8): nligne_AsystU(mem_syst)
      jd9 = jd8 + mem_syst
c
c
c
c
c
c      IWORK(jd9): mem_amg pour App_amg et JA_amg
      mem_amg = IWORK(jd9)
      jd10 = jd9 + 1

c      IWORK(jd10): mem_IA_amg pour IA_amg
      mem_IA_amg = IWORK(jd10)
      jd11 = jd10 + 1
c
c      IWORK(jd11): IApp_amg(mem_IA_amg)
      jd12 = jd11 + mem_IA_amg
c
c      IWORK(jd12): JApp_amg(mem_amg)
      jd13 = jd12 +  mem_amg
c
c      IWORK(jd13): IG_amg
      jd14 = jd13 + mem_IA_amg*2
c
c      IWORK(jd14): ind_ligne_Asyst(nsolve+1) 
      jd15 = jd14 + nsolve + 1
c
c
c
cccccccccccccccccc
c
c
c
c
c
        id1 = 1
c
c
c     pour ILU0 COMBINATIVE 
c
c        RWORK(id1): AsystL
        id2 = id1 + mem_syst

c        RWORK(id2): AsystU
        id3 = id2 + mem_syst

c        RWORK(id3): AsystinvDiag
        id4 = id3 + nsolve

c        RWORK(id4): RRR
        id5 = id4 + nsolve

c        RWORK(id5): ZZZ
        id6 = id5 + nsolve
c
c
c
c        RWORK(id6): App_amg(mem_amg)
        id7 = id6 + mem_amg

c        RWORK(id7): U_amg(mem_IA_amg)
        id8 = id7 + mem_IA_amg

c        RWORK(id8): F_amg(mem_IA_amg)
        id9 = id8 + mem_IA_amg

c     diagonale de App pour decouplage ATHOS 
c      (supression du scaling)
c
c        RWORK(id9): diagApp
c
        id10 = id9 + ncell 
c
c
c
ccccccccccccccccccccccccccccccccccccccccc



ccccccccccccccccccccccccccccccccccccccccc

c
c     PRECONDITIONNEMENT 
c
c     ILU0 sur le systeme  
c
c
c     Precond ILU0 avec inconnues  
c     ordonnee par blocs mailles
c
c
ccccccccccccccccccccccccccccccccccccccccc
c
c     preconditionnement ILU0 du systeme Z = ILU0(A) R
c
      call pilu0_syst_bloc(nsolve,ncell,ninc,mem_syst,
     &    RWORK(id1),IWORK(jd5),IWORK(jd6),RWORK(id3),
     &    RWORK(id2),IWORK(jd7),IWORK(jd8),
     &    R,Z)
c
c     si pas de preconditionnement
c
c        do i=1,nsolve
c           Z(i) = R(i)
c        enddo
c
c     si AMG sur la pression ou pression
c
c     ou si IPREC=0 --> ILU0 seul 
c

c      write(*,*)' IPREC ',IPREC

      if (IPREC.ne.0) then 


c
c     calcul du residu des equations de pression
c
c     RWORK(id5) = (1,h=0,l=0,puits=0)(R-AZ)
c

      call residu_pression_p(
     &               nsolve,ncell,mem_syst,ninc,
     &               IWORK(jd14),IA,JA,A,
     &               R,Z,RWORK(id5))


      if (IPRECP.eq.1) then 
         do k=1,ncell
            RWORK(id5+k-1) = RWORK(id5+k-1)*RWORK(id9 + k -1)
         enddo
      endif
c
c
c     preconditionnement AMG de la pression
c
c     RWORK(id4) = AMG RWORK(id5) 
c

      call  AMG_App(ncell,mem_amg,mem_IA_amg,
     &       RWORK(id6),IWORK(jd11),IWORK(jd12),
     &       RWORK(id7),RWORK(id8),IWORK(jd13),
     &       RWORK(id5),RWORK(id4))

c
c
c     correction de la pression reservoir
c
      do k=1,ncell
         Z(k) = Z(k) + RWORK(id4+k-1)
      enddo
c
c
c     fin if AMG 
c
      endif
c
ccccccccccccccccccccccccccccccccccccccccc
      RETURN
      END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     En entree: systemes par blocs Aij et puits 
c
c
c
c     En sortie: Factorisation 
c     du preconditionnement ILU0 du systeme 
c
c     Stockage du systeme Asyst par blocs d'inconnues  
c
c
c     PUITS: ON MET LA DIAGONALE SUR LA PARTIE PUITS 
c           
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine init_prec_ilu0_bloc(
     &       nsolve,ncell,ninc,mem_bloc,mem_syst,
     &       nligne_Aij,ncol_Aij,ind_ligne_Aij,Aij,
     &       Asyst,nligne_Asyst,ncol_Asyst,
     &       AsystL,ncol_AsystL,nligne_AsystL,AsystinvDiag,
     &       AsystU,ncol_AsystU,nligne_AsystU,
     &       nrow,ncol)
c
c
ccccccccccccccccccccccccccccccccccccccccc
c
	implicit double precision (a-h,o-z)
	implicit integer (i-n)
c
c
c
c
        dimension Aij(mem_bloc,ninc,ninc)
        dimension ncol_Aij(mem_bloc)
        dimension nligne_Aij(mem_bloc)
        dimension ind_ligne_Aij(ncell+1)
c
c
c     FACTORISATION ILU0 de Asyst 
c
c     Asyst 
c
        dimension Asyst(mem_syst)
        integer ncol_Asyst(mem_syst),nligne_Asyst(mem_syst)
c
c     factorisation AsystL et AsystU et AsystinvDiag
c
        integer ncol_AsystL(mem_syst)
        integer nligne_AsystL(mem_syst)
        integer ncol_AsystU(mem_syst)
        integer nligne_AsystU(mem_syst)

        dimension AsystL(mem_syst)
        dimension AsystU(mem_syst)

        dimension AsystinvDiag(nsolve)
c
c     workspace pour fact ILU0
c
        integer nrow(nsolve), ncol(nsolve)
c
c       
        integer NL,NU,ISYM
c
c
c
c
c
c
c
c     stockage creux du SYSTEME 
c     Numerotation par inconnues pour chaque maille 
c     pour reduire la largeur de bande 
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        nsolve_res = ncell*ninc


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     stockage CSC pour ILU0 
c
c
c
c     comptage du nb d'elt par colonne ncol
c
        do jj=1,nsolve
           ncol(jj) = 0
        enddo


c
c     diagonale first
c      
        do i=1,ninc
        do k=1,ncell
        
         ii = k + (i-1)*ncell   
      
           j=i
           m=ind_ligne_Aij(k)+1

             kv = ncol_Aij(m)
             jj = kv + (j-1)*ncell
             ncol(jj) = ncol(jj)+1
        enddo
        enddo



        do i=1,ninc
        do k=1,ncell
        
         ii = k + (i-1)*ncell   
      
           j=i
           do m=ind_ligne_Aij(k)+2,ind_ligne_Aij(k+1)

             kv = ncol_Aij(m)

             jj = kv + (j-1)*ncell
             ncol(jj) = ncol(jj)+1

           enddo

           do j=1,ninc
            if (j.ne.i) then 
               do m=ind_ligne_Aij(k)+1,ind_ligne_Aij(k+1)

                  kv = ncol_Aij(m)                  
                  jj = kv + (j-1)*ncell
                  ncol(jj) = ncol(jj)+1
                  
               enddo
            endif
           enddo
  
         enddo
         enddo


cccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     remplissage de ncol_Asyst(1,...,nsolve+1): 
c
c     colonne jj Asyst de ncol_Asyst(jj),...,ncol_Asyst(jj+1)-1
c
c
         ncol_Asyst(1) = 1
         do jj=1,nsolve
            ncol_Asyst(jj+1) = ncol_Asyst(jj) + ncol(jj)
         enddo
c
c
c     remplissage de Asyst
c
        do jj=1,nsolve
           ncol(jj) = 0
        enddo


        do i=1,ninc
        do k=1,ncell
        
         ii = k + (i-1)*ncell   
      
           j=i
           m=ind_ligne_Aij(k)+1

             kv = ncol_Aij(m)
             jj = kv + (j-1)*ncell
             ncol(jj) = ncol(jj)+1
             l = ncol_Asyst(jj)+ncol(jj)-1
             Asyst(l) =  Aij(m,i,j)
             nligne_Asyst(l) = ii
        enddo
        enddo

        do i=1,ninc
        do k=1,ncell
        
         ii = k + (i-1)*ncell   
      
           j=i
           do m=ind_ligne_Aij(k)+2,ind_ligne_Aij(k+1)

             kv = ncol_Aij(m)

             jj = kv + (j-1)*ncell
             ncol(jj) = ncol(jj)+1
             l = ncol_Asyst(jj)+ncol(jj)-1
             Asyst(l) =  Aij(m,i,j)
             nligne_Asyst(l) = ii
           enddo

           do j=1,ninc
            if (j.ne.i) then 
               do m=ind_ligne_Aij(k)+1,ind_ligne_Aij(k+1)

                  kv = ncol_Aij(m)                  
                  jj = kv + (j-1)*ncell
                  ncol(jj) = ncol(jj)+1
                  l = ncol_Asyst(jj)+ncol(jj)-1
                  Asyst(l) =  Aij(m,i,j)
                  nligne_Asyst(l) = ii                 
               enddo
            endif
           enddo
  
         enddo
         enddo

c
c
c
c
c
         mem_ilu = ncol_Asyst(nsolve+1)-1

c         write(*,*)' mem_ilu ',mem_ilu
c         write(*,*)' mem_syst ',mem_syst
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     NB: mem_syst est la dimension du tableau passe 
c     a DSILUS, pas le nb d'elts non nuls
c
c
      ISYM = 0
      call DSILUS (nsolve,mem_syst,nligne_Asyst,ncol_Asyst,Asyst, 
     &   ISYM,NL,nligne_AsystL,ncol_AsystL,AsystL,AsystinvDiag,
     &   NU,nligne_AsystU,ncol_AsystU,AsystU,nrow,ncol)

c      call CPU_TIME(TIME2_init) 

c      write(*,*)' ilu0 ',TIME2_init-TIME1_init


c      write(*,*)' NL NU ',NL,NU,ISYM

c      do i=1,NL
c         write(*,*)' L ',i,nligne_AsystL(i),ncol_AsystL(i),AsystL(i)
c      enddo
 
c 
c      write(*,*)' fin init ILU0 ',mem_syst,nsolve,ncell,ninc,mem_bloc
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      subroutine pilu0_syst_bloc(nsolve,ncell,ninc,mem_syst,
     &       AsystL,ncol_AsystL,nligne_AsystL,AsystinvDiag,
     &       AsystU,ncol_AsystU,nligne_AsystU,
     &       R,Z)
c
c
ccccccccccccccccccccccccccccccccccccccccc
c
	implicit double precision (a-h,o-z)
	implicit integer (i-n)
c
c
        dimension R(nsolve), Z(nsolve)
c
c
c
c     FACTORISATION ILU0 de Asyst 
c
c     factorisation AsystL et AsystU
c
        integer ncol_AsystL(mem_syst)
        integer nligne_AsystL(mem_syst)
        integer ncol_AsystU(mem_syst)
        integer nligne_AsystU(mem_syst)

        dimension AsystL(mem_syst)
        dimension AsystU(mem_syst)

        dimension AsystinvDiag(nsolve)
c
c
cccccccccccccccccccccccccccccccccccccccccccc

c
c     Z = ILU0(Asyst)^{-1} (R)
c
      call DSLUI2 (nsolve, R, Z, 
     &    nligne_AsystL,ncol_AsystL,AsystL,AsystinvDiag,
     &    nligne_AsystU,ncol_AsystU,AsystU)
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      return
      end
c
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     INITIALISATION DU PRECONDITIONNEUR AMG (App)
c  
c     * SETUP AMG App
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine init_prec_AMG_blocp(ncell,ninc,mem_bloc,
     &       mem_amg,mem_IA_amg,   
     &       Aij,ind_ligne_Aij,ncol_Aij,
     &       App_amg,IApp_amg,JApp_amg,U,F,IG,
     &       diagApp,IPRECP,
     &       newt)

	implicit double precision (a-h,o-z)
	implicit integer (i-n)

c
c     en entree: 
c
c
c     blocs Aij
c
        dimension Aij(mem_bloc,ninc,ninc)

        integer ncol_Aij(mem_bloc)
        integer ind_ligne_Aij(ncell+1)
c
c     sorties: setup AMG: App_amg IApp_amg JApp_amg
c      

c
c
        integer IApp_amg(mem_IA_amg)
        integer JApp_amg(mem_amg)        
        dimension App_amg(mem_amg)
c
        integer IG(2*mem_IA_amg)
c
c
c
c     sol init (zero) et finale VCYCLEAMG 
c       (ne sert pas: setup AMG slt)
c
        dimension U(mem_IA_amg)
c
c
c
c     RHS (ne sert pas car setup AMG slt)
c
        dimension F(mem_IA_amg)
c
c
c
c     Diagonale du bloc App ds le cas du decouplage 
c     ATHOS pour supression du scaling slt pour 
c     IPRECP = 3
c
c
        dimension diagApp(ncell)


cccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
        if (IPRECP.eq.1) then 
         do k=1,ncell
         do l=ind_ligne_Aij(k) + 1,ind_ligne_Aij(k+1)
           App_amg(l) = Aij(l,1,1)*diagApp(k)
           JApp_amg(l) = ncol_Aij(l)
         enddo
         enddo

        else 

         do k=1,ncell
         do l=ind_ligne_Aij(k) + 1,ind_ligne_Aij(k+1)
           App_amg(l) = Aij(l,1,1)
           JApp_amg(l) = ncol_Aij(l)
         enddo
         enddo

        endif

        do k=1,ncell+1
           IApp_amg(k) = ind_ligne_Aij(k) + 1
        enddo

c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccc
c
c        do k=1,ncell
c           write(*,*)' ligne ',k
c           do m=IApp_amg(k),IApp_amg(k+1)-1
c             write(*,*)m,JApp_amg(m),App_amg(m)
c           enddo  
c           write(*,*)' '         
c        enddo
        


c
c
c     test de diagonale dominance des eqs reservoir 
c
          nl = 0
          nlf = 0
          do k=1,ncell
             m = IApp_amg(k)
             d = App_amg(m)
c             write(*,*)' diag App ',k,d
             s = 0.d0
             do m=IApp_amg(k)+1,IApp_amg(k+1)-1
                s = s + dabs(App_amg(m))
             enddo

             if (s.gt.d*(1.+1.0E-8)) then
c                write(30,*)newt,(s-d)/d
                nl = nl + 1
             endif

             if (s.gt.d*1.5) then
                nlf = nlf + 1
             endif
c
c     on corrige la non diag dom
c
             s = 0.d0
             do m=IApp_amg(k)+1,IApp_amg(k+1)-1
                s = s - App_amg(m)
             enddo
             if (s.gt.d) then 
c                App_amg(IApp_amg(k)) = s
             endif
c
c     fin correction
c
          enddo
          nbligne_nondiagdom = nl
          nbligne_fort_nondiagdom = nlf

c
c     test de negativite des termes hors diagonaux 
c           eqs reservoir 
c
          nl = 0
          nlf = 0
          do k=1,ncell
             md = IApp_amg(k)
             d = App_amg(md)
             amax = 0.0
             do m=IApp_amg(k)+1,IApp_amg(k+1)-1
                amax = dmax1(amax,dabs(App_amg(m)))
             enddo

             do m=IApp_amg(k)+1,IApp_amg(k+1)-1
               if (App_amg(m).gt.0.0) then
c                 write(31,*)newt,App_amg(m)/d
                 nl = nl + 1
               endif
               if (App_amg(m).gt.0.5*amax) then
                nlf = nlf + 1
               endif
             enddo
          enddo
          nb_horsdiag_positif = nl
          nb_horsdiag_fort_positif = nlf

c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     test de positivite de la diagonale
c
          nl = 0
          do k=1,ncell
             m = IApp_amg(k)
             d = App_amg(m)
             if (d.lt.0.0) then
                nl = nl + 1

c                write(*,*)' diag App <0 ',k,d

c                stop
             endif
          enddo
          nb_diag_negatif = nl


c          write(*,*)' nb ligne non diag dom (et fort) bloc P = ',
c     &            nbligne_nondiagdom,
c     &            nbligne_fort_nondiagdom

          


c
c
c
c          write(*,*)' nb hors diag positif (et fort) bloc P = ',
c     &                             nb_horsdiag_positif,
c     &                             nb_horsdiag_fort_positif



c          write(*,*)' nb diag negatif bloc P = ',nb_diag_negatif
c
c
c
c
c     Correction du bloc pression pour supprimer les 
c     termes hors diagonaux positifs : 
c
c     on les condense sur la diagonale et on les annule 
c                         OU 
c     on les annule sans condenser (choix A FAIRE!)
c
c
c
c
      do k=1,ncell
         md = IApp_amg(k)

c         write(*,*)' diagonale App ',App_amg(md)

c         amax = 0.0
c         do m=ind_ligne_Aij(k)+2,ind_ligne_Aij(k+1)
c            amax = dmax1(amax,dabs(App_amg(m)))
c         enddo

         do m=IApp_amg(k)+1,IApp_amg(k+1)-1
c            if (App_amg(m).gt.0.5*amax) then
            if (App_amg(m).gt.0.0) then
c             App_amg(md) = App_amg(md) + App_amg(m)
c             App_amg(m) = 0.d0
            endif
         enddo
      enddo

c
c     on corrige la non diagonale dominance des eqs de reservoir
c
      do k=1,ncell
             m = IApp_amg(k)
             d = App_amg(m)
             s = 0.d0
             do m=IApp_amg(k)+1,IApp_amg(k+1)-1
                s = s - App_amg(m)
             enddo            
             if (s.gt.d) then
               m = IApp_amg(k)
c               App_amg(m) = App_amg(m) + (s-d)
             endif
      enddo
c
c
c
cccccccccccccccccccccccccccccccccc
c
c
c     parametre AMG: setup
c
c
          ISWTCH = 4
c          IOUT   = 12 default 11 used
          IOUT = 10
          IPRINT = 10606 
c          IPRINT = -1
C
          LEVELX = 20
c          IFIRST = 13 default
          IFIRST = 10
c          NCYC   = 10110 default
          NCYC   = 1010
c
c     NO CYCLING (4th digit=0) = SETUP
c
          EPS    = 1.D-12
          MADAPT = 27
          NRD    = 1131
          NSOLCO = 110
          NRU    = 1131
C
C
c     ECG1: si je comprends biens definit 
c       les pts isoles si sum |off diag| <= ECG1*diag
c
          ECG1   = 0.
c
c     ECG2: define strong connections (alpha)
c
c     default 0.25
c          
          ECG2   = 0.25

c
c     define strong dependence on a A set (beta)
c
c     default = 0.35
c
          EWT2   = 0.35
c
c 
          NWT    = 2
          NTR    = 0
c
c
          MATRIX = 22

          
          NNU = ncell
          

          NDA = mem_amg
          NDJA = mem_amg
          NDIA = mem_IA_amg
          NDU = mem_IA_amg
          NDF = mem_IA_amg
          NDIG = 2*mem_IA_amg

c        write(*,*)' DEBUT INIT AMG '


c      write(*,*)NDA,NDIA,NDJA,NDU,NDF,NDIG,NNU,MATRIX
c      write(*,*)ISWTCH,IOUT,IPRINT
c      write(*,*)LEVELX,IFIRST,NCYC,EPS,MADAPT,NRD,NSOLCO,NRU
c      write(*,*)ECG1,ECG2,EWT2,NWT,NTR


      call AMG1R5(App_amg,IApp_amg,JApp_amg,U,F,IG,
     +                  NDA,NDIA,NDJA,NDU,NDF,NDIG,NNU,MATRIX,
     +                  ISWTCH,IOUT,IPRINT,
     +                  LEVELX,IFIRST,NCYC,EPS,MADAPT,NRD,NSOLCO,NRU,
     +                  ECG1,ECG2,EWT2,NWT,NTR,
     +                  IERR)


c      write(*,*)' fin init prec amg1r5 '

c        if (newt.eq.337) then 
c           stop
c        endif
c
c
c      write(*,*)' IERR ',IERR
c
c
c
cccccccccccccccccccccccccccccccccc
c
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     PRECONDTIONNEMENT DU BLOC App par AMG1R5 (1 VCYCLE)
c
c     EN ENTREE: matrice App_amg, IApp_amg, IApp_amg 
c         issu de la setup phase de AMG1R5 
c
c     EN ENTREE smembre R(1,...,ncell)
c
c     EN SORTIE Z(1,...,ncell) = VCYCLEAMG R
c
c
c
      subroutine  AMG_App(ncell,mem_amg,mem_IA_amg,
     &       App_amg,IApp_amg,JApp_amg,U,F,IG,
     &       R,Z)

	implicit double precision (a-h,o-z)
	implicit integer (i-n)

c
c     AMG
c
c

c
c
        integer IApp_amg(mem_IA_amg)
        integer JApp_amg(mem_amg)        
        dimension App_amg(mem_amg)
c
        integer IG(2*mem_IA_amg)
c
c     sol init (zero) et finale VCYCLEAMG
c
        dimension U(mem_IA_amg)
c
c     RHS
c
        dimension F(mem_IA_amg)
c
c
c
        dimension R(ncell), Z(ncell)
c
cccccccccccccccccccccc
        do k=1,ncell
           U(k) = 0.d0
           F(k) = R(k)
        enddo

c
c     parametre AMG: solve (setup deja effectue)
c
c
          ISWTCH = 2
c          IOUT   = 12 default
c          IOUT = 11: affiche les Vcycles, 10: affiche rien
          IOUT = 10
c          IPRINT = 10606 default
          IPRINT = 10606
C
          LEVELX = 20
c          IFIRST = 13 default
          IFIRST = 10
c          NCYC   = 10110 default
c
c     ici 1 VCYCLE
c
          NCYC   = 1011
c
c     ici 5 seul VCYCLE
c
c          NCYC   = 1015
c
          EPS    = 1.D-12
          MADAPT = 27
          NRD    = 1131
          NSOLCO = 110
          NRU    = 1131
C
c     ECG1: si je comprends biens definit 
c       les pts isoles si sum |off diag| <= ECG1*diag
c
          ECG1   = 0.
c
c     ECG2: define strong connections (alpha) 
c            ij strong iff |Aij| > alpha Sup_{l.ne.j} |Ail| 
c 
c     default = 0.25
c         
          ECG2   = 0.25
c    
c
c     define strong dependence on a A set (beta)
c
c     default = 0.35
c
          EWT2   = 0.35

c
          NWT    = 2
          NTR    = 0
c
c
          MATRIX = 22

          NNU = ncell
          NDA = mem_amg
          NDJA = mem_amg
          NDIA = mem_IA_amg
          NDU = mem_IA_amg
          NDF = mem_IA_amg
          NDIG = 2*mem_IA_amg


      call AMG1R5(App_amg,IApp_amg,JApp_amg,U,F,IG,
     +                  NDA,NDIA,NDJA,NDU,NDF,NDIG,NNU,MATRIX,
     +                  ISWTCH,IOUT,IPRINT,
     +                  LEVELX,IFIRST,NCYC,EPS,MADAPT,NRD,NSOLCO,NRU,
     +                  ECG1,ECG2,EWT2,NWT,NTR,
     +                  IERR)
c
c
      if (IERR.le.0) then 
c         write(*,*)' IERR ',IERR
        do k=1,ncell
           Z(k) = U(k)
        enddo
      else
         write(*,*)' IERR ',IERR
         stop
      endif

c
ccccccccccccccccccccccc
c
        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     CALCUL DU RESIDU DE LA PRESSION
c     RP = (1,h=0,l=0,puits=0)(R-AZ)
c
c
c
c     stockage creux de A: A_IA(l),JA(l) = A(l) 
c
c     NELT: nb d'elts non nuls de A 
c
c     N: dimension de la matrice A 
c
c
             subroutine residu_pression_p(
     &               nsolve,ncell,mem_syst,ninc,
     &               INDLA,IA,JA,A,
     &               R,Z,RP)

c
	implicit double precision (a-h,o-z)
	implicit integer (i-n)
c
        dimension IA(mem_syst), JA(mem_syst), A(mem_syst)
        dimension INDLA(nsolve+1)
        dimension R(nsolve), Z(nsolve)
        dimension RP(ncell)
c
c     ---------------------------------
c

        do k=1,ncell
           RP(k) = R(k)
           do m = INDLA(k)+1,INDLA(k+1)
              RP(k) = RP(k) - A(m)*Z(JA(m))
           enddo
        enddo


      RETURN
      END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     CALCUL DU RESIDU Du SYSTEME
c     Res = (R-AZ)
c
c
c
c     stockage creux de A: A_IA(l),JA(l) = A(l) 
c
c     NELT: nb d'elts non nuls de A 
c
c     N: dimension de la matrice A 
c
c
             subroutine residu_syst(
     &               nsolve,mem_syst,
     &               IA,JA,A,
     &               R,Z,Res)

c
	implicit double precision (a-h,o-z)
	implicit integer (i-n)
c
        dimension IA(mem_syst), JA(mem_syst), A(mem_syst)
        dimension R(nsolve), Z(nsolve)
        dimension Res(nsolve)
c
c     ---------------------------------
c
      DO I=1,nsolve
         Res(I) = R(I)
      ENDDO

      DO m=1,mem_syst
         Res(IA(m)) = Res(IA(m)) - A(m)*Z(JA(m))
      ENDDO


      RETURN
      END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccc

      FUNCTION D1MACH(I)

      DOUBLE PRECISION D1MACH
      INTEGER I

      D1MACH = 1.0E-14

      RETURN
      END
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccc




