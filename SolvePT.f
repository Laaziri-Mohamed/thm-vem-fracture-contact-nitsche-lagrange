ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Calcul de la pression 
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
      subroutine SolvePT(
     &     Deltat, 
     &     NbCell,NbFaceFrac,NbFace,NbArete,NbCV,
     &     NbdofMecaContact, 
     &     VolCell,VolFaceFrac,SurfaceFace,XInc,     
     &     NumFacebyCell,NumCellbyFace,      
     &     NumFaceVersFaceFrac,IndFaceFrac,
     &     NbAretebyFace,NumAretebyFace,      
     &     NumIncCell,NumIncFaceFrac, 
     &     IndIncDir,IndIncNeu,IndIncDirT,IndIncNeuT,NumIncFace,  
     &     NbInterfacebyCell,NumInterfacebyCell,
     &     NbInterfacebyFaceFrac,NumInterfacebyFaceFrac,
     &     NumIncCellbyFaceFrac,NumFaceFracVersFace,             
     &     TransCell,TransFaceFrac,TransCellbyFaceFrac,
     &     TransCellFourier,TransFaceFracFourier,
     &     TransCellbyFaceFracFourier,            
     &     dfcontact,relaxff,relaxmm,      
     &     PorobyCell,DerPorobyCell,PorobyCell_nm1,
     &     dfbyFaceFrac,dfbyFaceFrac_nm1,
     &     AccP,AccP_nm1,
     &     AccT,AccT_nm1,  
     &     SolU_nm1,SolU_prev,SolU0,
     &     IndContact_prev,Sautt_Proj_nm1,       
     &     SolP,SolP_nm1,SolP_prev,SolPDir,SolP0,
     &     SolT,SolT_nm1,SolT_prev,SolTDir,SolT0,           
     &     NbIncGlobalbyCell,NumIncGlobalbyCell,GradCellVEM,     
     &     NbIncGlobalbyFaceFrac,NumIncGlobalbyFaceFrac,SautbyFaceFrac,
     &     VecNormalbyFace,VecTanFace1,VecTanFace2)
c
c
	implicit double precision (a-h,o-z)
	implicit integer (i-n)
c
        include 'include_parameter'
c
c
ccccccccccccccccccccccccccccccccccccccc
c
c     ENTREES
c
c
        dimension NumFacebyCell(NbCell,NbFacebyCellMax)
        dimension NumCellbyFace(NbFace,2)

        dimension NumFaceVersFaceFrac(NbFace)
        dimension IndFaceFrac(NbFace)

        dimension NbAretebyFace(NbFace)
        dimension NumAretebyFace(NbFace,NbAretebyFaceMax)        
c
c
c     coordonnee des inconnues
c
        dimension XInc(NbCV,NbDim) 
c
c     Volume des mailles 
c     num cell  >  VolCell 
c
        dimension VolCell(NbCell)
c
c     Volume des faces frac (surface face * epaisseur de la face fracture)
c     num face frac  >  VolFaceFrac 
c
        dimension VolFaceFrac(NbFaceFrac)
c
c     Surface by Face 
c
        dimension SurfaceFace(NbFace)
c        
c
c     Numero des inconnues mailles et faces fractures 
c
c     num cell       >  num inc cell
c     num face frac  >  num inc face frac  
c
        dimension NumIncCell(NbCell)
        dimension NumIncFaceFrac(NbFaceFrac)

        
c
c     pour HFV (toutes les faces y compris les faces frac)
c     num face  > num inc 
c
        dimension NumIncFace(NbFace)
c        
c     pour HFV discontinu 
c     num face frac  > num inc Ksigma et Lsigma 
c        
        dimension NumIncCellbyFaceFrac(NbFaceFrac,2)              
c        
c
c     Indice inc de bord Dirichlet ou non 
c
c     O si Int ou Neumann, 1 si Dir 
c
c     num inc > indice 0 ou 1 
c
        dimension IndIncDir(NbCV)
        dimension IndIncDirT(NbCV)        
c
c     num inc > indice 0 ou 1 (neu homog) ou 2 (neu non homog) 
c
        dimension IndIncNeu(NbCV)
        dimension IndIncNeuT(NbCV)            
c
c
c     Numero des inconnues d'interfaces IK par maille et Isig par face fracture
c
c     num cell x  num local interface cell             > num inc interface 
c     num face frac x  num local interface face frac   > num inc interface 
c
        dimension NbInterfacebyCell(NbCell)
        dimension NumInterfacebyCell(NbCell,NbInterfacebyCellMax)

        dimension NbInterfacebyFaceFrac(NbFaceFrac)
        dimension NumInterfacebyFaceFrac(NbFaceFrac,
     &                              NbInterfacebyFaceFracMax)

c
c     Transmissivites TK et Tsig 
c     num cell x num local interface cell x num local interface cell                 > TKij 
c     num face frac x num local interface face frac x num local interface face frac  > Tsigij 
c
        dimension TransCell(NbCell,
     &        NbInterfacebyCellMax,NbInterfacebyCellMax)

        dimension TransFaceFrac(NbFaceFrac,
     &       NbInterfacebyFaceFracMax,NbInterfacebyFaceFracMax)


        dimension TransCellFourier(NbCell,
     &        NbInterfacebyCellMax,NbInterfacebyCellMax)

        dimension TransFaceFracFourier(NbFaceFrac,
     &       NbInterfacebyFaceFracMax,NbInterfacebyFaceFracMax)
        
        

c
c     Pour HFV discontinu: TransCellbyFaceFrac : TKsigma et TLsigma 
c
c     num des cells K et L donnes par NumCellbyFace(numface,1:2) 
c
c        
        dimension TransCellbyFaceFrac(NbFaceFrac,2)
        
        dimension TransCellbyFaceFracFourier(NbFaceFrac,2)            
c        
c
c       num face frac > num face 
        dimension NumFaceFracVersFace(NbFaceFrac)        
c
c     en entree 
c
c     CL Dirichlet aux inconnues interfaces de bord DIR 
c
        dimension SolPDir(NbCV)

        dimension SolTDir(NbCV)
c        
c
c     porosite courante (non convergee !!)  
c     donc a recalculer par appel a couplage
c       
       dimension PorobyCell(NbCell)

c     porosite nm1 et der de la porosite courante / P
c       
      dimension PorobyCell_nm1(NbCell)
      dimension DerPorobyCell(NbCell,2)

c
c     aperture au contact par face frac 
c      
      dimension dfcontact(NbFaceFrac)      
c
c     aperture courante et nm1 
c      
      dimension dfbyFaceFrac(NbFaceFrac)
      dimension dfbyFaceFrac_nm1(NbFaceFrac)
c
c     Acc Darcy nm1 
c      
      dimension AccP_nm1(NbCV)

c
c     Acc Fourier nm1 
c      
      dimension AccT_nm1(NbCV)
c
      
c
c
c     U au temps initial, nm1 et a l'iteration prev (previous) 
c
      dimension SolU0(NbdofMecaContact,NbDim)
      dimension SolU_nm1(NbdofMecaContact,NbDim)
      dimension SolU_prev(NbdofMecaContact,NbDim)

c
c     Indice contact 
c
      dimension IndContact_prev(NbFaceFrac)
c
c     Saut tangentiel moyen projete a nm1 
c
      dimension Sautt_Proj_nm1(NbFaceFrac,NbDim-1)      
c
c     Pression au temps nm1 et a l'iteration prev du point fixe 
c
      dimension SolP_nm1(NbCV)
      dimension SolP_prev(NbCV)
c
c     Pression a t=0 
c      
      dimension SolP0(NbCV)


c
c    Temperature au temps nm1 et a l'iteration prev du point fixe 
c
      dimension SolT_nm1(NbCV)
      dimension SolT_prev(NbCV)
c
c     Temperature a t=0 
c      
      dimension SolT0(NbCV)
      
c
c        
c     Numerotation locale par cell ds l'ordre NodebyCell inc Ks
c               puis les inc face frac Ksigma si  k=NumCellbyFace(n,1) (cote maille 1 de la face) 
c     Local -> global 
c        
        dimension NbIncGlobalbyCell(NbCell)
        dimension NumIncGlobalbyCell(NbCell,NbIncbyCellMax)             
c
c   operateur gradient scalaire constant par maille 
c
      dimension GradCellVEM(NbCell,NbdofbyCellMax,NbDim)  
c
c     Operateur saut par face frac 
c
        dimension NbIncGlobalbyFaceFrac(NbFaceFrac)
        dimension NumIncGlobalbyFaceFrac(NbFaceFrac,2*NbNodebyFaceMax+2)           
        dimension SautbyFaceFrac(NbFaceFrac,2*NbNodebyFaceMax+2)
c
c     vecteur normal par face sortant maille 1 de NumCellbyFace 
c     
        dimension VecNormalbyFace(NbFace,NbDim)

c
c     Vecteurs unitaires tangents par face format un repere orthonorme avec le vec normal
c
        dimension VecTanFace1(NbFace,NbDim)
        dimension VecTanFace2(NbFace,NbDim)                  
c
c     coeff de relaxation ff
c
        dimension relaxff(NbFaceFrac)
c
c     coeff de relaxation mm
c
        dimension relaxmm(NbCell)
c
c        
cccccccccccccccccccccccccccccccccccccccc
c
c     SORTIES 
c
cccccccccccccccccccccccccccccccccccccccc
c
c     Solution En ENTREE ET EN SORTIE  
c
        dimension SolP(NbCV)

        dimension SolT(NbCV)
c
c
c     Acc Darcy et fourier  avec la solution courante (non convergee !!)  
c     donc a recalculer par appel a couplage et a JacSm
c       
      dimension AccP(NbCV)

c
      dimension AccT(NbCV)
c
c      
ccccccccccccccccccccccccccccccccccccccc
c
c     WorkSpaces 
c
c
c      Indices m CSR des j (indice col) des lignes K et I    
c
c       dimension NumCSRK(NbCV), NumCSRI(NbCV)
c
c     
c     
c     Indices m CSR des inc colonnes des lignes K et I ou Sig et I    
c     
      integer, dimension(:), allocatable :: NumCSRK
      integer, dimension(:), allocatable :: NumCSRI

c     
c     Jacobienne 
c     
      double precision, dimension(:,:,:), allocatable :: AA

c     
c     Jacobienne CSR mailles + noeuds 
c     On met l'elt diagonal en premier 

      integer, dimension(:), allocatable :: IndLigneAA
      integer, dimension(:), allocatable :: NLigneAA
      integer, dimension(:), allocatable :: NColAA

      integer, dimension(:), allocatable :: NbNzbyLine
      integer, dimension(:,:), allocatable :: NumNzbyLine
      
c     
c     Second membre 
c
      double precision, dimension(:,:), allocatable :: Sm
c
c     
c     increment 
c
      double precision, dimension(:), allocatable :: dSolP

      double precision, dimension(:), allocatable :: dSolT     
     
      double precision, dimension(:,:), allocatable :: dSolPT
c    
c
ccccccccccc
c      
       dimension X(NbDim)
c
c     par arete, nb et numeros des faces frac connectees 
c       
      integer, dimension(:), allocatable :: NbFaceFracbyArete
      integer, dimension(:,:), allocatable :: NumFaceFracbyArete
c
c
c       
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c     face frac voisine d'une arete 
c
       allocate(NbFaceFracbyArete(NbArete))
       allocate(NumFaceFracbyArete(NbArete,NbFacebyAreteMax))

       NbFaceFracbyArete(:) = 0
       NumFaceFracbyArete(:,:) = 0
       
       do k=1,NbFaceFrac

          nf = NumFaceFracVersFace(k)

          do ia = 1,NbAretebyFace(nf)
             na = NumAretebyFace(nf,ia)

             NbFaceFracbyArete(na) = NbFaceFracbyArete(na) + 1
             if (NbFaceFracbyArete(na).gt.NbFacebyAreteMax) then
                write(*,*)' trop de face frac par arete redim ',
     &               NbFaceFracbyArete(na),NbFacebyAreteMax
                stop
             endif

             NumFaceFracbyArete(na,NbFaceFracbyArete(na)) = k 
             
             
          enddo
          
       enddo       
c
c
c       do na=1,NbArete
c          if ( NbFaceFracbyArete(na).ne.0) then
c             write(*,*)NbFaceFracbyArete(na),na,
c     &            NumFaceFracbyArete(na,1:NbFaceFracbyArete(na))             
c          endif          
c       enddo
c       stop       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Calcul de la structure creuse CSR:  IndLigneAA,NColAA,NLigneAA
c     pour le FLOW
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
c     
       allocate(NbNzbyLine(NbCV))
       allocate(NumNzbyLine(NbCV,NbNzbyLineMax))


       call StructureCreuse(
     &      NbCell,NbFace,NbArete,NbFaceFrac,NbCV,
     &      NumFacebyCell,NumCellbyFace,
     &      NumFaceVersFaceFrac,IndFaceFrac,         
     &      NumIncCell,NumIncFaceFrac, 
     &      NbInterfacebyCell,NumInterfacebyCell,
     &      NbInterfacebyFaceFrac,NumInterfacebyFaceFrac,
     &      NumIncCellbyFaceFrac,
     &      NbFaceFracbyArete,NumFaceFracbyArete,
     &      NumFaceFracVersFace,NumAretebyFace,       
     &      NbNzbyLine,NumNzbyLine)
c     
c     IndLigneAA, NColAA, NLigneAA  
c     

       allocate(IndLigneAA(NbCV+1))
       IndLigneAA(1) = 0
       do i=1,NbCV
          IndLigneAA(i+1) = IndLigneAA(i) + NbNzbyLine(i)          
       enddo
       MemBloc = IndLigneAA(NbCV+1)
       
       allocate(NLigneAA(MemBloc))
       allocate(NColAA(MemBloc))
       do i=1,NbCV
          do m=IndLigneAA(i)+1,IndLigneAA(i+1)
             n = m-IndLigneAA(i)
             NColAA(m) = NumNzbyLine(i,n)
             NLigneAA(m) = i
          enddo
       enddo
c     
c     Permutation pour mettre l'elt diagonal en premier sur la ligne 
c     
       do i=1,NbCV
          md = IndLigneAA(i)+1
          jd = NColAA(md)
          do m = md,IndLigneAA(i+1)
             j = NColAA(m)
             if (j.eq.i) then 
                m1 = m
             endif
          enddo
          NColAA(md) = i
          NColAA(m1) = jd
       enddo

      


c      do i=1,NbCV+1
c      write(*,*)' IndLigneAA ',i,IndLigneAA(i)
c      enddo
c      write(*,*)
c      do i=1,NbCV
c      write(*,*)' NColAA ',i,(NColAA(m),
c     &              m=IndLigneAA(i)+1,IndLigneAA(i+1))
c      enddo
      
       
c
ccccccccccccccccccccccccccccccc

       
      allocate(AA(MemBloc,NbIncPrim,NbIncPrim))

       
      allocate(Sm(NbCV,NbIncPrim))
      allocate(dSolP(NbCV))
      allocate(dSolT(NbCV))
      allocate(NumCSRK(NbCV))
      allocate(NumCSRI(NbCV))


      allocate(dSolPT(NbCV,NbIncPrim))
      
c     
c
c
ccccccccc
c
c     calcul du residu initial 
c  
c      
c     calcul de la porosite et de l'aperture 
c

      
      call Compute_Porosity_Aperture(
     &     NbCell,NbFaceFrac,NbFace,NbCV,
     &     NbdofMecaContact,
     &     dfcontact,
     &     PorobyCell,DerPorobyCell,PorobyCell_nm1,
     &     dfbyFaceFrac, 
     &     SolU_nm1,SolU_prev,SolU0,
     &     IndContact_prev,Sautt_Proj_nm1,       
     &     SolP,SolP_nm1,SolP_prev,
     &     SolT,SolT_nm1,SolT_prev,
     &     NbIncGlobalbyCell,NumIncGlobalbyCell,GradCellVEM,     
     &     NbIncGlobalbyFaceFrac,NumIncGlobalbyFaceFrac,SautbyFaceFrac,
     &     NumFaceFracVersFace,
     &     VecNormalbyFace,VecTanFace1,VecTanFace2,
     &     NumIncCell)    
             
      call JacSmPT(
     &     Deltat,       
     &     NbCell,NbFaceFrac,NbFace,NbArete,NbCV,
     &     VolCell,VolFaceFrac,SurfaceFace,XInc,
     &     NumIncCell,NumIncFaceFrac, 
     &     IndIncDir,IndIncNeu,IndIncDirT,IndIncNeuT,NumIncFace,   
     &     NbInterfacebyCell,NumInterfacebyCell,
     &     NbInterfacebyFaceFrac,NumInterfacebyFaceFrac,
     &     NumIncCellbyFaceFrac,NumFaceFracVersFace,      
     &     NumFacebyCell,NumCellbyFace,
     &     NumFaceVersFaceFrac,IndFaceFrac,
     &     NbFaceFracbyArete,NumFaceFracbyArete,
     &     NumAretebyFace,     
     &     IndLigneAA,NColAA,NLigneAA,MemBloc, 
     &     TransCell,TransFaceFrac,TransCellbyFaceFrac,
     &     TransCellFourier,TransFaceFracFourier,
     &     TransCellbyFaceFracFourier,      
     &     PorobyCell,DerPorobyCell,PorobyCell_nm1,
     &     dfbyFaceFrac,dfbyFaceFrac_nm1,
     &     relaxff,relaxmm,
     &     AccP,AccP_nm1,       
     &     SolPDir,SolP,SolP_nm1,SolP_prev,SolP0,
     &     AccT,AccT_nm1,       
     &     SolTDir,SolT,SolT_nm1,SolT_prev,SolT0,    
     &     AA,Sm,
     &     NumCSRK,NumCSRI)

cccccccccccccccccccccccccccccccccccccccccccc
c
c     test solution 
c
c$$$      SolT(:) = 1.d0 
c$$$
c$$$      do i=1,NbCV
c$$$         s = 0.d0 
c$$$         do m=IndLigneAA(i)+1,IndLigneAA(i+1)
c$$$
c$$$            ncol = NColAA(m)
c$$$            
c$$$            s = s + AA(m,2,2)*SolT(ncol)             
c$$$            
c$$$         enddo
c$$$         if (dabs(s).ge.1.d-10) then
c$$$            write(*,*)
c$$$            write(*,*)
c$$$            write(*,*)' A22 Sol T ',i,IndIncDirT(i),IndIncNeuT(i),s
c$$$            write(*,*)
c$$$            write(*,*)' ligne i ',
c$$$     &           AA(IndLigneAA(i)+1:IndLigneAA(i+1),2,2) 
c$$$            write(*,*) 
c$$$         endif
c$$$      enddo
c$$$      stop
cccccccccccccccccccccccccccccccccccccccccccc
c      
c     calcul de la norme du residu initial 
c
      sp = 0.d0 
      do i=1,NbCV
         sp = sp + dabs(Sm(i,1))**2
      enddo
      sp = dsqrt(sp)
      resp0 = sp

      sT = 0.d0 
      do i=1,NbCV
         sT = sT + dabs(Sm(i,2))**2
      enddo
      sT = dsqrt(sT)
      resT0 = sT

c      res0 = resp0 + resT0
      
      write(*,*)
      write(*,*)' norme residu init ',resp0,resT0
c      write(*,*)      
c
c     criteres d'arret lineaire et non lineaire 
c      
      nit = 0
      nitmax = 100
      
      critere_arret_newton = 1.0d-7
c      critere_arret_newton = 1.0d-9

      critere_arret_gmres = 1.0d-6
      
      residu_relatif = 1.d0
      dXmax = 1.d0 
c      
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Boucle Newton 
c      
      do while ( (residu_relatif.ge.critere_arret_newton)
     &     .and.(nit.le.nitmax-1).and.(dXmax.ge.critere_arret_newton) )
c
c
c      do while (nit.le.1)         
c
c         
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


         
         nit = nit + 1
      
     
ccccccccccccccccccccccccccccc
c
c     Resolution du systeme lineaire 
c
cccccccccccccccccccccccccccc

      dSolPT= 0.d0 
                 
      ITER = 0
      NbIncSys = NbIncPrim
      call solveurSuperLU(NbCV,MemBloc,NbIncSys,
     &     AA,NLigneAA,NColAA,IndLigneAA,
     &     Sm,dSolPT)      


c      MAXL = 400
c      IPREC = 1
c      NbIncSys = NbIncPrim      
c      call solveur_iteratif(NbCV,MemBloc,NbIncSys,
c     &       AA,NLigneAA,NColAA,IndLigneAA,
c     &       Sm,dSolPT,critere_arret_gmres,MAXL,IPREC,
c     &       ITER,IERR,ERR)
      


      dSolP(:) = dSolPT(:,1)
      dSolT(:) = dSolPT(:,2)
  
      
ccccccccccccc      
c
c     calcul de l'increment max pour critere d'arret 
c
      dPmax = 0.d0
      dTmax = 0.d0       
      do i=1,NbCV
c         write(*,*)' dsolT ',i,dSolT(i)
         dPmax = dmax1(dabs(dSolP(i)),dPmax)
         dTmax = dmax1(dabs(dSolT(i)),dTmax)

c         if (dabs(dSolT(i)).gt.1000.d0) then
c            write(*,*)' dsolT ',i,dSolT(i),SolT(i)+dSolT(i)
c         endif
      enddo
c      write(*,*)' dTmax ',dTmax 


      
ccccccccccccccc  TMP TMP cccccccccc
      
      dXmax = dPmax/1.d+5 + dTmax/100.d0


c      dXmax = dPmax + dTmax ! TMP TRACEUR !!!!!       
c
cccccccccccccccc
c      
c
ccccccccccccccccccccccccccccccccccccc
c      
c      
c     Incrementation du Newton   
c

c      dTobj = 0.1d0
c      theta = dmin1(1.d0,dTobj/dTmax)
      theta = 1.d0
      
      SolP(:) = SolP(:) + dSolP(:)

      SolT(:) = SolT(:) + dSolT(:)*theta 
c
c
c      
cccccccccccccccccccccccccccccccccccc
c
c     Calcul du residu et de la Jacobienne 
c
ccccccccccccccccccccccccccccccccccc
c      
c     calcul de la porosite et de l'aperture 
c           
      call Compute_Porosity_Aperture(
     &     NbCell,NbFaceFrac,NbFace,NbCV,
     &     NbdofMecaContact,
     &     dfcontact,
     &     PorobyCell,DerPorobyCell,PorobyCell_nm1,
     &     dfbyFaceFrac, 
     &     SolU_nm1,SolU_prev,SolU0,
     &     IndContact_prev,Sautt_Proj_nm1,       
     &     SolP,SolP_nm1,SolP_prev,
     &     SolT,SolT_nm1,SolT_prev,  
     &     NbIncGlobalbyCell,NumIncGlobalbyCell,GradCellVEM,     
     &     NbIncGlobalbyFaceFrac,NumIncGlobalbyFaceFrac,SautbyFaceFrac,
     &     NumFaceFracVersFace,
     &     VecNormalbyFace,VecTanFace1,VecTanFace2,
     &     NumIncCell)    
          
       
      call JacSmPT(
     &     Deltat,       
     &     NbCell,NbFaceFrac,NbFace,NbArete,NbCV,
     &     VolCell,VolFaceFrac,SurfaceFace,XInc,
     &     NumIncCell,NumIncFaceFrac, 
     &     IndIncDir,IndIncNeu,IndIncDirT,IndIncNeuT,NumIncFace,  
     &     NbInterfacebyCell,NumInterfacebyCell,
     &     NbInterfacebyFaceFrac,NumInterfacebyFaceFrac,
     &     NumIncCellbyFaceFrac,NumFaceFracVersFace,
     &     NumFacebyCell,NumCellbyFace,
     &     NumFaceVersFaceFrac,IndFaceFrac,
     &     NbFaceFracbyArete,NumFaceFracbyArete,
     &     NumAretebyFace,     
     &     IndLigneAA,NColAA,NLigneAA,MemBloc, 
     &     TransCell,TransFaceFrac,TransCellbyFaceFrac,
     &     TransCellFourier,TransFaceFracFourier,
     &     TransCellbyFaceFracFourier,      
     &     PorobyCell,DerPorobyCell,PorobyCell_nm1,
     &     dfbyFaceFrac,dfbyFaceFrac_nm1,
     &     relaxff,relaxmm,
     &     AccP,AccP_nm1,       
     &     SolPDir,SolP,SolP_nm1,SolP_prev,SolP0,
     &     AccT,AccT_nm1,       
     &     SolTDir,SolT,SolT_nm1,SolT_prev,SolT0,    
     &     AA,Sm,
     &     NumCSRK,NumCSRI)
c
c     calcul de la norme du residu et residu relatif 
c
      sp = 0.d0 
      do i=1,NbCV
         sp = sp + dabs(Sm(i,1))**2
      enddo
      sp = dsqrt(sp)
      resp = sp

      sT = 0.d0 
      do i=1,NbCV
c         write(*,*)' T ',i,SolT(i),Sm(i,2)
         sT = sT + dabs(Sm(i,2))**2
      enddo
      sT = dsqrt(sT)
      resT = sT



      residu_relatif = resp/resp0 + resT/resT0
      
c      write(*,*)
       write(*,*)' res relatif dpTmax ITER ',
     &  nit,residu_relatif,dXmax,dPmax,dTmax,resp,resT,ITER 
c      write(*,*)

       
ccccccccccccccccccccccccccccccccccccc
      enddo ! fin boucle while  Newton 
cccccccccccccccccccccccccccccccccccc

c      stop
      
c      call CPU_TIME(time_jacsolve2)

c      write(*,*)' CPU Jac assembly + solve ',
c     &     time_jacsolve2-time_jacsolve1

ccccccccccccccccccccccccccccc


      deallocate(AA)
      deallocate(Sm)
      deallocate(dSolP)
      deallocate(dSolT)
      deallocate(NumCSRK)
      deallocate(NumCSRI)

      deallocate(dSolPT)

      deallocate(NbNzbyLine)
      deallocate(NumNzbyLine)      

      deallocate(IndLigneAA,NLigneAA,NColAA)
      
       deallocate(NbFaceFracbyArete)
       deallocate(NumFaceFracbyArete)      
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       return
       end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
