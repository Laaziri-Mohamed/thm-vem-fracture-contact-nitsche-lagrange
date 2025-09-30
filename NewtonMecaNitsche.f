
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     Modele d'elasticite lineaire isotrope avec fractures + contact sans frottement 
c
c     Champ u discretise avec VEM de degre 1 
c
c     Formulation Nitsche 
c
c
c
c     Elasticite + contact sans frottement 
c
c
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c      
      subroutine ElasticiteNitscheFrac(
     &     NbNode,NbCell,NbFace,NbFaceFrac,
     &     NbdofMeca,NbdofMecaContact,MemBlocU, 
     &     IndDirNodeMeca,IndFace,             
     &     NbNodebyCell,NumNodebyCell,
     &     NbNodebyFace,NumNodebyFace,
     &     NbFacebyCell,NumFacebyCell,
     &     NumCellbyFace,NumFaceFracVersFace,
     &     NumIncFaceFrac,NumIncCell,NumIncCellbyFaceFrac, 
     &     VecNormalbyFace,PoidsXCellCG,      
     &     XNode,XCellCG,XFaceCG,SurfaceFace,VolCell,      
     &     NbIncGlobalbyCell,NumIncGlobalbyCell,
     &     NumNodeFaceLocalbyCell,
     &     NbIncGlobalbyFaceFrac,NumIncGlobalbyFaceFrac,
     &     NbCV,SolP,
     &     SolT,
     &     SautbyFaceFrac,
     &     ACell,
     &     VecTanFace1,VecTanFace2,
     &     Sigmanbyff,
     &     PoidsXFaceCG,      
     &     IndBulleNumCellbyFace,
     &     GradTangbyFace,GradCellVEM,                       
     &     SolUL,SolUL_nm1,IndContact,nit)      
c
c
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
c
c
      include 'include_parameter'


c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     DECLARATIONS
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ENTREES       
c
c
c     coordonnees des noeuds XS 
c
      dimension XNode(NbNode,NbDim)
c
c     centre de maille 
c
      dimension XCellCG(NbCell,NbDim)
c
c     coordonnee du centre de gravite de la Face     
c
        dimension XFaceCG(NbFace,NbDim)
c
c     surface by face 
c        
        dimension SurfaceFace(NbFace)        
c
c     Volume des mailles 
c
        dimension VolCell(NbCell)
c
c     Numero des inconnues mailles et faces fractures 
c
c     num cell       >  num inc cell
c     num face frac  >  num inc face frac  
c
        dimension NumIncCell(NbCell)
        dimension NumIncFaceFrac(NbFaceFrac)

        
c        
c     Mailles par les nodes 
c
        dimension NbNodebyCell(NbCell)
        dimension NumNodebyCell(NbCell,NbNodebyCellMax)
c        
c     pour HFV discontinu 
c     num face frac  > num inc Ksigma et Lsigma 
c        
        dimension NumIncCellbyFaceFrac(NbFaceFrac,2)                
c
c     Maille par les faces 
c
        dimension NbFacebyCell(NbCell)
        dimension NumFacebyCell(NbCell,NbFacebyCellMax)       
c
c     Mailles voisines des faces 
c
        dimension NumCellbyFace(NbFace,2)        
c
c     Faces par les nodes 
c        
        dimension NbNodebyFace(NbFace)
        dimension NumNodebyFace(NbFace,NbNodebyFaceMax)      
c
c     numero de la face pour une face frac donnee  
c        
        dimension NumFaceFracVersFace(NbFaceFrac)
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
c     Poids aux nodes de la maille du CG de maille 
c        
        dimension PoidsXCellCG(NbCell,NbNodebyCellMax)

c
c      Poids barycentriques du CG de la face ds l'ordre des nodes by face 
        dimension PoidsXFaceCG(NbFace,NbNodebyFaceMax)
              
c
c
c        indice bulle dans la structure NumCellbyFace (cote maille 1 / 2) 
c
      dimension IndBulleNumCellbyFace(NbFace,2)

c
c     Operateur gradient tangentiel constant par face 
c
c     Grad(nf)_l = sum_{is \in V_nf} GradTangbyFace(nf,is,l)*U(inc)
c     inc = num global du node is de la face nf
c     l = 1,Nbdim 
c        
c        
      dimension GradTangbyFace(NbFace,NbNodebyFaceMax,NbDim)



cc   vecteur gradiant pour chaque fonction de bases scalaire v^dof (dof=1,...,NbdofbyCell(K)) pour chaque cellule K

      dimension GradCellVEM(NbCell,NbdofbyCellMax,NbDim) 

c     pressions dans la matrice et dans la fracture

      dimension SolP(NbCV)


c     temperature  dans la matrice 

      dimension SolT(NbCV)


c      
c     saut moyen 
c        

      dimension  SautbyFaceFrac(NbFaceFrac,2*NbNodebyFaceMax+2)

      
ccccccc
c        
c     Indices Dirichlet par node (1 si Dir, 0 sinon) 
c        
      dimension IndDirNodeMeca(NbNode)
c        
c     indice aux faces (1,2,,..,6 des 6 cotes) pour CL Neumann non homogene
c
      dimension IndFace(NbFace)
c
c      
c        
c
ccccccccc
c        
c     Numerotation locale par cell ds l'ordre NodebyCell inc Ks
c               puis les inc face frac Ksigma si  k=NumCellbyFace(n,1) (cote maille 1 de la face) 
c     Local -> global 
c        
        dimension NbIncGlobalbyCell(NbCell)
        dimension NumIncGlobalbyCell(NbCell,NbIncbyCellMax)        
c
c     Numero local du noeud ou de la face par inc cell 
c
        dimension NumNodeFaceLocalbyCell(NbCell,NbIncbyCellMax)
c
c     Structure Operateur saut par face frac 
c
        dimension NbIncGlobalbyFaceFrac(NbFaceFrac)
        dimension NumIncGlobalbyFaceFrac(NbFaceFrac,2*NbNodebyFaceMax+2)           

c
c     Matrice locale par maille VEM 
c
        dimension ACell(NbCell,NbdofbyCellMax,NbdofbyCellMax,
     &       NbDim,NbDim)

c
c     operateur sigman by face frac (tractions normale et tangentielles) 
c     sigman_nf(U,:) = sum_{dof,l} Sigmanbyff(nf,dof,l,:) U^dof_l
c
      dimension Sigmanbyff(NbFaceFrac,NbdofbyCellMax,NbDim,NbDim)       
c
c   
ccccccccccccccccccccccccccccccccc
c
c     Sorties 
c
c
c     Solution UL (init en entree et calcul en sortie)
c     contient le deplacement puis le multiplicateur lambda (post-processing)
c
      dimension SolUL(NbdofMecaContact,NbDim)
c
c     valeur au pas de temps precedent pour Coulomb 
c      
      dimension SolUL_nm1(NbdofMecaContact,NbDim)

      
c       
c     Indice de contact par face frac
c
       dimension IndContact(NbFaceFrac)     
c
c
c       
cccccccccccccccccccccccccccccccccccccccccccccc
c
c     Workspaces, declarations locales 
c
       dimension XXk(NbDim)

c
c
c
       dimension XXs(NbDim)
c       
c     SolU restreint aux dof U 
c
       double precision, dimension(:,:), allocatable :: SolU

c       
c     SolU_nm1 restreint aux dof U 
c
       double precision, dimension(:,:), allocatable :: SolU_nm1        
c
c     SolLambda
c
       double precision, dimension(:,:), allocatable :: SolLambda
c       
c
c     Jacobienne VEM complete 
c
c
c    structure creuse mecanique 
c
c        
       integer, dimension(:), allocatable :: IndLigneAU
       integer, dimension(:), allocatable :: NLigneAU
       integer, dimension(:), allocatable :: NColAU

       integer, dimension(:), allocatable :: NbNzbyLineMeca
       integer, dimension(:,:), allocatable :: NumNzbyLineMeca
       
       double precision, dimension(:,:,:), allocatable :: AU
       double precision, dimension(:,:,:), allocatable :: AU22 
       

c
c     Jacobienne VEM avec partie elastique + sigman*sigman (fixe) 
c
       double precision, dimension(:,:,:), allocatable :: AU0  
       
c
c     Second membre du systeme Meca 
c
       double precision, dimension(:,:), allocatable :: SmU

c
c     Residu du systeme Meca 
c
       double precision, dimension(:,:), allocatable :: ResU
       double precision, dimension(:,:), allocatable :: ResU2
       

c
c     Increment du Newton 
c
       double precision, dimension(:,:), allocatable :: dSolU
       double precision, dimension(:,:), allocatable :: dSolU2
       

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Workspaces 
c
       
c
c      Operateur Saut        
c
       dimension EvalXSautbyff(2*NbNodebyFaceMax+2)

c
c
c     numero global des inconnues de l'operateur
c     On numerote les inconnues locales dans l'ordre suivant:
c        inc maille k1       
c        nodes face sigma cote k2 puis bulle (si presente) cote k2        
c      
c      NbIncGlobalPThetaBeta
      
      dimension NumIncGlobalPThetaBeta(NbdofbyCellMax+NbNodebyFaceMax+1)      
c
c      
      dimension EvalXPthetabeta(
     &     NbdofbyCellMax+NbNodebyFaceMax+1,NbDim,NbDim)
c
c      
ccccccccccccccccccccccccccccccccccccccc
       
       dimension Vec1(NbDim),Vec2(NbDim),Vec3(NbDim),Xf(NbDim)
c
c     produit scalaire (ei,fj) ou ei base canonique, fj base locale (vecnormal,tau1,tau2)
c       
       dimension BaseFace(NbDim,NbDim) 
c
c     parametre de penalisation de Nitsche: betant(nf) = beta0nt/sqrt(surfaceFace(numf))
c
       double precision, dimension(:,:), allocatable :: BetabyFaceFrac

       double precision, dimension(:,:), allocatable :: SolU1
c
c
c       
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
c
c
 
c        theta = 1.d0
c        theta = 0.d0
         theta = -1.d0
       
       beta0n = f_beta0n()      ! parameter beta0 normal de penalisation de Nitsche
       beta0t = f_beta0t()      ! parameter beta0 tangent de penalisation de Nitsche

        write(*,*)' beta0 ',beta0n,beta0t
        
c
       allocate(BetabyFaceFrac(NbFaceFrac,2))

       do k=1,NbFaceFrac
          nf = NumFaceFracVersFace(k)
          BetabyFaceFrac(k,1) = beta0n/dsqrt(SurfaceFace(nf))
          BetabyFaceFrac(k,2) = beta0t/dsqrt(SurfaceFace(nf))
          
       enddo

      
cccccccccccccccccccccccccccccccccccccccc
c
c     Structure creuse de la Jacobienne Meca 
c
cccccccccccccccccccccccccccccccccccccccc
c
c       
      allocate(NbNzbyLineMeca(NbdofMeca))
      allocate(NumNzbyLineMeca(NbdofMeca,NbNzbyLineMecaMax))

      
c
c     init avec la diagonale first 
c      
      do k=1,NbdofMeca
         NbNzbyLineMeca(k) = 1
         NumNzbyLineMeca(k,NbNzbyLineMeca(k)) = k 
      enddo
c
c      
c     boucle sur les mailles et dof meca des mailles 
c
c      
      do k=1,NbCell
         do i=1,NbIncGlobalbyCell(k)
            ni = NumIncGlobalbyCell(k,i)
            do j=1,NbIncGlobalbyCell(k)
               nj = NumIncGlobalbyCell(k,j)
               indij = -1 ! test si ni,nj deja existant ds la ligne ni 
               do l=1, NbNzbyLineMeca(ni)
                  n = NumNzbyLineMeca(ni,l)
                  if (n.eq.nj) then
                     indij = 1 
                  endif
               enddo

               if (indij.eq.-1) then ! on rajoute l'elt ni nj 
                  NbNzbyLineMeca(ni) = NbNzbyLineMeca(ni)+1
                  if ( NbNzbyLineMeca(ni).gt.NbNzbyLineMecaMax) then
                     write(*,*)' increase NbNzbyLineMecaMax ',
     &                 NbNzbyLineMecaMax,NbNzbyLineMeca(ni)
                  endif
                  
                  NumNzbyLineMeca(ni,NbNzbyLineMeca(ni)) = nj                   
               endif

               
            enddo
         enddo
      enddo

c
c
c     On rajoute les termes de sauts [u][v] aux faces frac 
c
c
      do i=1,NbFaceFrac
         
         nf = NumFaceFracVersFace(i)

         numfacefrac = i 

         XXk = 0.d0
         theta0 = 0.d0
         betan = 0.d0  
         betat = 0.d0  

         call EvalXPThetaBetabyFaceFrac(
     &        XXk,numfacefrac,theta0,betan,betat,
     &        NbFaceFrac,NbFace,NbCell,
     &        NbNodebyFace,NumNodebyFace,
     &        NbNodebyCell,NumNodebyCell,       
     &        NbIncGlobalbyCell,
     &        NumIncGlobalbyCell,
     &        NumNodeFaceLocalbyCell,
     &        NumFaceFracVersFace,
     &        NumFacebyCell,NumCellbyFace,
     &        IndBulleNumCellbyFace,            
     &        PoidsXFaceCG,XFaceCG,
     &        GradTangbyFace, 
     &        NbIncGlobalbyFaceFrac,
     &        NumIncGlobalbyFaceFrac,
     &        Sigmanbyff,
     &        VecNormalbyFace,
     &        VecTanFace1,
     &        VecTanFace2,             
     &        EvalXPThetaBeta,
     &        NbIncGlobalPThetaBeta,
     &        NumIncGlobalPThetaBeta) 
         
 
         do ik=1,NbIncGlobalPThetaBeta       
            ni = NumIncGlobalPThetaBeta(ik)

            do jk=1,NbIncGlobalPThetaBeta       
               nj = NumIncGlobalPThetaBeta(jk)

               indij = -1       ! test si ni,nj deja existant ds la ligne ni 
               do l=1, NbNzbyLineMeca(ni)
                  n = NumNzbyLineMeca(ni,l)
                  if (n.eq.nj) then
                     indij = 1 
                  endif
               enddo

               if (indij.eq.-1) then ! on rajoute l'elt ni nj 
                  NbNzbyLineMeca(ni) = NbNzbyLineMeca(ni)+1
                  if ( NbNzbyLineMeca(ni).gt.NbNzbyLineMecaMax) then
                     write(*,*)' increase NbNzbyLineMecaMax ',
     &                 NbNzbyLineMecaMax,NbNzbyLineMeca(ni)
                  endif
                  
                  NumNzbyLineMeca(ni,NbNzbyLineMeca(ni)) = nj                   
               endif              
               
               
            enddo           
         enddo            
 
         
      enddo            
c
ccccccccccccc      

c      
      write(*,*)'NbdofMeca ',NbdofMeca
      write(*,*)'Nbfacefrac ',NbFaceFrac     

      
      allocate(IndLigneAU(NbdofMeca+1))
      IndLigneAU(1) = 0
      do i=1,NbdofMeca
         IndLigneAU(i+1) = IndLigneAU(i) + NbNzbyLineMeca(i)          
      enddo
      MemBlocU = IndLigneAU(NbdofMeca+1)


      write(*,*)' MemBlocU ',MemBlocU 
      
      allocate(NLigneAU(MemBlocU))
      allocate(NColAU(MemBlocU))
      

      do i=1,NbdofMeca
         do m=IndLigneAU(i)+1,IndLigneAU(i+1)
            n = m-IndLigneAU(i)
            NColAU(m) = NumNzbyLineMeca(i,n)
            NLigneAU(m) = i
         enddo
      enddo
         

      deallocate(NbNzbyLineMeca)
      deallocate(NumNzbyLineMeca)      
      
cccccccccccccccccccccccccccccccccccc
c
c     fin structure creuse meca 
c      
ccccccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccccccccccccccccccccc      
c
c     SolU init 
c
      allocate(SolU(NbdofMeca,NbDim))
c
c     copie de SolUL dans SolU (valeur init) 
c      
      SolU(:,:) = SolUL(1:NbdofMeca,:)
c
c     SolLambda, mise a zero 
c      
      allocate(SolLambda(NbFaceFrac,NbDim))      

      SolLambda = 0.d0 
c
c     SolU_nm1 
c
      allocate(SolU_nm1(NbdofMeca,NbDim))
      SolU_nm1(:,:) = SolUL_nm1(1:NbdofMeca,:)
c
c
c      
cccccccccccccccccccccccccccccccccccccccccccc      
c
c     Assemblage Jacobienne VEM \sum_K  aK  
c
cccccccccccccccccccccccccccccccccccccccccccc

      
      allocate(AU(MemBlocU,NbDim,NbDim))
      allocate(AU0(MemBlocU,NbDim,NbDim))

      AU = 0.d0 ! on met AU a zero 
      AU0 = 0.d0 ! on met AU0 a zero 
              
ccccccccccc
      do k = 1,NbCell
cccccccccc
      
         do i=1,NbIncGlobalbyCell(k)
            ni = NumIncGlobalbyCell(k,i)
            
            do j=1,NbIncGlobalbyCell(k)
               nj = NumIncGlobalbyCell(k,j)
            
               mij = 0          ! recherche elt ni nj dans la ligne ni du CSR 
               do m=IndLigneAU(ni)+1,IndLigneAU(ni+1)
                  n = NColAU(m)
                  if (n.eq.nj) then
                     mij = m 
                  endif
               enddo
c
               
               if (mij.eq.0) then
                  write(*,*)' elt ni nj non trouve dans CSR ',k,ni,nj
                  write(*,*)' ni nj ',ni,nj
                  write(*,*)' m1 m2 ',IndLigneAU(ni),IndLigneAU(ni+1)
                  do m=IndLigneAU(ni)+1,IndLigneAU(ni+1)
                     n = NColAU(m)
                     write(*,*)' ni n m ',ni,n,m
                  enddo                  
                  stop
               else
                 
                  AU0(mij,:,:) = AU0(mij,:,:) + ACell(k,i,j,:,:)
                                    
               endif

            
            enddo
         enddo
         
ccccccccccccc         
      enddo                     ! fin boucle mailles


ccccccccccccccccccccccccccccccccccc old ccccccccccccccccccccccccccccccc      
ccccccccccccc
c      
c     assemblage du terme - \int_Gamma theta/beta sigman(U) sigman(V) dans la Jacobienne 
c
       do nff = 1,NbFaceFrac
          nf = NumFaceFracVersFace(nff)
          k1 = NumCellbyFace(nf,1)
          k2 = NumCellbyFace(nf,2)           

 
          do ik = 1,NbIncGlobalbyCell(k1)
             ni = NumIncGlobalbyCell(k1,ik)
             do jk = 1,NbIncGlobalbyCell(k1)
                nj = NumIncGlobalbyCell(k1,jk)


                mij = 0         ! recherche elt ni nj dans la ligne ni du CSR 
                do m=IndLigneAU(ni)+1,IndLigneAU(ni+1)
                   n = NColAU(m)
                   if (n.eq.nj) then
                      mij = m 
                   endif
                enddo
c
               
                if (mij.eq.0) then
                   write(*,*)' elt ni nj non trouve dans CSR ',k1,ni,nj
                   write(*,*)' ni nj ',ni,nj
                   write(*,*)' m1 m2 ',IndLigneAU(ni),IndLigneAU(ni+1)
                   do m=IndLigneAU(ni)+1,IndLigneAU(ni+1)
                      n = NColAU(m)
                      write(*,*)' ni n m ',ni,n,m
                   enddo                  
                   stop
                else
                   
                   do id = 1,Nbdim
                      do jd = 1,Nbdim
                         
                         sn = Sigmanbyff(nff,ik,id,1)
     &                        *Sigmanbyff(nff,jk,jd,1)
                         
                         sn = - sn*SurfaceFace(nf)
     &                        *theta/BetabyFaceFrac(nff,1)


                         st = Sigmanbyff(nff,ik,id,2)
     &                        *Sigmanbyff(nff,jk,jd,2)
                         st = st + Sigmanbyff(nff,ik,id,3)
     &                        *Sigmanbyff(nff,jk,jd,3)
                         
                         
                         st = - st*SurfaceFace(nf)
     &                        *theta/BetabyFaceFrac(nff,2)
                                                  
                         AU0(mij,id,jd) = AU0(mij,id,jd) + sn + st  
                      enddo
                   enddo
                   
                endif
                
               
                
             enddo
          enddo
          
          
       enddo
c
c
c      
cccccccccccccccccccccccccccccccccccc      
c
c     Calcul du second membre SmU 
c
ccccccccccccccccccccccccccccccccccc
c
      allocate(SmU(NbdofMeca,NbDim))

      
      SmU = 0.d0 
c      
c     Quadrature a un point au centre de gravite:
c
c       \int_K f \PiK v ~ |K| f(xCGK) * sum_s PoidsXCellCG(k,is)*Vks  
c
      do k=1,NbCell

         XXk(:) = XCellCG(k,:)
         
         Vec1 = 0.d0           
         Vec1(3) = - f_rhosf()*f_gravite()         

         do ik=1,NbNodebyCell(k)
            inc = NumIncGlobalbyCell(k,ik)
            SmU(inc,:) = SmU(inc,:)
     &          + PoidsXCellCG(k,ik)*VolCell(k)*Vec1(:)
         enddo
         
      enddo
c
c
c
c
cccccccccccccccccccccccccc      
c
c calcul du terme (dans le second membre) \sum_{sigma} \int_{sigma} pf(sigma)*saut(v)_n
c
cccccccccccccccccccccccccc
c
c
c      
      do ifrac = 1,NbFaceFrac
         nf = NumFaceFracVersFace(ifrac)
         inck = NumIncFaceFrac(ifrac)
     
         do ik=1,NbIncGlobalbyFaceFrac(ifrac)
            inc = NumIncGlobalbyFaceFrac(ifrac,ik)

              SmU(inc,:) = SmU(inc,:)
     &              - SolP(inck)*SurfaceFace(nf)
     &              * SautbyFaceFrac(ifrac,ik)*VecNormalbyFace(nf,:)       
          
                        
         enddo         
         
      enddo
c
cccccccccccccccccccccc      
c
c  calcul du terme (dans le second membre) \sum_k int_{K} b*pm div(v) 
c
ccccccccccccccccccc
c
c
c      
      b_coef = f_biot()

      alpha0 = f_alpha0()

      rK0 =  f_K0()

      Tref = f_TREF()


        do k=1,NbCell

           inck = NumIncCell(k) 

           do ik=1,NbIncGlobalbyCell(k)

               inc = NumIncGlobalbyCell(k,ik)

                 SmU(inc,:) = SmU(inc,:)
     &              + VolCell(k)* ( b_coef*SolP(inck)
     &              + alpha0 * rK0 * ( SolT(inck) - Tref) )
     &               *GradCellVEM(k,ik,:)

            enddo
         
         enddo


      
ccccccccccccccccccccccccccccccccccccccccccccc
c
c     Rajout des termes sources de type Neumann non homogene 
c      
ccccccccccccccccccccccccccccccccccccccccccccc
c
c 
c      
      do i=1,NbFace

         Xf(:) = XFaceCG(i,:)

ccccccccc  cote z=zmax 
         
         if (IndFace(i).eq.6) then 
            k = NumCellbyFace(i,1)
            nnk = NbNodebyCell(k)
          
            do in=1,NbNodebyFace(i)
               n = NumNodebyFace(i,in)
               
                icell = - 1 
                do is=1,nnk
                   if(NumNodebyCell(k,is).eq.n) then
                      icell = is 
                   endif  
                enddo
                if (icell.eq.-1) then
                   write(*,*)' on ne trouve pas le node
     &                         de la face ds GradvK ',
     &                  k,i 
                   stop
                endif               

                inc = NumIncGlobalbyCell(k,icell)           

                poids =  - f_rhosf()*f_gravite()*f_htop()  
                

                SmU(inc,:) = SmU(inc,:)
     &               + poids*SurfaceFace(i)*PoidsXFaceCG(i,in)
     &                      *VecNormalbyFace(i,:)               


c                write(*,*)' poids ',SmU(inc,:)
                
            enddo            
         endif

cccccccccc   cote x = xmax 

         if (IndFace(i).eq.2) then 
            k = NumCellbyFace(i,1)
            nnk = NbNodebyCell(k)

            
            do in=1,NbNodebyFace(i)
               n = NumNodebyFace(i,in)
               
                icell = - 1 
                do is=1,nnk
                   if(NumNodebyCell(k,is).eq.n) then
                      icell = is 
                   endif  
                enddo
                if (icell.eq.-1) then
                   write(*,*)' on ne trouve pas le node
     &                         de la face ds GradvK ',
     &                  k,i 
                   stop
                endif               

                inc = NumIncGlobalbyCell(k,icell)           


                ratiohv = f_ratiohv() ! ration contrainte h sur v                
                
                poids =  - ratiohv*f_rhosf()*f_gravite()
     &                     *(-Xf(3) + f_htop() + 2000.d0 ) 

                
                SmU(inc,:) = SmU(inc,:)
     &               + poids*SurfaceFace(i)*PoidsXFaceCG(i,in)
     &                      *VecNormalbyFace(i,:)               


c                write(*,*)' poids ',SmU(inc,:)
                
            enddo            
         endif

ccccccccc  cote y = ymax          
c$$$
c$$$         if (IndFace(i).eq.4) then 
c$$$            k = NumCellbyFace(i,1)
c$$$            nnk = NbNodebyCell(k)
c$$$
c$$$            
c$$$            do in=1,NbNodebyFace(i)
c$$$               n = NumNodebyFace(i,in)
c$$$               
c$$$                icell = - 1 
c$$$                do is=1,nnk
c$$$                   if(NumNodebyCell(k,is).eq.n) then
c$$$                      icell = is 
c$$$                   endif  
c$$$                enddo
c$$$                if (icell.eq.-1) then
c$$$                   write(*,*)' on ne trouve pas le node
c$$$     &                         de la face ds GradvK ',
c$$$     &                  k,i 
c$$$                   stop
c$$$                endif               
c$$$
c$$$                inc = NumIncGlobalbyCell(k,icell)           
c$$$
c$$$
c$$$                rationhv = 0.7d0
c$$$                poids =  - rationhv*f_rhosf()*f_gravite()
c$$$     &                     *(-Xf(3) + f_htop() + 2000.d0 ) 
c$$$
c$$$                
c$$$                SmU(inc,:) = SmU(inc,:)
c$$$     &               + poids*SurfaceFace(i)*PoidsXFaceCG(i,in)
c$$$     &                      *VecNormalbyFace(i,:)               
c$$$
c$$$
c$$$c                write(*,*)' poids ',SmU(inc,:)
c$$$                
c$$$            enddo            
c$$$         endif
         
         
      enddo
c      stop
c      
c         
c
cccccccccccccccccccccc         
c         
c   Calcul du terme (dans le scecond membre) \sum_{sigma} \int_{sigma}theta/beta (pf-b pm)*(sigma_n(v) )
c
ccccccccccccccccccccc
c         
        do ifrac = 1,NbFaceFrac
           nf = NumFaceFracVersFace(ifrac)
           inck = NumIncFaceFrac(ifrac)

           kcell = NumCellbyFace(nf,1)
           inccell = NumIncCell(kcell)

          incik1 = NumIncCellbyFaceFrac(ifrac,1)           
c
c        on prend la valeur de  sur l'interface K1sigma 
c           
          pTinterface = - b_coef*SolP(incik1) + SolP(inck)
     &     - alpha0 * rK0 * ( SolT(incik1) - Tref)     
c           pinterface = - b_coef*SolP(inccell) + SolP(inck)  
           
           do ik=1,NbIncGlobalbyFaceFrac(ifrac)
             inc = NumIncGlobalbyFaceFrac(ifrac,ik)

              SmU(inc,:) = SmU(inc,:)
     &           +  pTinterface*SurfaceFace(nf)
     &          *Sigmanbyff(ifrac,ik,:,1)*theta/BetabyFaceFrac(ifrac,1)

                     
           enddo                  
        enddo           
c
c
c
ccccccccccccccccccccccccccccccc
c
c     Residu initial + Jacobienne
c
ccccccccccccccccccccccccccccccc
c      
      allocate(ResU(NbdofMeca,NbDim))
c
cccccccccccccccccccccc      
c      
      ResU = 0.d0 
      AU = AU0


      call ResiduJacContact(
     &     theta,BetabyFaceFrac,SurfaceFace,      
     &     NbFaceFrac,NbFace,NbCell,NbNode,
     &     XNode,XCellCG,           
     &     NbNodebyFace,NumNodebyFace,
     &     NbNodebyCell,NumNodebyCell,       
     &     NbIncGlobalbyCell,
     &     NumIncGlobalbyCell,
     &     NumNodeFaceLocalbyCell,
     &     NumFacebyCell,NumCellbyFace,
     &     NumIncFaceFrac,NumIncCell,NumIncCellbyFaceFrac,  
     &     IndBulleNumCellbyFace,
     &     IndDirNodeMeca,      
     &     PoidsXFaceCG,
     &     GradTangbyFace,       
     &     Sigmanbyff,      
     &     NbdofMeca,MemBlocU,  
     &     NbIncGlobalbyFaceFrac,NumIncGlobalbyFaceFrac,
     &     XFaceCG,NumFaceFracVersFace,
     &     VecNormalbyFace,
     &     VecTanFace1,VecTanFace2,
     &     SolP, NbCV,
     &     SolT, 
     &     IndLigneAU,NLigneAU,NColAU,AU,       
     &     SolU,SolU_nm1,SmU,
     &     IndContact,ResU,SolLambda)      



ccccccccccccccccccccccccccccccccccccccccccccccccc     
c
c     norme l2 du residu init 
c
      s = 0.d0 
      do inc=1,NbdofMeca
         s = s + prodscal(ResU(inc,:),ResU(inc,:))
      enddo
      resinit = dsqrt(s) 

      write(*,*)' Res init ',resinit
      
      residu_relatif = 1.d0 


      
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Boucle semi-smooth Newton 
c

      nit = 0
      nitmax = 200
      critere_arret_newton = 1.0E-9
      critere_arret_dUmax = 1.0E-8      
      
      dUmax = 1.d0 

      
      allocate(dSolU(NbdofMeca,NbDim))
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
      do while (
c     &     (residu_relatif.ge.critere_arret_newton)
     &     (dUmax.ge.critere_arret_dUmax)           
     &     .and.(nit.le.nitmax-1) )
c         
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         nit = nit + 1
         
ccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccc
c     On resoud le systeme lineaire avec -ResU en second membre -> dSolU 
ccccccccccccccccccccccccccccccccc
c      
c      critere_arret_gmres = 1.0E-10

      
      ResU = - ResU
     
c$$$
c$$$      ITER = 0
c$$$      NbIncMeca = NbDim 
c$$$      call solveurSuperLU(NbdofMeca,MemBlocU,NbIncMeca,
c$$$     &     AU,NLigneAU,NColAU,IndLigneAU,
c$$$     &     ResU,dSolU)
c$$$


c      do k=1,NbdofMeca
c         write(*,*)' dSolU ',k,dSolU(k)
c      enddo
c      stop      
c

c
c
cccccccccccccccccccccccccccccccc      
c
c     reduction du systeme (2,2) pour cas 2D -> avec eq et inc 1 et 3 ici  
c      
      allocate(ResU2(NbdofMeca,2))
      allocate(dSolU2(NbdofMeca,2))

      ResU2(:,1) = ResU(:,1)
      ResU2(:,2) = ResU(:,3)

      allocate(AU22(MemblocU,2,2))

      AU22(:,1,1) = AU(:,1,1)
      AU22(:,1,2) = AU(:,1,3)
      AU22(:,2,1) = AU(:,3,1)
      AU22(:,2,2) = AU(:,3,3)

      ITER = 0
      NbIncMeca = NbDim - 1
      call solveurSuperLU(NbdofMeca,MemBlocU,NbIncMeca,
     &     AU22,NLigneAU,NColAU,IndLigneAU,
     &     ResU2,dSolU2)

      dSolU(:,1) = dSolU2(:,1)
      dSolU(:,2) = 0.d0 
      dSolU(:,3) = dSolU2(:,2)

      deallocate(ResU2,dSolU2,AU22)
c
cccccccccccccccccccccccccccccccccccc
      
cccccccccccccccccccccccccccccccc      
c
c     incrementation de Newton      
c
c     relaxation 
c
      dUmax = 0.d0
      do k=1,NbdofMeca
         dUmax = dmax1(dSolU(k,1),dUmax)
         dUmax = dmax1(dSolU(k,2),dUmax)
         dUmax = dmax1(dSolU(k,3),dUmax)        
      enddo

      write(*,*)' dUmax ',dUmax
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
c         irelax = 1
c          irelax = 2
           irelax = 3


      if (irelax.eq.1) then 
c
ccccccccc relaxation avec dUobj fixe ccccccccccccccccccccccccccc
c
      
c         dUobj = 5.d-6         
c         relax = dmin1(dUobj/dUmax,1.d0)
         relax = 1.d0         
      
         write(*,*)'dUmax ',dUmax,' relax ',relax
         if (nit.eq.1) then
            relax = 1.d0 
         endif
         
         SolU = SolU + dSolU*relax 


      else if (irelax.eq.2) then 

ccccccccccc line search basique ccccccccccccccccccccccccccccccccc
      
      allocate(SolU1(NbdofMeca,NbDim))
      SolU1 = 0.d0

c      write(*,*)
      smin = 1.d+20 
      do ils = 1,11
c      do ils = 1,4

         relax1 = 2.d0**(1-ils)
         
         SolU1 = SolU + dSolU*relax1 
c
ccccccccccccccccccccccccccccccc      
c
c     calcul du residu + Jacobienne   
c

         AU = AU0
     

      call ResiduJacContact(
     &     theta,BetabyFaceFrac,SurfaceFace,      
     &     NbFaceFrac,NbFace,NbCell,NbNode,
     &     XNode,XCellCG,           
     &     NbNodebyFace,NumNodebyFace,
     &     NbNodebyCell,NumNodebyCell,       
     &     NbIncGlobalbyCell,
     &     NumIncGlobalbyCell,
     &     NumNodeFaceLocalbyCell,
     &     NumFacebyCell,NumCellbyFace,
     &     NumIncFaceFrac,NumIncCell,NumIncCellbyFaceFrac,  
     &     IndBulleNumCellbyFace,
     &     IndDirNodeMeca,      
     &     PoidsXFaceCG,
     &     GradTangbyFace,       
     &     Sigmanbyff,      
     &     NbdofMeca,MemBlocU,  
     &     NbIncGlobalbyFaceFrac,NumIncGlobalbyFaceFrac,
     &     XFaceCG,NumFaceFracVersFace,
     &     VecNormalbyFace,
     &     VecTanFace1,VecTanFace2,
     &     SolP, NbCV,
     &     SolT,
     &     IndLigneAU,NLigneAU,NColAU,AU,       
     &     SolU1,SolU_nm1,SmU,
     &     IndContact,ResU,SolLambda)      
         

CCCCCCCCCCCCC      
      ResU(:,2) = 0.d0          ! TMP 2D 
CCCCCCCCCCCCC
      
         s = 0.d0 
         do inc=1,NbdofMeca
            s = s + prodscal(ResU(inc,:),ResU(inc,:))
         enddo
         s = dsqrt(s)
         
c         write(*,*)relax1,s

         if (s.le.smin) then
            smin = s
            relax = relax1 
         endif

      enddo

      if (nit.eq.1) then
         relax = 1.d0 
      endif
      
c      write(*,*)
      write(*,*)' relax ',relax 


      
      deallocate(SolU1)


      SolU = SolU + dSolU*relax

cccccccccccccccccccccccccccccccc fin line search basique ccccccccccc     
c


      
      else if (irelax.eq.3) then
         
cccccccccccccccccccccccccccccccc backtracking line search  ccccccccccc     


         allocate(SolU1(NbdofMeca,NbDim))

         SolU1 = 0.d0


         relax1 = 1.d0

c     parameters alpha /in (0,1/2) and beta /in (0,1)


        rAlpha = 1.d-5
c        rBeta= 0.67d0
c     rBeta= 0.5d0

        rBeta= 0.6d0

CCCCCCCCCCC        
      ResU(:,2) = 0.d0          ! TMP 2D 
CCCCCCCCCCC
      
        sInit = 0.d0
        do inc=1,NbdofMeca
           sInit = sInit + prodscal(ResU(inc,:),ResU(inc,:))
        enddo



        rcondition = sInit - 2*rAlpha*relax1*sInit

        SolU1 = SolU + dSolU*relax1



ccccccccccccccccccccccccccccccc      
c
c     calcul du residu + Jacobienne   
c

       AU = AU0  



      call ResiduJacContact(
     &     theta,BetabyFaceFrac,SurfaceFace,      
     &     NbFaceFrac,NbFace,NbCell,NbNode,
     &     XNode,XCellCG,           
     &     NbNodebyFace,NumNodebyFace,
     &     NbNodebyCell,NumNodebyCell,       
     &     NbIncGlobalbyCell,
     &     NumIncGlobalbyCell,
     &     NumNodeFaceLocalbyCell,
     &     NumFacebyCell,NumCellbyFace,
     &     NumIncFaceFrac,NumIncCell,NumIncCellbyFaceFrac,  
     &     IndBulleNumCellbyFace,
     &     IndDirNodeMeca,      
     &     PoidsXFaceCG,
     &     GradTangbyFace,       
     &     Sigmanbyff,      
     &     NbdofMeca,MemBlocU,  
     &     NbIncGlobalbyFaceFrac,NumIncGlobalbyFaceFrac,
     &     XFaceCG,NumFaceFracVersFace,
     &     VecNormalbyFace,
     &     VecTanFace1,VecTanFace2,
     &     SolP, NbCV,
     &     SolT,
     &     IndLigneAU,NLigneAU,NColAU,AU,       
     &     SolU1,SolU_nm1,SmU,
     &     IndContact,ResU,SolLambda) 

CCCCCCCCCCCC
      ResU(:,2) = 0.d0 ! TMP 2D 
CCCCCCCCCCCC
      
      s = 0.d0
      do inc=1,NbdofMeca
         s = s + prodscal(ResU(inc,:),ResU(inc,:))
      enddo



          relax = relax1



         iBsteps = 0

c         iBstepsmax = 5

         iBstepsmax = 10





c      Backtracking line search loop

         do while((s.gt.rcondition).and.(iBsteps.lt.iBstepsmax)
     &         .and.(nit.ne.1))




            iBsteps= iBsteps +1

c        nbiBsteps= nbiBsteps + iBsteps

         relax = rBeta*relax


         rcondition = sInit - 2*rAlpha*relax*sInit


         SolU1 = SolU + dSolU*relax




ccccccccccccccccccccccccccccccc      
c
c     calcul du residu + Jacobienne   
c

      AU = AU0   


       call ResiduJacContact(
     &     theta,BetabyFaceFrac,SurfaceFace,      
     &     NbFaceFrac,NbFace,NbCell,NbNode,
     &     XNode,XCellCG,           
     &     NbNodebyFace,NumNodebyFace,
     &     NbNodebyCell,NumNodebyCell,       
     &     NbIncGlobalbyCell,
     &     NumIncGlobalbyCell,
     &     NumNodeFaceLocalbyCell,
     &     NumFacebyCell,NumCellbyFace,
     &     NumIncFaceFrac,NumIncCell,NumIncCellbyFaceFrac,  
     &     IndBulleNumCellbyFace,
     &     IndDirNodeMeca,      
     &     PoidsXFaceCG,
     &     GradTangbyFace,       
     &     Sigmanbyff,      
     &     NbdofMeca,MemBlocU,  
     &     NbIncGlobalbyFaceFrac,NumIncGlobalbyFaceFrac,
     &     XFaceCG,NumFaceFracVersFace,
     &     VecNormalbyFace,
     &     VecTanFace1,VecTanFace2,
     &     SolP, NbCV,
     &     SolT,
     &     IndLigneAU,NLigneAU,NColAU,AU,       
     &     SolU1,SolU_nm1,SmU,
     &     IndContact,ResU,SolLambda)      

CCCCCCCCCCCCCC
      ResU(:,2) = 0.d0 ! TMP 2D       
CCCCCCCCCCCCC
      
      s = 0.d0
      do inc=1,NbdofMeca
         s = s + prodscal(ResU(inc,:),ResU(inc,:))
      enddo

      enddo




        if (nit.eq.1) then
           relax = 1.d0 
        endif




c      write(*,*)
       write(*,*)' relax ',relax

       write(*,*)' nb Backtracking steps ',iBsteps




         deallocate(SolU1)


        SolU = SolU + dSolU*relax


         

cccccccccccccccccccccccccccccccc fin backtracking line search  ccccccccccc     
      
      else

         write(*,*)' relax option adoit etre 1 or 2 '
         write(*,*)' irelax ',irelax
         stop

      endif ! fin if relax = 1 ou 2       
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c     calcul du residu + Jacobienne 
c

         AU = AU0
      
         call ResiduJacContact(
     &        theta,BetabyFaceFrac,SurfaceFace,      
     &        NbFaceFrac,NbFace,NbCell,NbNode,
     &        XNode,XCellCG,            
     &        NbNodebyFace,NumNodebyFace,
     &        NbNodebyCell,NumNodebyCell,       
     &        NbIncGlobalbyCell,
     &        NumIncGlobalbyCell,
     &        NumNodeFaceLocalbyCell,
     &        NumFacebyCell,NumCellbyFace,
     &        NumIncFaceFrac,NumIncCell,NumIncCellbyFaceFrac, 
     &        IndBulleNumCellbyFace,
     &        IndDirNodeMeca,       
     &        PoidsXFaceCG,
     &        GradTangbyFace,       
     &        Sigmanbyff,      
     &        NbdofMeca,MemBlocU,  
     &        NbIncGlobalbyFaceFrac,NumIncGlobalbyFaceFrac,
     &        XFaceCG,NumFaceFracVersFace,
     &        VecNormalbyFace,
     &        VecTanFace1,VecTanFace2,
     &        SolP,NbCV,
     &        SolT,
     &        IndLigneAU,NLigneAU,NColAU,AU,       
     &        SolU,SolU_nm1,SmU,
     &        IndContact,ResU,SolLambda)            
c
ccccccccccccccccccccccccccccccccccccccccccccccccc           
c
c     norme l2 du residu
c
      s = 0.d0 
      do inc=1,NbdofMeca 
         s = s + prodscal(ResU(inc,:),ResU(inc,:))
c         write(*,*)' res ',inc,ResU(inc,:)
      enddo
      res = dsqrt(s) 

      residu_relatif = res/resinit
      
c      write(*,*)
      write(*,*)' it newton ',nit,res,resinit,residu_relatif,ITER
c      write(*,*)       
      
ccccccccccccccccccccccccccccccccccccc
      enddo ! fin boucle while  Newton 
cccccccccccccccccccccccccccccccccccc
c
c
c
c      
ccccccccccccccccccccccccccccccccccccc
c      
c     Copie SolU et SolLambda dans SolUL
c
ccccccccccccccccccccccccccccccccccccc
c      
      SolUL(1:NbdofMeca,:) = SolU(:,:)

      SolUL(NbdofMeca+1:NbdofMecaContact,:) = SolLambda(:,:)

      deallocate(SolU,SolU_nm1,SolLambda)
c
c
c
c      
c      
ccccccccccccccccccccccccccccccccccc

      deallocate(AU,AU0,SmU,ResU,dSolU)
      deallocate(BetabyFaceFrac)      

      deallocate(NLigneAU,NColAU,IndLigneAU)
c       
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     Elasticite + contact sans frottement 
c
c     partie lineaire sur i=1,NbdofMeca calculee par produit matrice vecteur 
c
c
c     Termes de contact a rajouter 
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c      
      subroutine ResiduJacContact(
     &     theta,BetabyFaceFrac,SurfaceFace,      
     &     NbFaceFrac,NbFace,NbCell,NbNode,
     &     XNode,XCellCG,      
     &     NbNodebyFace,NumNodebyFace,
     &     NbNodebyCell,NumNodebyCell,       
     &     NbIncGlobalbyCell,
     &     NumIncGlobalbyCell,
     &     NumNodeFaceLocalbyCell,
     &     NumFacebyCell,NumCellbyFace,
     &     NumIncFaceFrac,NumIncCell,NumIncCellbyFaceFrac, 
     &     IndBulleNumCellbyFace,
     &     IndDirNodeMeca,
     &     PoidsXFaceCG,
     &     GradTangbyFace,       
     &     Sigmanbyff,      
     &     NbdofMeca,MemBlocU,  
     &     NbIncGlobalbyFaceFrac,NumIncGlobalbyFaceFrac,
     &     XFaceCG,NumFaceFracVersFace,
     &     VecNormalbyFace,
     &     VecTanFace1,VecTanFace2,
     &     SolP,NbCV,
     &     SolT,
     &     IndLigneAU,NLigneAU,NColAU,AU,       
     &     SolU,SolU_nm1,SmU,
     &     IndContact,ResU,SolLambda)      
c
c
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
c
c
      include 'include_parameter'


c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     DECLARATIONS
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ENTREES       
c
ccccccccc
c
c     noeuds 
c
      dimension XNode(NbNode,Nbdim)

c
c     centre de maille 
c
      dimension XCellCG(NbCell,NbDim)      
c        
c     nb de noeuds par maille K 
c
        dimension NbNodebyCell(NbCell)
c
c     num des noeuds de chaque maille K dans l'ordre cyclique 
c
        dimension NumNodebyCell(NbCell,NbNodebyCellMax)
c
c      le nombre de noeuds de chaque face
        dimension NbNodebyFace(NbFace)
        dimension NumNodebyFace(NbFace,NbNodebyFaceMax)       
c
c     numero des faces par maille 
c        
        dimension NumFacebyCell(NbCell,NbFacebyCellMax)
c        
c
c      Poids barycentriques du CG de la face ds l'ordre des nodes by face 
        dimension PoidsXFaceCG(NbFace,NbNodebyFaceMax)
      

c
c     Operateur gradient tangentiel constant par face 
c
c     Grad(nf)_l = sum_{is \in V_nf} GradTangbyFace(nf,is,l)*U(inc)
c     inc = num global du node is de la face nf
c     l = 1,Nbdim 
c        
c        
        dimension GradTangbyFace(NbFace,NbNodebyFaceMax,NbDim)        

        
c     Numerotation locale par cell ds l'ordre NodebyCell inc Ks
c               puis les inc face frac Ksigma si  k=NumCellbyFace(n,1) (cote maille 1 de la face) 
c     Local -> global 
c        
        dimension NbIncGlobalbyCell(NbCell)
        dimension NumIncGlobalbyCell(NbCell,NbIncbyCellMax)        
c
c     Numero local du noeud ou de la face par inc cell 
c
        dimension NumNodeFaceLocalbyCell(NbCell,NbIncbyCellMax)

c
c     mailles voisines d'une face 
c        
        dimension NumCellbyFace(NbFace,2)
c
c
c
        dimension NumIncCell(NbCell)        
        dimension NumIncFaceFrac(NbFaceFrac)
c        
c     pour HFV discontinu 
c     num face frac  > num inc Ksigma et Lsigma 
c        
        dimension NumIncCellbyFaceFrac(NbFaceFrac,2)        
c
c        indice bulle dans la structure NumCellbyFace (cote maille 1 / 2) 
c
        dimension IndBulleNumCellbyFace(NbFace,2)
        
c
c     operateur sigman by face frac (traction normale et tangentielle) 
c     sigman_nf(U,:) = sum_{dof,l} Sigmanbyff(nf,dof,l,:) U^dof_l    
c
c        
        dimension Sigmanbyff(NbFaceFrac,NbdofbyCellMax,NbDim,NbDim)


cccccccc      
c
c     structure operateur saut par face frac 
c
        dimension NbIncGlobalbyFaceFrac(NbFaceFrac)
        dimension NumIncGlobalbyFaceFrac(NbFaceFrac,2*NbNodebyFaceMax+2)           


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
c     coordonnee du centre de gravite de la Face     
c
        dimension XFaceCG(NbFace,NbDim)

c
c     numero de la face pour une face frac donnee  
c        
        dimension NumFaceFracVersFace(NbFaceFrac)        
c
c     CSR matrice meca et matrice sans les eqs du mulriplicateur 
c
       dimension IndLigneAU(NbdofMeca+1)
       dimension NLigneAU(MemBlocU)
       dimension NColAU(MemBlocU)

       dimension AU(MemBlocU,NbDim,NbDim)          
c
c
c     Solution U courante du Newton 
c
       dimension SolU(NbdofMeca,NbDim)

c
c     Solution U nm1 
c
       dimension SolU_nm1(NbdofMeca,NbDim)       
c
c
c     Solution Pression
c
       dimension SolP(NbCV)


c
c     Solution temperature
c
       dimension SolT(NbCV)
       
       


       
c
c     second membre de l'equation (or equation de "contact") 
c
       dimension SmU(NbdofMeca,NbDim)
c
c     penalisation de Nitsche 
c       
       dimension BetabyFaceFrac(NbFaceFrac,2)
c
c     surface by face 
c        
       dimension SurfaceFace(NbFace)

ccccccc
c        
c     Indices Dirichlet par node et par face (1 si Dir, 0 sinon) 
c        
      dimension IndDirNodeMeca(NbNode)
c       
c       
ccccccccccccccccccccccccccccccccc
c
c     Sortie 
c
c
c     Residu Residu 
c
       dimension ResU(NbdofMeca,NbDim)
c
c
c     Indice de contact par face frac
c
       dimension IndContact(NbFaceFrac)
c
c     SolLambda: on ne remplit que la composante normale ds le cas sans friction 
c
c     Moyenne sur la face de [sigmann(U) - beta [U]n ]_R-
c       
       dimension SolLambda(NbFaceFrac,NbDim)
c
c
c       
cccccccccccccccccccccccccccccc
c
c     workspaces 
c
       dimension XXk(NbDim)    
       
       dimension Vec1(NbDim), Vec2(NbDim)

       dimension Saut(NbDim), Sautnt(NbDim)  
c
c
c     numero global des inconnues de l'operateur
c     On numerote les inconnues locales dans l'ordre suivant:
c        inc maille k1       
c        nodes face sigma cote k2 puis bulle (si presente) cote k2        
c      
c      NbIncGlobalPThetaBeta
      
      dimension NumIncGlobalPThetaBeta(NbdofbyCellMax+NbNodebyFaceMax+1)      
c
c      
c
c     theta*Sigmanbyff - beta*EvalXSautbyff
c
c     On numerote les inconnues locales dans l'ordre suivant:
c        inc maille k1       
c        nodes face sigma cote k2 puis bulle (si presente) cote k2 
c      
      dimension EvalXPbeta(
     &     NbdofbyCellMax+NbNodebyFaceMax+1,NbDim,NbDim)

      
      dimension EvalXPthetabeta(
     &     NbdofbyCellMax+NbNodebyFaceMax+1,NbDim,NbDim)       
c
c     NbPoints = 2 x nb de nodes de la face 
c        
        dimension PoidsPoints(2*NbNodebyFaceMax)        
        dimension XPoints(2*NbNodebyFaceMax,NbDim)        
c
c     point de quadrature 
c
        dimension XXq(Nbdim)
c
c     PbetaUt
c
        dimension sPbetaUt(2)
c
c
c     
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        
        
cccccccccccccccccccccccccccccccccccccccc
c     Calcul du residu AU*SolU - SmU hors termes de contact 
cccccccccccccccccccccccccccccccccccccccc
c       
       do inc=1,NbdofMeca

          ResU(inc,:) = - SmU(inc,:)
          
          do m=IndLigneAU(inc)+1,IndLigneAU(inc+1)

             call ProdMatVec(AU(m,:,:),SolU(NColAU(m),:),Vec1)             

             ResU(inc,:) = ResU(inc,:) + Vec1(:)
             
          enddo


c          write(*,*)' res U AU0 U ',inc,ResU(inc,:)
          
       enddo
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       
c     terme de Nitsche int_Gamma [Pbeta(U)]_R- * Pthetabeta(V) 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
c       
c
c     Indice de contact par face frac: contact = 1 si contact en au moins un point de quadrature 
c       
       IndContact(:) = 0.d0 


       theta1 = 1.d0            ! 1.d0   

       
       Pst1 = 0.d0
       Pst2 = 0.d0
       sNormDerive = 0.d0

       TermTau1 = 0.d0
       TermTau2 = 0.d0

       SeuilSlip = 0.d0

       friction =  f_Frottement()

       b_coef = f_biot()

       alpha0 = f_alpha0()

       rK0 =  f_K0()

       Tref = f_TREF()
       

       
       
cccccccccccccccccccccccccccccc       
       do nff = 1,NbFaceFrac
cccccccccccccccccccccccccccccc
          
          inck = NumIncFaceFrac(nff)
          
          nf = NumFaceFracVersFace(nff)

          kcell = NumCellbyFace(nf,1)
          inccell = NumIncCell(kcell)
               
          incik1 = NumIncCellbyFaceFrac(nff,1)
c
c        on prend la valeur de  sur l'interface K1sigma 
c 
c          pinterface = - b_coef*SolP(incik1) + SolP(inck)  
c     pinterface = - b_coef*SolP(inccell) + SolP(inck)

         pTinterface = - b_coef*SolP(incik1) + SolP(inck)
     &     - alpha0 * rK0 * ( SolT(incik1) - Tref)  
                  

          betan = BetabyFaceFrac(nff,1)
          betat = BetabyFaceFrac(nff,2)

ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     normal stiffness
c
       unsurkn = f_UnSurkn()
       
c       
c 
       qn = 1.d0/( 1.d0 + betan*unsurkn )

          call QuadFace(
     &         nf,
     &         NbNode,NbFace,
     &         XNode,
     &         NbNodebyFace,NumNodebyFace,
     &         XFaceCG,SurfaceFace,
     &         NbPoints,XPoints,PoidsPoints) 

cccccccccccccccccccccccccccccccc
          do iq = 1,NbPoints  ! boucle sur les points de quadrature 
ccccccccccccccccccccccccccccccc

             
             XXq(:) = XPoints(iq,:)
             wq = PoidsPoints(iq)

c          do iq = 1,1   ! un seul pt de quadrature au CG                 
c             XXq(:) = XFaceCG(nf,:)
c             wq = SurfaceFace(nf)           
c

             theta0 = 0.d0 
             
             call EvalXPThetaBetabyFaceFrac(
     &            XXq,nff,theta0,betan,betat,
     &            NbFaceFrac,NbFace,NbCell,
     &            NbNodebyFace,NumNodebyFace,
     &            NbNodebyCell,NumNodebyCell,       
     &            NbIncGlobalbyCell,
     &            NumIncGlobalbyCell,
     &            NumNodeFaceLocalbyCell,
     &            NumFaceFracVersFace,
     &            NumFacebyCell,NumCellbyFace,
     &            IndBulleNumCellbyFace,            
     &            PoidsXFaceCG,XFaceCG,
     &            GradTangbyFace, 
     &            NbIncGlobalbyFaceFrac,
     &            NumIncGlobalbyFaceFrac,
     &            Sigmanbyff,
     &            VecNormalbyFace,
     &            VecTanFace1,
     &            VecTanFace2,                  
     &            EvalXPbeta,
     &            NbIncGlobalPThetaBeta,
     &            NumIncGlobalPThetaBeta)


             sPbetaUt(1) = 0.d0 
             do jk=1,NbIncGlobalPThetaBeta
                incj = NumIncGlobalPThetaBeta(jk)

                Vec1(:) = - SolU_nm1(incj,:)
                
                sPbetaUt(1) = sPbetaUt(1)
     &               + prodscal(EvalXPbeta(jk,:,2),Vec1)
             enddo
             
             
            sPbetaUt(2) = 0.d0 
             do jk=1,NbIncGlobalPThetaBeta
                incj = NumIncGlobalPThetaBeta(jk)
           
                Vec1(:) = - SolU_nm1(incj,:)                
                
                sPbetaUt(2) = sPbetaUt(2)
     &               + prodscal(EvalXPbeta(jk,:,3),Vec1)
             enddo
             
             call EvalXPThetaBetabyFaceFrac(
     &            XXq,nff,theta1,betan,betat,
     &            NbFaceFrac,NbFace,NbCell,
     &            NbNodebyFace,NumNodebyFace,
     &            NbNodebyCell,NumNodebyCell,       
     &            NbIncGlobalbyCell,
     &            NumIncGlobalbyCell,
     &            NumNodeFaceLocalbyCell,
     &            NumFaceFracVersFace,
     &            NumFacebyCell,NumCellbyFace,
     &            IndBulleNumCellbyFace,            
     &            PoidsXFaceCG,XFaceCG,
     &            GradTangbyFace, 
     &            NbIncGlobalbyFaceFrac,
     &            NumIncGlobalbyFaceFrac,
     &            Sigmanbyff,
     &            VecNormalbyFace,
     &            VecTanFace1,
     &            VecTanFace2,                  
     &            EvalXPbeta,
     &            NbIncGlobalPThetaBeta,
     &            NumIncGlobalPThetaBeta) 
          

             sPbetaUn = 0.d0 
             do jk=1,NbIncGlobalPThetaBeta
                incj = NumIncGlobalPThetaBeta(jk)
                sPbetaUn = sPbetaUn
     &               + prodscal(EvalXPbeta(jk,:,1),SolU(incj,:))
             enddo

             sPbetaUn = sPbetaUn + pTinterface

             


             do jk=1,NbIncGlobalPThetaBeta
                incj = NumIncGlobalPThetaBeta(jk)

                Vec1(:) = SolU(incj,:) 
                
                sPbetaUt(1) = sPbetaUt(1)
     &               + prodscal(EvalXPbeta(jk,:,2),Vec1)
             enddo
             
             
             do jk=1,NbIncGlobalPThetaBeta
                incj = NumIncGlobalPThetaBeta(jk)

                Vec1(:) = SolU(incj,:)       
                
                sPbetaUt(2) = sPbetaUt(2)
     &               + prodscal(EvalXPbeta(jk,:,3),Vec1)
             enddo


             sPbetat_norme = dsqrt(sPbetaUt(1)**2 + sPbetaUt(2)**2)

             SeuilSlip = -friction*qn*dmin1(sPbetaUn,0.d0)


                
                call EvalXPThetaBetabyFaceFrac(
     &               XXq,nff,theta,betan,betat,
     &               NbFaceFrac,NbFace,NbCell,
     &               NbNodebyFace,NumNodebyFace,
     &               NbNodebyCell,NumNodebyCell,       
     &               NbIncGlobalbyCell,
     &               NumIncGlobalbyCell,
     &               NumNodeFaceLocalbyCell,
     &               NumFaceFracVersFace,
     &               NumFacebyCell,NumCellbyFace,
     &               IndBulleNumCellbyFace,            
     &               PoidsXFaceCG,XFaceCG,
     &               GradTangbyFace, 
     &               NbIncGlobalbyFaceFrac,
     &               NumIncGlobalbyFaceFrac,
     &               Sigmanbyff,
     &               VecNormalbyFace,
     &               VecTanFace1,
     &               VecTanFace2,                     
     &               EvalXPthetabeta,
     &               NbIncGlobalPThetaBeta,
     &               NumIncGlobalPThetaBeta)

             
             
c
c             
c             
c             write(*,*)' sPbetaUn ',nff,sPbetaUn             
c
               if (sPbetaUn.le.0.d0) then
c               if (sPbetaUn.lt.0.d0) then
                   
                   if (sPbetat_norme.le.SeuilSlip) then

                      IndContact(nff) = 1 ! stick contact

                      Pst1 =  sPbetaUt(1)
                      Pst2 =  sPbetaUt(2)

                  else

                      IndContact(nff) = 2 !slip contact
                    
                      Pst1 =  SeuilSlip*sPbetaUt(1)/sPbetat_norme 

                      Pst2 =  SeuilSlip*sPbetaUt(2)/sPbetat_norme

                 endif



                do ik=1,NbIncGlobalPThetaBeta
                   ni = NumIncGlobalPThetaBeta(ik)

                   ResU(ni,:) = ResU(ni,:)
     &                  + wq*qn/betan
     &                  *sPbetaUn
     &                  *EvalXPthetabeta(ik,:,1)

                   ResU(ni,:) = ResU(ni,:)
     &                  + wq/betat
     &                  *Pst1
     &                  *EvalXPthetabeta(ik,:,2)


                   
                   ResU(ni,:) = ResU(ni,:)
     &                  + wq/betat
     &                  *Pst2
     &                  *EvalXPthetabeta(ik,:,3)



                   
                   do jk=1,NbIncGlobalPThetaBeta
                      nj = NumIncGlobalPThetaBeta(jk)
                     
                      mij = 0   ! recherche elt ni nj dans la ligne ni du CSR 
                      do m=IndLigneAU(ni)+1,IndLigneAU(ni+1)
                         n = NColAU(m)
                         if (n.eq.nj) then
                            mij = m 
                         endif
                      enddo
c     
                      
                      if (mij.eq.0) then
                         write(*,*)
     &                        ' elt ni nj non trouve dans CSR Nitsche ',
     &                        ni,nj
                         write(*,*)' ni nj ',ni,nj
                         write(*,*)' m1 m2 ',
     &                        IndLigneAU(ni),IndLigneAU(ni+1)
                         do m=IndLigneAU(ni)+1,IndLigneAU(ni+1)
                            n = NColAU(m)
                            write(*,*)' ni n  m ',ni,n,m
                         enddo                  
                         stop
                      else

                         do id = 1,Nbdim
                            do jd = 1,Nbdim
                               
                               AU(mij,id,jd) = AU(mij,id,jd)
     &                              + wq*qn/betan
     &                              *EvalXPbeta(jk,jd,1)
     &                              *EvalXPthetabeta(ik,id,1)


                               
c                               if (sPbetat_norme.le.SeuilSlip) then
                               if (IndContact(nff).eq.1) then 
                                
                                  AU(mij,id,jd) = AU(mij,id,jd)
     &                                 + wq/betat
     &                                 *EvalXPbeta(jk,jd,2)
     &                                 *EvalXPthetabeta(ik,id,2)
                                  

                                  AU(mij,id,jd) = AU(mij,id,jd)
     &                                 + wq/betat
     &                                 *EvalXPbeta(jk,jd,3)
     &                                 *EvalXPthetabeta(ik,id,3)
                                  
                                else
                            
                                  sNormDerive =
     &                              ( sPbetaUt(1)*EvalXPbeta(jk,jd,2)
     &                              + sPbetaUt(2)*EvalXPbeta(jk,jd,3) )
     &                                  /sPbetat_norme**3

                                  TermTau1=
     &                            ( EvalXPbeta(jk,jd,1)*sPbetaUt(1)
     &                             + sPbetaUn*EvalXPbeta(jk,jd,2) )
     &                             /sPbetat_norme
     &                             - sPbetaUn*sPbetaUt(1)*sNormDerive

                                  TermTau2=
     &                            ( EvalXPbeta(jk,jd,1)*sPbetaUt(2)
     &                             + sPbetaUn*EvalXPbeta(jk,jd,3) )
     &                             /sPbetat_norme
     &                             - sPbetaUn*sPbetaUt(2)*sNormDerive
                                                                 

                                  AU(mij,id,jd) = AU(mij,id,jd)
     &                               - friction*wq*qn/betat
     &                              *TermTau1
     &                              *EvalXPthetabeta(ik,id,2)
                                  
                                  AU(mij,id,jd) = AU(mij,id,jd)
     &                               -friction*wq*qn/betat
     &                               *TermTau2                           
     &                               *EvalXPthetabeta(ik,id,3)
                                  
                                   
                                endif 

                            enddo
                         enddo
                         
                      endif
                      

                   enddo
                   
                enddo                
                
                
             endif






c             write(*,*) 'friction', friction

c             stop

             
ccccccccccccccccccccccccccc
          enddo ! fin boucle points de quadrature 
ccccccccccccccccccccccccccc



c      
cccccccccccccccccccccccccccccccccccc          
c
c     calcul de rlambdan et IndContact au CG 
c          
          XXq(:) = XFaceCG(nf,:)


             theta0 = 0.d0 
             
             call EvalXPThetaBetabyFaceFrac(
     &            XXq,nff,theta0,betan,betat,
     &            NbFaceFrac,NbFace,NbCell,
     &            NbNodebyFace,NumNodebyFace,
     &            NbNodebyCell,NumNodebyCell,       
     &            NbIncGlobalbyCell,
     &            NumIncGlobalbyCell,
     &            NumNodeFaceLocalbyCell,
     &            NumFaceFracVersFace,
     &            NumFacebyCell,NumCellbyFace,
     &            IndBulleNumCellbyFace,            
     &            PoidsXFaceCG,XFaceCG,
     &            GradTangbyFace, 
     &            NbIncGlobalbyFaceFrac,
     &            NumIncGlobalbyFaceFrac,
     &            Sigmanbyff,
     &            VecNormalbyFace,
     &            VecTanFace1,
     &            VecTanFace2,                  
     &            EvalXPbeta,
     &            NbIncGlobalPThetaBeta,
     &            NumIncGlobalPThetaBeta)


             sPbetaUt(1) = 0.d0 
             do jk=1,NbIncGlobalPThetaBeta
                incj = NumIncGlobalPThetaBeta(jk)

                Vec1(:) = - SolU_nm1(incj,:)
                
                sPbetaUt(1) = sPbetaUt(1)
     &               + prodscal(EvalXPbeta(jk,:,2),Vec1)
             enddo
             
             
            sPbetaUt(2) = 0.d0 
             do jk=1,NbIncGlobalPThetaBeta
                incj = NumIncGlobalPThetaBeta(jk)
           
                Vec1(:) = - SolU_nm1(incj,:)                
                
                sPbetaUt(2) = sPbetaUt(2)
     &               + prodscal(EvalXPbeta(jk,:,3),Vec1)
             enddo

          
             call EvalXPThetaBetabyFaceFrac(
     &            XXq,nff,theta1,betan,betat,
     &            NbFaceFrac,NbFace,NbCell,
     &            NbNodebyFace,NumNodebyFace,
     &            NbNodebyCell,NumNodebyCell,       
     &            NbIncGlobalbyCell,
     &            NumIncGlobalbyCell,
     &            NumNodeFaceLocalbyCell,
     &            NumFaceFracVersFace,
     &            NumFacebyCell,NumCellbyFace,
     &            IndBulleNumCellbyFace,            
     &            PoidsXFaceCG,XFaceCG,
     &            GradTangbyFace, 
     &            NbIncGlobalbyFaceFrac,
     &            NumIncGlobalbyFaceFrac,
     &            Sigmanbyff,
     &            VecNormalbyFace,
     &            VecTanFace1,
     &            VecTanFace2,              
     &            EvalXPbeta,
     &            NbIncGlobalPThetaBeta,
     &            NumIncGlobalPThetaBeta) 
          

             sPbetaUn = 0.d0 
             do jk=1,NbIncGlobalPThetaBeta
                incj = NumIncGlobalPThetaBeta(jk)
                sPbetaUn = sPbetaUn
     &               + prodscal(EvalXPbeta(jk,:,1),SolU(incj,:))
             enddo

             
             sPbetaUn = sPbetaUn + pTinterface

             
             do jk=1,NbIncGlobalPThetaBeta
                incj = NumIncGlobalPThetaBeta(jk)

                Vec1(:) = SolU(incj,:) 
                
                sPbetaUt(1) = sPbetaUt(1)
     &               + prodscal(EvalXPbeta(jk,:,2),Vec1)
             enddo             

             do jk=1,NbIncGlobalPThetaBeta
                incj = NumIncGlobalPThetaBeta(jk)

                Vec1(:) = SolU(incj,:) 
                
                sPbetaUt(2) = sPbetaUt(2)
     &               + prodscal(EvalXPbeta(jk,:,3),Vec1)
             enddo

             
             rlambdan = -qn*dmin1(0.d0,sPbetaUn)

             sautn = -(qn*rlambdan + sPbetaUn)/(betan*qn)

             SeuilSlip = friction*rlambdan

             

             sPbetat_norme = dsqrt(sPbetaUt(1)**2 + sPbetaUt(2)**2)

c             write(*,*)' contact ',nff,pinterface,
c     &          sPbetat_norme-SeuilSlip

             

             if(sPbetat_norme.le.SeuilSlip) then

                Pst1=  sPbetaUt(1)
                Pst2=  sPbetaUt(2)

             else
                
                Pst1= seuilSlip*sPbetaUt(1)/sPbetat_norme

                Pst2= seuilSlip*sPbetaUt(2)/sPbetat_norme
               
             endif

             dsautt1 = -(Pst1 - sPbetaUt(1))/betat 
             dsautt2 = -(Pst2 - sPbetaUt(2))/betat

             
             SolLambda(nff,1) = rlambdan
             SolLambda(nff,2) = - Pst1
             SolLambda(nff,3) = - Pst2 
c             
c     Indice contact au CG 
c             
             IndContact(nff) = 0             
             if (sPbetaUn.le.0.d0) then                 
                if (sPbetat_norme.le.SeuilSlip) then
                   IndContact(nff) = 1 ! stick contact
                else
                   IndContact(nff) = 2 !slip contact
                endif                   
             endif

c             if (IndContact(nff).ne.1) then 
c                write(*,*)' saut ',nff,sautn,dsautt1,dsautt2,
c     &               IndContact(nff)
c             endif
             
c          write(*,*)' rlambdan ',nff,rlambdan

c             write(*,*)' ddsautt ',nff,dsautt1,dsautt2,IndContact(nff)

             
cccccccccccccccccccccccccccccc
       enddo  ! fin boucle faces frac 
cccccccccccccccccccccccccccccc
c
c
c       
c       
ccccccccccccccccccccccccccccccccccccccccccc      
c
c     CL de Dirichlet aux bord, correction AU et ResU
c
cccccccccccccccccccccccccccccccccccccccccc
c      
      do k = 1,NbCell

         XXk(:) = XCellCG(k,:)
         
         do i=1,NbNodebyCell(k)
            inc = NumIncGlobalbyCell(k,i)
            
            n = NumNodebyCell(k,i) ! numero global du node
            
            if (IndDirNodeMeca(n).eq.5) then

               Vec1(3) = 0.d0               
               ResU(inc,3) = SolU(inc,3) - Vec1(3) 
               
               do m = IndLigneAU(inc)+1,IndLigneAU(inc+1)
                  AU(m,3,:) = 0.d0 
               enddo
               do id = 3,3
                  AU(IndLigneAU(inc)+1,id,id) = 1.d0 ! mise a Id de la ligne inc
               enddo
            endif


            if (IndDirNodeMeca(n).eq.3) then

               Vec1(2) = 0.d0               
               ResU(inc,2) = SolU(inc,2) - Vec1(2) 
               
               do m = IndLigneAU(inc)+1,IndLigneAU(inc+1)
                  AU(m,2,:) = 0.d0 
               enddo
               do id = 2,2
                  AU(IndLigneAU(inc)+1,id,id) = 1.d0 ! mise a Id de la ligne inc
               enddo
            endif            


            if (IndDirNodeMeca(n).eq.4) then

               Vec1(2) = 0.d0               
               ResU(inc,2) = SolU(inc,2) - Vec1(2) 
               
               do m = IndLigneAU(inc)+1,IndLigneAU(inc+1)
                  AU(m,2,:) = 0.d0 
               enddo
               do id = 2,2
                  AU(IndLigneAU(inc)+1,id,id) = 1.d0 ! mise a Id de la ligne inc
               enddo
            endif


            if (IndDirNodeMeca(n).eq.1) then

               Vec1(1) = 0.d0               
               ResU(inc,1) = SolU(inc,1) - Vec1(1) 
               
               do m = IndLigneAU(inc)+1,IndLigneAU(inc+1)
                  AU(m,1,:) = 0.d0 
               enddo
               do id = 1,1
                  AU(IndLigneAU(inc)+1,id,id) = 1.d0 ! mise a Id de la ligne inc
               enddo
            endif            

            
            
         enddo       
                           
      enddo
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccc       
       
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       
