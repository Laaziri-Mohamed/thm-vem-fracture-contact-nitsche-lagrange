
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     Modele d'elasticite lineaire isotrope avec fractures + contact avec frottement 
c
c     Champ u discretise avec VEM de degre 1 + Bulles aux faces (frac a priori) 
c
c     Formulation mixte avec multiplicateurs constants par face frac 
c
c
c
c     Elasticite + contact avec frottement (coefficient de frottement constant par face)
c
c     prise en compte de la NORMAL STIFFNESS (cas avec penetration)  f_kn()
c
c     avec la fonction f definie par
c
c                     f(s) = 0 si s <=0, et f(s) = kn*s sinon  
c
c     Prise en compte de la Cohesion f_cohesion()
c      
c
c     Calcul du Residu Contact VEM Bulle avec fractures 
c
c     partie lineaire sur i=1,NbdofMeca calculee par produit matrice vecteur 
c
c     eq du contact: [u]n = (uKsigma - uLsigma). nKsigma <= 0 
c
c                    rlambdan = LambdaKsigma_n >= 0
c
c                    rlambdan = f([u]n)   
c
c                    rlambda_t = (rlambdat1,rlambdat2) ds un repere orthonorme sur la face frac
c
c                    |rlambdat| <= F rlambdan + Cohesion 
c
c                    rlambdat.[u]t = (F rlambdan + Cohesion) |[u]t| 
c
c      
c     semi-smooth Newton sur la fonction (avec beta parametre > 0) 
c
c     rlambdan - f([u]n) = 0
c
c
c      <=>   rlambdan = qn [ bn ]+   avec qn = 1/(1 + beta/kn)  et bn = rlambdan + beta [u]n 
c
c      
c     rlambdat = (rlambdat1,rlambdat2) 
c
c     rlambdat = [rlambdat + beta*[u]t ]_{F rlamndan  + Cohesion}
c
c      
c     Se simplifie en une methode de type Active Set
c      avec (eps=1.d-10 pour eviter une division par zero en mode slip) 
c
c     1. si rlamndan + beta* [u]n <= 0  : NO CONTACT  
c      
c
c        rlambdan = 0
c        rlambdat = 0
c
c     2. sinon on pose qn = kn/(kn+beta) = 1/( 1 + beta*1/kn )
c
c          rlambdan - kn*[u]n = 0  ou encore    beta*1/kn*rlambdan - beta [u]n = 0 
c
c        2.1 si |rlambdat + beta* [u]t| <= F *( qn*( rlamndan + beta* [u]n ) + Cohesion)  CONTACT STICK 
c
c          - beta*[u]t = 0 
c
c
c        2.2 sinon CONTACT SLIP 
c
c
c          rlambdat -  F (rlamndat + beta* [u]t)/|rlamndat + beta* [u]t|*( qn( rlamndan + beta* [u]n ) + Cohesion) =0
c
ccccccc
c
c     On retrouve la condition de non penetration pour unsurkn = 0. 
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c      
      subroutine ElasticiteVEMBulleFrac(Temps, 
     &     NbNode,NbCell,NbFace,NbFaceFrac,
     &     NbdofMeca,NbdofMecaContact,
     &     rF, 
     &     IndDirNodeMeca,IndFace,          
     &     NbNodebyCell,NumNodebyCell,
     &     NbNodebyFace,NumNodebyFace,
     &     NbFacebyCell,NumFacebyCell,
     &     NumCellbyFace,NumFaceFracVersFace,NumFaceVersFaceFrac,      
     &     VecNormalbyFace,PoidsXCellCG,PoidsXFaceCG,      
     &     XNode,XCellCG,XFaceCG,SurfaceFace,VolCell,      
     &     NbIncGlobalbyCell,NumIncGlobalbyCell,
     &     NumNodeFaceLocalbyCell,
     &     NbIncGlobalbyFaceFrac,NumIncGlobalbyFaceFrac,
     &     ACell,SautbyFaceFrac,
     &     VecTanFace1,VecTanFace2, 
     &     SolU,SolU_nm1,IndContact,nit,
     &     SolT, 
     &     SolP,NbCV,NumIncFaceFrac,NumIncCell,
     &     GradCellVEM)      
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

c     pressions dans la matrice et dans la fracture

        dimension SolP(NbCV)


c     Temperature dans la matrice 

        dimension SolT(NbCV)        
c

c     num cell       >  num inc cell
c     num face frac  >  num inc face frac  

        dimension NumIncFaceFrac(NbFaceFrac)
        dimension NumIncCell(NbCell)

c        
c     Mailles par les nodes 
c
        dimension NbNodebyCell(NbCell)
        dimension NumNodebyCell(NbCell,NbNodebyCellMax)
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
c     numero de la face frac pour une face donnee (0 si non face frac) 
c        
        dimension NumFaceVersFaceFrac(NbFace)        
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
c     poids barycentrique aux faces 
c        
       dimension PoidsXFaceCG(NbFace,NbNodebyFaceMax)
        

cc   vecteur gradiant pour chaque fonction de bases scalaire v^dof (dof=1,...,NbdofbyCell(K)) pour chaque cellule K

      dimension GradCellVEM(NbCell,NbdofbyCellMax,NbDim) 
c
c      
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
c     Operateur saut par face frac 
c
        dimension NbIncGlobalbyFaceFrac(NbFaceFrac)
        dimension NumIncGlobalbyFaceFrac(NbFaceFrac,2*NbNodebyFaceMax+2)           
        dimension SautbyFaceFrac(NbFaceFrac,2*NbNodebyFaceMax+2)

c
c     Matrice locale par maille VEM 
c
       dimension ACell(NbCell,NbdofbyCellMax,NbdofbyCellMax,NbDim,NbDim)       
c
c
cccccccccccccccccccccccccccccccc
c
c     coefficient de frottement par face frac 
c
       dimension rF(NbFaceFrac)
c
c
c       
ccccccccccccccccccccccccccccccccc
c
c     Sorties 
c
c
c     Solution U (init en entree et calcul en sortie) 
c
       dimension SolU(NbdofMecaContact,NbDim)
c
c     Solution U a tnm1 
c
       dimension SolU_nm1(NbdofMecaContact,NbDim)       
c
c
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
c     Jacobienne VEM 
c
c
c    structure creuse mecanique 
c
c        
       integer, dimension(:), allocatable :: IndLigneAU
       integer, dimension(:), allocatable :: NLigneAU
       integer, dimension(:), allocatable :: NColAU
       
       double precision, dimension(:,:,:), allocatable :: AU 
       double precision, dimension(:,:,:), allocatable :: AU22 

       integer, dimension(:), allocatable :: NbNzbyLineMeca
       integer, dimension(:,:), allocatable :: NumNzbyLineMeca
       
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
       dimension Vec1(NbDim),Vec2(NbDim),Vec3(NbDim),Xf(NbDim)
c
c     produit scalaire (ei,fj) ou ei base canonique, fj base locale (vecnormal,tau1,tau2)
c       
       dimension BaseFace(NbDim,NbDim) 

c     Vecteur unitaire       
      dimension vecUnit(NbDim)
c
      dimension vecbtau(NbDim)

      dimension Saut(NbDim), Sautnt(NbDim)  

      dimension bnfrac(NbDim)
c
c              
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
c
c
       beta = f_beta_semi_smooth() ! parameter du semi-smooth Newton pour le contact 
c
c       
ccccccccccccccccccccccccccccccccccccccccccc


       
      
cccccccccccccccccccccccccccccccccccccccc
c
c     Structure creuse de la Jacobienne Meca 
c
cccccccccccccccccccccccccccccccccccccccc
c
c       
      allocate(NbNzbyLineMeca(NbdofMecaContact))
      allocate(NumNzbyLineMeca(NbdofMecaContact,NbNzbyLineMecaMax))      
      
c
c     init avec la diagonale first 
c      
      do k=1,NbdofMecaContact 
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

               if (indij.eq.-1) then ! on rajoute l'elt ij 
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


      do ifrac = 1,NbFaceFrac

         incfrac = NbdofMeca + ifrac 
         
         do ik=1,NbIncGlobalbyFaceFrac(ifrac)
            inc = NumIncGlobalbyFaceFrac(ifrac,ik)
c
c     element (inc,incfrac)
c            
            indij = -1          ! test si (inc,incfrac) deja existant ds la ligne inc 
            do l=1, NbNzbyLineMeca(inc)
               n = NumNzbyLineMeca(inc,l)
               if (n.eq.incfrac) then
                  indij = 1 
               endif
            enddo            

            if (indij.eq.-1) then ! on rajoute l'elt (inc,incfrac) 
               NbNzbyLineMeca(inc) = NbNzbyLineMeca(inc)+1
               if ( NbNzbyLineMeca(inc).gt.NbNzbyLineMecaMax) then
                  write(*,*)' increase NbNzbyLineMecaMax ',
     &                 NbNzbyLineMecaMax,NbNzbyLineMeca(inc)
                  stop
               endif               
               NumNzbyLineMeca(inc,NbNzbyLineMeca(inc)) = incfrac                  
            endif            
c
c     element (incfrac,inc)
c
c            
            indij = -1          ! test si (incfrac,inc) deja existant ds la ligne incfrac 
            do l=1, NbNzbyLineMeca(incfrac)
               n = NumNzbyLineMeca(incfrac,l)
               if (n.eq.inc) then
                  indij = 1 
               endif
            enddo            

            if (indij.eq.-1) then ! on rajoute l'elt (incfrac,inc)             
               NbNzbyLineMeca(incfrac) = NbNzbyLineMeca(incfrac)+1
               if ( NbNzbyLineMeca(incfrac).gt.NbNzbyLineMecaMax) then
                  write(*,*)' increase NbNzbyLineMecaMax ',
     &                 NbNzbyLineMecaMax,NbNzbyLineMeca(incfrac)
                  stop
               endif               
               NumNzbyLineMeca(incfrac,NbNzbyLineMeca(incfrac)) = inc      
            endif
            
         enddo         
c
c
c
c         
      enddo
      
      allocate(IndLigneAU(NbdofMecaContact+1))
      IndLigneAU(1) = 0
      do i=1,NbdofMecaContact
         IndLigneAU(i+1) = IndLigneAU(i) + NbNzbyLineMeca(i)          
      enddo
      MemBlocU = IndLigneAU(NbdofMecaContact+1)


c      write(*,*)' MemBlocU ',MemBlocU 
      
      allocate(NLigneAU(MemBlocU))
      allocate(NColAU(MemBlocU))
      

      do i=1,NbdofMecaContact
         do m=IndLigneAU(i)+1,IndLigneAU(i+1)
            n = m-IndLigneAU(i)
            NColAU(m) = NumNzbyLineMeca(i,n)
            NLigneAU(m) = i
         enddo
      enddo
         

      deallocate(NbNzbyLineMeca)
      deallocate(NumNzbyLineMeca)      

c      do i=1,NbdofMecaContact 
c         write(*,*)' ligne i ',i
c         do m=IndLigneAU(i)+1,IndLigneAU(i+1)
c            write(*,*)' m i j ',m,NLigneAU(m),NColAU(m)
c         enddo
c         write(*,*)
c     enddo
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       
cccccccccccccccccccccccccccccccccccccccccccc      
c
c     Assemblage Jacobienne VEM \sum_K  aK  
c
cccccccccccccccccccccccccccccccccccccccccccc

      
      allocate(AU(MemBlocU,NbDim,NbDim))

      AU = 0.d0 ! on met AU a zero 
              
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
c                     write(*,*)' n nj m ',n,nj,m
                  enddo                  
                  stop
               else
                 
                  AU(mij,:,:) = AU(mij,:,:) + ACell(k,i,j,:,:)
                                    
               endif

            
            enddo
         enddo
         
ccccccccccccc         
      enddo  ! fin boucle mailles 
ccccccccccccc
c

      
ccccccccccccccccc
c      
c      terme   sum_sig lambda_sig (v+ - v-)_sig  
c
      do ifrac = 1,NbFaceFrac
         nf = NumFaceFracVersFace(ifrac)

         do id=1,NbDim
            BaseFace(id,2) = VecTanFace1(nf,id)
            BaseFace(id,3) = VecTanFace2(nf,id)
            BaseFace(id,1) = VecNormalbyFace(nf,id)            
         enddo
      
         incfrac = NbdofMeca + ifrac 
         
         do ik=1,NbIncGlobalbyFaceFrac(ifrac)
            inc = NumIncGlobalbyFaceFrac(ifrac,ik)

            
            mij = 0             ! recherche elt (inc,incfrac) dans la ligne inc du CSR 
            do m=IndLigneAU(inc)+1,IndLigneAU(inc+1)
               n = NColAU(m)
               if (n.eq.incfrac) then
                  mij = m 
               endif
            enddo
               
            if (mij.eq.0) then
               write(*,*)' elt inc,incfrac non trouve dans CSR ',
     &                  ifrac,inc
               stop
            else


               do id=1,NbDim
                  do jd=1,NbDim
                     AU(mij,id,jd) = AU(mij,id,jd)
     &                    + SautbyFaceFrac(ifrac,ik)*SurfaceFace(nf)
     &                      *BaseFace(id,jd)          
                  enddo
               enddo
               
            endif

         
         enddo         
         
      enddo
c      
c
ccccccccccccccccccccccccccccccccccccccccccc      
c
c     CL de Dirichlet aux bord -> mise a Id des lignes sur la Jacobienne 
c
cccccccccccccccccccccccccccccccccccccccccc
c      
      do k = 1,NbCell
         do i=1,NbNodebyCell(k)
            inc = NumIncGlobalbyCell(k,i)
            
            n = NumNodebyCell(k,i) ! numero global du node
            
            if (IndDirNodeMeca(n).eq.5) then               
               do m = IndLigneAU(inc)+1,IndLigneAU(inc+1)
                  AU(m,3,:) = 0.d0 
               enddo
              
               do id = 3,3
                  AU(IndLigneAU(inc)+1,id,id) = 1.d0 ! mise a Id de la ligne inc
               enddo
               
            endif

            if (IndDirNodeMeca(n).eq.3) then               
               do m = IndLigneAU(inc)+1,IndLigneAU(inc+1)
                  AU(m,2,:) = 0.d0 
               enddo
              
               do id = 2,2
                  AU(IndLigneAU(inc)+1,id,id) = 1.d0 ! mise a Id de la ligne inc
               enddo
               
            endif

            if (IndDirNodeMeca(n).eq.4) then               
               do m = IndLigneAU(inc)+1,IndLigneAU(inc+1)
                  AU(m,2,:) = 0.d0 
               enddo
              
               do id = 2,2
                  AU(IndLigneAU(inc)+1,id,id) = 1.d0 ! mise a Id de la ligne inc
               enddo
               
            endif
            

            if (IndDirNodeMeca(n).eq.1) then               
               do m = IndLigneAU(inc)+1,IndLigneAU(inc+1)
                  AU(m,1,:) = 0.d0 
               enddo
              
               do id = 1,1
                  AU(IndLigneAU(inc)+1,id,id) = 1.d0 ! mise a Id de la ligne inc
               enddo
               
            endif                  


c$$$            if (IndDirNodeMeca(n).eq.6) then               
c$$$c               write(*,*)' node Dir ',inc,n,XNode(n,:) 
c$$$               do m = IndLigneAU(inc)+1,IndLigneAU(inc+1)
c$$$                  AU(m,1:3,:) = 0.d0 
c$$$               enddo
c$$$              
c$$$               do id = 1,3
c$$$                  AU(IndLigneAU(inc)+1,id,id) = 1.d0 ! mise a Id de la ligne inc
c$$$               enddo
c$$$               
c$$$            endif            


         enddo
         
         do i=NbNodebyCell(k)+1,NbIncGlobalbyCell(k)
            inc = NumIncGlobalbyCell(k,i)           

            nk = NumNodeFaceLocalbyCell(k,i) ! numero local de la face 
            n = NumFacebyCell(k,nk)            
            if (IndFace(n).ne.0) then

               write(*,*)' face Dir non permis ',inc,n,XFaceCG(n,:)
               stop

               
c               write(*,*)' face Dir ',inc,n,XFaceCG(n,:)                
               do m = IndLigneAU(inc)+1,IndLigneAU(inc+1)
                  AU(m,:,:) = 0.d0 
               enddo
               do id = 1,NbDim              
                  AU(IndLigneAU(inc)+1,id,id) = 1.d0 ! mise a Id de la ligne inc
               enddo
            endif            
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
      allocate(SmU(NbdofMecaContact,NbDim))

      
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






cc  Calcul du terme (dans le second membre) \sum_{sigma} \int_{sigma} pf(sigma)*saut(v)_n

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
c
cc  calcul du terme (dans le second membre) \sum_k int_{K} b*pm div(v) 
c
c

      b_coef = f_biot()
      
      alpha0 = f_alpha0()

      rK0 =  f_K0()

      Tref = f_TREF()

      do k=1,NbCell

         inck = NumIncCell(k)
         
         bpk =  b_coef*SolP(inck) + alpha0 * rK0 * ( SolT(inck) - Tref) 

cccccc  TMP on prend la p d'interface si maille adjacente a une face frac cccccc         
c         s = 0.d0
c         ss = 0.d0
c         ind = 0
c         do ik = 1,NbFacebyCell(k)
c            nfk = NumFacebyCell(k,ik)
c            nfrac = NumFaceVersFaceFrac(nfk)
c            if(nfrac.ge.1) then
c               incik = NumIncFaceFrac(nfrac)
c               s = s + SolP(incik)
c               ss = ss + 1.d0
c               ind = 1
c            endif
c         enddo
c         if (ind.eq.1) then
c            bpk = b_coef*s/ss
c         endif
cccccccccccccccccccccccccccccccccccc         
         

         do ik=1,NbIncGlobalbyCell(k)

            inc = NumIncGlobalbyCell(k,ik)


            SmU(inc,:) = SmU(inc,:)
     &           + VolCell(k)*bpk*GradCellVEM(k,ik,:)

         enddo
         
      enddo

      
c
c
c     Rajout des termes sources de type Neumann non homogene 
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
ccccccccccccccccccccccccccccccccccc      
c
c     Correction du SmU par les CL de Dirichlet au bord 
c
cccccccccccccccccccccccccccccccccc
c

c      TU = 10.d0 

      
      do k = 1,NbCell

         XXk(:) = XCellCG(k,:)

         
         do i=1,NbNodebyCell(k)
            inc = NumIncGlobalbyCell(k,i)
            
            n = NumNodebyCell(k,i) ! numero global du node

            if (IndDirNodeMeca(n).eq.5) then
               SmU(inc,3) = 0.d0                       
            endif
            if (IndDirNodeMeca(n).eq.3) then
               SmU(inc,2) = 0.d0                       
            endif
            if (IndDirNodeMeca(n).eq.4) then
               SmU(inc,2) = 0.d0                       
            endif            
            if (IndDirNodeMeca(n).eq.1) then
               SmU(inc,1) = 0.d0                       
            endif                   

c$$$            if (IndDirNodeMeca(n).eq.6) then
c$$$
c$$$               if (Temps.le.TU) then 
c$$$   
c$$$                  
c$$$                  SmU(inc,1) = 4.d0*Temps/TU
c$$$                  SmU(inc,2) = 4.d0*Temps/TU
c$$$                  SmU(inc,3) = - 4.d0*Temps/TU
c$$$                  
c$$$            else
c$$$                 
c$$$
c$$$                  SmU(inc,1) = 4.d0
c$$$                  SmU(inc,2) = 4.d0
c$$$                  SmU(inc,3) = - 4.d0
c$$$                  
c$$$               endif
c$$$            endif

            
         enddo

         do i=NbNodebyCell(k)+1,NbIncGlobalbyCell(k)
            inc = NumIncGlobalbyCell(k,i)           

            nk = NumNodeFaceLocalbyCell(k,i) ! numero local de la face 
            n = NumFacebyCell(k,nk)            
            if (IndFace(n).ne.0) then

               write(*,*)' bulle sur le bord '
               stop

               SmU(inc,:) = 0.d0 ! ne sert pas car pas de bulle de bord 
               
            endif            
         enddo            
                           
      enddo
c
c      
c      
ccccccccccccccccccccccccccccccc
c
c     Residu initial et init IndContact 
c
ccccccccccccccccccccccccccccccc
c      
      allocate(ResU(NbdofMecaContact,NbDim))

      ResU = 0.d0 
c
c     calcul du residu init ResU et de IndContact init 
c      
      call ResiduContact(
     &     beta,rF, 
     &     NbFaceFrac,NbFace,
     &     NbdofMeca,NbdofMecaContact,MemBlocU,  
     &     NbIncGlobalbyFaceFrac,NumIncGlobalbyFaceFrac,
     &     SautbyFaceFrac,    
     &     XFaceCG,NumFaceFracVersFace,
     &     VecNormalbyFace,
     &     VecTanFace1,VecTanFace2,          
     &     IndLigneAU,NLigneAU,NColAU,AU,       
     &     SolU,SolU_nm1,SmU,
     &     IndContact,ResU)

c
c     TMP 2D ccccccccccccc
c     
      ResU(:,2) = 0.d0 
c
cccccccccccccccccccccccccc
      
c
c     norme l2 du residu init 
c
      s = 0.d0 
      do inc=1,NbdofMecaContact 
         s = s + prodscal(ResU(inc,:),ResU(inc,:))
      enddo
      resinit = dsqrt(s) 
      
      residu_relatif = 1.d0 
      
      write(*,*)' res init ',resinit 
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Boucle semi-smooth Newton 
c

      nit = 0
      nitmax = 40
      critere_arret_newton = 1.0E-9
      critere_arret_dUmax = 1.0E-8      
      
      dUmax = 1.d0 
      
      allocate(dSolU(NbdofMecaContact,NbDim))
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c      do while (
c     &         (residu_relatif.ge.critere_arret_newton)
c     &    .and.(dUmax.ge.critere_arret_dUmax)         
c     &     .and.(nit.le.nitmax-1)
c     &     )


      do while (
c     &         (residu_relatif.ge.critere_arret_newton)
     &    (dUmax.ge.critere_arret_dUmax)         
     &     .and.(nit.le.nitmax-1)
     &        )         
c         
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         nit = nit + 1
c
c     on rajoute les termes lies aux equations du contact dans la Jacobienne 
c     Il faut SolU et IndContact en entree 
c         
         call JACContact(
     &        beta,rF, 
     &        NbFaceFrac,NbFace,
     &        NbdofMeca,NbdofMecaContact,MemBlocU,  
     &        NbIncGlobalbyFaceFrac,NumIncGlobalbyFaceFrac,
     &        SautbyFaceFrac,
     &        NumFaceFracVersFace,
     &        VecNormalbyFace,
     &        VecTanFace1,VecTanFace2,      
     &        IndLigneAU,NColAU,AU,SolU,SolU_nm1,       
     &        IndContact)      
c
c
c         
cccccccccccccccccccccccccccccccccc
c     On resoud le systeme lineaire avec -ResU en second membre -> dSolU 
ccccccccccccccccccccccccccccccccc
c      
      critere_arret_gmres = 1.0E-10

      
      ResU = - ResU

      
c      NbIncMeca = NbDim
c      IPREC = 0  ! ILU0 
c      MAXL = 1000
c      call solveur_iteratif(
c     &       NbdofMecaContact,MemBlocU,NbIncMeca,
c     &       AU,NLigneAU,NColAU,IndLigneAU,
c     &       ResU,dSolU,critere_arret_gmres,MAXL,IPREC,
c     &       ITER,IERR,ERR)      
c
c
cccccccccccccccccccccccccccccccc      
c
c     reduction du systeme (2,2) pour cas 2D -> avec eq et inc 1 et 3 ici  
c      
      allocate(ResU2(NbdofMecaContact,2))
      allocate(dSolU2(NbdofMecaContact,2))

      ResU2(:,1) = ResU(:,1)
      ResU2(:,2) = ResU(:,3)

      allocate(AU22(MemblocU,2,2))

      AU22(:,1,1) = AU(:,1,1)
      AU22(:,1,2) = AU(:,1,3)
      AU22(:,2,1) = AU(:,3,1)
      AU22(:,2,2) = AU(:,3,3)

      ITER = 0
      NbIncMeca = NbDim - 1
      call solveurSuperLU(NbdofMecaContact,MemBlocU,NbIncMeca,
     &     AU22,NLigneAU,NColAU,IndLigneAU,
     &     ResU2,dSolU2)

      dSolU(:,1) = dSolU2(:,1)
      dSolU(:,2) = 0.d0 
      dSolU(:,3) = dSolU2(:,2)

      deallocate(ResU2,dSolU2,AU22)
c
cccccccccccccccccccccccccccccccccccc
c      
c      
c      write(*,*)' debut solveur SuperLU '
c$$$      ITER = 0
c$$$      NbIncMeca = NbDim 
c$$$      call solveurSuperLU(NbdofMecaContact,MemBlocU,NbIncMeca,
c$$$     &     AU,NLigneAU,NColAU,IndLigneAU,
c$$$     &     ResU,dSolU)
c
c      
c      
c      call solveurDirect(NbdofMecaContact,MemBlocU,NbIncMeca,
c     &     AU,NLigneAU,NColAU,IndLigneAU,
c     &     ResU,dSolU)      

c      write(*,*)' fin solveur SuperLU '      
c
c     Increment max de U 
c      
      dUmax = 0.d0 
      do k=1,NbdofMeca
         dUmax = dmax1(dUmax,dabs(dSolU(k,1)))
         dUmax = dmax1(dUmax,dabs(dSolU(k,2)))
         dUmax = dmax1(dUmax,dabs(dSolU(k,3)))         
      enddo
c
c
c      
cccccccccccccccccccccccccccccccc      
c
c     incrementation de Newton      
c      
      SolU = SolU + dSolU
c
ccccccccccccccccccccccccccccccc      
c
c     calcul du residu ResU et de IndContact 
c
      call ResiduContact(
     &     beta,rF, 
     &     NbFaceFrac,NbFace,
     &     NbdofMeca,NbdofMecaContact,MemBlocU,  
     &     NbIncGlobalbyFaceFrac,NumIncGlobalbyFaceFrac,
     &     SautbyFaceFrac,      
     &     XFaceCG,NumFaceFracVersFace,
     &     VecNormalbyFace,
     &     VecTanFace1,VecTanFace2,         
     &     IndLigneAU,NLigneAU,NColAU,AU,       
     &     SolU,SolU_nm1,SmU,
     &     IndContact,ResU)    
c
c     TMP 2D ccccccccccccc
c     
      ResU(:,2) = 0.d0 
c
cccccccccccccccccccccccccc
c      
c     norme l2 du residu
c
      s = 0.d0 
      do inc=1,NbdofMecaContact 
         s = s + prodscal(ResU(inc,:),ResU(inc,:))
c         write(*,*)' res ',inc,ResU(inc,:)
      enddo
      res = dsqrt(s) 

      residu_relatif = res/resinit
      
c      write(*,*)
      write(*,*)' it newton Meca ',nit,res,resinit,residu_relatif,dUmax
c      write(*,*)      
      
ccccccccccccccccccccccccccccccccccccc
      enddo ! fin boucle while  Newton 
cccccccccccccccccccccccccccccccccccc

      
cccccccccccccccccccccccccccccccc

      deallocate(AU,SmU,ResU,dSolU)
      deallocate(IndLigneAU,NLigneAU,NColAU)

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
c
c
c      
      subroutine ResiduContact(
     &     beta,rF, 
     &     NbFaceFrac,NbFace,
     &     NbdofMeca,NbdofMecaContact,MemBlocU,  
     &     NbIncGlobalbyFaceFrac,NumIncGlobalbyFaceFrac,
     &     SautbyFaceFrac,
     &     XFaceCG,NumFaceFracVersFace,
     &     VecNormalbyFace,
     &     VecTanFace1,VecTanFace2,      
     &     IndLigneAU,NLigneAU,NColAU,AU,       
     &     SolU,SolU_nm1,SmU,
     &     IndContact,ResU)      
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
       dimension IndLigneAU(NbdofMecaContact+1)
       dimension NLigneAU(MemBlocU)
       dimension NColAU(MemBlocU)

       dimension AU(MemBlocU,NbDim,NbDim)          
c
c
c     Solution U courante du Newton 
c
       dimension SolU(NbdofMecaContact,NbDim)
c
c
c     Solution U a tnm1 
c
       dimension SolU_nm1(NbdofMecaContact,NbDim)
c
c       
c
c     second membre de l'equation (or equation de "contact") 
c
       dimension SmU(NbdofMecaContact,NbDim)
c
c
c
c     coefficient de frottement par face frac 
c
       dimension rF(NbFaceFrac)
c
c       
ccccccccccccccccccccccccccccccccc
c
c     Sortie 
c
c
c     Residu Residu 
c
       dimension ResU(NbdofMecaContact,NbDim)
c
c
c     Indice de contact par face frac
c
       dimension IndContact(NbFaceFrac)
c
cccccccccccccccccccccccccccccc
c
c     workspaces 
c
       dimension Vec1(NbDim), Vec2(NbDim)

       dimension Saut(NbDim), Sautnt(NbDim)  

       dimension bnfrac(NbDim)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     normal stiffness
c
       unsurkn = f_UnSurkn()
       

       eps = 1.d-10       ! pour éviter la division par 0
c
c
c     rlambdan = qn*bn si contact avec 
c       
c
       qn = 1.d0/( 1.d0 + beta*unsurkn )
c
c     Cohesion 
c
       cohesion = f_cohesion()

       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       
cccccccccccccccccccccccccccccccccccccccc
c     Calcul du residu AU*SolU - SmU hors equations du contact 
cccccccccccccccccccccccccccccccccccccccc
c     
       do inc=1,NbdofMeca
          
          ResU(inc,:) = - SmU(inc,:)
          
          do m=IndLigneAU(inc)+1,IndLigneAU(inc+1)
             
             call ProdMatVec(AU(m,:,:),SolU(NColAU(m),:),Vec1)             
             
             ResU(inc,:) = ResU(inc,:) + Vec1(:)
             
          enddo
          
       enddo
c     
c     
cccccccccccccccccccccccccccccccccccccccccccccc
c     on rajoute les equation du  contact
cccccccccccccccccccccccccccccccccccccccccccccc
c     
       
       do ifrac = 1,NbFaceFrac

          
          Fcoef = rF(ifrac)
          
          
          nf = NumFaceFracVersFace(ifrac)
          
          incfrac = NbdofMeca + ifrac
c     
c     tractions normales et tangentielles 
c     
          rlambdan  = SolU(incfrac,1)
          rlambdat1 = SolU(incfrac,2)
          rlambdat2 = SolU(incfrac,3)
                      
c     
c     
c     calcul du saut normal et tangentiel
c     -> Sautnt(1) composante normale de [U],
c     Sautnt(2,3) composantes tangentielles [U-Unm1]   
c
c     Saut [U] 
c          
          Saut(:) = 0.d0 
          do ik=1,NbIncGlobalbyFaceFrac(ifrac)
             inc = NumIncGlobalbyFaceFrac(ifrac,ik)
             Saut(:) = Saut(:) + SautbyFaceFrac(ifrac,ik)*SolU(inc,:)
          enddo         
          
          Sautnt(1) = prodscal(Saut,VecNormalbyFace(nf,:))  ! [U]n 
c
c     Saut [U-Unm1]
c
          do ik=1,NbIncGlobalbyFaceFrac(ifrac)
             inc = NumIncGlobalbyFaceFrac(ifrac,ik)
             Saut(:) = Saut(:)
     &            - SautbyFaceFrac(ifrac,ik)*SolU_nm1(inc,:)
          enddo
          
          Sautnt(2) = prodscal(Saut,VecTanFace1(nf,:)) ! [U-Unm1]_t1 
          Sautnt(3) = prodscal(Saut,VecTanFace2(nf,:)) ! [U-Unm1]_t2 
          
          
          bnfrac(1) = rlambdan  + beta*Sautnt(1)
          bnfrac(2) = rlambdat1 + beta*Sautnt(2)
          bnfrac(3) = rlambdat2 + beta*Sautnt(3)  
          
          btaunorme = dsqrt( bnfrac(2)**2 + bnfrac(3)**2 )  
          
          
          
          if (bnfrac(1).gt.0.0) then ! on a du contact
             
             
             if (btaunorme .le. Fcoef*bnfrac(1)*qn + cohesion + eps)
     &            then
                
c                write(*,*)ifrac,' contact stick ',bnfrac(:)                
                IndContact(ifrac) = 1 ! contact with stick
                
                
                ResU(incfrac,1) = beta*unsurkn*rlambdan - beta*Sautnt(1)
                
                ResU(incfrac,2) = - beta*Sautnt(2) 
                
                ResU(incfrac,3) = - beta*Sautnt(3) 

                
             else

c                write(*,*)ifrac,' contact slip '
                IndContact(ifrac) = 2 ! contact with slip
                
                ResU(incfrac,1) = beta*unsurkn*rlambdan - beta*Sautnt(1)
                
                ResU(incfrac,2) =  rlambdat1 - bnfrac(2)/
     &               btaunorme*(Fcoef*bnfrac(1)*qn + cohesion) 


                ResU(incfrac,3) =  rlambdat2 - bnfrac(3)/
     &               btaunorme*(Fcoef*bnfrac(1)*qn + cohesion) 
                
   
             endif
             
             
          else

c             write(*,*)ifrac,' no contact '             
             IndContact(ifrac) = 0 ! pas de contact 
             
             ResU(incfrac,1) = rlambdan
             ResU(incfrac,2) = rlambdat1
             ResU(incfrac,3) = rlambdat2    
             
          endif
          
          
       enddo
       
       
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       return
       end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       
c
c
c      
      subroutine JACContact(
     &     beta,rF, 
     &     NbFaceFrac,NbFace,
     &     NbdofMeca,NbdofMecaContact,MemBlocU,  
     &     NbIncGlobalbyFaceFrac,NumIncGlobalbyFaceFrac,
     &     SautbyFaceFrac,
     &     NumFaceFracVersFace,
     &     VecNormalbyFace,
     &     VecTanFace1,VecTanFace2,      
     &     IndLigneAU,NColAU,AU,SolU,SolU_nm1,       
     &     IndContact)      
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
c     numero de la face pour une face frac donnee  
c        
        dimension NumFaceFracVersFace(NbFaceFrac)        
c
c     CSR matrice meca et matrice sans les eqs du mulriplicateur 
c
       dimension IndLigneAU(NbdofMecaContact+1)
       dimension NColAU(MemBlocU)   
c
c     coefficient de frottement par face frac 
c
       dimension rF(NbFaceFrac)
c
c
c
c     Solution U courante du Newton 
c
       dimension SolU(NbdofMecaContact,NbDim)
c
c
c
c    U a tnm1
c
       dimension SolU_nm1(NbdofMecaContact,NbDim)
c
c       
c
c     Indice de contact par face frac
c
       dimension IndContact(NbFaceFrac)
c
c
c       
ccccccccccccccccccccccccccccccccc
c
c     Sortie 
c
c     Jacobienne avec termes lies aux equations du contact 
c
       dimension AU(MemBlocU,NbDim,NbDim)      
c
cccccccccccccccccccccccccccccc
c
c     workspaces 
c
       dimension Saut(NbDim), Sautnt(NbDim)  

       dimension bnfrac(NbDim)
c
c       
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     normal stiffness
c
       unsurkn = f_UnSurkn()
       

       eps = 1.d-10       ! pour éviter la division par 0
c
c
c     rlambdan = qn*bn si contact avec 
c       
c
       qn = 1.d0/( 1.d0 + beta*unsurkn )
c
c     Cohesion 
c
       cohesion = f_cohesion()
c
c        
ccccccccccccccccccccccccccc
c     
c     remise a zero des lignes du contact dans la jacobienne 
c         
         do ifrac=1,NbFaceFrac
            incfrac = NbdofMeca + ifrac
            do m=IndLigneAU(incfrac)+1,IndLigneAU(incfrac+1)
               AU(m,:,:) = 0.d0 
            enddo            
         enddo
c     
c
c     Rajout des equations du contact dans la Jacobienne AU 
c
c
         do ifrac = 1,NbFaceFrac

            
            Fcoef = rF(ifrac) ! coefficient de frottement par face frac 

            
            incfrac = NbdofMeca + ifrac
            nf = NumFaceFracVersFace(ifrac)
                        
            if (IndContact(ifrac).eq.0) then ! on n'a pas de contact
               
               mii = IndLigneAU(incfrac)+1
               
               AU(mii,1,1) = 1.d0 ! rlambdan = 0
               AU(mii,2,2) = 1.d0 ! rlambdantau1 = 0
               AU(mii,3,3) = 1.d0 ! rlambdantau2 = 0
               
               
            else if (IndContact(ifrac).eq.1) then ! on a du sticky contact


               mii = IndLigneAU(incfrac)+1
               
               AU(mii,1,1) = beta*unsurkn ! unsurkn*beta*rlambdan  - beta*sautn(u) = 0  der / rlambdan               
                  
                  do ik=1,NbIncGlobalbyFaceFrac(ifrac)
                     inc = NumIncGlobalbyFaceFrac(ifrac,ik)
                     
                     mij = 0    ! recherche elt (incfrac,inc) dans la ligne inc du CSR 
                     do m=IndLigneAU(incfrac)+1,IndLigneAU(incfrac+1)
                        n = NColAU(m)
                        if (n.eq.inc) then
                           mij = m 
                        endif
                     enddo
                     
                     if (mij.eq.0) then
                        write(*,*)'elt incincfrac non trouve dans CSR',
     &                    ifrac,inc
                        stop
                     endif  

ccccccccccccccc pour l'équation unsurkn*beta*rlambdan - beta*sautn(u) = 0  der / u             
                        
                        AU(mij,1,:) = AU(mij,1,:) 
     &                       - beta*SautbyFaceFrac(ifrac,ik)
     &                       *VecNormalbyFace(nf,:)
                        
                        

cccccccccccccccpour l'équation - beta*sauttau1(u) = 0  
                        
                        AU(mij,2,:) = AU(mij,2,:) 
     &                       - beta*SautbyFaceFrac(ifrac,ik)
     &                       *VecTanFace1(nf,:)
                        

                        
ccccccccccccccc pour l'équation - beta*sauttau2(u) = 0                
c
                        AU(mij,3,:) = AU(mij,3,:) 
     &                       - beta*SautbyFaceFrac(ifrac,ik)
     &                       *VecTanFace2(nf,:)
                                                
                     
                  enddo  ! end boucle dof
                  
                  
                  
               else             ! on a du slippy contact


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccc   Calcul de bn, btau1, btau2, lamdan, lambdatau1, lambdatau2, sautn(u), sauttau1(u) et sauttau2(u) etc....  cccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


                  rlambdan  = SolU(incfrac,1)
                  rlambdat1 = SolU(incfrac,2)
                  rlambdat2 = SolU(incfrac,3)
                  
                     
c
c     Saut [U] 
c                       
                  Saut(:) = 0.d0 
                  do ikk=1,NbIncGlobalbyFaceFrac(ifrac)
                     incc = NumIncGlobalbyFaceFrac(ifrac,ikk)
                     Saut(:) = Saut(:) +
     &                    SautbyFaceFrac(ifrac,ikk)*SolU(incc,:)
                  enddo         
                     
                  Sautnt(1) = prodscal(Saut,VecNormalbyFace(nf,:)) ! sautn(u)


c
c     Saut [U-Unm1]
c
                  do ikk=1,NbIncGlobalbyFaceFrac(ifrac)
                     incc = NumIncGlobalbyFaceFrac(ifrac,ikk)
                     Saut(:) = Saut(:)
     &                    - SautbyFaceFrac(ifrac,ikk)*SolU_nm1(incc,:)
                  enddo

                  
                  Sautnt(2) = prodscal(Saut,VecTanFace1(nf,:)) ! sauttau1(u-unm1)
                  Sautnt(3) = prodscal(Saut,VecTanFace2(nf,:)) ! sauttau2(u-unm1)
                  
                     
                  bnfrac(1) = rlambdan  + beta*Sautnt(1) ! bn
                  bnfrac(2) = rlambdat1 + beta*Sautnt(2) ! btau1
                  bnfrac(3) = rlambdat2 + beta*Sautnt(3) ! btau2
                  
                  btaunormecarre = bnfrac(2)**2 + bnfrac(3)**2 ! || btau ||^2 
                  



                  mii = IndLigneAU(incfrac)+1

               
                  AU(mii,1,1) = beta*unsurkn ! eq beta*unsurkn*rlambdan  - beta*sautn(u) = 0  der / rlambdan    
                  
                        
ccccccccccccccc pour l'équation lambdatau1 - (btau1/normebtau)*(F*bn*qn + cohesion) = 0  der / rlambda 




                  AU(mii,2,1) = - qn*Fcoef*bnfrac(2)
     &                 /(btaunormecarre**(0.5d0))   
                     


                  AU(mii,2,2) = 1.d0 - (qn*Fcoef*bnfrac(1)+cohesion)
     &                 * (bnfrac(3)**2)
     &                 / (btaunormecarre**(1.5d0))   


                  AU(mii,2,3) =
     &                 (qn*Fcoef*bnfrac(1)+cohesion)*bnfrac(2)*bnfrac(3)
     &                 / (btaunormecarre**(1.5d0)) 





ccccccccccccccc pour l'équation lambdatau2 - (btau2/normebteau)(F*bn*qn+coh)  = 0   der / rlambda 


                  AU(mii,3,1) = - qn*Fcoef*bnfrac(3)
     &                 /(btaunormecarre**(0.5d0))   


                  AU(mii,3,2) =
     &                 (qn*Fcoef*bnfrac(1)+cohesion)*bnfrac(2)*bnfrac(3)
     &                 / (btaunormecarre**(1.5d0))   
                  
                     

                  AU(mii,3,3) = 1.d0 - (qn*Fcoef*bnfrac(1)+cohesion)
     &                 *(bnfrac(2)**2)
     &                 /(btaunormecarre**(1.5d0))   
                  

                 

                  

   
                  do ik=1,NbIncGlobalbyFaceFrac(ifrac)
                     inc = NumIncGlobalbyFaceFrac(ifrac,ik)
                     
                     mij = 0    ! recherche elt (incfrac,inc) dans la ligne inc du CSR 
                     do m=IndLigneAU(incfrac)+1,IndLigneAU(incfrac+1)
                        n = NColAU(m)
                        if (n.eq.inc) then
                           mij = m 
                        endif
                     enddo
                     
                     if (mij.eq.0) then
                        write(*,*)'elt inc,incfrac non trouve dans CSR',
     &                     ifrac,inc
                        stop
                     endif  
                     
                           
                           

ccccccccccccccc pour l'équation beta*unsurkn*rlambdan - beta*sautn(u) = 0  terme de sautn             
                      



                
                        AU(mij,1,:) = AU(mij,1,:) !  beta*unsurkn*rlambdan -beta*Sautn = 0  der / u 
     &                       - beta*SautbyFaceFrac(ifrac,ik)
     &                       *VecNormalbyFace(nf,:)
                     



ccccccccccccccc pour l'équation lambdatau1 - (btau1/normebtau)(F*bn*qn + coh) = 0

                        AU(mij,2,:) = AU(mij,2,:) 
     &                       - beta*qn*Fcoef*SautbyFaceFrac(ifrac,ik)                                        
     &                         *bnfrac(2)/(btaunormecarre**(0.5d0))
     &                         *VecNormalbyFace(nf,:)                        
                        
                         AU(mij,2,:) = AU(mij,2,:) 
     &                       + beta*SautbyFaceFrac(ifrac,ik)
     &                         *(qn*Fcoef*bnfrac(1) + cohesion)*(                  
     &                       - bnfrac(3)**2
     &                       /(btaunormecarre**(1.5d0))
     &                       * VecTanFace1(nf,:)                                
     &                       +  bnfrac(2)*bnfrac(3)
     &                       /(btaunormecarre**(1.5d0))
     &                       *VecTanFace2(nf,:)
     &                                                )
                        
                        



ccccccccccccccc pour l'équation lambdatau2 - (btau2/normebtau)F*bn*qn  = 0



                       AU(mij,3,:) = AU(mij,3,:) 
     &                       - beta*qn*Fcoef*SautbyFaceFrac(ifrac,ik)                                     
     &                        *bnfrac(3)/(btaunormecarre**(0.5d0)) 
     &                        *VecNormalbyFace(nf,:)                       

                        
                        


                       AU(mij,3,:) = AU(mij,3,:) 
     &                      + beta*SautbyFaceFrac(ifrac,ik)
     &                         *(qn*Fcoef*bnfrac(1) + cohesion)*(                          
     &                       - bnfrac(2)**2
     &                       /(btaunormecarre**(1.5d0))
     &                       *VecTanFace2(nf,:)                                       
     &                       + bnfrac(2)*bnfrac(3)
     &                       /(btaunormecarre**(1.5d0))
     &                       * VecTanFace1(nf,:)
     &                                                )
                        
                        




                     
                  enddo  ! end boucle dof
                                 
            endif   ! en if indcontact 
            
               
         enddo ! end boucle sur les fractures
            



       
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       return
       end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       
c
