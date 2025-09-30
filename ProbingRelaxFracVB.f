
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Calcul du coefficient de relaxation diagonal sur les fractures par PROBING avec dpf constant 
c     
c
c     En entree SolU, SolU_nm1, SolP, IndContact t.q. F(U,Lambda,SolP, SolT) = 0 
c
c     (SolT n'intervient pas ici car T est fixee et sa dependance est lineaire) 
c      
c     En sorties relaxff par probing
c
c     relaxmm non calcule ici (on prend juste le Crm du fixed stress) 
c
c
c     calcul de dSolU tq  dF/dX dX + dF/dP dP + dF/dT dT = 0 avec X=(U,Lambda) 
c
c     on deduit ensuite la variation linearisee ddf de df -> relaxff = ddf/dpf 
c      
c
c     RQ: la dependance en P,T de F(X,P,T) est lineaire en formulation mixte 
c
c     RQ: T etant fixee et n'intervenant qu'au Sm, elle n'intervient pas dans le calcul 
c
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c      
      subroutine ProbingRelaxFracVB( 
     &     NbNode,NbCell,NbFace,NbFaceFrac,
     &     NbdofMeca,NbdofMecaContact,
     &     rF, 
     &     IndDirNodeMeca,IndFace,          
     &     NbNodebyCell,NumNodebyCell,
     &     NbNodebyFace,NumNodebyFace,
     &     NbFacebyCell,NumFacebyCell,
     &     NumCellbyFace,NumFaceFracVersFace,      
     &     VecNormalbyFace,PoidsXCellCG,PoidsXFaceCG,      
     &     XNode,XCellCG,XFaceCG,SurfaceFace,VolCell,      
     &     NbIncGlobalbyCell,NumIncGlobalbyCell,
     &     NumNodeFaceLocalbyCell,
     &     NbIncGlobalbyFaceFrac,NumIncGlobalbyFaceFrac,
     &     ACell,SautbyFaceFrac,GradCellVEM,
     &     VecTanFace1,VecTanFace2, 
     &     SolU,SolU_nm1,SolU0,SolP,IndContact,
     &     NbCV,NumIncFaceFrac,NumIncCell,
     &     relaxff,relaxmm)      
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
c     Solution U 
c
       dimension SolU(NbdofMecaContact,NbDim)
c
c     Solution U a tnm1 
c
       dimension SolU_nm1(NbdofMecaContact,NbDim)

c
c     Solution U init 
c
       dimension SolU0(NbdofMecaContact,NbDim)       
c
c
       dimension IndContact(NbFaceFrac)
c
c       
cccccccccccccccccccccccccccccc
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
cc   vecteur gradiant pour chaque fonction de bases scalaire v^dof (dof=1,...,NbdofbyCell(K)) pour chaque cellule K

      dimension GradCellVEM(NbCell,NbdofbyCellMax,NbDim) 
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
       dimension relaxff(NbFaceFrac)
       dimension relaxmm(NbCell)       
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
       beta = 1.d+6 ! parameter du semi-smooth Newton pour le contact 
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
ccccccccccccccccccccccccccccccccc
c      
c     dpf = 1.d+5 ici 
c
cc  Calcul du terme (dans le second membre) \sum_{sigma} \int_{sigma} pf(sigma)*saut(v)_n

      dpf = 1.d+5 ! constant le long des frac 
      
      do ifrac = 1,NbFaceFrac
         nf = NumFaceFracVersFace(ifrac)
         inck = NumIncFaceFrac(ifrac)
     
         do ik=1,NbIncGlobalbyFaceFrac(ifrac)
            inc = NumIncGlobalbyFaceFrac(ifrac,ik)


              SmU(inc,:) = SmU(inc,:)
     &              - dpf*SurfaceFace(nf)
     &              * SautbyFaceFrac(ifrac,ik)*VecNormalbyFace(nf,:)       
          
               
         
         enddo         
         
      enddo
c
c
c     
c
cc  calcul du terme (dans le second membre) \sum_k int_{K} b* dpf * div(v) 
c
c
c$$$
c$$$      b_coef = f_biot()
c$$$
c$$$
c$$$      do k=1,NbCell
c$$$
c$$$         do ik=1,NbIncGlobalbyCell(k)
c$$$
c$$$            inc = NumIncGlobalbyCell(k,ik)
c$$$
c$$$
c$$$            SmU(inc,:) = SmU(inc,:)
c$$$     &           + VolCell(k)*b_coef*dpf 
c$$$     &           *GradCellVEM(k,ik,:)
c$$$
c$$$         enddo
c$$$         
c$$$      enddo
c$$$
c$$$            
c
c      
ccccccccccccccccccccccccccccccccccc      
c
c     Correction du SmU par les CL de Dirichlet au bord 
c
cccccccccccccccccccccccccccccccccc
c
      
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

c
c     calcul du residu ResU et de IndContact
c      
c      call ResiduContact(
c     &     beta,rF, 
c     &     NbFaceFrac,NbFace,
c     &     NbdofMeca,NbdofMecaContact,MemBlocU,  
c     &     NbIncGlobalbyFaceFrac,NumIncGlobalbyFaceFrac,
c     &     SautbyFaceFrac,    
c     &     XFaceCG,NumFaceFracVersFace,
c     &     VecNormalbyFace,
c     &     VecTanFace1,VecTanFace2,          
c     &     IndLigneAU,NLigneAU,NColAU,AU,       
c     &     SolU,SolU_nm1,SmU,
c     &     IndContact,ResU)     

      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     calcul de dSolU lie au dpf 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      
      allocate(dSolU(NbdofMecaContact,NbDim))

      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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

      
      ResU = SmU 

      
c      NbIncMeca = NbDim
c      IPREC = 0  ! ILU0 
c      MAXL = 1000
c      call solveur_iteratif(
c     &       NbdofMecaContact,MemBlocU,NbIncMeca,
c     &       AU,NLigneAU,NColAU,IndLigneAU,
c     &       ResU,dSolU,critere_arret_gmres,MAXL,IPREC,
c     &       ITER,IERR,ERR)      

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
      

c      write(*,*)' debut solveur SuperLU '
c$$$      ITER = 0
c$$$      NbIncMeca = NbDim 
c$$$      call solveurSuperLU(NbdofMecaContact,MemBlocU,NbIncMeca,
c$$$     &     AU,NLigneAU,NColAU,IndLigneAU,
c$$$     &     ResU,dSolU)

      
cccccccccccccccccccccccccccccccc
c
c     calcul des sauts moyens SautU (vec1) et dSautU (vec2)par face frac 
c
c     puis du coeff de relaxation fracture par face frac 
c


      SDilation = f_sheardilation()

      
      do ifrac = 1,NbFaceFrac
         nf = NumFaceFracVersFace(ifrac)
         
         Vec1(:) = 0.d0
         Vec2(:) = 0.d0
         Vec3(:) = 0.d0
         
         do ik=1,NbIncGlobalbyFaceFrac(ifrac)
            inc = NumIncGlobalbyFaceFrac(ifrac,ik)
            
            Vec1(:) = Vec1(:) + SautbyFaceFrac(ifrac,ik)*SolU(inc,:)
            
            Vec3(:) = Vec3(:) + SautbyFaceFrac(ifrac,ik)
     &           *(SolU(inc,:)-SolU0(inc,:))
            
            Vec2(:) = Vec2(:) + SautbyFaceFrac(ifrac,ik)*dSolU(inc,:)
            
         enddo        

         dsn = prodscal(Vec2(:),VecNormalbyFace(nf,:)) 


         ddf = 0.d0 

c         if (IndContact(ifrac).eq.0) then 
            ddf = - dsn 
c         endif

         if (IndContact(ifrac).ne.1) then 
c         
c     saut et dsaut tangentiels 
c
            st1 = prodscal(Vec3(:),VecTanFace1(nf,:))
            st2 = prodscal(Vec3(:),VecTanFace2(nf,:))         

            dst1 = prodscal(Vec2(:),VecTanFace1(nf,:))
            dst2 = prodscal(Vec2(:),VecTanFace2(nf,:))         

            dnormest = ( st1*dst1 + st2*dst2 )/dsqrt( st1**2 + st2**2)

            ddf = ddf + f_sheardilation()*dnormest
         
         
c            write(*,*)' relaxff ',ifrac,relaxff(ifrac),dsn,st2,dst2 
            
            
         endif

         relaxff(ifrac) = dabs(ddf)/dpf   


         
      enddo
      
ccccccccccccccccccc
c
c     calcul de div U et relaxmm 
c
c$$$      do k=1,NbCell
c$$$ 
c$$$         divk = 0.d0 
c$$$         do ik=1,NbIncGlobalbyCell(k)
c$$$            inc = NumIncGlobalbyCell(k,ik)
c$$$
c$$$            do i = 1,NbDim
c$$$               divk = divk +  SolU(inc,i)
c$$$     &                       *GradCellVEM(k,ik,i) ! divergence de U cell k 
c$$$            enddo       
c$$$         enddo
c$$$
c$$$
c$$$
c$$$         relaxmm(k) = - divk*b_coef/dpf 
c$$$
c$$$         
c$$$
c$$$      enddo
      
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
