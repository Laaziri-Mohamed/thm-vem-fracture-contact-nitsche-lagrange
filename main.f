ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     Modele: Elasticite lineaire avec fracture + contact frictionnel de type Coulomb 
c
c
c     Schema: VEM d'ordre 1 avec bulle(s) aux faces frac 
c
c     Orientation de la normale locale a la face (sortant cote maille 1): VecNormalbyFace
c      
c     En revanche les deux vecteurs tangents ont une orientation intrinseque VecTanFace1 et VecTanFace2  
c
c     -> pour la visu il faut corriger les deux composantes tangentielles du saut u^+-u^- et de sa variable duale Lambda (multiplicateur)
c        pour avoir une orientation intrinseque       : correction par multiplication avec OrientationbyFaceFrac 
c
c     -> mais par definition la composante normale du saut sautn = (u^+ - u^-).n^+
c                           et celle lambdan de la variable duale lambda du saut sont intrinseques 
c      
cccccccccccccccccccccccccccccccccc
c      
c     fluide incompressible discretise avec le schema HFV, modele a pression discontinue 
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     
c     
      program main
c     
c     
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c     
c     
c     
      include 'include_parameter'
c     
c     affichage paraview 
c     
      parameter ( iaffich = 1)


c         parameter ( iMecaVBulle = 1 ) ! VEM Bulle ->  on met au moins une bulle sur les faces fra
         parameter ( iMecaVBulle = 0 ) ! Nitsche   ->  pas de bulle 
      
c
c     Permeabilite matrice rock type 1,2,3 
c
         parameter ( Permeabilitem1 = 1.d-15 )
         parameter ( Permeabilitem2 = 1.d-19 )
c         parameter ( Permeabilitem2 = 1.d-18 )                  
c          parameter ( Permeabilitem2 = 1.d-17 )         
c         parameter ( Permeabilitem3 = 1.d-13 )
         parameter ( Permeabilitem3 = 6.d-13 )         
c
c     Porosite init 
c         
         parameter ( Phiminit1 = 0.1d0 )
         parameter ( Phiminit2 = 0.01d0 )         
         parameter ( Phiminit3 = 0.1d0 )         

c
c     Viscosite 
c
         parameter ( Viscof = 1.d-3 )  
c
c     Conductivite thermique 
c
c         parameter ( CondThermique = 2.d-3 ) 
         parameter ( CondThermique = 2.d0 )
         
c
c     Epaisseur fracture au contact 
c
c          parameter ( dfmin = 1.d-4 )

          parameter ( dfmin = 1.d-5 )

c         parameter ( dfmin = 5.d-6 )                
c
c      
ccccccccccc
c
c     Temps de simulation
c
         parameter ( UneHeure = 3600.d0 )
         parameter ( UnJour = 24.d0*UneHeure )
         parameter ( UnAn = UnJour*365.d0 )

cccccccccccccccccccc
c      
c        
c     Temps final et marche en temps 
c



         parameter ( TempsFinal = UnJour*10000.d0 )

         
         parameter ( Deltat1  = UnJour*0.25d0  )
         parameter ( Deltat2  = UnJour*1.d0  )         
         parameter ( Deltat3  = UnJour*10.d0  )
         parameter ( Deltat4  = UnJour*50.d0  )            
         parameter ( Deltatinit = Deltat1  ) 
         parameter ( DeltatMin  = Deltat2  )


c
ccccccccccccccccccc        
c
c     affichage ts les naffichpasdetemps 
c
      parameter ( naffichpasdetemps = 1 )
c      parameter ( naffichpasdetemps = 10 )
         
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     DECLARATIONS
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     LECTURE DU FORMAT BENCHMARK FVCA6 3D 
c     
ccccccccccccccccccccccccccccccccccccccc
c     
c     coordonnees des noeuds XS 
c     
      double precision, dimension(:,:), allocatable :: XNode
c     
c     nb de noeuds S par maille K 
c     
      integer, dimension(:), allocatable :: NbNodebyCell
c     
c     num des noeuds S de chaque maille K dans l'ordre cyclique 
c     
      integer, dimension(:,:), allocatable :: NumNodebyCell
c     
c     Maille par les faces 
c     
      integer, dimension(:), allocatable :: NbFacebyCell
      integer, dimension(:,:), allocatable :: NumFacebyCell
c     
c     Faces par les aretes 
c     
      integer, dimension(:), allocatable :: NbAretebyFace
      integer, dimension(:,:), allocatable :: NumAretebyFace
c     
c     Faces par les noeuds 
c     
      integer, dimension(:), allocatable :: NbNodebyFace
      integer, dimension(:,:), allocatable :: NumNodebyFace
c     
c     Mailles voisines des faces 
c     
      integer, dimension(:,:), allocatable :: NumCellbyFace
c     
c     2 noeuds de chaque arete 
c     
      integer, dimension(:,:), allocatable :: NumNodebyArete   
c     
c     Label = 1 si face frac 
c     
      integer, dimension(:), allocatable :: LabelbyFace

c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Ensemble des mailles voisines d'un noeud 
c
c        dimension NbCellbyNode(NbNode)
c        dimension NumCellbyNode(NbNode,NbCellbyNodeMax)
c
c     
      integer, dimension(:), allocatable :: NbCellbyNode    
      integer, dimension(:,:), allocatable :: NumCellbyNode


c      
c     classe d'equivalence Ks autour du noeud -> un label par classe 
c     on stocke le numero du label pour chaque maille autour du noeud : LabelCellbyNode     
c      
c        dimension LabelCellbyNode(NbNode,NbCellbyNodeMax)
c        dimension NbLabelbyNode(NbNode)

      integer, dimension(:), allocatable :: NbLabelbyNode    
      integer, dimension(:,:), allocatable :: LabelCellbyNode        
c
c
c
c     Indice 1 si Node frac, 0 sinon 
c        
c     dimension IndNodeFrac(NbNode)
      integer, dimension(:), allocatable :: IndNodeFrac          
c
c
cccccccccc
c
c
c     Numerotation locale par cell ds l'ordre NodebyCell inc Ks
c               puis les inc face frac Ksigma si  k=NumCellbyFace(n,1) (cote maille 1 de la face) 
c     Local -> global 
c        
c        dimension NbIncGlobalbyCell(NbCell)
c        dimension NumIncGlobalbyCell(NbCell,NbIncbyCellMax)        
c
      integer, dimension(:), allocatable :: NbIncGlobalbyCell
      integer, dimension(:,:), allocatable :: NumIncGlobalbyCell          
c
c
c     Numero local du noeud ou de la face par inc cell 
c
c        dimension NumNodeFaceLocalbyCell(NbCell,NbIncbyCellMax)                
c
      integer, dimension(:,:), allocatable :: NumNodeFaceLocalbyCell           
c
cccccccccc
c
c     Numerotation globale des inc de l'operateur saut pour chaque face frac pour l'ordre local suivant:
c        - nodes de la face cote maille 1 
c        - nodes de la face cote maille 2 
c        - face cote maille 1 si bulle ( cf si IndBulleNumCellbyface(n,1) = 1 )
c        - face cote maille 2 si bulle ( cf si IndBulleNumCellbyface(n,2) = 1 )
c
c        dimension NbIncGlobalbyFaceFrac(NbFaceFrac)
c        dimension NumIncGlobalbyFaceFrac(NbFaceFrac,2*NbNodebyFaceMax+2)              
c
      integer, dimension(:), allocatable :: NbIncGlobalbyFaceFrac
      integer, dimension(:,:), allocatable :: NumIncGlobalbyFaceFrac          
c
c      
ccccccccccccccccccccccccccccccccccccc
c     
c     Structure Mailles K, Faces fractures sig, IK, Isig
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Identification des faces fractures parmi l'ensemble des faces 
c     
c     num face frac > num face  
      integer, dimension(:), allocatable :: NumFaceFracVersFace     
c     
c     num face  > num face frac
      integer, dimension(:), allocatable :: NumFaceVersFaceFrac    
c     
c     num face > ind 1 si face frac 0 sinon
      integer, dimension(:), allocatable :: IndFaceFrac        
c     
c     Choix de la CL par face pour P,T (flow)
c
c     Pour la pression     
c     
      integer, dimension(:), allocatable :: IndDirFace ! Dir 1, 0 sinon (neumann ou int)
      integer, dimension(:), allocatable :: IndNeuFace ! Neumann homogene 1, Neumann non homogene 2
c
c     Pour la Temperature 
c
      integer, dimension(:), allocatable :: IndDirFaceT ! Dir 1, 0 sinon 
      integer, dimension(:), allocatable :: IndNeuFaceT ! Neumann homogene 1, 0 sinon       
c
c
c      
c     Node Dirichlet pour la Meca 
c      
      integer, dimension(:), allocatable :: IndDirNodeMeca       
c      
c     
c     Numerotation des inconnues 
c     
c     Inc Interfaces Dir 
c     Inc Mailles 
c     Inc Face Frac (on suppose qu'elles ne sont pas sur le bord)
c     Autres Inc Interfaces (faces ou noeuds ou aretes) 
c     
c     
c     
c     
c     Numero des inconnues mailles et faces fractures 
c     
c     num cell       >  num inc cell
c     num face frac  >  num inc face frac  
c     
      integer, dimension(:), allocatable :: NumIncCell      
      integer, dimension(:), allocatable :: NumIncFaceFrac       
c     
c     Numero des inconnues interfaces 
c     
c     pour VAG (tous les noeuds) 
c     num node  > num inc  
c     
      integer, dimension(:), allocatable :: NumIncNode       
c     
c     pour HFV (toutes les faces y compris les faces frac)
c     num face  > num inc 
c     
      integer, dimension(:), allocatable :: NumIncFace 
c     
c     pour HFV (toutes les aretes fracture)
c     
c     Ind 11 arete frac DIR
c     Ind 10 arete frac non DIR 
c     Ind 0 non frac non DIR
c     Ind 1 non frac DIR 
c     
      integer, dimension(:), allocatable :: IndIncArete ! pression
      integer, dimension(:), allocatable :: IndIncAreteT ! temperature        
c     
c     num arete  > num inc 
      integer, dimension(:), allocatable :: NumIncArete
c     
c   
c     
c     pour HFV discontinu 
c     num face frac  > num inc Ksigma et Lsigma 
c     
      integer, dimension(:,:), allocatable :: NumIncCellbyFaceFrac       
c
c
c     Ind 1 si Inc frac (arete ou face frac), 0 sinon
c      
      integer, dimension(:), allocatable :: IndIncFaceFrac
c     
c
c      
ccccccccccccc
c     
c     Numero des inconnues d'interfaces IK par maille et Isig par face fracture
c     
c     num cell x  num local interface cell             > num inc interface 
c     num face frac x  num local interface face frac   > num inc interface 
c     
c     
      integer, dimension(:), allocatable :: NbInterfacebyCell
      integer, dimension(:,:), allocatable :: NumInterfacebyCell

      integer, dimension(:), allocatable :: NbInterfacebyFaceFrac
      integer, dimension(:,:), allocatable :: NumInterfacebyFaceFrac

c     
c     Transmissivites TK et Tsig 
c     num cell x num local interface cell x num local interface cell                 > TKij 
c     num face frac x num local interface face frac x num local interface face frac  > Tsigij 
c     
      double precision, dimension(:,:,:), allocatable :: TransCell

      double precision, dimension(:,:,:), allocatable :: 
     &     TransFaceFrac


      double precision, dimension(:,:,:), allocatable ::
     & TransCellFourier

      double precision, dimension(:,:,:), allocatable :: 
     &     TransFaceFracFourier
c
c
c     Pour HFV discontinu: TransCellbyFaceFrac : TKsigma et TLsigma 
c
      double precision, dimension(:,:), allocatable :: 
     &     TransCellbyFaceFrac

      double precision, dimension(:,:), allocatable :: 
     &     TransCellbyFaceFracFourier    
      
C     
C     
C     Gradient par maille : g = \sum_i\in IK (uk-ui)*GCell(k,i) 
c     
      double precision, dimension(:,:,:), allocatable :: GCell
c     dimension GCell(NbCell,NbInterfacebyCellMax,NbDim)
c     
c     Gradient par face frac : g = \sum_i\in IS (uk-ui)*GFaceFrac(k,i) 
c     
      double precision, dimension(:,:,:), allocatable :: GFaceFrac
c     dimension GFaceFrac(NbFaceFrac,
c     &         NbInterfacebyFaceFracMax,NbDim)
c     
c     
c     Indice inc de bord Dirichlet, Neumann pour la pression et la temperature 
c     
      integer, dimension(:), allocatable :: IndIncDir ! pression 1 Dir, 0 sinon 
      integer, dimension(:), allocatable :: IndIncNeu ! pression 1 neumann homgene, 2 non homogene, 0 sinon 
      integer, dimension(:), allocatable :: IndIncDirT ! temperature 1 Dir, 0 sinon
      integer, dimension(:), allocatable :: IndIncNeuT ! temperature 1 Neumann homogene fourier, 0 sinon 

      
c
c     rock tyoe par maille 
c
      integer, dimension(:), allocatable :: IndRockTypebyCell
c
c
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Isobarycentre des nodes de la maille 
c      
      double precision, dimension(:,:), allocatable :: XCell

c
c     centre de gravite de la maille 
c      
      double precision, dimension(:,:), allocatable :: XCellCG

c
c     Poids barycentriques du centre de gravite de la maille / aux nodes de la maille 
c      
      double precision, dimension(:,:), allocatable :: PoidsXCellCG           
c     
c     coordonnee du centre de Face Frac
c     num face frac  >  XFaceFrac 
c     
      double precision, dimension(:,:), allocatable :: XFaceFrac      
c     
c     coordonnee du centre de Face     
c     
      double precision, dimension(:,:), allocatable :: XFace
c     
c     
c     coordonnee des inconnues (XK, Xsig, Xs, Xa, ...) 
c     
      double precision, dimension(:,:), allocatable :: XInc 
c     
c     Volume des mailles 
c     num cell  >  VolCell 
c     
      double precision, dimension(:), allocatable :: VolCell
c     
c     Volume des faces frac (surface face * epaisseur de la face fracture)
c     num face frac  >  VolFaceFrac 
c     
      double precision, dimension(:), allocatable :: VolFaceFrac
c
c
c
c     coordonnee du centre de gravite de Face     
c
c        dimension XFaceCG(NbFace,NbDim)
c
      double precision, dimension(:,:), allocatable :: XFaceCG
      
c     surface de la face 
c
c        dimension SurfaceFace(NbFace)
c
      double precision, dimension(:), allocatable :: SurfaceFace
      
c     vecteur normal VecKsig oriente sortant de K
c
c        dimension VecNormalKSigma(NbCell,NbFacebyCellMax)
c
      double precision, dimension(:,:,:), allocatable :: VecNormalKSigma

c     vecteur normal par face oriente sortant de la maille 1 de NumCellbyFace 
c
      double precision, dimension(:,:), allocatable :: VecNormalbyFace
c
c     orientation by face frac pour signorini (et pour visu meca)    
c      
      double precision, dimension(:), allocatable ::
     &    OrientationbyFaceFrac       
c      
c     Longueur arete 
c
c        dimension SizeArete(NbArete)
c
      double precision, dimension(:), allocatable :: SizeArete


c     la base orthonormale pour chaque face, une base par face

      double precision, dimension(:,:), allocatable :: VecTanFace1
      double precision, dimension(:,:), allocatable :: VecTanFace2

      double precision, dimension(:), allocatable :: diamK


      
c     Vecteur normal VecSigmaArete oriente sortant de sigma 
c
c       dimension VecNormalSigmaArete(NbFace,NbAretebyFaceMax)
c
      double precision, dimension(:,:,:), allocatable ::
     &                                     VecNormalSigmaArete

c
c
c     Poids barycentriques du CG de la face ds l'ordre des nodes by face 
c
c       dimension PoidsXFaceCG(NbFace,NbNodebyFaceMax)

      double precision, dimension(:,:), allocatable ::
     &                                    PoidsXFaceCG
       
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c      
c     Permeabilite matrice par maille 
c     num cell  >  perm
c
c      
      double precision, dimension(:,:,:), allocatable :: perm        
c     
c     Permeabilite longitudinale fracture par face fracture (permf * epaisseur df )
c     supposee isotrope
c     num face frac > permfdf 
c     
      double precision, dimension(:), allocatable :: permfdf
c
c     Pour HFV discontinu 
c      
c     Permeabilite normale des fracture par face fracture 
c     
      double precision, dimension(:), allocatable :: permfn    
c      
c     
c     
c     Identifiant Face: bords 1,2,3,4,5,6 (boite) ou 0 si face int
c     Centre face 
c     
c     num face  >  indface , xface 
c     
      integer, dimension(:), allocatable :: IndFace
c     
c 
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     solution pression et Temperature
c     
c     
      double precision, dimension(:), allocatable :: SolP
      double precision, dimension(:), allocatable :: SolP_nm1
      double precision, dimension(:), allocatable :: SolP_nm2      
      double precision, dimension(:), allocatable :: SolP_prev


      double precision, dimension(:), allocatable :: SolT
      double precision, dimension(:), allocatable :: SolT_nm1
      double precision, dimension(:), allocatable :: SolT_nm2      
      double precision, dimension(:), allocatable :: SolT_prev
c
c
c     Porosite by cell et der / a Pk Tk 
c
      double precision, dimension(:), allocatable :: PorobyCell
      double precision, dimension(:), allocatable :: PorobyCell_nm1 
      double precision, dimension(:,:), allocatable :: DerPorobyCell
c      
c     Ouverture hydraulique 
c
      double precision, dimension(:), allocatable :: dfbyFaceFrac
      double precision, dimension(:), allocatable :: dfbyFaceFrac_nm1        
c
c     Epaisseur frac au contact by face frac 
c
      double precision, dimension(:), allocatable :: dfcontact         
c
c     Accumulation Darcy et fourier (NbCV) 
c
      double precision, dimension(:), allocatable :: AccP       
      double precision, dimension(:), allocatable :: AccP_nm1

      double precision, dimension(:), allocatable :: AccT       
      double precision, dimension(:), allocatable :: AccT_nm1

      


      
c     
c     
c     CL Dirichlet aux inconnues interfaces de bord DIR 
c     
      double precision, dimension(:), allocatable :: SolPDir
      double precision, dimension(:), allocatable :: SolTDir
c
c     pression et temperature init 
c
      double precision, dimension(:), allocatable :: SolP0

      double precision, dimension(:), allocatable :: SolT0       
c
c     pour visu valeurs de P et T by cell et face frac 
c      
      double precision, dimension(:), allocatable :: PCell
      double precision, dimension(:), allocatable :: PFaceFrac


      double precision, dimension(:), allocatable :: TCell
      double precision, dimension(:), allocatable :: TFaceFrac
c
c
c     pour visu valeurs de traceur by cell et face frac 
c      
      double precision, dimension(:), allocatable :: CCell
      double precision, dimension(:), allocatable :: CFaceFrac      
c
c
c     test relaxation fracture 
c        
      double precision, dimension(:), allocatable :: relaxff

c
c
c     test relaxation matrice
c        
      double precision, dimension(:), allocatable :: relaxmm

      
cccccccccccccccccc
c
c     solution traceur 
c
      double precision, dimension(:), allocatable :: SolC
      double precision, dimension(:), allocatable :: SolC_nm1
c
c
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     VEM 
c
cccccccccccccccccccc
c
c     Indice bulle stocke dans la structure NumCellbyFace : 1 si bulle sur la face, 0 sinon 
c       ->   dans le cas d'une face frac on peut mettre ou non une bulle selon le cote de la frac
c       ->   dans le cas d'une face int non frac, l'indice doit etre le meme des deux cotes 
c      
      integer, dimension(:,:), allocatable :: IndBulleNumCellbyFace      
      
c     Operateur gradient local par maille VEM (avec tous les nodes et toutes les faces) 
c
c     Grad v = sum_s\in VK  vKs*GradCellVEM(k,is,:) + sum_sig\in FK  vKsig*GradCellVEM(k,nnodecell+isig,:) 
c      
c
      double precision, dimension(:,:,:), allocatable :: GradCellAll
c
c     sous operateur avec les seules bulles selectionnes dans IndBulleNumCellbyFace 
c      
      double precision, dimension(:,:,:), allocatable :: GradCellVEM      
c
c     matrice carre de taille nbnodebycell(k) x (nbnodebycell(k) + nbfacebycell(k))
c     telle que dof_nodes (PiK U) = DofPiCell(k,:,:)*U 
c     stockage sur toutes les mailles k 
c
c     All: avec ttes les bulles
c     VEM: sous operateur avec les seules bulles selectionnes dans IndBulleNumCellbyFace 
c      
      double precision, dimension(:,:,:), allocatable :: DofPiCellAll
      double precision, dimension(:,:,:), allocatable :: DofPiCellVEM

c
c     projecteur piK pour une maille k donnee et evalue en un point x
c     matrice ligne de taille nbnodebycell(k) + nbfacebycell(k) telle que piK U(x) = PiK*U  
c     
c            
      double precision, dimension(:), allocatable :: piK       
c      
c     Operateur saut aux faces frac
c      
c     1/|sig|\int_sig (v^+ - v^-) = sum_i=1^n SautbyFaceFrac(numfrac,i)*Vinc(inc) 
c     cote + = cote maille 1 de NumCellbyFace
c     cote - = cote maille 2 
c     inc = NumIncGlobalbyFaceFrac(numfrac,i)
c     n = NbIncGlobalbyFaceFrac((numfrac)
c      
      double precision, dimension(:,:), allocatable :: SautbyFaceFrac     

c
c     Operateur GradTangbyFace calculant le gradient tangentiel nodal constant par face 
c
c     G_l(nf) = sum_is GradTangbyFace(nf,is,l) U(inc)
c
c     is: numero local a la face du noeud, numerote dans l'ordre des noeuds de la face  
c     inc = inc globale du node is de la face 
c      
      double precision, dimension(:,:,:), allocatable :: GradTangbyFace     

      
c     matrice locale ACell de taille (nbnodebycell(k) + nbfacebycell(k)), stockage
c     pour toutes les mailles k (premier indice) 
c
      double precision, dimension(:,:,:,:,:), allocatable :: ACell     
c
c
c
c     operateur sigman by face frac (tractions normale et tangentielles) 
c     sigman_nf(U,:) = sum_{dof,l} Sigmanbyff(nf,dof,l,:) U^dof_l    
c
c
      double precision, dimension(:,:,:,:), allocatable :: Sigmanbyff      
c        
c
c     Sol du systeme Meca 
c
        double precision, dimension(:,:), allocatable :: SolU
        double precision, dimension(:,:), allocatable :: SolU0        
        double precision, dimension(:,:), allocatable :: SolU_nm1
        double precision, dimension(:,:), allocatable :: SolU_nm2        
        double precision, dimension(:,:), allocatable :: SolU_prev


c
c     Indice contact par face frac 
c
        integer, dimension(:), allocatable :: IndContact          
c
c     pour visu 
c
        double precision, dimension(:,:), allocatable :: UCell     
        double precision, dimension(:,:), allocatable :: SolSaut       
        double precision, dimension(:,:), allocatable :: SolLambda 
        double precision, dimension(:,:), allocatable :: solSigma    

c     
c     Sautn de U et Sautt de  U-U0 avec projection pour Nitsche 
c        
      double precision, dimension(:), allocatable :: Sautn_Proj     
      double precision, dimension(:,:), allocatable :: Sautt_Proj
      double precision, dimension(:,:), allocatable :: Sautt_Proj_nm1 
      
c
c
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     WorkSpaces 
c
c
c
c     

        
      dimension EvalXPThetaBetabyff(
     &     NbdofbyCellMax+NbNodebyFaceMax+1,NbDim,NbDim)    
      dimension NumIncGlobalPThetaBeta(NbdofbyCellMax+NbNodebyFaceMax+1)


c      
      dimension EvalXPbeta(
     &     NbdofbyCellMax+NbNodebyFaceMax+1,NbDim,NbDim)

       double precision, dimension(:,:), allocatable :: BetabyFaceFrac
c
c     PbetaUt
c
        dimension sPbetaUt(2)
c
c      
      
      dimension PoidsPoints(2*NbNodebyFaceMax)        
      dimension XPoints(2*NbNodebyFaceMax,NbDim)            
                     
                
      dimension X(NbDim), Vec(NbDim), Diff(NbDim,NbDim)
      dimension Grad(NbDim,NbDim)
      dimension GradEx(NbDim,NbDim)
      dimension deltaU(NbDim)        
c
c     
c
      dimension evalx(NbDim),XXk(NbDim)
c
c      
      dimension BBK(NbdofbyCellMax,NbdofbyCellMax)
      dimension CCK(NbdofbyCellMax,NbdofbyCellMax)

      dimension AId(NbDim,NbDim)
      dimension tensor1(NbDim,NbDim), tensor2(NbDim,NbDim) 
      dimension tensor3(NbDim,NbDim)


      dimension Uex(NbDim),Uex1(NbDim),Uex2(NbDim)
      dimension Vec_s(NbDim), Vec2(NbDim)
      dimension sautex(NbDim)
c
c     produit scalaire (ei,fj) ou ei base canonique, et fj base locale j=(vecnormal,tau1,tau2)
c        
      dimension BaseFace(NbDim,NbDim)
c
c
c      
      dimension Resultat(NbDim)
c
c
      dimension vecsaut_1(NbDim), Vec1(NbDim)
      dimension vecsaut_2(NbDim)
      dimension vecsaut_f(NbDim)
c
c     sorties frac gnuplot
c
        dimension sortiesfrac(10000,15)
        
        dimension numsortiesfrac(10000)
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
ccc   coefficient de frottement

      double precision, dimension(:), allocatable :: rF

      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     cas test a une frac x=0  ou deux fractures en x=0 et z=0 
c
c     pour z=0 il faut permuter y et z
c     ne pas permuter pour les test a une fracture 
c      

            itest = 1                 ! une fracture
c           itest = 2                 ! deux fractures 

c             itest = 3                 ! pas de fractures 
      
ccccccccccccccccccccccc
c
c     Boucle sur les maillages 
c

          i3D = 1
          i2D5 = 0
          itet = 0
          icart = 0
          irand = 0 ! exa random avec decoupage des faces non planes en deux triangles 

            iflecturelitho = 1           
c             iflecturelitho = 0    ! pas de litho ds le fichier maillage pour hexa rand           


          if ( (icart.eq.1).and.(iflecturelitho.eq.0) ) then              
             write(*,*)' il faut ilecturelitho = 1 '
             stop
          endif


          if ( (itet.eq.1).and.(iflecturelitho.eq.0) ) then              
             write(*,*)' il faut ilecturelitho = 1 '
             stop
          endif


          if ( (irand.eq.1).and.(iflecturelitho.eq.1) ) then              
             write(*,*)' il faut ilecturelitho = 0 '
             stop
          endif

          
             
         do imesh=3,3     
c        do imesh=2,2
c        do imesh = 1,10

         write(*,*)
         write(*,*)' mesh ',imesh
         write(*,*)
c
c         
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Lecture du Fichier maillage au format Benchmark 3D FVCA6 
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
c     
      num = 81
c

      if (i3D.eq.1) then

c         open(unit=num,file =
c     &        'cpg.msh',status = 'old')

c         open(unit=num,file =
c     &        'meshfine_2D5.msh',status = 'old')

c         open(unit=num,file =
c     &        'meshmedium_2D5.msh',status = 'old') 

c         open(unit=num,file =
c     &        'meshmedium_raff1_2D5.msh',status = 'old')  ! 122 faces frac   

c         open(unit=num,file =
c     &        'meshcoarse_2D5.msh',status = 'old')

         open(unit=num,file =
     &        'meshraff_2D5.msh',status = 'old') ! 1000 faces frac 


c         open(unit=num,file =
c     &        'mesh244_2D5.msh',status = 'old')    ! 244 faces frac      

         
      else if (i2D5.eq.1) then 



         open(unit=num,file =
     &    '../../../../MAILLAGES/mesh0_2D5.msh',status = 'old')

c         open(unit=num,file =
c     &    '../../../../MAILLAGES/mesh1_2D5.msh',status = 'old')

c         open(unit=num,file =
c     &    '../../../../MAILLAGES/mesh2_2D5.msh',status = 'old')

c         open(unit=num,file =
c     &    '../../../../MAILLAGES/mesh3_2D5.msh',status = 'old')

c         open(unit=num,file =
c     &    '../../../../MAILLAGES/mesh4_2D5.msh',status = 'old')
 
         
cccccccccccc
c
c     maillages tetra avec plans de fractures en x,y,z=0.5  
c      
c
      else if (itet.eq.1) then 
      
      if (imesh.eq.1) then 
         open(unit=num,file =
     &     '../../MAILLAGES/Tetra/tet_1th.msh',status = 'old')
      else if (imesh.eq.2) then
         open(unit=num,file =
     &    '../../MAILLAGES/Tetra/tet_10th.msh',status = 'old')
      else if (imesh.eq.3) then
         open(unit=num,file =
     &    '../../MAILLAGES/Tetra/tet_100th.msh',status = 'old')
      else if (imesh.eq.4) then
         open(unit=num,file =
     &    '../../MAILLAGES/Tetra/tet_200th.msh',status = 'old')
      else if (imesh.eq.5) then
         open(unit=num,file =
     &        '../../MAILLAGES/Tetra/tet1.msh',status = 'old')
      else if (imesh.eq.6) then          
         open(unit=num,file =
     &        '../../MAILLAGES/Tetra/tet2.msh',status = 'old')
      else if (imesh.eq.7) then          
         open(unit=num,file =
     &        '../../MAILLAGES/Tetra/tet3.msh',status = 'old')
      else if (imesh.eq.8) then          
         open(unit=num,file =
     &        '../../MAILLAGES/Tetra/tet4.msh',status = 'old')
      else if (imesh.eq.9) then          
         open(unit=num,file =
     &        '../../MAILLAGES/Tetra/tet5.msh',status = 'old')
      else if (imesh.eq.10) then      
         open(unit=num,file =
     &        '../../MAILLAGES/Tetra/tet6.msh',status = 'old')
      else
         write(*,*)' mesh tet > 10 not found ',imesh
         stop
      endif
      
c
ccccccccccccccccc      
c     Voronoi (pas de plans de frac)
c
c$$$      if (imesh.eq.1) then          
c$$$         open(unit=num,file = '../../Vor/vmesh_1.msh',status = 'old')
c$$$      else if (imesh.eq.2) then           
c$$$         open(unit=num,file = '../../Vor/vmesh_2.msh',status = 'old')
c$$$      else if (imesh.eq.3) then           
c$$$         open(unit=num,file = '../../Vor/vmesh_3.msh',status = 'old')
c$$$      else if (imesh.eq.4) then           
c$$$         open(unit=num,file = '../../Vor/vmesh_4.msh',status = 'old')
c$$$      else if (imesh.eq.5) then           
c$$$         open(unit=num,file = '../../Vor/vmesh_5.msh',status = 'old')
c$$$      else
c$$$         write(*,*)' mesh cart > 5 not found ',imesh
c$$$         stop
c$$$      endif         
c
ccccccccccccc      
c     Kershaw (pas de plan de frac)
c      
c         open(unit=num,file = '../../Kershaw/dkershaw08.msh',status = 'old')             
c
cccccccccccccccc      
c     Cartesien 
c      
c
      else if (icart.eq.1) then 
      
      if (imesh.eq.1) then       
         open(unit=num,file =
     &   '../../MAILLAGES/Cart/hexa_4x4x4.msh',status = 'old')
      else if (imesh.eq.2) then       
         open(unit=num,file =
     &   '../../MAILLAGES/Cart/hexa_8x8x8.msh', status = 'old')
      else if (imesh.eq.3) then       
         open(unit=num,file=
     &    '../../MAILLAGES/Cart/hexa_16x16x16.msh',status = 'old')
      else if (imesh.eq.4) then       
         open(unit=num,file=
     &   '../../MAILLAGES/Cart/hexa_32x32x32.msh',status = 'old')
      else if (imesh.eq.5) then       
         open(unit=num,file =
     &    '../../MAILLAGES/Cart/hexa_64x64x64.msh',status = 'old')
      else if (imesh.eq.6) then       
         open(unit=num,file =
     &  '../../MAILLAGES/Cart/hexa_128x128x128.msh',status = 'old') 
      else
         write(*,*)' mesh cart > 6 not found ',imesh
         stop
      endif
c     
c
      else if (irand.eq.1) then 
      
      if (imesh.eq.1) then       
         open(unit=num,file =
     &    '../../MAILLAGES/Cubic-Cells-transformed/gcube_2x2x2.msh',
     &        status = 'old')
      else if (imesh.eq.2) then       
         open(unit=num,file =
     &     '../../MAILLAGES/Cubic-Cells-transformed/gcube_4x4x4.msh',
     &        status = 'old')
      else if (imesh.eq.3) then       
         open(unit=num,file=
     &    '../../MAILLAGES/Cubic-Cells-transformed/gcube_8x8x8.msh',
     &        status = 'old')
      else if (imesh.eq.4) then       
         open(unit=num,file=
     &    '../../MAILLAGES/Cubic-Cells-transformed/gcube_16x16x16.msh',
     &        status = 'old')
      else if (imesh.eq.5) then       
         open(unit=num,file =
     &    '../../MAILLAGES/Cubic-Cells-transformed/gcube_32x32x32.msh',
     &        status = 'old')
      else if (imesh.eq.6) then       
         open(unit=num,file =
     &    '../../MAILLAGES/Cubic-Cells-transformed/gcube_48x48x48.msh',
     &        status = 'old')
      else if (imesh.eq.7) then       
         open(unit=num,file =
     &     '../../MAILLAGES/Cubic-Cells-transformed/gcube_64x64x64.msh',
     &        status = 'old')          
      else
         write(*,*)' mesh cart > 7 not found ',imesh
         stop
      endif
         
      
      endif
c     
ccccccccccccccccccccccccccccccccccccccc
c     
c     lecture fichier de maillage 3D
c     
ccccccccccccccccccccccccccccccccccccccc
c     
c     
c     lecture des dimensions 
c     

      read(num,*)
      read(num,*)
      read(num,*)
      read(num,*)
      read(num,*)
      read(num,*)
      read(num,*)
      read(num,*)

      read(num,*)
      read(num,*)NbNode
      write(*,*)' NbNode ',NbNode 

      read(num,*)
      read(num,*)NbCell
      write(*,*)' NbCell ',NbCell 

      read(num,*)
      read(num,*)NbFace
      write(*,*)' NbFace ',NbFace

      read(num,*)
      read(num,*)NbArete
      write(*,*)' NbArete ',NbArete
c     
c     allocations 
c     
      allocate(NumNodebyArete(NbArete,2))

      allocate(NbAretebyFace(NbFace))
      allocate(NumAretebyFace(NbFace,NbAretebyFaceMax))

      allocate(NumCellbyFace(NbFace,2))

      allocate(NbFacebyCell(NbCell))
      allocate(NumFacebyCell(NbCell,NbFacebyCellMax))

      allocate(NbNodebyCell(NbCell))
      allocate(NumNodebyCell(NbCell,NbNodebyCellMax))

      allocate(NbNodebyFace(NbFace))
      allocate(NumNodebyFace(NbFace,NbNodebyFaceMax))

      allocate(XNode(NbNode,NbDim))

      allocate(LabelbyFace(NbFace))


      allocate (VecTanFace1(NbFace,NbDim))
      allocate (VecTanFace2(NbFace,NbDim))

c     Diametre de chaque cellule
      allocate (diamK(NbCell))

      allocate(IndRockTypebyCell(NbCell))
c     
c     lecture des tableaux 
c     
c      iflabelfrac = 0
       iflabelfrac = 1

      call lectureBench3D(num,iflabelfrac,iflecturelitho,
     &     NbCell,NbNode,NbFace,NbArete,
     &     XNode,NbNodebyCell,NumNodebyCell,
     &     NbFacebyCell,NumFacebyCell,
     &     NbAretebyFace,NumAretebyFace,
     &     NbNodebyFace,NumNodebyFace,
     &     NumCellbyFace,NumNodebyArete,
     &     IndRockTypebyCell,LabelbyFace)

	do n=1,NbNode
c
c           if (irand.ne.1) then  ! maillage hexa random deja dans (-1,1)^3 

c              XXk(:) = XNode(n,:)
c              XNode(n,:) = 2*XXk(:) - 1.d0  ! mapping (0,1)^3 -> (-1,1)^3 pour maillages tet et cart 
              
c              if (itest.eq.2) then ! on permutte y et z pour maillages tet et cart 
              
            XXk(:) = XNode(n,:)

            XNode(n,2) = XXk(2)/10.d0
           
c           XNode(n,2) = XXk(2)/2.d0           
      
c           XNode(n,1) = XXk(1)/2000.d0 
c           XNode(n,2) = XXk(2)/2000.d0
c           XNode(n,3) = XXk(3)/2000.d0
                 
c              endif
              
c           endif
           
	enddo
c     
c
c
c     On teste l'ordre cyclique des noeuds des faces 
c

      do k=1,NbFace

         n1 = NumNodebyFace(k,1)
         do in = 2,NbNodebyFace(k)+1

            if (in.eq.NbNodebyFace(k)+1) then 
               n2 = NumNodebyFace(k,1)
            else
               n2 = NumNodebyFace(k,in)               
            endif
c
c     on cherche si l'arete n1 n2 existe 
c
            isearch = 0
            do ia = 1,NbAretebyFace(k)

               na = NumAretebyFace(k,ia)
               m1 = NumNodebyArete(na,1)
               m2 = NumNodebyArete(na,2)

               if ((m1.eq.n1).and.(m2.eq.n2)) then
                  isearch = 1
               endif

               if ((m1.eq.n2).and.(m2.eq.n1)) then
                  isearch = 1
               endif

               
            enddo

            if (isearch.eq.1) then
c               write(*,*)' arete n1 n2 exists ',k,in-1,n1,n2               
            else             
               write(*,*)' arete n1 n2 does not exist ',k,in-1,n1,n2
               stop
            endif          
            
            n1 = n2 
            
         enddo
         
         
      enddo
c
c      stop
c     
cccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Identification des faces fractures 
c     
cccccccccccccccccccccccccccccccccccccccccccccc
c


      
      if (iflabelfrac.eq.1) then

cccccccccc
c         
c     On peut enlever les fractures ici 
c         
c         LabelbyFace = 0 ! TMP TMP TMP
cccccccccc
c         
c     
c     lecture label faces 
c     
c     
         allocate(IndFaceFrac(NbFace))
         allocate(NumFaceVersFaceFrac(NbFace))
         allocate(IndFace(NbFace))

         
         NumFaceVersFaceFrac(:) = 0       
c     
c     Comptage 
c     
         NbFaceFrac = 0
         do i=1,NbFace
            if (LabelbyFace(i).eq.1) then 
               NbFaceFrac = NbFaceFrac + 1
            endif
         enddo


         
c     
c     allocation
c     
         allocate(NumFaceFracVersFace(NbFaceFrac))
c     
c     remplissage 
c     
         NbFaceFrac = 0
         do i=1,NbFace
            if (LabelbyFace(i).eq.1) then 
               NbFaceFrac = NbFaceFrac + 1
               NumFaceFracVersFace(NbFaceFrac) = i
               NumFaceVersFaceFrac(i) = NbFaceFrac
               IndFaceFrac(i) = 1
            else 
               IndFaceFrac(i) = 0
            endif
         enddo

         write(*,*)' NbFaceFrac ',NbFaceFrac
c     do i=1,NbFaceFrac
c     write(*,*)' num face frac ',i,NumFaceFracVersFace(i)
c     enddo
c     do i=1,NbFace
c     write(*,*)' IndFaceFrac ',i,IndFaceFrac(i)
c     enddo
         

      else 
c     
cccccccccccccccccccccccccccccccc
c     
c     Identification geometrique des faces fractures 
c     
c     comptage 
c     
         epsi = 1.0d-6
         NbFaceFrac = 0
         do i=1,NbFace
            nsf = NbNodebyFace(i)
            xf = 0.d0
            yf = 0.d0
            zf = 0.d0
            do is=1,nsf
               xf = xf + XNode(NumNodebyFace(i,is),1)
               yf = yf + XNode(NumNodebyFace(i,is),2)
               zf = zf + XNode(NumNodebyFace(i,is),3)
            enddo
            xf = xf/nsf
            yf = yf/nsf
            zf = zf/nsf
            
            erx = abs(xf)       ! ?? une fracture x = 0 ??           
            ery = abs(yf) ! ?? une fracture y = 0 ??
            erz = abs(zf)   ! ?? une fracture z = 0 ??              
c            erz = 1.d0 

            if (      (erz.lt.epsi)
     &           .and.(erx.le.0.5d0)
     &           .and.(ery.le.0.5d0) ) then ! un plan de frac en z = 0
               
               NbFaceFrac = NbFaceFrac + 1
               
               endif


            
         enddo

         write(*,*)' NbFaceFrac ',NbFaceFrac 
c     
c     allocation 
c     
         allocate(IndFaceFrac(NbFace))
         allocate(NumFaceVersFaceFrac(NbFace))
         allocate(IndFace(NbFace))
         allocate(NumFaceFracVersFace(NbFaceFrac))
    
c     
c     remplissage 
c
         
         IndFaceFrac(:) = 0 
         
         epsi = 1.0E-6
         NbFaceFrac = 0
         do i=1,NbFace
            nsf = NbNodebyFace(i)
            xf = 0.d0
            yf = 0.d0
            zf = 0.d0
            do is=1,nsf
               xf = xf + XNode(NumNodebyFace(i,is),1)
               yf = yf + XNode(NumNodebyFace(i,is),2)
               zf = zf + XNode(NumNodebyFace(i,is),3)
            enddo
            xf = xf/nsf
            yf = yf/nsf
            zf = zf/nsf
c     write(*,*)' Face i xf = ',i,xf,yf,zf
            
            erx = abs(xf)       ! ?? une fracture x = 0 ??         
            ery = abs(yf) ! ?? une fracture y = 0 ??            
            erz = abs(zf)   ! ?? une fracture z = 0 ??              
c            erz = 1.d0 


            if (      (erz.lt.epsi)
     &           .and.(erx.le.0.5d0)
     &           .and.(ery.le.0.5d0) ) then ! un plan de frac en z = 0
               
               NbFaceFrac = NbFaceFrac + 1
               NumFaceFracVersFace(NbFaceFrac) = i
               NumFaceVersFaceFrac(i) = NbFaceFrac
               IndFaceFrac(i) = 1
               
            endif

              
            
         enddo

c         stop
c     
ccccccccccccccccccccccccccccccccccccccccccccc
c     
      endif
c
c      do i=1,NbFace
c         write(*,*)' face ',i,IndFaceFrac(i)
c      enddo
      
cccccccccccccccccccccccccccccccccccccccccccc
      
      
      allocate(IndNodeFrac(NbNode))
c
c     on repere les nodes frac 
c
      IndNodeFrac(:) = 0
      do k=1,NbFaceFrac
         numf = NumFaceFracVersFace(k)
         do nk=1,NbNodebyFace(numf)
            n = NumNodebyFace(numf,nk)
            IndNodeFrac(n) = 1
         enddo
      enddo
c
c
c
cccccccccccccccccccccccccccccccccccccc
c     
c     Identification des cotes de la boite supposee rectangle 
c     
cccccccccccccccccccccccccccccccccccccc
c     
      xmax = -1.E+10
      xmin = 1.E+10
      ymax = -1.E+10
      ymin = 1.E+10
      zmax = -1.E+10
      zmin = 1.E+10
      do i=1,NbNode
         xi = XNode(i,1)
         yi = XNode(i,2)
         zi = XNode(i,3)

         xmax = dmax1(xmax,xi)
         ymax = dmax1(ymax,yi)
         zmax = dmax1(zmax,zi)

         xmin = dmin1(xmin,xi)
         ymin = dmin1(ymin,yi)
         zmin = dmin1(zmin,zi)

      enddo
      write(*,*)' xmin xmax ',xmin,xmax
      write(*,*)' ymin ymax ',ymin,ymax
      write(*,*)' zmin zmax ',zmin,zmax

c     stop
c     
c     on identifie les faces de bord pour chaque cote 
c     
c     1 = x = xmin 
c     2 = x = xmax 
c     3 = y = ymin 
c     4 = y = ymax 
c     5 = z = zmin 
c     6 = z = zmax 
c     

      eps = 1.0E-6
      do i=1,NbFace
         IndFace(i) = 0

         if (NumCellbyFace(i,2).eq.-1) then



            
            nsf = NbNodebyFace(i)
            xi = 0.d0
            yi = 0.d0
            zi = 0.d0
            do is=1,nsf
               xi = xi + XNode(NumNodebyFace(i,is),1)
               yi = yi + XNode(NumNodebyFace(i,is),2)
               zi = zi + XNode(NumNodebyFace(i,is),3)
            enddo
            xi = xi/nsf
            yi = yi/nsf
            zi = zi/nsf

            if (IndFaceFrac(i).eq.1) then
               write(*,*)' face frac de bord ',i,xi,yi,zi
               stop
            endif

            if (xi.lt.xmin+eps) then 
               IndFace(i) = 1
            else if (xi.gt.xmax-eps) then 
               IndFace(i) = 2
            else if (yi.gt.ymax-eps) then 
               IndFace(i) = 4
            else if (yi.lt.ymin+eps) then              
               IndFace(i) = 3

            else if (zi.gt.zmax-eps) then 
               IndFace(i) = 6
            else if (zi.lt.zmin+eps) then              
               IndFace(i) = 5
            endif

c     write(*,*)' Face bord ',i,IndFace(i),xi,yi,zi
         endif   
      enddo


      

cccccccccccccccccccccccccccccccccccccc


      allocate(IndIncArete(NbArete))
      allocate(IndIncAreteT(NbArete))      
      allocate(NumIncArete(NbArete))
      allocate(NumIncNode(NbNode))
      allocate(NumIncFace(NbFace))

      allocate(IndDirNodeMeca(NbNode))      
      allocate(IndDirFace(NbFace))
      allocate(IndNeuFace(NbFace))
      allocate(IndDirFaceT(NbFace))
      allocate(IndNeuFaceT(NbFace))      

      allocate(XCell(NbCell,NbDim))
      allocate(XFace(NbFace,NbDim))
      allocate(XFaceFrac(NbFaceFrac,NbDim))



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Choix CL Dirichlet non homogene ou Neumann homogene 
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     Choix CL pour la Meca 
c
c     0 si Neumann 
c
c     > 0 si Dir ou Dir/Neumann 
c      
      IndDirNodeMeca = 0 ! neumann par defaut 
            
      
c
c     Dir selon la normale + neumann homogene tangentiel 
c

      
      do i=1,NbFace         
         if (IndFace(i).eq.5) then ! z = zmin 

            do ik = 1,NbNodebyFace(i)
               n = NumNodebyFace(i,ik)
               if ( IndDirNodeMeca(n).eq.0 ) then 
                  IndDirNodeMeca(n) = 5
               endif
            enddo
            
         endif
      enddo
      

      do i=1,NbFace         
         if (IndFace(i).eq.1) then ! x = xmin 

            do ik = 1,NbNodebyFace(i)
               n = NumNodebyFace(i,ik)
               if ( IndDirNodeMeca(n).eq.0 ) then 
                  IndDirNodeMeca(n) = 1
               endif
            enddo
            
         endif
      enddo



      do i=1,NbFace         
         if (IndFace(i).eq.3) then ! y = ymin 

            do ik = 1,NbNodebyFace(i)
               n = NumNodebyFace(i,ik)
               if ( IndDirNodeMeca(n).eq.0 ) then 
                  IndDirNodeMeca(n) = 3
               endif
            enddo
            
         endif
      enddo


      do i=1,NbFace         
         if (IndFace(i).eq.4) then ! y = ymin 

            do ik = 1,NbNodebyFace(i)
               n = NumNodebyFace(i,ik)
               if ( IndDirNodeMeca(n).eq.0 ) then 
                  IndDirNodeMeca(n) = 4
               endif
            enddo
            
         endif
      enddo        
   
c
c
c      
c      do i=1,NbFace      
c         if (IndFace(i).eq.5) then ! z = zmin 
!     dir en z
c            do ik = 1,NbNodebyFace(i)
c               n = NumNodebyFace(i,ik)
c               if ( IndDirNodeMeca(n).eq.0 ) then                
c                  IndDirNodeMeca(n) = 5
c               endif
c            enddo            
c         endif
c      enddo

c      do i=1,NbFace      
c         if (IndFace(i).eq.6) then ! z = zmin 
!     dir en z
c            do ik = 1,NbNodebyFace(i)
c               n = NumNodebyFace(i,ik)
c               if ( IndDirNodeMeca(n).eq.0 ) then                
c                  IndDirNodeMeca(n) = 5
c               endif
c            enddo            
c         endif
c      enddo      
c         
c
c
c
c      do n=1,NbNode
c         if (IndDirNodeMeca(n).eq.4) then 
c            write(*,*)' node ',n,IndDirNodeMeca(n),XNode(n,:)
c         endif
c      enddo
c      stop
c      
c      
cccccccccccccccccccccccccc
c      
c     Choix CL par face pour P et T 
c
ccccccccccccccccccccccccc
c      
c     Pression 
c      
c     IndDirFace -> 1 Dir, 0 sinon 
c     IndNeuFace -> 1 Neumann homogene, 2 Neumann non homogene 
c     On deduit IndIncArete -> 11 Dir, 10 Neumann
c     On deduit ensuite dans NumerotationInc -> IndIncDir (0,1), IndIncNeu (0,1,2)
c
c     Temperature 
c
c     IndDirFaceT (0,1)
c     IndNeuFaceT -> 1 Neumann homogene fourier, 0 sinon    
c     On deduit IndIncAreteT -> 11 Dir, 10 Neumann
c     On deduit ensuite dans NumerotationInc -> IndIncDirT  (0,1)et IndIncNeuT  (0,1)
c     Si IndIncDirT = 1 -> inc Dir
c     Si IndIncNeuT = 1  -> inc Neumann homogene =flux diffusif Fourier nul 
c      
c      
ccccccccccccccccccccccccc
c     
      IndDirFace(:) = 0
      IndNeuFace(:) = 0
      IndDirFaceT(:) = 0
      IndNeuFaceT(:) = 0      
     
c     
c     Choix par cote 
c     
      do i=1,NbFace

         X(:) = 0.d0 
         do ik = 1,NbNodebyFace(i)
            n = NumNodebyFace(i,ik)
            X(:) = X(:) + XNode(n,:)
         enddo
         X(:) = X(:)/dfloat(NbNodebyFace(i))
         
         if (IndFace(i).eq.1) then 
c            IndDirFace(i) = 1
            IndDirFace(i) = 0
            IndNeuFace(i) = 1
            if ((X(3).le.1050.d0).and.(X(3).ge.950.d0)) then
c            if ((X(2).le.1200.d0).and.(X(2).ge.800.d0)) then

               
               IndNeuFace(i) = 2 ! Neumann non homogene
               write(*,*)' neumann ',i,X(:)
c            endif
            endif

            IndDirFaceT(i) = 1
            IndNeuFaceT(i) = 0            
            
         endif
         if (IndFace(i).eq.2) then
c           IndDirFace(i) = 0            
            IndDirFace(i) = 1
            IndNeuFace(i) = 0

            IndDirFaceT(i) = 0
            IndNeuFaceT(i) = 1    
            
         endif
         if (IndFace(i).eq.3) then 
c            IndDirFace(i) = 1
            IndDirFace(i) = 0
            IndNeuFace(i) = 1

            IndDirFaceT(i) = 0
            IndNeuFaceT(i) = 1    
            
         endif
         if (IndFace(i).eq.4) then 
c           IndDirFace(i) = 1
            IndDirFace(i) = 0
            IndNeuFace(i) = 1
            
            IndDirFaceT(i) = 0
            IndNeuFaceT(i) = 1    
            
          endif
         if (IndFace(i).eq.5) then 
c            IndDirFace(i) = 1 
            IndDirFace(i) = 0
            IndNeuFace(i) = 1

            IndDirFaceT(i) = 0
            IndNeuFaceT(i) = 1    
            
         endif
         if (IndFace(i).eq.6) then
c           IndDirFace(i) = 1  
            IndDirFace(i) = 0
            IndNeuFace(i) = 1

            IndDirFaceT(i) = 0
            IndNeuFaceT(i) = 1    
            
         endif
         
      enddo
      
c     write(*,*)
c     do i=1,NbFace
c     write(*,*)' IndDirFace ',i,IndDirFace(i),IndDirFaceT(i)
c     enddo
c
c      do i=1,NbFace
c         write(*,*)' IndNeuFace ',i,IndNeuFace(i),IndNeuFaceT(i)
c      enddo
c      stop      
cccccccccccccccccccccc
c     Choix CL par arete pour HFV
c     arete DIR si elle appartient a au moins une face DIR 
cccccccccccccccccccccc
c     
      do i=1,NbArete
         IndIncArete(i) = 0
         IndIncAreteT(i) = 0         
      enddo
c     
c     Arete de bord DIR (toutes pas slt les FRAC) > Indice 1 
c     
      do i=1,NbFace
         if (IndDirFace(i).eq.1) then 
            do is=1,NbAretebyFace(i)
               js = NumAretebyFace(i,is)
               IndIncArete(js) = 1
            enddo
         endif
      enddo

      do i=1,NbFace
         if (IndDirFaceT(i).eq.1) then 
            do is=1,NbAretebyFace(i)
               js = NumAretebyFace(i,is)
               IndIncAreteT(js) = 1
            enddo
         endif
      enddo      
c
c     On force les aretes frac du bord 5 a etre en Dirichlet 
c
c      do i=1,NbFace
c         if (IndFaceFrac(i).eq.1) then 
c            do is=1,NbAretebyFace(i)
c               js = NumAretebyFace(i,is)
c               n1 = NumNodebyArete(js,1)
c               n2 = NumNodebyArete(js,2)

c               za = (XNode(n1,3)+XNode(n2,3))/2.d0 

c               if (za.le.-1.d0+1.d-3) then 
               
c                  IndIncArete(js) = 1
c                  write(*,*)' arete dir ',za
c               endif
c            enddo
c         endif
c      enddo      
c
c
c     On met les aretes frac sur une face de bord Dirichlet en mode Dir 
c     Les autres aretes de bord seront en mode Neumann homogene
c
c     Il faudrait completer si on veut mettre du neumann homogene en Fourier par ex !! 
c      
c     
c     arete frac DIR      > Indice 11
c     arete frac non DIR  > Indice 10 
c     
      do i=1,NbFace
         if (IndFaceFrac(i).eq.1) then 
            do is=1,NbAretebyFace(i)
               js = NumAretebyFace(i,is)
               if (IndIncArete(js).eq.1) then 
                  IndIncArete(js) = 11
               else 
                  IndIncArete(js) = 10                  
               endif
            enddo
         endif
      enddo

      do i=1,NbFace
         if (IndFaceFrac(i).eq.1) then 
            do is=1,NbAretebyFace(i)
               js = NumAretebyFace(i,is)
               if (IndIncAreteT(js).eq.1) then 
                  IndIncAreteT(js) = 11
               else 
                  IndIncAreteT(js) = 10                  
               endif
            enddo
         endif
      enddo
      
c     
c     Fin choix des CLs pour pression et temperature 
c
ccccccccccccccccccccccccccccc
c
c     transformation affine pour maillage Cartesien -> hexa a faces planes 
c
ccccccccccccccccccccccccccc
c      
c      do n=1,NbNode
c         xn = XNode(n,1)
c         yn = XNode(n,2)
c         zn = XNode(n,3)

c         XNode(n,1) = (xn + 2.d0*yn)/3.d0 
c         XNode(n,2) = (3.d0*zn + yn)/4.d0 
c         XNode(n,3) = (zn + xn)/2.d0 
         
c      enddo
c
c
cccccccccccccccccccccccccccccccccccc
c     
c     Visu Paraview: maillages matrice et fracture 
c     
      call InitVisuParaview(
     &     iaffich,
     &     NbCell,NbFaceFrac,NbNode,NbFace,XNode,
     &     NbNodebyCell,NumNodebyCell,
     &     NumFaceFracVersFace,NbNodebyFace,NumNodebyFace,
     &     NbFacebyCell,NumFacebyCell,            
     &     Ntetra4,Nhexa8,Nnfaced,Nquad4,Ntria3)

      write(*,*)' cells ',Ntetra4,Nhexa8,Nnfaced
      
c     
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     
c     Numerotation des inconnues fct du schema VAG/HFV et des CLs 
c     
c     COMPTAGE POUR ALLOCATION 
c


cccccccccccccccccccccccccccccccccc      


      NbCV = NbCell + NbFace + 2*NbFaceFrac 
           
      do i=1,NbArete ! comptage arete frac 
         if (IndIncArete(i).eq.11) then 
            NbCV = NbCV + 1           
         endif         
         if (IndIncArete(i).eq.10) then 
            NbCV = NbCV + 1
         endif
      enddo   



      
      write(*,*)
      write(*,*)' Nb Control Volume ',NbCV
      write(*,*)

      
c     
c     ALLOCATION 
c     
      allocate(IndIncDir(NbCV))
      allocate(IndIncNeu(NbCV))
      allocate(IndIncDirT(NbCV))
      allocate(IndIncNeuT(NbCV))      
      
      allocate(NumIncCell(NbCell))
      allocate(NumIncFaceFrac(NbFaceFrac))
      allocate(NumIncCellbyFaceFrac(NbFaceFrac,2))
      allocate(IndIncFaceFrac(NbCV))


c      call NumerotationInc(
c     &     NbCell,NbFaceFrac,NbFace,NbNode,NbArete,
c     &     NbNodebyFace,NumNodebyFace,
c     &     NbAretebyFace,NumAretebyFace,
c     &     IndFaceFrac,NumFaceVersFaceFrac,NumFaceFracVersFace,
c     &     IndIncArete,IndDirFace,IndNeuFace,
c     &     NbCV,
c     &     NumIncCell,NumIncFaceFrac, 
c     &     NumIncArete,NumIncFace,NumIncCellbyFaceFrac,  
c     &     IndIncDir,IndIncNeu,IndIncFaceFrac)
      

      
      call NumerotationInc(
     &   NbCell,NbFaceFrac,NbFace,NbNode,NbArete,
     &   NbNodebyFace,NumNodebyFace,
     &   NbAretebyFace,NumAretebyFace,
     &   IndFaceFrac,NumFaceVersFaceFrac,NumFaceFracVersFace,
     &   IndIncArete,IndDirFace,IndNeuFace,
     &   IndIncAreteT,IndDirFaceT,IndNeuFaceT,      
     &   NbCV,
     &   NumIncCell,NumIncFaceFrac, 
     &   NumIncArete,NumIncFace,NumIncCellbyFaceFrac, 
     &   IndIncDir,IndIncNeu,
     &   IndIncDirT,IndIncNeuT,      
     &   IndIncFaceFrac)



      
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Flux de Darcy pour les schemas VAG ou HFV 
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Permeabilites des mailles 
c     
      allocate(perm(NbCell,NbDim,NbDim))
c  
c      
      do k=1,NbCell

         irtk = IndRockTypebyCell(k)
c     
c     calcul de XCell isobarycentre des noeuds (recalcule dans le schema) 
c     
         nsk = NbNodebyCell(k)
         do m=1,NbDim
            XCell(k,m) = 0.d0
            do i=1,nsk
               is = NumNodebyCell(k,i)
               XCell(k,m) = XCell(k,m) + XNode(is,m)/dfloat(nsk) 
            enddo
         enddo

         xk = XCell(k,1)
         zk = XCell(k,3)

         pi = 4.d0*datan(1.d0)
         tethaf = 10.d0*pi/180.d0
         
         xtk = xk - dtan(tethaf)*(zk-1000.d0) - 500.d0 
         
         perm(k,:,:) = 0.d0 
         do m=1,NbDim
            if (irtk.eq.1) then               
               perm(k,m,m) = Permeabilitem1/Viscof
            else if (irtk.eq.2) then
               
               perm(k,m,m) = Permeabilitem2/Viscof

c               if (dabs(xtk).le.50.d0) then ! zone endommagee maillage meshraff !!!!!
c               if (dabs(xtk).le.25.d0) then ! zone endommagee maillage meshraff !!!!!                  
c                  write(*,*)xk,xtk,zk

c                   perm(k,m,m) = 1000.d0*Permeabilitem2/Viscof ! zone endommagee barriere 10-16
c                   perm(k,m,m) = 100.d0*Permeabilitem2/Viscof ! zone endommagee barriere  10-17 

                  
c               endif
c               
            else if (irtk.eq.3) then
               perm(k,m,m) = Permeabilitem3/Viscof               
            else
               write(*,*)' rock ne 1 2 ou 3 ',k,irtk
               stop
            endif            
         enddo

      enddo

c      stop
c     
c     Permeabilite des faces fractures 
c
c     
c     epaisseur des faces frac 
c     
      allocate(dfbyFaceFrac(NbFaceFrac))

      do i=1,NbFaceFrac
         dfbyFaceFrac(i) = dfmin 
      enddo
c
c
c     
      allocate(permfdf(NbFaceFrac))
      allocate(permfn(NbFaceFrac))

      do i=1,NbFaceFrac

         nf = NumFaceFracVersFace(i)
         k1 = NumCellbyFace(nf,1)
         k2 = NumCellbyFace(nf,2)
         
         permk1 = perm(k1,1,1)
         permk2 = perm(k2,1,1)
         
         rKfn = (permk1 + permk2)/2.d0  
         
c          rKfn = 1.d-6 

              
         rKf = 1.d0/12.d0/Viscof         
         permfdf(i) = rKf 
c         permfdf(i) = rKf*(dfbyFaceFrac(i)**3)
                  

         permfn(i)  = rKfn      ! perm moy des mailles k1 et k2  
         
c         write(*,*)' permfdf permfn ',i,permfdf(i),permfn(i)
      enddo
      
c     
c     
c     
      write(*,*)' NbCell ',NbCell
      write(*,*)' NbFaceFrac ',NbFaceFrac
      write(*,*)' NbCV ',NbCV
c     
c     
c     
      allocate(TransCell(NbCell,
     &     NbInterfacebyCellMax,NbInterfacebyCellMax))

      allocate(TransFaceFrac(NbFaceFrac,
     &     NbInterfacebyFaceFracMax,NbInterfacebyFaceFracMax))


      allocate(TransCellbyFaceFrac(NbFaceFrac,2))



         
      allocate(TransCellFourier(NbCell,
     &     NbInterfacebyCellMax,NbInterfacebyCellMax))

      allocate(TransFaceFracFourier(NbFaceFrac,
     &     NbInterfacebyFaceFracMax,NbInterfacebyFaceFracMax))


      allocate(TransCellbyFaceFracFourier(NbFaceFrac,2))

      
      allocate(GCell(NbCell,NbInterfacebyCellMax,NbDim))

      allocate(GFaceFrac(NbFaceFrac,NbInterfacebyFaceFracMax,NbDim))

      allocate(NbInterfacebyCell(NbCell))
      allocate(NumInterfacebyCell(NbCell,NbInterfacebyCellMax))

      allocate(NbInterfacebyFaceFrac(NbFaceFrac))
      allocate(NumInterfacebyFaceFrac(NbFaceFrac,
     &     NbInterfacebyFaceFracMax))


      allocate(XInc(NbCV,NbDim))
      allocate(VolCell(NbCell))
      allocate(VolFaceFrac(NbFaceFrac))


c     
c     schema HFV 
c     
      call HFV3D (CondThermique, 
     &     NbCell,NbNode,NbFace,NbArete,
     &     NbFaceFrac,NbCV, 
     &     XNode,NbNodebyCell,NumNodebyCell,
     &     NbFacebyCell,NumFacebyCell,
     &     NbAretebyFace,NumAretebyFace,
     &     NbNodebyFace,NumNodebyFace,
     &     NumCellbyFace,NumNodebyArete,
     &     NumFaceFracVersFace,NumFaceVersFaceFrac,IndFaceFrac,
     &     NumIncCell,NumIncFaceFrac,NumIncFace,
     &     IndIncArete,NumIncArete,NumIncCellbyFaceFrac,  
     &     perm,permfdf,permfn,dfbyFaceFrac, 
     &     XCell,XFace,XFaceFrac,XInc,
     &     VolCell,VolFaceFrac,
     &     NbInterfacebyCell,NumInterfacebyCell,
     &     NbInterfacebyFaceFrac,NumInterfacebyFaceFrac,
     &     TransCell,TransFaceFrac,TransCellbyFaceFrac,
     &     TransCellFourier,TransFaceFracFourier,
     &     TransCellbyFaceFracFourier,
     &     GCell,GFaceFrac)


     

      do inc=1,NbCV
         if (IndIncNeu(inc).eq.2) then
            write(*,*)' inc neumann non homog Pression ',inc,XInc(inc,:) 
         endif
      enddo
c      stop

      do inc=1,NbCV
         if (IndIncDir(inc).eq.1) then
            write(*,*)' Dir P ',inc,XInc(inc,:)
         endif
      enddo
      write(*,*)
      do inc=1,NbCV
         if (IndIncDirT(inc).eq.1) then
            write(*,*)' Dir T ',inc,XInc(inc,:)
         endif
      enddo      
      write(*,*)
      do inc=1,NbCV
         if (IndIncNeu(inc).eq.1) then
            if (dabs(XInc(inc,2)-2.5d0).le.1.d-3) then 
               write(*,*)' Neu P ',inc,XInc(inc,:)
            endif
         endif
      enddo
      write(*,*)
      do inc=1,NbCV
         if (IndIncNeuT(inc).eq.1) then
            if (dabs(XInc(inc,2)-2.5d0).le.1.d-3) then 
               write(*,*)' Neu T ',inc,XInc(inc,:)
            endif
         endif
      enddo    
      
c      call CPU_TIME(time_scheme2)
c
c      write(*,*)' CPU scheme ',
c     &     time_scheme2-time_scheme1


c$$$      do i=1,NbFaceFrac
c$$$c         write(*,*)' TKsig TKsig ',i,TransCellbyFaceFrac(i,:)
c$$$
c$$$         incf = NumIncFaceFrac(i)
c$$$
c$$$         
c$$$         numf = NumFaceFracVersFace(i)
c$$$
c$$$         k1 = NumCellbyFace(numf,1)
c$$$         k2 = NumCellbyFace(numf,2)
c$$$         inck1 = NumIncCell(k1)
c$$$         inck2 = NumIncCell(k2)
c$$$         
c$$$         incfk1 = NumIncCellbyFaceFrac(i,1)
c$$$         incfk2 = NumIncCellbyFaceFrac(i,2)
c$$$
c$$$         write(*,*)' Xf ',XInc(incf,:)
c$$$         write(*,*)' Xf ',( XInc(incfk1,:)+XInc(incfk2,:) )/2.d0          
c$$$         write(*,*)' Xfk1 ',XInc(incfk1,:)
c$$$         write(*,*)' Xfk2 ',XInc(incfk2,:)
c$$$         write(*,*)' Xk1 ',XInc(inck1,:)
c$$$         write(*,*)' Xk2 ',XInc(inck2,:)         
c$$$         write(*,*)
c$$$         
c$$$      enddo


c$$$      do k=1,NbCell
c$$$
c$$$         inck = NumIncCell(k)
c$$$         write(*,*)' Xk ',XInc(inck,:)
c$$$         
c$$$         do ik=1,NbInterfacebyCell(k)
c$$$
c$$$            incik = NumInterfacebyCell(k,ik)
c$$$            write(*,*)' Xik ',ik,XInc(incik,:)            
c$$$            
c$$$         enddo
c$$$         write(*,*)
c$$$
c$$$         
c$$$      enddo
      

c      stop
c     
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     FIN COMPTAGE ET STRUCTURE CREUSE DE LA JACOBIENNE POUR LE FLOW 
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c
c
c      
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Structure de donnees et inconnues de la MECANIQUE DU CONTACT  
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c
c     On CHOISIT ICI OU ON MET LES BULLES 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
      allocate(IndBulleNumCellbyFace(NbFace,2))

      
      IndBulleNumCellbyFace(:,:) = 0 
c      IndBulleNumCellbyFace(:,:) = 1 

      
      do i=1,NbFace
c
c     on met des bulles sur ttes les face pour commencer 
c

         
c          IndBulleNumCellbyFace(i,:) = 1         
c         IndBulleNumCellbyFace(i,:) = 0         


         if (NumCellbyFace(i,2).le.0) then
c            IndBulleNumCellbyFace(i,:) = 0                  
         else

c             IndBulleNumCellbyFace(i,:) = 1

            if (IndFaceFrac(i).eq.1) then


c               ifrac = NumFaceVersFaceFrac(i)
c               write(*,*)' ifrac ',i,ifrac

               if (iMecaVBulle.eq.1) then 
               
                 IndBulleNumCellbyFace(i,1) = 1 ! bulle du cote maille 1 sur la face frac
c               IndBulleNumCellbyFace(i,2) = 1 ! bulle du cote maille 2 sur la face frac

              endif

              
            endif
            
            
c            yf = XFace(i,2)                        
c            if (dabs(yf-0.5d0).le.1.d-6) then 
            
c               IndBulleNumCellbyFace(i,:) = 1
c               write(*,*)' face bulle ',i,yf
c            endif
               
         endif
         
      enddo
c      stop
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c     bulles aux faces non bord dans les barrieres 
c      
c$$$      do k=1,NbCell
c$$$
c$$$         irtk = IndRockTypebyCell(k)
c$$$
c$$$
c$$$         if (irtk.eq.2) then 
c$$$            do ik=1,NbFacebyCell(k)
c$$$
c$$$               nf = NumFacebyCell(k,ik)
c$$$
c$$$               if (NumCellbyFace(nf,2).gt.0) then
c$$$
c$$$                  IndBulleNumCellbyFace(nf,:) = 1               
c$$$               
c$$$               endif
c$$$               
c$$$            enddo
c$$$         endif
c$$$         
c$$$      enddo
c$$$      
c$$$      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      
cccccccccccccccccccccccccccccccccccc      
c
c     connectivite: ensemble des mailles connectees a chaque node 
c
cccccccccccccccccccccccccccccccccccc
c      
      allocate(NbCellbyNode(NbNode))
      allocate(NumCellbyNode(NbNode,NbCellbyNodeMax))

       do n=1,NbNode
          NbCellbyNode(n) = 0
       enddo
       
       do k=1,NbCell
          do in=1,NbNodebyCell(k)
             n = NumNodebyCell(k,in)
             NbCellbyNode(n) = NbCellbyNode(n) + 1
             if (NbCellbyNode(n).gt.NbCellbyNodeMax) then
                write(*,*)' redim NbCellbyNodeMax ',
     &               NbCellbyNode(n),NbCellbyNodeMax
                stop
             endif
             NumCellbyNode(n,NbCellbyNode(n)) = k 
          enddo             
       enddo

c       do n=1,NbNode
c          write(*,*)'CellbyNode ',n,NbCellbyNode(n),
c     &   NumCellbyNode(n,1:NbCellbyNode(n))
c       enddo
c       stop
c       
c             
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     classes d'equivalences pour les ensembles cellbynode 
c
c     deux mailles sont dans la meme classe si elles partagent une face non frac 
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       
       allocate(LabelCellbyNode(NbNode,NbCellbyNodeMax))
       allocate(NbLabelbyNode(NbNode))

       call Calcul_LabelCellbyNode(NbNode,NbCell,NbFace,
     &      NbCellbyNode,NumCellbyNode,
     &      NbFacebyCell,NumFacebyCell,NumCellbyFace,
     &      IndFaceFrac,IndNodeFrac,      
     &      NbLabelbyNode,LabelCellbyNode) 


c       do n=1,NbNode
c          write(*,*)'Nblabel ',n,IndNodeFrac(n),
c     &      NbLabelbyNode(n),LabelCellbyNode(n,1:NbCellbyNode(n))
c       enddo
c
c
c       stop
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       
c     Numerotation locale -> globale des inconnues Meca (Ks + Ksigma) 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       
        allocate(NbIncGlobalbyCell(NbCell))
        allocate(NumIncGlobalbyCell(NbCell,NbIncbyCellMax))       
c
c     numero local node puis face 
c        
        allocate(NumNodeFaceLocalbyCell(NbCell,NbIncbyCellMax))       



        call NumIncMeca(NbNode,NbCell,NbFace,NbFaceFrac, 
     &     NbCellbyNode,NumCellbyNode,
     &     NbNodebyCell,NumNodebyCell,
     &     NbFacebyCell,NumFacebyCell,NumCellbyFace,
     &     IndFaceFrac,IndBulleNumCellbyFace,
     &     NbLabelbyNode,LabelCellbyNode,
     &     NumIncGlobalbyCell,NbIncGlobalbyCell,
     &     NumNodeFaceLocalbyCell,      
     &     NbdofMeca)        
c
c      do k=1,NbCell         
c         nk = NbIncGlobalbyCell(k)
c         nkk = NbNodebyCell(k)         
c         write(*,*)' cell k ',k,XCell(k,:)
c         write(*,*)' nb inc meca by cell ',k,nk,nk-nkk         
c         do ik = 1,nkk 
c            n = NumNodebyCell(k,ik)
c            write(*,*)' num inc Ks et numero node ',
c     &           NumIncGlobalbyCell(k,ik),
c     &           NumNodeFaceLocalbyCell(k,ik),n,         
c     &           XNode(n,:)
c         enddo         
c         do ik=nkk+1,nk
c            n = NumNodeFaceLocalbyCell(k,ik)
c            write(*,*)' num inc Ksigma et numero face '
c     &           ,NumIncGlobalbyCell(k,ik),n
c     &           ,XFace(n,:)       
c         enddo         
c         write(*,*)
c      enddo
c      stop
c

        if (iMecaVBulle.eq.1) then 
c
c     pour chaque face frac on choisit une inconnue bulle associee 
c        
        do k=1,NbFaceFrac
           nf = NumFaceFracVersFace(k)
           k1 = NumCellbyFace(nf,1)
           k2 = NumCellbyFace(nf,2)


c           write(*,*)' ind bulle ',k,nf,IndBulleNumCellbyFace(nf,:)
           
           if (IndBulleNumCellbyFace(nf,1).eq.1) then
              kf = k1
           else if (IndBulleNumCellbyFace(nf,2).eq.1) then
              kf = k2
           else
              write(*,*)' pas de bulle pour la face frac ',
     &                   k,nf,XFace(nf,:)
              stop
           endif

           ind = -1
           do i = NbNodebyCell(kf)+1,NbIncGlobalbyCell(kf)
              inc = NumIncGlobalbyCell(kf,i)
              ilocalface = NumNodeFaceLocalbyCell(kf,i)
              iface = NumFacebyCell(kf,ilocalface)
              if (iface.eq.nf) then
                 ind = 1
                 incbulle = inc 
              endif
           enddo

           if (ind.eq.-1) then
              write(*,*)' on ne trouve pas l inc bulle ',k,nf
              stop
           else
c              write(*,*)' inc bulle by face frac ',k,incbulle
           endif
           
        enddo

        endif

c        stop
c     
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Geometrie pour VEM 
c
c     Centre de gravite face
c     Surface face
c     Vecteur normal sortant nKsigma 
c
c     Longueur arete 
c     Vecteur normal sortant nSigma_arete
c
ccccc      
c      
c     VolCell, XFace, XCell sont calcules dans les schemas HFV ou VAG
c     Warning: XCell calcule dans HFV ou VAG n'est pas necessaire le centre de gravite
c     XFace est le centre de gravite si schema HFV, mais pas necessairement si schema VAG 
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      
      allocate(XFaceCG(NbFace,NbDim))
      allocate(SurfaceFace(NbFace))
      allocate(VecNormalKSigma(NbCell,NbFacebyCellMax,NbDim))
      allocate(PoidsXFaceCG(NbFace,NbNodebyFaceMax))
      allocate(VecNormalbyFace(NbFace,NbDim))
      
      allocate(SizeArete(NbArete))      
      allocate(VecNormalSigmaArete(NbFace,NbAretebyFaceMax,NbDim))


      allocate(XCellCG(NbCell,NbDim))
      allocate(PoidsXCellCG(NbCell,NbNodebyCellMax))

      
      
      call Geometry(NbCell,NbNode,NbFace,NbArete,
     &             XNode,NbNodebyCell,NumNodebyCell,
     &             NbFacebyCell,NumFacebyCell,
     &             NbAretebyFace,NumAretebyFace,
     &             NbNodebyFace,NumNodebyFace,
     &             NumCellbyFace,NumNodebyArete,
     &             XCell,VolCell,
     &             XCellCG,PoidsXCellCG,       
     &             XFaceCG,PoidsXFaceCG,SurfaceFace,SizeArete,
     &             VecNormalKSigma,VecNormalSigmaArete,
     &             VecNormalbyFace)

ccccc
ccccc Calcul du dimaetre de chaque cellule
ccccc

c
c     TMP 2D !!!!!!!!!!!!!!!! calcul du diametre xz 
c
      
      call ComputeDiambyCell(NbCell,NbNode,
     &             XNode,NbNodebyCell,
     &             NumNodebyCell,
     &             diamK)

      diam_max = 0.d0
      diam_min = 1.d+10
      do k=1,NbCell
         diam_max = dmax1(diamK(k),diam_max)
         diam_min = dmin1(diamK(k),diam_min)         
      enddo
      write(*,*)
      write(*,*)' diam cell max ',diam_max
      write(*,*)' diam cell min ',diam_min       
      write(*,*)' warning  TMP 2D calcul du diametre xz '
      write(*,*) 
c    
c       stop
c
c
c     calcul d'une base orthonormale avec les vecteurs unitaires 1 et 2 tangents a la face
c     suivi du vecteur normal a la face VecNormalbyFace
c
c     on calcule aussi le facteur OrientationbyFaceFrac qui donne une orientation intrinseque de la normale
c
c     VecNormalbyFace(iface,:)*OrientationbyFaceFrac(ifrac) 
c      
c     (VecNormalbyFace a une orientation locale) 
c
c
c      
      allocate(OrientationbyFaceFrac(NbFaceFrac))
      
      call RepereOrthoByFace (
     &     NbFace,NbFaceFrac,VecNormalbyFace,
     &     VecTanFace1,VecTanFace2,OrientationbyFaceFrac,
     &     NumFaceFracVersFace)

c
c
c      do ifrac = 1,NbFaceFrac
c         write(*,*)' orientation ',ifrac,OrientationbyFaceFrac(ifrac)         
c      enddo
c      stop
c      
c     Tester la subroutine "RepereOrthoByFace"c
c
c
c$$$      do iface = 1,NbFace
c$$$c
c$$$         if (IndFaceFrac(iface).eq.1) then
c$$$
c$$$            ifrac = NumFaceVersFaceFrac(iface)
c$$$            
c$$$         write(*,*)'ifrac',ifrac
c$$$         write(*,*)'VecTan1',VecTanFace1(iface,:)
c$$$         write(*,*)'VecTan2,',VecTanFace2(iface,:)
c$$$         write(*,*)'Vecnorm ',VecNormalbyFace(iface,:)
c$$$         write(*,*)
c$$$
c$$$c         call ProdVectoriel(VecTanFace1(iface,:),VecTanFace2(iface,:),
c$$$c     &                       Resultat)
c$$$c         write(*,*)'erreur1',
c$$$c     &     vecnorme(Resultat(:) - VecNormalbyFace(iface,:))
c$$$         
c$$$         write(*,*)'prodscal',
c$$$     &   prodscal(VecTanFace1(iface,:),VecNormalbyFace(iface,:))
c$$$         write(*,*)'prodscal',
c$$$     &   prodscal(VecTanFace2(iface,:),VecNormalbyFace(iface,:))
c$$$         write(*,*)'prodsca3',
c$$$     &    prodscal(VecTanFace1(iface,:),VecTanFace2(iface,:))
c$$$         write(*,*)
c$$$
c$$$         endif
c$$$         
c$$$      enddo
c$$$      STOP

      
c      do i=1,NbFace
c         write(*,*)' XCG face  ',i,XFaceCG(i,:)
c         write(*,*)' surface face ',i,SurfaceFace(i)        
c         do ia=1,NbAretebyFace(i)
c            na = NumAretebyFace(i,ia)
c            n1 = NumNodebyArete(na,1)
c            n2 = NumNodebyArete(na,2)
c            write(*,*)' longueur arete ',na,SizeArete(na)
c            write(*,*)' X node 1 ',n1,XNode(n1,:)
c            write(*,*)' X node 2 ',n2,XNode(n2,:)
            
c            write(*,*)' vec normal sortant arete ',
c     &                       VecNormalSigmaArete(i,ia,:)
c            write(*,*)' Poids CG face ',
c     &           PoidsXFaceCG(i,1:NbNodebyFace(i))
c            write(*,*)
c         enddo
c         write(*,*)
c      enddo


c      vol = 0.d0 
c      do k=1,NbCell
c
c         write(*,*)' Vol maille k ',k,VolCell(k)
c
c         vol = vol + VolCell(k)
c         
c         write(*,*)' Xk maille ',k,XCell(k,:)
c         do ik = 1,NbFacebyCell(k)
c            numf = NumFacebyCell(k,ik)
c            write(*,*)' XCG face ',numf,XFaceCG(numf,:)
c            write(*,*)' vec normal Ksigma ',VecNormalKSigma(k,ik,:)
c            write(*,*)
c         enddo
c         write(*,*)
c      enddo
c      write(*,*)' vol ',vol 
c
c      
ccccccccccccccccccccccccccccccccccccccccc
c
c      
ccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Gradient maille All  (ici on met toutes les faces 'All')  
c
cccccccccccccccccccccccccccccccccccccccccccccccccc      
c
      write(*,*)' test grad cell VEM '

      allocate(GradCellAll(NbCell,NbdofbyCellMax,NbDim))
      
      allocate (piK(NbdofbyCellMax))

      allocate(DofPiCellAll(NbCell,NbNodebyCellMax,NbdofbyCellMax))

      
      
c      call ComputegradvK(NbCell,NbFace,NbNode,
c     &     SurfaceFace,XNode,XFaceCG,
c     &     NbNodebyFace,NumNodebyFace,
c     &     NbNodebyCell,NumNodebyCell,       
c     &     NbFacebyCell,NumFacebyCell,
c     &     PoidsXFaceCG,
c     &     VecNormalKSigma,
c     &     VolCell,
c     &     GradCellAll)

     
      call ComputegradvKnonplanar(
     &     NbCell,NbFace,NbNode,
     &     SurfaceFace,XNode,NbArete,
     &     NbNodebyFace,NumNodebyFace,
     &     NbNodebyCell,NumNodebyCell,       
     &     NbFacebyCell,NumFacebyCell,
     &     PoidsXFaceCG,
     &     VecNormalKSigma,
     &     VolCell,
     &     GradCellAll,
     &     NbAretebyFace,
     &     NumAretebyFace,
     &     XCellCG,
     &     NumNodebyArete)    
c
c     
c      
cccccccccccccccccccccccccccccccccccccccccccccccccc      
c
c     Calcul de l'operateur DofPiCellAll: (ici on met toutes les faces 'All')
c
c     matrice rectangulaire telle que DofPiCell*U = dof PiK(U)
c      
c     cf on en met pas les lignes nulles pour les dof de face 
c     
c      
ccccccccccccccccccccccccccccccccccccccccccccccccc     
c
cc Note: DofPiCell(ncell,:,:) est une matrice de "NbnodebyCellMax" ligne
cc et "NbdofbyCellMax" colonne
cc tels ques le ligne "s" de "DofPiCell" est: piK(nCell,i)(xs), i = 1,...,NbnodebyCell(nCell) + NbFacebyCell(nCell)
cc avec "xs" est le s \E8me noeud de K(nCell)
        
      
      do nCell=1,NbCell   

         nnK = NbnodebyCell(nCell)

         do nnode = 1,nnK
            
            n = NumNodebyCell(nCell,nnode)
            evalx = XNode(n,:)

           
            call ComputepiK(NbCell,nCell,NbNode, 
     &           NbNodebyCell,NumNodebyCell,           
     &           NbFacebyCell,XNode,
     &           XCellCG,PoidsXCellCG,
     &           GradCellAll,
     &           evalx,
     &           piK)


            DofPiCellAll(nCell,nnode,:) = piK(:)
            
         enddo
         
      enddo

      

cccccccccccccccccccccccccccccccccccccccccccccccc
c
c     On extrait les sous operateur en fct des bulles selectionnees  
c     par IndBulleNumCellbyFace 
c      
cccccccccccccccccccccccccccccccccccccccccccccccc

      allocate(GradCellVEM(NbCell,NbdofbyCellMax,NbDim))      

      allocate(DofPiCellVEM(NbCell,NbNodebyCellMax,NbdofbyCellMax))


      
      do k=1,NbCell


         nnk = NbNodebyCell(k)
         
         GradCellVEM(k,1:nnk,:) = GradCellAll(k,1:nnk,:)

         DofPiCellVEM(k,1:nnk,1:nnk) = DofPiCellAll(k,1:nnk,1:nnk)

         
         do i=nnk+1,NbIncGlobalbyCell(k)

            j = NumNodeFaceLocalbyCell(k,i) + nnk 
            
            GradCellVEM(k,i,:) = GradCellAll(k,j,:)

            DofPiCellVEM(k,1:nnk,i) = DofPiCellAll(k,1:nnk,j)
            
         enddo

         
      enddo
         
      deallocate(GradCellAll)
      deallocate(DofPiCellAll)
      

      
      
cccccccccccccccccccccccccccccccccccccccccccccccccccc      
c
c     Matrice locale ACell = |K|(gradK,gradK) + alpha*diam(K)* t(I-DofPiCellVEM)*(I-DofPiCellVEM) 
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      
      allocate(ACell(NbCell,NbdofbyCellMax,NbdofbyCellMax,NbDim,NbDim))      
c
c
c     Matrice identite 3 X 3
c      
         AId(:,:) = 0.d0
         do i=1,NbDim
            AId(i,i) = 1.d0
         enddo
c
c

         rmu = f_Lame_mu()

         rlambda = f_Lame_lambda()

         
         write(*,*)' rmu ',rmu
         write(*,*)' rlambda ',rlambda
 

         
cccccccccccccccccccccc      
      do k=1,NbCell
cccccccccccccccccccccc
c
c         write(*,*)' cell ',k
c         
c
c     BBK = I - DofPiCellVEM
c         
         BBK = 0.d0 
         do i=1,NbNodebyCell(k)
            do j=1,NbIncGlobalbyCell(k)
               BBK(i,j) = - DofPiCellVEM(k,i,j)                
            enddo
            BBK(i,i) = 1.d0 + BBK(i,i)
         enddo
         do i=NbNodebyCell(k)+1,NbIncGlobalbyCell(k)
            BBK(i,i) = 1.d0
         enddo
         
c
c     CCK = t^BBK*BBK 
c         
         do i=1,NbIncGlobalbyCell(k)
            do j=1,NbIncGlobalbyCell(k)
               s = 0.d0               
               do l=1,NbIncGlobalbyCell(k)
                  s = s + BBK(l,i)*BBK(l,j)
               enddo
               CCK(i,j) = s 
            enddo
         enddo


cccccccccccc
         
         alpha = 2.d0*rmu + rlambda ! coefficient de stabilisation 

        
         do ialp = 1,NbIncGlobalbyCell(k)
         do ibet = 1,NbIncGlobalbyCell(k)


            ss = prodscal(GradCellVEM(k,ialp,:)
     &                          ,GradCellVEM(k,ibet,:))

            call prodtensor(GradCellVEM(k,ialp,:)
     &              ,GradCellVEM(k,ibet,:),tensor1)
               
            call prodtensor(GradCellVEM(k,ibet,:)
     &              ,GradCellVEM(k,ialp,:),tensor2)   
               

c           ACell(k,ialp,ibet,:,:) =  VolCell(k)*ss*AId(:,:)

               
c            ACell(k,ialp,ibet,:,:) =
c     &           0.5d0*rmu*VolCell(k)*(
c     &           2.d0*ss*AId(:,:)
c     &           + tensor1(:,:) + tensor2(:,:)   )             
c     &           + VolCell(k)*rlambda*tensor1(:,:) 
            
                
c            ACell(k,ialp,ibet,:,:) =
c     &           0.5d0*rmu*VolCell(k)*(
c     &           2.d0*ss*AId(:,:)
c     &           + tensor1(:,:) + tensor2(:,:)   )             
c     &           + VolCell(k)*rlambda*(tensor1(:,:) + tensor2(:,:))/2.d0
            
               
            ACell(k,ialp,ibet,:,:) =
     &           0.5d0*rmu*VolCell(k)*(
     &           2.d0*ss*AId(:,:)
     &           + 2.d0*tensor2(:,:) )             
     &           + VolCell(k)*rlambda*tensor1(:,:) 
            
               
               
            do i=1,NbDim
               
               ACell(k,ialp,ibet,i,i) =
     &              ACell(k,ialp,ibet,i,i)
     &              +   alpha*diamK(k)*CCK(ialp,ibet)
               
            enddo

               
c           write(*,*)' ACell ',k,ibet,ialp,Acell(k,ibet,ialp,:,:)

         enddo
         enddo
cccccccccccccccc


      enddo ! fin boucle sur les mailles 

      


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
c
c     Operateur saut aux faces frac: SautbyFaceFrac
c
c     1/|sig|\int_sig (v^+ - v^-) = sum_i=1^n SautbyFaceFrac(numfrac,i)*Vinc(inc) 
c     cote + = cote maille 1 de NumCellbyFace
c     inc = NumIncGlobalbyFaceFrac(numfrac,i)
c     n = NbIncGlobalbyFaceFrac((numfrac)
c
c        
c     Numerotation globale des inc de l'op de saut pour l'ordre local suivant:
c        - nodes de la face cote maille 1 
c        - nodes de la face cote maille 2 
c        - face cote maille 1 si bulle ( cf si IndBulleNumCellbyface(n,1) = 1 )
c        - face cote maille 2 si bulle ( cf si IndBulleNumCellbyface(n,2) = 1 )
c
c<
        allocate(NbIncGlobalbyFaceFrac(NbFaceFrac))
        allocate(NumIncGlobalbyFaceFrac(NbFaceFrac,2*NbNodebyFaceMax+2))            

        allocate(SautbyFaceFrac(NbFaceFrac,2*NbNodebyFaceMax+2))
        

       call ComputeSautbyFaceFrac(
     &     NbFaceFrac,NbFace,NbCell,
     &     NbNodebyFace,NumNodebyFace,
     &     NbNodebyCell,NumNodebyCell,       
     &     NbIncGlobalbyCell,
     &     NumIncGlobalbyCell,
     &     NumNodeFaceLocalbyCell,
     &     NumFaceFracVersFace,
     &     NumFacebyCell,NumCellbyFace,
     &     IndBulleNumCellbyFace,            
     &     PoidsXFaceCG,      
     &     NbIncGlobalbyFaceFrac,
     &     NumIncGlobalbyFaceFrac,
     &     SautbyFaceFrac) 

       

cccccccccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccccccccccccccccccc
c
c     Operateur SigmanbyFaceFrac 
c
c      
      allocate(Sigmanbyff(NbFaceFrac,NbdofbyCellMax,Nbdim,NbDim))
       
      call SigmanbyFaceFrac(
     &     NbCell,NbFace,NbFaceFrac,
     &     NumCellbyFace,
     &     NumFaceFracVersFace,      
     &     GradCellVEM,
     &     VecNormalbyFace,
     &     VecTanFace1,
     &     VecTanFace2,           
     &     NbIncGlobalbyCell,
     &     Sigmanbyff)       
c
c     test operateur SigmanbyFaceFrac (composante normale) 
c
c
c       do i=1,NbFaceFrac
c          nf = NumFaceFracVersFace(i)
c          k1 = NumCellbyFace(nf,1)
c          k2 = NumCellbyFace(nf,2)
c
c          write(*,*)' vec n1 ',VecNormalbyFace(nf,:)          
c          write(*,*)' vec t1 ',VecTanFace1(nf,:)
c          write(*,*)' vec t2 ',VecTanFace2(nf,:)
c
c          Vec2(:) = 0.d0 
c          do ik = 1,NbIncGlobalbyCell(k1)
c             inc = NumIncGlobalbyCell(k1,ik)
c             do l=1,NbDim
c                Vec2(:) = Vec2(:) + Sigmanbyff(i,ik,l,:)*SolUEx(inc,l)
c             enddo             
c          enddo
c
c          write(*,*)' sigman ',i,Vec2(:)
c          write(*,*)
c         
c       enddo
c
c       stop
c
c
cccccccccccccccccccccccccccccccc
c
c     Gradient tangentiel nodal constant par face 
c
       allocate(GradTangbyFace(NbFace,NbNodebyFaceMax,NbDim))        


      
      call GradTangentbyFace(NbNode,NbFace,
     &     XNode,
     &     NbNodebyFace,NumNodebyFace,
     &     XFaceCG,SurfaceFace,
     &     GradTangbyFace)

c
c      
c      do k=1,NbFace
c         Vec = 0.d0 
c         nsf = NbNodebyFace(k)
c         do is = 1,nsf
c            ns = NumNodebyFace(k,is)
c            X(:) = XNode(ns,:)
c            fs = f_test(X)
c            Vec(:) = Vec(:) + GradTangbyFace(k,is,:)*fs            
c         enddo
c         write(*,*)' grad tang ',k,Vec(:)
         
c      enddo
c      stop
c
cccccccccccccccccccccccccccccc
c      
c     Test Pthetabeta sur fonction affine f_SolU 
c     On utilise le calcul de Sigmann avec Sigmanbyff pour valider le test 
c
c
c$$$
c$$$      theta = 1.d0
c$$$      betan  = 1.d0
c$$$      betat  = 1.d0       
c$$$           
c$$$      
c$$$      do k = 1,NbCell
c$$$
c$$$         XXk(:) = XCell(k,:) ! sert a positionner le sous domaine 
c$$$         
c$$$c         write(*,*)' cell ',k
c$$$         
c$$$         do i=1,NbNodebyCell(k)
c$$$            inc = NumIncGlobalbyCell(k,i)
c$$$            
c$$$            n = NumNodebyCell(k,i) ! numero global du node
c$$$
c$$$            call f_SolU(XNode(n,:),XXk,Vec)
c$$$            
c$$$            SolUEx(inc,:) =  Vec(:)
c$$$
c$$$
c$$$         enddo
c$$$
c$$$         do i=NbNodebyCell(k)+1,NbIncGlobalbyCell(k)
c$$$            inc = NumIncGlobalbyCell(k,i)           
c$$$
c$$$            nk = NumNodeFaceLocalbyCell(k,i) ! numero local de la face 
c$$$            n = NumFacebyCell(k,nk)
c$$$    
c$$$            SolUEx(inc,:) = 0.d0  
c$$$           
c$$$         enddo            
c$$$c         write(*,*)
c$$$      enddo
c$$$
c$$$      
c$$$      
c$$$      
c$$$      do i=1,NbFaceFrac
c$$$
c$$$         numfacefrac = i 
c$$$       
c$$$         nf = NumFaceFracVersFace(i)
c$$$         k1 = NumCellbyFace(nf,1)
c$$$         k2 = NumCellbyFace(nf,2)
c$$$
c$$$
c$$$         s = 0.d0 
c$$$         do ik = 1,NbIncGlobalbyCell(k1)
c$$$            inc = NumIncGlobalbyCell(k1,ik)
c$$$   
c$$$            s = s + prodscal(Sigmanbyff(i,ik,:,1),SolUEx(inc,:))
c$$$              
c$$$         enddo         
c$$$
c$$$         sigmann = s
c$$$c         write(*,*)' sigmann ',sigmann 
c$$$         
c$$$         nsf = NbNodebyFace(nf)
c$$$         XXk(:) = 0.d0 
c$$$         do is = 1,nsf
c$$$            ns = NumNodebyFace(nf,is)
c$$$            XXk(:) = XXk(:)  + 2.d0*XNode(ns,:)/dfloat(nsf)        
c$$$         enddo
c$$$c         write(*,*)' XXk ',XXk 
c$$$         
c$$$
c$$$      call EvalXPThetaBetabyFaceFrac(
c$$$     &     XXk,numfacefrac,theta,betan,betat,
c$$$     &     NbFaceFrac,NbFace,NbCell,
c$$$     &     NbNodebyFace,NumNodebyFace,
c$$$     &     NbNodebyCell,NumNodebyCell,       
c$$$     &     NbIncGlobalbyCell,
c$$$     &     NumIncGlobalbyCell,
c$$$     &     NumNodeFaceLocalbyCell,
c$$$     &     NumFaceFracVersFace,
c$$$     &     NumFacebyCell,NumCellbyFace,
c$$$     &     IndBulleNumCellbyFace,            
c$$$     &     PoidsXFaceCG,XFaceCG,
c$$$     &     GradTangbyFace, 
c$$$     &     NbIncGlobalbyFaceFrac,
c$$$     &     NumIncGlobalbyFaceFrac,
c$$$     &     Sigmanbyff,
c$$$     &     VecNormalbyFace,
c$$$     &     VecTanFace1,
c$$$     &     VecTanFace2,         
c$$$     &     EvalXPThetaBetabyff,
c$$$     &     NbIncGlobalPThetaBeta,
c$$$     &     NumIncGlobalPThetaBeta) 
c$$$           
c$$$           
c$$$         call f_SolU(XXk,XCell(k1,:),Vec)
c$$$         call f_SolU(XXk,XCell(k2,:),Vec2)
c$$$
c$$$         sautex(:) = Vec(:) - Vec2(:)
c$$$c         write(*,*)' saut ex ',sautex(:)
c$$$
c$$$         sautn = prodscal(sautex(:),VecNormalbyFace(nf,:))
c$$$
c$$$         s = 0.d0 
c$$$         do ik=1,NbIncGlobalPThetaBeta
c$$$            inc = NumIncGlobalPThetaBeta(ik)
c$$$            s = s + prodscal(EvalXPThetaBetabyff(ik,:,1),SolUEx(inc,:))
c$$$         enddo
c$$$
c$$$         sc = theta*sigmann - beta*sautn
c$$$
c$$$         write(*,*)' PThetaBeta ',s,sc,s-sc
c$$$c         write(*,*)
c$$$      enddo
c$$$      stop 
c$$$
c$$$      

c
c     tests quadrature face
c
c$$$      sint = 0.d0 
c$$$      do nff = 1,NbFaceFrac 
c$$$         numf = NumFaceFracVersFace(nff)
c$$$         call QuadFace(
c$$$     &        numf,
c$$$     &        NbNode,NbFace,
c$$$     &        XNode,
c$$$     &        NbNodebyFace,NumNodebyFace,
c$$$     &        XFaceCG,SurfaceFace,
c$$$     &        NbPoints,XPoints,PoidsPoints)       
c$$$         
c$$$
c$$$         s = 0.d0 
c$$$         do iq = 1,NbPoints
c$$$            
c$$$            fq = XPoints(iq,2)*XPoints(iq,3)
c$$$     &           + XPoints(iq,2)**2 + XPoints(iq,3)**2
c$$$     &           - 2.d0 + XPoints(iq,2) + XPoints(iq,3)
c$$$     &           + XPoints(iq,1)        
c$$$            
c$$$            s = s + PoidsPoints(iq)*fq  
c$$$            
c$$$         enddo
c$$$         write(*,*)' s ',s
c$$$
c$$$         sint = sint + s 
c$$$         
c$$$      enddo
c$$$
c$$$      sexact = ( (1.d0)**3/3.d0 - (-1.d0)**3/3.d0 )*4.d0 -2.d0*4.d0 
c$$$      write(*,*)' integral ',sint,sexact
c$$$      
c$$$      stop
c    
c       
c
cccccccccccccccccccccccccccccccccccccccc
c
c     nb de dof en rajoutant les multiplicateurs 
c

       
      NbdofMecaContact = NbdofMeca + NbFaceFrac

      
      write(*,*)'NbdofMeca ',NbdofMeca
      write(*,*)'Nbfacefrac ',NbFaceFrac     
      write(*,*)'NbdofMecaContact ',NbdofMecaContact

ccccccccccccccccccccccccccccccccccccccc
c
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Resolution d'un modele Poro-elastique fracture avec contact de type Coulomb 
c
c     Formulation mixte 
c      
c     U: VEM Bulle (aux faces frac) 
c
c     Multiplicateur: constant par face frac 
c
c     Darcy = HFV, modele a pressions discontinues 
c      
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
c
      open(unit=101,file = 'dup_fixedpoint.val',status = 'unknown')
      open(unit=102,file = 'phimoydfmoy.val',status = 'unknown')
      open(unit=103,file = 'pmmoy_pfmoy.val',status = 'unknown')

      open(unit=104,file = 'dfpf.val',status = 'unknown')

      open(unit=105,file = 'dt.val',status = 'unknown')


      open(unit=106,file = 'SautntNitsche.val',status = 'unknown')


      open(unit=107,file = 'totalnewtonmeca.val',status = 'unknown')
      open(unit=108,file = 'dUtan.val',status = 'unknown')

      open(unit=109,file = 'dunodes.val',status = 'unknown')

      
      open(unit=119,file = 'tempsaffich.val',status = 'unknown')
      
c
c
c
c
c      
c
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Coefficient de frottement de Coulomb 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

      allocate(rF(NbFaceFrac))

      
      rF(:) = f_Frottement()
      write(*,*)' rF ',rF(1)

     
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c      
c
c     
cccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     CL Dirichlet Darcy 
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
      
      allocate(SolPDir(NbCV))

      do i=1,NbCV
         
         X(:) = XInc(i,:)

        
         SolPDir(i) = f_RHOREF()*f_gravite()
     &        *( -X(3) + 2000.d0 + f_htop() ) 

c         write(*,*)' SolPDir MPa ',i,X(3),SolPDir(i)/1.d+6
      enddo
c     stop


cccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     CL Dirichlet Fourier ( to review !!!!!!!!!!!!!)
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
      
      allocate(SolTDir(NbCV))

      do i=1,NbCV
         
         X(:) = XInc(i,:)
        
c         SolTDir(i) = 345.d0 

         SolTDir(i) = 285.d0 + 4.d-2*( 2000.d0 - X(3) + f_htop() )

c         write(*,*)' T ',i,X(3),SolTDir(i)
      enddo
c      stop

      
c     
c 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     aperture au contact 
c     
      allocate(dfcontact(NbFaceFrac))

      dfcontact(:) = dfmin 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      allocate(SolP(NbCV))        
      allocate(SolP_nm1(NbCV))
      allocate(SolP_nm2(NbCV))      
      allocate(SolP_prev(NbCV))

      allocate(SolP0(NbCV))


      allocate(SolT(NbCV))        
      allocate(SolT_nm1(NbCV))
      allocate(SolT_nm2(NbCV))      
      allocate(SolT_prev(NbCV))

      allocate(SolT0(NbCV))         

      allocate(PorobyCell(NbCell))
      allocate(PorobyCell_nm1(NbCell))
      allocate(DerPorobyCell(NbCell,2))

      allocate(dfbyFaceFrac_nm1(NbFaceFrac))

      allocate(AccP(NbCV))
      allocate(AccP_nm1(NbCV))


      allocate(AccT(NbCV))
      allocate(AccT_nm1(NbCV))


cccccccccccccccccccccccccccccccccccccccccccccccccc      
c
c     Init Pression et temperature
c

      
      SolP_nm1(:) = SolPDir(:) 
      SolP = SolP_nm1
      SolP_nm2 = SolP_nm1

      SolP0 = SolP              ! on memorise la pression init


            
      SolT_nm1(:) = SolTDir(:) 
      SolT = SolT_nm1
      SolT_nm2 = SolT_nm1

      SolT0 = SolT ! on memorise la Temperature init 
c      
c
c
cccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c      
ccccccccccccccccccccccccccccccccccccccccc      
c
c     Calcul de U init : solve meca a P donne par Pinit 
c
ccccccccccccccccccccccccccccccccccccccccc
      
      allocate(SolU(NbdofMecaContact,NbDim))
      allocate(SolU0(NbdofMecaContact,NbDim))      
      allocate(SolU_nm1(NbdofMecaContact,NbDim))
      allocate(SolU_nm2(NbdofMecaContact,NbDim))      
      allocate(SolU_prev(NbdofMecaContact,NbDim))


      
      allocate(IndContact(NbFaceFrac))
c
c     init SolU_nm1 pour calcul de SolU init 
c
      SolU_nm1 = 0.d0      
c
c     init SolU pour Newton 
c     
      SolU = 0.d0               ! init a zero de SolU ici
      do i=1,NbFaceFrac
         inc = i + NbdofMeca
c         SolU(inc,1) = 100.d0 ! on force le contact a l'etat init  
      enddo      



      call CPU_TIME(time_VEM_solve1)
   
      
      Temps = 0.d0
c
c
c      
      if (iMecaVBulle.eq.0) then


      call ElasticiteNitscheFrac(
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
     &     SolU,SolU_nm1,IndContact,nitnewton)      
         
      else 

      
      call ElasticiteVEMBulleFrac(Temps,
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
     &     SolU,SolU_nm1,IndContact,nitnewton,
     &     SolT,
     &     SolP,NbCV,NumIncFaceFrac,NumIncCell,
     &     GradCellVEM)

      
      endif 
c
c     init SolU_nm1 et SolU_nm2 (pour schema fixed stress sequentiel) 
c
      SolU0 = SolU   ! on memorise le deplacement initial pour la shear dilation    
      SolU_nm1 = SolU
      SolU_nm2 = SolU_nm1 

      call CPU_TIME(time_VEM_solve2)



      CPU_VEM_solve = time_VEM_solve2 - time_VEM_solve1

      write(*,*)
      write(*,*)' CPU MECA CONTACT SOLVE ',CPU_VEM_solve 
      write(*,*)      

      do k = 1,NbCell

         
         do i=1,NbNodebyCell(k)
            inc = NumIncGlobalbyCell(k,i)
            
            n = NumNodebyCell(k,i) ! numero global du node

            if (IndDirNodeMeca(n).eq.4) then
c               write(*,*)' Uy node bord 4 ',n,XNode(n,:),SolU(inc,2)                      
            endif            
                

            
         enddo                                  
      enddo


      
      
ccccccccccccccccccccccccccccccccc
c
c     Calcul de l'epaisseur init df = dfcontact -[u]_n 
c
ccccccccccccccccccccccccccccccccc      

      do ifrac = 1,NbFaceFrac
         nf = NumFaceFracVersFace(ifrac)

         sn = 0.d0 
         do ik=1,NbIncGlobalbyFaceFrac(ifrac)
            inc = NumIncGlobalbyFaceFrac(ifrac,ik)

            do i=1,NbDim
               sn = sn + SolU(inc,i)*SautbyFaceFrac(ifrac,ik)
     &                               *VecNormalbyFace(nf,i)       
            enddo

         enddo         
c
c     projection pour Nitsche 
c
c         if (IndContact(ifrac).ne.0) then
c            sn = 0.d0 
c         endif

         
         dfcontact(ifrac) = dfcontact(ifrac) + sn
         
         dfbyFaceFrac_nm1(ifrac) = dfcontact(ifrac) - sn



c            dfbyFaceFrac_nm1(ifrac) = dfcontact(ifrac)


            write(*,*)' saut n et ind contact init ',
     &           ifrac,sn,dfbyFaceFrac_nm1(ifrac),
     &           SolU(NbdofMeca+ifrac,:),IndContact(ifrac)
         
      enddo
c      stop
ccccccccccccccccccccccccccccccccc
c
c     sautn de U et sautt de U-U0 avec projection pour Nitsche
c      
      allocate(Sautt_Proj_nm1(NbFaceFrac,NbDim-1))
      allocate(Sautt_Proj(NbFaceFrac,NbDim-1))      
      allocate(Sautn_Proj(NbFaceFrac))

      Sautt_Proj = 0.d0
      Sautt_Proj_nm1 = 0.d0      
      Sautn_Proj = 0.d0      
c
c        
ccccccccccccccccccccccccccccccccccc      
c      
c     Porosite init 
c
cccccccccccccccccccccccccccccccccc      
c

      do k=1,NbCell

         irtk = IndRockTypebyCell(k)

         if (irtk.eq.1) then              
            PorobyCell_nm1(k) = Phiminit1
         else if (irtk.eq.2) then
            PorobyCell_nm1(k) = Phiminit2                  
         else if (irtk.eq.3) then
            PorobyCell_nm1(k) = Phiminit3                           
         else
            write(*,*)' rock ne 1 2 ou 3 ',k,irtk
            stop
         endif

         
         enddo
c
c
c     accumulations initiales aux cells et face frac
c
c       
c     Chaleur volumique drainee du squelette  J.m^-3 K^-1 
c       
         C_Roche = f_C0()         
         
      do k=1,NbCell
         inck = NumIncCell(k)

         pk = SolP(inck)
         Tk = SolT(inck)         
         rhof = f_rho(pk,Tk)
         eik = f_ef(pk,Tk)
         
         AccP_nm1(inck) = VolCell(k)*PorobyCell_nm1(k)*rhof

         AccT_nm1(inck) = AccP_nm1(inck) * eik 
     &            + VolCell(k) * C_Roche * Tk    
         
      enddo

      do k=1,NbFaceFrac
         inck = NumIncFaceFrac(k)

         pk = SolP(inck)
         Tk = SolT(inck)         
         rhof = f_rho(pk,Tk)
         eik = f_ef(pk,Tk)
         
         
         nf = NumFaceFracVersFace(k)
         AccP_nm1(inck) = SurfaceFace(nf)*dfbyFaceFrac_nm1(k)*rhof

         AccT_nm1(inck) = AccP_nm1(inck) * eik 

         
c         write(*,*)' acc ',k,AccP_nm1(inck)
      enddo      
c
c
c      stop
c
ccccccccccccccccccccccccccccccccccccccccccc
c
c     init traceur 
c
      allocate(SolC(NbCell+NbFaceFrac))
      allocate(SolC_nm1(NbCell+NbFaceFrac))
      

      SolC_nm1 = 0.d0
      SolC = SolC_nm1 
      
ccccccccccccccccccccccccccccccccccccccccccccc
c
c      
c     Visu de la solution init 
c
c      
ccccccccccccccccccccccccccccccccccccccccccccc
c
      allocate(CCell(NbCell))
      allocate(PCell(NbCell))
      allocate(PFaceFrac(NbFaceFrac))
      allocate(TCell(NbCell))
      allocate(TFaceFrac(NbFaceFrac))
      allocate(CFaceFrac(NbFaceFrac))      
      allocate(UCell(NbCell,NbDim))
      allocate(SolSaut(NbFaceFrac,NbDim+1))
      allocate(SolLambda(NbFaceFrac,NbDim))
      allocate(solSigma(NbFaceFrac,NbDim))
c
c     P et T mat by Cell et P et T frac by FaceFrac 
c      
      do k=1,NbCell
         inck = NumIncCell(k)
c         PCell(k) = SolP(inck) - SolPDir(inck)
         PCell(k) = SolP(inck) - SolP0(inck)

         TCell(k) = SolT(inck) - SolT0(inck)
         
          CCell(k) = 0.d0

c          CCell(k) = dfloat(IndRockTypebyCell(k))

c          CCell(k) = perm(k,1,1)*1.d+12*Viscof ! mettre la visu en log ! 

         
      enddo

      do k=1,NbFaceFrac
         inck = NumIncFaceFrac(k)
         PFaceFrac(k) = SolP(inck) - SolP0(inck)

         TFaceFrac(k) = SolT(inck)  - SolT0(inck)
c        PFaceFrac(k) = SolP(inck)- SolPDir(inck)
         CFaceFrac(k) = 0.d0          
      enddo

c
c
c     U by Cell, du Saut by face frac, et de Lambda by face frac  
c      
      do k=1,NbCell


         
         Vec_s = 0.d0
         do ik=1,NbNodebyCell(k)
            inc = NumIncGlobalbyCell(k,ik)
            Vec_s = Vec_s + PoidsXCellCG(k,ik)*SolU(inc,:)
         enddo

          UCell(k,:) = Vec_s

      enddo

      do i=1,NbFaceFrac

           nf = NumFaceFracVersFace(i)
           k1 = NumCellbyFace(nf,1)
           k2 = NumCellbyFace(nf,2)
           

           
           Vec_s = 0.d0
           do ik=1,NbIncGlobalbyFaceFrac(i)
              inc = NumIncGlobalbyFaceFrac(i,ik)
              Vec_s = Vec_s + SautbyFaceFrac(i,ik)*SolU(inc,:)
           enddo

           SolSaut(i,1) = prodscal(Vec_s,VecNormalbyFace(nf,:)) ! intrinseque (pas de correction d'orientation) ( u+ - u-).n+ 
           SolSaut(i,2) = prodscal(Vec_s,VecTanFace1(nf,:))
     &          *OrientationbyFaceFrac(i) ! car le vec_s = ( u+ - u-) a une orientation locale liee a la normale VecNormalyFace 
           SolSaut(i,3) = prodscal(Vec_s,VecTanFace2(nf,:))
     &          *OrientationbyFaceFrac(i) ! car le vec_s = ( u+ - u-) a une orientation locale liee a la normale VecNormalyFace 

           
           SolSaut(i,4) = IndContact(i)
           
c
c     multiplicateur lambda : variable duale des composantes du saut -> meme correction d'orientation a effectuer 
c           
           Vec_s = SolU(NbdofMeca+i,:)
           Vec_s(2:3) = Vec_s(2:3)*OrientationbyFaceFrac(i) ! correction d'orientation des composantes tangentielle car lambda est en dualite avec les composantes du saut 
c        
         
           SolLambda(i,:) = Vec_s    



c calcul de la valeur moyenne du stresse sur chaque face fracture: (sigma(u+)n+ + sigma(u-)n+)/2, composante normale (scalaire n) et tangentielle (scalaire tau_1) !!! A jouter la composante tau2
c cette subroutine calcul aussi les stress sur les non-faces fractures (ou on a sigma(u+)n+ = sigma(u-)n+)

           StressFracLoc_n  = 0.d0
           StressFracLoc_t1 = 0.d0
           StressFracLoc_t2 = 0.d0
               
           numface = nf
               
           call StressFracMoy(NbCell,
     &          NbFace,
     &          NbdofMecaContact,
     &          NumCellbyFace,
     &          SolU,
     &          GradCellVEM,
     &          VecNormalbyFace,
     &          VecTanFace1,
     &          VecTanFace2,
     &          NbIncGlobalbyCell,
     &          NumIncGlobalbyCell,
     &          numface,
     &          StressFracLoc_n,
     &          StressFracLoc_t1,
     &          StressFracLoc_t2)
           
           
           solSigma(i,1) = StressFracLoc_n
           solSigma(i,2) = StressFracLoc_t1
           solSigma(i,3) = StressFracLoc_t2
                     
        enddo
       
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Affichage init: Visu de la solution P et U, [U], lambda, PorobyCell, dfbyFaceFrac   
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      Ndt = 1 ! Comptage nb de pas de temps avec affichage paraview 

      write(119,*)' Ndt ',Ndt,Temps,Temps/UnJour
      
      

      call VisuSolUP(
     &     iaffich,Ndt,
     &     Ntetra4,Nhexa8,Nnfaced,Nquad4,Ntria3,
     &     NbCell,NbFaceFrac,NbFace,
     &     NbNodebyCell,
     &     NbNodebyFace,NumFaceFracVersFace,
     &     UCell,SolSaut,SolLambda,PCell,PFaceFrac,
     &     TCell,TFaceFrac,
     &     CCell,CFaceFrac, 
     &     PorobyCell_nm1,dfbyFaceFrac_nm1)       
      
      call FichierEnseightCase(iaffich,Ndt)
c     
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc        
c
c
ccccccccccccccc
c      
c     affichage 1D frac 
c
      do k=1,NbFaceFrac
         inck = NumIncFaceFrac(k)
         nf = NumFaceFracVersFace(k)

         k1 = NumCellbyFace(nf,1)
         k2 = NumCellbyFace(nf,2)         
         inck1 = NumIncCell(k1)
         inck2 = NumIncCell(k2)

         incik1 = NumIncCellbyFaceFrac(k,1)         
         incik2 = NumIncCellbyFaceFrac(k,2)         
         
c
c     saut normal de U 
c         
         sn = 0.d0 
         do ik=1,NbIncGlobalbyFaceFrac(k)
            inc = NumIncGlobalbyFaceFrac(k,ik)
            sn = sn + SautbyFaceFrac(k,ik)   
     &           *prodscal(SolU(inc,:),VecNormalbyFace(nf,:))                        
         enddo         
c
c     saut tangentiel 
c
         sndu = 0.d0 
         do ik=1,NbIncGlobalbyFaceFrac(k)
            inc = NumIncGlobalbyFaceFrac(k,ik)
            Vec2(:) = SolU(inc,:)-SolU0(inc,:)            
            pscal = prodscal(Vec2,VecNormalbyFace(nf,:))
            sndu = sndu + pscal*SautbyFaceFrac(k,ik)       
         enddo       
         
         Vec2(:) = 0.d0
         do ik=1,NbIncGlobalbyFaceFrac(k)
            inc = NumIncGlobalbyFaceFrac(k,ik)
            Vec2(:) = Vec2(:) + SautbyFaceFrac(k,ik)
     &           *(SolU(inc,:)-SolU0(inc,:))     
         enddo
         Vec2(:) = Vec2(:) - sndu*VecNormalbyFace(nf,:)
         dUtau = dsqrt(prodscal(Vec2,Vec2))         

         sortiesfrac(k,1) = XFaceCG(nf,3)
         sortiesfrac(k,2) = dfbyFaceFrac(k)
c         sortiesfrac(k,3) = SolP(inck)-SolPDir(inck)
         sortiesfrac(k,3) = SolP(inck)         
         sortiesfrac(k,4) = SolC(NbCell+k)
         sortiesfrac(k,5) = IndContact(k)
         sortiesfrac(k,6) = sn 
         sortiesfrac(k,7) = dUtau
         sortiesfrac(k,8) = dabs(SolU(k + NbdofMeca,1))*rF(k) 
         sortiesfrac(k,9) = SolP(inck1)
         sortiesfrac(k,10) = dabs(SolU(k + NbdofMeca,3))      
         sortiesfrac(k,11) = -f_biot()*SolP(inck1) + SolP(inck)
         sortiesfrac(k,12) = SolT(inck) 
      enddo

      call  QSORTI(numsortiesfrac,NbFaceFrac,sortiesfrac(:,1))      


      do i=1,NbFaceFrac
         
         k =  numsortiesfrac(i)
         
         write(104,*)sortiesfrac(k,1:12)
         
         
      enddo      
      write(104,*)
c
cccccccccccc


c      stop

      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     BOUCLE EN TEMPS  
c      
      write(*,*)' boucle en temps '

      call CPU_TIME(TIME1_simu)
       
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       
      Temps = 0.d0      
      IndEndTimeLoop = 0

      NbDeltat = 0
      NbEchecPasDeTemps = 0

      NbTotalFixedPoint = 0
      NbTotalNewtonMeca = 0
      
c
c
      Deltat = Deltatinit
      Deltatnm1 = Deltat
c
c     indicateur d'echec de pas de temps 
c
      ireprisedt = 0
c
c     nb it boucle pt fixe ! si 1 -> algo sequentiel 
c
        nit_boucle_ptfixe = 10000
c       nit_boucle_ptfixe = 1
      
c
c     nb max d'it du pt fixe avant reprise du dt 
c      
      nit_max_ptfixe = 4000
c      
c     critere d'arret du pt fixe sur dUmax/Uref + dPmax/pref 
c      
         critere_arret_ptfixe = 1.d-5
c        critere_arret_ptfixe = 1.d-6
c        critere_arret_ptfixe = 1.d-7          
c
c     p et U ref pour calcul des increments relatifs dUmax/Uref et dPmax/pref
c     sert de critere d'arret du pt fixe       
c      
      pref = 1.d+6
      Uref = 1.d-2 


      allocate(relaxff(NbFaceFrac))
      relaxff = 0.d0 

      allocate(relaxmm(NbCell))
      relaxmm = 0.d0 
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
      do IterationTemps = 1,10000

c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
	   if ( IndEndTimeLoop.eq.1 ) then 
c
c       fin de la simulation Temps = TempsFinal
c
	      goto 1999
	   endif 
c
c
ccccccccccccccccccccccccc
c
c     reprise pas de temps 
c
 1234      continue
ccccccccccccccccccccccccc


	   if ( Temps+Deltat.gt.TempsFinal ) then 
	      Deltat = TempsFinal - Temps
              if (Deltat.lt.1.0d-4*DeltatInit) then 
                 goto 1999
              endif
	      IndEndTimeLoop = 1
	   endif 
	   Temps = Temps + Deltat 
           NbDeltat = NbDeltat + 1

           write(*,*)'ccccccccccccccccccccccccccccccccccccccccccc '
           write(*,*)'ccccccccccccccccccccccccccccccccccccccccccc '
           write(*,*)' '
           write(*,*)'PAS DE TEMPS ',NbDeltat,' t =',Temps/UnJour,
     &               '  dt =',Deltat/UnJour,
     &               ' echec dt ',NbEchecPasDeTemps
           write(*,*)' '
           write(*,*)'ccccccccccccccccccccccccccccccccccccccccccc '
           write(*,*)'ccccccccccccccccccccccccccccccccccccccccccc '


        
cccccccccccccccccccccccccccccccccccccccccccccccccccccc            
c
c
c           
ccccccccccccccccccccccccccccccccccc           
c
c     Mise a jour du Schema HFV avec dfbyFaceFrac_nm1 (potentiellement aussi avec PorobyCell_nm1) 
c
cccccccccccccccccccccccccccccccccc
           
c     perm matrice / Viscof 
c           
           do k=1,NbCell
              irtk = IndRockTypebyCell(k)
              perm(k,:,:) = 0.d0
              do m=1,NbDim
                 if (irtk.eq.1) then               
                    perm(k,m,m) = Permeabilitem1/Viscof
                 else if (irtk.eq.2) then
                    perm(k,m,m) = Permeabilitem2/Viscof               
                 else if (irtk.eq.3) then
                    perm(k,m,m) = Permeabilitem3/Viscof               
                 else
                    write(*,*)' rock ne 1 2 ou 3 ',k,irtk
                    stop
                 endif            
              enddo          
           enddo
c
c     Schema HFV 
c           
c     permeabilites normale et tangentielle*df fractures / Viscof
c           
c     approximations  EXPLICITES en temps (prises au temps nm1) 
c           
           do i=1,NbFaceFrac

              nf = NumFaceFracVersFace(i)
              k1 = NumCellbyFace(nf,1)
              k2 = NumCellbyFace(nf,2)

              permk1 = perm(k1,1,1)
              permk2 = perm(k2,1,1)

              rKfn = (permk1 + permk2)/2.d0 
                           
              rKf = 1.d0/12.d0/Viscof
              
              permfdf(i) = rKf               
              permfn(i)  = rKfn  ! perm moy des mailles k1 et k2  
              
c             write(*,*)' permfdf permfn ',i,permfdf(i),permfn(i)
           enddo

           call  HFV3D (CondThermique, 
     &      NbCell,NbNode,NbFace,NbArete,
     &     NbFaceFrac,NbCV, 
     &     XNode,NbNodebyCell,NumNodebyCell,
     &     NbFacebyCell,NumFacebyCell,
     &     NbAretebyFace,NumAretebyFace,
     &     NbNodebyFace,NumNodebyFace,
     &     NumCellbyFace,NumNodebyArete,
     &     NumFaceFracVersFace,NumFaceVersFaceFrac,IndFaceFrac,
     &     NumIncCell,NumIncFaceFrac,NumIncFace,
     &     IndIncArete,NumIncArete,NumIncCellbyFaceFrac,  
     &     perm,permfdf,permfn,dfbyFaceFrac,  
     &     XCell,XFace,XFaceFrac,XInc,
     &     VolCell,VolFaceFrac,
     &     NbInterfacebyCell,NumInterfacebyCell,
     &     NbInterfacebyFaceFrac,NumInterfacebyFaceFrac,
     &     TransCell,TransFaceFrac,TransCellbyFaceFrac,
     &     TransCellFourier,TransFaceFracFourier,
     &     TransCellbyFaceFracFourier,
     &     GCell,GFaceFrac)



         
cccccccccccccccccccccccccccccccccccccc
c
c     init SolU_prev et SolP_prev = celle de l'algo sequentiel en temps 
c
c         
         SolU_prev = SolU_nm1 + (SolU_nm1 - SolU_nm2)*Deltat/Deltatnm1

         SolP_prev = SolP_nm1 + (SolP_nm1 - SolP_nm2)*Deltat/Deltatnm1
         
         SolT_prev = SolT_nm1 + (SolT_nm1 - SolT_nm2)*Deltat/Deltatnm1
         
c
c         
c
cccccccccccccccccccccccccccccccccccc         
c
c     choix algo couplage iteratif ou sequentiel selon le pas de temps 
c
cccccccccccccccccccccccccccccccccccc

         
c         if (Ndt.eq.1) then
c            nit_boucle_ptfixe = 100
c         elseif (Ndt.eq.11) then
c         elseif (Ndt.eq.101) then            
c            nit_boucle_ptfixe = 100            
c         else
c            nit_boucle_ptfixe = 1  ! sequentiel             
c         endif
c
c         
ccccccccccccccccccccccccccccccccccccccc
c
c     Boucle de l'algorithme de point fixe 
c
c         
c
         
         
         nit_ptfixe = 0 

        
         do it_fixed_point = 1,nit_boucle_ptfixe
 
c          do it_fixed_point = 1,1
 


            
ccccccccccccccccccccccccccccccccccccccc

            nit_ptfixe = nit_ptfixe + 1


cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     calcul des coeff de relaxation Crm et Crf 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
ccccccccccccccccccccc
c
c     Cas couple 
c
ccccccccccccccccccccc
c
c
            
      relaxmm(:) = f_Crm() ! relaxation donnee par le fixed stress Crm dans Lois.f 

      ! warning, le calcul ds Probing est desactive !!!!!!!!!!

      
c$$$      if (iMecaVBulle.eq.0) then 
c$$$
c$$$
c$$$         
c$$$      call ProbingRelaxFracNitsche(
c$$$     &     NbNode,NbCell,NbFace,NbFaceFrac,
c$$$     &     NbdofMeca,NbdofMecaContact, 
c$$$     &     IndDirNodeMeca,IndFace,             
c$$$     &     NbNodebyCell,NumNodebyCell,
c$$$     &     NbNodebyFace,NumNodebyFace,
c$$$     &     NbFacebyCell,NumFacebyCell,
c$$$     &     NumCellbyFace,NumFaceFracVersFace,
c$$$     &     NumIncFaceFrac,NumIncCell,NumIncCellbyFaceFrac, 
c$$$     &     VecNormalbyFace,PoidsXCellCG,      
c$$$     &     XNode,XCellCG,XFaceCG,SurfaceFace,VolCell,      
c$$$     &     NbIncGlobalbyCell,NumIncGlobalbyCell,
c$$$     &     NumNodeFaceLocalbyCell,
c$$$     &     NbIncGlobalbyFaceFrac,NumIncGlobalbyFaceFrac,
c$$$     &     NbCV,
c$$$     &     SautbyFaceFrac,
c$$$     &     ACell,
c$$$     &     VecTanFace1,VecTanFace2,
c$$$     &     Sigmanbyff,
c$$$     &     PoidsXFaceCG,      
c$$$     &     IndBulleNumCellbyFace,
c$$$     &     GradTangbyFace,GradCellVEM,                       
c$$$     &     SolU,SolU_nm1,SolP,SolT,
c$$$     &     IndContact,Sautt_Proj_nm1,
c$$$     &     relaxff,relaxmm)        
c$$$
c$$$         
c$$$      else 


         call ProbingRelaxFracVB( 
     &        NbNode,NbCell,NbFace,NbFaceFrac,
     &        NbdofMeca,NbdofMecaContact,
     &        rF, 
     &        IndDirNodeMeca,IndFace,          
     &        NbNodebyCell,NumNodebyCell,
     &        NbNodebyFace,NumNodebyFace,
     &        NbFacebyCell,NumFacebyCell,
     &        NumCellbyFace,NumFaceFracVersFace,      
     &        VecNormalbyFace,PoidsXCellCG,PoidsXFaceCG,      
     &        XNode,XCellCG,XFaceCG,SurfaceFace,VolCell,      
     &        NbIncGlobalbyCell,NumIncGlobalbyCell,
     &        NumNodeFaceLocalbyCell,
     &        NbIncGlobalbyFaceFrac,NumIncGlobalbyFaceFrac,
     &        ACell,SautbyFaceFrac,GradCellVEM,
     &        VecTanFace1,VecTanFace2, 
     &        SolU,SolU_nm1,SolU0,SolP,IndContact,
     &        NbCV,NumIncFaceFrac,NumIncCell,
     &        relaxff,relaxmm)            

c$$$
c$$$      endif
         
c
c      relaxff(:) = f_Crf()
c
ccccccccccccccccccccccccc            
c
c     En decouple 
c            
ccccccccccccccccccccccccc
c
c            
c            relaxff(:) = 0.d0
c            relaxmm(:) = 0.d0 
c
c                

ccccccccccccccccccccccccccccccccccccc
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc                 
            
         
c           
c     Calcul de SolP a SolU_prev donne 
c
            
         call SolvePT(
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
     &     IndContact,Sautt_Proj_nm1,       
     &     SolP,SolP_nm1,SolP_prev,SolPDir,SolP0,
     &     SolT,SolT_nm1,SolT_prev,SolTDir,SolT0,           
     &     NbIncGlobalbyCell,NumIncGlobalbyCell,GradCellVEM,     
     &     NbIncGlobalbyFaceFrac,NumIncGlobalbyFaceFrac,SautbyFaceFrac,
     &     VecNormalbyFace,VecTanFace1,VecTanFace2)
c
c     Calcul de la norme max de l'increment SolP - SolP_prev 
c         
         s = 0.d0
         do k=1,NbCV
            s = dmax1(s,dabs(SolP(k)-SolP_prev(k)))        
         enddo
         dPMax_ptfixe = s/pref 

c         do k=1,NbFaceFrac
c            inck = NumIncFaceFrac(k)
c            write(*,*)' P ',k,SolP(inck),SolP_prev(inck),
c     &           SolP(inck)-SolP_prev(inck)            
c         enddo

         
         write(*,*)' dp max ',dPMax_ptfixe,pref

c         stop
c         
c
c     Calcul de SolU a SolP donnee 
c

      if (iMecaVBulle.eq.0) then


      call ElasticiteNitscheFrac(
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
     &     SolU,SolU_nm1,IndContact,nitnewton)      
         
      else 
         
      call ElasticiteVEMBulleFrac(Temps, 
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
     &     SolU,SolU_nm1,IndContact,nitnewton,
     &     SolT,
     &     SolP,NbCV,NumIncFaceFrac,NumIncCell,
     &     GradCellVEM)      


      endif

         
         
      NbTotalNewtonMeca = NbTotalNewtonMeca + nitnewton
      write(*,*)' cumul newton meca ',NbTotalNewtonMeca
      
ccccccccc
c
c     test reprise du pas de temps si Echec Newton    !!!!!! ireprisedt A COMPLETER Dans solveP et MECA  !!!!!!!!
c        
      if (ireprisedt.eq.1) then

         Temps = Temps - Deltat
         Deltat = Deltat/2.d0 
         
         NbEchecPasDeTemps = NbEchecPasDeTemps + 1
         IndEndTimeLoop = 0
         
         write(*,*)
         write(*,*)' REPRISE DU PAS DE TEMPS ',Deltat
         write(*,*)


         SolU = SolU_nm1
         SolP = SolP_nm1
         SolT = SolT_nm1
         
         AccP = AccP_nm1
         AccT = AccT_nm1          
         SolC = SolC_nm1 
         
         goto 1234         
         
      endif         
c         
c     Calcul de la norme max de l'increment SolU - SolU_prev 
c
      s = 0.d0 
      do k=1,NbdofMeca
         s = dmax1(s,dabs(SolU_prev(k,1)-SolU(k,1)))
         s = dmax1(s,dabs(SolU_prev(k,2)-SolU(k,2)))
         s = dmax1(s,dabs(SolU_prev(k,3)-SolU(k,3)))          
      enddo
      dUMax_ptfixe = s/Uref 

      
      write(*,*)
      write(*,*)' cccccccccccccccccccccccccccccccccccccccc'      
      write(*,*)' dU et dP max pt fixe ',nit_ptfixe,'dU ',dUMax_ptfixe,
     &                                        '   dP ',dPMax_ptfixe  
      write(*,*)' cccccccccccccccccccccccccccccccccccccccc'


      NbTotalFixedPoint = NbTotalFixedPoint + 1

      write(101,*)NbTotalFixedPoint,dUMax_ptfixe,dPMax_ptfixe
      
c
c     Critere d'arret du point fixe 
c      
      if (dUMax_ptfixe+dPMax_ptfixe.le.critere_arret_ptfixe) then  
         goto 2143
      else if (nit_ptfixe.ge.nit_max_ptfixe) then

         Temps = Temps - Deltat
         Deltat = Deltat/2.d0 
         
         NbEchecPasDeTemps = NbEchecPasDeTemps + 1
         IndEndTimeLoop = 0
         
         write(*,*)
         write(*,*)' REPRISE DU PAS DE TEMPS ',Deltat
         write(*,*)
         stop


         SolU = SolU_nm1
         SolP = SolP_nm1
         SolT = SolT_nm1
         
         AccP = AccP_nm1
         AccT = AccT_nm1             
         SolC = SolC_nm1 
         
         goto 1234        

         
      endif
c
c     maj des variables U et P du point fixe 
c      
      SolU_prev = SolU
      SolP_prev = SolP
      SolT_prev = SolT
      
c
c
c
ccccccccccccccccccccccccc



      
ccccccccccccccccccccccccc
      enddo  ! fin it du point fixe 
ccccccccccccccccccccccccc      

 2143 continue 

ccccccccccccccccccccccccc
c
c     Restart si glissement [Un-Unm1]_tau > dUtau_obj 
c
ccccccccccccccccccccccccc      
c
      dUtau_obj = 1.d-2 
c      
c     affichage 1D frac 
c
      dUtau_max = 0.d0 
      do k=1,NbFaceFrac
         inck = NumIncFaceFrac(k)
         nf = NumFaceFracVersFace(k)
c         
c     saut tangentiel 
c
         sndu = 0.d0 
         do ik=1,NbIncGlobalbyFaceFrac(k)
            inc = NumIncGlobalbyFaceFrac(k,ik)
            Vec2(:) = SolU(inc,:)-SolU_nm1(inc,:)            
            pscal = prodscal(Vec2,VecNormalbyFace(nf,:))
            sndu = sndu + pscal*SautbyFaceFrac(k,ik)       
         enddo       
         
         Vec2(:) = 0.d0
         do ik=1,NbIncGlobalbyFaceFrac(k)
            inc = NumIncGlobalbyFaceFrac(k,ik)
            Vec2(:) = Vec2(:) + SautbyFaceFrac(k,ik)
     &           *(SolU(inc,:)-SolU_nm1(inc,:))     
         enddo
         Vec2(:) = Vec2(:) - sndu*VecNormalbyFace(nf,:)
         dUtau = dsqrt(prodscal(Vec2,Vec2))         
         dUtau_max = dmax1(dUtau_max,dUtau)
      enddo
      write(*,*)' slip max ',dUtau_max 
      

c$$$      if ((dUtau_max.gt.dUtau_obj*1.2d0).and.(Deltat.gt.DeltatMin))
c$$$     &       then
c$$$
c$$$         Temps = Temps - Deltat
c$$$         
c$$$         Deltat = Deltat*dUtau_obj/dUtau_max
c$$$         Deltat = dmax1(DeltatMin,Deltat) 
c$$$         
c$$$         NbEchecPasDeTemps = NbEchecPasDeTemps + 1
c$$$         IndEndTimeLoop = 0
c$$$         
c$$$         write(*,*)
c$$$         write(*,*)' pas de temps adaptatif ',dUtau_max,dUtau_obj,
c$$$     &                      Deltat/UnJour
c$$$         write(*,*)
c$$$         
c$$$
c$$$         SolU = SolU_nm1
c$$$         SolP = SolP_nm1
c$$$         AccP = AccP_nm1
c$$$         SolC = SolC_nm1 
c$$$         
c$$$         goto 1234        
c$$$
c$$$         
c$$$      endif

cccccccccccccccccccccccc      
c
c     Calcul de PorobyCell et dfbyFaceFrac avec le dernier SolU = SolU_prev
c     et SolP = SolP_prev et mise a jour de l'accumulation AccP avec les nouveaux PorobyCell et dfbyFaceFrac
c      
c
      SolU_prev = SolU
      SolP_prev = SolP
      SolT_prev = SolT
      
      call Compute_Porosity_Aperture(
     &     NbCell,NbFaceFrac,NbFace,NbCV,
     &     NbdofMecaContact,
     &     dfcontact,
     &     PorobyCell,DerPorobyCell,PorobyCell_nm1,
     &     dfbyFaceFrac, 
     &     SolU_nm1,SolU_prev,SolU0,
     &     IndContact,Sautt_Proj_nm1,       
     &     SolP,SolP_nm1,SolP_prev,
     &     SolT,SolT_nm1,SolT_prev,   
     &     NbIncGlobalbyCell,NumIncGlobalbyCell,GradCellVEM,     
     &     NbIncGlobalbyFaceFrac,NumIncGlobalbyFaceFrac,SautbyFaceFrac,
     &     NumFaceFracVersFace,
     &     VecNormalbyFace,VecTanFace1,VecTanFace2,
     &     NumIncCell)          

       
       do k=1,NbCell
          inck = NumIncCell(k)

          pk = SolP(inck)
          Tk = SolT(inck)         
          rhof = f_rho(pk,Tk)
          eik = f_ef(pk,Tk)

          
          AccP(inck) = VolCell(k) * PorobyCell(k) * rhof

          AccT(inck) = AccP(inck) * eik 
     &         + VolCell(k) * C_Roche * Tk              
       enddo
       do k=1,NbFaceFrac
          inck = NumIncFaceFrac(k)
          nf = NumFaceFracVersFace(k)

          pk = SolP(inck)
          Tk = SolT(inck)         
          rhof = f_rho(pk,Tk)
          eik = f_ef(pk,Tk)
          
          AccP(inck) = SurfaceFace(nf) * dfbyFaceFrac(k) * rhof

          AccT(inck) = AccP(inck) * eik 
          
       enddo
       
ccccccccccccccccccccccccc     
c
c     calcul du traceur 
c
      call JacSmTraceur(
     &     Deltat,       
     &     NbCell,NbFaceFrac,NbFace,NbArete,NbCV,
     &     VolCell,VolFaceFrac,SurfaceFace,XInc,
     &     NumIncCell,NumIncFaceFrac, 
     &     IndIncDir,IndIncNeu,NumIncFace, 
     &     NbInterfacebyCell,NumInterfacebyCell,
     &     NumFacebyCell,NumCellbyFace,
     &     NbAretebyFace,NumAretebyFace,
     &     IndFaceFrac,NumFaceFracVersFace,NumFaceVersFaceFrac, 
     &     NbInterfacebyFaceFrac,NumInterfacebyFaceFrac,
     &     NumIncCellbyFaceFrac,            
     &     TransCell,TransFaceFrac,TransCellbyFaceFrac,
     &     AccP,AccP_nm1,dfbyFaceFrac,dfbyFaceFrac_nm1,       
     &     SolPDir,SolP,SolP0,SolT,
     &     SolC,SolC_nm1) 
c
ccccccccccccccccccccccccccc
c
c     Sautn de U et Sautt de  U-U0 avec projection pour Nitsche 
c
cccccccccccccccccccccccccc


      sauttmoy = 0.d0
      surff = 0.d0 
      
      
cccccccccccccccccccccccccccccc       
       do nff = 1,NbFaceFrac
cccccccccccccccccccccccccccccc
          
          inck = NumIncFaceFrac(nff)
          
          nf = NumFaceFracVersFace(nff)
c
c     saut normal de U 
c         
         sn = 0.d0 
         do ik=1,NbIncGlobalbyFaceFrac(nff)
            inc = NumIncGlobalbyFaceFrac(nff,ik)
            sn = sn + SautbyFaceFrac(nff,ik)   
     &           *prodscal(SolU(inc,:),VecNormalbyFace(nf,:))                        
         enddo         
c
c     saut tangentiel de SolU - SolU_nm1 
c         
         st1 = 0.d0
         do ik=1,NbIncGlobalbyFaceFrac(nff)
            inc = NumIncGlobalbyFaceFrac(nff,ik)
            Vec2(:) = SolU(inc,:)-SolU_nm1(inc,:)
            st1 = st1 + SautbyFaceFrac(nff,ik)   
     &           *prodscal(Vec2(:),VecTanFace1(nf,:))                    
         enddo

         st2 = 0.d0
         do ik=1,NbIncGlobalbyFaceFrac(nff)
            inc = NumIncGlobalbyFaceFrac(nff,ik)
            Vec2(:) = SolU(inc,:)-SolU_nm1(inc,:)
            st2 = st2 + SautbyFaceFrac(nff,ik)   
     &           *prodscal(Vec2(:),VecTanFace2(nf,:))                    
         enddo
c
c
         if (IndContact(nff).eq.1) then
            st1 = 0.d0
            st2 = 0.d0 
         endif

c         if (IndContact(nff).ne.0) then
c            sn = 0.d0
c         endif         
         

         Sautn_Proj(nff) = sn  
         Sautt_Proj(nff,1) = Sautt_Proj_nm1(nff,1) + st1  
         Sautt_Proj(nff,2) = Sautt_Proj_nm1(nff,2) + st2  
         
         dUTau = dsqrt( Sautt_Proj(nff,1)**2
     &        +  Sautt_Proj(nff,2)**2 )


         sauttmoy = sauttmoy + dUTau*SurfaceFace(nf)
         surff = surff + SurfaceFace(nf)
             
         sortiesfrac(nff,1) = XFaceCG(nf,3)             
         sortiesfrac(nff,2) = Sautn_Proj(nff)        
         sortiesfrac(nff,3) = dUTau      

cccccccccccccccccccccccccccccc
       enddo  ! fin boucle faces frac 
cccccccccccccccccccccccccccccc
c             
c
       sauttmoy = sauttmoy/surff
       
       
       call  QSORTI(numsortiesfrac,NbFaceFrac,sortiesfrac(:,1))      


      do i=1,NbFaceFrac
         
         k =  numsortiesfrac(i)
         
         write(106,*)sortiesfrac(k,1:3)


      enddo             
      write(106,*)
c
c
      write(108,*)Temps/UnJour,sauttmoy 

      
ccccccccccccccccccccccccccc
      
      write(101,*)
c
c     calcul de phim moyen et df moyen 
c      
      s = 0.d0
      vol = 0.d0
      do k=1,NbCell
         inck = NumIncCell(k)
         s = s + PorobyCell(k)*VolCell(k)
         vol = vol + VolCell(k)
      enddo
      phimoy = s/vol

      s = 0.d0
      surf = 0.d0
      do k=1,NbFaceFrac
         nf = NumFaceFracVersFace(k)
         inck = NumIncFaceFrac(k)
         s = s + dfbyFaceFrac(k)*SurfaceFace(nf)
         surf = surf + SurfaceFace(nf)
      enddo
      dfmoy = s/surf

      write(102,*)Temps/UnJour,phimoy,dfmoy 

c
c     calcul de pm moyen et pf moyen 
c
      s = 0.d0
      vol = 0.d0
      do k=1,NbCell
         inck = NumIncCell(k)
         s = s + SolP(inck)*VolCell(k)
         vol = vol + VolCell(k)
      enddo
      pmmoy = s/vol

      s = 0.d0
      surf = 0.d0
      do k=1,NbFaceFrac
         nf = NumFaceFracVersFace(k)
         inck = NumIncFaceFrac(k)
         s = s + SolP(inck)*SurfaceFace(nf)
         surf = surf + SurfaceFace(nf)
      enddo
      pfmoy = s/surf

      write(103,*)Temps/UnJour,pmmoy,pfmoy 
c
c
c
c
cccccccccccccccccccccccccccccccccccccc
c      
c     nb cumule d'it de Newton Meca 
c     
      write(107,*)Temps/UnJour,NbTotalNewtonMeca

      write(*,*)
      write(*,*)'nb cumule Newton Meca ',NbTotalNewtonMeca
      write(*,*) 
c      
c
c      
ccccccccccccccc
c      
c     affichage 1D frac 
c
      do k=1,NbFaceFrac
         inck = NumIncFaceFrac(k)
         nf = NumFaceFracVersFace(k)

         k1 = NumCellbyFace(nf,1)
         k2 = NumCellbyFace(nf,2)         
         inck1 = NumIncCell(k1)
         inck2 = NumIncCell(k2)


         incik1 = NumIncCellbyFaceFrac(k,1)         
         incik2 = NumIncCellbyFaceFrac(k,2)         
         
c
c     saut normal de U 
c         
         sn = 0.d0 
         do ik=1,NbIncGlobalbyFaceFrac(k)
            inc = NumIncGlobalbyFaceFrac(k,ik)
            sn = sn + SautbyFaceFrac(k,ik)   
     &           *prodscal(SolU(inc,:),VecNormalbyFace(nf,:))                        
         enddo         
c
c     saut tangentiel 
c
         sndu = 0.d0 
         do ik=1,NbIncGlobalbyFaceFrac(k)
            inc = NumIncGlobalbyFaceFrac(k,ik)
            Vec2(:) = SolU(inc,:)-SolU0(inc,:)            
            pscal = prodscal(Vec2,VecNormalbyFace(nf,:))
            sndu = sndu + pscal*SautbyFaceFrac(k,ik)       
         enddo       
         
         Vec2(:) = 0.d0
         do ik=1,NbIncGlobalbyFaceFrac(k)
            inc = NumIncGlobalbyFaceFrac(k,ik)
            Vec2(:) = Vec2(:) + SautbyFaceFrac(k,ik)
     &           *(SolU(inc,:)-SolU0(inc,:))     
         enddo
         Vec2(:) = Vec2(:) - sndu*VecNormalbyFace(nf,:)
         dUtau = dsqrt(prodscal(Vec2,Vec2))         
c
c     coeff de relaxation ff 
c         
c         relaxff(k) = dabs((dfbyFaceFrac(k)-dfbyFaceFrac_nm1(k))
c     &        /(SolP(inck)-SolP_nm1(inck)))   

         
         sortiesfrac(k,1) = XFaceCG(nf,3)
         sortiesfrac(k,2) = dfbyFaceFrac(k)
         sortiesfrac(k,3) = SolP(inck)
c         sortiesfrac(k,3) = SolP(inck)-SolPDir(inck)         
         sortiesfrac(k,4) = SolC(NbCell+k)
         sortiesfrac(k,5) = IndContact(k)
         sortiesfrac(k,6) = sn 
         sortiesfrac(k,7) = dUtau
         sortiesfrac(k,8) = dabs(SolU(k + NbdofMeca,1))*rF(k)  
         sortiesfrac(k,9) = SolP(inck1)         
         sortiesfrac(k,10) = dabs(SolU(k + NbdofMeca,3))   
c         sortiesfrac(k,11) =  relaxff(k)
         sortiesfrac(k,11) = -f_biot()*SolP(inck1) + SolP(inck)
         sortiesfrac(k,12) = SolT(inck)         

         
      enddo

      call  QSORTI(numsortiesfrac,NbFaceFrac,sortiesfrac(:,1))      


      do i=1,NbFaceFrac
         
         k =  numsortiesfrac(i)
         
         write(104,*)sortiesfrac(k,1:12)
         
         
      enddo            
      write(104,*)
c
c     saut aux nodes 
c

      nn = 0
      do i=1,NbFaceFrac
           nf = NumFaceFracVersFace(i)
           k1 = NumCellbyFace(nf,1)
           k2 = NumCellbyFace(nf,2)
           nsf  = NbNodebyFace(nf)

           do isf = 1,nsf
              ns = NumNodebyFace(nf,isf)
              incks1 = NumIncGlobalbyFaceFrac(i,isf) 
              incks2 = NumIncGlobalbyFaceFrac(i,isf+nsf) 
              
              Vec1(:) = SolU(incks1,:) - SolU(incks2,:)
             
              sautn = prodscal(Vec1(:),VecNormalbyFace(nf,:))
              sautt1 = prodscal(Vec1(:),VecTanFace1(nf,:))
              sautt2 = prodscal(Vec1(:),VecTanFace2(nf,:))
              sautt = dsqrt(sautt1**2 + sautt2**2)

c              write(109,*)XNode(ns,3),sautn,sautt

              nn = nn + 1
              sortiesfrac(nn,1) = XNode(ns,3)
              sortiesfrac(nn,2) = sautn
              sortiesfrac(nn,3) = sautt
              
              
           enddo                   
        enddo

      call  QSORTI(numsortiesfrac,nn,sortiesfrac(:,1))      


      do i=1,nn
         
         k =  numsortiesfrac(i)
         
         write(109,*)sortiesfrac(k,1:3)
         
         
      enddo            
      write(109,*)

        
c        write(109,*)
        
      
ccccccccccccccccccccccccccccccccccc
c
c     glissement max de Un - Unm1 
c      
      dUtau_max = 0.d0 
      do k=1,NbFaceFrac
         inck = NumIncFaceFrac(k)
         nf = NumFaceFracVersFace(k)
c         
c     saut tangentiel 
c
         sndu = 0.d0 
         do ik=1,NbIncGlobalbyFaceFrac(k)
            inc = NumIncGlobalbyFaceFrac(k,ik)
            Vec2(:) = SolU(inc,:)-SolU_nm1(inc,:)            
            pscal = prodscal(Vec2,VecNormalbyFace(nf,:))
            sndu = sndu + pscal*SautbyFaceFrac(k,ik)       
         enddo       
         
         Vec2(:) = 0.d0
         do ik=1,NbIncGlobalbyFaceFrac(k)
            inc = NumIncGlobalbyFaceFrac(k,ik)
            Vec2(:) = Vec2(:) + SautbyFaceFrac(k,ik)
     &           *(SolU(inc,:)-SolU_nm1(inc,:))     
         enddo
         Vec2(:) = Vec2(:) - sndu*VecNormalbyFace(nf,:)
         dUtau = dsqrt(prodscal(Vec2,Vec2))         
         dUtau_max = dmax1(dUtau_max,dUtau)
      enddo
      write(*,*)' slip max [Un - Unm1]_tau  ',dUtau_max       
ccccccccccccccccccccccccccccccccccccccc
c
c     Mise a jour des variables a nm1 et nm2 
c      
ccccccccccccccccccccccccccccccccccccccc
c

      SolP_nm2 = SolP_nm1
      
      SolP_nm1 = SolP


      SolT_nm2 = SolT_nm1
      
      SolT_nm1 = SolT

            

      SolU_nm2 = SolU_nm1
      
      SolU_nm1 = SolU


      Sautt_Proj_nm1 = Sautt_Proj 
      
cccccccccccccccccccccccccccccccccccccc      
c
c     Mise a jour des accumulations df et phim a tnm1 
c
      AccP_nm1 = AccP

      AccT_nm1 = AccT      
c
      PorobyCell_nm1 = PorobyCell

      dfbyFaceFrac_nm1 = dfbyFaceFrac 
c
cccccccccccccccccccccccccccccccccccccc
c
c     mise a jour du traceur 
c
      SolC_nm1 = SolC 
c
c      
c      
ccccccccccccccccccccccccccccccccccccccc   
c
c     maj de dtnm1 
c         
         Deltatnm1 = Deltat 
c
c
c
c              
ccccccccccccccccccccccccccccccccccccc      
c
c
c     marche en temps 
c
ccccccccccccccccccccccccc
c
c     Pas de temps adaptatif fct du glissement 
c      
c
      dUtau_obj = 3.d-3 

c      if (dUtau_max.gt.dUtau_obj) then
c         Deltat = Deltat*dUtau_obj/dUtau_max
c         Deltat = dmax1(DeltatMin,Deltat)          
c      else                  
c        alpha = 1.1d0      
c        Deltat = Deltat*alpha       
c        Deltat = dmin1(Deltat,DeltatMax)              
c      endif

      

c       if (Temps.le.36.5d0*UnJour) then                   
c         Deltat = Deltat1
c       else if (Temps.le.Unjour*100.d0) then     
c         Deltat = Deltat2
c      else
c         Deltat = Deltat3 
c      endif


c       if (Temps.le.Temps12) then                   
c         Deltat = Deltat1
c       else 
c         Deltat = Deltat2    
c      endif      

c$$$
c$$$       if (Temps.le.39.5d0*UnJour) then                   
c$$$         Deltat = Deltat1
c$$$       else if (Temps.le.Unjour*41.45d0) then     
c$$$         Deltat = Deltat2
c$$$      else if (Temps.le.Unjour*41.8d0) then     
c$$$         Deltat = Deltat3
c$$$      else
c$$$         Deltat = Deltat4         
c$$$      endif

c$$$      
c$$$       if (Temps.le.59.5d0*UnJour) then                   
c$$$         Deltat = Deltat1
c$$$       else if (Temps.le.Unjour*77.5d0) then     
c$$$         Deltat = Deltat2
c$$$      else if (Temps.le.Unjour*81.d0) then     
c$$$         Deltat = Deltat3
c$$$      else
c$$$         Deltat = Deltat4         
c$$$      endif

      

c       if (Temps.le.89.d0*UnJour) then                   
c         Deltat = Deltat1
c       else if (Temps.le.Unjour*150.d0) then     
c         Deltat = Deltat2
c      else
c         Deltat = Deltat3 
c      endif      



        
       if (Temps.le.Unjour*10.d0) then                   
          Deltat = Deltat1
       else if (Temps.le.Unjour*50.d0) then
          Deltat = Deltat2          
       else if (Temps.le.Unjour*500.d0) then
          Deltat = Deltat3          
       else 
          Deltat = Deltat4   
      endif 
      
      
      write(105,*)Temps/UnJour,Deltat/UnJour
c
ccccccccccccccccccccccccccccccccccccc

c      do ifrac = 1,NbFaceFrac

c         write(*,*)' df et ind contact ',
c     &        ifrac,dfbyFaceFrac(ifrac),IndContact(ifrac)
         
c      enddo
c
cccccccccccccccccccccccccccccccccccc
c
c     erreur solT sans frac Neumann non homogene a gauche Dirichlet a droite 
c
c      ermax = 0.d0
c      do k = 1,NbCV
         
c         xk = XInc(k,1)

c         erk = SolT(k)-SolT0(k)-10.d0*(1.d0- xk/2000)
c         ermax = dmax1(dabs(erk),ermax)
         
c      enddo
c      write(*,*)
c      write(*,*)'erreur max T ',ermax
c      write(*,*)
c
c
c      
ccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Visu de la solution P et U, [U], lambda, Poro et df   
c     
ccccccccccccccccccccccccccccccccccccccccccccc
c

        
      do k=1,NbCell
         inck = NumIncCell(k)
         PCell(k) = SolP(inck) - SolP0(inck)

         TCell(k) = SolT(inck)  - SolT0(inck)
c     PCell(k) = IndRockTypebyCell(k)

         if (dabs(SolC(k)).gt.1.d-10) then          
            CCell(k) = SolC(k)
         else
            CCell(k) = 0.d0 
         endif
c         write(*,*)' error cell ',k,PCell(k)
      enddo

      do k=1,NbFaceFrac
         inck = NumIncFaceFrac(k)
         PFaceFrac(k) = SolP(inck) - SolP0(inck)

         TFaceFrac(k) = SolT(inck)  - SolT0(inck)

         if (dabs(SolC(k+NbCell)).gt.1.d-10) then                  
            CFaceFrac(k) = SolC(k+NbCell)
         else
            CFaceFrac(k) = 0.d0             
         endif
c         write(*,*)' error ff ',k,PFaceFrac(k)
      enddo


      do k=1,NbCell
        
         Vec_s = 0.d0
         do ik=1,NbNodebyCell(k)
            inc = NumIncGlobalbyCell(k,ik)
            Vec_s = Vec_s + PoidsXCellCG(k,ik)*SolU(inc,:)
         enddo

          UCell(k,:) = Vec_s

      enddo



      
      do i=1,NbFaceFrac

           nf = NumFaceFracVersFace(i)
           k1 = NumCellbyFace(nf,1)
           k2 = NumCellbyFace(nf,2)
                    
           Vec_s = 0.d0
           do ik=1,NbIncGlobalbyFaceFrac(i)
              inc = NumIncGlobalbyFaceFrac(i,ik)
              Vec_s = Vec_s + SautbyFaceFrac(i,ik)*SolU(inc,:)
           enddo

           SolSaut(i,1) = prodscal(Vec_s,VecNormalbyFace(nf,:)) 
           SolSaut(i,2) = prodscal(Vec_s,VecTanFace1(nf,:))
     &          *OrientationbyFaceFrac(i) 
           SolSaut(i,3) = prodscal(Vec_s,VecTanFace2(nf,:))
     &          *OrientationbyFaceFrac(i) 
           
           SolSaut(i,4) = IndContact(i)
           
c
c     multiplicateur lambda 
c           
           Vec_s = SolU(NbdofMeca+i,:)
           Vec_s(2:3) = Vec_s(2:3)*OrientationbyFaceFrac(i)        
         
           SolLambda(i,:) = Vec_s    
           
           
        enddo


cccccccccccccccccc

        epst = Deltat2*1.d-2

        nnat = naffichpasdetemps

c        if ( (NbDeltat.eq.2).or.((NbDeltat/nnat)*nnat.eq.NbDeltat)
c     &       .or.(Temps.ge.TempsFinal-epst) ) then 

        
           Ndt = Ndt + 1

           write(119,*)' Ndt ',Ndt,Temps,Temps/UnJour
           
           call VisuSolUP(
     &          iaffich,Ndt,
     &          Ntetra4,Nhexa8,Nnfaced,Nquad4,Ntria3,
     &          NbCell,NbFaceFrac,NbFace,
     &          NbNodebyCell,
     &          NbNodebyFace,NumFaceFracVersFace,
     &          UCell,SolSaut,SolLambda,PCell,PFaceFrac,
     &          TCell,TFaceFrac,
     &          CCell,CFaceFrac,            
     &          PorobyCell,dfbyFaceFrac)       
           
           
           call FichierEnseightCase(iaffich,Ndt)

c        endif
c     
cccccccccccccccc
c        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      call CPU_TIME(TIME2_simu)

      cputime = TIME2_simu-TIME1_simu

       write(*,*)' CPU Time simu ',cputime         
       write(*,*)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Fin boucle en temps 
c
      enddo
c
 1999 continue 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      close(101)
      close(102)
      close(103)      
c
c
c
ccccccccccccccccccccccccccccccccc
c     
c     
      write(*,*)' MemBloc ',MemBloc
      write(*,*)' NbCell ',NbCell
      write(*,*)' NbFaceFrac ',NbFaceFrac
      write(*,*)' NbCV ',NbCV
c
c        
c
c      
      deallocate(UCell)
      deallocate(SolSaut)
      deallocate(SolLambda)
      deallocate(solSigma)

      deallocate(PCell,PFaceFrac)

      deallocate(TCell,TFaceFrac)

      deallocate(NbNodebyCell)
      deallocate(NumNodebyCell)
      
      deallocate(NbNodebyFace)
      deallocate(NumNodebyFace)
      
      deallocate(XNode)
      
      deallocate(TransCell)
      deallocate(TransFaceFrac)
      deallocate(TransCellbyFaceFrac)

      deallocate(TransCellFourier)
      deallocate(TransFaceFracFourier)
      deallocate(TransCellbyFaceFracFourier)
      
      deallocate(GCell)
      deallocate(GFaceFrac)
             
      deallocate(NbInterfacebyCell)
      deallocate(NumInterfacebyCell)
      deallocate(NbInterfacebyFaceFrac)
      deallocate(NumInterfacebyFaceFrac)

      deallocate(NumIncCell)
      deallocate(NumIncFaceFrac)
      deallocate(NumIncCellbyFaceFrac)


      deallocate(SolP)
      
      deallocate(SolT)

      deallocate(XInc)


      deallocate(NumFaceFracVersFace)

      deallocate(SolPDir)
      deallocate(permfdf,permfn)
c
      deallocate(NbCellbyNode)
      deallocate(NumCellbyNode)

      deallocate(LabelCellbyNode)
      deallocate(NbLabelbyNode)

      deallocate(IndNodeFrac)

      deallocate(NbIncGlobalbyCell)
      deallocate(NumIncGlobalbyCell)
      deallocate(NumNodeFaceLocalbyCell)       

      deallocate(NbIncGlobalbyFaceFrac)
      deallocate(NumIncGlobalbyFaceFrac)            
      deallocate(SautbyFaceFrac)



      deallocate(perm)
      deallocate(NumNodebyArete)
      deallocate(NbAretebyFace)
      deallocate(NumAretebyFace)
      deallocate(NumCellbyFace)
      deallocate(NbFacebyCell)
      deallocate(NumFacebyCell)
      

      deallocate(IndIncArete)
      deallocate(NumIncArete)
      deallocate(NumIncNode)
      deallocate(NumIncFace)

      deallocate(IndFaceFrac)
      deallocate(NumFaceVersFaceFrac)

      deallocate(XCell)
      deallocate(XFace)
      deallocate(XFaceFrac)

      deallocate(IndFace)

      deallocate(XFaceCG)
      deallocate(SurfaceFace)
      deallocate(VecNormalKSigma)

      deallocate(SizeArete)            
      deallocate(VecNormalSigmaArete)
      deallocate(PoidsXFaceCG)


      deallocate(GradCellVEM)

      deallocate(ACell) 


      deallocate (VecTanFace1)
      deallocate (VecTanFace2)
     
      deallocate(diamK)

      deallocate(IndDirFace)
      deallocate(IndDirNodeMeca)            
      
      deallocate(SolU,SolU0)
      deallocate(SolU_nm1)
      deallocate(SolU_nm2)
      deallocate(SolU_prev)
 
      deallocate(XCellCG)
      deallocate(PoidsXCellCG)


      deallocate(piK)

      deallocate(DofPiCellVEM)

      deallocate(IndBulleNumCellbyFace) 
      deallocate(LabelbyFace)
      deallocate(IndIncDir)
      deallocate(VolCell)
      deallocate(VolFaceFrac)
      deallocate(VecNormalbyFace)

      deallocate(IndContact)

      deallocate(OrientationbyFaceFrac)
c      
c
      deallocate(rF)


      deallocate(dfcontact)
      
      deallocate(SolP_nm1)
      deallocate(SolP_prev)
      deallocate(SolP0)


      deallocate(SolT_nm1)
      deallocate(SolT_prev)
      deallocate(SolT0)

      deallocate(PorobyCell)
      deallocate(PorobyCell_nm1)
      deallocate(DerPorobyCell)

      deallocate(dfbyFaceFrac)
      deallocate(dfbyFaceFrac_nm1)

      deallocate(AccP)
      deallocate(AccP_nm1)

      deallocate(AccT)
      deallocate(AccT_nm1)

      deallocate(IndRockTypebyCell)

      deallocate(IndIncNeu,IndNeuFace) 

      deallocate(SolC)
      deallocate(SolC_nm1)

      deallocate(CCell,CFaceFrac)


      deallocate(GradTangbyFace)
      deallocate(Sigmanbyff)

      
ccccccccccccccccccc      
c

      enddo ! fin boucle imesh 
c
cccccccccccccccccc      
c     

      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     
