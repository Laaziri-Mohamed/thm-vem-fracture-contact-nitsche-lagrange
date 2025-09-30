cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     Numerotation des Inconnues pour le schema HFV  
c
c     fct  des CLs 
c
c     HFV: modele discontinu 
c
c         Aretes frac DIR, Faces DIR, Mailles, faces frac , 
c         faces non DIR et non frac, Ksigma et Lsigma aux sigma frac, aretes frac non DIR 
c
c
c     Hyp: les faces frac ne sont pas des aretes de bord 
c
ccccccccccccccccccc
c
c     EN ENTREE
c
c
c     - Faces frac
c
c         * IndFaceFrac = 1 si face frac, 0 sinon 
c
c         * NumFaceVersFaceFrac: numface > numfacefrac
c     
c     - CL Dirichlet 
c
c
c          * pour HFV: IndIncArete = 11 frac DIR, 10 frac non DIR, < 10 si non frac
c
c                        IndDirface = 1 si DIR, 0 sinon 
c
cccccccccccccccccc
c
c                  
c
c     EN SORTIE
c
c     NbCV, 
c
c     IndIncDir = 1 si inconnue DIR, 0 sinon
c
c     IndIncFaceFrac -> 1 si inc frac (pour HFV = arete ou face frac)
c
c     Pour HFV: NumIncArete, NumIncFace, NumIncCell, NumIncFaceFrac, NumIncCellbyFaceFrac 
c
c                 NumIncArete donne le num de l'inc arete si arete frac 
c                 (ie si IndIncArete >=10)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
      subroutine NumerotationInc(
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
c       num face  > num face frac
        dimension NumFaceVersFaceFrac(NbFace)

c       num facefrac  > num face 
        dimension NumFaceFracVersFace(NbFaceFrac)
c
c
c     Faces par les noeuds 
c       
        dimension NbNodebyFace(NbFace)
        dimension NumNodebyFace(NbFace,NbNodebyFaceMax)
c
c     Faces par les aretes 
c
        dimension NbAretebyFace(NbFace)
        dimension NumAretebyFace(NbFace,NbAretebyFaceMax)
c
c       num face > ind 1 si face frac 0 sinon 
        dimension IndFaceFrac(NbFace)

c
c     Choix de la CL par face pour P,T (flow) sur les faces 
c        
        dimension IndDirFace(NbFace)
        dimension IndNeuFace(NbFace)
        dimension IndDirFaceT(NbFace)
        dimension IndNeuFaceT(NbFace)              
c
c     pour HFV (toutes les aretes fracture)
c
c     Ind 11 arete frac DIR
c     Ind 10 arete frac non DIR 
c     Ind 0 non frac non DIR
c     Ind 1 non frac DIR 
c
        dimension IndIncArete(NbArete)
        dimension IndIncAreteT(NbArete)        
c
cccccccccccccccccccccccccccccccccccccccc
c
c     SORTIES 
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
c     Numero des inconnues interfaces 
c
c
c     pour HFV (toutes les faces y compris les faces frac)
c     num face  > num inc 
c
        dimension NumIncFace(NbFace)
c
c        
c     pour HFV discontinu 
c     num face frac  > num inc Ksigma et Lsigma 
c        
        dimension NumIncCellbyFaceFrac(NbFaceFrac,2)        
c
c
c     num arete  > num inc 
        dimension NumIncArete(NbArete)
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
c
        dimension IndIncFaceFrac(NbCV)
c
c
c     num inc Neumann > indice 0 ou 1 ou 2 
c
        dimension IndIncNeu(NbCV)
        dimension IndIncNeuT(NbCV)
        
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      write(*,*)' numerotation flow '
        
c
cccccccccccccccccccccccccccccccccccccc
c
c     Numerotation des inconues pour HFV 
c
cccccccccccccccccccccccccccccccccccccc
c     
c      Numerotation des inconnues HFV:
c      Arete Frac Dir, FaceDir, Cell, Face Frac,
c      autres faces non frac non dir, interfaces Ksigma,Lsigma, aretes frac non dir
c

         
      IndIncFaceFrac = 0    

      IndIncDir = 0
      IndIncNeu = 0

     
      IndIncDirT = 0
      IndIncNeuT = 0       
     
      NumIncArete = 0         


cccccccccccccccccccc

      NbCV1 = 0
      
      do i=1,NbCell ! cells 
         NbCV1 = NbCV1 + 1
         NumIncCell(i) = NbCV1  
      enddo

      write(*,*)' fin cell ',NbCV1 


      
      do i=1,NbFace ! faces 

         NbCV1 = NbCV1 + 1
         NumIncFace(i) = NbCV1

         if (IndFaceFrac(i).eq.1) then             
            j = NumFaceVersFaceFrac(i)
            NumIncFaceFrac(j) = NbCV1 
            IndIncFaceFrac(NbCV1) = 1 ! inc face frac 
         endif

         
         if (IndDirFace(i).eq.1) then
            IndIncDir(NbCV1) = 1  ! inc dir pression           
         endif
         
         if (IndDirFaceT(i).eq.1) then
            IndIncDirT(NbCV1) = 1  ! inc dir temp           
         endif 
         
         


         if ( (IndDirFace(i).eq.0).and.(IndFaceFrac(i).eq.0) ) then         
            if (IndNeuFace(i).eq.1) then 
               IndIncNeu(NbCV1) = 1 ! inc neumann homogene pression 
            else if (IndNeuFace(i).eq.2) then
               IndIncNeu(NbCV1) = 2  ! inc neumann non homogene pression              
            endif
         endif


         if ( (IndDirFaceT(i).eq.0).and.(IndFaceFrac(i).eq.0) ) then         
            if (IndNeuFaceT(i).eq.1) then 
               IndIncNeuT(NbCV1) = 1  ! inc neumann homogene fourier temp          
            endif
         endif         

            
      enddo

      write(*,*)' fin face  ',NbCV1 

      
      do i=1,NbFaceFrac
         
         NbCV1 = NbCV1 + 1        ! inc Ksigma
         NumIncCellbyFaceFrac(i,1) = NbCV1

         NbCV1 = NbCV1 + 1        ! inc Lsigma
         NumIncCellbyFaceFrac(i,2) = NbCV1

         
      enddo

      write(*,*)' fin mf ',NbCV1 

c      
c     on ne distingue par les arete frac de bord neumann (neumann homogene par defaut si non dir) 
c     a completer plus tard !!
c
c      
      do i=1,NbArete
         
         if (IndIncArete(i).ge.10) then
            
            NbCV1 = NbCV1 + 1
            NumIncArete(i) = NbCV1 ! arete frac bord dir
            IndIncFaceFrac(NbCV1) = 1

            if (IndIncArete(i).eq.11) then ! arete frac dir pression             
               IndIncDir(NbCV1) = 1 ! inc dir pression 
            endif

            if (IndIncAreteT(i).eq.11) then  ! arete frac dir temperature            
               IndIncDirT(NbCV1) = 1 ! inc dir temp 
            endif
            
         endif
                
      enddo   

      write(*,*)' fin aretes ',NbCV1 
c
c      
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       return
       end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
