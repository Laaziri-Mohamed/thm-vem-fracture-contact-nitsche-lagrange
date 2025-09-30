ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
      subroutine Calcul_LabelCellbyNode(NbNode,NbCell,NbFace,
     &     NbCellbyNode,NumCellbyNode,
     &     NbFacebyCell,NumFacebyCell,NumCellbyFace,
     &     IndFaceFrac,IndNodeFrac,      
     &     NbLabelbyNode,LabelCellbyNode) 
c
c
	implicit double precision (a-h,o-z)
	implicit integer (i-n)
c
c
        include 'include_parameter'
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     DECLARATIONS
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ENTREES 
c
c     NbNode
c
c     mailles voisines de chaque node 
c
        dimension NbCellbyNode(NbNode)
        dimension NumCellbyNode(NbNode,NbCellbyNodeMax)

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
c       num face > ind 1 si face frac 0 sinon 
        dimension IndFaceFrac(NbFace)
c
c     indnodefrac = 1 si nodefrac, 0 sinon 
c        
        dimension IndNodeFrac(NbNode)        
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     SORTIES
c
c
c     Label des composantes connexes sur les ensembles de mailles par node coupees par les faces frac 
c        
        dimension NbLabelbyNode(NbNode)
        dimension LabelCellbyNode(NbNode,NbCellbyNodeMax)             
c
c
c        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        
c     WorkSpaces 
c
c
c        dimension IndLabel(NbCellbyNodeMax)
c        dimension NumLabel(NbCellbyNodeMax)
c
      integer, dimension(:), allocatable :: IndLabel,NumLabel

        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      allocate(IndLabel(NbCellbyNodeMax))
      allocate(NumLabel(NbCellbyNodeMax))
      

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     classes d'equivalences pour les ensembles cellbynode 
c
c     deux mailles sont dans la meme classe si elles partagent une face non frac 
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call CPU_TIME(TIME1_classenode)
      
c
cccccccccccccccccccccc      
      do n=1,NbNode


c
c     label a 1 pour toutes les mailles si non frac
c

         ncell = NbCellbyNode(n)

         
         do nk = 1,ncell
            LabelCellbyNode(n,nk) = 1          
         enddo

         NbLabelbyNode(n) = 1          
         
      if (IndNodeFrac(n).eq.1) then 
         
ccccccccccccccccccccc
       
 

         do nk = 1,ncell
            LabelCellbyNode(n,nk) = nk          
         enddo
      
         lblm1 = 0
         do nk = 1,ncell

            lbl = LabelCellbyNode(n,nk)
            if (lbl.ne.lblm1) then 
c
c     propagation du label lbl de nk 
c
            do mk = 1,ncell
            do nk1 = 1,ncell
               k1 = NumCellbyNode(n,nk1)


               if (LabelCellbyNode(n,nk1).eq.lbl) then 
            
               do nk2 = 1,ncell
                  k2 = NumCellbyNode(n,nk2)

                  
c            
c     on teste si k1 et k2 appartiennent a la meme classe (ie ils partage une face non frac) 
c               
                  do iface2 = 1,NbFacebyCell(k2)
                     nface2 = NumFacebyCell(k2,iface2)
                     
                     if (IndFaceFrac(nface2).eq.0) then
                        
                        k21 = NumCellbyFace(nface2,1)
                        k22 = NumCellbyFace(nface2,2)
                        
                        if ( (k21.eq.k1).or.(k22.eq.k1) ) then
                           LabelCellbyNode(n,nk2) = lbl
c                           write(*,*)' on propage le label lbl ',nk2,lbl
                        endif
                        
                     endif
                  enddo              
            
               enddo
               endif

            enddo
            enddo

            lblm1 = lbl

            
         endif 
         enddo



         

         IndLabel(:) = 0
         do l=1,ncell
            IndLabel(LabelCellbyNode(n,l)) = 1
         enddo
         nblabel = 0
         do l=1,ncell
            if (IndLabel(l).eq.1) then
               nblabel = nblabel + 1
               do nk=1,ncell
                  if (LabelCellbyNode(n,nk).eq.l) then
                     NumLabel(nk) = nblabel 
                  endif
               enddo
            endif                        
         enddo
         NbLabelbyNode(n) = nblabel

         if (nblabel.gt.NbLabelbyNodeMax) then
            write(*,*)' nblabel >  NbLabelbyNodeMax ',
     &           nblabel,NbLabelbyNodeMax
            stop
         endif
         
         do l=1,ncell
            LabelCellbyNode(n,l) = NumLabel(l)
         enddo


c         if (IndNodeFrac(n).eq.1) then
c            write(*,*)n,NbLabelbyNode(n)
c            do l=1,ncell
c               write(*,*)l,NumCellbyNode(n,l),LabelCellbyNode(n,l)
c            enddo
c            write(*,*)
c         endif
         
ccccccccccccccccc
      endif
      enddo
ccccccccccccccccc
c      
c
      call CPU_TIME(TIME2_classenode)

      write(*,*)' CPU classe node ',TIME2_classenode-TIME1_classenode
c
c      stop

cccccccccccccc

      deallocate(IndLabel)
      deallocate(NumLabel)
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       return
       end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
      subroutine NumIncMeca(NbNode,NbCell,NbFace,NbFaceFrac, 
     &     NbCellbyNode,NumCellbyNode,
     &     NbNodebyCell,NumNodebyCell,
     &     NbFacebyCell,NumFacebyCell,NumCellbyFace,
     &     IndFaceFrac,IndBulleNumCellbyFace,
     &     NbLabelbyNode,LabelCellbyNode,
     &     NumIncGlobalbyCell,NbIncGlobalbyCell,
     &     NumNodeFaceLocalbyCell,
     &     NbdofMeca)    
c
c
	implicit double precision (a-h,o-z)
	implicit integer (i-n)
c
c
        include 'include_parameter'
c
c
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
c     mailles voisines de chaque node 
c
        dimension NbCellbyNode(NbNode)
        dimension NumCellbyNode(NbNode,NbCellbyNodeMax)
c
c     nb de noeuds S par maille K 
c
        dimension NbNodebyCell(NbCell)
c
c     num des noeuds S de chaque maille K dans l'ordre cyclique 
c
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
c       num face > ind 1 si face frac 0 sinon 
        dimension IndFaceFrac(NbFace)   
c
c     Label des composantes connexes sur les ensembles de mailles par node coupees par les faces frac 
c        
        dimension NbLabelbyNode(NbNode)
        dimension LabelCellbyNode(NbNode,NbCellbyNodeMax)         

c
c     Indice bulle stocke dans la structure NumCellbyFace 
c
        dimension IndBulleNumCellbyFace(NbFace,2)      
        
c
c        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     SORTIES
c
c       NbdofMeca
c
c     Numero global d'inconnue mecanique CellbyNode
c
c     dimension NumIncGlobalCellbyNode(NbNode,NbCellbyNodeMax)

      integer, dimension(:,:), allocatable :: NumIncGlobalCellbyNode          
c
c     Numero global d'inconnue mecanique CellbyFace 
c
c      dimension NumIncGlobalCellbyFace(NbFace,2)

      integer, dimension(:,:), allocatable :: NumIncGlobalCellbyFace        
c
c     Numerotation locale par cell k ds l'ordre
c               NumNodebyCell(k,i), i=1,...,NbNodebyCell(k) :  inc Ks
c               puis les inc face frac Ksigma   si  k = NumCellbyFace(n,1) (cote maille 1 de la face) 
c     Local -> global 
c        
        dimension NbIncGlobalbyCell(NbCell)
        dimension NumIncGlobalbyCell(NbCell,NbIncbyCellMax)        
c
c     Numero local du noeud ou de la face par inc cell 
c
        dimension NumNodeFaceLocalbyCell(NbCell,NbIncbyCellMax)                
c
c
c        
cccccccccccccccccccccccccccccccccccccccccccccccccc      
c
c
c     Numerotation globale et locale pour la mecanique: VEM ordre 1 sur Omega/Gamma 
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccc

        allocate(NumIncGlobalCellbyNode(NbNode,NbCellbyNodeMax))
        allocate(NumIncGlobalCellbyFace(NbFace,2))
c
c     numerotation globale des \overline Ks
c      
      numinc = 0
      do n=1,NbNode
         ncell = NbCellbyNode(n)       
         
         do nlabel = 1,NbLabelbyNode(n)
            
            numinc = numinc + 1

            
            do nk=1,ncell
               nl = LabelCellbyNode(n,nk)
               if (nl.eq.nlabel) then 
                  NumIncGlobalCellbyNode(n,nk) = numinc 
               endif
            enddo
            
         enddo
         
      enddo

      NumIncKs = numinc
      write(*,*)
      write(*,*)' nb inc meca Ks ',NumincKs 
      

c      do n=1,NbNode
c         ncell = NbCellbyNode(n)       
c         write(*,*)' node ',n,NumIncGlobalCellbyNode(n,1:ncell)
c      enddo
c      stop
c
c     numerotation globale des faces \sigma aux faces non frac et K\sigma aux faces frac  
c
      do n=1,NbFace
         if (NumCellbyFace(n,2).le.0) then
            ncell = 1
         else
            ncell = 2
         endif

         if (IndFaceFrac(n).eq.1) then
            if (ncell.ne.2) then
               write(*,*)' face frac sur le bord ',n
               stop
            endif

            if (IndBulleNumCellbyFace(n,1).eq.1) then 
               numinc = numinc + 1
               NumIncGlobalCellbyFace(n,1) = numinc
            endif
            
            if (IndBulleNumCellbyFace(n,2).eq.1) then             
               numinc = numinc + 1
               NumIncGlobalCellbyFace(n,2) = numinc
            endif
            
         else

            if (ncell.eq.2) then
               if ( IndBulleNumCellbyFace(n,1).ne.
     &              IndBulleNumCellbyFace(n,2) ) then
                  write(*,*)' les indices bulle doit etre egaux
     &            sur une face int non frac ',n
                  stop
               endif
            endif
            
            if (IndBulleNumCellbyFace(n,1).eq.1) then             
               numinc = numinc + 1            
               do nc=1,ncell
                  NumIncGlobalCellbyFace(n,nc) = numinc
               enddo
            endif
            
         endif

         
      enddo

c      do n=1,NbFace
c         if (NumCellbyFace(n,2).le.0) then
c            ncell = 1
c         else
c            ncell = 2
c         endif
c         write(*,*)' inc face ',NumIncGlobalCellbyFace(n,1:ncell)         
c      enddo

      
      NbdofMeca = numinc

      write(*,*)' nb inc meca Ksig ',NbdofMeca-NumIncKs        
      write(*,*)' nb inc meca ',NbdofMeca
     
      write(*,*)

c      stop
      
c      if ((NbdofMeca-NumincKs).ne.NbFaceFrac) then
c         write(*,*)' NbdofMeca ne nb Ks + nb face frac '
c         stop
c      endif

c
c     
c     Correspondance num locale -> globale par maille 
c
c      
      do k=1,NbCell
         

c     Inc Ks pour les noeuds de k
         
         do nk=1,NbNodebyCell(k)
            n = NumNodebyCell(k,nk)
            ncell = NbCellbyNode(n)
            nkn = -1 
            do nc = 1,ncell
               kp = NumCellbyNode(n,nc)               
               if (kp.eq.k) then
                  nkn = nc 
               endif
            enddo

            if (nkn.eq.-1) then
               write(*,*)' cell autour du node non trouvee ',k,nk
               stop
            else
               NumIncGlobalbyCell(k,nk) = NumIncGlobalCellbyNode(n,nkn)
            endif

            NumNodeFaceLocalbyCell(k,nk) = nk ! numero local du node 
            
            
         enddo


         NbIncGlobalbyCell(k) = NbNodebyCell(k)


         
c     Inc Ksigma pour les faces frac n et sigma pour les faces non frac 


         do nk=1,NbFacebyCell(k)

            n = NumFacebyCell(k,nk)


            if (NumCellbyFace(n,2).le.0) then ! face de bord 

               k1 = NumCellbyFace(n,1)
               if (k1.eq.k) then


                  if ( IndBulleNumCellbyFace(n,1).eq.1) then

                     NbIncGlobalbyCell(k) = NbIncGlobalbyCell(k)+1
                  
                     NumIncGlobalbyCell(k,NbIncGlobalbyCell(k)) =
     &                    NumIncGlobalCellbyFace(n,1)

                     NumNodeFaceLocalbyCell(k,NbIncGlobalbyCell(k))= nk ! numero local de la face 
                  endif
                  
               else
                  write(*,*)' on ne trouve pas la maille ',k,n
                  stop
               endif
               
            else

               k1 = NumCellbyFace(n,1)
               k2 = NumCellbyFace(n,2)
               
               if (k1.eq.k) then

                  if (IndBulleNumCellbyFace(n,1).eq.1) then

                    NbIncGlobalbyCell(k) = NbIncGlobalbyCell(k)+1
                      
                    NumIncGlobalbyCell(k,NbIncGlobalbyCell(k)) =
     &                   NumIncGlobalCellbyFace(n,1)

                    NumNodeFaceLocalbyCell(k,NbIncGlobalbyCell(k))= nk ! numero local de la face
                  
                  endif
                  
               else if (k2.eq.k) then

                  if (IndBulleNumCellbyFace(n,2).eq.1) then

                     NbIncGlobalbyCell(k) = NbIncGlobalbyCell(k)+1                     
                  
                     NumIncGlobalbyCell(k,NbIncGlobalbyCell(k)) =
     &                    NumIncGlobalCellbyFace(n,2)
                     
                     NumNodeFaceLocalbyCell(k,NbIncGlobalbyCell(k))= nk ! numero local de la face 

                  endif
                  
               else
                  write(*,*)' on ne trouve pas la maille ',k,n
                  stop                  
               endif
               
               
            endif

            
         enddo

                     
         
         
      enddo
c
c

        deallocate(NumIncGlobalCellbyNode)
        deallocate(NumIncGlobalCellbyFace)
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       return
       end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Operateur gradient local par maille VEM (avec tous les nodes et toutes les faces) 
c
c     Grad v = sum_s\in VK  vKs*Gradvk(k,is,:) + sum_sig\in FK  vKsig*Gradvk(k,nnodecell+isig,:) 
c
c     On met ici toutes les dof faces, on selectionnera plus tard 
c
c
c
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      
      subroutine ComputegradvKnonplanar(
     &     NbCell,NbFace,NbNode,
     &     SurfaceFace,XNode,NbArete,
     &     NbNodebyFace,NumNodebyFace,
     &     NbNodebyCell,NumNodebyCell,       
     &     NbFacebyCell,NumFacebyCell,
     &     PoidsXFaceCG,
     &     VecNormalKSigma,
     &     VolCell,
     &     gradvK,
     &     NbAretebyFace,
     &     NumAretebyFace,
     &     XCellCG,
     &     NumNodebyArete)


c
c
	implicit double precision (a-h,o-z)
	implicit integer (i-n)
c
c
        include 'include_parameter'
c
c
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
c
c     coordonnees des noeuds XS 
c
        dimension XNode(NbNode,NbDim)          
c
c     nb de noeuds par maille K 
c
        dimension NbNodebyCell(NbCell)
c
c     num des noeuds de chaque maille K dans l'ordre cyclique 
c
        dimension NumNodebyCell(NbCell,NbNodebyCellMax)
c
c
c        
c      La surface de chaque face
       dimension SurfaceFace(NbFace)

c      le nombre de noeuds de chaque face
       dimension NbNodebyFace(NbFace)
       dimension NumNodebyFace(NbFace,NbNodebyFaceMax)       


c      le nombre de faces de chaque cellule
       dimension NbFacebyCell(NbCell)

c      les numeros globales des faces de chaque cellule
       dimension NumFacebyCell(NbCell,NbFacebyCellMax)

c      Poids barycentriques du CG de la face ds l'ordre des nodes by face 
       dimension PoidsXFaceCG(NbFace,NbNodebyFaceMax)

c      vecteur normal VecKsigma oriente sortant de K
       dimension VecNormalKSigma(NbCell,NbFacebyCellMax,NbDim)

c      Volume des mailles  
       dimension VolCell(NbCell)

c      Nombre d'aretes par face
       dimension NbAretebyFace(NbFace)

c      Numero des aretes par face
       dimension NumAretebyFace(NbFace,NbAretebyFaceMax)

c     coordonnee du centre de gravite de la maille     
c
        dimension XCellCG(NbCell,NbDim)

c     2 noeuds de chaque arete 
c
        dimension NumNodebyArete(NbArete,2)





cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     SORTIES

       dimension gradvK(NbCell,NbdofbyCellMax,NbDim)             
c    
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
c
c
c
cccccccccccccccccccccccccccccccccccc      
c
c     WORKSPACE
c
       dimension Gd(NbDim), Ed(NbDim) 

cc     vecteut normal faces triangulaires
       dimension VecnormT(NbDim)
c
c     isobarycentre des nodes pour la face 
c       
       dimension Xf(NbDim) 
c
ccccccccccccccccccccccccccc
       do nCell=1,NbCell
cccccccccccccccccccccccccc

       
          nfk = NbFacebyCell(nCell)
          nnk = NbNodebyCell(nCell)
         
          gradvK(nCell,:,:) = 0.d0 
         
cccccccccccccccc         
          do nFace = 1,nfk
ccccccccccccccc
             
             Numface = NumFacebyCell(nCell,nFace)
             nnf     = NbNodebyFace(Numface)

             naf = NbAretebyFace(Numface)
c
c     calcul de l'isobarycentre des nodes de la face: Xf 
c
             Xf(:) = 0.d0 
             do nsk = 1,nnf
                ns = NumNodebyFace(Numface,nsk)                
                Xf(:) = Xf(:) + XNode(ns,:)
             enddo 
             Xf(:) = Xf(:)/dfloat(nnf)
 
             

             do ia = 1,naf  !! boucle sur les aretes

                 numa = NumAretebyFace(Numface,ia)


                 numn1 = NumNodebyArete(numa,1) 
                 numn2 = NumNodebyArete(numa,2)

                 x1 = XNode(numn1,1)
                 y1 = XNode(numn1,2)
                 z1 = XNode(numn1,3)

                 x2 = XNode(numn2,1)
                 y2 = XNode(numn2,2)
                 z2 = XNode(numn2,3)

                 x3 = Xf(1)
                 y3 = Xf(2)
                 z3 = Xf(3)


                 x4 = XCellCG(nCell,1)
                 y4 = XCellCG(nCell,2)
                 z4 = XCellCG(nCell,3)



               call VecNormalT(x1,y1,z1,x2,y2,z2,x3,y3,z3,
     &                         x4,y4,z4,vx,vy,vz,surf)   

    
            
                VecnormT(1) = vx     
                VecnormT(2) = vy      
                VecnormT(3) = vz       





               do nnode=1,2

                 NumNode = NumNodebyArete(numa,nnode)  


                icell = - 1 
                do i=1,nnk
                   if(NumNodebyCell(nCell,i).eq.NumNode) then
                      icell = i 
                   endif  
                enddo
                if (icell.eq.-1) then
                   write(*,*)' on ne trouve pas le node
     &                         de la face ds GradvK ',
     &                  ncell,NumFace
                   stop
                endif



                gradvK(nCell,icell,:) = gradvK(nCell,icell,:) !!!! pour le terme (|T|/3|K|)*(us_1 + us_2)*n_{T}
     &               +(surf/3.d0)
     &               *(1.d0/VolCell(nCell))
     &               *VecnormT(:)     






         enddo  !! fin de la boucle sur les noeuds de l'arete


         do nnode=1,nnf !! boucel sur les noeuds de la face courante  

              NumNode = NumNodebyFace(Numface,nnode)


              icell = - 1 
                do i=1,nnk
                   if(NumNodebyCell(nCell,i).eq.NumNode) then
                      icell = i 
                   endif  
                enddo
                if (icell.eq.-1) then
                   write(*,*)' on ne trouve pas le node
     &                         de la face ds GradvK ',
     &                  ncell,NumFace
                   stop
                endif



                gradvK(nCell,icell,:) = gradvK(nCell,icell,:) !!!! pour le terme (|T|/3|K|)*u_{sigma}*n_{T}
     &               +(surf/3.d0)
     &               *(1.d0/VolCell(nCell))
     &               *(1.d0/dfloat(nnf))
     &               *VecnormT(:)     

             enddo



          enddo !! fin de la bouble sur les aretes



             gradvK(nCell,nFace+nnk,:) = gradvK(nCell,nFace+nnk,:)  !!! pour le le terme bulle      
     &            + SurfaceFace(Numface)/(VolCell(nCell))
     &            *VecNormalKSigma(nCell,nFace,:)


cccccccccccccc          
          enddo   ! fin boucle face de la maille nCell 
cccccccccccccc




c
c     Test du gradient sur une fonction affine ( vKsigma = 0 ici ) 
c         
c     v = a1 + a2*x + a3*y + a4*z 
c
          a1 = 2.d0
          a2 = -3.d0
          a3 = 4.d0
          a4 = -5.d0 

c      
c    
c       
       Gd(:) = 0.d0 
       do nnode = 1,nnk

          n = NumNodebyCell(nCell,nnode)

          vKs = a1 + a2*XNode(n,1)  + a3*XNode(n,2)  + a4*XNode(n,3) 
          
          Gd(:) = Gd(:) + gradvK(nCell,nnode,:)*vKs 
          
       enddo

       Ed(1) = Gd(1) - a2 
       Ed(2) = Gd(2) - a3 
       Ed(3) = Gd(3) - a4 

       erreur = sqrt(Ed(1)**2 + Ed(2)**2 + Ed(3)**2)

       if (erreur.gt.1.d-7) then 
          write(*,*)' erreur gradvk fct affine ',nCell,erreur,Gd(:)
          write(*,*)' nnode nface ',nnk,nfk 
          stop
       endif
cccccccccccccccccc
       enddo                    ! fin boucle cell 
cccccccccccccccc
  
       

ccccccccccccccccccccccccccccccccccccccccc      
       return
       end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Operateur gradient local par maille VEM (avec tous les nodes et toutes les faces) 
c
c     Grad v = sum_s\in VK  vKs*Gradvk(k,is,:) + sum_sig\in FK  vKsig*Gradvk(k,nnodecell+isig,:) 
c
c     On met ici toutes les dof faces, on selectionnera plus tard 
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      
      subroutine ComputegradvK(NbCell,NbFace,NbNode,
     &     SurfaceFace,XNode,XFaceCG,
     &     NbNodebyFace,NumNodebyFace,
     &     NbNodebyCell,NumNodebyCell,       
     &     NbFacebyCell,NumFacebyCell,
     &     PoidsXFaceCG,
     &     VecNormalKSigma,
     &     VolCell,
     &     gradvK)

c
c
	implicit double precision (a-h,o-z)
	implicit integer (i-n)
c
c
        include 'include_parameter'
c
c
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
c
c     coordonnees des noeuds XS 
c
        dimension XNode(NbNode,NbDim)        
c
c
c
c     coordonnee du centre de gravite de la Face     
c
        dimension XFaceCG(NbFace,NbDim)        
c        
c     nb de noeuds par maille K 
c
        dimension NbNodebyCell(NbCell)
c
c     num des noeuds de chaque maille K dans l'ordre cyclique 
c
        dimension NumNodebyCell(NbCell,NbNodebyCellMax)
c
c
c        
c      La surface de chaque face
       dimension SurfaceFace(NbFace)

c      le nombre de noeuds de chaque face
       dimension NbNodebyFace(NbFace)
       dimension NumNodebyFace(NbFace,NbNodebyFaceMax)       


c      le nombre de faces de chaque cellule
       dimension NbFacebyCell(NbCell)

c      les numeros globales des faces de chaque cellule
       dimension NumFacebyCell(NbCell,NbFacebyCellMax)

c      Poids barycentriques du CG de la face ds l'ordre des nodes by face 
       dimension PoidsXFaceCG(NbFace,NbNodebyFaceMax)

c      vecteur normal VecKsigma oriente sortant de K
       dimension VecNormalKSigma(NbCell,NbFacebyCellMax,NbDim)

c      Volume des mailles  
       dimension VolCell(NbCell)



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     SORTIES

       dimension gradvK(NbCell,NbdofbyCellMax,NbDim)             
c    
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
c
c
c
cccccccccccccccccccccccccccccccccccc      
c
c     WORKSPACE
c
       dimension Gd(NbDim), Ed(NbDim) 
c
c       
c
ccccccccccccccccccccccccccc
       do nCell=1,NbCell
cccccccccccccccccccccccccc

       
          nfk = NbFacebyCell(ncell)
          nnk = NbNodebyCell(nCell)


          if (nfk+nnk.gt.NbdofbyCellMax) then
             write(*,*)' redim NbdofbyCellMax a >=',nfk+nnk
             stop
          endif
          
          gradvK(nCell,:,:) = 0.d0 
         
cccccccccccccccc         
          do nFace = 1,nfk
ccccccccccccccc
             
             Numface = NumFacebyCell(nCell,nFace)
             nnf     = NbNodebyFace(Numface)


             do nnode=1,nnf

                NumNode = NumNodebyFace(Numface,nnode)

                icell = - 1 
                do i=1,nnk
                   if(NumNodebyCell(nCell,i).eq.NumNode) then
                      icell = i 
                   endif  
                enddo
                if (icell.eq.-1) then
                   write(*,*)' on ne trouve pas le node
     &                         de la face ds GradvK ',
     &                  ncell,NumFace
                   stop
                endif
               
                gradvK(nCell,icell,:) = gradvK(nCell,icell,:)
     &               + PoidsXFaceCG(Numface,nnode)
     &               *SurfaceFace(Numface)/(VolCell(nCell))
     &               *VecNormalKSigma(nCell,nFace,:)               
             enddo

            
 
             gradvK(nCell,nFace+nnk,:) = gradvK(nCell,nFace+nnk,:)            
     &            + SurfaceFace(Numface)/(VolCell(nCell))
     &            *VecNormalKSigma(nCell,nFace,:)
             
             
cccccccccccccc          
          enddo   ! fin boucle face de la maille nCell 
cccccccccccccc

c
c     Test du gradient sur une fonction affine ( vKsigma = 0 ici ) 
c         
c     v = a1 + a2*x + a3*y + a4*z 
c
          a1 = 2.d0
          a2 = -3.d0
          a3 = 4.d0
          a4 = -5.d0 

c      
c    
c       
       Gd(:) = 0.d0 
       do nnode = 1,nnk

          n = NumNodebyCell(nCell,nnode)

          vKs = a1 + a2*XNode(n,1)  + a3*XNode(n,2)  + a4*XNode(n,3) 
          
          Gd(:) = Gd(:) + gradvK(nCell,nnode,:)*vKs 
          
       enddo

       Ed(1) = Gd(1) - a2 
       Ed(2) = Gd(2) - a3 
       Ed(3) = Gd(3) - a4 

       erreur = sqrt(Ed(1)**2 + Ed(2)**2 + Ed(3)**2)

       if (erreur.gt.1.d-10) then 
          write(*,*)' erreur gradvk fct affine ',nCell,erreur,Gd(:)
          stop
       endif
cccccccccccccccccc
       enddo                    ! fin boucle cell 
cccccccccccccccc
  
       

ccccccccccccccccccccccccccccccccccccccccc      
       return
       end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccsccccccccccccccccccccc
c
c Calcul du projecteur piK en utilisant le resultat de la subroutine "ComputegradvK"
c     piK = gradvK*(x-xK) + vK
c
c     maille k donnee = nCell 
c
c     vK = sum_ s \in Vk wks * vks, les wks donnes dans PoidsXCellXG 
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      subroutine ComputepiK(NbCell,nCell,NbNode,  
     &     NbNodebyCell,NumNodebyCell,           
     &     NbFacebyCell,XNode,
     &     XCellCG,PoidsXCellCG,
     &     gradvK,
     &     evalx,
     &     piK)
      
c     
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c     
c     
      include 'include_parameter'
c     
c     
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
      
c     nb de noeuds par maille K 
c     
      dimension NbNodebyCell(NbCell)
c     
c     num des noeuds de chaque maille K dans l'ordre cyclique 
c     
      dimension NumNodebyCell(NbCell,NbNodebyCellMax)
c     
c     le nombre de faces de chaque cellule
      dimension NbFacebyCell(NbCell)
      
c     cordonnees de nodes
      dimension XNode(NbNode,NbDim)
      
c     les piK son évalué en "evalx"  
      dimension evalx(NbDim)
      
c     le centre de mailles (centre de gravité)       
      dimension XCellCG(NbCell,NbDim)
      
c     les poids associé aux centre de gravité pour chaque cellule K
      
      dimension PoidsXCellCG(NbCell,NbNodebyCellMax)
      
c     pour chaque cellule K, gradvK(nCell,:,:) est une matrice qui contient les valeures moyennes des gradients des fonctions bases locale associé aux dof locales de K 
      dimension gradvK(NbCell,NbdofbyCellMax,NbDim)             
      
      
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     SORTIE >>> pi^Kv(evalx) pour toutes les fonction des bases "v" dans K associé aux degres de libertés dans K = ncell donnee 
c
c     piK vecteur de taille NbNodebyCell(k) + NbFacebyCell(k) : ie on met toutes les faces 
c      
c     piK U(x) = sum_i \in dofk  PiK(i)*Ui  
c      
      dimension piK(NbdofbyCellMax)             
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   
         
      nfk = NbFacebyCell(nCell)
      nnk = NbNodebyCell(nCell)
         

      do nnode=1,nnk
            
            
         piK(nnode) = prodscal(gradvK(nCell,nnode,:),
     &        evalx - XCellCG(nCell,:))
     &        + PoidsXCellCG(nCell,nnode) 
             
            
      enddo
         
         
      do nface=1,nfk 
            
            
         piK(nnk + nface) =  prodscal(
     &        gradvK(nCell,nnk + nface,:),
     &        evalx - XCellCG(nCell,:)) 
         
         
      enddo
               

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c     Test du piK sur une fonction affine (vKsigma = 0 ici)
c     
c     v = a1 + a2*x + a3*y + a4*z 
c              
      a1 = 2.d0
      a2 = -3.d0
      a3 = 4.d0
      a4 = -5.d0      
c     
c     
      piKv = 0.d0
      
      do nnode = 1,nnk
         
         n = NumNodebyCell(nCell,nnode)
         
         vKs = a1+a2*XNode(n,1)+a3*XNode(n,2)+a4*XNode(n,3)         
         
         piKv = piKv + piK(nnode)*vKs 
         
         
      enddo
               

      piKvex = a1+a2*evalx(1)+a3*evalx(2)+a4*evalx(3)         

      
      Ed = piKv - piKvex 
      
      
      erreur = dabs(Ed)
      if (erreur.ge.1.d-7) then 
         write(*,*)' erreur piK fct affine ',nCell,erreur
         stop
      endif
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      
      
      
ccccccccccccccccccccccccccccccccccccccccc
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c
c




cccc
cccc
cccc
cccc Calcul de pour chaque faces fracture on calcul 1/2(sigma u_1 + sigma u_2)
cccc Tels que sigma u_1 (respectivement sigma u_2) est une valeure constante de sigma dans k1 (respectivement k2)
cccc k1 et k2 sont les deux cellules (du coté 1 et 2 respectivement) qui partagent la face fracture numface

      subroutine StressFracMoy(NbCell,
     &     NbFace,
     &     NbdofMecaContact,
     &     NumCellbyFace,
     &     SolU,
     &     GradCellVEM,
     &     VecNormalbyFace,
     &     VecTanFace1,
     &     VecTanFace2,
     &     NbIncGlobalbyCell,
     &     NumIncGlobalbyCell,
     &     numface,
     &     StressFracLoc_n,
     &     StressFracLoc_t1,
     &     StressFracLoc_t2)
      
c     
c     
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
c
      include 'include_parameter'
c
c
c

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     DECLARATIONS
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      

cccc  ENTRE


      dimension SolU(NbdofMecaContact,NbDim)
      
      dimension GradCellVEM(NbCell,NbdofbyCellMax,NbDim) 

      dimension NumCellbyFace(NbFace,2)
      
      dimension VecNormalbyFace(NbFace,NbDim)

      dimension VecTanFace1(NbFace,NbDim)

      dimension VecTanFace2(NbFace,NbDim)

      dimension NbIncGlobalbyCell(NbCell)

      dimension NumIncGlobalbyCell(NbCell,NbIncbyCellMax)                
      

cccc  SORTI

cccc  StressFracLoc_n et StressFracLoc_t1,2
      
ccc   


ccc   workspace
      
      dimension sigma_k1(NbDim,NbDim)
      dimension sigma_k2(NbDim,NbDim) 
      
      dimension sigma_mean(NbDim,NbDim) 
      
      dimension Resultat(NbDim)



      rmu = f_Lame_mu()
      rlambda = f_Lame_lambda()



      k1 = NumCellbyFace(numface,1)
      k2 = NumCellbyFace(numface,2)

      
      
      sigma_k1(:,:) = 0.d0
      sigma_k2(:,:) = 0.d0

      Resultat(:) = 0.d0






cccc  calcul de 2 mu (1/2(gradu + gradu^T))

      do i=1,NbDim
         do j=1,NbDim
            
           

            do ik=1,NbIncGlobalbyCell(k1) 
               inc = NumIncGlobalbyCell(k1,ik)
               
               sigma_k1(i,j) = sigma_k1(i,j)
     &              + rmu*(SolU(inc,i)*GradCellVEM(k1,ik,j)    
     &              + SolU(inc,j)*GradCellVEM(k1,ik,i))

            enddo

            



            do ik=1,NbIncGlobalbyCell(k2)
               inc = NumIncGlobalbyCell(k2,ik)

               sigma_k2(i,j) = sigma_k2(i,j)
     &              + rmu*(SolU(inc,i)*GradCellVEM(k2,ik,j)
     &              + SolU(inc,j)*GradCellVEM(k2,ik,i))

            enddo
            
            
            
         enddo
      enddo




cccc  Ajouter le terme lambda*div(u)*Id

      do n=1,NbDim

         
         do i=1,NbDim
            
            do ik=1,NbIncGlobalbyCell(k1)
               inc = NumIncGlobalbyCell(k1,ik)
               
               
               sigma_k1(n,n) = sigma_k1(n,n) 
     &              + rlambda*SolU(inc,i)*GradCellVEM(k1,ik,i)  


            enddo
            
            
            do ik=1,NbIncGlobalbyCell(k2)
               inc = NumIncGlobalbyCell(k2,ik)


               sigma_k2(n,n) = sigma_k2(n,n) 
     &              + rlambda*SolU(inc,i)*GradCellVEM(k2,ik,i)  


            enddo
            
         enddo
         
         
         
      enddo



      sigma_mean(:,:) = 0.5d0*(sigma_k1(:,:)
     &                     + sigma_k2(:,:) )




      call ProdMatVec(sigma_mean,
     &     VecNormalbyFace(numface,:),Resultat)





      StressFracLoc_n = prodscal(Resultat,
     &     VecNormalbyFace(numface,:))


      StressFracLoc_t1 = prodscal(Resultat,
     &     VecTanFace1(numface,:))


      StressFracLoc_t2 = prodscal(Resultat,
     &     VecTanFace2(numface,:))



      return
      end



c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      
      subroutine ComputeSautbyFaceFrac(
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
c
c
	implicit double precision (a-h,o-z)
	implicit integer (i-n)
c
c
        include 'include_parameter'
c
c
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
c     numero de la face pour une face frac donnee  
c        
        dimension NumFaceFracVersFace(NbFaceFrac)
c
c     mailles voisines d'une face 
c        
        dimension NumCellbyFace(NbFace,2)
c
c        indice bulle dans la structure NumCellbyFace (cote maille 1 / 2) 
c
        dimension IndBulleNumCellbyFace(NbFace,2)

        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     SORTIES

        dimension NbIncGlobalbyFaceFrac(NbFaceFrac)
        dimension NumIncGlobalbyFaceFrac(NbFaceFrac,2*NbNodebyFaceMax+2)           

        dimension SautbyFaceFrac(NbFaceFrac,2*NbNodebyFaceMax+2)                 
c    
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
c
c
        do i=1,NbFaceFrac

           nf = NumFaceFracVersFace(i)

           k1 = NumCellbyFace(nf,1)
           k2 = NumCellbyFace(nf,2)

c           write(*,*)' face frac ',i,nf,k1,k2

           
           if (k2.le.0) then
              write(*,*)' face frac au bord ',i,nf,k1,k2
              stop
           endif

           nsk1 = NbNodebyCell(k1)
           nsk2 = NbNodebyCell(k2)
           nsf  = NbNodebyFace(nf)

c           write(*,*)' nb nodes ',nsf,nsk1,nsk2
           
           ii = -1
           do isf = 1,nsf 
              numsf = NumNodebyFace(nf,isf) ! nodes de la face nf 
              
              do ik=1,nsk1 
                 numik = NumNodebyCell(k1,ik)
                 if (numik.eq.numsf) then 
                    numsflocalcell1 = ik ! num local maille k1 du node numsf 
                    ii = 1
                 endif
              enddo
              if (ii.eq.-1) then 
                 write(*,*)'on ne trouve pas le noeud cote k1 ',k1,numsf
                 stop
              endif

              do ik=1,nsk2 
                 numik = NumNodebyCell(k2,ik)
                 if (numik.eq.numsf) then 
                    numsflocalcell2 = ik ! num local maille k2 du node numsf 
                    ii = 1
                 endif
              enddo
              if (ii.eq.-1) then 
                 write(*,*)'on ne trouve pas le noeud cote k2 ',k2,numsf
                 stop
              endif

              NumIncGlobalbyFaceFrac(i,isf) =
     &             NumIncGlobalbyCell(k1,numsflocalcell1)

              SautbyFaceFrac(i,isf) = PoidsXFaceCG(nf,isf)
              
              NumIncGlobalbyFaceFrac(i,isf+nsf) =
     &            NumIncGlobalbyCell(k2,numsflocalcell2)

              SautbyFaceFrac(i,isf+nsf) = - PoidsXFaceCG(nf,isf)
              
              
           enddo

          
           ibulle = 0
           
           if ( IndBulleNumCellbyFace(nf,1).eq.1) then
              ii = -1
              do ik=NbNodebyCell(k1)+1,NbIncGlobalbyCell(k1)
                 ilocalk1 = NumNodeFaceLocalbyCell(k1,ik)
                 nfacek1 = NumFacebyCell(k1,ilocalk1)
                 if (nf.eq.nfacek1) then
                    numnflocalk1 = ik
                    ii = 1
                 endif
              enddo
              if (ii.eq.-1) then 
                 write(*,*)' on ne trouve pas la face nf cote k1 ',k1,nf
                 stop
              endif
c             on rajoute la bulle cote k1          
              ibulle = ibulle + 1
              NumIncGlobalbyFaceFrac(i,2*nsf+ibulle) =
     &             NumIncGlobalbyCell(k1,numnflocalk1)

              SautbyFaceFrac(i,2*nsf+ibulle) = 1.d0 
              
           endif

           if ( IndBulleNumCellbyFace(nf,2).eq.1) then
              ii = -1
              do ik=NbNodebyCell(k2)+1,NbIncGlobalbyCell(k2)
                 ilocalk2 = NumNodeFaceLocalbyCell(k2,ik)
                 nfacek2 = NumFacebyCell(k2,ilocalk2)
                 if (nf.eq.nfacek2) then
                    numnflocalk2 = ik
                    ii = 1
                 endif
              enddo
              if (ii.eq.-1) then 
                 write(*,*)' on ne trouve pas la face nf cote k2 ',k2,nf
                 stop
              endif
c             on rajoute la bulle cote k2               
              ibulle = ibulle + 1
              NumIncGlobalbyFaceFrac(i,2*nsf+ibulle) =
     &             NumIncGlobalbyCell(k2,numnflocalk2)

              SautbyFaceFrac(i,2*nsf+ibulle) = - 1.d0 
                          
           endif

           NbIncGlobalbyFaceFrac(i) = 2*nsf + ibulle

c           write(*,*)' nb inc op saut ',i,nf,NbIncGlobalbyFaceFrac(i)
           
        enddo
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
       return
       end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       




ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     Pour chaque face fracture on calcule l'operateur Traction T
c     du cote de la maille k1 : (Sigmak1 nk1)
c     avec ses 3 composantes : (Sigmak1 nk1)nk1
c                              (Sigmak1 nk1)nt1
c                              (Sigmak1 nk1)nt2
c     
c     
c
c      
c     sigmanbyff(nff,U,:) = sum_{dof,l} Sigmanbyff(nff,dof,l,:) U^dof_l 
c         
c
c     les dof sont ceux de la maille k1 de la face 
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      
      subroutine SigmanbyFaceFrac(
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
c     
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
c
      include 'include_parameter'
c
c
c

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     DECLARATIONS
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      

cccc  ENTREES
      
      dimension GradCellVEM(NbCell,NbdofbyCellMax,NbDim) 

      dimension NumCellbyFace(NbFace,2)
      
      dimension VecNormalbyFace(NbFace,NbDim)
      dimension VecTanFace1(NbFace,NbDim)
      dimension VecTanFace2(NbFace,NbDim)      

      dimension NbIncGlobalbyCell(NbCell)
      
c
c     numero de la face pour une face frac donnee  
c        
      dimension NumFaceFracVersFace(NbFaceFrac)

      
cccc  SORTIE
c
c     traction T (composantes normale et tangentielles)  
c      
      dimension Sigmanbyff(NbFaceFrac,NbdofbyCellMax,Nbdim,NbDim)
c
c     sigman_nf(U,:) = sum_{dof,l} Sigmanbyff(nf,dof,l,:) U^dof_l 
c      
ccccccccccc  


ccc   workspace: operateur sigmak1 
      
      dimension sigma_k1(NbDim,NbDim,NbdofbyCellMax,NbDim)
            
c
c      Sigma_k1_ij = sum_{dof,l} sigma_k1(i,j,dof,l)*(U^dof)_l 
c      
cccccccccccccccccccccccccccccccccccccccccccccccccc

      rmu = f_Lame_mu()
      rlambda = f_Lame_lambda()


      Sigmanbyff(:,:,:,:) = 0.d0 
      
ccccccccccccccccccccccccccccc
      do nff = 1,NbFaceFrac ! boucle sur les faces frac 
cccccccccccccccccccccccccccc
         
         numface = NumFaceFracVersFace(nff)
      

         k1 = NumCellbyFace(numface,1) ! maille cote k1 de la face 
     
         sigma_k1(:,:,:,:) = 0.d0


cccc  Ajout du terme 2 mu *1/2*( grad(.) + grad(.)^T )

         do i=1,NbDim
            do j=1,NbDim
                      

               do ik=1,NbIncGlobalbyCell(k1) 
                  
                  sigma_k1(i,j,ik,i) = sigma_k1(i,j,ik,i)
     &                 + rmu*GradCellVEM(k1,ik,j)


                  sigma_k1(i,j,ik,j) = sigma_k1(i,j,ik,j)
     &                 + rmu*GradCellVEM(k1,ik,i)                  
                  
               enddo
            
            
            enddo
         enddo




cccc  Ajout du terme lambda*div(.)*Id

         do n=1,NbDim
         
            do i=1,NbDim              
               do ik=1,NbIncGlobalbyCell(k1)
                  
               
                  sigma_k1(n,n,ik,i) = sigma_k1(n,n,ik,i) 
     &                 + rlambda*GradCellVEM(k1,ik,i)  

               enddo                                 
            enddo
                  
         enddo
c
c     composante normale 
c
         do ik=1,NbIncGlobalbyCell(k1)
            do l=1,NbDim

               s = 0.d0 
               do i=1,NbDim
                  do j=1,NbDim

               
                     s = s + VecNormalbyFace(numface,i)
     &                    *VecNormalbyFace(numface,j)
     &                    *sigma_k1(i,j,ik,l)

                     
                  enddo
               enddo

               Sigmanbyff(nff,ik,l,1) = s 

               
            enddo
         enddo
c
c     composante tangentielle 1  
c
         do ik=1,NbIncGlobalbyCell(k1)
            do l=1,NbDim

               s = 0.d0 
               do i=1,NbDim
                  do j=1,NbDim

               
                     s = s + VecTanFace1(numface,i)
     &                    *VecNormalbyFace(numface,j)
     &                    *sigma_k1(i,j,ik,l)

                     
                  enddo
               enddo

               Sigmanbyff(nff,ik,l,2) = s 

               
            enddo
         enddo
c         
c     composante tangentielle 2  
c
         do ik=1,NbIncGlobalbyCell(k1)
            do l=1,NbDim

               s = 0.d0 
               do i=1,NbDim
                  do j=1,NbDim

               
                     s = s + VecTanFace2(numface,i)
     &                    *VecNormalbyFace(numface,j)
     &                    *sigma_k1(i,j,ik,l)

                     
                  enddo
               enddo

               Sigmanbyff(nff,ik,l,3) = s 

               
            enddo
         enddo
         
ccccccccccccccccccc
      enddo 
ccccccccccccccccccc      

      return
      end
c
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     




      

      
ccccccccccccccccccccccccccccccccccccc
ccccc fonction produit scalaire ccccc
ccccccccccccccccccccccccccccccccccccc



      function prodscal(vec1,vec2)
  
c
c
  
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c     
c     
      include 'include_parameter'
c     
c     
c     
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     DECLARATIONS
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     ENTREES 
c     

      dimension vec1(NbDim)
      dimension vec2(Nbdim)     
c     
      
      prodscal  = vec1(1)*vec2(1)  + vec1(2)*vec2(2) + vec1(3)*vec2(3)
               
      
ccccccccccccccccccccccccccccccccccccccccc
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c


      subroutine prodtensor(vec1,vec2,tensor)
  
c
c
  
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c     
c     
      include 'include_parameter'
c     
c     
c          
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     DECLARATIONS
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     ENTREES 
c     

      dimension vec1(NbDim)
      dimension vec2(Nbdim)     
c 
ccccccccccccccccccccccccccccccccccccc
    
      dimension tensor(NbDim,NbDim)     
     

       do i=1,NbDim 

         do j=1,NbDim 

           tensor(i,j) =  vec1(i)*vec2(j)          

         enddo
     
       enddo
               
      
ccccccccccccccccccccccccccccccccccccccccc
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c
c
      subroutine ProdMatVec(AA,Vec,Resultat)
  
c
c
  
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c     
c     
      include 'include_parameter'
c     
c     
c          
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     DECLARATIONS
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     ENTREES 
c     

      dimension AA(NbDim,NbDim)
      dimension Vec(Nbdim)     
c 
ccccccccccccccccccccccccccccccccccccc
    
c     SORTI

      dimension Resultat(NbDim)     
     

       do i=1,NbDim 
 
         Resultat(i) = 0.d0

         do j=1,NbDim 

           Resultat(i) =  Resultat(i) + AA(i,j)*Vec(j)          

         enddo
     
       enddo
               
      
ccccccccccccccccccccccccccccccccccccccccc
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccc fonction produit scalaire entre deux matrices ccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



      function prodscalMatrix(A1,A2)
  
c
c
  
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c     
c     
      include 'include_parameter'
c     
c     
c     
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     DECLARATIONS
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     ENTREES 
c     

      dimension A1(NbDim,NbDim)
      dimension A2(Nbdim,NbDim)     
c     
      prodscalMatrix = 0.d0

      
      do i = 1,NbDim
         do j = 1,NbDim
            prodscalMatrix = prodscalMatrix + A1(i,j)*A2(i,j)
         enddo
      enddo


      
ccccccccccccccccccccccccccccccccccccccccc
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
