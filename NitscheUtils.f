c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Operateur saut SCALAIRE evalue en point X pour la frac numfacefrac: EvalXSautbyFaceFrac
c     Saut avec les reconstructions affines nodales sur la face des cotes k1 et k2 
c     obtenues avec le gradient tangent constant par face et la valeur au CG par les poids  
c      
c     (v^+ - v^-)(X) = sum_idof=1^n EvalXSautbyff(idof)*Vinc(inc) 
c     cote + = cote maille 1 de NumCellbyFace
c     inc = NumIncGlobalbyFaceFrac(numfrac,idof)
c     n = NbIncGlobalbyFaceFrac((numfrac)
c
c
c
c      
c     Numerotation globale des inc de l'op de saut pour l'ordre local suivant:
c        - nodes de la face cote maille 1 
c        - nodes de la face cote maille 2 
c        - face cote maille 1 si bulle ( cf si IndBulleNumCellbyface(n,1) = 1 )
c        - face cote maille 2 si bulle ( cf si IndBulleNumCellbyface(n,2) = 1 )     
c
c
c     On met zero sur les inc bulle (operateur nodal) 
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      
      subroutine EvalXSautbyFaceFrac(
     &     XX,numfacefrac,
     &     NbFaceFrac,NbFace,NbCell,
     &     NbNodebyFace,NumNodebyFace,
     &     NbNodebyCell,NumNodebyCell,       
     &     NbIncGlobalbyCell,
     &     NumIncGlobalbyCell,
     &     NumNodeFaceLocalbyCell,
     &     NumFaceFracVersFace,
     &     NumFacebyCell,NumCellbyFace,
     &     IndBulleNumCellbyFace,            
     &     PoidsXFaceCG,XFaceCG,
     &     GradTangbyFace, 
     &     NbIncGlobalbyFaceFrac,
     &     NumIncGlobalbyFaceFrac,
     &     EvalXSautbyff) 
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
        dimension XX(NbDim) 
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

c
c     coordonnee du centre de gravite de la Face     
c
        dimension XFaceCG(NbFace,NbDim)        

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
c
c     inconnues locales et globales de l'operateur saut 
c     Numerotation globale des inc de l'op de saut pour l'ordre local suivant:
c        - nodes de la face cote maille 1 
c        - nodes de la face cote maille 2 
c        - face cote maille 1 si bulle ( cf si IndBulleNumCellbyface(n,1) = 1 )
c        - face cote maille 2 si bulle ( cf si IndBulleNumCellbyface(n,2) = 1 )     
        
        dimension NbIncGlobalbyFaceFrac(NbFaceFrac)
        dimension NumIncGlobalbyFaceFrac(NbFaceFrac,2*NbNodebyFaceMax+2)         
        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     SORTIES
  
        dimension EvalXSautbyff(2*NbNodebyFaceMax+2)                 
c    
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
c
c
        i = numfacefrac

        
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
     &          NumIncGlobalbyCell(k1,numsflocalcell1)
           
           EvalXSautbyff(isf) = PoidsXFaceCG(nf,isf)
           do l=1,NbDim
              EvalXSautbyff(isf) = EvalXSautbyff(isf)                 
     &             + GradTangbyFace(nf,isf,l)*(XX(l) - XFaceCG(nf,l))
           enddo

              
           NumIncGlobalbyFaceFrac(i,isf+nsf) =
     &          NumIncGlobalbyCell(k2,numsflocalcell2)
           
           EvalXSautbyff(isf+nsf) = - PoidsXFaceCG(nf,isf)
           do l=1,NbDim
              EvalXSautbyff(isf+nsf) = EvalXSautbyff(isf+nsf)                 
     &             - GradTangbyFace(nf,isf,l)*(XX(l) - XFaceCG(nf,l))
           enddo              
           
        enddo
c
c     on met zero sur les bulles 
c          
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
c     on rajoute la bulle cote k1          
           ibulle = ibulle + 1
           NumIncGlobalbyFaceFrac(i,2*nsf+ibulle) =
     &          NumIncGlobalbyCell(k1,numnflocalk1)

           EvalXSautbyff(2*nsf+ibulle) = 0.d0 
           
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
     &          NumIncGlobalbyCell(k2,numnflocalk2)
           
           EvalXSautbyff(2*nsf+ibulle) = 0.d0 
           
        endif
        
        NbIncGlobalbyFaceFrac(i) = 2*nsf + ibulle

c           write(*,*)' nb inc op saut ',i,nf,NbIncGlobalbyFaceFrac(i)
           
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
       return
       end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       




cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c
c     Operateur GradTangbyFace calculant le gradient tangentiel nodal constant par face 
c
c     G_l(nf) = sum_is GradTangbyFace(nf,is,l) U(inc)
c
c     is: numero local a la face du noeud, numerote dans l'ordre des noeuds de la face  
c     inc = inc globale du node is de la face 
c
c
c     On suppose l'ordre cyclique des noeuds de la face 
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c      
      subroutine GradTangentbyFace(NbNode,NbFace,
     &     XNode,
     &     NbNodebyFace,NumNodebyFace,
     &     XFaceCG,SurfaceFace,
     &     GradTangbyFace) 
     
c
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     FORMAT BENCHMARK 3D 
c
ccccccccccccccccccccccccccccccccccccccc
c
c     coordonnees des noeuds XS 
c
        dimension XNode(NbNode,NbDim)

c
c
c     Faces par les noeuds 
c       
        dimension NbNodebyFace(NbFace)
        dimension NumNodebyFace(NbFace,NbNodebyFaceMax)

c
c     coordonnee du centre de gravite de la Face     
c
        dimension XFaceCG(NbFace,NbDim)

c
c     surface by face 
c        
        dimension SurfaceFace(NbFace)     

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     SORTIES
        
c
c     Operateur gradient tangentiel constant par face 
c
c     Grad(nf)_l = sum_{is \in V_nf} GradTangbyFace(nf,is,l)*U(inc)
c     inc = num global du node is de la face nf
c     l = 1,Nbdim 
c        
c        
        dimension GradTangbyFace(NbFace,NbNodebyFaceMax,NbDim)        

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        dimension vecn(NbDim) 

        
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
c
c
        GradTangbyFace(:,:,:) = 0.d0 
        
        do numf = 1,NbFace

           xf = XFaceCG(numf,1)
           yf = XFaceCG(numf,2)
           zf = XFaceCG(numf,3)


           surf = SurfaceFace(numf)

           
           nsf = NbNodebyFace(numf)
           do is = 1,nsf

              is1 = is
              if (is.lt.nsf) then 
                 is2 = is+1
              else
                 is2 = 1
              endif


              ns1 = NumNodebyFace(numf,is1)
              ns2 = NumNodebyFace(numf,is2)


                 x1 = XNode(ns1,1)
                 y1 = XNode(ns1,2)
                 z1 = XNode(ns1,3)

                 x2 = XNode(ns2,1)
                 y2 = XNode(ns2,2)
                 z2 = XNode(ns2,3)              

c
c     vecteur normal (vx,vy,vz) a l'arete 12 ds le plan 12f sortant de la face 
c
c     n = 1-f + alpha*(2-1)  avec n.(2-1)=0  
c                 

                 s1f = (x1-xf)*(x2-x1) 
     &                + (y1-yf)*(y2-y1) + (z1-zf)*(z2-z1)
                 
                 s12 = (x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2

                 d12 = dsqrt(s12) ! longueur de l'arete s1s2

                 alpha = -s1f/s12

                 vx = x1-xf + alpha*(x2-x1)
                 vy = y1-yf + alpha*(y2-y1)
                 vz = z1-zf + alpha*(z2-z1)

                 s = dsqrt(vx**2+vy**2+vz**2) 
c
c     vecteur normal sortant 
c                 
                 vx = vx/s
                 vy = vy/s
                 vz = vz/s
                 
                 vecn(1) = vx
                 vecn(2) = vy
                 vecn(3) = vz 


                 GradTangbyFace(numf,is1,:) = GradTangbyFace(numf,is1,:) 
     &                + d12/surf/2.d0*vecn(:)
          

                 GradTangbyFace(numf,is2,:) = GradTangbyFace(numf,is2,:) 
     &                + d12/surf/2.d0*vecn(:)
          
                 
           enddo

           
        enddo



        
        
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
       return
       end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     En entree :
c
c      * theta et beta de la face fracture numfacefrac  
c
c      * point XX
c
c      * Sigmanbyff 
c
c     En sortie: operateur calculant le vecteur theta*sigman(U) - beta[U]
c     (composantes normale puis tangentielles) 
c
c     P_theta_beta(U,XX,:) = sum_{idof,l} EvalXPThetaBetabyff(idof,l,:)*U^dof_l 
c      
c     On numerote les inconnues locales dans l'ordre suivant:
c        inc maille k1       
c        nodes face sigma cote k2 puis bulle (si presente) cote k2       
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      
      subroutine EvalXPThetaBetabyFaceFrac(
     &     XX,numfacefrac,theta,betan,betat,
     &     NbFaceFrac,NbFace,NbCell,
     &     NbNodebyFace,NumNodebyFace,
     &     NbNodebyCell,NumNodebyCell,       
     &     NbIncGlobalbyCell,
     &     NumIncGlobalbyCell,
     &     NumNodeFaceLocalbyCell,
     &     NumFaceFracVersFace,
     &     NumFacebyCell,NumCellbyFace,
     &     IndBulleNumCellbyFace,            
     &     PoidsXFaceCG,XFaceCG,
     &     GradTangbyFace, 
     &     NbIncGlobalbyFaceFrac,
     &     NumIncGlobalbyFaceFrac,
     &     Sigmanbyff,
     &     VecNormalbyFace,
     &     VecTanFace1,
     &     VecTanFace2,          
     &     EvalXPThetaBetabyff,
     &     NbIncGlobalPThetaBeta,
     &     NumIncGlobalPThetaBeta) 
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
        dimension XX(NbDim) 
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

c
c     coordonnee du centre de gravite de la Face     
c
        dimension XFaceCG(NbFace,NbDim)        

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

c
c     inconnues locales et globales de l'operateur saut 
c     Numerotation globale des inc de l'op de saut pour l'ordre local suivant:
c        - nodes de la face cote maille 1 
c        - nodes de la face cote maille 2 
c        - face cote maille 1 si bulle ( cf si IndBulleNumCellbyface(n,1) = 1 )
c        - face cote maille 2 si bulle ( cf si IndBulleNumCellbyface(n,2) = 1 )     
        
        dimension NbIncGlobalbyFaceFrac(NbFaceFrac)
        dimension NumIncGlobalbyFaceFrac(NbFaceFrac,2*NbNodebyFaceMax+2)
c
c
c     operateur sigman by face frac (tractions normale et tangentielles) 
c     sigman_nf(U,:) = sum_{dof,l} Sigmanbyff(nf,dof,l,:) U^dof_l    
c
c        
      dimension Sigmanbyff(NbFaceFrac,NbdofbyCellMax,NbDim,NbDim)

c
c     vecteurs normal et tangents par face (normale sortante maille 1 de NumCellbyFace) 
c     
      dimension VecNormalbyFace(NbFace,NbDim)
      dimension VecTanFace1(NbFace,NbDim)
      dimension VecTanFace2(NbFace,NbDim)      
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     SORTIES
c
c     theta*Sigmanbyff - beta*EvalXSautbyff
c
c     On numerote les inconnues locales dans l'ordre suivant:
c        inc maille k1       
c        nodes face sigma cote k2 puis bulle (si presente) cote k2 
c      
      dimension EvalXPThetaBetabyff(
     &     NbdofbyCellMax+NbNodebyFaceMax+1,NbDim,NbDim)        
c
c     numero global des inconnues de l'operateur
c     On numerote les inconnues locales dans l'ordre suivant:
c        inc maille k1       
c        nodes face sigma cote k2 puis bulle (si presente) cote k2        
c      
c      NbIncGlobalPThetaBeta
      
      dimension NumIncGlobalPThetaBeta(NbdofbyCellMax+NbNodebyFaceMax+1)      
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc         
c
c     Workspace 
c
        dimension EvalXSautbyff(2*NbNodebyFaceMax+2)        

        
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
c
c
      call EvalXSautbyFaceFrac(
     &     XX,numfacefrac,
     &     NbFaceFrac,NbFace,NbCell,
     &     NbNodebyFace,NumNodebyFace,
     &     NbNodebyCell,NumNodebyCell,       
     &     NbIncGlobalbyCell,
     &     NumIncGlobalbyCell,
     &     NumNodeFaceLocalbyCell,
     &     NumFaceFracVersFace,
     &     NumFacebyCell,NumCellbyFace,
     &     IndBulleNumCellbyFace,            
     &     PoidsXFaceCG,XFaceCG,
     &     GradTangbyFace, 
     &     NbIncGlobalbyFaceFrac,
     &     NumIncGlobalbyFaceFrac,
     &     EvalXSautbyff) 

      

          EvalXPThetaBetabyff(:,:,:) = 0.d0
      
 
          nf = NumFaceFracVersFace(numfacefrac)
          k1 = NumCellbyFace(nf,1)
          k2 = NumCellbyFace(nf,2)

          do ik = 1,NbIncGlobalbyCell(k1)

             EvalXPThetaBetabyff(ik,:,:) =
     &            theta*Sigmanbyff(numfacefrac,ik,:,:)

             NumIncGlobalPThetaBeta(ik) = NumIncGlobalbyCell(k1,ik)
             
          enddo


          i = 0
          do ik=1,NbIncGlobalbyFaceFrac(numfacefrac)             
            inc = NumIncGlobalbyFaceFrac(numfacefrac,ik)

            ind = -1
            do ik1 = 1,NbIncGlobalbyCell(k1)
               inck1 = NumIncGlobalbyCell(k1,ik1)
               if (inck1.eq.inc) then
                  ind = 1
                  ikcell = ik1
               endif
            enddo

            

            if (ind.eq.1) then

               EvalXPThetaBetabyff(ikcell,:,1) =
     &              EvalXPThetaBetabyff(ikcell,:,1)
     &              - betan*EvalXSautbyff(ik)*VecNormalbyFace(nf,:)

               EvalXPThetaBetabyff(ikcell,:,2) =
     &              EvalXPThetaBetabyff(ikcell,:,2)
     &              - betat*EvalXSautbyff(ik)*VecTanFace1(nf,:)               


               EvalXPThetaBetabyff(ikcell,:,3) =
     &              EvalXPThetaBetabyff(ikcell,:,3)
     &              - betat*EvalXSautbyff(ik)*VecTanFace2(nf,:)                         
            else

               i = i + 1
               EvalXPThetaBetabyff(NbIncGlobalbyCell(k1)+i,:,1) =
     &              EvalXPThetaBetabyff(NbIncGlobalbyCell(k1)+i,:,1)
     &              - betan*EvalXSautbyff(ik)*VecNormalbyFace(nf,:)               


               EvalXPThetaBetabyff(NbIncGlobalbyCell(k1)+i,:,2) =
     &              EvalXPThetaBetabyff(NbIncGlobalbyCell(k1)+i,:,2)
     &              - betat*EvalXSautbyff(ik)*VecTanFace1(nf,:)    


               EvalXPThetaBetabyff(NbIncGlobalbyCell(k1)+i,:,3) =
     &              EvalXPThetaBetabyff(NbIncGlobalbyCell(k1)+i,:,3)
     &              - betat*EvalXSautbyff(ik)*VecTanFace2(nf,:)    
               
               NumIncGlobalPThetaBeta(NbIncGlobalbyCell(k1)+i) = inc 
               
            endif 
               
            


            
         enddo

         NbIncGlobalPThetaBeta = NbIncGlobalbyCell(k1) + i

         
c         write(*,*)' nb inc Pthetabeta ',NbIncGlobalPThetaBeta
c         do ik=1,NbIncGlobalPThetaBeta
c            inc = NumIncGlobalPThetaBeta(ik)
c            write(*,*)' P ',ik,inc,EvalXPThetaBetabyff(ik,:,:)
c         enddo
         
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
       return
       end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c
c     formule de quadrature sur une face numf
c
c     On coupe la face en triangles T = (arete12 - CG face)
c     et on prend les points de quadrature aux milieux
c     des aretes des triangles T avec les poids |T|/3
c     Formule exacte sur les polynomes d'ordre 2 
c      
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c      
      subroutine QuadFace(
     &     numf,
     &     NbNode,NbFace,
     &     XNode,
     &     NbNodebyFace,NumNodebyFace,
     &     XFaceCG,SurfaceFace,
     &     NbPoints,XPoints,PoidsPoints) 
     
c
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     coordonnees des noeuds XS 
c
        dimension XNode(NbNode,NbDim)

c
c
c     Faces par les noeuds 
c       
        dimension NbNodebyFace(NbFace)
        dimension NumNodebyFace(NbFace,NbNodebyFaceMax)

c
c     coordonnee du centre de gravite de la Face     
c
        dimension XFaceCG(NbFace,NbDim)

c
c     surface by face 
c        
        dimension SurfaceFace(NbFace)     

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     SORTIES
c       
c
c     NbPoints = 2 x nb de nodes de la face 
c        
        dimension PoidsPoints(2*NbNodebyFaceMax)        
        dimension XPoints(2*NbNodebyFaceMax,NbDim)        
c
c
c        
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
c
c

        PoidsPoints(:) = 0.d0 
        XPoints(:,:) = 0.d0 
        
        xf = XFaceCG(numf,1)
        yf = XFaceCG(numf,2)
        zf = XFaceCG(numf,3)
        

        surf = SurfaceFace(numf)


c        write(*,*)' surface face ',surf 
           
        nsf = NbNodebyFace(numf)


        NbPoints = 2*nsf        ! nb de points de quadrature 

        do is = 1,nsf

           is1 = is
           if (is.lt.nsf) then 
              is2 = is+1
           else
              is2 = 1
           endif

           ns1 = NumNodebyFace(numf,is1)
           ns2 = NumNodebyFace(numf,is2)


           x1 = XNode(ns1,1)
           y1 = XNode(ns1,2)
           z1 = XNode(ns1,3)

           x2 = XNode(ns2,1)
           y2 = XNode(ns2,2)
           z2 = XNode(ns2,3)

c           write(*,*)' X node ',x1,y1,z1
c
c     point de quad milieu de l'arete de bord 12 
c
           XPoints(is,1) = (x1 + x2)/2.d0 
           XPoints(is,2) = (y1 + y2)/2.d0 
           XPoints(is,3) = (z1 + z2)/2.d0 
c
c     point de quad milieu de l'arete 1f de T 
c           
           XPoints(nsf+is,1) = (x1 + xf)/2.d0 
           XPoints(nsf+is,2) = (y1 + yf)/2.d0 
           XPoints(nsf+is,3) = (z1 + zf)/2.d0            

           xT = ( x1 + x2 + xf )/3.d0 
           yT = ( y1 + y2 + yf )/3.d0 
           zT = ( z1 + z2 + zf )/3.d0 

           call VecNormalT(x1,y1,z1,x2,y2,z2,xf,yf,zf,xT,yT,zT,
     &          vx,vy,vz,surfT)
c
c     poids sur les 3 aretes 12, 1f et 2f du triangle T 
c
           PoidsPoints(is) = surfT/3.d0

           PoidsPoints(nsf+is1) = PoidsPoints(nsf+is1) + surfT/3.d0

           PoidsPoints(nsf+is2) = PoidsPoints(nsf+is2) + surfT/3.d0
          
        enddo

        

c        do is = 1,nsf
c           write(*,*)' points quad bord ',XPoints(is,:),
c     &          ' poids ',PoidsPoints(is)
c           write(*,*)' points quad int ',XPoints(nsf+is,:),
c     &          'poids ',PoidsPoints(nsf+is)                      
c        enddo
        
           

        
        
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
       return
       end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
c
c
