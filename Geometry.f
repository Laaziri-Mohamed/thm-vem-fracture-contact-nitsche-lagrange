c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                                                                                                       
c                         
c          Calcul geometriques pour le schema VEM bulle 
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
      subroutine Geometry(NbCell,NbNode,NbFace,NbArete,
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
c     ENTREES : on suppose XCell et XFace deja calcules par le schema HFV ou VAG
c
c     ce ne sont pas necessairement les centres de gravite !!! 
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
c     Faces par les aretes 
c
        dimension NbAretebyFace(NbFace)
        dimension NumAretebyFace(NbFace,NbAretebyFaceMax)
c
c     Faces par les noeuds 
c       
        dimension NbNodebyFace(NbFace)
        dimension NumNodebyFace(NbFace,NbNodebyFaceMax)
c
c     Mailles voisines des faces 
c
        dimension NumCellbyFace(NbFace,2)
c
c     2 noeuds de chaque arete 
c
        dimension NumNodebyArete(NbArete,2)
c
c     coordonnee du centre de maille XK (maille supposee etoilee par rapport a XK)
c     num cell  >  XCell
c
        dimension XCell(NbCell,NbDim)
c
c     Vol de la maille pour test 
c
        dimension VolCell(NbCell) 


        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     SORTIES 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c     coordonnee du centre de gravite de la maille     
c
        dimension XCellCG(NbCell,NbDim)
c
c     Poids barycentriques du CG de la maille ds l'ordre des nodes by cell  
c
       dimension PoidsXCellCG(NbCell,NbNodebyCellMax)        
c
c     coordonnee du centre de gravite de Face     
c
        dimension XFaceCG(NbFace,NbDim)
c
c
c     Poids barycentriques du CG de la face ds l'ordre des nodes by face 
c
       dimension PoidsXFaceCG(NbFace,NbNodebyFaceMax)

c        
c     surface de la face 
c
        dimension SurfaceFace(NbFace)
c
c     vecteur normal VecKsig oriente sortant de K
c
        dimension VecNormalKSigma(NbCell,NbFacebyCellMax,NbDim)
c
c     Longueur Arete 
c
        dimension SizeArete(NbArete)
c        
c
c     Vecteur normal VecSigmaArete oriente sortant de sigma 
c
       dimension VecNormalSigmaArete(NbFace,NbAretebyFaceMax,NbDim)

c
c     vecteur normal a la face sigma sortant de la maille 1 de NumCellbyFace 
c
       dimension VecNormalbyFace(NbFace,NbDim)
c
c       
cccccccccccccccccccccccc
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     WorkSpaces locaux 
c
ccccccccccccccccccccccccccccccccccccccccccc
c
c
c     Par maille: 
c     Surfi:  surface des faces, 
c     Xfi = XFace 
c     VecNormalKi: normale unitaire sortante  
c

        dimension Surfi(NbFacebyCellMax)
        dimension Xfi(NbFacebyCellMax,NbDim)
        dimension VecNormalKi(NbFacebyCellMax,NbDim)
c
c     par face 
c
c     VecNormalKi: normale unitaire sortante dans le plan du triangle 
c 
 
        dimension VecNormalSi(NbAretebyFaceMax,NbDim)

c
        dimension XCG(NbDim) 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
ccccccccccccccccccccccccccccc
c
c     BOUCLE SUR LES MAILLES 
c
        do k=1,NbCell
ccccccccccccccccccccccccccccc
c
c     nfk: nb de faces de la maille 
c
           nsk = NbNodebyCell(k)           
           nfk = NbFacebyCell(k)
c
c
cccccccccccccccccccccccccccc

           xk = XCell(k,1)
           yk = XCell(k,2)
           zk = XCell(k,3)


           Volk = 0.d0 
           XCellCG(k,:) = 0.d0 
           PoidsXCellCG(k,:) = 0.d0 
c
c
ccccccccccccccccccccccccccccc
c
c          boucle sur les faces de la maille k
c
           do jf = 1,nfk 
              numf = NumFacebyCell(k,jf)

              PoidsXFaceCG(numf,:) = 0.d0 ! poids du CG de la face numf 
              
              xf = 0.d0
              yf = 0.d0
              zf = 0.d0

              do nk = 1,NbNodebyFace(numf)
                 n = NumNodebyFace(numf,nk) 
                 xf = xf + XNode(n,1)
                 yf = yf + XNode(n,2)
                 zf = zf + XNode(n,3)                  
              enddo
                            
              xf = xf/dfloat(NbNodebyFace(numf))
              yf = yf/dfloat(NbNodebyFace(numf))
              zf = zf/dfloat(NbNodebyFace(numf))

              Surfi(jf) = 0.d0
              si = 0.d0
              do m=1,NbDim
                 Xfi(jf,m) = 0.d0
                 VecNormalKi(jf,m) = 0.d0
              enddo

c
c             boucle sur les aretes de la face numf 
c
              nsf = NbAretebyFace(numf)
              do ia = 1,nsf 
                 numa = NumAretebyFace(numf,ia)
c
cccccccccccccccccccccccccccc
c
c                noeuds de l'arete numa 
c
                 numn1 = NumNodebyArete(numa,1)
                 numn2 = NumNodebyArete(numa,2)

                 x1 = XNode(numn1,1)
                 y1 = XNode(numn1,2)
                 z1 = XNode(numn1,3)

                 x2 = XNode(numn2,1)
                 y2 = XNode(numn2,2)
                 z2 = XNode(numn2,3)

c
c     volume du tetra Xk,Xf,Arete X1X2
c
                 VolT = VolTetra(x1,y1,z1,x2,y2,z2,
     &                                  xf,yf,zf,xk,yk,zk)

                 
c
c     numero local des deux noeuds dans la liste des noeuds de la face 
c
                 ii = -1 
                 do is=1,nsf
                    numis = NumNodebyFace(numf,is)
                    if (numis.eq.numn1) then 
                       numn1localface = is 
                       ii = 1
                    endif
                 enddo
                 if (ii.eq.-1) then 
                    write(*,*)' on ne trouve pas le noeud ',numf,numn1
                    stop
                 endif

                 ii = -1 
                 do is=1,nsf
                    numis = NumNodebyFace(numf,is)
                    if (numis.eq.numn2) then 
                       numn2localface = is 
                       ii = 1
                    endif
                 enddo
                 if (ii.eq.-1) then 
                    write(*,*)' on ne trouve pas le noeud ',numf,numn2
                    stop
                 endif


c                numero local du noeud numn1 = num du noeud ds liste noeuds maille 
                 ii = -1 
                 do ik=1,nsk
                    numik = NumNodebyCell(k,ik)
                    if (numik.eq.numn1) then 
                       numn1localcell = ik 
                       ii = 1
                    endif
                 enddo
                 if (ii.eq.-1) then 
                    write(*,*)' on ne trouve pas le noeud ',k,numn1
                    stop
                 endif
c                 write(*,*)' numn1localcell ',numn1localcell

c                numero local du noeud numn2 = num du noeud ds liste noeuds maille                   
                 ii = -1 
                 do ik=1,nsk 
                    numik = NumNodebyCell(k,ik)
                    if (numik.eq.numn2) then 
                       numn2localcell = ik 
                       ii = 1
                    endif
                 enddo
                 if (ii.eq.-1) then 
                    write(*,*)' on ne trouve pas le noeud ',k,numn2
                    stop
                 endif
                 
c
c     centre du tetra X             
c
                 x = (x1+x2+xk+xf)/4.d0
                 y = (y1+y2+yk+yf)/4.d0
                 z = (z1+z2+zk+zf)/4.d0
c
c     X CG Cell 
c                 
                 Volk = Volk + VolT
                 XCellCG(k,1) = XCellCG(k,1) + VolT*x
                 XCellCG(k,2) = XCellCG(k,2) + VolT*y
                 XCellCG(k,3) = XCellCG(k,3) + VolT*z
c
c     Poids du CG Cell / aux nodes de la maille 
c                 
                 PoidsXCellCG(k,numn1localcell) =
     &                 PoidsXCellCG(k,numn1localcell)
     &                + VolT/4.d0

                 PoidsXCellCG(k,numn2localcell) =
     &                 PoidsXCellCG(k,numn2localcell)
     &                + VolT/4.d0                  

                 do isk=1,nsk
                    PoidsXCellCG(k,isk) =
     &                   PoidsXCellCG(k,isk)
     &                   + VolT/4.d0/dfloat(nsk)                    
                 enddo


                 do isf = 1,nsf

                    numsf = NumNodebyFace(numf,isf)
                    ii = -1 
                    do ik=1,nsk 
                       numik = NumNodebyCell(k,ik)
                       if (numik.eq.numsf) then 
                          numsflocalcell = ik 
                          ii = 1
                       endif
                    enddo
                    if (ii.eq.-1) then 
                       write(*,*)' on ne trouve pas le noeud ',k,numsf
                       stop
                    endif

                    PoidsXCellCG(k,numsflocalcell) =
     &                 PoidsXCellCG(k,numsflocalcell)
     &                + VolT/4.d0/dfloat(nsf)
                    
                    
                 enddo
                 
c                 
c     vecteur normal a 12f sortant de K et de norme la surface de 12f  
c
                 call VecNormalT(x1,y1,z1,x2,y2,z2,xf,yf,zf,x,y,z,
     &                      vx,vy,vz,surf)

                 v12fx = vx*surf
                 v12fy = vy*surf
                 v12fz = vz*surf

                 xi = (x1+x2+xf)/3.d0
                 yi = (y1+y2+yf)/3.d0
                 zi = (z1+z2+zf)/3.d0

                 Xfi(jf,1) = Xfi(jf,1) + xi*surf
                 Xfi(jf,2) = Xfi(jf,2) + yi*surf
                 Xfi(jf,3) = Xfi(jf,3) + zi*surf 
                 si = si + surf


                 PoidsXFaceCG(numf,numn1localface) =
     &             PoidsXFaceCG(numf,numn1localface) + surf/3.d0
                 PoidsXFaceCG(numf,numn2localface) =
     &             PoidsXFaceCG(numf,numn2localface) + surf/3.d0

                 do is=1,nsf
                    PoidsXFaceCG(numf,is) =
     &                   PoidsXFaceCG(numf,is) + surf/3.d0/dfloat(nsf)                    
                 enddo

                 VecNormalKi(jf,1) = VecNormalKi(jf,1) + v12fx
                 VecNormalKi(jf,2) = VecNormalKi(jf,2) + v12fy
                 VecNormalKi(jf,3) = VecNormalKi(jf,3) + v12fz

              enddo ! fin boucle sur les aretes de la face jf 


              
              do is=1,nsf
                 PoidsXFaceCG(numf,is) =
     &                PoidsXFaceCG(numf,is)/si                   
              enddo              

              Surfi(jf) = dsqrt( VecNormalKi(jf,1)**2 + 
     &              VecNormalKi(jf,2)**2 + VecNormalKi(jf,3)**2 )

              VecNormalKi(jf,1) = VecNormalKi(jf,1)/Surfi(jf)
              VecNormalKi(jf,2) = VecNormalKi(jf,2)/Surfi(jf)
              VecNormalKi(jf,3) = VecNormalKi(jf,3)/Surfi(jf)


              VecNormalKSigma(k,jf,:) = VecNormalKi(jf,:)

              Xfi(jf,1) = Xfi(jf,1)/si
              Xfi(jf,2) = Xfi(jf,2)/si
              Xfi(jf,3) = Xfi(jf,3)/si 

c              write(*,*)' Surfi ',jf,Surfi(jf),si
c              write(*,*)' Xfi ',(Xfi(jf,m),m=1,NbDim)
c              write(*,*)' vec normal ',(VecNormalKi(jf,m),m=1,NbDim)              
c
c
c
c     centre de gravite et surface de la face 
c
              do m=1,NbDim
                 XFaceCG(numf,m) = Xfi(jf,m)
              enddo
        
               SurfaceFace(numf) = Surfi(jf)


c
c     test poids face 
c
           s = 0.d0 
           XCG(:) = 0.d0 
           do ik=1,nsf
              n = NumNodebyFace(numf,ik)
              XCG(:) = XCG(:) + PoidsXFaceCG(numf,ik)*XNode(n,:)
              s = s + PoidsXFaceCG(k,ik)
           enddo


c           write(*,*)' XFace CG ',XFaceCG(numf,:)
c           write(*,*)' XFace CG ',XCG(:)       
c           write(*,*)' Poids CG face ',PoidsXFaceCG(numf,1:nsf)
c           write(*,*)' somme poids face ',s                
               

           enddo ! fin boucle sur les faces 

           PoidsXCellCG(k,:) =
     &          PoidsXCellCG(k,:)/Volk 
           

           XCellCG(k,:) = XCellCG(k,:)/Volk  
c
c     test poids cell 
c
           s = 0.d0 
           XCG(:) = 0.d0 
           do ik=1,nsk
              n = NumNodebyCell(k,ik)
              XCG(:) = XCG(:) + PoidsXCellCG(k,ik)*XNode(n,:)
              s = s + PoidsXCellCG(k,ik)
           enddo

           
c           write(*,*)' XCell CG ',XCellCG(k,:)
c           write(*,*)' XCell CG ',XCG(:)
c           write(*,*)' XCell IS ',XCell(k,:)           
c           write(*,*)' VolCell ',Volk,VolCell(k) 
c           write(*,*)' Poids CG Cell ',PoidsXCellCG(k,1:nsk)
c           write(*,*)' somme poids cell ',s 
c
c          
cccccccccccccccccccccccccccccc
c
c
c     Fin boucle sur les mailles 
c
        enddo
c
ccccccccccccccccccccccccccccccccc
c
c
c     Stockage du vecteur normal par face 
c
        do i=1,NbFace

           k1 = NumCellbyFace(i,1)
           ii = -1
           do ik=1,NbFacebyCell(k1)
              n = NumFacebyCell(k1,ik)
              if (n.eq.i) then
                 iik = ik
                 ii = 1
              endif
           enddo
           if (ii.eq.-1) then
              write(*,*)' face non trouvee ',i,k1
              stop
           endif
           
           VecNormalbyFace(i,:) = VecNormalKSigma(k1,iik,:)
           
        enddo
c        
c
c        
ccccccccccccccccccccccccccccccccc
c
c     boucle sur les faces 
c
ccccccccccccccccccccccccccccccc
c
c        write(*,*)' NbFaceFrac ',NbFaceFrac
c
c
c
        do numf = 1,NbFace

c
c     naf : nb d'aretes par face 
c
           naf = NbAretebyFace(numf)


           xf = XFaceCG(numf,1)
           yf = XFaceCG(numf,2)
           zf = XFaceCG(numf,3)
c
c     boucle sur les aretes de la face frac
c

              do ia = 1,naf 
                 numa = NumAretebyFace(numf,ia)
c
c
cccccccccccccccccccccccccccc
c
c                noeuds de l'arete numa 
c
                 numn1 = NumNodebyArete(numa,1)
                 numn2 = NumNodebyArete(numa,2)

                 x1 = XNode(numn1,1)
                 y1 = XNode(numn1,2)
                 z1 = XNode(numn1,3)

                 x2 = XNode(numn2,1)
                 y2 = XNode(numn2,2)
                 z2 = XNode(numn2,3)

c
c     Vecteur normal a l'arete 12 ds le plan 12f 
c
c     n = 1-f + alpha*(2-1)  avec n.(2-1)=0  
c                 

                 s1f = (x1-xf)*(x2-x1) 
     &              + (y1-yf)*(y2-y1) + (z1-zf)*(z2-z1)
                 s12 = (x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2

                 alpha = -s1f/s12

                 vx = x1-xf + alpha*(x2-x1)
                 vy = y1-yf + alpha*(y2-y1)
                 vz = z1-zf + alpha*(z2-z1)

                 s = dsqrt(vx**2+vy**2+vz**2) ! pour normalisation du vecteur 

                 vx = vx/s
                 vy = vy/s
                 vz = vz/s

                 d12 = dsqrt(s12)


                 VecNormalSi(ia,1) = vx 
                 VecNormalSi(ia,2) = vy 
                 VecNormalSi(ia,3) = vz


                 VecNormalSigmaArete(numf,ia,:) = VecNormalSi(ia,:)

                 SizeArete(numa) = d12
                 
c
c     fin boucle sur les aretes de la face frac 
c
              enddo

c
c    ccccccccccccccccccccccccccccccc
c
c     fin boucle sur les faces 
c
        enddo

cccccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccccccc
       return
       end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccc Une subroutine qui calcul le repère orthonormé "Directe" pur chaque face du maillage ccccccccc
ccccccc pour chaque "iface" (le numéro globale de la face), on associe (Veci(iface,:),Vecj(iface,:),Vecnormalbyface(iface,:))
ccccccc un repère othonormé tels que Vecnormalbyface(iface,:) est le vecteur normal de la face "iface".
cccccc  
cccccc     

c     on calcule aussi le facteur OrientationbyFaceFrac qui donne une orientation intrinseque de la normale
c
c     VecNormalbyFace(iface,:)*OrientationbyFaceFrac(ifrac) 
c      
c     (VecNormalbyFace a une orientation locale) 
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      
      
      subroutine RepereOrthoByFace (
     &     NbFace,NbFaceFrac,VecNormalbyFace,
     &     Veci,Vecj,OrientationbyFaceFrac,
     &     NumFaceFracVersFace)
  
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     DECLARATIONS
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
ccccccccccccccccccccccccccccccccc  ENTREE ccccccccccccccccccccccccccccccccccccccc
      
      dimension VecNormalbyFace(NbFace,NbDim)

      dimension NumFaceFracVersFace(NbFaceFrac)
      
ccccccccccccccccccccccccccccccccc  SORTIES cccccccccccccccccccccccccccccccccccccccc
      
      dimension Veci(NbFace,NbDim)
      dimension Vecj(NbFace,NbDim)


      dimension OrientationbyFaceFrac(NbFaceFrac)

      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c     VARIABLES LOCALES 
c      
c     Définition de la base canonique

      dimension Vec_e1(NbDim)
      dimension Vec_e2(NbDim)    
      dimension Vec_e3(NbDim)

      dimension Vec_z(NbDim,NbDim)
      dimension Vecteur_z(NbDim)


      dimension Vec_max_norme(NbDim)
      
      dimension produit_vect(NbDim)
            
ccccccccccccccccccccccccccccccccccccccccccccc

      
      Vec_e1(:) = 0.d0
      Vec_e2(:) = 0.d0
      Vec_e3(:) = 0.d0
      
      Vec_e1(1) = 1.d0
      Vec_e2(2) = 1.d0
      Vec_e3(3) = 1.d0
      
      do iface=1,NbFace
         

         Vec_z(1,:) = Vec_e1(:) -
     &        prodscal(Vec_e1,VecNormalbyFace(iface,:))
     &        *VecNormalbyFace(iface,:)



         Vec_z(2,:) = Vec_e2(:) -
     &        prodscal(Vec_e2,VecNormalbyFace(iface,:))
     &        *VecNormalbyFace(iface,:)


         Vec_z(3,:) = Vec_e3(:) -
     &        prodscal(Vec_e3,VecNormalbyFace(iface,:))
     &        *VecNormalbyFace(iface,:)





         
         value_max = 0.d0

         in = 0
         do i=1,NbDim
            
            if(vecnorme(Vec_z(i,:)) .ge. value_max) then

               Vec_max_norme(:) = Vec_z(i,:)
               value_max = vecnorme(Vec_z(i,:))
               in = i 

            endif

         enddo

         do i=1,NbDim
            if (i.ne.in) then
               i1 = i 
            endif
         enddo
         do i=1,NbDim
            if ( (i.ne.in).and.(i.ne.i1) ) then
               i2 = i 
            endif
         enddo         
         
         Veci(iface,:) = (1.d0/vecnorme(Vec_max_norme))*Vec_max_norme(:)


         psi1 = prodscal(Veci(iface,:),vec_z(i1,:))
         psi2 = prodscal(Veci(iface,:),vec_z(i2,:))

         call ProdVectoriel(Vec_z(i1,:),Veci(iface,:),
     &        produit_vect)

         pvec1 = vecnorme(produit_vect)


         call ProdVectoriel(Vec_z(i2,:),Veci(iface,:),
     &        produit_vect)

         pvec2 = vecnorme(produit_vect)         

         if (pvec1.gt.pvec2) then
            
            ip = i1
            psip = psi1
            
         else

            ip = i2
            psip = psi2     
            
         endif
c
c     vecj = alpha Veci(iface,:) + vec_z(ip,:) on calcule alpha tel que (vecj,veci) = 0
c
c     alpha = - psip 
c         
         alpha = - psip 

         Vecj(iface,:) = alpha*Veci(iface,:) + vec_z(ip,:)

         rnorme = vecnorme(Vecj(iface,:))

         Vecj(iface,:) = Vecj(iface,:)/rnorme
c
c     le produit vectoriel donne pour Vecj une orientation qui depend de l'orientation de VecNormalbyFace donc non fixe !!!!
c         
c         call ProdVectoriel(VecNormalbyface(iface,:),Veci(iface,:),
c     &                           produit_vect)
         
c         Vecj(iface,:) = produit_vect(:)
c
c         
        
      enddo
c
c
c      
c     calcul de l'orientation de la normale aux faces frac 
c 
      do ifrac=1,NbFaceFrac

         iface = NumFaceFracVersFace(ifrac)

         call ProdVectoriel(Veci(iface,:),Vecj(iface,:),
     &        produit_vect)         


         OrientationbyFaceFrac(ifrac) =
     &       prodscal(produit_vect,VecNormalbyFace(iface,:))
         
      enddo

      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
ccccccccccccccccc Subroutine qui calcul la norme L2 du vecteur Vect ccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      
      function vecnorme(Vect)



      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
c
      include 'include_parameter'
c
      

      dimension Vect(NbDim)
      
      
      vecnorme = (Vect(1)**2 + Vect(2)**2 + Vect(3)**2)
     &     **(0.5d0)

      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccc Subroutine qui calcul le produit vectoriel, vect1 vectoriel vect2 ccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  

      subroutine ProdVectoriel(Vect1,Vect2,Resultat)



c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
c
      include 'include_parameter'
c
c

cccccccccccccccccccccccccccccccccccccc  ENTRE ccccccccccccccccccccccccccccccccccccccc
      
      dimension Vect1(NbDim)
      dimension Vect2(NbDim)
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
cccccccccccccccccccccccccccccccccc  SORTIE cccccccccccccccccccccccccccccccccccccccccc

      dimension Resultat(NbDim)
      
      Resultat(1) = Vect1(2)*Vect2(3) - Vect2(2)*Vect1(3)
      
      Resultat(2) = -(Vect1(1)*Vect2(3) - Vect2(1)*Vect1(3))
      
      Resultat(3) = Vect1(1)*Vect2(2) - Vect2(1)*Vect1(2)
      
      
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     diametre des mailles
c
c     Warning: TMP 2D xz !!!!!!!!!!!!!!!!!!
c
ccccccccccccccccccccccccccccccccccccccc
c      
      subroutine ComputeDiambyCell(NbCell,NbNode,
     &             XNode,NbNodebyCell,
     &             NumNodebyCell,
     &             diamK)
c
c
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c     
c     
      include 'include_parameter'
      
c
ccccc ENTRE     

c     coordonnees des noeuds XS 
c
      dimension XNode(NbNode,NbDim)

c     nb de noeuds S par maille K 
c
      dimension NbNodebyCell(NbCell)
c     
c     num des noeuds S de chaque maille K dans l'ordre cyclique 
c
      dimension NumNodebyCell(NbCell,NbNodebyCellMax)
      
ccccc SORTI

      dimension diamK(NbCell)

ccccc WORKSPACE

      dimension vecij(NbDim)
      
      
ccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccc



      do ik = 1,NbCell

         nbnk = NbNodebyCell(ik)

         dist_max = 0.d0

         do i = 1,nbnk
            
            do j = 1,nbnk


               numnodei = NumNodebyCell(ik,i)
               numnodej = NumNodebyCell(ik,j)                        
               
               vecij(:) = XNode(numnodei,:)
     &              -XNode(numnodej,:)
 
c               distij = vecnorme(vecij)  ! TMP 2D xz 

               distij = dsqrt(vecij(1)**2 + vecij(3)**2) 
               
               if(distij .ge. dist_max) then
                  dist_max = distij 
               endif


            enddo

         enddo

         diamK(ik) = dist_max
         
      enddo
         


      return
      end


