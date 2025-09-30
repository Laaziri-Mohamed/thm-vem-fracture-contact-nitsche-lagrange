
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                                                     
c
c     
c         LECTURE DU FORMAT Benchmark 3D 
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
      subroutine lectureBench3D(num,iflabelfrac,iflecturelitho,
     &             NbCell,NbNode,NbFace,NbArete,
     &             XNode,NbNodebyCell,NumNodebyCell,
     &             NbFacebyCell,NumFacebyCell,
     &             NbAretebyFace,NumAretebyFace,
     &             NbNodebyFace,NumNodebyFace,
     &             NumCellbyFace,NumNodebyArete,
     &             IndRockTypebyCell,LabelbyFace)
      

c
c
	implicit double precision (a-h,o-z)
	implicit integer (i-n)
c
c
        include 'include_parameter'
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     LECTURE DU FORMAT BENCHMARK 3D 
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
c
c
c     Label = 1 si face frac 
c
        dimension LabelbyFace(NbFace)
c
c     rock type by cell
c
        dimension IndRockTypebyCell(NbCell)
c
c
c        
ccccccccccccccccccccccccccccccccccccc
c
c
c

        read(num,*)
        do n = 1,NbNode
           read(num,*)(XNode(n,i),i=1,NbDim)
c           write(*,*)(XNode(n,i),i=1,NbDim)
        enddo


        read(num,*)
        do k = 1,NbCell
           read(num,*)NbFacebyCell(k),
     &               (NumFacebyCell(k,i),i=1,NbFacebyCell(k))
           if (NbFacebyCell(k).gt.NbFacebyCellMax) then 
           write(*,*)' NbFacebyCell > Max ',
     &                 k,NbFacebyCell(k),NbFacebyCellMax
           stop
           endif
        enddo


        read(num,*)
        do k = 1,NbCell
           read(num,*)NbNodebyCell(k),
     &                (NumNodebyCell(k,i),i=1,NbNodebyCell(k))
           if (NbNodebyCell(k).gt.NbNodebyCellMax) then 
           write(*,*)' NbNodebyCell > Max ',
     &                 k,NbNodebyCell(k),NbNodebyCellMax
           stop
           endif
        enddo


        read(num,*)
        do k = 1,NbFace
        read(num,*)NbAretebyFace(k),
     &             (NumAretebyFace(k,i),i=1,NbAretebyFace(k))
           if (NbAretebyFace(k).gt.NbAretebyFaceMax) then 
           write(*,*)' NbAretebyFace > Max ',k,
     &                    NbAretebyFace(k),NbAretebyFaceMax
           stop
           endif
        enddo


        read(num,*)
        do k = 1,NbFace
           read(num,*)NbNodebyFace(k),
     &                (NumNodebyFace(k,i),i=1,NbNodebyFace(k))
           if (NbNodebyFace(k).gt.NbNodebyFaceMax) then 
           write(*,*)' NbNodebyFace > Max ',
     &                 k,NbNodebyFace(k),NbNodebyFaceMax
           stop
           endif
        enddo


        read(num,*)
        do k = 1,NbFace
           read(num,*)NumCellbyFace(k,1),NumCellbyFace(k,2)
        enddo


        read(num,*)
        do k = 1,NbArete     
           read(num,*)NumNodebyArete(k,1),NumNodebyArete(k,2)
        enddo

        if (iflecturelitho.eq.1) then 
           read(num,*)
           do k=1,NbCell
c              LITHO PAS MAILLE 
              read(num,*)iiii,IndRockTypebyCell(k)
           enddo
        endif
        
        if (iflabelfrac.eq.1) then 
           read(num,*)
           do k = 1,NbFace
              read(num,*)LabelbyFace(k)
           enddo    
        endif


        close(num)

        iaffich = 0
c        iaffich = 1
        if (iaffich.eq.1) then 
        write(*,*)
        write(*,*)' NbNode ',NbNode 
        write(*,*)' NbCell ',NbCell
        write(*,*)' NbFace ',NbFace
        write(*,*)' NbArete ',NbArete
        write(*,*)
        do n = 1,NbNode
           write(*,*)'Node n = ',n,(XNode(n,i),i=1,NbDim)
        enddo
        write(*,*)
        write(*,*)' volumes par les faces '
        do k = 1,NbCell
           write(*,*)'maille ',k,
     &        NbFacebyCell(k),(NumFacebyCell(k,i),i=1,NbFacebyCell(k))
        enddo
        write(*,*)
        write(*,*)' volumes par les noeuds '
        do k = 1,NbCell
           write(*,*)'maille ',k,
     &            NbNodebyCell(k),
     &            (NumNodebyCell(k,i),i=1,NbNodebyCell(k))
        enddo
        write(*,*)
        write(*,*)' faces par les aretes '
        do k = 1,NbFace
         write(*,*)NbAretebyFace(k),
     &             (NumAretebyFace(k,i),i=1,NbAretebyFace(k))
        enddo
        write(*,*)
        write(*,*)' faces par les noeuds '
        do k = 1,NbFace
           write(*,*)NbNodebyFace(k),
     &                (NumNodebyFace(k,i),i=1,NbNodebyFace(k))
        enddo
        write(*,*)
        write(*,*)' volumes adjacents aux faces '
        do k = 1,NbFace
           write(*,*)NumCellbyFace(k,1),NumCellbyFace(k,2)
        enddo
        write(*,*)
        write(*,*)' Aretes par les noeuds '
        do k = 1,NbArete
           write(*,*)NumNodebyArete(k,1),NumNodebyArete(k,2)
        enddo
        write(*,*)

        endif
c
c
c
ccccccccccccccccccccccccccccccccccccc

        return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                                                                                                       
c                         
c          Schema VAG 3D SUR MAILLAGE POLYHEDRIQUE 
c
c            NbDim = 3 (include_parameter)
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
      subroutine vag3D(NbCell,NbNode,NbFace,NbArete,
     &             NbFaceFrac,NbCV, 
     &             XNode,NbNodebyCell,NumNodebyCell,
     &             NbFacebyCell,NumFacebyCell,
     &             NbAretebyFace,NumAretebyFace,
     &             NbNodebyFace,NumNodebyFace,
     &             NumCellbyFace,NumNodebyArete,
     &             NumFaceFracVersFace,NumFaceVersFaceFrac,IndFaceFrac,
     &             NumIncCell,NumIncFaceFrac,NumIncNode,
     &             perm,permfdf,epaisseurf, 
     &             XCell,XFace,XFaceFrac,XInc,
     &             VolCell,VolFaceFrac,
     &             NbInterfacebyCell,NumInterfacebyCell,
     &             NbInterfacebyFaceFrac,NumInterfacebyFaceFrac,
     &             TransCell,TransFaceFrac,GCell,GFaceFrac)
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
c
ccccccccccccccccccccccccc
c
c     Identification des faces fractures parmi l'ensemble des faces 
c
c       num face frac > num face 
        dimension NumFaceFracVersFace(NbFaceFrac)
c
c       num face  > num face frac
        dimension NumFaceVersFaceFrac(NbFace)
c
c
c       num face > ind 1 si face frac 0 sinon 
        dimension IndFaceFrac(NbFace)
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
c     pour VAG (tous les noeuds) 
c     num node  > num inc  
c        
        dimension NumIncNode(NbNode)
c
c
ccccccccccccccccccccccccc
c
c     Permeabilite matrice par maille 
c     num cell  >  perm
c
        dimension perm(NbCell,NbDim,NbDim)
c
c     Permeabilite longitudinale fracture par face fracture (permf * epaisseur df )
c     supposee isotrope
c     num face frac > permfdf 
c
        dimension permfdf(NbFaceFrac) 
c
c     epaisseur fracture
c
        dimension epaisseurf(NbFaceFrac)

cccccccccccccccccccccccc



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     SORTIES 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     coordonnee du centre de maille XK (maille supposee etoilee par rapport a XK)
c     num cell  >  XCell
c
        dimension XCell(NbCell,NbDim)
c
c     coordonnee du centre de Face Frac
c     num face frac  >  XFaceFrac 
c
        dimension XFaceFrac(NbFaceFrac,NbDim)
c
c     coordonnee du centre de Face     
c
        dimension XFace(NbFace,NbDim)
c
c
c     coordonnee des inconnues (XK, Xsig, Xs, Xa, ...) 
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
cccccccccccccccccccccccc
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
     &         NbInterfacebyFaceFracMax,NbInterfacebyFaceFracMax)
c
c
c     Gradient par maille : g = \sum_i\in IK (uk-ui)*GCell(k,i) 
c
        dimension GCell(NbCell,NbInterfacebyCellMax,NbDim)
c
c     Gradient par face frac : g = \sum_i\in IS (uk-ui)*GFaceFrac(k,i) 
c
        dimension GFaceFrac(NbFaceFrac,
     &         NbInterfacebyFaceFracMax,NbDim)
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     WorkSpaces locaux 
c
ccccccccccccccccccccccccccccccccccccccccccc
c
c
c     Par maille: gradients par tetra xK,xsigma,arete de la face 
c
c     les inconnues aux faces frac sont conservées 
c     les inconnues aux faces non frac sont interpolées 
c
c     grad_T_iface,iarete = \sum_(n = noeuds maille+face frac) GKij(iface,iarete,n)(un-uk)
c
c     le num du noeud est celui dans la numerotation locale des noeuds 
c     de la maille suivi des numeros de face frac 
c
        dimension GKij(NbFacebyCellMax,NbAretebyFaceMax,
     &                 NbNodebyCellMax+NbFacebyCellMax,NbDim)
c
c     Par maille: volume par tetra T_iface,iarete 
c
        dimension VolT(NbFacebyCellMax,NbAretebyFaceMax)
c
c     Gradient tangentiel fct affine sur le triangle plonge en 3D 
c     Grad_sig_iarete = \sum_{noeuds ns de la face } GTtgt(iarete,ns)(uf-us)
c
        dimension GTtgt(NbAretebyFaceMax,NbNodebyFaceMax,NbDim)
c
c     Surface du triangle (arete,centre de face)
c
        dimension SurfT(NbAretebyFaceMax)
c
c     Inversion matrice 3-3
c
        dimension AA(3,3)
        dimension BB(3,3)
        dimension IPVT(4)

c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     centre des mailles = isobarycentre des nodes 
c
        do k = 1,NbCell
           do m = 1,NbDim
              XCell(k,m) = 0.d0
           enddo

           do nk = 1,NbNodebyCell(k)
              n = NumNodebyCell(k,nk)
              do m = 1,NbDim
                 XCell(k,m) = XCell(k,m) + XNode(n,m)
              enddo
           enddo

           do m = 1,NbDim
              XCell(k,m) = XCell(k,m)/dfloat(NbNodebyCell(k))
              XInc(NumIncCell(k),m) = XCell(k,m)
           enddo
c           write(*,*)' Centre maille ',k,(XCell(k,m),m=1,NbDim)



        enddo
c
c
c
c
c     Centre des faces  = isobarycentre des nodes 
c
        do k = 1,NbFace
           do m = 1,NbDim
              XFace(k,m) = 0.d0
           enddo

           do nk = 1,NbNodebyFace(k)
              n = NumNodebyFace(k,nk)
              do m = 1,NbDim
                 XFace(k,m) = XFace(k,m) + XNode(n,m)
              enddo
           enddo

           do m = 1,NbDim
              XFace(k,m) = XFace(k,m)/dfloat(NbNodebyFace(k))
           enddo

           if (IndFaceFrac(k).eq.1) then 
              kfrac = NumFaceVersFaceFrac(k)
              inc = NumIncFaceFrac(kfrac)
              do m = 1,NbDim
                 XFaceFrac(kfrac,m) = XFace(k,m)
                 XInc(inc,m) = XFace(k,m)
              enddo
           endif

c           write(*,*)' Centre Face ',k,(XFace(k,m),m=1,NbDim)
        enddo
c
c     XInc Nodes
c
        do i=1,NbNode
           do m = 1,NbDim
              XInc(NumIncNode(i),m) = XNode(i,m)
           enddo           
        enddo
        
c        write(*,*)
c        do i=1,NbCV
c           write(*,*)' XInc ',i,(XInc(i,m),m=1,NbDim)
c        enddo

c        stop
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
c     nsk: nb de noeuds de la maille 
c     nfk: nb de faces de la maille 
c     nffk: nb de faces frac de la maille 
c     nik = nsk + nffk 
c
           nsk = NbNodebyCell(k)
           nfk = NbFacebyCell(k)
c
c     nb de faces frac de la maille 
c
           nffk = 0
           do j = 1,nfk 
              numf = NumFacebyCell(k,j)
              if (IndFaceFrac(numf).eq.1) then 
                 nffk = nffk + 1
              endif
           enddo
c
c     NbInterfacebyCell
c
           nik = nsk + nffk 
           NbInterfacebyCell(k) = nik 
c
c     NumInterfacebyCell
c

           do in = 1,nsk 
              ns = NumNodebyCell(k,in)
              numinc = NumIncNode(ns)
              NumInterfacebyCell(k,in) = numinc
           enddo

           iffk = 0
           do j = 1,nfk
              numf = NumFacebyCell(k,j)
              if (IndFaceFrac(numf).eq.1) then 
                 iffk = iffk + 1
                 numfrac = NumFaceVersFaceFrac(numf)
                 numinc = NumIncFaceFrac(numfrac)
                 NumInterfacebyCell(k,iffk + nsk) = numinc
              endif
           enddo

c           write(*,*)
c           write(*,*)' IK ',k,nik,(NumInterfacebyCell(k,i),i=1,nik)

           

cccccccccccccccccccccccccccc

           xk = XCell(k,1)
           yk = XCell(k,2)
           zk = XCell(k,3)

           VolCell(k) = 0.d0
c
c
ccccccccccccccccccccccccccccc
c
c          boucle sur les faces de la maille k
c
           do jf = 1,NbFacebyCell(k)
              numf = NumFacebyCell(k,jf)
              xf = XFace(numf,1)
              yf = XFace(numf,2)
              zf = XFace(numf,3)

              if (IndFaceFrac(numf).eq.1) then 
                 numfrac = NumFaceVersFaceFrac(numf)
c
c                numero local de la face frac = numf, numfrac 
c                 
                 iffk = 0
                 numfraclocalcell = 0
                 do j = 1,nfk
                    numf1 = NumFacebyCell(k,j)
                    if (IndFaceFrac(numf1).eq.1) then 
                       iffk = iffk + 1
                       if (numf1.eq.numf) then 
                          numfraclocalcell = iffk + nsk 
                       endif
                    endif
                 enddo 
                 if (numfraclocalcell.eq.0) then 
                    write(*,*)'pb num local face frac ',numf,numfrac
                    stop
                 endif

c                 write(*,*)' numfraclocalcell ',numfraclocalcell


              endif 

c
c             boucle sur les aretes de la face numf 
c
              do ia = 1,NbAretebyFace(numf)
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
                 VolT(jf,ia) = VolTetra(x1,y1,z1,x2,y2,z2,
     &                                  xf,yf,zf,xk,yk,zk)
c
c     Volume de la maille = somme du volume des tetras 
c
                 VolCell(k) = VolCell(k) + VolT(jf,ia)


c                 write(*,*)' vol T ',k,jf,ia,VolT(jf,ia)
c
c     centre du tetra X             
c
                 x = (x1+x2+xk+xf)/4.d0
                 y = (y1+y2+yk+yf)/4.d0
                 z = (z1+z2+zk+zf)/4.d0
c                 
c     vecteur normal a 12k pointant vers X et de norme la surface de 12k  
c
                 call VecNormalT(x1,y1,z1,x2,y2,z2,xk,yk,zk,x,y,z,
     &                      vx,vy,vz,surf)

                 v12kx = vx*surf
                 v12ky = vy*surf
                 v12kz = vz*surf

c                 write(*,*)' v12k ',v12kx,v12kz,v12kz
c                 
c     vecteur normal a 1fk pointant vers X et de norme la surface de 1fk  
c
                 call VecNormalT(x1,y1,z1,xf,yf,zf,xk,yk,zk,x,y,z,
     &                      vx,vy,vz,surf)

                 v1fkx = vx*surf
                 v1fky = vy*surf
                 v1fkz = vz*surf

c                 write(*,*)' v1fk ',v1fkx,v1fkz,v1fkz
c                 
c     vecteur normal a f2k pointant vers X et de norme la surface de f2k  
c
                 call VecNormalT(xf,yf,zf,x2,y2,z2,xk,yk,zk,x,y,z,
     &                      vx,vy,vz,surf)

                 vf2kx = vx*surf
                 vf2ky = vy*surf
                 vf2kz = vz*surf

c                 write(*,*)' vf2k ',vf2kx,vf2kz,vf2kz
c
c ?? GKij ??

c     Vecteurs GKij : grad_T_jf,ia = \sum_(n=noeuds maille+faces frac) GKij(jf,ia,n)(uk-un)
c
c     le num du noeud est celui dans la numerotation locale des noeuds 
c     de la maille !!! 
c
c     Grad T = 1/(3*VolT) V12k (uk-uf) + V1fk (uk-u2) + Vf2k (uk-u1) 
c
c     avec uf = 1/NbNodebyFace /sum_(n= noeuds de la face) un   si face non frac 
c
c          uf si face frac  
c
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
c                 write(*,*)' numn2localcell ',numn2localcell
c
c                numero local de la face si face fracture 
c
                 
c
c                init a zero de GKij                  
                 do j=1,nik 
                    do n=1,NbDim
                       GKij(jf,ia,j,n) = 0.d0
                    enddo
                 enddo

                 if (IndFaceFrac(numf).eq.0) then 
c
c     numf face non frac: interpolation isobarycentrique 
c
c
c     boucle sur les noeuds de la face pour le terme 1/(3*VolT) V12k (uk-uf)
c
                 ss = 1.d0/(3.d0*VolT(jf,ia)*dfloat(NbNodebyFace(numf)))

                 do j=1,NbNodebyFace(numf)
                    numj = NumNodebyFace(numf,j)
c
c     numero du noeud numj dans la numerotation locale des noeuds de la maille k 
c
                    ii = -1 
                    do ik=1,NbNodebyCell(k)
                      numik = NumNodebyCell(k,ik)
                      if (numik.eq.numj) then 
                         numjlocalcell = ik 
                         ii = 1
                      endif
                    enddo
                    if (ii.eq.-1) then 
                    write(*,*)' on ne trouve pas le noeud ',k,numj
                    stop
                    endif
c                    write(*,*)' numjlocalcell ',j,numjlocalcell

                    GKij(jf,ia,numjlocalcell,1) = v12kx*ss
                    GKij(jf,ia,numjlocalcell,2) = v12ky*ss
                    GKij(jf,ia,numjlocalcell,3) = v12kz*ss
                    
                 enddo

                 ss = 1.d0/(3.d0*VolT(jf,ia))
          
                 GKij(jf,ia,numn1localcell,1) = 
     &                           GKij(jf,ia,numn1localcell,1) + vf2kx*ss
                 GKij(jf,ia,numn1localcell,2) = 
     &                           GKij(jf,ia,numn1localcell,2) + vf2ky*ss
                 GKij(jf,ia,numn1localcell,3) = 
     &                           GKij(jf,ia,numn1localcell,3) + vf2kz*ss


                 GKij(jf,ia,numn2localcell,1) = 
     &                           GKij(jf,ia,numn2localcell,1) + v1fkx*ss
                 GKij(jf,ia,numn2localcell,2) = 
     &                           GKij(jf,ia,numn2localcell,2) + v1fky*ss
                 GKij(jf,ia,numn2localcell,3) = 
     &                           GKij(jf,ia,numn2localcell,3) + v1fkz*ss

                 else 
c
c     numf face fracture numfrac , numfraclocalcell 
c
                 ss = 1.d0/(3.d0*VolT(jf,ia))

                 GKij(jf,ia,numfraclocalcell,1) = v12kx*ss
                 GKij(jf,ia,numfraclocalcell,2) = v12ky*ss
                 GKij(jf,ia,numfraclocalcell,3) = v12kz*ss
          
                 GKij(jf,ia,numn1localcell,1) =  vf2kx*ss
                 GKij(jf,ia,numn1localcell,2) =  vf2ky*ss
                 GKij(jf,ia,numn1localcell,3) =  vf2kz*ss

                 GKij(jf,ia,numn2localcell,1) =  v1fkx*ss
                 GKij(jf,ia,numn2localcell,2) =  v1fky*ss
                 GKij(jf,ia,numn2localcell,3) =  v1fkz*ss


                 endif 

c              write(*,*)   
c              do j=1,nik 
c                 write(*,*)' GKij ',k,jf,ia,j,GKij(jf,ia,j,1),
c     &                              GKij(jf,ia,j,2),GKij(jf,ia,j,3)
c              enddo
c              write(*,*)
c
cccccccccccccccccccccccccccc
c
c     test du gradient GKij 
c
              uk = xk + 2*yk + 3*zk + 1
              sx = 1
              sy = 2
              sz = 3
              do j=1,nik 
                 nsj =  NumInterfacebyCell(k,j)
                 xj = XInc(nsj,1)
                 yj = XInc(nsj,2)
                 zj = XInc(nsj,3)
c                 write(*,*)nik,nsj,xj,yj,zj
                 uj = xj + 2*yj + 3*zj + 1
                 sx = sx - GKij(jf,ia,j,1)*(uk-uj)
                 sy = sy - GKij(jf,ia,j,2)*(uk-uj)
                 sz = sz - GKij(jf,ia,j,3)*(uk-uj)
              enddo

              if (dabs(sx)+dabs(sy)+dabs(sz).gt.1.0E-10) then 
                 write(*,*) 
                 write(*,*)' test grad ',k,jf,ia,sx,sy,sz
                 write(*,*) 
                 stop
              endif

cccccccccccccccccccccccccccc
c
c
c
ccccccccccccccccccccccccccccc
c
c     fin boucle aretes de la face numf 
c
              enddo
c
c     fin boucle face de la maille k 
c
           enddo
cccccccccccccccccccccccccccccc
c          
c     Calcul des TKij = \sum_jf,ia VolT(jf,ia) (\Perm_K GKij(jf,ia,i,.) GKij(jf,ia,j,.)  ) 
c
c     i,j=1,...,nik 
c

           do i=1,nik
              inci = NumInterfacebyCell(k,i)
           do j=1,nik
              incj = NumInterfacebyCell(k,j)
              s = 0.d0
              do jf = 1,NbFacebyCell(k)
              numf = NumFacebyCell(k,jf)
              do ia = 1,NbAretebyFace(numf)            
                 do m=1,NbDim
                 do n=1,NbDim
                   s = s + VolT(jf,ia)*perm(k,m,n)
     &                     *GKij(jf,ia,i,m)*GKij(jf,ia,j,n)
                 enddo
                 enddo
              enddo
              enddo

              TransCell(k,i,j) = s 

c            write(*,*)' TransCell ',k,i,j,inci,incj,TransCell(k,i,j)
           enddo
           enddo
c           write(*,*)
c
          
cccccccccccccccccccccccccccccc


          
cccccccccccccccccccccccccccccc
c
c     calcul du gradient par maille : moyenne au pro rata des volumes des volT 
c
           do ni = 1,nik
              do m=1,NbDim
                 GCell(k,ni,m) = 0.d0
              enddo
           enddo
c
           s = 0.d0
           do jf = 1,NbFacebyCell(k)
           numf = NumFacebyCell(k,jf)
           do ia = 1,NbAretebyFace(numf)   
              s = s + VolT(jf,ia)
              do m = 1,NbDim
              do nj = 1,nik
               GCell(k,nj,m) = GCell(k,nj,m) 
     &                + GKij(jf,ia,nj,m)*VolT(jf,ia)
              enddo     
              enddo
           enddo   
           enddo
c      
           do m = 1,NbDim
              do nj = 1,nik
                 GCell(k,nj,m) = GCell(k,nj,m)/s
              enddo     
           enddo
c

c         if (NbInterfacebyCell(k).lt.0) then 
c            write(*,*)' cell k vag ',k,NbInterfacebyCell(k)
c            stop
c          endif


cccccccccccccccccccccccccccccc
c
c     Fin boucle sur les mailles 
c
        enddo
c
cccccccccccccccccccccccccc
c
c
c     boucle sur les face frac 
c
c        write(*,*)' NbFaceFrac ',NbFaceFrac
c
c
c
        do ifrac = 1,NbFaceFrac
           numf = NumFaceFracVersFace(ifrac)
           inc = NumIncFaceFrac(ifrac)
c
c     nsf : nb de noeuds de la face = nb d'aretes 
c
           nsf = NbNodebyFace(numf)


           xf = XFace(numf,1)
           yf = XFace(numf,2)
           zf = XFace(numf,3)

c           xfrac = XFaceFrac(ifrac,1)
c           yfrac = XFaceFrac(ifrac,2)
c           zfrac = XFaceFrac(ifrac,3)
c           xinco = XInc(inc,1)
c           yinco = XInc(inc,2)
c           zinco = XInc(inc,3)
c           write(*,*)' xf = xfrac = xinc ',xf,xfrac,xinco 
c           write(*,*)' yf = yfrac = yinc ',yf,yfrac,yinco 
c           write(*,*)' zf = zfrac = zinc ',zf,zfrac,zinco

c
c     NbInterfacebyFaceFrac
c  

           NbInterfacebyFaceFrac(ifrac) = nsf
c
c     NumInterfacebyFaceFrac: dans l'ordre des noeuds de la face 
c 
           do is = 1,nsf
              nums = NumNodebyFace(numf,is)
              numinc = NumIncNode(nums)
              NumInterfacebyFaceFrac(ifrac,is) = numinc 
           enddo

c           write(*,*)
c           write(*,*)' Ifrac ',ifrac,nsf,
c     &       (NumInterfacebyFaceFrac(ifrac,is),is=1,nsf)


           volf = 0.d0 

c
c     boucle sur les aretes de la face frac
c

              do ia = 1,nsf 
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
c     numero local des deux noeuds dans la liste des noeuds de la face 
c
                 ii = -1 
                 do is=1,nsf
                    numis = NumNodebyFace(numf,is)
                    if (numis.eq.numn1) then 
                       numn1localcell = is 
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
                       numn2localcell = is 
                       ii = 1
                    endif
                 enddo
                 if (ii.eq.-1) then 
                    write(*,*)' on ne trouve pas le noeud ',numf,numn2
                    stop
                 endif
c
c     Surface du triangle (arete ia, centre de face) et vecteur normal unitaire 
c
                 x = (x1+x2+xf)/3.d0
                 y = (y1+y2+yf)/3.d0
                 z = (z1+z2+zf)/3.d0
                 call VecNormalT(x1,y1,z1,x2,y2,z2,xf,yf,zf,x,y,z,
     &                      vx,vy,vz,surf)

                 SurfT(ia) = surf

                 volf = volf + surf 
      
c                 write(*,*)' normale T ',vx,vy,vz
c                 write(*,*)' xf -x1 ',xf-x1,yf-y1,zf-z1
c                 write(*,*)' xf -x2 ',xf-x2,yf-y2,zf-z2

c
c     Gradient tangentiel 
c
c     Grad_sig_iarete = \sum_{noeuds ns de la face } GTtgt(ia,ns)(uf-us)
c
                 AA(1,1) = -x1+xf
                 AA(1,2) = -y1+yf
                 AA(1,3) = -z1+zf

                 AA(2,1) = -x2+xf
                 AA(2,2) = -y2+yf
                 AA(2,3) = -z2+zf

                 AA(3,1) = vx
                 AA(3,2) = vy                 
                 AA(3,3) = vz                

                 BB(1,1) = 1.d0
                 BB(1,2) = 0.d0
                 BB(1,3) = 0.d0

                 BB(2,1) = 0.d0
                 BB(2,2) = 1.d0
                 BB(2,3) = 0.d0

                 BB(3,1) = 0.d0
                 BB(3,2) = 0.d0
                 BB(3,3) = 1.d0

                 nsm = 3
                 nn  = 3
                 mm1 = 3
                 mm2 = mm1
                 info = 0
c                 call DGESV(nn,nsm,AA,mm1,IPVT,BB,mm2,info)
c                write(*,*)' info ',info
c                 if (info.ne.0) then 
c                    write(*,*)' pb DGESV AA ',info
c                    stop
c                 endif

                 call InvA(AA,BB)
c                 stop

                 do is = 1,nsf
                 do m = 1,NbDim
                    GTtgt(ia,is,m) = 0.d0 
                 enddo
                 enddo

                 GTtgt(ia,numn1localcell,1) = BB(1,1)
                 GTtgt(ia,numn1localcell,2) = BB(2,1)
                 GTtgt(ia,numn1localcell,3) = BB(3,1)

                 GTtgt(ia,numn2localcell,1) = BB(1,2)
                 GTtgt(ia,numn2localcell,2) = BB(2,2)
                 GTtgt(ia,numn2localcell,3) = BB(3,2)


                 s1 = BB(1,1)*vx + BB(2,1)*vy + BB(3,1)*vz
                 s2 = BB(1,2)*vx + BB(2,2)*vy + BB(3,2)*vz

c                 write(*,*)' check ps ',s1,s2
c
c     check gradient tangentiel 
c
                 u1 = 3*x1 + 2*y1 + z1 - 1
                 u2 = 3*x2 + 2*y2 + z2 - 1
                 uf = 3*xf + 2*yf + zf - 1
                 gx = 0.d0
                 gy = 0.d0
                 gz = 0.d0
                 do is=1,nsf
                    inc = NumInterfacebyFaceFrac(ifrac,is)
                    xi = XInc(inc,1)
                    yi = XInc(inc,2)
                    zi = XInc(inc,3)
                    ui = 3*xi + 2*yi + zi -1
                    gx = gx + GTtgt(ia,is,1)*(uf-ui)
                    gy = gy + GTtgt(ia,is,2)*(uf-ui)
                    gz = gz + GTtgt(ia,is,3)*(uf-ui) 
                 enddo
c               write(*,*)' grad tan ',gx,gy,gz
c               write(*,*)' g * n ',gx*vx + gy*vy + gz*vz
c               write(*,*)' g*(xf-x1)-(uf-u1) ',gx*(xf-x1)+gy*(yf-y1)
c     &                    + gz*(zf-z1)-(uf-u1)
c               write(*,*)' g*(xf-x2)-(uf-u2) ',gx*(xf-x2)+gy*(yf-y2) 
c     &                    + gz*(zf-z2)-(uf-u2)
c              write(*,*)         
c
c     fin boucle sur les aretes de la face frac 
c
              enddo
c
c
c
              VolFaceFrac(ifrac) = volf*epaisseurf(ifrac)

c              write(*,*)' vol face frac ',ifrac,VolFaceFrac(ifrac)


c          
c     Calcul des TSigij = \sum_ia SurfT(ia) (\Permf_sig GTtgt(ia,i,.) GTtgt(ia,j,.)  ) 
c
c     i,j=1,...,nsf 
c

           do i=1,nsf
           do j=1,nsf
              s = 0.d0
              do ia = 1,nsf          
                 do m=1,NbDim
                   s = s + SurfT(ia)*permfdf(ifrac)
     &                     *GTtgt(ia,i,m)*GTtgt(ia,j,m)
                 enddo
              enddo

              TransFaceFrac(ifrac,i,j) = s 

c            write(*,*)' TransFaceFrac ',ifrac,i,j,
c     &                    TransFaceFrac(ifrac,i,j)
           enddo
           enddo
c           write(*,*)
cccccccccccccccccccccccccccccc
c
c     calcul du gradient par face frac : 
c         moyenne au pro rata des SurfT 
c
           do ni = 1,nsf
              do m=1,NbDim
                 GFaceFrac(ifrac,ni,m) = 0.d0
              enddo
           enddo
c
           s = 0.d0
           do ia = 1,nsf 
              s = s + SurfT(ia)
              do m = 1,NbDim
              do nj = 1,nsf
               GFaceFrac(ifrac,nj,m) = GFaceFrac(ifrac,nj,m) 
     &            + GTtgt(ia,nj,m)*SurfT(ia)
              enddo     
              enddo
           enddo   
c      
           do m = 1,NbDim
              do nj = 1,nsf
                 GFaceFrac(ifrac,nj,m) = GFaceFrac(ifrac,nj,m)/s
              enddo     
           enddo
c
cccccccccccccccccccccccccccccc
c
c
c    
c
c     fin boucle sur les faces frac 
c
        enddo

cccccccccccccccccccccccccccc
       return
       end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     calcul de l'aire d'un triangle 
c
c
      subroutine aireT(x1,y1,x2,y2,x3,y3,surf)
c
c
      implicit double precision (a-h,o-z)

      vz = (y1-y3)*(x1-x2) - (x1-x3)*(y1-y2)
      
      surf = dabs(vz)/2.d0 

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     volume du tetra defini par ses4 points 
c
c
      function VolTetra(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4)
c
c
c
      implicit double precision (a-h,o-z)
c
c     normale a 123 orientee vers 4 et surface 123
c
      call VecNormalT(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &                      vx,vy,vz,surf)
c
c     hauteur 
c     
      s = vx*(x1-x4)+vy*(y1-y4)+vz*(z1-z4)
      s = dabs(s)

      VolTetra = s*surf/3.d0

      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     calcul du vecteur normal unitaire d'un triangle defini par les 
c     coordonnees de ses trois points sortant par rapport 
c     au point x,y,z > vx,vy,vz 
c
c     + surface du triangle = surf 
c
c
      subroutine VecNormalT(x1,y1,z1,x2,y2,z2,x3,y3,z3,x,y,z,
     &                      vx,vy,vz,surf)
c
c
      implicit double precision (a-h,o-z)

      xt = (x1+x2+x3)/3.d0
      yt = (y1+y2+y3)/3.d0
      zt = (z1+z2+z3)/3.d0


      vx = (z1-z3)*(y1-y2) - (y1-y3)*(z1-z2)
      vy = (x1-x3)*(z1-z2) - (z1-z3)*(x1-x2)
      vz = (y1-y3)*(x1-x2) - (x1-x3)*(y1-y2)

      s = dsqrt(vx**2+vy**2+vz**2)
      
      surf = s/2.d0 

      vx = vx/s
      vy = vy/s
      vz = vz/s

      s = (xt-x)*vx+(yt-y)*vy+(zt-z)*vz

      if (s.lt.0.0) then 
         vx = - vx
         vy = - vy
         vz = - vz
      endif   


      return
      end       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c
c     Inverse une matrice (3,3) par 1/det CoMatrice^t 
c
c
      subroutine InvA(A,Ainv)
c
c
      implicit double precision (a-h,o-z)

      dimension A(3,3)
      dimension Ainv(3,3), coA(3,3)
c
ccccccccccccccccccccc
c
      det = A(1,1)*( A(3,3)*A(2,2)-A(3,2)*A(2,3) )
      det = det - A(1,2)*( A(3,3)*A(2,1)-A(3,1)*A(2,3) ) 
      det = det + A(1,3)*( A(3,2)*A(2,1)-A(3,1)*A(2,2) ) 

      coA(1,1) = + ( A(3,3)*A(2,2)-A(3,2)*A(2,3) )/det 
      coA(1,2) = - ( A(3,3)*A(2,1)-A(3,1)*A(2,3) )/det 
      coA(1,3) = + ( A(3,2)*A(2,1)-A(3,1)*A(2,2) )/det  

      coA(2,1) = - ( A(3,3)*A(1,2)-A(3,2)*A(1,3) )/det 
      coA(2,2) = + ( A(3,3)*A(1,1)-A(3,1)*A(1,3) )/det 
      coA(2,3) = - ( A(3,2)*A(1,1)-A(3,1)*A(1,2) )/det  

      coA(3,1) = + ( A(2,3)*A(1,2)-A(2,2)*A(1,3) )/det 
      coA(3,2) = - ( A(2,3)*A(1,1)-A(2,1)*A(1,3) )/det 
      coA(3,3) = + ( A(2,2)*A(1,1)-A(2,1)*A(1,2) )/det  


      Ainv(1,1) = coA(1,1)
      Ainv(2,2) = coA(2,2)
      Ainv(3,3) = coA(3,3)

      Ainv(1,2) = coA(2,1)
      Ainv(2,1) = coA(1,2)

      Ainv(1,3) = coA(3,1)
      Ainv(3,1) = coA(1,3)

      Ainv(3,2) = coA(2,3)
      Ainv(2,3) = coA(3,2)
c
c     test 
c
c      do i=1,3
c      do j=1,3
c         s = 0.d0
c         do l=1,3
c            s = s + A(i,l)*Ainv(l,j)
c         enddo
c         write(*,*)' i j ',i,j,s
c      enddo
c      enddo
c
c
c
cccccccccccccccccccccc
      return
      end       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
