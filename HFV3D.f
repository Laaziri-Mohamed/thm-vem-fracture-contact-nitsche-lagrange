c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                                                                                                       
c                         
c          Schema HFV 3D SUR MAILLAGE POLYHEDRIQUE, MODELE A PRESSION DISCONTINUE 
c
c            NbDim = 3 (include_parameter)
c
c
c     En entree:
c
c              perm: permabilite matrice / Viscof 
c      
c
c               epaisseurf 
c
c               permfdf: permeabilite tangentielle de la fracture * epaisseurf/Viscof 
c
c               permn: permeabilite normale de la fracture
c
c               Correction du Z interface mf avec l'epaisseur frac pour la gravite 
c
cccccccccccccccccccccc
c
c     En sortie 
c
c     TransCell, TransFaceFrac et TransCellbyFaceFrac 
c
c     XCell, XFace, XFaceFrac, XInc
c
c     VolCell, VolFaceFrac 
c
c     NumInterfacebyCell, NumInterfacebyFaceFrac 
c      
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
      subroutine HFV3D(CondThermique,
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
     &     perm,permfdf,permfn,epaisseurf, 
     &     XCell,XFace,XFaceFrac,XInc,
     &     VolCell,VolFaceFrac,
     &     NbInterfacebyCell,NumInterfacebyCell,
     &     NbInterfacebyFaceFrac,NumInterfacebyFaceFrac,
     &     TransCell,TransFaceFrac,TransCellbyFaceFrac,
     &     TransCellFourier,TransFaceFracFourier,
     &     TransCellbyFaceFracFourier,
     &     GCell,GFaceFrac)
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
c
c     pour HFV (toutes les faces y compris les faces frac)
c     num face  > num inc 
c
        dimension NumIncFace(NbFace)
c
c     pour HFV (toutes les aretes fracture)
c
c     Ind 11 arete frac DIR
c     Ind 10 arete frac non DIR 
c     Ind 0 non frac non DIR
c     Ind 1 non frac DIR 
c
        dimension IndIncArete(NbArete)

c     num arete  > num inc 
        dimension NumIncArete(NbArete)
c
c
c        
c     pour HFV discontinu 
c     num face frac  > num inc Ksigma et Lsigma 
c        
        dimension NumIncCellbyFaceFrac(NbFaceFrac,2)              
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
c     permeabilite normale pour HFV discontinu 
c        
        dimension permfdf(NbFaceFrac)
       dimension permfn(NbFaceFrac)         
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     WorkSpaces locaux 
c
ccccccccccccccccccccccccccccccccccccccccccc
c
c
c     Par maille: 
c     VolKi : volume du cone
c     Surfi:  surface des faces, 
c     XKi = Xfi-XK 
c     Xfi = XFace 
c     VecNormalKi: normale unitaire sortante  
c
        dimension VolKi(NbFacebyCellMax)
        dimension Surfi(NbFacebyCellMax)
        dimension XKi(NbFacebyCellMax,NbDim)
        dimension Xfi(NbFacebyCellMax,NbDim)
        dimension VecNormalKi(NbFacebyCellMax,NbDim)

        dimension BKij(NbFacebyCellMax,NbFacebyCellMax)

        double precision, dimension(:), allocatable :: EtaK 

c
c     Gradient par cone K,sig 
c
c     grad_Ksig_iface = \sum_(j = face K) GKij(iface,j)(uj-uk)
c
c     le num de la face est celui dans la numerotation locale des faces  
c     de la maille suivi des numeros de face frac 
c
        dimension GKij(NbFacebyCellMax,NbFacebyCellMax,NbDim)
c
c
c     Par face frac: 
c     SurfSi : surface du triangle 
c     Sizei:  longueur des aretes
c     XSi = Xai - XS 
c     Xai = centre de l'arete  
c     VecNormalKi: normale unitaire sortante dans le plan du triangle 
c 
        dimension SurfSi(NbAretebyFaceMax)
        dimension Sizei(NbAretebyFaceMax)
        dimension XSi(NbAretebyFaceMax,NbDim)
        dimension Xai(NbAretebyFaceMax,NbDim)
        dimension VecNormalSi(NbAretebyFaceMax,NbDim)

        dimension BSij(NbAretebyFaceMax,NbAretebyFaceMax)


       double precision, dimension(:), allocatable :: EtaS

c
c     Gradient tangentiel par cone sig, arete i (triangle) 
c     Grad_sigi = \sum_{arete i de la face } GSij(iarete,j)(uf-uj)
c
        dimension GSij(NbAretebyFaceMax,NbAretebyFaceMax,NbDim)
c
c     workspace 
c
        dimension UU(NbAretebyFaceMax)
        dimension X(NbDim)
c
c
       double precision, dimension(:), allocatable :: SurfaceFaceFrac        
c
c       
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        allocate(EtaK(NbCell))
        allocate(EtaS(NbFaceFrac))

        allocate(SurfaceFaceFrac(NbFaceFrac))

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

           inc = NumIncFace(k)
           do m = 1,NbDim
              XFace(k,m) = XFace(k,m)/dfloat(NbNodebyFace(k))
              XInc(inc,m) = XFace(k,m)
           enddo


           if (IndFaceFrac(k).eq.1) then 
              kfrac = NumFaceVersFaceFrac(k)
              inc = NumIncFaceFrac(kfrac)
              do m = 1,NbDim
                 XFaceFrac(kfrac,m) = XFace(k,m)
                 XInc(inc,m) = XFace(k,m)
              enddo

              incKsig = NumIncCellbyFaceFrac(kfrac,1)
              incLsig = NumIncCellbyFaceFrac(kfrac,2)              
              XInc(incKsig,:) = XFace(k,:)
              XInc(incLsig,:) = XFace(k,:)
              
           endif

c           write(*,*)' Centre Face ',k,(XFace(k,m),m=1,NbDim)
        enddo
c
c     XInc Aretes 
c
        do i=1,NbArete
           is1 = NumNodebyArete(i,1)
           is2 = NumNodebyArete(i,2)
           if (IndIncArete(i).ge.10) then 
              do m = 1,NbDim
               XInc(NumIncArete(i),m) = (XNode(is1,m)+XNode(is2,m))/2.d0
              enddo           
            endif
        enddo
        
c        write(*,*)
c        do i=1,NbCV
c           write(*,*)' XInc ',i,(XInc(i,m),m=1,NbDim)
c        enddo
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
           nfk = NbFacebyCell(k)
c
c     NbInterfacebyCell
c
           nik = nfk 
           NbInterfacebyCell(k) = nik 
c
c     NumInterfacebyCell
c

           do in = 1,nfk 
              nf = NumFacebyCell(k,in)
              
              if (IndFaceFrac(nf).eq.1) then
                 nfrac = NumFaceVersFaceFrac(nf)
                 k1 = NumCellbyFace(nf,1)
                 k2 = NumCellbyFace(nf,2)
                 if (k1.eq.k) then
                    numinc = NumIncCellbyFaceFrac(nfrac,1)
                 elseif (k2.eq.k) then
                    numinc = NumIncCellbyFaceFrac(nfrac,2)                    
                 else
                    write(*,*)' cellbyface non trouvee ',k,k1,k2
                    stop
                 endif
              else 
                 numinc = NumIncFace(nf)
              endif

              NumInterfacebyCell(k,in) = numinc              
              
           enddo

c           write(*,*)
c           write(*,*)' IK ',k,nik,(NumInterfacebyCell(k,i),i=1,nik)

           
c           EtaK(k) = 1.d0/dfloat(NbDim)
           EtaK(k) = 1.d0/dsqrt(dfloat(NbDim))

cccccccccccccccccccccccccccc

           xk = XCell(k,1)
           yk = XCell(k,2)
           zk = XCell(k,3)
c
c
ccccccccccccccccccccccccccccc
c
c          boucle sur les faces de la maille k
c
           do jf = 1,nfk 
              numf = NumFacebyCell(k,jf)
              xf = XFace(numf,1)
              yf = XFace(numf,2)
              zf = XFace(numf,3)

              Surfi(jf) = 0.d0
              si = 0.d0
              do m=1,NbDim
                 Xfi(jf,m) = 0.d0
                 VecNormalKi(jf,m) = 0.d0
              enddo

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
c     centre du tetra X             
c
                 x = (x1+x2+xk+xf)/4.d0
                 y = (y1+y2+yk+yf)/4.d0
                 z = (z1+z2+zk+zf)/4.d0
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

                 VecNormalKi(jf,1) = VecNormalKi(jf,1) + v12fx
                 VecNormalKi(jf,2) = VecNormalKi(jf,2) + v12fy
                 VecNormalKi(jf,3) = VecNormalKi(jf,3) + v12fz

              enddo

              Surfi(jf) = dsqrt( VecNormalKi(jf,1)**2 + 
     &              VecNormalKi(jf,2)**2 + VecNormalKi(jf,3)**2 )

              VecNormalKi(jf,1) = VecNormalKi(jf,1)/Surfi(jf)
              VecNormalKi(jf,2) = VecNormalKi(jf,2)/Surfi(jf)
              VecNormalKi(jf,3) = VecNormalKi(jf,3)/Surfi(jf)

              Xfi(jf,1) = Xfi(jf,1)/si
              Xfi(jf,2) = Xfi(jf,2)/si
              Xfi(jf,3) = Xfi(jf,3)/si 

              XKi(jf,1) = Xfi(jf,1) - xk
              XKi(jf,2) = Xfi(jf,2) - yk
              XKi(jf,3) = Xfi(jf,3) - zk

              s = 0.d0
              do m = 1,NbDim
                 s = s + XKi(jf,m)*VecNormalKi(jf,m)
              enddo
              VolKi(jf) = dabs(s)*Surfi(jf)/dfloat(NbDim)

c              write(*,*)' VolKi ',jf,VolKi(jf)
c              write(*,*)' Surfi ',jf,Surfi(jf)
c              write(*,*)' Xfi ',(Xfi(jf,m),m=1,NbDim)
c              write(*,*)' vec normal ',(VecNormalKi(jf,m),m=1,NbDim)
c
c
c
c     on corrige le centre des faces et la surface des faces frac 
c
              incf = NumIncFace(numf) ! si frac il s'agit de l'inconnue face frac 
              do m=1,NbDim
                 XFace(numf,m) = Xfi(jf,m)
                 XInc(incf,m) = Xfi(jf,m)
              enddo
c
c     Correction des XInc pour les incKsig de la maille k ici (si sig = frac) 
c
              if (IndFaceFrac(numf).eq.1) then

                 numffrac = NumFaceVersFaceFrac(numf)

                 k1 = NumCellbyFace(numf,1)
                 k2 = NumCellbyFace(numf,2)

                 inck1sig = NumIncCellbyFaceFrac(numffrac,1)
                 inck2sig = NumIncCellbyFaceFrac(numffrac,2)                      


                 if (k1.eq.k) then
                    incksig = inck1sig
                 elseif (k2.eq.k) then
                    incksig = inck2sig           
                 else
                    write(*,*)' cellbyface non trouvee ',k,k1,k2
                    stop
                 endif
                 
                 XInc(incksig,:) = XInc(incf,:)              
     &                - VecNormalKi(jf,:)*epaisseurf(numffrac)/2.d0  ! correction pour la gravite 
                 
              endif
              
              if (IndFaceFrac(numf).eq.1) then 
               numffrac = NumFaceVersFaceFrac(numf)
               do m=1,NbDim
                  XFaceFrac(numffrac,m) = Xfi(jf,m)
               enddo               
               VolFaceFrac(numffrac) = Surfi(jf)*epaisseurf(numffrac)
               SurfaceFaceFrac(numffrac) = Surfi(jf)
c              write(*,*)' vol face frac ',numffrac,VolFaceFrac(numffrac)               
              endif

           enddo

           s = 0.d0
           do jf = 1,nfk
              s = s + VolKi(jf)
           enddo
           VolCell(k) = s

c           write(*,*)
c           write(*,*)' VolCell ',k,VolCell(k)
c
c
c
c
c
c     calcul de BKIJ = m_i delta_i_j - m_i m_j/m_k * n_k_j.(x_i-x_k)
c
c
c
              do ni = 1,nfk
                 do nj = 1,nfk
                    s = 0.d0
                    do m = 1,NbDim 
                       s = s + XKi(ni,m)*VecNormalKi(nj,m)
                    enddo
                    BKij(ni,nj) = - s*Surfi(ni)*Surfi(nj)/VolCell(k)
                 enddo
                 BKij(ni,ni) = BKij(ni,ni) + Surfi(ni)
              enddo
c
c
c
c     calcul de GKij = m_K_i/mk * m_j * n_K_j + eta_k *B_k_i_j * n_k_i
c
c
c
              do ni = 1,nfk
                 do nj = 1,nfk
                    do m = 1,NbDim
                       GKij(ni,nj,m) = VolKi(ni)*Surfi(nj)/VolCell(k)
     &                      *VecNormalKi(nj,m) 
     &                      + EtaK(k)*BKij(ni,nj)*VecNormalKi(ni,m)
c                    write(*,*)' GKIJ ',ni,nj,m,GKij(ni,nj,m)
                    enddo
                 enddo
              enddo
c     
c
c
c     test schema GKIJ: gradient exact sur sotions affines 
c
c     u = X1 + 2*X2 + 3*X3  
c 
              do ni = 1,nfk       
c                write(*,*)' maille k face ni ',k,ni
                 do m = 1,NbDim
                    s = 0.d0
                    do nj = 1,nfk
                       f = 0.d0
                       do n = 1,NbDim
                          f = f + XKi(nj,n)*dfloat(n)
                       enddo
                       s = s + GKij(ni,nj,m)*f/VolKi(ni)
                    enddo     
                    if (dabs(s-dfloat(m)).gt.1.0E-4) then 
                     write(*,*)' maille k face ni ',k,ni
                     write(*,*)' grad = 1 ',m,s-dfloat(m) 
                     write(*,*)' erreur gradient vol '
                     stop
                    endif
                 enddo
c                 write(*,*)
              enddo
c
c
c
cccccccccccccccccccccccccccccc

           do ni = 1,nfk
           do nj = 1,nfk
              sp = 0.d0
              do nl = 1,nfk
                 do n = 1,NbDim
                 do m = 1,NbDim
                  sp = sp + perm(k,n,m)
     &                   *GKij(nl,nj,n)*GKij(nl,ni,m)/VolKi(nl)             
                 enddo
                 enddo
              enddo
              
              TransCell(k,ni,nj) = sp
                         
           enddo
        enddo



           do ni = 1,nfk
           do nj = 1,nfk
              st = 0.d0
              do nl = 1,nfk
                 do n = 1,NbDim

                  st = st +        !thermal conductivity to be added !!!!!
     &                GKij(nl,nj,n)*GKij(nl,ni,n)/VolKi(nl)                  

                 enddo
              enddo
              
              TransCellFourier(k,ni,nj) = st * CondThermique
              
           enddo
           enddo
        
            
c           write(*,*)
c           do ni = 1,nfk
c              inci = NumInterfacebyCell(k,ni)
c           do nj = 1,nfk
c              incj = NumInterfacebyCell(k,nj)
c              write(*,*)' TransCell ',
c     &          k,ni,nj,inci,incj,TransCell(k,ni,nj) 
c           enddo
c           enddo


c
c     lumping ccccccccccccccccccccccccccc
c           
c$$$           do ni = 1,nfk
c$$$              s = 0.d0 
c$$$              do nj = 1,nfk
c$$$
c$$$                 s = s + dabs(TransCell(k,ni,nj))
c$$$                 
c$$$              enddo
c$$$              
c$$$              TransCell(k,ni,:) = 0.d0
c$$$              TransCell(k,ni,ni) = s
c$$$              
c$$$              
c$$$           enddo

cccccccccccccccccccccccccccccc           

           
cccccccccccccccccccccccccccccc
c
c     calcul du gradient par maille : moyenne au pro rata des volumes des cones 
c
           do nj = 1,nfk
              do m=1,NbDim
                 GCell(k,nj,m) = 0.d0
              enddo
           enddo
c
           s = 0.d0
           do ni = 1,nfk
              s = s + VolKi(ni)
              do m = 1,NbDim
              do nj = 1,nfk
                 GCell(k,nj,m) = GCell(k,nj,m) 
     &                          - GKij(ni,nj,m)
              enddo     
              enddo
           enddo   
c      
           do m = 1,NbDim
              do nj = 1,nfk
                 GCell(k,nj,m) = GCell(k,nj,m)/s
              enddo     
           enddo
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
c
c        
ccccccccccccccccccccccccccccccccc
c
c     boucle sur les face frac 
c
ccccccccccccccccccccccccccccccc
c
c        write(*,*)' NbFaceFrac ',NbFaceFrac
c
c
c
        do ifrac = 1,NbFaceFrac
           numf = NumFaceFracVersFace(ifrac)
           inc = NumIncFaceFrac(ifrac)
c
c     naf : nb d'aretes par face 
c
           naf = NbAretebyFace(numf)


           xf = XFace(numf,1)
           yf = XFace(numf,2)
           zf = XFace(numf,3)
c
c     NbInterfacebyFaceFrac
c  

           NbInterfacebyFaceFrac(ifrac) = naf
c
c     NumInterfacebyFaceFrac: dans l'ordre des aretes de la face 
c 
           do is = 1,naf
              nums = NumAretebyFace(numf,is)
              numinc = NumIncArete(nums)
              NumInterfacebyFaceFrac(ifrac,is) = numinc 
           enddo

c           write(*,*)
c           write(*,*)' Ifrac ',ifrac,naf,
c     &       (NumInterfacebyFaceFrac(ifrac,is),is=1,naf)


c           EtaS(ifrac) = 1.d0/(dfloat(NbDim-1))
           EtaS(ifrac) = 1.d0/dsqrt(dfloat(NbDim-1))

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
c     Surface du triangle 12f et vecteur normal a l'arete 12 ds le plan 12f 
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

                 s = dsqrt(vx**2+vy**2+vz**2) 

                 vx = vx/s
                 vy = vy/s
                 vz = vz/s

                 d12f = (xf-x1)*vx + (yf-y1)*vy + (zf-z1)*vz 
                 d12 = dsqrt(s12)

                 surfSi(ia) = abs(d12f)*d12/dfloat(NbDim-1)
                 Sizei(ia) = d12

                 VecNormalSi(ia,1) = vx 
                 VecNormalSi(ia,2) = vy 
                 VecNormalSi(ia,3) = vz 

                 Xai(ia,1) = (x1+x2)/2.d0
                 Xai(ia,2) = (y1+y2)/2.d0
                 Xai(ia,3) = (z1+z2)/2.d0

                 XSi(ia,1) = Xai(ia,1) - xf 
                 XSi(ia,2) = Xai(ia,2) - yf 
                 XSi(ia,3) = Xai(ia,3) - zf 

c              write(*,*)' SurfSi ',ia,SurfSi(ia)
c              write(*,*)' Sizei ',jf,Sizei(ia)
c              write(*,*)' Xai ',(Xai(ia,m),m=1,NbDim)
c              write(*,*)' vec normal ',(VecNormalSi(ia,m),m=1,NbDim) 
c
c     fin boucle sur les aretes de la face frac 
c
              enddo
c
c
           s = 0.d0
           do ia = 1,naf
              s = s + SurfSi(ia)
           enddo
           surf = s
c
c
c
c
c     calcul de BSIJ = m_i delta_i_j - m_i m_j/m_k * n_k_j.(x_i-x_k)
c
c
c
              do ni = 1,naf
                 do nj = 1,naf
                    s = 0.d0
                    do m = 1,NbDim 
                       s = s + XSi(ni,m)*VecNormalSi(nj,m)
                    enddo
                    BSij(ni,nj) = - s*Sizei(ni)*Sizei(nj)/surf
                 enddo
                 BSij(ni,ni) = BSij(ni,ni) + Sizei(ni)
              enddo
c
c
c
c     calcul de GSij = m_K_i/mk * m_j * n_K_j + eta_k *B_k_i_j * n_k_i
c
c
c
              do ni = 1,naf
c                 write(*,*)' ni = ',ni
                 do nj = 1,naf
                    do m = 1,NbDim
                       GSij(ni,nj,m) = SurfSi(ni)*Sizei(nj)/surf
     &                      *VecNormalSi(nj,m) 
     &                      + EtaS(ifrac)*BSij(ni,nj)*VecNormalSi(ni,m)
c                    write(*,*)' GSIJ sur SurfSi ',
c     &                 nj,m,GSij(ni,nj,m)/SurfSi(ni)
                    enddo
                 enddo
c                 write(*,*)
              enddo
c     
c
c
c     test schema GSIJ: gradient tangentiel 
c     on teste sur des vecteurs du plan de la face : XSi = Xai - Xf 
c
c     grad_tau u.XSi = uai - uf 
c
c     avec 
c
c     grad_tau u = 1/SurfSi(ni) \sum_j GSij (uaj-uf)
c
c 
              uf = 2*xf + yf - zf + 1.d0
              do ni = 1,naf
                 UU(ni) = 2*Xai(ni,1) + Xai(ni,2) - Xai(ni,3) + 1.d0
              enddo

              do ni = 1,naf
                 s = 0.d0
                 do nj = 1,naf
                    do m = 1,NbDim
                     s = s + XSi(ni,m)*GSij(ni,nj,m)
     &                       *(UU(nj)-uf)/SurfSi(ni)
                    enddo
                 enddo
                 if (abs(s-UU(ni)+uf).gt.1.0E-6) then 
                  write(*,*)' test grad tan ',i,s,UU(ni)-uf,s-UU(ni)+uf
                 stop
                 endif
              enddo
c
c     calcul des Tsij transmissivites (permfdf isotrope) 
c                 
              do ni = 1,naf
                 do nj = 1,naf
                    sp = 0.d0
                    st = 0.d0
                    do nl = 1,naf
                       do n = 1,NbDim
                          sp = sp + permfdf(ifrac)
     &                         *GSij(nl,nj,n)*GSij(nl,ni,n)/SurfSi(nl)

                          st = st +                  
     &                         GSij(nl,nj,n)*GSij(nl,ni,n)/SurfSi(nl)                          
                       enddo
                    enddo
                    TransFaceFrac(ifrac,ni,nj) = sp

                    TransFaceFracFourier(ifrac,ni,nj) = st
     &                          * CondThermique
                    
                 enddo
              enddo

c
c     lumping ccccccccccccccccccc
c           
c$$$           do ni = 1,naf
c$$$              s = 0.d0 
c$$$              do nj = 1,naf
c$$$
c$$$                 s = s + dabs(TransFaceFrac(ifrac,ni,nj))
c$$$                 
c$$$              enddo
c$$$              
c$$$              TransFaceFrac(ifrac,ni,:) = 0.d0
c$$$              TransFaceFrac(ifrac,ni,ni) = s
c$$$              
c$$$              
c$$$           enddo

cccccccccccccccccccccccccccccccc
           
              
c              write(*,*)
c              do ni = 1,naf
c                 inci = NumInterfacebyFaceFrac(ifrac,ni)
c                 do nj = 1,naf
c                    incj = NumInterfacebyFaceFrac(ifrac,nj)
c                    write(*,*)' TransFaceFrac ',
c     &               ifrac,ni,nj,inci,incj,TransFaceFrac(ifrac,ni,nj) 
c                 enddo
c              enddo
c
c              stop
c
c
c     calcul du gradient par face frac : 
c         moyenne au pro rata des SurfSi 
c
           do ni = 1,naf
              do m=1,NbDim
                 GFaceFrac(ifrac,ni,m) = 0.d0
              enddo
           enddo
c
           s = 0.d0
           do ia = 1,naf 
              s = s + SurfSi(ia)
              do m = 1,NbDim
              do nj = 1,naf
               GFaceFrac(ifrac,nj,m) = GFaceFrac(ifrac,nj,m) 
     &            - GSij(ia,nj,m)
              enddo     
              enddo
           enddo   
c      
           do m = 1,NbDim
              do nj = 1,naf
                 GFaceFrac(ifrac,nj,m) = GFaceFrac(ifrac,nj,m)/s
              enddo     
           enddo

c
c     Transmissivites TKsigma et TLSigma 
c
           k1 = NumCellbyFace(numf,1)
           k2 = NumCellbyFace(numf,2)

           inck1sig = NumIncCellbyFaceFrac(ifrac,1)
           inck2sig = NumIncCellbyFaceFrac(ifrac,2)

       

           sfp =   SurfaceFaceFrac(ifrac)
     &          *permfn(ifrac)*2.d0/epaisseurf(ifrac)

           sft=  SurfaceFaceFrac(ifrac)  ! fracture normal thermal condutivity to be added
     &         *2.d0/epaisseurf(ifrac)
           
           TransCellbyFaceFrac(ifrac,1) = sfp
           
           TransCellbyFaceFrac(ifrac,2) = sfp

           TransCellbyFaceFracFourier(ifrac,1) = sft * CondThermique
           
           TransCellbyFaceFracFourier(ifrac,2) = sft * CondThermique

           

           
c
c    ccccccccccccccccccccccccccccccc
c
c     fin boucle sur les faces frac 
c
        enddo

cccccccccccccccccccccccccccccccccccc

        write(*,*)' fin HFV '

        deallocate(EtaK)
        deallocate(EtaS)
        deallocate(SurfaceFaceFrac)

cccccccccccccccccccccccccccc
       return
       end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
