cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                                                              
c     Assemblage de la Jacobienne et du Second Membre 
c          
c     En entree: solution courante SolP, SolT 
c
c     En sortie: Jacobienne: AA
c                Second Membre: Sm = - residu fct de la solution courante SolP et SolT  
c
cccccccccccccccccccc     
c
c     Traitement des intersections > = 3 faces 
c
c      -> pour la convection thermique on ecrit le flux Fourier + convection, on decentre entre la face et l'arete 
c
c      -> pour le flux Darcy: on commence par le cas 1, ok pour un liquide peu compressible 
c
c        1. si schema centre (rhof arete), on peut garder la continuite du flux Darcy avec rho0  
c
c        2. si schema decentre ( rhof decentre entre face et arete), il faut ecrire la conservation du flux Darcy avec rhof  
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
      subroutine JacSmPT(
     &     Deltat,       
     &     NbCell,NbFaceFrac,NbFace,NbArete,NbCV,
     &     VolCell,VolFaceFrac,SurfaceFace,XInc,
     &     NumIncCell,NumIncFaceFrac, 
     &     IndIncDir,IndIncNeu,IndIncDirT,IndIncNeuT,NumIncFace, 
     &     NbInterfacebyCell,NumInterfacebyCell,
     &     NbInterfacebyFaceFrac,NumInterfacebyFaceFrac,
     &     NumIncCellbyFaceFrac,NumFaceFracVersFace,     
     &     NumFacebyCell,NumCellbyFace,
     &     NumFaceVersFaceFrac,IndFaceFrac,
     &     NbFaceFracbyArete,NumFaceFracbyArete,
     &     NumAretebyFace,    
     &     IndLigneAA,NColAA,NLigneAA,MemBloc, 
     &     TransCell,TransFaceFrac,TransCellbyFaceFrac,
     &     TransCellFourier,TransFaceFracFourier,
     &     TransCellbyFaceFracFourier,      
     &     PorobyCell,DerPorobyCell,PorobyCell_nm1,
     &     dfbyFaceFrac,dfbyFaceFrac_nm1,
     &     relaxff,relaxmm,
     &     AccP,AccP_nm1,       
     &     SolPDir,SolP,SolP_nm1,SolP_prev,SolP0,
     &     AccT,AccT_nm1,       
     &     SolTDir,SolT,SolT_nm1,SolT_prev,SolT0,    
     &     AA,Sm,
     &     NumCSRK,NumCSRI)
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
        dimension NumFacebyCell(NbCell,NbFacebyCellMax)
        dimension NumCellbyFace(NbFace,2)
        
        dimension NumFaceVersFaceFrac(NbFace)
        dimension IndFaceFrac(NbFace)

        dimension NbFaceFracbyArete(NbArete)
        dimension NumFaceFracbyArete(NbArete,NbFacebyAreteMax)

        dimension NumAretebyFace(NbFace,NbAretebyFaceMax)        
        
c     Solution courante 
c
        dimension SolP(NbCV)

        dimension SolP_nm1(NbCV)

        dimension SolP_prev(NbCV)

        dimension SolP0(NbCV)

        dimension SolT(NbCV)

        dimension SolT_nm1(NbCV)

        dimension SolT_prev(NbCV)

        dimension SolT0(NbCV)           

        
c        
c
c
c     porosite courante (non convergee !!)  
c     donc a recalculer par appel a couplage
c       
       dimension PorobyCell(NbCell)

c     porosite nm1 et der de la porosite courante / P
c       
      dimension PorobyCell_nm1(NbCell)
      dimension DerPorobyCell(NbCell,2)

c
c     aperture au contact par face frac 
c      
      dimension dfcontact(NbFaceFrac)      
c
c     aperture courante et nm1 
c      
      dimension dfbyFaceFrac(NbFaceFrac)
      dimension dfbyFaceFrac_nm1(NbFaceFrac)
c
c     Acc Darcy nm1 
c      
      dimension AccP_nm1(NbCV)

c
c     Acc Fourier nm1 
c      
      dimension AccT_nm1(NbCV)  
c
c      
ccccccccccccccccccccccccccccccccccccc
c      
c     coordonnee des inconnues
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
c     Surface by Face 
c
        dimension SurfaceFace(NbFace)
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
c     pour HFV (toutes les faces y compris les faces frac)
c     num face  > num inc 
c
        dimension NumIncFace(NbFace)        
c        
c     pour HFV discontinu 
c     num face frac  > num inc Ksigma et Lsigma 
c        
        dimension NumIncCellbyFaceFrac(NbFaceFrac,2)              
c        
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
c     num inc > indice 0 ou 1 (neu homog) ou 2 (neu non homog) 
c
        dimension IndIncNeu(NbCV)
        dimension IndIncNeuT(NbCV)           
c
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
c       num face frac > num face 
        dimension NumFaceFracVersFace(NbFaceFrac)        
c
c
c     Structure CSR 
c     On met l'elt diagonal en premier 
c
        dimension IndLigneAA(NbCV+1)
        dimension NLigneAA(MemBloc)
        dimension NColAA(MemBloc)
c
c
c
c     en entree 
c
c     CL Dirichlet aux inconnues interfaces de bord DIR 
c
        dimension SolPDir(NbCV)

        dimension SolTDir(NbCV)
c
c
c
c     coeff de relaxation ff et mm 
c
       dimension relaxff(NbFaceFrac)
       
       dimension relaxmm(NbCell)
c
c       
cccccccccccccccccccccccccccccccccccccccc
c
c     SORTIES 
c
cccccccccccccccccccccccccccccccccccccccc
c
c     Second Membre du systeme et solution 
c
       dimension Sm(NbCV,NbIncPrim)
c
c     Jacobienne Mailles + Noeuds 
c
        dimension AA(MemBloc,NbIncPrim,NbIncPrim)
c
c
c     Acc Darcy et Fourier avec la solution courante (non convergee !!)  
c     donc a recalculer par appel a couplage et a JacSm
c       
        dimension AccP(NbCV)


        dimension AccT(NbCV)      
c
ccccccccccccccccccccccccccccccccccccccc
c
c     WorkSpaces 
c
c
c      Indices m CSR des j (indice col) des lignes K et I    
c
       dimension NumCSRK(NbCV), NumCSRI(NbCV)
c
c
c
       dimension X(NbDim)
c
c     Densite et Derivees de la densite par rapport a P , T 
c       
      double precision, dimension(:), allocatable :: Rho       
      double precision, dimension(:,:), allocatable :: DerRho
      dimension Drhof(2)



      double precision, dimension(:), allocatable :: Rho_nm1 
c
c
c     Energie Interne et Derivees par rapport a P , T 
c       
      double precision, dimension(:), allocatable :: EInt       
      double precision, dimension(:,:), allocatable :: DerEInt
      dimension DEf(2)

      double precision, dimension(:), allocatable :: EInt_prev 
c
c
c
c     Enthalpie et Derivees par rapport a P , T 
c       
      double precision, dimension(:), allocatable :: Enth       
      double precision, dimension(:,:), allocatable :: DerEnth
      dimension DHf(2)       
c
c            
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Remplissage Jacobienne et second membre 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     gravite avec densite constante 
c       
       grav = f_gravite()
       rho0  = f_RHOREF()


c
c        Coefficients
c       

      b_coef = f_biot()
      
      alpha0 = f_alpha0()

      rK0 =  f_K0()
       
c
c
c     coefficients de relaxation matrice et fracture  
c      
c      Crf = f_Crf()
c      Crm = f_Crm()
c
c
c      
ccccccccccccccccccccccccccc       
c
c     mise a zero de AA et Sm 
c
 
       AA(:,:,:) = 0.d0
 
       Sm(:,:) = 0.d0

       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Lois thermo et derivees des lois thermo / p,T 
c
       allocate(Rho(NbCV))
       Rho = 0.d0

       allocate(Rho_nm1(NbCV))
       Rho_nm1 = 0.d0        


       allocate(DerRho(NbCV,2))
       DerRho = 0.d0 
       
       do k=1,NbCV
          pk = SolP(k)
          Tk = SolT(k)
          Rho(k) = f_rho(pk,Tk)
          DerRho(k,1) = f_dprho(pk,Tk)
          DerRho(k,2) = f_dTrho(pk,Tk)

          pkm1 = SolP_nm1(k)
          Tkm1 = SolT_nm1(k)
          
          Rho_nm1(k) = f_rho(pkm1,Tkm1)
          
       enddo


       
       allocate(EInt(NbCV))
       EInt = 0.d0        

       allocate(DerEInt(NbCV,2))
       DerEInt = 0.d0

       allocate(EInt_prev(NbCV))
       EInt_prev = 0.d0
       
       
       do k=1,NbCV
          pk = SolP(k)
          Tk = SolT(k)
          EInt(k) = f_ef(pk,Tk)
          DerEInt(k,1) = f_dpef(pk,Tk)
          DerEInt(k,2) = f_dTef(pk,Tk)

          pkm1 = SolP_nm1(k)
          Tkm1 = SolT_nm1(k)
          EInt_prev(k) = f_ef(pkm1,Tkm1)


         
       enddo

         
       
       allocate(Enth(NbCV))
       Enth = 0.d0        

       allocate(DerEnth(NbCV,2))
       DerEnth = 0.d0 
       
       do k=1,NbCV
          pk = SolP(k)
          Tk = SolT(k)
          Enth(k) = f_hf(pk,Tk)
          DerEnth(k,1) = f_dphf(pk,Tk)
          DerEnth(k,2) = f_dThf(pk,Tk)          
       enddo       
       
c       
c     Chaleur volumique drainee du squelette  J.m^-3 K^-1 
c       
       C_Roche = f_C0()
c       write(*,*)' CR ',C_Roche 
c       
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Assemblage des termes d'accumulation aux mailles et aux faces fracture 
c
ccccccccccccccccccccccccccccccccccccc       
c
c    

         
       
       
       do k=1,NbCell
          inck = NumIncCell(k)

         rhof = Rho(inck)
         Drhof(:) = DerRho(inck,:)
         DEf(:) = DerEInt(inck,:)


         EInt_m1 = EInt_prev(inck)

         rhof_nm1 = Rho_nm1(inck)
         
          AccP(inck) = VolCell(k)*PorobyCell(k)*rhof   ! acc |K|phi*rhof 
 
          AccT(inck) = AccP(inck) * EInt(inck) ! acc |K|phi*rhof * EI  
     &             + VolCell(k) * C_Roche * SolT(inck)     
         
          Sm(inck,1) = Sm(inck,1)
     &         + ( AccP(inck) - AccP_nm1(inck) )/Deltat ! terme en |K| d_t (phi*rhof)

         
          Sm(inck,2) = Sm(inck,2)
     &         + ( AccT(inck) - AccT_nm1(inck) )/Deltat ! terme en |K| d_t (phi rhof * T ) 


cccccccccccccccccccc          
c
c     on rajoute la relaxation matrice 
c          
c
          A =  alpha0 *rK0 / b_coef
          ss = VolCell(k)*relaxmm(k) * rhof_nm1 *
     &      ( ( SolP(inck) - SolP_prev(inck) )
     &         + A * ( SolT(inck) - SolT_prev(inck) ) ) ! relaxation

          
          Sm(inck,1) = Sm(inck,1) + ss/Deltat

          Sm(inck,2) = Sm(inck,2) +   EInt_m1 * ss /Deltat  

ccccccccccccccccccccc          
          
          mkk = IndLigneAA(inck) + 1

          AA(mkk,1,:) = AA(mkk,1,:)
     &          + VolCell(k)*DerPorobyCell(k,:)*rhof/Deltat 
          
          AA(mkk,1,:) = AA(mkk,1,:)
     &         + VolCell(k)*PorobyCell(k)*Drhof(:)/Deltat 
          
                   
          AA(mkk,1,1) = AA(mkk,1,1)
     &         + VolCell(k)*relaxmm(k)* rhof_nm1/Deltat ! relaxation, der / P

          AA(mkk,1,2) = AA(mkk,1,2)
     &         + VolCell(k)*relaxmm(k)* rhof_nm1 *A/Deltat ! relaxation, der / T
          

          

          AA(mkk,2,:) = AA(mkk,2,:)
     &         + VolCell(k)*DerPorobyCell(k,:)*rhof/Deltat ! der AccT ds energie porosite 
     &           *  EInt(inck)
          
          AA(mkk,2,:) = AA(mkk,2,:)
     &         + VolCell(k)*PorobyCell(k)*Drhof(:)/Deltat  ! der / P,T ds energie rhof         
     &           *  EInt(inck)

         

          AA(mkk,2,:) = AA(mkk,2,:) + DEf(:) * AccP(inck)/Deltat ! der / P,T ds energie EI
          
          AA(mkk,2,2) = AA(mkk,2,2) + VolCell(k) * C_Roche / Deltat ! der / T ds energie C_Roche
          

           AA(mkk,2,1) = AA(mkk,2,1)
     &         + VolCell(k)*relaxmm(k) *  rhof_nm1 * EInt_m1/Deltat ! relaxation, der / P


           
           AA(mkk,2,2) = AA(mkk,2,2)
     &         + VolCell(k)*relaxmm(k)*  rhof_nm1 * A * EInt_m1/Deltat ! relaxation, der / T
          
       enddo


       do k=1,NbFaceFrac
          inck = NumIncFaceFrac(k)

          nf = NumFaceFracVersFace(k)

          rhof = Rho(inck)
          Drhof(:) = DerRho(inck,:)
          DEf(:)   = DerEInt(inck,:)


           rhof_nm1 = Rho_nm1(inck)


          
          AccP(inck) = SurfaceFace(nf)*dfbyFaceFrac(k)*rhof ! acc |sigma|df *rhof


          AccT(inck) = AccP(inck) *  EInt(inck) ! acc |sigma| df*rhof * EI 

          
          
          Sm(inck,1) = Sm(inck,1)
     &         + ( AccP(inck) - AccP_nm1(inck) )/Deltat ! |sigma| d_t (df*rhof )

          Sm(inck,2) = Sm(inck,2)
     &         + ( AccT(inck) - AccT_nm1(inck) )/Deltat ! terme en |sigma| d_t ( df*rhof * T) 

c
cccccccccccccccccccc          
c
c     on rajoute la relaxation 
c          

          
           ss = rhof_nm1*relaxff(k)*( SolP(inck) - SolP_prev(inck) )
     
          
          Sm(inck,1) = Sm(inck,1) + SurfaceFace(nf)*ss/Deltat

c
cccccccccccccccccccc

          

          mkk = IndLigneAA(inck) + 1

          AA(mkk,1,:) = AA(mkk,1,:)
     &         + SurfaceFace(nf)*dfbyFaceFrac(k)*Drhof(:)/Deltat
  
          
          AA(mkk,1,1) = AA(mkk,1,1)
     &         + SurfaceFace(nf)*rhof_nm1*relaxff(k)/Deltat ! relaxation 

          

          AA(mkk,2,:) = AA(mkk,2,:)
     &         + SurfaceFace(nf)*dfbyFaceFrac(k)*Drhof(:)/Deltat ! der / P,T ds energie 
     &         * EInt(inck)
          
          AA(mkk,2,:) = AA(mkk,2,:) + DEf(:) * AccP(inck)/Deltat ! der / T ds energie 
          

          
c          write(*,*)' Sm Acc ',k,Sm(inck)
          
       enddo       



       
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Assemblage des flux 
c
cccccccccccccccccccccccccccccccccccccc
c     boucle sur les mailles 
c
       do k=1,NbCell
          nik = NbInterfacebyCell(k) 
          inck = NumIncCell(k)
c          
c         on stocke les indices CSR des ind de col de la ligne inck  
c
          do m = IndLigneAA(inck)+1,IndLigneAA(inck+1)
             NumCSRK(NColAA(m)) = m
          enddo

          do ik = 1,nik
             incik = NumInterfacebyCell(k,ik)
             nf = NumFacebyCell(k,ik) ! numero de la face ik 


             
              ! on centre la densite sur l'interface incik 
              rhof = Rho(incik)
              Drhof(:) = DerRho(incik,:)
             
c
c     Calcul du flux de Darcy et du flux de Fourier 
c
             fik = 0.d0
             gik = 0.d0
             do jk = 1,nik
                
                incjk = NumInterfacebyCell(k,jk)
                
                fik = fik + TransCell(k,ik,jk)
     &               *( SolP(inck) - SolP(incjk)
     &               + grav * Rho_nm1(incik)
     &                     * ( XInc(inck,3) - XInc(incjk,3) )  )

                 gik = gik + TransCellFourier(k,ik,jk)
     &               *( SolT(inck) - SolT(incjk) )
                               
              enddo

              ! decentrage pour la convection thermique 
              if (fik.ge.0.d0) then 
                 
                 incup = inck 

              else

                 if (IndFaceFrac(nf).eq.1) then ! face frac 
                    
                    ifrac = NumFaceVersFaceFrac(nf)
                    incup = NumIncFaceFrac(ifrac)

                 else if (NumCellbyFace(nf,2).le.0) then! face de bord 

                    incup = incik

                 else

                   if (NumCellbyFace(nf,1).eq.k) then
                      kup = NumCellbyFace(nf,2)
                   else
                      kup = NumCellbyFace(nf,1)
                   endif
                   
                   incup = NumIncCell(kup)
                    
                 endif
                 
              endif


              Hf = Enth(incup)  ! enthalpie amont 
              DHf(:) = DerEnth(incup,:) ! derivee de l'enthalpie incup 
              
c
c     assemblage du flux fik 
c
c             
c     eq de conservation de la masse 
c             
             Sm(inck,1) = Sm(inck,1) + fik*rhof ! flux darcy avec rhof  

             if (IndIncNeu(incik).eq.2) then
                
                Sm(incik,1) = Sm(incik,1) - fik*rhof ! flux complet de Darcy avec rhof si Neumann non homogene
                
             else 
                
                Sm(incik,1) = Sm(incik,1) - fik*rho0 ! eq de continuite de flux Darcy avec scaling rho0              
                
             endif
c             
c     eq de conservation de l'energie 
c             
             Sm(inck,2) = Sm(inck,2) + gik
     &                               + fik * rhof * Hf ! flux complet Fourier + convection 

             Sm(incik,2) = Sm(incik,2) - gik ! eq de continuite du flux Fourier                
                           
c             
c          
c          on stocke les indices CSR des inc de col de la ligne incik  
c
           do m = IndLigneAA(incik)+1,IndLigneAA(incik+1)
             NumCSRI(NColAA(m)) = m
           enddo

           
           do jk = 1,nik ! boucle jk 
              incjk = NumInterfacebyCell(k,jk)

c          eq inck col inck Darcy 
              mkk = NumCSRK(inck)
              AA(mkk,1,1) = AA(mkk,1,1) + TransCell(k,ik,jk)*rhof

c          eq inck col inck convection thermique              
              AA(mkk,2,1) = AA(mkk,2,1) + TransCell(k,ik,jk)*rhof ! eq inck energie der fik / Pinck 
     &                                    * Hf
c          eq inck col inck Fourier         
              AA(mkk,2,2) = AA(mkk,2,2) +  TransCellFourier(k,ik,jk)

              
c           eq incik col inck  Darcy         
              if (IndIncNeu(incik).eq.2) then                
                 mik = NumCSRI(inck)
                 AA(mik,1,1) = AA(mik,1,1) - TransCell(k,ik,jk)*rhof                                               
              else                  
                 mik = NumCSRI(inck)
                 AA(mik,1,1) = AA(mik,1,1) - TransCell(k,ik,jk)*rho0                 
              endif

                
c           eq incik col inck Fourier 
              mik = NumCSRI(inck)                 
              AA(mik,2,2) = AA(mik,2,2) -  TransCellFourier(k,ik,jk)
                             
c           eq inck col incjk Darcy 
              mkj = NumCSRK(incjk)
              AA(mkj,1,1) = AA(mkj,1,1) - TransCell(k,ik,jk)*rhof

c           eq inck col incjk Convection thermique               
              AA(mkj,2,1) = AA(mkj,2,1) - TransCell(k,ik,jk)*rhof ! eq inck energie der fik / Pincjk 
     &             * Hf

c           eq inck col incjk Fourier          
              AA(mkj,2,2) = AA(mkj,2,2) -  TransCellFourier(k,ik,jk)

                 
c     eq incik col incjk Darcy 
              if (IndIncNeu(incik).eq.2) then                 
                 mij = NumCSRI(incjk)
                 AA(mij,1,1) = AA(mij,1,1) + TransCell(k,ik,jk)*rhof                       
              else                 
                 mij = NumCSRI(incjk)
                 AA(mij,1,1) = AA(mij,1,1) + TransCell(k,ik,jk)*rho0                
              endif
                 
c           eq incik col incjk  Fourier 
              mij = NumCSRI(incjk)
              AA(mij,2,2) = AA(mij,2,2) + TransCellFourier(k,ik,jk)              
              
           enddo ! fin boucle jk 

c          eq inck Darcy, terme rhof 
           mki = NumCSRK(incik)
           AA(mki,1,:) = AA(mki,1,:) + fik * Drhof(:)
           
           AA(mki,2,:) = AA(mki,2,:)
     &          + fik * Drhof(:) * Hf ! eq inck energie, der / P,T ds rhof 
           
c           eq incik Darcy terme rhof si Neumann non homogene 
           if (IndIncNeu(incik).eq.2) then
              mii = NumCSRI(incik)
              AA(mii,1,:) = AA(mii,1,:) - fik * Drhof(:)              
           endif

c          eq inck convection thermique, terme fik rhof Enth(incup) 
           mkup = NumCSRK(incup)
           AA(mkup,2,:) = AA(mkup,2,:) + fik * rhof * DHf(:)  ! eq inck der / p,T Enth incup            
           
                   
          enddo ! fin boucle ik 
       enddo  ! fin boucle mailles k        
c
cccccccccccccccccccccccccccccccccccccc
c     boucle sur les faces frac 
c
c     Il faudra rajouter la distinction entre les aretes avec <= 2 faces et celles >= 3 faces
c     dans ce dernier cas il faut du volume de controle sur l'arete et des flux complets P et T 
c
       
       do k=1,NbFaceFrac
          nik = NbInterfacebyFaceFrac(k) 
          inck = NumIncFaceFrac(k)
          nf = NumFaceFracVersFace(k)


c          Condf = dfbyFaceFrac(k)**3 ! dependence en df**3 (iteration prev)        
          Condf = dfbyFaceFrac_nm1(k)**3 ! dependence en df**3 explicite
          
          Condf_th = dfbyFaceFrac_nm1(k)          
c          
c         on stocke les indices CSR des ind de col de la ligne inck  
c
          do m = IndLigneAA(inck)+1,IndLigneAA(inck+1)
             NumCSRK(NColAA(m)) = m
          enddo

          do ik = 1,nik
             incik = NumInterfacebyFaceFrac(k,ik)
             
             na = NumAretebyFace(nf,ik)

             ! on centre la densite sur l'interface incik 
             rhof = Rho(incik)
             Drhof(:) = DerRho(incik,:)
             
c
c     Calcul du flux de Darcy et du flux de Fourier 
c
             fik = 0.d0
             gik = 0.d0            
             do jk = 1,nik
                incjk = NumInterfacebyFaceFrac(k,jk)
                
                fik = fik + TransFaceFrac(k,ik,jk)*Condf 
     &               *( SolP(inck)  - SolP(incjk)
     &               + grav * Rho_nm1(incik)
     &                 * ( XInc(inck,3) - XInc(incjk,3) )  )

                                
                gik = gik + TransFaceFracFourier(k,ik,jk)
     &               *( SolT(inck)  - SolT(incjk)) * Condf_th 
                              
             enddo
c
c     decentrage pour la convection thermique 
c
             if (fik.ge.0.d0) then

                incup = inck

             else

                if (IndIncDirT(incik).eq.1) then ! dirichlet T

                   incup = incik                    
                
                else if (IndIncNeuT(incik).eq.1) then ! Neumann T 

                   incup = incik          
                   
                else ! arete interne 
                   
                   if (NbFaceFracbyArete(na).eq.1) then ! tip 

                      incup = incik 
                   
                   else if (NbFaceFracbyArete(na).eq.2) then ! 2 faces 
                      
                      if (NumFaceFracbyArete(na,1).eq.k) then
                         ifrac = NumFaceFracbyArete(na,2)
                      else
                         ifrac = NumFaceFracbyArete(na,1) 
                      endif
                      
                      incup = NumIncFaceFrac(ifrac) 
                      
                   else ! au moins 3 faces 
                      
                      incup = incik        
                                            
                   endif
                   
                endif
                
             endif


             Hf = Enth(incup) 
             DHf(:) = DerEnth(incup,:) ! derivee de l'enthalpie incup 
             
c
c     assemblage du flux fik 
c

             
             
             Sm(inck,1) = Sm(inck,1) + fik*rhof ! flux complet darcy avec rhof 

             Sm(incik,1) = Sm(incik,1) - fik*rho0 ! equation de continuite du flux darcy avec scaling rho0
            
             Sm(inck,2) = Sm(inck,2) + gik
     &            + fik * rhof * Hf  !  flux complet fourier + convection 
                                           
                
             if ( (IndIncNeuT(incik).eq.0)
     &            .and.(NbFaceFracbyArete(na).ge.3) ) then ! >= 3 faces et pas sur le bord       

                Sm(incik,2) = Sm(incik,2) - gik
     &                  - fik * rhof * Hf ! flux complet fourier + convection 
                           
             else  ! Neuman ou <= 2 faces 

                Sm(incik,2) = Sm(incik,2) - gik ! flux fourier seul   

             endif                                                             
c          
c          on stocke les indices CSR des inc de col de la ligne incik  
c
             do m = IndLigneAA(incik)+1,IndLigneAA(incik+1)
                NumCSRI(NColAA(m)) = m
             enddo
c
             do jk = 1,nik
                incjk = NumInterfacebyFaceFrac(k,jk)

           ! eq inck col inck Darcy 
                mkk = NumCSRK(inck)
                AA(mkk,1,1) = AA(mkk,1,1)
     &               + TransFaceFrac(k,ik,jk)*Condf*rhof

           ! eq inck col inck convection thermique 
                AA(mkk,2,1) = AA(mkk,2,1)
     &               + TransFaceFrac(k,ik,jk)*Condf*rhof
     &               * Hf 

          ! eq inck col inck Fourier             
                AA(mkk,2,2) = AA(mkk,2,2)
     &               + TransFaceFracFourier(k,ik,jk)
     &               * Condf_th
           
c           eq incik col inck Darcy 
                mik = NumCSRI(inck)
                AA(mik,1,1) = AA(mik,1,1)
     &               -TransFaceFrac(k,ik,jk)*Condf*rho0
 

c           eq incik col inck Fourier 
                mik = NumCSRI(inck)
                AA(mik,2,2) = AA(mik,2,2)
     &               - TransFaceFracFourier(k,ik,jk)
     &               * Condf_th            

           
c          eq inck col incjk Darcy 
                mkj = NumCSRK(incjk)
                AA(mkj,1,1) = AA(mkj,1,1)
     &               - TransFaceFrac(k,ik,jk)*Condf*rhof

c          eq inck col incjk convection thermique            
                AA(mkj,2,1) = AA(mkj,2,1)
     &               - TransFaceFrac(k,ik,jk)*Condf*rhof
     &               * Hf ! eq inck energie, der / Pincjk ds fik 

c          eq inck col incjk Fourier            
                AA(mkj,2,2) = AA(mkj,2,2)
     &               - TransFaceFracFourier(k,ik,jk)
     &               * Condf_th
           
c           eq incik col incjk Darcy 
                mij = NumCSRI(incjk)
                AA(mij,1,1) = AA(mij,1,1)
     &               + TransFaceFrac(k,ik,jk)*Condf*rho0

c           eq incik col incjk Fourier 
                mij = NumCSRI(incjk)
                AA(mij,2,2) = AA(mij,2,2)
     &               + TransFaceFracFourier(k,ik,jk)
     &               * Condf_th            

        
                if ( (IndIncNeuT(incik).eq.0) 
     &               .and.(NbFaceFracbyArete(na).ge.3) ) then ! >= 3 faces et pas sur le bord, on rajoute la conv thermique

                   ! eq incik col inck convection thermique 
                   mik = NumCSRI(inck)
                   AA(mik,2,1) = AA(mik,2,1)
     &                  - TransFaceFrac(k,ik,jk)
     &                  * Condf * rhof * Hf              
                   
                   ! eq incik col incjk convection thermique 
                   mij = NumCSRI(incjk)
                   AA(mij,2,1) = AA(mij,2,1) + TransFaceFrac(k,ik,jk)
     &                  * Condf * rhof * Hf
              
                endif

           
             enddo              ! fin boucle jk 


             ! eq inck col incik Darcy, rhof 
             mki = NumCSRK(incik)
             AA(mki,1,:) = AA(mki,1,:) + fik * Drhof(:)

             ! eq inck col incik convection thermique rhof             
             AA(mki,2,:) = AA(mki,2,:) + fik * Drhof(:) * Hf ! eq inck energie, der / P,T ds rhof

             ! eq inck col incup, convection thermique 
             mkup = NumCSRK(incup)
             AA(mkup,2,:) = AA(mkup,2,:) + fik * rhof * DHf(:)! eq inck energie, der / H incup            

    
             if ( (IndIncNeuT(incik).eq.0)   
     &            .and.(NbFaceFracbyArete(na).ge.3) ) then ! >= 3 faces et pas sur le bord, on rajoute la conv thermique

                ! eq incik col incik convection thermique rhof 
                mii = NumCSRI(incik)
                AA(mii,2,:) = AA(mii,2,:)
     &               - fik * Drhof(:) * Hf
                
                ! eq incik col incup convection thermique, der / solT(incup)  
                miup = NumCSRI(incup)
                AA(miup,2,2) = AA(miup,2,2) - fik * rhof 
              
             endif
          
          enddo ! fin boucle ik 
       enddo   ! fin boucle faces frac k 
c
c
c       
ccccccccccccccccccccccccccccccccc
c
c     Assemblage des flux d'interface mf 
c
       do nf = 1,NbFaceFrac
          
          nface = NumFaceFracVersFace(nf)
          k1 = NumCellbyFace(nface,1)
          k2 = NumCellbyFace(nface,2)
          inck1 = NumIncCell(k1)
          inck2 = NumIncCell(k2)
          
c
c     flux f-ik1 
c          
          incf = NumIncFaceFrac(nf)
          incik = NumIncCellbyFaceFrac(nf,1)
          
          
          Tfik = TransCellbyFaceFrac(nf,1)
          Tgik = TransCellbyFaceFracFourier(nf,1)
c
c     Assemblage du flux TPFA 
c          
          fik = Tfik*( SolP(incf) - SolP(incik)
     &         + grav * Rho_nm1(incik)
     &         * ( XInc(incf,3) - XInc(incik,3) )  )

          gik = Tgik*( SolT(incf) - SolT(incik))


          if (fik.ge.0.d0) then

             incup = incf 
             
          else

             incup = inck1 
             
          endif
          

          Hf = Enth(incup) 
          DHf(:) = DerEnth(incup,:) ! derivee de l'enthalpie incup 
          
          ! on centre la densite sur l'interface incik 
          rhof = Rho(incik)
          Drhof(:) = DerRho(incik,:)                            

          
          Sm(incf,1) = Sm(incf,1)   + fik * rhof ! flux darcy avec rhof  
          Sm(incik,1) = Sm(incik,1) - fik * rho0 ! flux darcy avec scaling rho0 

          Sm(incf,2) = Sm(incf,2)   + gik ! flux fourier + convection
     &         + fik * rhof * Hf 

          
          Sm(incik,2) = Sm(incik,2) - gik ! flux fourier seul  
          
          
          
          do m = IndLigneAA(incf)+1,IndLigneAA(incf+1)
             NumCSRK(NColAA(m)) = m
          enddo          

           do m = IndLigneAA(incik)+1,IndLigneAA(incik+1)
             NumCSRI(NColAA(m)) = m
           enddo


           mff = NumCSRK(incf)
           mfi = NumCSRK(incik)
           
           mii = NumCSRI(incik)
           mif = NumCSRI(incf)

           AA(mff,1,1) = AA(mff,1,1) + Tfik*rhof 
           AA(mfi,1,1) = AA(mfi,1,1) - Tfik*rhof 

           AA(mii,1,1) = AA(mii,1,1) + Tfik*rho0  
           AA(mif,1,1) = AA(mif,1,1) - Tfik*rho0 

           AA(mfi,1,:) = AA(mfi,1,:) + fik * Drhof(:)
           

           AA(mff,2,1) = AA(mff,2,1) + Tfik*rhof * Hf ! eq incf energie, der / Pincf ds fik 
           AA(mfi,2,1) = AA(mfi,2,1) - Tfik*rhof * Hf ! eq incf energie, der / Pincik ds fik 
           
           AA(mff,2,2) = AA(mff,2,2) + Tgik 
           AA(mfi,2,2) = AA(mfi,2,2) - Tgik 

           AA(mii,2,2) = AA(mii,2,2) + Tgik  
           AA(mif,2,2) = AA(mif,2,2) - Tgik

           AA(mfi,2,:) = AA(mfi,2,:)
     &          + fik * Drhof(:) * Hf ! eq incf energie, der / P,T ds rhof 

           
           mfup = NumCSRK(incup)           
           AA(mfup,2,:) = AA(mfup,2,:) + fik * rhof * DHf(:) ! eq incf energie, der / Tincup 
           
c
c     flux f-ik2 
c          
          incf = NumIncFaceFrac(nf)
          incik = NumIncCellbyFaceFrac(nf,2)
          Tfik = TransCellbyFaceFrac(nf,2)
          Tgik = TransCellbyFaceFracFourier(nf,2)
c
c     Assemblage du flux TPFA 
c 
          fik = Tfik*( SolP(incf) - SolP(incik)
     &         + grav * Rho_nm1(incik)
     &           * ( XInc(incf,3) - XInc(incik,3) )  )

          gik = Tgik*( SolT(incf) - SolT(incik))

          
          if (fik.ge.0.d0) then

             incup = incf 
             
          else

             incup = inck2 
             
          endif


          Hf = Enth(incup) 
          DHf(:) = DerEnth(incup,:) ! derivee de l'enthalpie incup           

          ! on centre la densite sur l'interface incik 
          rhof = Rho(incik)
          Drhof(:) = DerRho(incik,:)
         
          
          
          Sm(incf,1) = Sm(incf,1)   + fik * rhof ! flux darcy avec rhof        
          Sm(incik,1) = Sm(incik,1) - fik * rho0 ! flux darcy avec scaling rho0 

          Sm(incf,2) = Sm(incf,2)   + gik ! flux fourier + convection
     &         + fik * rhof * Hf

          
          Sm(incik,2) = Sm(incik,2) - gik ! flux fourier seul 

          

          
          do m = IndLigneAA(incf)+1,IndLigneAA(incf+1)
             NumCSRK(NColAA(m)) = m
          enddo          

           do m = IndLigneAA(incik)+1,IndLigneAA(incik+1)
             NumCSRI(NColAA(m)) = m
           enddo


           mff = NumCSRK(incf)
           mfi = NumCSRK(incik)
           
           mii = NumCSRI(incik)
           mif = NumCSRI(incf)
           

           AA(mff,1,1) = AA(mff,1,1) + Tfik*rhof 
           AA(mfi,1,1) = AA(mfi,1,1) - Tfik*rhof 

           AA(mii,1,1) = AA(mii,1,1) + Tfik*rho0  
           AA(mif,1,1) = AA(mif,1,1) - Tfik*rho0 

           AA(mfi,1,:) = AA(mfi,1,:) + fik * Drhof(:)


           AA(mff,2,1) = AA(mff,2,1) + Tfik*rhof * Hf ! eq incf energie, der / Pincf ds fik 
           AA(mfi,2,1) = AA(mfi,2,1) - Tfik*rhof * Hf ! eq incf energie, der / Pincik ds fik 

           
           AA(mff,2,2) = AA(mff,2,2) + Tgik 
           AA(mfi,2,2) = AA(mfi,2,2) - Tgik 

           AA(mii,2,2) = AA(mii,2,2) + Tgik  
           AA(mif,2,2) = AA(mif,2,2) - Tgik

           AA(mfi,2,:) = AA(mfi,2,:)
     &          + fik * Drhof(:) * Hf ! eq incf energie, der / P,T ds rhof 

           mfup = NumCSRK(incup)           
           AA(mfup,2,:) = AA(mfup,2,:) + fik * rhof * DHf(:) ! eq incf energie, der / Tincup 
           
           
           
       enddo


       
ccccccccccccccccccccccccccccccccc        
c
c     Equations et Sm Interfaces Dirichlet 
c
       do i=1,NbCV
          
          if (IndIncDir(i).eq.1) then 

             
             m = IndLigneAA(i)+1
             AA(m:IndLigneAA(i+1),1,1) = 0.d0 
             AA(m,1,1) = AA(m,1,1) + 1.d0
             Sm(i,1) = SolP(i) - SolPDir(i)
            
          endif


          if (IndIncDirT(i).eq.1) then 

             m = IndLigneAA(i)+1
             AA(m:IndLigneAA(i+1),2,2) = 0.d0              
             AA(m,2,2) = AA(m,2,2) + 1.d0

             
             if (IndIncNeu(i).eq.2) then
                
                Sm(i,2) = SolT(i) - 285.d0
                
             else

                Sm(i,2) = SolT(i) - SolTDir(i)
                
             endif
             
                          
          endif
        
       enddo
c
c
ccccccccccccccccccccccccccccccccccc
c
c     condition de Neumann non homogene (IndIncNeu(inc) = 2 )
c
        surfT = 50.d0*100.d0 ! le facteur 50 ne sert pas car il est compense par le volT !!!!!!!!!
       
        volT = surfT*2000.d0
       
        Temps = f_temps_debit()
        
       
       do i=1,NbFace
          surf = SurfaceFace(i)
          inc = NumIncFace(i)

          flux = - (surf/surfT)*volT/Temps
          
          
          if (IndIncNeu(inc).eq.2) then

             
                         
             Sm(inc,1) = Sm(inc,1) + flux*rho0             
!            correspond au debit = 2000.d0*rho0/Temps  en Kg /m^2/s
                          
          endif

c          if (IndIncNeuT(inc).eq.1) then
c
c             if (XInc(inc,1).le.1.d-3) then 
c                Sm(inc,2) = Sm(inc,2) - flux*rho0*300.d0 
c             endif
c          endif
                    
       enddo

cccccccccccccccccccccccccccc
c
c     Second membre = - residu 
c
       Sm(:,:) = - Sm(:,:) 
       
       
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


       deallocate(Rho,DerRho)
       deallocate(EInt,DerEInt)
       deallocate(Enth,DerEnth)
       deallocate(EInt_prev)
       deallocate(Rho_nm1)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       return
       end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
