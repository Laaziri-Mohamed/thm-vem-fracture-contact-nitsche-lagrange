cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Calcul du traceur au pas de temps n                                                          
c
c     Inconnues: traceur aux mailles et aux faces frac 
c
c     On suppose qu'il n'y a pas d'intersection de fractures car sinon il faudrait rajouter
c     les inconnues aux intersections avec au moins 3 faces frac 
c        
ccccccccccccccc
c      
c     En entree: SolP, AccP, AccP_nm1 : on recalcule les flux de darcy qui sont conservatifs 
c
c     En entree: SolC_nm1 
c
cccccccccc
c      
c     En sortie: SolC 
c
c
cccccccccc
c    
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
      subroutine JacSmTraceur(
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
      dimension NbAretebyFace(NbFace)
      dimension NumAretebyFace(NbFace,NbAretebyFaceMax)
        
      dimension NumFacebyCell(NbCell,NbFacebyCellMax)
      dimension NumCellbyFace(NbFace,2)

      dimension IndFaceFrac(NbFace)
      
c
c     num face frac > num face et en sens inverse
c      
      dimension NumFaceFracVersFace(NbFaceFrac)

      dimension NumFaceVersFaceFrac(NbFace)
      
ccccccccccccccccccccc
c        
c     Pression P et presion init P0 
c
        dimension SolP(NbCV)

        dimension SolT(NbCV)
        
        dimension SolP0(NbCV)                
c
c     Acc Darcy n et nm1 
c      
        dimension AccP(NbCV)       
        dimension AccP_nm1(NbCV)       
c
c     epaisseur face frac 
c
        dimension dfbyFaceFrac(NbFaceFrac)
        dimension dfbyFaceFrac_nm1(NbFaceFrac)
        
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
c
c     num inc > indice 0 ou 1 (neu homog) ou 2 (neu non homog) 
c
        dimension IndIncNeu(NbCV)                
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

c
c     Pour HFV discontinu: TransCellbyFaceFrac : TKsigma et TLsigma 
c
c     num des cells K et L donnes par NumCellbyFace(numface,1:2) 
c
c        
        dimension TransCellbyFaceFrac(NbFaceFrac,2)       
c
c     en entree 
c
c     CL Dirichlet aux inconnues interfaces de bord DIR 
c
       dimension SolPDir(NbCV)
c
c
ccccccccccc       
c
c     traceur au pas de temps precedent 
c       
       dimension SolC_nm1(NbCell + NbFaceFrac)

       
cccccccccccccccccccccccccccccccccccccccc
c
c     SORTIES 
c
cccccccccccccccccccccccccccccccccccccccc
c
c     traceur au nouveau pas de temps 
c       
       dimension SolC(NbCell + NbFaceFrac)

       
c
ccccccccccccccccccccccccccccccccccccccc
c
c     WorkSpaces 
c
c
c
c     
c     Jacobienne CSR mailles + noeuds 
c     On met l'elt diagonal en premier 

      integer, dimension(:), allocatable :: IndLigneAA
      integer, dimension(:), allocatable :: NLigneAA
      integer, dimension(:), allocatable :: NColAA

      double precision, dimension(:), allocatable :: AA
      
      integer, dimension(:), allocatable :: NbNzbyLine
      integer, dimension(:,:), allocatable :: NumNzbyLine

c     
c     Indices m CSR des inc colonnes des lignes K et I ou Sig et I    
c     
      integer, dimension(:), allocatable :: NumCSRK
      integer, dimension(:), allocatable :: NumCSRI

      
      integer, dimension(:), allocatable :: NbFaceFracbyArete
      integer, dimension(:,:), allocatable :: NumFaceFracbyArete


      

      double precision, dimension(:), allocatable :: Sm      

c     
c     increment 
c
      double precision, dimension(:), allocatable :: dSolC     
      
       dimension X(NbDim)
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     face frac voisine d'une arete 
c
       allocate(NbFaceFracbyArete(NbArete))
       allocate(NumFaceFracbyArete(NbArete,NbFacebyAreteMax))

       NbFaceFracbyArete(:) = 0
       NumFaceFracbyArete(:,:) = 0
       
       do k=1,NbFaceFrac

          nf = NumFaceFracVersFace(k)

          do ia = 1,NbAretebyFace(nf)
             na = NumAretebyFace(nf,ia)

             NbFaceFracbyArete(na) = NbFaceFracbyArete(na) + 1
             if (NbFaceFracbyArete(na).gt.NbFacebyAreteMax) then
                write(*,*)' trop de face frac par arete redim ',
     &               NbFaceFracbyArete(na),NbFacebyAreteMax
                stop
             endif

             NumFaceFracbyArete(na,NbFaceFracbyArete(na)) = k 
             
             
          enddo
          
       enddo



c       do na=1,NbArete
c          if ( NbFaceFracbyArete(na).ne.0) then
c             write(*,*)NbFaceFracbyArete(na),na,
c     &            NumFaceFracbyArete(na,1:NbFaceFracbyArete(na))             
c          endif          
c       enddo
c       stop
ccccccccccccccccccccccccccccccccccccc
c
c     structure creuse 
c       
       NCVT = NbCell + NbFaceFrac

       write(*,*)' NCVT ',NCVT 
       
       allocate(NbNzbyLine(NCVT))
       allocate(NumNzbyLine(NCVT,NbNzbyLineMax))


       NbNzbyLine = 0
       NumNzbyLine = 0
c
c     Boucle mailles 
c
       do k=1,NbCell

          incc = k

c         rajout incc eq ligne incc
          NbNzbyLine(incc) = NbNzbyLine(incc) + 1
          if (NbNzbyLine(incc).gt.NbNzbyLineMax) then 
           write(*,*)'NbNzbyLine > NbNzbyLineMax ',incc,
     &            NbNzbyLine(incc),NbNzbyLineMax
           stop
          endif
          NumNzbyLine(incc,NbNzbyLine(incc)) =  incc

c             write(*,*) 'nz ',incc,NbNzbyLine(incc),
c     &             NumNzbyLine(incc,NbNzbyLine(incc))           
c
c     boucle interfaces ik
c
          nik = NbInterfacebyCell(k)
          do ik = 1,nik
             nf = NumFacebyCell(k,ik) ! numero de la face ik 
             

             if (IndFaceFrac(nf).eq.1) then
                ifrac = NumFaceVersFaceFrac(nf)
                   
                incup = ifrac + NbCell ! inc frac traceur 


c         rajout incup eq ligne incc
                NbNzbyLine(incc) = NbNzbyLine(incc) + 1
                if (NbNzbyLine(incc).gt.NbNzbyLineMax) then 
                   write(*,*)'NbNzbyLine > NbNzbyLineMax ',incc,
     &                  NbNzbyLine(incc),NbNzbyLineMax
                   stop
                endif
                NumNzbyLine(incc,NbNzbyLine(incc)) =  incup   
                

                incf = ifrac + NbCell ! on rajoute la maille k (incc) sur l 'eq frac incf 

c         rajout incup eq ligne incc
                NbNzbyLine(incf) = NbNzbyLine(incf) + 1
                if (NbNzbyLine(incf).gt.NbNzbyLineMax) then 
                   write(*,*)'NbNzbyLine > NbNzbyLineMax ',incf,
     &                  NbNzbyLine(incf),NbNzbyLineMax
                   stop
                endif
                NumNzbyLine(incf,NbNzbyLine(incf)) =  incc                    

                
             else if (NumCellbyFace(nf,2).ge.1) then 
                   
                if (NumCellbyFace(nf,1).eq.k) then
                   incup = NumCellbyFace(nf,2)
                else
                   incup = NumCellbyFace(nf,1)
                endif

c         rajout incup eq ligne incc
                NbNzbyLine(incc) = NbNzbyLine(incc) + 1
                if (NbNzbyLine(incc).gt.NbNzbyLineMax) then 
                   write(*,*)'NbNzbyLine > NbNzbyLineMax ',incc,
     &                  NbNzbyLine(incc),NbNzbyLineMax
                   stop
                endif
                NumNzbyLine(incc,NbNzbyLine(incc)) =  incup   
                

             endif


c             write(*,*) 'nz ',incc,NbNzbyLine(incc),
c     &             NumNzbyLine(incc,NbNzbyLine(incc)) 
             
          enddo
       enddo
c
c     boucle face frac 
c       
       do k=1,NbFaceFrac
          nik = NbInterfacebyFaceFrac(k) 
          nf = NumFaceFracVersFace(k)
          
          incc = k + NbCell ! inc face frac traceur 

c         rajout incc eq ligne incc
          NbNzbyLine(incc) = NbNzbyLine(incc) + 1
          if (NbNzbyLine(incc).gt.NbNzbyLineMax) then 
           write(*,*)'NbNzbyLine > NbNzbyLineMax ',incc,
     &            NbNzbyLine(incc),NbNzbyLineMax
           stop
          endif
          NumNzbyLine(incc,NbNzbyLine(incc)) =  incc            

          
          do ik = 1,nik
             na = NumAretebyFace(nf,ik)


             if (NbFaceFracbyArete(na).eq.1) then
                
! on suppose flux nul ici donc rien a assembler
                   
             else if (NbFaceFracbyArete(na).eq.2) then
                if (NumFaceFracbyArete(na,1).eq.k) then
                   incup = NumFaceFracbyArete(na,2) + NbCell
                else
                   incup = NumFaceFracbyArete(na,1) + NbCell 
                endif

c                write(*,*) 'incup frac ',incup - NbCell
                   

c         rajout incup eq ligne incc
                NbNzbyLine(incc) = NbNzbyLine(incc) + 1
                if (NbNzbyLine(incc).gt.NbNzbyLineMax) then 
                   write(*,*)'NbNzbyLine > NbNzbyLineMax ',incc,
     &                  NbNzbyLine(incc),NbNzbyLineMax
                   stop
                endif
                NumNzbyLine(incc,NbNzbyLine(incc)) =  incup      

                   
             else                   
                write(*,*)' on suppose 2 face frac vois max ',
     &               NbFaceFracbyArete(na)
                stop
             endif


             
             
          enddo
       enddo
c
c
c
c     
c     IndLigneAA, NColAA, NLigneAA  
c     

       allocate(IndLigneAA(NCVT+1))
       IndLigneAA(1) = 0
       do i=1,NCVT
          IndLigneAA(i+1) = IndLigneAA(i) + NbNzbyLine(i)          
       enddo
       MemBloc = IndLigneAA(NCVT+1)

       write(*,*)' MemBloc ',MemBloc 
       
       allocate(NLigneAA(MemBloc))
       allocate(NColAA(MemBloc))
       do i=1,NCVT
          do m=IndLigneAA(i)+1,IndLigneAA(i+1)
             n = m-IndLigneAA(i)
             NColAA(m) = NumNzbyLine(i,n)
             NLigneAA(m) = i
          enddo
       enddo
c     
c     Permutation pour mettre l'elt diagonal en premier sur la ligne 
c     
       do i=1,NCVT
          md = IndLigneAA(i)+1
          jd = NColAA(md)
          do m = md,IndLigneAA(i+1)
             j = NColAA(m)
             if (j.eq.i) then 
                m1 = m
             endif
          enddo
          NColAA(md) = i
          NColAA(m1) = jd
       enddo
c
ccccccccccccccccc
c
c       do k=1,NCVT
c       do k=1,NbCell   
c          do m=IndLigneAA(k)+1,IndLigneAA(k+1)
c             if (NColAA(m).gt.NbCell) then 
c                write(*,*)' ligne k ',k,m,NColAA(m)
c             endif
c          enddo
c          write(*,*)
c       enddo
c
c       stop
c
ccccccccccccccc


       
       
       allocate(NumCSRK(MemBloc))
       allocate(NumCSRI(MemBloc))       
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
ccccccccccccccccccccccccccc

       allocate(AA(MemBloc))
       allocate(Sm(NCVT))
       
       
c
c     mise a zero de AA et Sm 
c
 
       AA(:) = 0.d0
 
       Sm(:) = 0.d0

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Assemblage des termes d'accumulation aux mailles et aux faces fracture 
c
ccccccccccccccccccccccccccccccccccccc       
c
c
       
       do k=1,NbCell
          inck = NumIncCell(k) ! inc maille pression 

          incc = k ! inc maille traceur 
         
          Sm(incc) = Sm(incc) + ( AccP(inck)*SolC(incc)
     &         - AccP_nm1(inck)*SolC_nm1(incc) )/Deltat 
   
          
          mkk = IndLigneAA(incc) + 1

          AA(mkk) = AA(mkk) + AccP(inck)/Deltat 

                    
       enddo


       do k=1,NbFaceFrac
          inck = NumIncFaceFrac(k) ! inc face frac pression 

          incc = k + NbCell ! inc face frac traceur 

          Sm(incc) = Sm(incc) + ( AccP(inck)*SolC(incc)
     &         - AccP_nm1(inck)*SolC_nm1(incc) )/Deltat 
   
          
          mkk = IndLigneAA(incc) + 1

          AA(mkk) = AA(mkk) + AccP(inck)/Deltat 
          
       enddo       



c       write(*,*)' fin acc '
       
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Assemblage des flux 
c
cccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccc
c       
c     boucle sur les mailles 
c
ccccccccccccccccccccccccccccc
c       
       do k=1,NbCell
          nik = NbInterfacebyCell(k) 
          inck = NumIncCell(k) ! num inc pression de la maille k 

          incc = k ! num inc traceur maille k 
c          
c         on stocke les indices CSR des ind de col de la ligne incc  
c
          do m = IndLigneAA(incc)+1,IndLigneAA(incc+1)
             NumCSRK(NColAA(m)) = m
          enddo

          do ik = 1,nik
             incik = NumInterfacebyCell(k,ik) ! num inc pression de la face ik 
             nf = NumFacebyCell(k,ik) ! numero de la face ik 

c
c     Calcul du flux Darcy
c
             fik = 0.d0 
             do jk = 1,nik
                incjk = NumInterfacebyCell(k,jk)
                
                fik = fik + TransCell(k,ik,jk)
     &               *( SolP(inck) - SolP(incjk)
     &                + grav*rho0*( XInc(inck,3) - XInc(incjk,3) )  )
                
             enddo

             pik = SolP(incik)
             Tik = SolT(incik)             
             rhof = f_rho(pik,Tik)

             fik = fik*rhof ! flux de Darcy 

             if (fik.ge.0.d0) then 

                incup = incc
                Sm(incc) = Sm(incc) + SolC(incup)*fik ! flux sur inc maille incc 

                mkup = NumCSRK(incup) 
                AA(mkup) = AA(mkup) + fik                  

                if (IndFaceFrac(nf).eq.1) then
 
                   ifrac = NumFaceVersFaceFrac(nf)                   
                   incf = ifrac + NbCell ! inc frac traceur                    

                   do m = IndLigneAA(incf)+1,IndLigneAA(incf+1)
                      NumCSRI(NColAA(m)) = m
                   enddo                   

                   Sm(incf) = Sm(incf) - SolC(incup)*fik ! si face frac on rajoute -flux sur inc frac incf 
                   
                   mfup = NumCSRI(incup) 
                   AA(mfup) = AA(mfup) - fik      
                   
                   
                endif
                
                
             else

                if (IndIncNeu(incik).eq.2) then

                   Sm(incc) = Sm(incc) + fik ! bord entrant on injecte C=1 
                   
                else if ( (IndIncNeu(incik).eq.1)
     &                  .or.(IndIncDir(incik).eq.1) ) then 

                      ! flux darcy nul ou bord sortant (on injecte 0 sinon) 
                      
                else if (IndFaceFrac(nf).eq.1) then
                   ifrac = NumFaceVersFaceFrac(nf)
                   
                   incup = ifrac + NbCell ! inc face frac traceur inc amont ici 
                   incf  = ifrac + NbCell ! inc face frac traceur 
                   
                   Sm(incc) = Sm(incc) + fik*SolC(incup) ! flux sur inc maille incc 
                   
                   mkup = NumCSRK(incup) 
                   AA(mkup) = AA(mkup) + fik

                   do m = IndLigneAA(incf)+1,IndLigneAA(incf+1)
                      NumCSRI(NColAA(m)) = m
                   enddo                   

                   Sm(incf) = Sm(incf) - SolC(incup)*fik ! on rajoute -flux sur inc frac incf 
                   
                   mfup = NumCSRI(incup) 
                   AA(mfup) = AA(mfup) - fik    
                   
                   
                else

                   if (NumCellbyFace(nf,1).eq.k) then
                      incup = NumCellbyFace(nf,2)
                   else
                      incup = NumCellbyFace(nf,1)
                   endif

                   Sm(incc) = Sm(incc) + fik*SolC(incup)
                   
                   mkup = NumCSRK(incup) 
                   AA(mkup) = AA(mkup) + fik      
                   
                   
                endif
                      
                
             endif
             
                   
          enddo
       enddo   


c       write(*,*)' fin boucle maille '
c
cccccccccccccccccccccccccccccccccccccc
c
c       
c     boucle sur les faces frac 
c       
cccccccccccccccccccccccccccccccccccccc
c       
       do k=1,NbFaceFrac
          nik = NbInterfacebyFaceFrac(k) 
          inck = NumIncFaceFrac(k) ! inc face frac pression 

          nf = NumFaceFracVersFace(k)
          
          incc = k + NbCell ! inc face frac traceur 
          
c          Condf = dfbyFaceFrac(k)**3 ! conductivite frac
        Condf = dfbyFaceFrac_nm1(k)**3 ! dependence en df**3 explicite         
          
c          
c         on stocke les indices CSR des ind de col de la ligne incc  
c
          do m = IndLigneAA(incc)+1,IndLigneAA(incc+1)
             NumCSRK(NColAA(m)) = m
          enddo

          do ik = 1,nik
             incik = NumInterfacebyFaceFrac(k,ik)


             na = NumAretebyFace(nf,ik)
c
c     Calcul du flux de Darcy 
c
             fik = 0.d0 
             do jk = 1,nik
                incjk = NumInterfacebyFaceFrac(k,jk)
                
                fik = fik + TransFaceFrac(k,ik,jk)*Condf 
     &               *( SolP(inck)  - SolP(incjk)
     &                + grav*rho0*( XInc(inck,3) - XInc(incjk,3) )  )
                
                
             enddo

             pik = SolP(incik)
             Tik = SolT(incik)             
             rhof = f_rho(pik,Tik)             
                          
             fik = fik*rhof  ! flux de Darcy         


             if (fik.ge.0.d0) then

                incup = k + NbCell

                Sm(incc) = Sm(incc) + fik*SolC(incup)

                mkup = NumCSRK(incup) 
                AA(mkup) = AA(mkup) + fik       

             else
                
                if (NbFaceFracbyArete(na).eq.1) then

                ! on suppose flux nul ici donc rien a assembler 
                
                else if (NbFaceFracbyArete(na).eq.2) then

                   if (NumFaceFracbyArete(na,1).eq.k) then
                      incup = NumFaceFracbyArete(na,2) + NbCell
                   else
                      incup = NumFaceFracbyArete(na,1) + NbCell 
                   endif

                   Sm(incc) = Sm(incc) + fik*SolC(incup)

                   mkup = NumCSRK(incup) 
                   AA(mkup) = AA(mkup) + fik       
                   
                else
                   
                   write(*,*)' on suppose 2 face frac vois max ',
     &                  NbFaceFracbyArete(na)
                   stop
                endif
                
             endif


             
          enddo
       enddo   
c
c

c       write(*,*)' fin boucle face frac '
       
cccccccccccccccccccccccccc
c
c       do k=1,NCVT
c          write(*,*)' Sm ',k,Sm(k)
c       enddo
c
c       
cccccccccccccccccccccccccccc
c
c     Second membre = - residu 
c
       Sm(:) = - Sm(:)


       allocate(dSolC(NCVT))

       dSolC = 0.d0 

       NbIncMeca = 1
       ITER = 0
       call solveurSuperLU(NCVT,MemBloc,NbIncMeca,
     &      AA,NLigneAA,NColAA,IndLigneAA,
     &      Sm,dSolC)      

c      critere_arret_gmres = 1.0d-8
c      call solveurILU0(
c     &       NbCVT,MemBloc,
c     &       AA,NLigneAA,NColAA,IndLigneAA,
c     &       Sm,dSolC,critere_arret_gmres,
c     &       ITER,IERR,ERR)
       
       SolC = SolC + dSolC 

       
       deallocate(dSolC)

cccccccccccccccccccccccccc
c
c       do k=1,NCVT
c          write(*,*)' SolC ',k,SolC(k)
c       enddo
c
c     
c       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       

       deallocate(AA,IndLigneAA,NColAA,NLigneAA,Sm)
       deallocate(NumCSRK,NumCSRI,NbNzbyLine,NumNzbyLine)

       deallocate(NbFaceFracbyArete)
       deallocate(NumFaceFracbyArete)       

c       stop
       
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       return
       end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
