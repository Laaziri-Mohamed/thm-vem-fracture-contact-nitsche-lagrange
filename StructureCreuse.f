cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Structure creuse CSR: IndLigneAA, NColAA, NLigneAA, MemBloc 
c
c      Ligne i: m=IndLigneAA(i)+1, ..., IndLigneAA(i+1)
c
c
c     On permute a la fin pour mettre l'elt diagonal en premier sur la ligne
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
      subroutine StructureCreuse(
     &   NbCell,NbFace,NbArete,NbFaceFrac,NbCV,
     &   NumFacebyCell,NumCellbyFace,
     &   NumFaceVersFaceFrac,IndFaceFrac,  
     &   NumIncCell,NumIncFaceFrac, 
     &   NbInterfacebyCell,NumInterfacebyCell,
     &   NbInterfacebyFaceFrac,NumInterfacebyFaceFrac,
     &   NumIncCellbyFaceFrac,
     &   NbFaceFracbyArete,NumFaceFracbyArete,
     &   NumFaceFracVersFace,NumAretebyFace,
     &   NbNzbyLine,NumNzbyLine)
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

        dimension NumFaceFracVersFace(NbFaceFrac)
        dimension NumAretebyFace(NbFace,NbAretebyFaceMax)        
c
c     Numero des inconnues mailles et faces fractures 
c       
c     num cell       >  num inc cell
c     num face frac  >  num inc face frac  
c
        dimension NumIncCell(NbCell)
        dimension NumIncFaceFrac(NbFaceFrac)
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
c     
c     pour HFV discontinu 
c     num face frac  > num inc Ksigma et Lsigma 
c     
        dimension NumIncCellbyFaceFrac(NbFaceFrac,2)    
c          
c
cccccccccccccccccccccccccccccccccccccccc
c
c     SORTIES 
c
c
        dimension NbNzbyLine(NbCV)
        dimension NumNzbyLine(NbCV,NbNzbyLineMax)
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     Comptage pour le stockage CSR 
c
       do i=1,NbCV
         NbNzbyLine(i) = 0
       enddo


cccccccccccccccccccccccccccc
c
c     Boucle maille k \in M 
c
       do k=1,NbCell

          inck = NumIncCell(k)
c          write(*,*)' k inc k ',k,inck

c         rajout inc Uk (inck) eq ligne inck
          NbNzbyLine(inck) = NbNzbyLine(inck) + 1
          if (NbNzbyLine(inck).gt.NbNzbyLineMax) then 
           write(*,*)'NbNzbyLine > NbNzbyLineMax ',
     &                NbNzbyLine(inck),NbNzbyLineMax
          endif
          NumNzbyLine(inck,NbNzbyLine(inck)) =  inck  
c
c     boucle interfaces ik \in IK
c
          nik = NbInterfacebyCell(k)
          do ik = 1,nik 

          incik = NumInterfacebyCell(k,ik)
c          write(*,*)' incik ',incik 
          
c         terme FKi = \sum_{j\in IK} (Uk-Uj)  
c
c      
c          rajout inc Uk (inck) eq ligne incik (pas besoin de test )
           NbNzbyLine(incik) = NbNzbyLine(incik) + 1  
           if (NbNzbyLine(incik).gt.NbNzbyLineMax) then 
              write(*,*)'NbNzbyLine > NbNzbyLineMax ',
     &             NbNzbyLine(incik),NbNzbyLineMax
           endif                      
           NumNzbyLine(incik,NbNzbyLine(incik)) =  inck  

          do jk = 1,nik 
             incjk = NumInterfacebyCell(k,jk)
c             write(*,*)' incjk ',incjk 
          
c            test ajout inc Uj (incjk) eq ligne inck    

             iajout = 1
             do m=1,NbNzbyLine(inck)
                if (NumNzbyLine(inck,m).eq.incjk) then 
                   iajout = 0
                endif
             enddo
             if (iajout.eq.1) then 
              NbNzbyLine(inck) = NbNzbyLine(inck) + 1      
              NumNzbyLine(inck,NbNzbyLine(inck)) =  incjk  
              if (NbNzbyLine(inck).gt.NbNzbyLineMax) then 
               write(*,*)' NbNzbyLine > NbNzbyLineMax ',
     &                   NbNzbyLine(inck),NbNzbyLineMax
               stop
              endif                   
             endif           

c            test ajout inc Uj (incjk) eq ligne incik  

             iajout = 1
             do m=1,NbNzbyLine(incik)
                if (NumNzbyLine(incik,m).eq.incjk) then 
                   iajout = 0
                endif
             enddo
             if (iajout.eq.1) then 
              NbNzbyLine(incik) = NbNzbyLine(incik) + 1      
              NumNzbyLine(incik,NbNzbyLine(incik)) =  incjk  
              if (NbNzbyLine(incik).gt.NbNzbyLineMax) then 
               write(*,*)' NbNzbyLine > NbNzbyLineMax ',
     &                   NbNzbyLine(incik),NbNzbyLineMax
               stop
              endif                   
             endif 


          enddo

c
c     on rajoute a l'eq inck l'inc de la maille voisine de k ou l'inc face frac si ik = face frac 
c
          nf = NumFacebyCell(k,ik) 
          if (NumCellbyFace(nf,2).ge.1) then

             
             if (NumCellbyFace(nf,1).eq.k) then
                kvois = NumCellbyFace(nf,2)
             else
                kvois = NumCellbyFace(nf,1)
             endif

c             write(*,*)' kvois ',k,ik,kvois,nf, NumCellbyFace(nf,:)

             if ( IndFaceFrac(nf).ne.1 ) then
                
                NbNzbyLine(inck) = NbNzbyLine(inck) + 1
                if (NbNzbyLine(inck).gt.NbNzbyLineMax) then 
                   write(*,*)' NbNzbyLine > NbNzbyLineMax ',
     &                  NbNzbyLine(inck),NbNzbyLineMax
                   stop
                endif                                   
                NumNzbyLine(inck,NbNzbyLine(inck)) =  NumIncCell(kvois)  

             else ! face frac 

                ifrac = NumFaceVersFaceFrac(nf)
                incf = NumIncFaceFrac(ifrac)
                
                NbNzbyLine(inck) = NbNzbyLine(inck) + 1
                if (NbNzbyLine(inck).gt.NbNzbyLineMax) then 
                   write(*,*)' NbNzbyLine > NbNzbyLineMax ',
     &                  NbNzbyLine(inck),NbNzbyLineMax
                   stop
                endif                                   
                NumNzbyLine(inck,NbNzbyLine(inck)) =  incf  ! terme inck - incf 

                NbNzbyLine(incf) = NbNzbyLine(incf) + 1
                if (NbNzbyLine(incf).gt.NbNzbyLineMax) then 
                   write(*,*)' NbNzbyLine > NbNzbyLineMax ',
     &                  NbNzbyLine(incf),NbNzbyLineMax
                   stop
                endif                                   
                NumNzbyLine(incf,NbNzbyLine(incf)) =  inck  ! terme incf - inck 
                
             endif
                
          endif
             
c
c     fin boucle IK 
c
          enddo

c
c     fin boucle maille 
c
       enddo
c
cccccccccccccccccccccccccccc
c
c
c       write(*,*)
c       do i = 1,NbCV
c          write(*,*)' Line i ',i,NbNzbyLine(i)
c       enddo


c       write(*,*)
c       do i = 1,NbCV
c          write(*,*)' Line i ',i,(NumNzbyLine(i,m),m=1,NbNzbyLine(i))
c       enddo
c

c       stop
c
cccccccccccccccccccccccccccc
c
c     Boucle face frac k \in F_G 
c
       do k=1,NbFaceFrac

          inck = NumIncFaceFrac(k)

          nf = NumFaceFracVersFace(k)
          

c            test ajout inc Uk (inck) eq ligne inck    
             iajout = 1
             do m=1,NbNzbyLine(inck)
                if (NumNzbyLine(inck,m).eq.inck) then 
                   iajout = 0
                endif
             enddo
             if (iajout.eq.1) then 
              NbNzbyLine(inck) = NbNzbyLine(inck) + 1        
              NumNzbyLine(inck,NbNzbyLine(inck)) =  inck  
              if (NbNzbyLine(inck).gt.NbNzbyLineMax) then 
               write(*,*)' NbNzbyLine > NbNzbyLineMax ',
     &                   NbNzbyLine(inck),NbNzbyLineMax
               stop
              endif                   
             endif   
c
c     boucle interfaces ik \in Isig
c
          nik = NbInterfacebyFaceFrac(k)
          do ik = 1,nik

          na = NumAretebyFace(nf,ik)
             

          incik = NumInterfacebyFaceFrac(k,ik)

c         terme FKi = \sum_{j\in IK} (Uk-Uj)  
c
c            test ajout inc Uk (inck) eq ligne incik  
c

             iajout = 1
             do m=1,NbNzbyLine(incik)
                if (NumNzbyLine(incik,m).eq.inck) then 
                   iajout = 0
                endif
             enddo
             if (iajout.eq.1) then 
              NbNzbyLine(incik) = NbNzbyLine(incik) + 1        
              NumNzbyLine(incik,NbNzbyLine(incik)) =  inck  
              if (NbNzbyLine(incik).gt.NbNzbyLineMax) then 
               write(*,*)' NbNzbyLine > NbNzbyLineMax ',
     &                   NbNzbyLine(incik),NbNzbyLineMax
               stop
              endif                   
             endif 
             

          do jk = 1,nik 
             incjk = NumInterfacebyFaceFrac(k,jk)
          
c            test ajout inc Uj (incjk) eq ligne inck    

             iajout = 1
             do m=1,NbNzbyLine(inck)
                if (NumNzbyLine(inck,m).eq.incjk) then 
                   iajout = 0
                endif
             enddo
             if (iajout.eq.1) then 
              NbNzbyLine(inck) = NbNzbyLine(inck) + 1        
              NumNzbyLine(inck,NbNzbyLine(inck)) =  incjk  
              if (NbNzbyLine(inck).gt.NbNzbyLineMax) then 
               write(*,*)' NbNzbyLine > NbNzbyLineMax ',
     &                   NbNzbyLine(inck),NbNzbyLineMax
               stop
              endif                   
             endif           

c            test ajout inc Uj (incjk) eq ligne incik  


             iajout = 1
             do m=1,NbNzbyLine(incik)
                if (NumNzbyLine(incik,m).eq.incjk) then 
                   iajout = 0
                endif
             enddo
             if (iajout.eq.1) then 
              NbNzbyLine(incik) = NbNzbyLine(incik) + 1        
              NumNzbyLine(incik,NbNzbyLine(incik)) =  incjk  
              if (NbNzbyLine(incik).gt.NbNzbyLineMax) then 
               write(*,*)' NbNzbyLine > NbNzbyLineMax ',
     &                   NbNzbyLine(incik),NbNzbyLineMax
               stop
              endif                   
             endif 


          enddo

c
c     rajout de la face frac voisine de l'arete (si 2 faces frac voisines)  
c

          if (NbFaceFracbyArete(na).eq.2) then  

             if (NumFaceFracbyArete(na,1).eq.k) then
                kvois = NumFaceFracbyArete(na,2)
             else
                kvois = NumFaceFracbyArete(na,1) 
             endif

             inckvois = NumIncFaceFrac(kvois) 

             NbNzbyLine(inck) = NbNzbyLine(inck) + 1
             if (NbNzbyLine(inck).gt.NbNzbyLineMax) then 
                write(*,*)' NbNzbyLine > NbNzbyLineMax ',
     &               NbNzbyLine(inck),NbNzbyLineMax
                stop
             endif                                   
             NumNzbyLine(inck,NbNzbyLine(inck)) =  inckvois
             
          endif
c
c     fin boucle Isig
c
          enddo

c
c     fin boucle Face Frac 
c
       enddo
c
cccccccccccccccccc
ccccccccccccccccccc
c
c     Boucle sur les flux TPFA d'interface 
c     On rajoute les termes de couplage hors diagonaux incf <-> interface k1 et k2 
c       
       do nf = 1,NbFaceFrac

          incf = NumIncFaceFrac(nf)

          incik1 = NumIncCellbyFaceFrac(nf,1)
          
          incik2 = NumIncCellbyFaceFrac(nf,2)

          
          NbNzbyLine(incf) = NbNzbyLine(incf) + 1        
          if (NbNzbyLine(incf).gt.NbNzbyLineMax) then 
             write(*,*)' NbNzbyLine > NbNzbyLineMax ',
     &            NbNzbyLine(incf),NbNzbyLineMax
             stop
          endif                   
          NumNzbyLine(incf,NbNzbyLine(incf)) =  incik1  ! eq incf col ik1

          NbNzbyLine(incf) = NbNzbyLine(incf) + 1        
          if (NbNzbyLine(incf).gt.NbNzbyLineMax) then 
             write(*,*)' NbNzbyLine > NbNzbyLineMax ',
     &            NbNzbyLine(incf),NbNzbyLineMax
             stop
          endif                   
          NumNzbyLine(incf,NbNzbyLine(incf)) =  incik2  ! eq incf col ik2

          NbNzbyLine(incik1) = NbNzbyLine(incik1) + 1        
          if (NbNzbyLine(incik1).gt.NbNzbyLineMax) then 
             write(*,*)' NbNzbyLine > NbNzbyLineMax ',
     &            NbNzbyLine(incik1),NbNzbyLineMax
             stop
          endif                   
          NumNzbyLine(incik1,NbNzbyLine(incik1)) =  incf  ! eq ik1 col incf 

          NbNzbyLine(incik2) = NbNzbyLine(incik2) + 1        
          if (NbNzbyLine(incik2).gt.NbNzbyLineMax) then 
             write(*,*)' NbNzbyLine > NbNzbyLineMax ',
     &            NbNzbyLine(incik2),NbNzbyLineMax
             stop
          endif                   
          NumNzbyLine(incik2,NbNzbyLine(incik2)) =  incf  ! eq ik2 col incf 
          
       enddo
c
c       
c 
ccccccccccccccccccccc       
c
c     Fin comptage 
c
c       write(*,*)
c       do i = 1,NbCV
c          write(*,*)' Line i ',i,NbNzbyLine(i)
c       enddo

c       write(*,*)
c       do i = 1,NbCV
c          write(*,*)' Line i ',i,(NumNzbyLine(i,m),m=1,NbNzbyLine(i))
c       enddo
c
c       stop

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       return
       end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
