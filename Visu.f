cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Visu Paraview : init maillage matrice et maillage frac c
c
c     On suppose des mailles Hexa ou Tetra pour la matrice 
c
c     et des faces quad ou triangle pour les fractures 
c
c
c     En sorties 
c
c     fichiers maillages ensight1.geo000000 matrice 
c
c                     et ensight2.geo000000 fracture 
c
c     Ntetra4,Nhexa8
c
c     Nquad4,Ntria3
c     
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
      subroutine InitVisuParaview(
     &   iaffich,
     &   NbCell,NbFaceFrac,NbNode,NbFace,XNode,
     &   NbNodebyCell,NumNodebyCell,
     &   NumFaceFracVersFace,NbNodebyFace,NumNodebyFace,
     &   NbFacebyCell,NumFacebyCell,      
     &   Ntetra4,Nhexa8,Nnfaced,Nquad4,Ntria3)

c
c
	implicit double precision (a-h,o-z)
	implicit integer (i-n)
c
        include 'include_parameter'
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ENTREES
c
c
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
c
c     Faces par les noeuds 
c       
        dimension NbNodebyFace(NbFace)
        dimension NumNodebyFace(NbFace,NbNodebyFaceMax)
c
c
c
c     Maille par les faces 
c
        dimension NbFacebyCell(NbCell)
        dimension NumFacebyCell(NbCell,NbFacebyCellMax)
c
c        
c       num face frac > num face 
        dimension NumFaceFracVersFace(NbFaceFrac)
c
c
        character*30 xstring        
c
cccccccccccccccccccccccccccccccccccccccccccccccccc

       if (iaffich.eq.1) then 
c
c     Affichage PARAVIEW 
c
c     Creation du fichier maillage au Format Enseight Gold 
c

c    ----- creation du fichier maillage au bon format ------
      open(unit=13,file='ensight1.geo000000',status='unknown')
c   entete du fichier
      write(13,'(A)')'Code frac '
      write(13,'(A)')'Maillage format ensight'
      write(13,'(A)')'node id assign'
      write(13,'(A)')'element id assign'
      write(13,'(A)')'part'
      write(13,'(A)')'         1'
c  coord des noeuds
      write(13,'(A)')'AllCells'
      write(13,'(A)')'coordinates'
      write(13,'(I10)')NbNode
      do i=1,NbNode
         write(13,'(E12.5)')XNode(i,1)
      enddo
      do i=1,NbNode
         write(13,'(E12.5)')XNode(i,2)
      enddo
      do i=1,NbNode
         write(13,'(E12.5)')XNode(i,3)
      enddo
c  cpt type de maille
      Ntetra4 = 0
      Nhexa8 = 0
      Nnfaced = 0      
      do k=1,NbCell
         nsk = NbNodebyCell(k)
c         if (nsk.eq.4) then 
c            Ntetra4 = Ntetra4 + 1
c         elseif (nsk.eq.8) then
c            Nhexa8 = Nhexa8 + 1
c         else if ((nsk.ne.4).and.(nsk.ne.8)) then
            Nnfaced = Nnfaced + 1            
c         endif
      enddo 
c  description des mailles par type
c      if (Ntetra4.gt.0) then
c         write(13,'(A)')'tetra4'
c         write(13,'(I10)')Ntetra4
c         do k=1,NbCell
c            nsk = NbNodebyCell(k) 
c            if (nsk.eq.4) then
c	       write(13,'(4I10)')(NumNodebyCell(k,i),i=1,nsk)
c            endif
c         enddo
c      endif
c      if (Nhexa8.gt.0) then
c         write(13,'(A)')'hexa8'
c         write(13,'(I10)')Nhexa8
c         do k=1,NbCell
c            nsk = NbNodebyCell(k) 
c            if (nsk.eq.8) then
c	       write(13,'(8I10)')(NumNodebyCell(k,i),i=1,nsk)
c            endif
c         enddo
c      endif


      if (Nnfaced.gt.0) then
         write(13,'(A)')'nfaced'
         write(13,'(I10)')Nnfaced
         do k=1,NbCell
            nsk = NbNodebyCell(k) 
c            if ((nsk.ne.4).and.(nsk.ne.8)) then
               write(13,'(I10)')NbFacebyCell(k)              
c            endif
         enddo
         do k=1,NbCell
            nsk = NbNodebyCell(k) 
c            if ((nsk.ne.4).and.(nsk.ne.8)) then
               do ifk = 1,NbFacebyCell(k)    
                  nfk = NumFacebyCell(k,ifk)   
                  write(13,'(I10)')NbNodebyFace(nfk)              
               enddo
c            endif
         enddo
         do k=1,NbCell
            nsk = NbNodebyCell(k) 
c            if ((nsk.ne.4).and.(nsk.ne.8)) then
               do ifk = 1,NbFacebyCell(k)    
                  nfk = NumFacebyCell(k,ifk)  
                  numnode = NbNodebyFace(nfk)
                  if (numnode.lt.10) then 
	       xstring(1:1) = "("
	       write(xstring(2:2), '(I1)' )numnode ! pr format on suppose numnode < 10
	       xstring(3:6)="I10)"
	       write(13,xstring(1:6))(NumNodebyFace(nfk,i),i=1,numnode)
                 else if (numnode.lt.100) then 
	       xstring(1:1) = "("
	       write(xstring(2:3), '(I2)' )numnode ! pr format on suppose numnode < 100
	       xstring(4:7)="I10)"
	       write(13,xstring(1:7))(NumNodebyFace(nfk,i),i=1,numnode)
            else
	       xstring(1:1) = "("
	       write(xstring(2:4), '(I3)' )numnode ! pr format on suppose numnode < 1000
	       xstring(5:8)="I10)"
	       write(13,xstring(1:8))(NumNodebyFace(nfk,i),i=1,numnode)
               
                 endif
               enddo
c            endif
         enddo
      endif

      
      close(13)
      
c    ----------- creation des directory de sortie -----------------
      call system('mkdir output1')
      call system('cd output1 ; rm -rf *')
      call system('cd output1 ; mkdir ensight.geo')
      call system('cd output1 ; mkdir Ux')
      call system('cd output1 ; mkdir Uy')
      call system('cd output1 ; mkdir Uz')
      call system('cd output1 ; mkdir Pm')
      call system('cd output1 ; mkdir Tm')
      call system('cd output1 ; mkdir Cm')      
      call system('cd output1 ; mkdir Ph')      
      call system('mv ensight1.geo000000 output1/ensight.geo')
      
c
c     maillage frac 
c

c    ----- creation du fichier maillage au bon format ------
      open(unit=23,file='ensight2.geo000000',status='unknown')
c   entete du fichier
      write(23,'(A)')'Code frac '
      write(23,'(A)')'Maillage format ensight'
      write(23,'(A)')'node id assign'
      write(23,'(A)')'element id assign'
      write(23,'(A)')'part'
      write(23,'(A)')'         1'
c  coord des noeuds
      write(23,'(A)')'AllCells'
      write(23,'(A)')'coordinates'
      write(23,'(I10)')NbNode
      do i=1,NbNode
         write(23,'(E12.5)')XNode(i,1)
      enddo
      do i=1,NbNode
         write(23,'(E12.5)')XNode(i,2)
      enddo
      do i=1,NbNode
         write(23,'(E12.5)')XNode(i,3)
      enddo
c  cpt type de maille
      Nquad4 = 0
      Ntria3 = 0
      do k=1,NbFaceFrac
         numf = NumFaceFracVersFace(k)
         nsk = NbNodebyFace(numf)
         if (nsk.eq.4) then 
            Nquad4 = Nquad4 + 1
         elseif (nsk.eq.3) then
            Ntria3 = Ntria3 + 1
         else
            write(*,*)'pb pr visu ds nb noeud face frac ',k
            write(*,*)'nsk ne 3 ou 4 : ',nsk
            stop
         endif
      enddo      
c  description des mailles par type
      if (Nquad4.gt.0) then
         write(23,'(A)')'quad4'
         write(23,'(I10)')Nquad4
         do k=1,NbFaceFrac
            numf = NumFaceFracVersFace(k)
            nsk = NbNodebyFace(numf)
            if (nsk.eq.4) then
	       write(23,'(4I10)')(NumNodebyFace(numf,i),i=1,nsk)
            endif
         enddo
      endif
      if (Ntria3.gt.0) then
         write(23,'(A)')'tria3'
         write(23,'(I10)')Ntria3
         do k=1,NbFaceFrac
            numf = NumFaceFracVersFace(k)
            nsk = NbNodebyFace(numf)
            if (nsk.eq.3) then
	       write(23,'(3I10)')(NumNodebyFace(numf,i),i=1,nsk)
            endif
         enddo
      endif
      close(23)
      
c    ----------- creation des directory de sortie -----------------
      call system('mkdir output2')
      call system('cd output2 ; rm -rf *')
      call system('cd output2 ; mkdir ensight.geo')
      call system('cd output2 ; mkdir Sn')
      call system('cd output2 ; mkdir S2')
      call system('cd output2 ; mkdir S3')         
      call system('cd output2 ; mkdir Ln')
      call system('cd output2 ; mkdir L2')
      call system('cd output2 ; mkdir L3')
      call system('cd output2 ; mkdir CS')
      call system('cd output2 ; mkdir Pf')
      call system('cd output2 ; mkdir Tf')
      call system('cd output2 ; mkdir Cf')      
      call system('cd output2 ; mkdir df')      
      call system('mv ensight2.geo000000 output2/ensight.geo')
      
c
c     FIn fichier maillage format Enseight Gold pour PARAVIEW 
c
      endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       return
       end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    Affichage de la solution SolU pour les mailles matrice et du saut aux faces fractures 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
      subroutine VisuSolUP(
     &   iaffich,Ndt,
     &   Ntetra4,Nhexa8,Nnfaced,Nquad4,Ntria3, 
     &   NbCell,NbFaceFrac,NbFace,
     &   NbNodebyCell,
     &   NbNodebyFace,NumFaceFracVersFace,
     &   SolCell,SolSaut,SolLambda,PCell,PFaceFrac,
     &   TCell,TFaceFrac, 
     &   CCell,CFaceFrac, 
     &   PorobyCell,dfbyFaceFrac) 


c
c
	implicit double precision (a-h,o-z)
	implicit integer (i-n)
c
        include 'include_parameter'
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ENTREES
c
c
c     Solution 
c
        dimension SolCell(NbCell,NbDim), SolSaut(NbFaceFrac,NbDim+1)
        dimension SolLambda(NbFaceFrac,NbDim)
        dimension PCell(NbCell)
        dimension PFaceFrac(NbFaceFrac)
        dimension TCell(NbCell)
        dimension TFaceFrac(NbFaceFrac)
        dimension CCell(NbCell)
        dimension CFaceFrac(NbFaceFrac)        
        dimension PorobyCell(NbCell)
        dimension dfbyFaceFrac(NbFaceFrac)
c
c
c     nb de noeuds S par maille K 
c
        dimension NbNodebyCell(NbCell)
c
c     Faces par les noeuds 
c       
        dimension NbNodebyFace(NbFace)
c
c
c       num face frac > num face 
        dimension NumFaceFracVersFace(NbFaceFrac)
c
c
c
c     VISU PARAVIEW 
c
        character*30 xstring
        character*60 commande
        character*6 endF
c        
c
cccccccccccccccccccccccccccccccccccccccccccccccccc

       if (iaffich.eq.1) then 
c
c    --------------------- ecriture de la solution   ----------------------
c
          if (Ndt<10) then
             endF(1:5)='00000'
             write( endF(6:6),'(I1)') Ndt
          elseif (Ndt<100) then
             endF(1:4)='0000'
             write( endF(5:6),'(I2)') Ndt
          elseif (Ndt<1000) then
             endF(1:3)='0000'
             write( endF(4:6),'(I3)') Ndt
          else
             write(*,*)'pb nb pas de temps ds post pro'
             stop
          endif
          
          xstring(3:8) = endF(1:6)

          xstring(1:2) = 'Ux'
          open(unit=11,file =xstring(1:8),status = 'unknown')
          
          xstring(1:2) = 'Uy'
          open(unit=21,file =xstring(1:8),status = 'unknown')
          
          xstring(1:2) = 'Uz'
          open(unit=31,file =xstring(1:8),status = 'unknown')
          
          xstring(1:2) = 'Pm'
          open(unit=41,file =xstring(1:8),status = 'unknown')

          xstring(1:2) = 'Ph'
          open(unit=51,file =xstring(1:8),status = 'unknown')                 

          xstring(1:2) = 'Cm'
          open(unit=61,file =xstring(1:8),status = 'unknown')          

          xstring(1:2) = 'Tm'
          open(unit=71,file =xstring(1:8),status = 'unknown')

           

          xstring(1:2) = 'Sn'
          open(unit=12,file =xstring(1:8),status = 'unknown')
          
          xstring(1:2) = 'S2'
          open(unit=22,file =xstring(1:8),status = 'unknown')
          
          xstring(1:2) = 'S3'
          open(unit=32,file =xstring(1:8),status = 'unknown')          

          xstring(1:2) = 'Ln'
          open(unit=13,file =xstring(1:8),status = 'unknown')
          
          xstring(1:2) = 'L2'
          open(unit=23,file =xstring(1:8),status = 'unknown')
          
          xstring(1:2) = 'L3'
          open(unit=33,file =xstring(1:8),status = 'unknown')
          
          xstring(1:2) = 'CS'
          open(unit=34,file =xstring(1:8),status = 'unknown')
          
          xstring(1:2) = 'Pf'
          open(unit=35,file =xstring(1:8),status = 'unknown')
          
          xstring(1:2) = 'df'
          open(unit=36,file =xstring(1:8),status = 'unknown')

          xstring(1:2) = 'Cf'
          open(unit=37,file =xstring(1:8),status = 'unknown')

          xstring(1:2) = 'Tf'
          open(unit=38,file =xstring(1:8),status = 'unknown')

      
      write(11,'(A)')'Ux'
      write(21,'(A)')'Uy'
      write(31,'(A)')'Uz'
      write(41,'(A)')'Pm'
      write(51,'(A)')'Ph'
      write(61,'(A)')'Cm'
      write(71,'(A)')'Tm'
      
      write(12,'(A)')'Sn'
      write(22,'(A)')'S2'
      write(32,'(A)')'S3'
      
      write(13,'(A)')'Ln'
      write(23,'(A)')'L2'
      write(33,'(A)')'L3'
      
      write(34,'(A)')'CS'
      write(35,'(A)')'Pf'
      write(36,'(A)')'df'
      write(37,'(A)')'Cf'
      write(38,'(A)')'Tf'

      

      write(11,'(A)')'part'
      write(11,'(A)')'         1'
      
      write(12,'(A)')'part'
      write(12,'(A)')'         1'

      write(13,'(A)')'part'
      write(13,'(A)')'         1'
      
      write(21,'(A)')'part'
      write(21,'(A)')'         1'

      write(22,'(A)')'part'
      write(22,'(A)')'         1'

      write(23,'(A)')'part'
      write(23,'(A)')'         1'

      write(31,'(A)')'part'
      write(31,'(A)')'         1'

      write(41,'(A)')'part'
      write(41,'(A)')'         1'


      write(51,'(A)')'part'
      write(51,'(A)')'         1'

      write(61,'(A)')'part'
      write(61,'(A)')'         1'

      write(71,'(A)')'part'
      write(71,'(A)')'         1'    

      
      write(32,'(A)')'part'
      write(32,'(A)')'         1'

      write(33,'(A)')'part'
      write(33,'(A)')'         1'          

      write(34,'(A)')'part'
      write(34,'(A)')'         1'        
      
      write(35,'(A)')'part'
      write(35,'(A)')'         1'        
      
      write(36,'(A)')'part'
      write(36,'(A)')'         1'
      
      write(37,'(A)')'part'
      write(37,'(A)')'         1'

      write(38,'(A)')'part'
      write(38,'(A)')'         1'  

      
c      if (Ntetra4.gt.0) then
c         write(11,'(A)')'tetra4'
c         do k=1,NbCell
c            nsk = NbNodebyCell(k) 
c            inck = NumIncCell(k)
c            if (nsk.eq.4) then
c	       write(11,'(E12.5)')Sol(inck)
c            endif
c         enddo
c      endif
c      if (Nhexa8.gt.0) then
c         write(11,'(A)')'hexa8'
c         do k=1,NbCell
c            nsk = NbNodebyCell(k) 
c            inck = NumIncCell(k)
c            if (nsk.eq.8) then
c	       write(11,'(E12.5)')Sol(inck)
c            endif
c         enddo
c     endif
      
      if (Nnfaced.gt.0) then
         write(11,'(A)')'nfaced'
         write(21,'(A)')'nfaced'
         write(31,'(A)')'nfaced'
         write(41,'(A)')'nfaced'
         write(51,'(A)')'nfaced'
         write(61,'(A)')'nfaced'
         write(71,'(A)')'nfaced' 
         do k=1,NbCell
            nsk = NbNodebyCell(k) 
c            if ((nsk.ne.4).and.(nsk.ne.8)) then                  
            write(11,'(E12.5)')SolCell(k,1)
            write(21,'(E12.5)')SolCell(k,2)
            write(31,'(E12.5)')SolCell(k,3)
            write(41,'(E12.5)')PCell(k)
            write(51,'(E12.5)')PorobyCell(k)
            write(61,'(E12.5)')CCell(k)
            write(71,'(E12.5)')TCell(k)
c            endif
         enddo
      endif
      
      if (Nquad4.gt.0) then
         write(12,'(A)')'quad4'
         write(22,'(A)')'quad4'
         write(32,'(A)')'quad4'
         write(13,'(A)')'quad4'
         write(23,'(A)')'quad4'
         write(33,'(A)')'quad4'
         write(34,'(A)')'quad4'
         write(35,'(A)')'quad4'
         write(36,'(A)')'quad4'
         write(37,'(A)')'quad4'
         write(38,'(A)')'quad4'
         
         do k=1,NbFaceFrac
            numf = NumFaceFracVersFace(k)
            nsk = NbNodebyFace(numf)
            if (nsk.eq.4) then
	       write(12,'(E12.5)')SolSaut(k,1)
	       write(22,'(E12.5)')SolSaut(k,2)
	       write(32,'(E12.5)')SolSaut(k,3)
               
	       write(34,'(E12.5)')SolSaut(k,4)
	       write(35,'(E12.5)')PFaceFrac(k)
	       write(36,'(E12.5)')dfbyFaceFrac(k)
	       write(37,'(E12.5)')CFaceFrac(k)
               write(38,'(E12.5)')TFaceFrac(k)

               
	       write(13,'(E12.5)')SolLambda(k,1)
	       write(23,'(E12.5)')SolLambda(k,2)
	       write(33,'(E12.5)')SolLambda(k,3)
               
            endif
         enddo
      endif
      if (Ntria3.gt.0) then
         write(12,'(A)')'tria3'
         write(22,'(A)')'tria3'
         write(32,'(A)')'tria3'
         
         write(13,'(A)')'tria3'
         write(23,'(A)')'tria3'
         write(33,'(A)')'tria3'
         write(34,'(A)')'tria3'
         write(35,'(A)')'tria3'
         write(36,'(A)')'tria3'
         write(37,'(A)')'tria3'
         write(38,'(A)')'tria3'
         
         do k=1,NbFaceFrac
            numf = NumFaceFracVersFace(k)
            nsk = NbNodebyFace(numf)
            if (nsk.eq.3) then
	       write(12,'(E12.5)')SolSaut(k,1)
	       write(22,'(E12.5)')SolSaut(k,2)
	       write(32,'(E12.5)')SolSaut(k,3)

	       write(34,'(E12.5)')SolSaut(k,4)
	       write(35,'(E12.5)')PFaceFrac(k)
	       write(36,'(E12.5)')dfbyFaceFrac(k)
	       write(37,'(E12.5)')CFaceFrac(k)
               write(38,'(E12.5)')TFaceFrac(k)
               
               
	       write(13,'(E12.5)')SolLambda(k,1)
	       write(23,'(E12.5)')SolLambda(k,2)
	       write(33,'(E12.5)')SolLambda(k,3)
               
            endif
         enddo
      endif
      
      close(11)
      close(21)
      close(31)
      close(41)
      close(51)
      close(61)
      close(71)
      
      close(12)
      close(22)
      close(32)
      
      close(13)
      close(23)
      close(33)
      
      close(34)
      close(35)
      close(36)
      close(37)
      close(38)
      


      
c     --------------------- move vers le bon directory -------------------------------

         xstring(1:2) = 'Ux'
         call system('mv '//xstring(1:8)//' output1/Ux') 

         xstring(1:2) = 'Uy'
         call system('mv '//xstring(1:8)//' output1/Uy') 

         xstring(1:2) = 'Uz'
         call system('mv '//xstring(1:8)//' output1/Uz') 

         xstring(1:2) = 'Pm'
         call system('mv '//xstring(1:8)//' output1/Pm')

         xstring(1:2) = 'Tm'
         call system('mv '//xstring(1:8)//' output1/Tm')
         
         xstring(1:2) = 'Ph'
         call system('mv '//xstring(1:8)//' output1/Ph')
         
         xstring(1:2) = 'Cm'
         call system('mv '//xstring(1:8)//' output1/Cm') 
          

         xstring(1:2) = 'Sn'
         call system('mv '//xstring(1:8)//' output2/Sn') 

         xstring(1:2) = 'S2'
         call system('mv '//xstring(1:8)//' output2/S2') 

         xstring(1:2) = 'S3'
         call system('mv '//xstring(1:8)//' output2/S3') 

         xstring(1:2) = 'Ln'
         call system('mv '//xstring(1:8)//' output2/Ln') 
         
         xstring(1:2) = 'L2'
         call system('mv '//xstring(1:8)//' output2/L2') 

         xstring(1:2) = 'L3'
         call system('mv '//xstring(1:8)//' output2/L3') 

         xstring(1:2) = 'CS'
         call system('mv '//xstring(1:8)//' output2/CS')
         
         xstring(1:2) = 'Pf'
         call system('mv '//xstring(1:8)//' output2/Pf')

         xstring(1:2) = 'Tf'
         call system('mv '//xstring(1:8)//' output2/Tf')
         
         xstring(1:2) = 'df'
         call system('mv '//xstring(1:8)//' output2/df') 

         xstring(1:2) = 'Cf'
         call system('mv '//xstring(1:8)//' output2/Cf') 

         
      

      
      endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       return
       end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    Affichage de la solution Sol pour les mailles matrice et les faces fractures 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
      subroutine FichierEnseightCase(iaffich,Ndt)


c
c
	implicit double precision (a-h,o-z)
	implicit integer (i-n)
c
        include 'include_parameter'
c
ccccccccccccccccccccccccccccccccccccccccccc
c

       if (iaffich.eq.1) then 

c    ! creation du fichier ensight.case
      open(unit=13,file ='ensight.case',status = 'unknown')
      write(13,'(A)')'FORMAT'
      write(13,'(A)')'type: ensight gold'
      write(13,*)
      write(13,'(A)')'GEOMETRY'
      write(13,'(A)')'model: 1 ensight.geo/ensight1.geo000000'
      write(13,*)
      write(13,'(A)')'VARIABLE'
      write(13,'(A)')'scalar per element:    1   Ux Ux/Ux******'
      write(13,'(A)')'scalar per element:    1   Uy Uy/Uy******'
      write(13,'(A)')'scalar per element:    1   Uz Uz/Uz******'
      write(13,'(A)')'scalar per element:    1   Pm Pm/Pm******'
      write(13,'(A)')'scalar per element:    1   Tm Tm/Tm******'
      write(13,'(A)')'scalar per element:    1   Ph Ph/Ph******'
      write(13,'(A)')'scalar per element:    1   Cm Cm/Cm******'          
      write(13,*)
      write(13,'(A)')'TIME'
      write(13,'(A,i6)')'time set:              ',1
      write(13,'(A,i6)')'number of steps:       ',Ndt
      write(13,'(A)')'filename start number: 1'
      write(13,'(A)')'filename increment:    1'
      write(13,'(A)')'time values:'
      do i = 1,Ndt
         write(13,*)i
      enddo
      close(13)
      call system('mv ensight.case output1/') 

c
c
c
c     Maillage frac 
c

c    ! creation du fichier ensight.case
      open(unit=23,file ='ensight.case',status = 'unknown')
      write(23,'(A)')'FORMAT'
      write(23,'(A)')'type: ensight gold'
      write(23,*)
      write(23,'(A)')'GEOMETRY'
      write(23,'(A)')'model: 1 ensight.geo/ensight2.geo000000'
      write(23,*)
      write(23,'(A)')'VARIABLE'
      write(23,'(A)')'scalar per element:    1   Sn Sn/Sn******'
      write(23,'(A)')'scalar per element:    1   S2 S2/S2******'
      write(23,'(A)')'scalar per element:    1   S3 S3/S3******'
      write(23,'(A)')'scalar per element:    1   Ln Ln/Ln******'
      write(23,'(A)')'scalar per element:    1   L2 L2/L2******'
      write(23,'(A)')'scalar per element:    1   L3 L3/L3******'
      write(23,'(A)')'scalar per element:    1   CS CS/CS******'
      write(23,'(A)')'scalar per element:    1   Pf Pf/Pf******'
      write(23,'(A)')'scalar per element:    1   Tf Tf/Tf******'
      write(23,'(A)')'scalar per element:    1   df df/df******'
      write(23,'(A)')'scalar per element:    1   Cf Cf/Cf******'     
      write(23,*)
      write(23,'(A)')'TIME'
      write(23,'(A,i6)')'time set:              ',1
      write(23,'(A,i6)')'number of steps:       ',Ndt
      write(23,'(A)')'filename start number: 1'
      write(23,'(A)')'filename increment:    1'
      write(23,'(A)')'time values:'
      do i = 1,Ndt
         write(23,*)i
      enddo
      close(23)
      call system('mv ensight.case output2/') 


      endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       return
       end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
