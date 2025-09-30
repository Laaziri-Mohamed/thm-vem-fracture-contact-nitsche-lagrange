ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c

      subroutine solveurSuperLU(NbCV,MemBloc,NbInc,
     &       Aij,nligne_Aij,ncol_Aij,ind_ligne_Aij,
     &       Smm,Soll)



	implicit double precision (a-h,o-z)
	implicit integer (i-n)

cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     JACOBIENNE PAR BLOCS 
c
c
c     structure bloc de la jacobienne: on stocke tous les zeros  
c     lies aux decentrage de facon a avoir la structure bloc maille complete 
c
        dimension Aij(MemBloc,NbInc,NbInc)
c
c
        dimension ncol_Aij(MemBloc)
        dimension nligne_Aij(MemBloc)
        dimension ind_ligne_Aij(NbCV+1)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Sm  
c
        dimension Smm(NbCV,NbInc)
c
c
c     Sol 
c
        dimension Soll(NbCV,NbInc)
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     VARIABLES LOCALES 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Systeme par blocs d'inconnues 
c
c
        double precision, dimension(:), allocatable :: Asyst
        integer, dimension(:), allocatable :: ind_col_Asyst
        integer, dimension(:), allocatable :: nligne_Asyst
        double precision, dimension(:), allocatable :: ncol

       double precision, dimension(:), allocatable :: Sol

c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c


cccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     dimension du vecteur d'inconnues 
c     et nb non zero de la matrice 
c

        
        nsolve = NbCV*NbInc
        mem_syst = MemBloc*NbInc**2 

c        write(*,*)' nsolve ',nsolve
c        write(*,*)' mem_syst ',mem_syst


cccccccccccccccc

        allocate(Sol(nsolve))
c
c
        allocate(Asyst(mem_syst))
        allocate(ind_col_Asyst(nsolve+1))
        allocate(nligne_Asyst(mem_syst))

        allocate(ncol(nsolve))
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     stockage CSC pour SUPERLU 
c
c
c
c     comptage du nb d'elt par colonne ncol
c
        do jj=1,nsolve
           ncol(jj) = 0
        enddo


c
c     diagonale first
c      
        do i=1,NbInc
        do k=1,NbCV
        
c         ii = k + (i-1)*NbCV   
          ii = i + (k-1)*NbInc    
      
           j=i
           m=ind_ligne_Aij(k)+1

             kv = ncol_Aij(m)
c             jj = kv + (j-1)*NbCV
             jj = j + (kv-1)*NbInc             
             ncol(jj) = ncol(jj)+1
        enddo
        enddo



        do i=1,NbInc
        do k=1,NbCV
        
c         ii = k + (i-1)*NbCV   
          ii = i + (k-1)*NbInc    
          
           j=i
           do m=ind_ligne_Aij(k)+2,ind_ligne_Aij(k+1)

             kv = ncol_Aij(m)

c             jj = kv + (j-1)*NbCV
             jj = j + (kv-1)*NbInc                  
             ncol(jj) = ncol(jj)+1

           enddo

           do j=1,NbInc
            if (j.ne.i) then 
               do m=ind_ligne_Aij(k)+1,ind_ligne_Aij(k+1)

                  kv = ncol_Aij(m)                  
c                  jj = kv + (j-1)*NbCV
                  jj = j + (kv-1)*NbInc                     
                  ncol(jj) = ncol(jj)+1
                  
               enddo
            endif
           enddo
  
         enddo
         enddo


cccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     remplissage de ind_col_Asyst(1,...,nsolve+1): 
c
c     colonne jj Asyst de ind_col_Asyst(jj),...,ind_col_Asyst(jj+1)-1
c
c
         ind_col_Asyst(1) = 1
         do jj=1,nsolve
            ind_col_Asyst(jj+1) = ind_col_Asyst(jj) + ncol(jj)
         enddo
c
c
c     remplissage de Asyst
c
        do jj=1,nsolve
           ncol(jj) = 0
        enddo


        do i=1,NbInc
        do k=1,NbCV
        
c         ii = k + (i-1)*NbCV   
          ii = i + (k-1)*NbInc    
      
           j=i
           m=ind_ligne_Aij(k)+1

             kv = ncol_Aij(m)
c             jj = kv + (j-1)*NbCV
             jj = j + (kv-1)*NbInc
             
             ncol(jj) = ncol(jj)+1
             l = ind_col_Asyst(jj)+ncol(jj)-1
             Asyst(l) =  Aij(m,i,j)
             nligne_Asyst(l) = ii
        enddo
        enddo

        do i=1,NbInc
        do k=1,NbCV
        
c         ii = k + (i-1)*NbCV   
          ii = i + (k-1)*NbInc    
      
           j=i
           do m=ind_ligne_Aij(k)+2,ind_ligne_Aij(k+1)

             kv = ncol_Aij(m)

c             jj = kv + (j-1)*NbCV
             jj = j + (kv-1)*NbInc
             
             ncol(jj) = ncol(jj)+1
             l = ind_col_Asyst(jj)+ncol(jj)-1
             Asyst(l) =  Aij(m,i,j)
             nligne_Asyst(l) = ii
           enddo

           do j=1,NbInc
            if (j.ne.i) then 
               do m=ind_ligne_Aij(k)+1,ind_ligne_Aij(k+1)

                  kv = ncol_Aij(m)                  
c                  jj = kv + (j-1)*NbCV
                  jj = j + (kv-1)*NbInc
                  
                  ncol(jj) = ncol(jj)+1
                  l = ind_col_Asyst(jj)+ncol(jj)-1
                  Asyst(l) =  Aij(m,i,j)
                  nligne_Asyst(l) = ii                 
               enddo
            endif
           enddo
  
         enddo
         enddo

c
c
c
c
c
         nnz = ind_col_Asyst(nsolve+1)-1
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc

        do k=1,NbCV
           do i=1,NbInc
c             ii = k + (i-1)*NbCV   
              ii = i + (k-1)*NbInc                  
              Sol(ii) = Smm(k,i)
           enddo
        enddo


ccccccccccccccccccccccccccccccccccccccc
c
c        write(*,*)' debut superlu'
      info = 0
      iopt = 1
      nrhs = 1
      ldb = nsolve
      call c_fortran_dgssv( iopt, nsolve, nnz, nrhs, 
     &                      Asyst, nligne_Asyst, ind_col_Asyst, 
     $                      Sol, ldb, cpufac, cpusolve,info)

c       write(*,*)' info ',info
      
c      write(*,*)' cpu fact solve ',cpufac,cpusolve


      
        do k=1,NbCV
           do i=1,NbInc
c             ii = k + (i-1)*NbCV   
              ii = i + (k-1)*NbInc                
             Soll(k,i) = Sol(ii)
           enddo
        enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccc


        deallocate(Asyst)
        deallocate(ind_col_Asyst)
        deallocate(ncol)
        deallocate(nligne_Asyst)
        deallocate(Sol)

ccccccccccccccccccccccccccccccccccccccc

        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

