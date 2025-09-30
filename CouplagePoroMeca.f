c
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c
c     calcul de la porosite par maille et de sa derivee par rapport a P 
c
c          PorobyCell et DerPorobyCell 
c      
c     calcul de l'epaisseur frac par face frac
c
c          dfbyFaceFrac 
c      
c     Inclut la relaxation de type fixed-stress dans la matrice seulement 
c
c     Pour le moment on ne met pas la relaxation dans la fracture 
c      
ccccccccccccccccccccccccccccccccccccccccc
c
c     phi = phinm1 + b div( u_prev - u_nm1 ) + 1/N (p - p_nm1) + Crm (p - p_prev ) 
c
c      
c     on suppose df_init = dfcontact -[u_init]_n   avec dfcontact l'ouverture au contact 
c
c     -> df = dfcontact -[u_prev]_n   qui donne moins d'erreurs d'arrondi que
c
c     df = df_nm1 - [u_prev - u_nm1]_n 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
      subroutine Compute_Porosity_Aperture(
     &     NbCell,NbFaceFrac,NbFace,NbCV,
     &     NbdofMecaContact,
     &     dfcontact,
     &     PorobyCell,DerPorobyCell,PorobyCell_nm1,
     &     dfbyFaceFrac, 
     &     SolU_nm1,SolU_prev,SolU0,
     &     IndContact_prev,Sautt_Proj_nm1,       
     &     SolP,SolP_nm1,SolP_prev,
     &     SolT,SolT_nm1,SolT_prev, 
     &     NbIncGlobalbyCell,NumIncGlobalbyCell,GradCellVEM,     
     &     NbIncGlobalbyFaceFrac,NumIncGlobalbyFaceFrac,SautbyFaceFrac,
     &     NumFaceFracVersFace,
     &     VecNormalbyFace,VecTanFace1,VecTanFace2,
     &     NumIncCell)    
      
     

      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      
      include 'include_parameter'


c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     DECLARATIONS
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ENTREES       
c
c     aperture au contact par face frac 
c      
      dimension dfcontact(NbFaceFrac)
c
c     Porosite par maille au temps nm1 
c
      dimension PorobyCell_nm1(NbCell)

c
c     U au temps initial, nm1 et a l'iteration prev (previous) 
c
      dimension SolU0(NbdofMecaContact,NbDim)
      dimension SolU_nm1(NbdofMecaContact,NbDim)
      dimension SolU_prev(NbdofMecaContact,NbDim)
c
c     Indice contact 
c
      dimension IndContact_prev(NbFaceFrac)
c
c     Saut tangentiel moyen projete a nm1 
c
      dimension Sautt_Proj_nm1(NbFaceFrac,NbDim-1)
c
c     Pression courante, au temps nm1 et a l'iteration prev 
c
      dimension SolP(NbCV)
      dimension SolP_nm1(NbCV)
      dimension SolP_prev(NbCV)


c
c     Temperature courante, au temps nm1 et a l'iteration prev 
c
      dimension SolT(NbCV)
      dimension SolT_nm1(NbCV)
      dimension SolT_prev(NbCV)      

c
c        
c     Numerotation locale par cell ds l'ordre NodebyCell inc Ks
c               puis les inc face frac Ksigma si  k=NumCellbyFace(n,1) (cote maille 1 de la face) 
c     Local -> global 
c        
        dimension NbIncGlobalbyCell(NbCell)
        dimension NumIncGlobalbyCell(NbCell,NbIncbyCellMax)             
c
c   operateur gradient scalaire constant par maille 
c
      dimension GradCellVEM(NbCell,NbdofbyCellMax,NbDim)  
c
c     Operateur saut par face frac 
c
        dimension NbIncGlobalbyFaceFrac(NbFaceFrac)
        dimension NumIncGlobalbyFaceFrac(NbFaceFrac,2*NbNodebyFaceMax+2)           
        dimension SautbyFaceFrac(NbFaceFrac,2*NbNodebyFaceMax+2)
c
c        
c     numero de la face pour une face frac donnee  
c        
        dimension NumFaceFracVersFace(NbFaceFrac)
c
c     vecteur normal par face sortant maille 1 de NumCellbyFace 
c     
        dimension VecNormalbyFace(NbFace,NbDim)

c
c     Vecteurs unitaires tangents par face format un repere orthonorme avec le vec normal
c
        dimension VecTanFace1(NbFace,NbDim)
        dimension VecTanFace2(NbFace,NbDim)                       
c
c
c     num cell       >  num inc cell
c
        dimension NumIncCell(NbCell)
c
c        
ccccccccccccccccccccccccccccccccc
c
c     Sorties 
c
c
c     Porosite par maille et derivee par rapport a P et T de la maille 
c
      dimension PorobyCell(NbCell)
      dimension DerPorobyCell(NbCell,2)
c     
c     Ouverture hydraulique de la fracture par face frac 
c
      dimension dfbyFaceFrac(NbFaceFrac)
c
cccccccccccccccccccccccccccccccc
c
c
      dimension Vec(NbDim),deltaU(NbDim)  
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      UnSurN = f_UnSurN() ! module de Biot 
      
      bp = f_biot() ! coefficient de Biot 

c      Crm = f_Crm()             ! coefficient de relaxation de la porosite

      SDilation = f_sheardilation()


       alpha_phi =  f_alphaphi()  ! the volumetric thermal dilation coefficient related to the porosity


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c
c     calcul de l'Aperture 
c
ccccccccccccccccccccccccccccccccc
c
      dfmoy = 0.d0 
      
      do ifrac = 1,NbFaceFrac
         nf = NumFaceFracVersFace(ifrac)
c
c     saut normal de U_prev  
c         
         sn = 0.d0 
         do ik=1,NbIncGlobalbyFaceFrac(ifrac)
            inc = NumIncGlobalbyFaceFrac(ifrac,ik)

c            deltaU(:) = SolU_prev(inc,:)-SolU0(inc,:)
            deltaU(:) = SolU_prev(inc,:)

            pscal = prodscal(deltaU(:),VecNormalbyFace(nf,:)) 
            
            sn = sn + SautbyFaceFrac(ifrac,ik)*pscal                
            
         enddo
c
c     projection (pour Nitsche)
c
c         if (IndContact_prev(ifrac).ne.0) then
c            sn = 0.d0 
c         endif
c         
c
c     calcul de la norme du saut tangentiel de du=U_nm1-U0 (explicite) ou du=U_prev -U0 (implicite)
c         
         st1 = 0.d0
         st2 = 0.d0 
         do ik=1,NbIncGlobalbyFaceFrac(ifrac)
            inc = NumIncGlobalbyFaceFrac(ifrac,ik)         

            deltaU(:) = SolU_prev(inc,:)-SolU_nm1(inc,:)

            pscal1 = prodscal(deltaU,VecTanFace1(nf,:))
            pscal2 = prodscal(deltaU,VecTanFace2(nf,:))
            
            
            st1 = st1 + SautbyFaceFrac(ifrac,ik)*pscal1  
            st2 = st2 + SautbyFaceFrac(ifrac,ik)*pscal2  
     
         enddo

c$$$         st1 = 0.d0
c$$$         st2 = 0.d0 
c$$$         do ik=1,NbIncGlobalbyFaceFrac(ifrac)
c$$$            inc = NumIncGlobalbyFaceFrac(ifrac,ik)
c$$$
c$$$            deltaU(:) = SolU_prev(inc,:)-SolU0(inc,:)            
c$$$c            deltaU(:) = SolU_nm1(inc,:) -SolU0(inc,:)             
c$$$
c$$$
c$$$            pscal1 = prodscal(deltaU,VecTanFace1(nf,:))
c$$$            pscal2 = prodscal(deltaU,VecTanFace2(nf,:))
c$$$                        
c$$$            st1 = st1 + SautbyFaceFrac(ifrac,ik)*pscal1  
c$$$            st2 = st2 + SautbyFaceFrac(ifrac,ik)*pscal2  
c$$$     
c$$$         enddo
c$$$         Vec(1) = st1
c$$$         Vec(2) = st2 
         
         
c
c     projection (pour Nitsche)
c
         if (IndContact_prev(ifrac).eq.1) then
            st1 = 0.d0
            st2 = 0.d0 
         endif
         
         Vec(1) = Sautt_Proj_nm1(ifrac,1) + st1
         Vec(2) = Sautt_Proj_nm1(ifrac,2) + st2


         
         dUtau = dsqrt( Vec(1)**2 + Vec(2)**2 )
         

c         if (dUtau.gt.1.d-10) then 
c            write(*,*) 'dUtau Inc ',ifrac,dUtau
c         endif         
         

         dfbyFaceFrac(ifrac) = dfcontact(ifrac) - sn + dUtau*SDilation ! CAS COUPLE 


c         dfbyFaceFrac(ifrac) = dfcontact(ifrac) ! CAS DECOUPLE 

         

         dfmoy = dfmoy + dfbyFaceFrac(ifrac) 
         
c          write(*,*)' sn df ',ifrac,sn,dfbyFaceFrac(ifrac) 
c          dfbyFaceFrac(ifrac) = dfcontact(ifrac)
         
      enddo

      dfmoy = dfmoy/dfloat(NbFaceFrac)
      

      
c      stop



      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c
c     calcul de la porosite et de sa derivee par rapport a P 
c
cccccccccccccccccccccccccccccc
c
      phimoy = 0.d0


      
      do k=1,NbCell

         inck = NumIncCell(k)
 
         divk = 0.d0 
         do ik=1,NbIncGlobalbyCell(k)
            inc = NumIncGlobalbyCell(k,ik)

            do i = 1,NbDim
               divk = divk + ( SolU_prev(inc,i) - SolU_nm1(inc,i) )
     &                       *GradCellVEM(k,ik,i) ! divergence de (U_prev - U_nm1) cell k 
            enddo       
         enddo
         

         PorobyCell(k) = PorobyCell_nm1(k) + bp*divk 
     &        + UnSurN*( SolP(inck) - SolP_nm1(inck) ) ! CAS COUPLE
     &        -  alpha_phi * ( SolT(inck) - SolT_nm1(inck) )      
         


c         PorobyCell(k) = PorobyCell_nm1(k)
c     &        + UnSurN*( SolP(inck) - SolP_nm1(inck) )  ! CAS DECOUPLE

         
     

c          DerPorobyCell(k) = UnSurN + Crm         
           DerPorobyCell(k,1) = UnSurN ! Der / p de porosite
c           DerPorobyCell(k,2) = 0.d0 ! Der / T de porosite
           DerPorobyCell(k,2) = - alpha_phi ! Der / T de porosite
         
c          DerPorobyCell(k) = 0.d0 


c         write(*,*)' divk ',k,PorobyCell(k),divk,bp 


         phimoy = phimoy + PorobyCell(k)
         
      enddo

      phimoy = phimoy/dfloat(NbCell) 

      write(*,*)' df phim moy ',dfmoy,phimoy 


      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       return
       end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      
