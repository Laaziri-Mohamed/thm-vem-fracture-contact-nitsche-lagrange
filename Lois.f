
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
      function f_Lame_mu()

      
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      
      E0 = 10.d+9
      pois0 = 0.25d0
      
         
      rmu = E0/2.d0/(1.d0+pois0)
      

      F_Lame_mu = rmu 
               
      
ccccccccccccccccccccccccccccccccccccccccc
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
c
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
      function f_Lame_lambda()

      
      implicit double precision (a-h,o-z)
      implicit integer (i-n)


      
      E0 = 10.d+9
      pois0 = 0.25d0
               

      rlambda = E0*pois0/(1.d0+pois0)/(1.d0-2.d0*pois0)
      
      

      F_Lame_lambda = rlambda 
               
      
ccccccccccccccccccccccccccccccccccccccccc
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c     coefficient de compression draine Pa
c
c
c      en dimension d:    K0 = 2*mu/d + lambda)
c      ici on est en dimension d = 2 
c      
      function f_K0()

      implicit double precision (a-h,o-z)
      implicit integer (i-n)


      rmu = f_Lame_mu()
         
      rlambda = f_Lame_lambda()
      

      f_K0 = rlambda + rmu ! en dimension d = 2 
      

      return
      end

c
c      
c         
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



      function f_biot()

      
      implicit double precision (a-h,o-z)
      implicit integer (i-n)


c        f_biot = 0.d0       
       
       f_biot = 0.8d0 
c      f_biot = 1.d0 
               
      
ccccccccccccccccccccccccccccccccccccccccc
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



      function f_UnSurN()

      
      implicit double precision (a-h,o-z)
      implicit integer (i-n)




      f_UnSurN = 1.d-10
               

      
ccccccccccccccccccccccccccccccccccccccccc
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc




ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     relaxation matrice en rho*Crm*dP/dt 
c
      function f_Crm()

      
      implicit double precision (a-h,o-z)
      implicit integer (i-n)


      rmu     = f_Lame_mu()
      rlambda = f_Lame_lambda()

      bp = f_biot()
c
c     coeff de relaxation du fixed stress 
c      
        f_Crm = 3.d0*bp**2/(2.d0*rmu + 3.d0*rlambda) ! CAS COUPLE 
               
c      f_Crm = 20.d0*bp**2/(2.d0*rmu + 3.d0*rlambda)  
               
c      f_Crm = 100.d0*bp**2/(2.d0*rmu + 3.d0*rlambda) 
               
c        f_Crm = 0.d0  ! CAS DECOUPLE 
      
ccccccccccccccccccccccccccccccccccccccccc
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc





ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     relaxation fracture en rho*Crf*df*DP/Dt  
c
      function f_Crf()

      
      implicit double precision (a-h,o-z)
      implicit integer (i-n)


      rmu     = f_Lame_mu()
      rlambda = f_Lame_lambda()

      bp = f_biot()
c
c     coeff de relaxation fracture 
c      

c       f_Crf = 3.d0*bp**2/(2.d0*rmu + 3.d0*rlambda)
c       f_Crf = f_Crf*2.d+5
c          f_Crf = f_Crf*3.d+5 ! CAS COUPLE 
c         f_Crf = f_Crf*5.d+5
c         f_Crf = f_Crf*1.d+6


        f_Crf = f_UnSurkn()
      
c         f_Crf = 0.d0 ! CAS DECOUPLE 
               
               
               
      
ccccccccccccccccccccccccccccccccccccccccc
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc




      


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c     unsurkn = 1/kn=  1/Stiffness : lambdan = kn [u]n   si [u]n >= 0 , kn en Pa 
c
c
c     1/kn = 0 -> condition de non penetration 
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
      function f_UnSurkn()

      
      implicit double precision (a-h,o-z)
      implicit integer (i-n)


      

      f_UnSurkn = 1.d0/50.d+9
      
c      f_UnSurkn = 0.d0   
               
      
ccccccccccccccccccccccccccccccccccccccccc
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c     Cohesion: critere de glissement modifie:  |lambdat| < = Cohesion + F lambdan 
c
c     en Pa 
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
      function f_cohesion()

      
      implicit double precision (a-h,o-z)
      implicit integer (i-n)


      

c      f_cohesion = 0.1d+6 
      f_cohesion = 0.d0 
               
      
ccccccccccccccccccccccccccccccccccccccccc
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c     Gravite 
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
      function f_gravite()

      
      implicit double precision (a-h,o-z)
      implicit integer (i-n)


      




      f_gravite = 10.d0 
               
      
ccccccccccccccccccccccccccccccccccccccccc
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c     hauteur du top a l'atmosphere 
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
      function f_htop()

      
      implicit double precision (a-h,o-z)
      implicit integer (i-n)


      


      f_htop = 500.d0 
               
      
ccccccccccccccccccccccccccccccccccccccccc
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c     densite moyenne fluide + roche au dessus du top 
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
      function f_rhosf()

      
      implicit double precision (a-h,o-z)
      implicit integer (i-n)


      


      f_rhosf = 2500.d0 
               
      
ccccccccccccccccccccccccccccccccccccccccc
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
           
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c     ratio des contraintes h sur v pour CL 
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
      function f_ratiohv()

      
      implicit double precision (a-h,o-z)
      implicit integer (i-n)


c      f_ratiohv = 0.d0      

        f_ratiohv = 0.7d0
c       f_ratiohv = 0.65d0
               
      
ccccccccccccccccccccccccccccccccccccccccc
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c     coefficient de frottement 
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
      function f_Frottement()

      
      implicit double precision (a-h,o-z)
      implicit integer (i-n)


      pi = 4.d0*datan(1.d0)
      rf = dtan(25.d0/180.d0*pi)      


      f_Frottement =  rf
               
      
ccccccccccccccccccccccccccccccccccccccccc
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c     coefficient de dilatation par cisaillement  
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
      function f_sheardilation()

      
      implicit double precision (a-h,o-z)
      implicit integer (i-n)


      pi = 4.d0*datan(1.d0)
      sd = dtan(20.d0/180.d0*pi)      


      f_sheardilation =  sd
               
      
ccccccccccccccccccccccccccccccccccccccccc
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c    Beta0n Nitsche 
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
      function f_beta0n()

      
      implicit double precision (a-h,o-z)
      implicit integer (i-n)


cc       f_beta0n = 1000.d0*f_Lame_mu()

      
      f_beta0n = 100.d0*f_Lame_mu()

       
      
ccccccccccccccccccccccccccccccccccccccccc
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c    Beta0t Nitsche 
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
      function f_beta0t()

      
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      
cc        f_beta0t = 100.d0*f_Lame_mu()

      
      f_beta0t = 10.d0*f_Lame_mu()       
       
      
ccccccccccccccccccccccccccccccccccccccccc
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      

            
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c     parametre du semi-smooth Newton pour le contact VEM BULLE  
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
      function f_beta_semi_smooth()

      
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      
      f_beta_semi_smooth = 1.d+6 
 
               
      
ccccccccccccccccccccccccccccccccccccccccc
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
            
      
            
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c     Temps pour le calcul du debit dans JacSm 
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
      function f_temps_debit()

      
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      
c      f_temps_debit = 1.d+9*5.d0/4.d0 
 
      f_temps_debit = 1.d+9 
 
               
      
ccccccccccccccccccccccccccccccccccccccccc
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c           
c      
c     LOIS THERMOS DU FLUIDE 
c
c      
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c     T reference  
c
      function f_TREF()

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      

        f_TREF = 300.d0  
      

      return
      end      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c     P reference 
c
      function f_PREF()

      implicit double precision (a-h,o-z)
      implicit integer (i-n)


      rho0 = f_RHOREF()
      p0 = rho0*f_gravite()*( 1000.d0 + f_htop() )   
      
      f_PREF = p0 
      

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
c      
c     rho reference 
c
      function f_RHOREF()

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      
       f_RHOREF = 1000.d0   
      

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c
c
c
c
c      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c     Module de compression isotherme du fluide: Kf en Pa
c     (eau 20 degres)
c      
      function f_Kf()

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      
       f_Kf = 2.18d+9   
      

      return
      end         
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c     coefficient de dilatation thermique volumique du fluide a pression constante: alphaf en K-1   
c     (eau 20 degres) 
c      
      function f_alphaf()

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      
       f_alphaf = 2.07d-4  
      

      return
      end      
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c    Compressibilite du fluide 
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
      function f_UnSurKf()

      
      implicit double precision (a-h,o-z)
      implicit integer (i-n)


      rKf = f_Kf()
       
      f_UnSurKf = 1.d0/rKf  
          
      
ccccccccccccccccccccccccccccccccccccccccc
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      

c     
c
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c     densite massique kg m^-3 
c
      function f_rho(p,T)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)


      pref = f_PREF()
      Tref = f_TREF()
      rhoref = f_RHOREF()
      
      alphaf = f_alphaf()
      rKf = f_Kf()

      
      ss = 1.d0 - (p-pref)/rKf + alphaf*(T-Tref) 
      f_rho = rhoref/ss 


c      f_rho = 1000.d0 

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c     derivee par rapport a p de la densite massique kg m^-3 
c
      function f_dprho(p,T)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)


      pref = f_PREF()
      Tref = f_TREF()
      rhoref = f_RHOREF()
      
      alphaf = f_alphaf()
      rKf = f_Kf()

      ss = 1.d0 - (p-pref)/rKf + alphaf*(T-Tref)
      dpss = -1.d0/rKf 
      f_dprho =   - rhoref*dpss/ss**2     


c      f_dprho = 0.d0 

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c     derivee par rapport a T de la densite massique kg m^-3 
c
      function f_dTrho(p,T)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      

      pref = f_PREF()
      Tref = f_TREF()
      rhoref = f_RHOREF()
      
      alphaf = f_alphaf()
      rKf = f_Kf()

      ss = 1.d0 - (p-pref)/rKf + alphaf*(T-Tref)
      dTss = alphaf  
      f_dTrho =   - rhoref*dTss/ss**2     


c      f_dTrho = 0.d0 

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c     Chaleur specifique du fluide a pression constante: Cf en J Kg^-1 K^-1    
c     (eau 20 degres) 
c      
      function f_Cf()

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      
       f_Cf = 4180.d0   
      

      return
      end   
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c     energie interne specifique g(p,T) J Kg^-1 K^-1  
c
      function f_ef(p,T)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

       


      pref = f_PREF()
      Tref = f_TREF()
      rhoref = f_RHOREF()

      Cf = f_Cf()
      alphaf = f_alphaf()
      unsKf = f_UnSurKf()

      gref = pref/rhoref + Cf*Tref
      sref = 0.d0 
      
      dp = p - pref
      dT = T - Tref 
      
c      ss0 =  (gref + sref*Tref - pref/rhoref - Cf*Tref) ! zero 
      ss0 = 0.d0 
      
      ss = ss0 + Cf*T  - alphaf/rhoref*Tref*dp  - alphaf/rhoref*p*dT   
      
      ss = ss + (p**2 - pref**2)/(2.d0*rhoref/unsKf) 
      
      f_ef = ss      

       f_ef = Cf*T  
c      f_ef = 0.d0   
      
      return
      end      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c     derivee / a T de l'Energie interne specifique 
c
      function f_dTef(p,T)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)


      pref = f_PREF()
      Tref = f_TREF()
      rhoref = f_RHOREF()

      Cf = f_Cf()
      alphaf = f_alphaf()
      unsKf = f_UnSurKf()      
      

      gref = pref/rhoref + Cf*Tref
      sref = 0.d0 
      
      dp = p - pref
      dT = T - Tref 
            
      ss = Cf - alphaf/rhoref*p  
           
      f_dTef = ss 


      f_dTef = Cf 
c      f_dTef = 0.d0

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c     derivee / a p de l'Energie interne specifique 
c
      function f_dpef(p,T)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)


      pref = f_PREF()
      Tref = f_TREF()
      rhoref = f_RHOREF()

      Cf = f_Cf()
      alphaf = f_alphaf()
      unsKf = f_UnSurKf()

      gref = pref/rhoref + Cf*Tref
      sref = 0.d0 
      
      dp = p - pref
      dT = T - Tref 
            
      ss = - alphaf/rhoref*Tref  - alphaf/rhoref*dT   
      
      ss = ss + p/(rhoref/unsKf) 
      
      
       f_dpef = ss  


      f_dpef = 0.d0 


      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
c      
c     enthalpie specifique g(p,T) J Kg^-1 K^-1  
c
      function f_hf(p,T)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)


      
         f_hf = f_ef(p,T) + p/f_rho(p,T)

      
c        f_hf = f_ef(p,T)  ! hf = ef 
      

      return
      end
c
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c     derivee / a T de l'enthalpie specifique g(p,T) J Kg^-1 K^-1  
c
      function f_dThf(p,T)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      rho = f_rho(p,T)      
      ss = f_dTef(p,T)  - f_dTrho(p,T)*p/rho**2      
      f_dThf = ss
      
      
c       f_dThf = f_dTef(p,T)  ! hf = ef 
      

      return
      end
c      
c
c
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c     derivee / a p de l'enthalpie specifique g(p,T) J Kg^-1 K^-1  
c
      function f_dphf(p,T)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      rho = f_rho(p,T)
      ss = f_dpef(p,T) - f_dprho(p,T)*p/rho**2
      ss = ss + 1.d0/rho
      f_dphf = ss
      

c      f_dphf = f_dpef(p,T) ! hf = ef 
      

      return
      end
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc                  
c      
c      
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c    Chaleur volumique drainee du squelette  J.m^-3 K^-1 
c
c     C0 = (1-phi0) * rhoS * Cs   
c      
      function f_C0()

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      


        f_C0 = 2.d+6    
      
c         f_C0 = 0.d0     
      

      return
      end

c
c      
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c    Coefficient de dilatation thermique volumique  K^-1 
c
      function f_alpha0()

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      

      f_alpha0 = 1.5d-5   

            

      return
      end

c
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c    Coefficient de dilatation isochore des pores J.m^-3 K^-1 
c
c     (b - phi0) alpha0 
c      
      function f_alphaphi()

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      

          f_alphaphi = f_alpha0()*( f_biot() - 0.1d0 ) 


      

      return
      end

c
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc            
