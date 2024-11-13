#You can run all with ctrl+alt+r



#dans cette version 10, on rajoute le code nécéssaire à la partie cross 
#validation. On rajoute aussi la fonction check (voir ligne 738)

#install.packages("pracma")
library(pracma)#utilis? pour calcul? la norme d'un vecteur
#install.packages("qlcMatrix")
library(qlcMatrix)
#install.packages("xtable")
library(xtable)
#install.packages("GoFKernel")
library(GoFKernel)



################################################################################
################################################################################
###In this part, we visualize the "Enriched Laplace Distribution". The Laguerre 
###function computes the Laguerre polynomials.
###

## ELD.pdf function computes th p.d.f. of the ELD.
## y corresponds to the set of points for which we are interested in evaluating
## the pdf.
## Tau, mu and theta are the parameters of the distribution.
##
##
ELD.pdf<-function(y,tau,mu,theta,theta_tilde){
  
  #we separe the observations larger than mu from the ones smaller than mu
  yu<-y[y>=mu]
  yd<-y[y<mu]
  
  #we compute n+ (and n-) the number of observations larger (smaller) than mu
  nplus<-length(yu)
  nminus<-length(yd)
  
  n=nplus+nminus
  
  #correspond to the m in our equation
  m <- length(theta)-1
  m_tilde <- length(theta_tilde)-1
  
  #In this part, we compute the Laguerre polynomial :
  # The Laguerre polynomial of order k evaluated in x is given by :
  #\sum_k=0^m (k choose m) (-x)^k / k!
  #We compute this in 3 separated parts. In each part, the code is repeated for
  #observations smaller and larger than mu.
  
  #Part 1 computes the binomial coefficients needed for the computation of the 
  #Laguerre Polynomial. These coefficients are stored in a (lower) triangular
  #matrix with ones on the diagonal. In other words, Part 1 computes k choose m
  #and store this in a matrix called B.
  B <- matrix(0,ncol=m+1,nrow=m+1)
  for(k in 1:(m+1)) {
    B[k,1:k] <- choose(k-1,0:(k-1))
  }
  
  B_tilde <- matrix(0,ncol=m_tilde+1,nrow=m_tilde+1)
  for(k in 1:(m_tilde+1)) {
    B_tilde[k,1:k] <- choose(k-1,0:(k-1))
  }
  
  
  #In part 2, we compute the part " (-x)^k / k!" of the Laguerre polynomial.  
  #We store these results in a matrix with m columns and n lines (where n is the
  #number of observations). So the element of the i-th row and j-th column of Z
  #corresponds to (-x_i ^(j-1))/j-1!. The j-th column is given by the j-1th 
  #column multiplied by -x/j.
  #We also note that the Laguerre polynomial is evauated in tau*(y-mu) for y >mu
  #and in -(1-tau)(y-mu) for y < mu.
  Zu <- matrix(1,ncol=m+1,nrow=length(yu))
  if(m+1>=2) {
    for(i in 2:(m+1)) {
      Zu[,i] <- Zu[,i-1]*(-tau*(yu-mu))/(i-1)
    }
  }
  
  Zd <- matrix(1,ncol=m_tilde+1,nrow=length(yd))
  if(m_tilde+1>=2) {
    for(i in 2:(m_tilde+1)) {
      Zd[,i] <- Zd[,i-1]*((1-tau)*(yd-mu))/(i-1)
    }
  }
  
  #Part 3 computes the coefficients of the Laguerre Polynomial. These coeffi-
  #cients are stored in a (lower) triangular matrix. The i+1th-row of the matrix
  #contains the k+1 coefficients of i-th order Laguerre polynomial. In this
  #part, the results of part 1 and 2 are put together and the matrician product 
  #computes the sum in the definition of the Laguerre polynomial.
  Lu <- matrix(0,ncol=m+1,nrow=m+1)
  Lu <- Zu%*%t(B)
  
  Ld <- matrix(0,ncol=m_tilde+1,nrow=m_tilde+1)
  Ld <- Zd%*%t(B_tilde)
  
  
  #Now that the Laguerre polynomial is computed, we can construct the pdf of the
  #Enriched Laplace distribution.
  
  f<-rep(0,n)
  f[y>=mu]<-exp(-tau*(yu-mu))/Norm(theta)^2*(Lu%*%theta)^2
  f[y<mu]<-exp((1-tau)*(yd-mu))/Norm(theta_tilde)^2*(Ld%*%theta_tilde)^2
  f<-f*(1-tau)*tau
  return(f)
}

################################################################################
################################################################################
###Our goal now is to compute the loglikelihood and to find the parameters that
###maximises this loglikelihood.

##The ELD.loglik function computes the likelihood for a vector of observations y
##given that the ELD has parameters tau, mu, theta, theta tilde. Here y 
##represents a vector of size n of observations.
##The beginning of the function is the same as for function ELD.pdf but of
##course the end is different.
ELD.loglik<-function(y,tau,mu,theta,theta_tilde){
  
  #we separe the observations larger than mu from the ones smaller than mu
  yu<-y[y>=mu]
  yd<-y[y<mu]
  
  #we compute n+ (and n-) the number of observations larger (smaller) than mu
  nplus<-length(yu)
  nminus<-length(yd)
  
  n=nplus+nminus
  
  #correspond to the m in our equation
  m <- length(theta)-1
  m_tilde <- length(theta_tilde)-1
  
  #In this part, we compute the Laguerre polynomial :
  # The Laguerre polynomial of order k evaluated in x is given by :
  #\sum_k=0^m (k choose m) (-x)^k / k!
  #We compute this in 3 separated parts. In each part, the code is repeated for
  #observations smaller and larger than mu.
  
  #Part 1 computes the binomial coefficients needed for the computation of the 
  #Laguerre Polynomial. These coefficients are stored in a (lower) triangular
  #matrix with ones on the diagonal. In other words, Part 1 computes k choose m
  #and store this in a matrix called B.
  B <- matrix(0,ncol=m+1,nrow=m+1)
  for(k in 1:(m+1)) {
    B[k,1:k] <- choose(k-1,0:(k-1))
  }
  
  B_tilde <- matrix(0,ncol=m_tilde+1,nrow=m_tilde+1)
  for(k in 1:(m_tilde+1)) {
    B_tilde[k,1:k] <- choose(k-1,0:(k-1))
  }
  
  
  #In part 2, we compute the part " (-x)^k / k!" of the Laguerre polynomial.  
  #We store these results in a matrix with m columns and n lines (where n is the
  #number of observations). So the element of the i-th row and j-th column of Z
  #corresponds to (-x_i ^(j-1))/j-1!. The j-th column is given by the j-1th 
  #column multiplied by -x/j.
  #We also note that the Laguerre polynomial is evauated in tau*(y-mu) for y >mu
  #and in -(1-tau)(y-mu) for y < mu.
  Zu <- matrix(1,ncol=m+1,nrow=length(yu))
  if(m+1>=2) {
    for(i in 2:(m+1)) {
      Zu[,i] <- Zu[,i-1]*(-tau*(yu-mu))/(i-1)
    }
  }
  
  Zd <- matrix(1,ncol=m_tilde+1,nrow=length(yd))
  if(m_tilde+1>=2) {
    for(i in 2:(m_tilde+1)) {
      Zd[,i] <- Zd[,i-1]*((1-tau)*(yd-mu))/(i-1)
    }
  }
  
  #Part 3 computes the coefficients of the Laguerre Polynomial. These coeffi-
  #cients are stored in a (lower) triangular matrix. The i+1th-row of the matrix
  #contains the k+1 coefficients of i-th order Laguerre polynomial. In this
  #part, the results of part 1 and 2 are put together and the matrician product 
  #computes the sum in the definition of the Laguerre polynomial.
  Lu <- matrix(0,ncol=m+1,nrow=m+1)
  Lu <- Zu%*%t(B)
  
  Ld <- matrix(0,ncol=m_tilde+1,nrow=m_tilde+1)
  Ld <- Zd%*%t(B_tilde)
  
  #Now that the Laguerre polynomial is calculated, we can compute the 
  #log-likelihood by using the formula that we derived
  
  
  return(
    n*log(tau*(1-tau)) +
      -tau*sum(yu-mu) -2*nplus*log(Norm(theta)) + sum(log((Lu%*%theta)^2))
    +(1-tau)*sum(yd-mu) -2*nminus*log(Norm(theta_tilde)) + sum(log((Ld%*%theta_tilde)^2))
  )
}

##Now we want to compute the MLE of the parameters of the distribution.

#On utilise optim pour optimiser la logvraisemblance. La fonction optim a en 
#entr?e une fonction et un point de d?part pour l'optimisation. L'entr?e qui est
#est une fonction est optimis?e par rapport ? tous ses param?tres.

#Puisqu'ici, on ne veut pas optimiser par rapport ? tous les param?tres (on
#optimise pas par rapport ? tau par exemple) on doit cr?er des sous fonctions 
#qui d?pendent juste des param?tres qu'on souhaite optimiser. Ces osus fonctions
#appelent exactement la fonction ELD.loglik.

#Par ailleurs la fonction optim ne sait apparemment que minimiser (et pas
#maximiser). Ceci explique la pr?sence du -1* plus bas).


ELD.loglik2<-function(param2){
  m=length(theta)-1;m_tilde=length(theta_tilde)-1
  mu=param2[1]
  return(-1*ELD.loglik(y,tau,mu,theta,theta_tilde))
}

#Apres la version2 on avait remarqu? un probl?me d'identifiabilit?.
#On avait d?cid? de fixer theta0 = theta tilde = 0 pour r?gler le probleme.

#Avant
ELD.loglik3<-function(param3){
  m=length(theta)-1;m_tilde=length(theta_tilde)-1
  mu=param3[1]
  
  if (m != 0 & m_tilde != 0) {
    theta=c(1,param3[2:(2+m-1)])
    theta_tilde=c(1,param3[(2+m):(2+m+m_tilde-1)])
    return(-1*ELD.loglik(y,tau,mu,theta,theta_tilde))
  } else if (m == 0 & m_tilde == 0) {
    theta=1
    theta_tilde=1
    return(-1*ELD.loglik(y,tau,mu,theta,theta_tilde))
  } else if (m == 0 & m_tilde != 0) {
    theta=1
    theta_tilde=c(1,param3[2:(1+m_tilde)])
    return(-1*ELD.loglik(y,tau,mu,theta,theta_tilde))
  } else if (m != 0 & m_tilde == 0) {
    theta_tilde=1
    theta=c(1,param3[2:(1+m)])
    return(-1*ELD.loglik(y,tau,mu,theta,theta_tilde))
  }  
  
}


################################################################################
################################################################################

"Dans cette partie, on construit la cdf de l'ELD. Cela sera n?c?ssaire pour 
simuler des donn?es ? partir de cette distribution. On fait ca en int?grant la 
pdf. Pour int?grer on a besoin d'une fonction d'un parametre d'o? la cr?ation
de la fonction ELD.pdf2"

ELD.pdf2<-function(x){
  ELD.pdf(x,tau,mu,theta,theta_tilde)
}

##x = valeurs pour lesquelles on calcule la cdf 
ELD.cdf<-function(x,tau,mu,theta,theta_tilde){
  y<-rep(NA,length(x))
  tau=tau;
  mu=mu;
  theta=theta;
  theta_tilde=theta_tilde
  for (i in 1:length(y)) {
    y[i]<-integrate(ELD.pdf2, lower = Inf, upper = x[i])$value    
  }
  return(y)
}


#In this part, we construct a function which simulates data from an "Enriched
#Laplace distribution" (ELD). The simulation method used is the rejection method.
#To do this, data are generated uniformly on a rectangle sufficiently large to
#cover a large part of the p.d.f. The simulated points situated beyond the
#curve are kept and the simulated points situated above the curve are excluded
#The n points that were kept give our simulated data from the ELD.
#This part uses the function ELD.pdf available earlier.

#On suppose mu est toujours le mode de la ELD. On commence par calculer la
#valeur de la pdf en ce point mu. Ce maximum est appell? C. 
#Pour la distribution de Laplace classique (non asym?trique), C vaut 0.25.
#On prend le maximum de f(?-) et f(?+) au cas o? il y a discontinuit? en mu.

#Ensuite, on d?termine les quantiles 0.001 et 0.999. J'ai choisi ces deux 
#quantiles de mani?re ? repr?senter la plus grande partie possible de la 
#distribution sans perdre trop de temps. Cela implique qu'on ne pourra pas 
#simuler de donn?es qui se trouvent en dehors de ces deux quantiles. Ces deux 
#quantiles sont calcul?s ? l'aide de qELD.

#Apr?s, on g?n?re des points uniform?ment dans le rectangle form?s par 
#{(x,y)| (x,y) in [q001 ; q999] X [0,C]}. Parmi tous ces points, on garde ceux
#qui sont situ?s en dessous de la pdf et on exclut les autres. Les abscisses des
#points simul?s situ?s en-dessous de la pdf donneront les donn?es simul?es.

#La fonction ELD.simul permet de g?n?rer des observations d'une ELD de parame-
#-tres mu, tau, theta, theta_tilde.
#n repr?sente le nombre d'observations que l'on souhaite simuler.
#Cette fonction fonctionne de la mani?re suivante :
#Tout d'abord on simule un point dans le rectangle pr?c?demment d?crit.
#Si ce point est situ? sous la courbe de la pdf, alors on le garde. Sinon on le
#met de c?t?. La variable count compte le nombre d'observations qui ont ?t? simu
#-l?es. Ce compte augmente de 1 ? chaque fois qu'un point simul? est gard?. La
#fonction s'interromp une fois que ce count atteint n.
#Les donn?es simul?es sont stock?es dans le vecteur appell? simul.
ELD.simul<-function(n,tau,mu,theta,theta_tilde){
  simul<-rep(NA,n)
  count<-0
  C<-max(ELD.pdf(y=mu-0.00001,tau,mu,theta,theta_tilde),ELD.pdf(y=mu+0.00001,tau,mu,theta,theta_tilde))
  q.005<-qELD(0.005,tau,mu,theta,theta_tilde)
  q.995<-qELD(0.995,tau,mu,theta,theta_tilde)
  while (count < n) {
    U1<-runif(1,q.005,q.995)
    U2<-runif(1,0,C)
    if(ELD.pdf(y=U1,tau,mu,theta,theta_tilde)>U2) {
      simul[count+1]<-U1
      count<-count+1
    }
    if (count==n) {return(simul)}
  }
}





################################################################################
################cdf en utilisant la (longue) formule explicite##################

#La fonction ELD.cdf est un peu lente. Puisqu'on va en avoir grandement besoin 
#pour maximiser la vraisemblance dans le cas de donn?es censur?es, j'ai pris le
#temps de calculer la forme explicite de la fonction de r?partition. Cette forme
#explicite est impl?ment?e ci dessous. La forme ?tant assez longue, je l'ai 
#construite petit ? petit.


ELD.cdf2<-function(y,tau,mu,theta,theta_tilde){
  m <- length(theta)-1
  m_tilde <- length(theta_tilde)-1
  
  ################################### D?but partie 1 #############################
  #Cette partie contient la partie o? on calcule tous les coefficients
  #Son temps de calcul n'est pas tres optimis? mais cette partie n'est ?x?cut?e
  #qu'une fois et les matrices sont souvent assez petites). Le code est r?p?t?
  #deux fois 
  Au <- matrix(0,ncol=m+1,nrow=m+1)
  for (k in 1:(m+1)) {
    Au[k,1:k] <- choose(k-1,0:(k-1)) *(-1)^(0:(k-1))/ factorial(0:(k-1))
  }
  
  Ad <- matrix(0,ncol=m_tilde+1,nrow=m_tilde+1)
  for (k in 1:(m_tilde+1)) {
    Ad[k,1:k] <- choose(k-1,0:(k-1)) *(-1)^(0:(k-1))/ factorial(0:(k-1))
  }
  
  
  #B=(j+j')! / l!
  Bu <- matrix(0,ncol=2*m+1,nrow=2*m+1)
  for(j.j_prime in 0:(2*m))
  {
    Bu[j.j_prime+1,(1:(j.j_prime+1))]<-factorial(j.j_prime)/factorial(0:j.j_prime)
  }
  #Bu
  
  
  Bd <- matrix(0,ncol=2*m_tilde+1,nrow=2*m_tilde+1)
  for(j.j_prime in 0:(2*m_tilde))
  {
    Bd[j.j_prime+1,(1:(j.j_prime+1))]<-factorial(j.j_prime)/factorial(0:j.j_prime)
  }
  #Bd
  
  
  #C contient toutes les combinaisons du produit de A j,k avec A jprime,kprime
  Cu <- array(0, dim=c(m+1,m+1,m+1,m+1))
  
  for (k in 0:m) {
    for (k_prime in 0:m) {
      for (j in 0:k) {
        for (j_prime in 0:k_prime) {
          Cu[j+1,j_prime+1,k+1,k_prime+1]<-Au[k+1,j+1] * Au[k_prime+1,j_prime+1]
        }
      }
    }
  }
  
  
  Cd <- array(0, dim=c(m_tilde+1,m_tilde+1,m_tilde+1,m_tilde+1))
  
  for (k in 0:m_tilde) {
    for (k_prime in 0:m_tilde) {
      for (j in 0:k) {
        for (j_prime in 0:k_prime) {
          Cd[j+1,j_prime+1,k+1,k_prime+1]<-Ad[k+1,j+1] * Ad[k_prime+1,j_prime+1]
        }
      }
    }
  }
  
  ###Cetta partie est un peu plus rapide que la pr?c?dente mais il y a un probleme avec la sym?trie. Bizarre :'(
  #  tic()
  
  #    for (j in 0:k) {
  #      for (j_prime in 0:k_prime) {
  #       my.array[j+1,j_prime+1,1:(m+1),1:(m+1)]<-A[(0:m)+1,j+1] * A[1:(m+1),j_prime+1]
  
  #      }
  #   }
  
  #toc()
  #my.array
  
  
  #Maintenant on veut connaitre le coefficient devant le terme en j* du produit
  #de vk et vkprime (avec j* in 0,1,...,k+k'). Pour ca on prend :
  #pour une valeur de k,k' la somme tel des ??ments tels que j+j' = j*
  
  #Dk kprime jpjprime = sum  Ajk*Ajprimekprime tel que j+jprime = jpjprime
  
  #Un as.mtrix a du etre rajout? car visiblement les fonctions row et col  ne 
  #marchent pas quand m et/ou mtilde = 1
  
  Du <- array(0, dim=c(m+1,m+1,2*m+1))
  for (k in 1:(m+1)) {
    for (k_prime in 1:(m+1)) {
      for (j_star in 1:(k+k_prime-1)) {
        Du[k,k_prime,j_star]<-
          sum(Cu[,,k,k_prime][(row(as.matrix(Cu[,,k,k_prime])) +col(as.matrix(Cu[,,k,k_prime])) -1)==j_star])
      }
    }
  }
  
  Dd <- array(0, dim=c(m_tilde+1,m_tilde+1,2*m_tilde+1))
  for (k in 1:(m_tilde+1)) {
    for (k_prime in 1:(m_tilde+1)) {
      for (j_star in 1:(k+k_prime-1)) {
        Dd[k,k_prime,j_star]<-
          sum(Cd[,,k,k_prime][(row(as.matrix(Cd[,,k,k_prime])) +col(as.matrix(Cd[,,k,k_prime])) -1)==j_star])
      }
    }
  }
  
  
  ################################### D?but partie 2 #############################
  
  index<-rank(y)
  y<-sort(y)
  
  yu <- y[y>=mu]
  su <- (yu-mu)*tau
  lenu<-length(su)
  resu <- rep(NA,lenu)#ce vecteur contiendra le r?sultat pour les valeurs 
  #demand?es sup?rieures ? ?.
  if(lenu == 0) {resu=numeric(0)} else
    if(lenu != 0) {
      #E contient -(s)^l * exp(-s) pour diff?rentes valeurs de l dans 0,...,2m
      Eu<-matrix(NA,2*m +1,lenu)
      for (l in 1:(2*m+1)) {
        Eu[l,] = (su)^(l-1) * exp(-su) #pas sur du moins ici, ? v?rifier.  
      }
      
      "Gj+jprime,l  contient - j+j!/l! s^l exp(-s) pour toutes les valeurs de l 
  #entre 0 et 2m et de j + j_prime entre 0 et 2m"
      Gu <- array(0,dim=c(2*m+1, 2*m+1,lenu))
      for (j.j_prime in 1:(2*m+1)) {
        for (l in 1:(j.j_prime)) { #-1 cark+kprime sont plus longs de 1
          Gu[j.j_prime,l,]<- Bu[j.j_prime,l] * Eu[l,]
        }
      }
      #Hj+jprime contient - sum l=0^j+jprime j+j!/l! s^l exp(-s)
      Hu<-array(0,dim=c(2*m+1,lenu))
      if(m!=0){
        for (i in 1:lenu) {
          Hu[,i] <- rowSums(Gu[,,i])
        }
      } else if (m==0) {  
        for (i in 1:lenu) {
          Hu[,i] <- Gu[,,i]
        }}
      #J contient 
      Ju<- array(0,dim=c(m+1, m+1, lenu))
      for (k in 1:(m+1)) {
        for (k_prime in 1:(m+1)) {
          Ju[k,k_prime,]<- (Du[k,k_prime,] %*% Hu)
        }
      }
      THETA <- (theta) %*% t(theta)
      
      #Le probleme avec J C'est que c'est difficile de faire des sommes 
      #sans utiliser de boucles for sur y. Ainsi, on transforme l'array en une
      #matrice avec m+1^2 colonnes et lenu lignes et on fait le sop?rations 
      #n?c?ssaires sur cette matrice K
      Ku <- matrix(Ju[,,1:lenu],ncol = (m+1)^2,nrow=lenu,byrow=T)*matrix(rep(THETA),byrow=T,ncol = (m+1)^2,nrow=lenu)
      
      if(m!=0){
        resu <- 1 + (1 -tau)* -rowSums(Ku)/Norm(theta)^2 
      } else if (m==0) {
        resu <- 1 + ((1-tau)* -Ku/Norm(theta)^2) 
      }
      
    }
  
  
  ##############################################################################
  #Pareil pour les observations < mu
  yd <- y[y<mu]
  sd <- -(yd-mu)*(1-tau)
  lend<-length(sd)
  if(lend == 0) {resd=numeric(0)} else
    if(lend != 0){
      
      resd<-rep(NA,lend)
      
      Ed<-matrix(NA,2*m_tilde +1,lend)
      for (l in 1:(2*m_tilde+1)) {
        Ed[l,] = (sd)^(l-1) * exp(-sd) #pas sur du moins ici, ? v?rifier.  
      }
      
      Gd <- array(0,dim=c(2*m_tilde+1, 2*m_tilde+1,lend))
      for (j.j_prime in 1:(2*m_tilde+1)) {
        for (l in 1:(j.j_prime)) { #-1 cark+kprime sont plus longs de 1
          Gd[j.j_prime,l,]<- Bd[j.j_prime,l] * Ed[l,]
        }
      }
      
      Hd<-array(0,dim=c(2*m_tilde+1,lend))
      if(m_tilde!=0){
        for (i in 1:lend) {
          Hd[,i] <- rowSums(Gd[,,i])
        }
      } else if (m_tilde==0) {  
        for (i in 1:lend) {
          Hd[,i] <- Gd[,,i]
        }}
      
      Jd<- array(0,dim=c(m_tilde+1, m_tilde+1, lend))
      for (k in 1:(m_tilde+1)) {
        for (k_prime in 1:(m_tilde+1)) {
          Jd[k,k_prime,]<- (Dd[k,k_prime,] %*% Hd)
        }
      }
      THETA_tilde <- (theta_tilde) %*% t(theta_tilde)
      
      Kd <- matrix(Jd[,,1:lend],ncol = (m_tilde+1)^2,nrow=lend,byrow=T)*matrix(rep(THETA_tilde),byrow=T,ncol = (m_tilde+1)^2,nrow=lend)
      
      if(m_tilde!=0){
        resd <- -tau* -rowSums(Kd)/Norm(theta_tilde)^2 
      } else if (m_tilde==0) {
        resd <- -tau* -Kd/Norm(theta_tilde)^2 
      }
      
    }
  res<-c(resd,resu)
  indexinfinite = is.nan(res) #Note : ces deux lignes ci sont présentes dans le 
  #cas où on introduit de la censure à droite. Si on a de la censure à droite,
  #le bord droit de l'intervalle est +infini. A l'étape du calcul de Eu cela 
  #génère des NaN car pour R exp(-infini) * infini^a (avec a >0) ne vaut pas 0 
  #mais vaut NaN. Note 2 cette modification ne convient pas si le bord gauche
  res[indexinfinite] = 1 #d'un intervalle est -infini
  res[index]
  
}


ELD.cdf3<-function(x){
  ELD.cdf2(x,tau,mu,theta,theta_tilde)
}

#Cr?ation de la fonction quantile. Ici il faudrait investiguer plus sur ce que
#fait exactement la fonction inverse 

qELD<-function(q,tau,mu,theta,theta_tilde){
  ELD.inversecdf<-inverse(ELD.cdf3)
  return(ELD.inversecdf(q))
}





rELD2<-function(n,tau=0.5,mu=0,theta=c(1),theta_tilde=c(1)){
  setwd("C:/Users/bdeketelaere/Documents/Sauvegardessimulations")
  m=length(theta)-1
  m_tilde=length(theta_tilde)-1
  if( (m==1) & (m_tilde==1) & (mu == 0) & (tau == 0.1) & (theta[1] ==1 ) & (theta[2] ==0.5 )  & (theta_tilde[1]==1) & (theta_tilde[2]==0.5) ){
    samp<-as.matrix(read.table("mu=0;tau=0,1;theta=c(1,0,5);theta_tilde=c(1,0,5)"))
    index<-runif(n,1,9999)
    return(samp[floor(index)]/2 + samp[ceiling(index)]/2)
  } else if( (m==1) & (m_tilde==1) & (mu == 0) & (tau == 0.9) & (theta[1] ==1 ) & (theta[2] ==0.5 )  & (theta_tilde[1]==1) & (theta_tilde[2]==0.5) ){
    samp<-as.matrix(read.table("mu=0;tau=0,9;theta=c(1,0,5);theta_tilde=c(1,0,5)"))
    index<-runif(n,1,9999)
    return(samp[floor(index)]/2 + samp[ceiling(index)]/2)
  } else if( (m==2) & (m_tilde==1) & (mu == 0) & (tau == 0.1) & (theta[1] ==1 ) & (theta[2] ==0.5 ) & (theta[3] ==0.25) & (theta_tilde[1]==1) & (theta_tilde[2]==1) ){
    samp<-as.matrix(read.table("mu=0;tau=0,1;theta=c(1,0,5,0,25);theta_tilde=c(1,1)"))
    index<-runif(n,1,9999)
    return(samp[floor(index)]/2 + samp[ceiling(index)]/2)
  } else if( (m==2) & (m_tilde==0) & (mu == 0) & (tau == 0.5) & (theta[1] ==1 ) & (theta[2] ==0.5 ) & (theta[3] ==0.25) & (theta_tilde[1]==1) ){
    samp<-as.matrix(read.table("mu=0;tau=0,5;theta=c(1,0,5,0,25);theta_tilde=1"))
    index<-runif(n,1,49999)
    return(samp[floor(index)]/2 + samp[ceiling(index)]/2)
  } else if( (m==2) & (m_tilde==0) & (mu == 0) & (tau == 0.1) & (theta[1] ==1 ) & (theta[2] ==0.5 ) & (theta[3] ==0.25) & (theta_tilde[1]==1) ){
    samp<-as.matrix(read.table("mu=0;tau=0,1;theta=c(1,0,5,0,25);theta_tilde=1"))
    index<-runif(n,1,49999)
    return(samp[floor(index)]/2 + samp[ceiling(index)]/2)
  } else if( (m==1) & (m_tilde==0) & (mu == 0) & (tau == 0.45) & (theta[1] ==1 ) & (theta[2] ==0.5 ) & (theta_tilde[1]==1) ){
    samp<-as.matrix(read.table("mu=0;tau=0,45;theta=c(1,0,5);theta_tilde=1"))
    index<-runif(n,50,24999)
    return(samp[floor(index)]/2 + samp[ceiling(index)]/2)
  }  else if( (m==0) & (m_tilde==0) & (mu == 0) & (tau == 0.5) & (theta[1] ==1 ) & (theta_tilde[1]==1) ){
    samp<-as.matrix(read.table("mu=0;tau=0,5;theta=1;theta_tilde=1"))
    index<-runif(n,1,24999)
    return(samp[floor(index)]/2 + samp[ceiling(index)]/2)
  }  else { 
    inverse.all<-rep(NA,1000)
    inverse.fun.all<-inverse(ELD.cdf3)
    for (i in 1:1000) {
      inverse.all[i]<-inverse.fun.all(i/(1000+1))
    }
    index<-runif(n,1,1000)
    floorindex<-floor(index)
    return(inverse.all[floorindex]/2 + inverse.all[floorindex+1]/2)  
  }
}


################################################################################
################################################################################
###In this part, the goal is to simulate Interval Censored data (IC).
###At the beginning, we will do this for original data coming from a normal
###distribution. Then, we will simulate IC data from original data coming from 
###an ELD. Finally, we will repeat this for data coming from other distributions

###The method used here comes "from Tutorial on methods for interval-censored 
###data and their implementation in R" from : Guadalupe Gomez, M Luz Calle and 
###Klaus Langohr published in Statistical Modelling 2009; 9(4): 259-297.

##The method consists in :
##Step 1 : Simulate T from a known distribution
##Step 2 : Simulate U1 and U2 both from a uniform distribution (min=0,max=c)
##Step 3 : Choose a value for c. 2c is the largest possible interval observable.
##The smaller it is, the easier it will be. 
##Step 4 : Construct L and R by the following formulas
## L=max(T+U2-c,T-U1) and R=min(T+U2,T-U1+c)

simul.IC.norm <- function(n,mu=0,sigma=1,c){
  X <- rnorm(n,mu,sigma)
  U1 <- runif(n,0,c) 
  U2 <- runif(n,0,c)
  
  L1<-X+U2-c;L2<-X-U1
  L12<-cbind(L1,L2)
  L<-as.vector(rowMax(L12))
  
  R1<-X+U2;R2<-X-U1+c
  R12<-cbind(R1,R2)
  R<-as.vector(rowMin(R12))
  return(cbind(L,R))
}

##Ensuite, on fait la m?me pour des donn?es originales issues de la distribution
##de Laplace. On utilise pour cela la fonction ELD.simul qui a ?t? d?crite plus
##haut. Cette fonction a pour but de simuler des donn?es d'une ELDistribution.

simul.IC.ELD <- function(n,c,tau=0.5,mu=0,theta=c(1),theta_tilde=c(1),...){
  X <- rELD2(n,tau,mu,theta,theta_tilde)
  u1 <- runif(n,0,c) 
  u2 <- runif(n,0,c)
  
  L1<-X+u2-c;L2<-X-u1
  L12<-cbind(L1,L2)
  L<-as.vector(rowMax(L12))
  
  R1<-X+u2;R2<-X-u1+c
  R12<-cbind(R1,R2)
  R<-as.vector(rowMin(R12))
  cbind(L,R)
}



################################################################################
#Now the same with a gamma distribution.
simul.IC.gamma <- function(n,shape,rate,c){
  X <- rgamma(n,shape,rate)
  U1 <- runif(n,0,c) 
  U2 <- runif(n,0,c)
  
  L1<-X+U2-c;L2<-X-U1
  L12<-cbind(L1,L2)
  L<-as.vector(rowMax(L12))
  
  R1<-X+U2;R2<-X-U1+c
  R12<-cbind(R1,R2)
  R<-as.vector(rowMin(R12))
  return(cbind(L,R))
}

################################################################################
#Now the same with a Weibull distribution.
simul.IC.weib <- function(n,shape=1.5,scale=2,c){
  X <- rweibull(n,shape=shape,scale=scale)
  U1 <- runif(n,0,c) 
  U2 <- runif(n,0,c)
  
  L1<-X+U2-c;L2<-X-U1
  L12<-cbind(L1,L2)
  L<-as.vector(rowMax(L12))
  
  R1<-X+U2;R2<-X-U1+c
  R12<-cbind(R1,R2)
  R<-as.vector(rowMin(R12))
  return(cbind(L,R))
}


################################################################################
#As for ELD.loglik versions 1 2 and 3, we compute the loglikelihood for Interval
#censored observations. Note that the contribution of observation i to the 
#likelihood is F(R_i)-F(L_i)

loglik.IC<-function(mu,theta,theta_tilde,y,tau)
{
  if(mu > max(y[,2])){return(-Inf)}
  n <- nrow(y)
  longy <- c(y[,1],y[,2])
  res<-ELD.cdf2(longy,tau,mu,theta,theta_tilde)
  sum(log(res[(n+1):(2*n)]-res[1:n]))
}


loglik.IC2<-function(param3){
  m=length(theta)-1;m_tilde=length(theta_tilde)-1
  mu=param3[1]
  if (m != 0 & m_tilde != 0) {
    theta=c(1,param3[2:(2+m-1)])
    theta_tilde=c(1,param3[(2+m):(2+m+m_tilde-1)])
    return(-1*loglik.IC(mu,theta,theta_tilde,y,tau))
  } else if (m == 0 & m_tilde == 0) {
    theta=1
    theta_tilde=1
    return(-1*loglik.IC(mu,theta,theta_tilde,y,tau))
  } else if (m == 0 & m_tilde != 0) {
    theta=1
    theta_tilde=c(1,param3[2:(1+m_tilde)])
    return(-1*loglik.IC(mu,theta,theta_tilde,y,tau))
  } else if (m != 0 & m_tilde == 0) {
    theta_tilde=1
    theta=c(1,param3[2:(1+m)])
    return(-1*loglik.IC(mu,theta,theta_tilde,y,tau))
  }  
}

loglik.IC3<-function(param2){
  m=length(theta)-1;m_tilde=length(theta_tilde)-1
  mu=param2[1]
  #theta=c(param2[2:(2+m)])
  #theta_tilde=param2[(2+m+1):(2+m+1+m_tilde)]
  return(-1*loglik.IC(mu,theta,theta_tilde,y,tau))
}

################################################################################
################################################################################
################Cross validation (K-fold and Leave-one-out)#####################
################################################################################
################################################################################


## Implementation of the check function
check_function <- function(z,tau) {
  return(z*(tau-as.numeric(z<=0)))
}


################################################################################
################################################################################
###################### Penalized Likelihood (via LASSO) ########################
################################################################################
################################################################################

#Dans cette partie  On tente l'approche de vraisemblance pénalisée avec LASSO.


PEN.loglik<-function(lambda,lambda_tilde,y,tau=0.5,mu=0,theta,theta_tilde){
  m=length(theta)-1;m_tilde=length(theta_tilde)-1
  n=length(y)
  
  loglik=ELD.loglik(y,tau,mu,theta,theta_tilde)
  return(loglik/n - lambda*Norm(theta[2:(m+1)],p=1)-lambda_tilde*Norm(theta_tilde[(2:(m_tilde+1))],p=1))
  
}



PEN.loglik3<-function(param3){
  m=length(theta)-1;m_tilde=length(theta_tilde)-1
  mu=param3[1]
  if (m != 0 & m_tilde != 0) {
    theta=c(1,param3[2:(2+m-1)])
    theta_tilde=c(1,param3[(2+m):(2+m+m_tilde-1)])
    return(-1*PEN.loglik(lambda,lambda_tilde,y,tau,mu,theta,theta_tilde))
  } else if (m == 0 & m_tilde == 0) {
    theta=1
    theta_tilde=1
    return(-1*PEN.loglik(lambda,lambda_tilde,y,tau,mu,theta,theta_tilde))
  } else if (m == 0 & m_tilde != 0) {
    theta=1
    theta_tilde=c(1,param3[2:(1+m_tilde)])
    return(-1*PEN.loglik(lambda,lambda_tilde,y,tau,mu,theta,theta_tilde))
  } else if (m != 0 & m_tilde == 0) {
    theta_tilde=1
    theta=c(1,param3[2:(1+m)])
    return(-1*PEN.loglik(lambda,lambda_tilde,y,tau,mu,theta,theta_tilde))
  }  
}


#Pareil qu'avant mais on ajoute la censure maintenant


PEN.loglik.IC<-function(lambda,lambda_tilde,y,tau=0.5,mu=0,theta,theta_tilde){
  m=length(theta)-1;m_tilde=length(theta_tilde)-1
  n=length(y)/2 ######################   Attention à ce divisé par 2
  
  loglik=loglik.IC(mu,theta,theta_tilde,y,tau)
  return(loglik/n - lambda*Norm(theta[2:(m+1)],p=1)-lambda_tilde*Norm(theta_tilde[(2:(m_tilde+1))],p=1))
  
}

#PEN.loglik.IC(0,0,y,tau,mu,theta,theta_tilde)
#loglik.IC(mu,theta,theta_tilde,y,tau)



PEN.loglik.IC3<-function(param3){
  m=length(theta)-1;m_tilde=length(theta_tilde)-1
  mu=param3[1]
  if (m != 0 & m_tilde != 0) {
    theta=c(1,param3[2:(2+m-1)])
    theta_tilde=c(1,param3[(2+m):(2+m+m_tilde-1)])
    return(-1*PEN.loglik.IC(lambda,lambda_tilde,y,tau,mu,theta,theta_tilde))
  } else if (m == 0 & m_tilde == 0) {
    theta=1
    theta_tilde=1
    return(-1*PEN.loglik.IC(lambda,lambda_tilde,y,tau,mu,theta,theta_tilde))
  } else if (m == 0 & m_tilde != 0) {
    theta=1
    theta_tilde=c(1,param3[2:(1+m_tilde)])
    return(-1*PEN.loglik.IC(lambda,lambda_tilde,y,tau,mu,theta,theta_tilde))
  } else if (m != 0 & m_tilde == 0) {
    theta_tilde=1
    theta=c(1,param3[2:(1+m)])
    return(-1*PEN.loglik.IC(lambda,lambda_tilde,y,tau,mu,theta,theta_tilde))
  }  
}

#lambda=0;lambda_tilde=0
#PEN.loglik.IC3(c(0,0.5,0.5))
#PEN.loglik.IC(0,0,y,tau,mu,theta,theta_tilde)






##Now we want to compute the MLE of the parameters of the distribution.

#On utilise optim pour optimiser la logvraisemblance. La fonction optim a en 
#entr?e une fonction et un point de d?part pour l'optimisation. L'entr?e qui est
#est une fonction est optimis?e par rapport ? tous ses param?tres.

#Puisqu'ici, on ne veut pas optimiser par rapport ? tous les param?tres (on
#optimise pas par rapport ? tau par exemple) on doit cr?er des sous fonctions 
#qui d?pendent juste des param?tres qu'on souhaite optimiser. Ces osus fonctions
#appelent exactement la fonction ELD.loglik.

#Par ailleurs la fonction optim ne sait apparemment que minimiser (et pas
#maximiser). Ceci explique la pr?sence du -1* plus bas).


ELD.loglik2<-function(param2){
  m=length(theta)-1;m_tilde=length(theta_tilde)-1
  mu=param2[1]
  return(-1*ELD.loglik(y,tau,mu,theta,theta_tilde))
}

#Apres la version2 on avait remarqu? un probl?me d'identifiabilit?.
#On avait d?cid? de fixer theta0 = theta tilde = 0 pour r?gler le probleme.

#Avant
ELD.loglik3<-function(param3){
  m=length(theta)-1;m_tilde=length(theta_tilde)-1
  mu=param3[1]
  
  if (m != 0 & m_tilde != 0) {
    theta=c(1,param3[2:(2+m-1)])
    theta_tilde=c(1,param3[(2+m):(2+m+m_tilde-1)])
    return(-1*ELD.loglik(y,tau,mu,theta,theta_tilde))
  } else if (m == 0 & m_tilde == 0) {
    theta=1
    theta_tilde=1
    return(-1*ELD.loglik(y,tau,mu,theta,theta_tilde))
  } else if (m == 0 & m_tilde != 0) {
    theta=1
    theta_tilde=c(1,param3[2:(1+m_tilde)])
    return(-1*ELD.loglik(y,tau,mu,theta,theta_tilde))
  } else if (m != 0 & m_tilde == 0) {
    theta_tilde=1
    theta=c(1,param3[2:(1+m)])
    return(-1*ELD.loglik(y,tau,mu,theta,theta_tilde))
  }  
  
}



################################################################################
################################################################################
############################### Juste pour moi #################################
################################################################################
################################################################################

setting <- function(lambda=F){
  
  if(lambda==T){cat("\n mu =", mu , "tau = ", tau, "sigma =", sigma , "\n theta =", theta, "thetatilde =", theta_tilde,
                    "\n m =", m , "mtilde =", m_tilde, "\n lambda =", lambda, "lambdatilde =", lambda_tilde, "\n alpha = ", alpha)
  } else {cat("\n mu =", mu , "tau = ", tau, "sigma =", sigma , "\n theta =", theta, "thetatilde =", theta_tilde,
              "\n m =", m , "mtilde =", m_tilde)
  }
}




#A remetrre dans un ordre plus clair apres :

myinverse <- function (f, lower = -Inf, upper = Inf) 
{
  if (!is.numeric(lower) || !is.numeric(upper) || lower >= 
      upper) 
    stop("lower < upper is not fulfilled")
  if (!is.finite(f(lower)) & is.finite(lower)) 
    lower <- lower + .Machine$double.eps
  if (!is.finite(f(upper)) & is.finite(upper)) 
    upper <- upper - .Machine$double.eps
  if (is.infinite(lower) & is.infinite(upper)) {
    function(y) optim(0, (function(x) abs(f(x) - y)), method = "Brent",lower = -200,upper=200)$par
  }
  else if (is.infinite(lower) & is.finite(upper)) {
    function(y) optim(upper, (function(x) abs(f(x) - y)), 
                      lower = lower, upper = upper, method = "L-BFGS-B")$par
  }
  else if (is.finite(lower) & is.infinite(upper)) {
    function(y) optim(lower, (function(x) abs(f(x) - y)), 
                      lower = lower, upper = upper, method = "L-BFGS-B")$par
  }
  else {
    if (f(lower) < f(upper)) {
      function(y) uniroot((function(x) f(x) - y), lower = lower, 
                          upper = upper)$root
    }
    else {
      function(y) optim((lower + upper)/2, (function(x) (f(x) - 
                                                           y)^2), lower = lower, upper = upper, method = "L-BFGS-B")$par
    }
  }
}


#q est une seule valeur et pas un vecteur
qELD2<-function(q,tau,mu,theta,theta_tilde){
  ELD.inversecdf<-myinverse(ELD.cdf3)
  return(ELD.inversecdf(q))
}


#Expliquer aussi le 1.05
#rejection method with q.001 and q.999 computed by my inverse

rELD4 <- function(n,tau=0.5,mu=0,theta=c(1),theta_tilde=c(1)){
  setwd("C:/Users/bdeketelaere/Documents/Sauvegardessimulations/Data Simulation 3")
  m=length(theta)-1
  m_tilde=length(theta_tilde)-1
  if( (m==0) & (m_tilde==0) ){
    return(ald::rALD(n,p=tau))
  } else if( (m==2) & (m_tilde==0) & (mu == 0) & (tau == 0.1) & (theta[1] ==1 ) & (theta[2] ==0.5 ) & (theta[3] ==0.25) & (theta_tilde[1]==1) ){
    samp<-as.matrix(read.table("mu=0;tau=0,1;theta=c(1,0,5,0,25);theta_tilde=1",nrows=1000000))
    index<-runif(n,1,1000000)
    return(samp[round(index)] )
  } else if( (m==1) & (m_tilde==1) & (mu == 0) & (tau == 0.1) & (theta[1] ==1 ) & (theta[2] ==0.5 ) & (theta_tilde[1]==1) & (theta_tilde[2]==0.5) ){
    samp<-as.matrix(read.table("mu=0;tau=0,1;theta=c(1,0,5);theta_tilde=c(1,0,5)",nrows=1000000))
    index<-runif(n,1,1000000)
    return(samp[round(index)] )
  } else if( (m==1) & (m_tilde==1) & (mu == 0) & (tau == 0.5) & (theta[1] ==1 ) & (theta[2] ==0.5 ) & (theta_tilde[1]==1) & (theta_tilde[2]==0.5) ){
    samp<-as.matrix(read.table("mu=0;tau=0,5;theta=c(1,0,5);theta_tilde=c(1,0,5)",nrows=1000000))
    index<-runif(n,1,1000000)
    return(samp[round(index)] )
  } else if( (m==1) & (m_tilde==0) & (mu == 0) & (tau == 0.5) & (theta[1] ==1 ) & (theta[2] ==0.5 ) & (theta_tilde[1]==1) ){
    samp<-as.matrix(read.table("mu=0;tau=0,5;theta=c(1,0,5);theta_tilde=c(1)",nrows=1000000))
    index<-runif(n,1,1000000)
    return(samp[round(index)] )
  } else if( (m==1) & (m_tilde==1) & (mu == 0) & (tau == 0.9) & (theta[1] ==1 ) & (theta[2] ==0.5 ) & (theta_tilde[1]==1) & (theta_tilde[2]==0.5) ){
    samp<-as.matrix(read.table("mu=0;tau=0,9;theta=c(1,0,5);theta_tilde=c(1,0,5)",nrows=1000000))
    index<-runif(n,1,1000000)
    return(samp[round(index)] )
  } else if( (m==1) & (m_tilde==1) & (mu == 0) & (tau == 0.5) & (theta[1] ==1 ) & (theta[2] ==1 ) & (theta_tilde[1]==1) & (theta_tilde[2]==1) ){
    samp<-as.matrix(read.table("mu=0;tau=0,5;theta=c(1,1);theta_tilde=c(1,1)",nrows=1000000))
    index<-runif(n,1,1000000)
    return(samp[round(index)] )
  } else if( (m==2) & (m_tilde==2) & (mu == 0) & (tau == 0.5) & (theta[1] ==1 ) & (theta[2] ==1 ) & (theta[3] ==1 ) & (theta_tilde[1]==1) & (theta_tilde[2]==1) & (theta_tilde[3]==1) ){
    samp<-as.matrix(read.table("mu=0;tau=0,5;theta=c(1,1,1);theta_tilde=c(1,1,1)",nrows=1000000))
    index<-runif(n,1,1000000)
    return(samp[round(index)] )
  } else if( (m==2) & (m_tilde==2) & (mu == 0) & (tau == 0.5) & (theta[1] ==1 ) & (theta[2] ==1 ) & (theta[3] ==1 ) & (theta_tilde[1]==1) & (theta_tilde[2]==1) & (theta_tilde[3]==1) ){
    samp<-as.matrix(read.table("mu=0;tau=0,5;theta=c(1,1,1);theta_tilde=c(1,1,1)",nrows=1000000))
    index<-runif(n,1,1000000)
    return(samp[round(index)] )
  } else { 
    simul<-rep(NA,n)
    count<-0
    q.001=qELD2(0.0001,tau,mu,theta,theta_tilde)
    q.999=qELD2(0.9999,tau,mu,theta,theta_tilde)
    C=1.01*max(ELD.pdf(mu-0.001,tau,mu,theta,theta_tilde),ELD.pdf(mu+0.001,tau,mu,theta,theta_tilde))
    while (count < n) {
      U1<-runif(1,q.001,q.999)
      U2<-runif(1,0,C)
      if(ELD.pdf(y=U1,tau,mu,theta,theta_tilde)>U2) {
        simul[count+1]<-U1
        count<-count+1
      }
      if(count%%500000==0){print(count)}
      if (count==n) {return(simul)}
    }
  }
}





#rejection method with q.001 and q.999 computed by my inverse

#rELD5 <- function(n,tau=0.5,mu=0,theta=c(1),theta_tilde=c(1)){
#  setwd("C:/Users/bdeketelaere/Documents/Sauvegardessimulations/Data Simulation2")
#  m=length(theta)-1
#  m_tilde=length(theta_tilde)-1
#  if( (m==0) & (m_tilde==0) ){
#    return(rALD(n,p=tau))
#  } else if( (m==2) & (m_tilde==0) & (mu == 0) & (tau == 0.1) & (theta[1] ==1 ) & (theta[2] ==0.5 ) & (theta[3] ==0.25) & (theta_tilde[1]==1) ){
#    samp<-as.matrix(read.table("mu=0;tau=0,1;theta=c(1,0,5,0,25);theta_tilde=1",nrows=1000000))
#    index<-runif(n,1,1000000)
#    return(samp[round(index)] )
#  } else if( (m==1) & (m_tilde==1) & (mu == 0) & (tau == 0.1) & (theta[1] ==1 ) & (theta[2] ==0.5 ) & (theta_tilde[1]==1) & (theta_tilde[2]==0.5) ){
#    samp<-as.matrix(read.table("mu=0;tau=0,1;theta=c(1,0,5);theta_tilde=c(1,0,5)",nrows=1000000))
#    index<-runif(n,1,1000000)
#    return(samp[round(index)] )
#  } else if( (m==1) & (m_tilde==1) & (mu == 0) & (tau == 0.5) & (theta[1] ==1 ) & (theta[2] ==0.5 ) & (theta_tilde[1]==1) & (theta_tilde[2]==0.5) ){
#    samp<-as.matrix(read.table("mu=0;tau=0,5;theta=c(1,0,5);theta_tilde=c(1,0,5)",nrows=1000000))
#    index<-runif(n,1,1000000)
#    return(samp[round(index)] )
#  } else if( (m==1) & (m_tilde==0) & (mu == 0) & (tau == 0.5) & (theta[1] ==1 ) & (theta[2] ==0.5 ) & (theta_tilde[1]==1) ){
#    samp<-as.matrix(read.table("mu=0;tau=0,5;theta=c(1,0,5);theta_tilde=c(1)",nrows=1000000))
#    index<-runif(n,1,1000000)
#    return(samp[round(index)] )
#  } else { 
#    simul<-rep(NA,n)
#    count<-0
#    q.001=qELD2(0.001,tau,mu,theta,theta_tilde)
#    q.999=qELD2(0.999,tau,mu,theta,theta_tilde)
#    C=1.05*max(ELD.pdf(mu-0.001,tau,mu,theta,theta_tilde),ELD.pdf(mu+0.001,tau,mu,theta,theta_tilde))
#    while (count < n) {
#      U1<-runif(1,q.001,q.999)
#      U2<-runif(1,0,C)
#      if(ELD.pdf(y=U1,tau,mu,theta,theta_tilde)>U2) {
#        simul[count+1]<-U1
#        count<-count+1
#      }
#      if (count==n) {return(simul)}
#    }
#  }
#}




#tic()# XXsec pour 100 qd débranché
#tau=0.1;
#theta_tilde=c(1,0.5,0.25);theta_tilde=c(1)
#ELD.inversecdf<-myinverse(ELD.cdf3)

#inverse.all<-rep(NA,9999)
#inverse.fun.all<-inverse(ELD.cdf3)
#for (i in 1:9999) {
#  inverse.all[i]<-inverse.fun.all(i/(10000))
#  if(i%%1000 == 0){print(i)}
#}


#test=matrix(NA,nrow=1000,ncol=3)
#Y=rELD4(100000,tau=0.9,theta = c(1,0.5),theta_tilde = c(1,0.5))
#tau=0.9
#theta=theta_tilde=c(1,0.5)
#for (j in 1:1000) {
#index<-runif(100,1,9999)
#floorindex<-floor(index)
#y=(inverse.all[floorindex]/2 + inverse.all[floorindex+1]/2)
#y=Y[( 1+(100*(j-1)) ): (100*j)]
#quantile((inverse.all[floorindex]/2 + inverse.all[floorindex+1]/2),0.1)
#write.table(inverse.all,"C:/Users/bdeketelaere/Documents/Sauvegardessimulations/Data Simulation 2/mu=0;tau=0,1;theta=c(1,0,5,0,25);theta_tilde=1")
#test[j,]=optim(par = c(0,0.5,0.5),fn=ELD.loglik3)$par
#if(j%%50 == 0){print(j)}
#}
#col_means(test)
#sd(test[,1]);sd(test[,2]);sd(test[,3])



#tic()# 1heure
#tau=0.1;
#theta_tilde=c(1,0.5,0.25);theta_tilde=c(1)
#fullpool=rELD4(2000000,tau=0.1,mu=0,theta=c(1,0.5,0.25),theta_tilde=1)
#quantile(fullpool,0.1)
#plot(density(fullpool))

#write.table(fullpool,"C:/Users/bdeketelaere/Documents/Sauvegardessimulations/Data Simulation 3/mu=0;tau=0,1;theta=c(1,0,5,0,25);theta_tilde=1")
#head(read.table("mu=0;tau=0,1;theta=c(1,0,5,0,25);theta_tilde=1",nrows=1000000))
#toc()



#tic()# 40min
#tau=0.1;
#theta=theta_tilde=c(1,0.5)
#fullpool2=rELD4(2000000,tau=0.1,mu=0,theta=c(1,0.5),theta_tilde=c(1,0.5))
#quantile(fullpool2,0.1)
#plot(density(fullpool2))
#toc()
#write.table(fullpool2,"C:/Users/bdeketelaere/Documents/Sauvegardessimulations/Data Simulation 3/mu=0;tau=0,1;theta=c(1,0,5);theta_tilde=c(1,0,5)")



#tic()
#tau=0.5;
#theta=theta_tilde=c(1,0.5)
#fullpool3=rELD4(2000000,tau=0.5,mu=0,theta=c(1,0.5),theta_tilde=c(1,0.5))
#toc()

#write.table(fullpool3,"C:/Users/bdeketelaere/Documents/Sauvegardessimulations/Data Simulation 3/mu=0;tau=0,5;theta=c(1,0,5);theta_tilde=c(1,0,5)")



#tic()#30 min
#tau=0.5;
#theta=c(1,0.5);theta_tilde=c(1)
#fullpool4=rELD4(2000000,tau=0.5,mu=0,theta=c(1,0.5),theta_tilde=c(1))
#toc()
#quantile(fullpool4,0.5)
#plot(density(fullpool4))

#write.table(fullpool4,"C:/Users/bdeketelaere/Documents/Sauvegardessimulations/Data Simulation 3/mu=0;tau=0,5;theta=c(1,0,5);theta_tilde=c(1)")



#tic()#1heure
#tau=0.9;
#theta=c(1,0.5);theta_tilde=c(1,0.5)
#fullpool5=rELD4(2000000,tau=0.9,mu=0,theta=c(1,0.5),theta_tilde=c(1,0.5))
#toc()


#write.table(fullpool5,"C:/Users/bdeketelaere/Documents/Sauvegardessimulations/Data Simulation 3/mu=0;tau=0,9;theta=c(1,0,5);theta_tilde=c(1,0,5)")



#tic()#1heure
#tau=0.5;
#theta=c(1,1);theta_tilde=c(1,1)
#fullpool6=rELD4(2000000,tau=0.5,mu=0,theta=c(1,1),theta_tilde=c(1,1))
#toc()


#write.table(fullpool6,"C:/Users/bdeketelaere/Documents/Sauvegardessimulations/Data Simulation 3/mu=0;tau=0,5;theta=c(1,1);theta_tilde=c(1,1)")



#tic()#2heure
#tau=0.5;
#theta=c(1,1,1);theta_tilde=c(1,1,1)
#fullpool7=rELD4(2000000,tau=0.5,mu=0,theta=c(1,1,1),theta_tilde=c(1,1,1))
#toc()


#write.table(fullpool7,"C:/Users/bdeketelaere/Documents/Sauvegardessimulations/Data Simulation 3/mu=0;tau=0,5;theta=c(1,1,1);theta_tilde=c(1,1,1)")



#tic()#1heure
#tau=0.1;
#theta=c(1,1,1);theta_tilde=c(1,1,1)
#fullpool8=rELD4(2000000,tau=0.1,mu=0,theta=c(1,1,1),theta_tilde=c(1,1,1))
#quantile(fullpool8,0.1)
#toc()

#write.table(fullpool8,"C:/Users/bdeketelaere/Documents/Sauvegardessimulations/Data Simulation 3/mu=0;tau=0,1;theta=c(1,1,1);theta_tilde=c(1,1,1)")



which.min2 <- function(x, last.index = T, ...){
  if(last.index) max(which(x == min(x, ...))) else which.min(x)
}
# which.min2(c(1,1,1))
# which.min2(c(1,0,1))

which.max2 <- function(x, last.index = T, ...){
  if(last.index) max(which(x == max(x, ...))) else which.max(x)
}
#which.max2(c(1,0,1))



#####################################
#A partir de maintenant introduction d une covariate
ELD.loglik.beta<-function(Y,tau,X,beta,theta,theta_tilde){
  
  n = length(Y)
  
  X = cbind(rep(1,n),X)
  
  #we separe the observations larger than mu from the ones smaller than mu
  indexup = Y>= X %*% beta
  indexdown = Y< X %*% beta
  
  
  Yu<-Y[indexup]
  Yd<-Y[indexdown]
  
  Xu<-X[indexup,]
  Xd<-X[indexdown,]
  
  
  #we compute n+ (and n-) the number of observations larger (smaller) than mu
  nplus<-length(Yu)
  nminus<-length(Yd)
  
  #correspond to the m in our equation
  m <- length(theta)-1
  m_tilde <- length(theta_tilde)-1
  
  #In this part, we compute the Laguerre polynomial :
  # The Laguerre polynomial of order k evaluated in x is given by :
  #\sum_k=0^m (k choose m) (-x)^k / k!
  #We compute this in 3 separated parts. In each part, the code is repeated for
  #observations smaller and larger than mu.
  
  #Part 1 computes the binomial coefficients needed for the computation of the 
  #Laguerre Polynomial. These coefficients are stored in a (lower) triangular
  #matrix with ones on the diagonal. In other words, Part 1 computes k choose m
  #and store this in a matrix called B.
  B <- matrix(0,ncol=m+1,nrow=m+1)
  for(k in 1:(m+1)) {
    B[k,1:k] <- choose(k-1,0:(k-1))
  }
  
  B_tilde <- matrix(0,ncol=m_tilde+1,nrow=m_tilde+1)
  for(k in 1:(m_tilde+1)) {
    B_tilde[k,1:k] <- choose(k-1,0:(k-1))
  }
  
  
  #In part 2, we compute the part " (-x)^k / k!" of the Laguerre polynomial.  
  #We store these results in a matrix with m columns and n lines (where n is the
  #number of observations). So the element of the i-th row and j-th column of Z
  #corresponds to (-x_i ^(j-1))/j-1!. The j-th column is given by the j-1th 
  #column multiplied by -x/j.
  #We also note that the Laguerre polynomial is evauated in tau*(y-mu) for y >mu
  #and in -(1-tau)(y-mu) for y < mu.
  Zu <- matrix(1,ncol=m+1,nrow=nplus)
  if(m+1>=2) {
    for(i in 2:(m+1)) {
      Zu[,i] <- Zu[,i-1]*(-tau*(Yu- (Xu %*% beta)))/(i-1)
    }
  }
  
  Zd <- matrix(1,ncol=m_tilde+1,nrow=nminus)
  if(m_tilde+1>=2) {
    for(i in 2:(m_tilde+1)) {
      Zd[,i] <- Zd[,i-1]*((1-tau)*(Yd-(Xd %*% beta)))/(i-1)
    }
  }
  
  #Part 3 computes the coefficients of the Laguerre Polynomial. These coeffi-
  #cients are stored in a (lower) triangular matrix. The i+1th-row of the matrix
  #contains the k+1 coefficients of i-th order Laguerre polynomial. In this
  #part, the results of part 1 and 2 are put together and the matrician product 
  #computes the sum in the definition of the Laguerre polynomial.
  Lu <- matrix(0,ncol=m+1,nrow=m+1)
  Lu <- Zu%*%t(B)
  
  Ld <- matrix(0,ncol=m_tilde+1,nrow=m_tilde+1)
  Ld <- Zd%*%t(B_tilde)
  
  #Now that the Laguerre polynomial is calculated, we can compute the 
  #log-likelihood by using the formula that we derived
  
  
  return(
    n*log(tau*(1-tau)) +
      -tau*sum(Yu- (Xu %*% beta)) -2*nplus*log(Norm(theta)) + sum(log((Lu%*%theta)^2))
    +(1-tau)*sum(Yd- (Xd %*% beta)) -2*nminus*log(Norm(theta_tilde)) + sum(log((Ld%*%theta_tilde)^2))
  )
}



ELD.loglik.beta2<-function(param2){
  m=length(theta)-1;m_tilde=length(theta_tilde)-1
  
  lenbeta=ncol(as.matrix(X))+1
  beta=param2[1:lenbeta]
  
  if (m != 0 & m_tilde != 0) {
    theta=c(1,param2[(lenbeta+1):(lenbeta+1+m-1)])
    theta_tilde=c(1,param2[(lenbeta+1+m):(lenbeta+1+m+m_tilde-1)])
    return(-1*ELD.loglik.beta(Y,tau,X,beta,theta,theta_tilde))
  } else if (m == 0 & m_tilde == 0) {
    theta=1
    theta_tilde=1
    return(-1*ELD.loglik.beta(Y,tau,X,beta,theta,theta_tilde))
  } else if (m == 0 & m_tilde != 0) {
    theta=1
    theta_tilde=c(1,param2[(lenbeta+1):(lenbeta+m_tilde)])
    return(-1*ELD.loglik.beta(Y,tau,X,beta,theta,theta_tilde))
  } else if (m != 0 & m_tilde == 0) {
    theta_tilde=1
    theta=c(1,param2[(lenbeta+1):(lenbeta+m)])
    return(-1*ELD.loglik.beta(Y,tau,X,beta,theta,theta_tilde))
  }  
  
}

#Passons maintenant au cas avec censure

loglik.IC.beta<-function(beta,theta,theta_tilde,Y,tau,X)
{
  n <- nrow(Y)
  X=cbind(rep(1,n),X)
  #  if(mu > max(Y[,2])){return(-Inf)}
  longY <- c(Y[,1],Y[,2])
  res=rep(NA,2*n)
  for (i in 1:(n-1)) {
    res[c(i,n+i)]<-ELD.cdf2(longY[c(i,n+i)],tau,mu= as.numeric(X[(i%%n),]%*%beta),theta,theta_tilde)
  }
  
  res[c(n,n+n)]<-ELD.cdf2(longY[c(n,n+n)],tau,as.numeric(X[(n),]%*%beta),theta,theta_tilde)
  sum(log(res[(n+1):(2*n)]-res[1:n]))
}




loglik.IC.beta2<-function(param2){
  m=length(theta)-1;m_tilde=length(theta_tilde)-1
  
  lenbeta=ncol(as.matrix(X))+1
  beta=param2[1:lenbeta]
  
  if (m != 0 & m_tilde != 0) {
    theta=c(1,param2[(lenbeta+1):(lenbeta+1+m-1)])
    theta_tilde=c(1,param2[(lenbeta+1+m):(lenbeta+1+m+m_tilde-1)])
    return(-1*loglik.IC.beta(beta,theta,theta_tilde,Y,tau,X))
  } else if (m == 0 & m_tilde == 0) {
    theta=1
    theta_tilde=1
    return(-1*loglik.IC.beta(beta,theta,theta_tilde,Y,tau,X))
  } else if (m == 0 & m_tilde != 0) {
    theta=1
    theta_tilde=c(1,param2[(lenbeta+1):(lenbeta+m_tilde)])
    return(-1*loglik.IC.beta(beta,theta,theta_tilde,Y,tau,X))
  } else if (m != 0 & m_tilde == 0) {
    theta_tilde=1
    theta=c(1,param2[(lenbeta+1):(lenbeta+m)])
    return(-1*loglik.IC.beta(beta,theta,theta_tilde,Y,tau,X))
  }  
  
}




##### On cree une fonction qui n'esat adéquate que quand il n'y a qu'une seule covariate discrete

loglik.IC.beta.1covdis<-function(beta,theta,theta_tilde,Y,tau,X)
{
  n <- nrow(Y)
  index1= X==1
  index0= X==0
  X=cbind(rep(1,n),X)
  Y0=Y[index0,]
  Y1=Y[index1,]
  #  if(mu > max(Y[,2])){return(-Inf)}
  longY0 <- c(Y0[,1],Y0[,2])
  longY1 <- c(Y1[,1],Y1[,2])
  
  res0 = ELD.cdf2(longY0,tau,beta[1],theta,theta_tilde)
  res1 = ELD.cdf2(longY1,tau,beta[1]+beta[2],theta,theta_tilde)
  
  res=rep(NA,2*n)
  res[index0]=res0
  res[index1]=res1
  
  sum(log(res[(n+1):(2*n)]-res[1:n]))
}




loglik.IC.beta.1covdis2<-function(param2){
  m=length(theta)-1;m_tilde=length(theta_tilde)-1
  
  lenbeta=ncol(as.matrix(X))+1
  beta=param2[1:lenbeta]
  
  if (m != 0 & m_tilde != 0) {
    theta=c(1,param2[(lenbeta+1):(lenbeta+1+m-1)])
    theta_tilde=c(1,param2[(lenbeta+1+m):(lenbeta+1+m+m_tilde-1)])
    return(-1*loglik.IC.beta.1covdis(beta,theta,theta_tilde,Y,tau,X))
  } else if (m == 0 & m_tilde == 0) {
    theta=1
    theta_tilde=1
    return(-1*loglik.IC.beta.1covdis(beta,theta,theta_tilde,Y,tau,X))
  } else if (m == 0 & m_tilde != 0) {
    theta=1
    theta_tilde=c(1,param2[(lenbeta+1):(lenbeta+m_tilde)])
    return(-1*loglik.IC.beta.1covdis(beta,theta,theta_tilde,Y,tau,X))
  } else if (m != 0 & m_tilde == 0) {
    theta_tilde=1
    theta=c(1,param2[(lenbeta+1):(lenbeta+m)])
    return(-1*loglik.IC.beta.1covdis(beta,theta,theta_tilde,Y,tau,X))
  }  
  
}


###On passe a 2 covariates discretes

loglik.IC.beta.2covdis<-function(beta,theta,theta_tilde,Y,tau,X)
{
  n <- nrow(Y)
  index11= (X[,1]==1 & X[,2]==1)
  index00= (X[,1]==0 & X[,2]==0)
  index10= (X[,1]==1 & X[,2]==0)
  index01= (X[,1]==0 & X[,2]==1)
  
  X=cbind(rep(1,n),X)
  
  Y00=Y[index00,]
  Y01=Y[index01,]
  Y10=Y[index10,]
  Y11=Y[index11,]
  
  #  if(mu > max(Y[,2])){return(-Inf)}
  longY00 <- c(Y00[,1],Y00[,2])
  longY01 <- c(Y01[,1],Y01[,2])
  longY10 <- c(Y10[,1],Y10[,2])
  longY11 <- c(Y11[,1],Y11[,2])
  
  res00 = ELD.cdf2(longY00,tau,beta[1],theta,theta_tilde)
  res01 = ELD.cdf2(longY01,tau,beta[1]+beta[3],theta,theta_tilde)
  res10 = ELD.cdf2(longY10,tau,beta[1]+beta[2],theta,theta_tilde)
  res11 = ELD.cdf2(longY11,tau,beta[1]+beta[2]+beta[3],theta,theta_tilde)
  
  res=rep(NA,2*n)
  res[index00]=res00
  res[index01]=res01
  res[index10]=res10
  res[index11]=res11
  
  sum(log(res[(n+1):(2*n)]-res[1:n]))
}



loglik.IC.beta.2covdis2<-function(param2){
  m=length(theta)-1;m_tilde=length(theta_tilde)-1
  
  lenbeta=ncol(as.matrix(X))+1
  beta=param2[1:lenbeta]
  
  if (m != 0 & m_tilde != 0) {
    theta=c(1,param2[(lenbeta+1):(lenbeta+1+m-1)])
    theta_tilde=c(1,param2[(lenbeta+1+m):(lenbeta+1+m+m_tilde-1)])
    return(-1*loglik.IC.beta.2covdis(beta,theta,theta_tilde,Y,tau,X))
  } else if (m == 0 & m_tilde == 0) {
    theta=1
    theta_tilde=1
    return(-1*loglik.IC.beta.2covdis(beta,theta,theta_tilde,Y,tau,X))
  } else if (m == 0 & m_tilde != 0) {
    theta=1
    theta_tilde=c(1,param2[(lenbeta+1):(lenbeta+m_tilde)])
    return(-1*loglik.IC.beta.2covdis(beta,theta,theta_tilde,Y,tau,X))
  } else if (m != 0 & m_tilde == 0) {
    theta_tilde=1
    theta=c(1,param2[(lenbeta+1):(lenbeta+m)])
    return(-1*loglik.IC.beta.2covdis(beta,theta,theta_tilde,Y,tau,X))
  }  
  
}

#tau=0.1
#theta=theta_tilde=c(1,1,1,1);
#optim(par=c(qnorm(tau),1.1,1.1,0.5,0,0,0.5,0,0),fn=loglik.IC.beta.2covdis2,control = list(maxit=10000))

#tau=0.1
#theta=c(1,1,1)
#theta_tilde=c(1,1);
#optim(par=c(qnorm(tau),1.1,1.1,0.5,0,0.5),fn=loglik.IC.beta.2covdis2,control = list(maxit=10000))



CVtobestm = function(CV.crit,mmax){
  indmin=which.min(CV.crit)
  if(indmin %% (mmax+1) == 0){
    bestm = mmax
    bestmtilde  = (indmin /(mmax+1))-1
  } else if(indmin %% (mmax+1) != 0) {
    bestm = (indmin %% (mmax+1))-1
    bestmtilde = floor(indmin/(mmax+1))
  }
  return(c(bestm,bestmtilde))
}

# CV.crit = matrix(1:12,nrow=4,ncol=3)
# 
# CVtobestm(CV.crit,3)
# 
# CV.crit[4,1]=0
# CVtobestm(CV.crit,3)
# 
# CV.crit[3,2]=-1
# CVtobestm(CV.crit,3)[1]



################################################################################
##################### Fonctions nécéssaires pour Zhou ##########################
################################################################################


## Implementation of the check function of Zhou
check_ZHOU <- function(Y,tau,beta,X) {
  n=nrow(Y)
  
  res=rep(NA,n)
  
  X=cbind(rep(1,n),X)
  xbeta = X %*% beta
  
  indexdown=Y[,1] >= xbeta
  indexup <- Y[,2] < xbeta
  
  res[indexdown]= tau*abs( Y[,1][indexdown] - xbeta[indexdown] )
  res[indexup]=  (1-tau) * abs( Y[,2][indexup] - xbeta[indexup] )
  res[is.na(res)] = 0
  
  return(res)
}
# 
# tau=0.5
# check_ZHOU(Y,tau,beta,X)
# tau=0.1
# check_ZHOU(Y,tau,beta,X)

optim_ZHOU <- function(beta){
  n=nrow(Y)
  sum(check_ZHOU(Y,tau,beta,X))/n
}

#tau=0.5
#optim(par=c(-1,1),fn=optim_ZHOU)


#BETA0 = seq(from=-2,to=2,by=0.005)
#test=rep(NA,length(BETA0))
#for (i in 1:length(BETA0)) {
#test[i]=optim_ZHOU(c(BETA0[i],1))  
#}
#plot(BETA0,test)


Biascor_ZHOU<-function(Y,tau,beta,X,betahat){
  n=nrow(Y)
  
  res=rep(NA,n)
  
  X=cbind(rep(1,n),X)
  xbeta = X %*% beta
  xbetahat = X %*% betahat
  
  indexdown=Y[,1] >= xbetahat
  indexup <- Y[,2] < xbetahat
  
  res[indexdown]= xbeta[indexdown] * -tau
  res[indexup]=  (1-tau) * xbeta[indexup]
  res[is.na(res)] = 0
  
  return(res) 
}

# tau=0.5
# Biascor_ZHOU(Y,tau,beta,X,betahat)
# tau=0.1
# Biascor_ZHOU(Y,tau,beta=c(0,1),X,betahat=c(0,1))


optim_ZHOU2 <- function(beta){
  n=nrow(Y)
  biascorrection=Biascor_ZHOU(Y,tau,beta,X,betahat)
  sum(check_ZHOU(Y,tau,beta,X)-biascorrection)
}





####################################

simul.CAM.norm<-function(n,tau,sd){
  L=R=rep(NA,n) 
  
  Y=rnorm(n,sd=sd)
  Y=Y-quantile(Y,tau)
  
  
  U=min(Y)-sd + runif(n,min=0,max=sd)
  
  #le 50 peut tout bloquer ici
  lCAM = runif(50,min=0,max=sd)
  lCAM[1] = 0
  
  
  for (i in 1:n) {
    j=1
    for (j in 1:50) {
      if ( !( (U[i] + sum(lCAM[1:j]) < Y[i])  &  (U[i] + sum(lCAM[1:(j+1)]) > Y[i]) ) ) {
        j=j+1
      } else {
        L[i] = U[i] + sum(lCAM[1:j])
        R[i] = U[i] + sum(lCAM[1:(j+1)])
        break
      }
    }  
  }
  
  cbind(L,R)
}

#set.seed(1)
#test=simul.CAM.norm(100,0.1,0.3)
#mean(test[,2])-mean(test[,1])




##### Une variable continue
###Différence entre ELD.cdf2 et ELD.cdf2.cont est que ce dernier admet comme 
###valeur pour mu un vecteur de la meme taille que y. C'est nécéssaire pour 
###l'introduction d'une covariate continue
ELD.cdf2.cont<-function(y,tau,mu,theta,theta_tilde){
  m <- length(theta)-1
  m_tilde <- length(theta_tilde)-1
  n <- length(y)
  
  if(length(mu) == 1){print("plutot utiliser ELD.cdf2") 
  } else if (length(mu) != length(y)){print("erreur de dim dans mu ou dans y");return()}
  
  ################################### D?but partie 1 #############################
  #Cette partie contient la partie o? on calcule tous les coefficients
  #Son temps de calcul n'est pas tres optimis? mais cette partie n'est ?x?cut?e
  #qu'une fois et les matrices sont souvent assez petites). Le code est r?p?t?
  #deux fois 
  Au <- matrix(0,ncol=m+1,nrow=m+1)
  for (k in 1:(m+1)) {
    Au[k,1:k] <- choose(k-1,0:(k-1)) *(-1)^(0:(k-1))/ factorial(0:(k-1))
  }
  
  Ad <- matrix(0,ncol=m_tilde+1,nrow=m_tilde+1)
  for (k in 1:(m_tilde+1)) {
    Ad[k,1:k] <- choose(k-1,0:(k-1)) *(-1)^(0:(k-1))/ factorial(0:(k-1))
  }
  
  
  #B=(j+j')! / l!
  Bu <- matrix(0,ncol=2*m+1,nrow=2*m+1)
  for(j.j_prime in 0:(2*m))
  {
    Bu[j.j_prime+1,(1:(j.j_prime+1))]<-factorial(j.j_prime)/factorial(0:j.j_prime)
  }
  #Bu
  
  
  Bd <- matrix(0,ncol=2*m_tilde+1,nrow=2*m_tilde+1)
  for(j.j_prime in 0:(2*m_tilde))
  {
    Bd[j.j_prime+1,(1:(j.j_prime+1))]<-factorial(j.j_prime)/factorial(0:j.j_prime)
  }
  #Bd
  
  
  #C contient toutes les combinaisons du produit de A j,k avec A jprime,kprime
  Cu <- array(0, dim=c(m+1,m+1,m+1,m+1))
  
  for (k in 0:m) {
    for (k_prime in 0:m) {
      for (j in 0:k) {
        for (j_prime in 0:k_prime) {
          Cu[j+1,j_prime+1,k+1,k_prime+1]<-Au[k+1,j+1] * Au[k_prime+1,j_prime+1]
        }
      }
    }
  }
  
  
  Cd <- array(0, dim=c(m_tilde+1,m_tilde+1,m_tilde+1,m_tilde+1))
  
  for (k in 0:m_tilde) {
    for (k_prime in 0:m_tilde) {
      for (j in 0:k) {
        for (j_prime in 0:k_prime) {
          Cd[j+1,j_prime+1,k+1,k_prime+1]<-Ad[k+1,j+1] * Ad[k_prime+1,j_prime+1]
        }
      }
    }
  }
  
  ###Cetta partie est un peu plus rapide que la pr?c?dente mais il y a un probleme avec la sym?trie. Bizarre :'(
  #  tic()
  
  #    for (j in 0:k) {
  #      for (j_prime in 0:k_prime) {
  #       my.array[j+1,j_prime+1,1:(m+1),1:(m+1)]<-A[(0:m)+1,j+1] * A[1:(m+1),j_prime+1]
  
  #      }
  #   }
  
  #toc()
  #my.array
  
  
  #Maintenant on veut connaitre le coefficient devant le terme en j* du produit
  #de vk et vkprime (avec j* in 0,1,...,k+k'). Pour ca on prend :
  #pour une valeur de k,k' la somme tel des ??ments tels que j+j' = j*
  
  #Dk kprime jpjprime = sum  Ajk*Ajprimekprime tel que j+jprime = jpjprime
  
  #Un as.mtrix a du etre rajout? car visiblement les fonctions row et col  ne 
  #marchent pas quand m et/ou mtilde = 1
  
  Du <- array(0, dim=c(m+1,m+1,2*m+1))
  for (k in 1:(m+1)) {
    for (k_prime in 1:(m+1)) {
      for (j_star in 1:(k+k_prime-1)) {
        Du[k,k_prime,j_star]<-
          sum(Cu[,,k,k_prime][(row(as.matrix(Cu[,,k,k_prime])) +col(as.matrix(Cu[,,k,k_prime])) -1)==j_star])
      }
    }
  }
  
  Dd <- array(0, dim=c(m_tilde+1,m_tilde+1,2*m_tilde+1))
  for (k in 1:(m_tilde+1)) {
    for (k_prime in 1:(m_tilde+1)) {
      for (j_star in 1:(k+k_prime-1)) {
        Dd[k,k_prime,j_star]<-
          sum(Cd[,,k,k_prime][(row(as.matrix(Cd[,,k,k_prime])) +col(as.matrix(Cd[,,k,k_prime])) -1)==j_star])
      }
    }
  }
  
  
  ################################### D?but partie 2 #############################
  
  #C'est un peu embetan de devoir passer par des data frame :( Je fais ca pour
  #pouvoir rendre le vecteur réponse dans l'ordre correspondant aux vecteurs 
  #d'entrées)
  
  data.temp=data.frame(num=1:n,y=y,mu=mu)
  
  data.up=subset(data.temp, y >= mu)
  
  num.up=data.up$num
  
  yu <- data.up$y
  muu <- data.up$mu
  su <- (yu-muu)*tau
  lenu<-length(su)
  resu <- rep(NA,lenu)#ce vecteur contiendra le r?sultat pour les valeurs 
  #demand?es sup?rieures ? ?.
  if(lenu == 0) {resu=numeric(0)} else
    if(lenu != 0) {
      #E contient -(s)^l * exp(-s) pour diff?rentes valeurs de l dans 0,...,2m
      Eu<-matrix(NA,2*m +1,lenu)
      for (l in 1:(2*m+1)) {
        Eu[l,] = (su)^(l-1) * exp(-su) #pas sur du moins ici, ? v?rifier.  
      }
      
      "Gj+jprime,l  contient - j+j!/l! s^l exp(-s) pour toutes les valeurs de l 
  #entre 0 et 2m et de j + j_prime entre 0 et 2m"
      Gu <- array(0,dim=c(2*m+1, 2*m+1,lenu))
      for (j.j_prime in 1:(2*m+1)) {
        for (l in 1:(j.j_prime)) { #-1 cark+kprime sont plus longs de 1
          Gu[j.j_prime,l,]<- Bu[j.j_prime,l] * Eu[l,]
        }
      }
      #Hj+jprime contient - sum l=0^j+jprime j+j!/l! s^l exp(-s)
      Hu<-array(0,dim=c(2*m+1,lenu))
      if(m!=0){
        for (i in 1:lenu) {
          Hu[,i] <- rowSums(Gu[,,i])
        }
      } else if (m==0) {  
        for (i in 1:lenu) {
          Hu[,i] <- Gu[,,i]
        }}
      #J contient 
      Ju<- array(0,dim=c(m+1, m+1, lenu))
      for (k in 1:(m+1)) {
        for (k_prime in 1:(m+1)) {
          Ju[k,k_prime,]<- (Du[k,k_prime,] %*% Hu)
        }
      }
      THETA <- (theta) %*% t(theta)
      
      #Le probleme avec J C'est que c'est difficile de faire des sommes 
      #sans utiliser de boucles for sur y. Ainsi, on transforme l'array en une
      #matrice avec m+1^2 colonnes et lenu lignes et on fait les operations 
      #necessaires sur cette matrice K
      Ku <- matrix(Ju[,,1:lenu],ncol = (m+1)^2,nrow=lenu,byrow=T)*matrix(rep(THETA),byrow=T,ncol = (m+1)^2,nrow=lenu)
      
      if(m!=0){
        resu <- 1 + (1 -tau)* -rowSums(Ku)/Norm(theta)^2 
      } else if (m==0) {
        resu <- 1 + ((1-tau)* -Ku/Norm(theta)^2) 
      }
      
    }
  
  
  ##############################################################################
  #Pareil pour les observations < mu
  
  data.down=subset(data.temp, y < mu)
  
  num.down=data.down$num
  
  yd <- data.down$y
  mud <- data.down$mu
  sd <- -(yd-mud)*(1-tau)
  lend<-length(sd)
  if(lend == 0) {resd=numeric(0)} else
    if(lend != 0){
      
      resd<-rep(NA,lend)
      
      Ed<-matrix(NA,2*m_tilde +1,lend)
      for (l in 1:(2*m_tilde+1)) {
        Ed[l,] = (sd)^(l-1) * exp(-sd) #pas sur du moins ici, ? v?rifier.  
      }
      
      Gd <- array(0,dim=c(2*m_tilde+1, 2*m_tilde+1,lend))
      for (j.j_prime in 1:(2*m_tilde+1)) {
        for (l in 1:(j.j_prime)) { #-1 cark+kprime sont plus longs de 1
          Gd[j.j_prime,l,]<- Bd[j.j_prime,l] * Ed[l,]
        }
      }
      
      Hd<-array(0,dim=c(2*m_tilde+1,lend))
      if(m_tilde!=0){
        for (i in 1:lend) {
          Hd[,i] <- rowSums(Gd[,,i])
        }
      } else if (m_tilde==0) {  
        for (i in 1:lend) {
          Hd[,i] <- Gd[,,i]
        }}
      
      Jd<- array(0,dim=c(m_tilde+1, m_tilde+1, lend))
      for (k in 1:(m_tilde+1)) {
        for (k_prime in 1:(m_tilde+1)) {
          Jd[k,k_prime,]<- (Dd[k,k_prime,] %*% Hd)
        }
      }
      THETA_tilde <- (theta_tilde) %*% t(theta_tilde)
      
      Kd <- matrix(Jd[,,1:lend],ncol = (m_tilde+1)^2,nrow=lend,byrow=T)*matrix(rep(THETA_tilde),byrow=T,ncol = (m_tilde+1)^2,nrow=lend)
      
      if(m_tilde!=0){
        resd <- -tau* -rowSums(Kd)/Norm(theta_tilde)^2 
      } else if (m_tilde==0) {
        resd <- -tau* -Kd/Norm(theta_tilde)^2 
      }
      
    }
  
  data.temp2<-data.frame(num=c(num.down,num.up),res=c(resd,resu))
  data.temp2[order(data.temp2$num),]$res
  
}



loglik.IC.beta.1covcont<-function(beta,theta,theta_tilde,Y,tau,X)
{
  n <- nrow(Y)
  X=cbind(rep(1,n),X)
  #  if(mu > max(Y[,2])){return(-Inf)}
  longY <- c(Y[,1],Y[,2])
  
  res = ELD.cdf2.cont(longY,tau,rep(X%*%beta,2),theta,theta_tilde)
  
  sum(log(res[(n+1):(2*n)]-res[1:n]))
}

#loglik.IC.beta.1covcont(beta,theta,theta_tilde,Y,tau,X)
#loglik.IC.beta.1covcont(beta=c(0,0),theta,theta_tilde,Y,tau,X)
#loglik.IC.beta.1covcont(beta=c(1,1),theta,theta_tilde,Y,tau,X)
#loglik.IC.beta.1covcont(beta=c(0,0.9),theta,theta_tilde,Y,tau,X)
#loglik.IC.beta.1covcont(beta=c(0,1.1),theta,theta_tilde,Y,tau,X)


loglik.IC.beta.1covcont2<-function(param2){
  m=length(theta)-1;m_tilde=length(theta_tilde)-1
  
  lenbeta=ncol(as.matrix(X))+1
  beta=param2[1:lenbeta]
  
  if (m != 0 & m_tilde != 0) {
    theta=c(1,param2[(lenbeta+1):(lenbeta+1+m-1)])
    theta_tilde=c(1,param2[(lenbeta+1+m):(lenbeta+1+m+m_tilde-1)])
    return(-1*loglik.IC.beta.1covcont(beta,theta,theta_tilde,Y,tau,X))
  } else if (m == 0 & m_tilde == 0) {
    theta=1
    theta_tilde=1
    return(-1*loglik.IC.beta.1covcont(beta,theta,theta_tilde,Y,tau,X))
  } else if (m == 0 & m_tilde != 0) {
    theta=1
    theta_tilde=c(1,param2[(lenbeta+1):(lenbeta+m_tilde)])
    return(-1*loglik.IC.beta.1covcont(beta,theta,theta_tilde,Y,tau,X))
  } else if (m != 0 & m_tilde == 0) {
    theta_tilde=1
    theta=c(1,param2[(lenbeta+1):(lenbeta+m)])
    return(-1*loglik.IC.beta.1covcont(beta,theta,theta_tilde,Y,tau,X))
  }  
  
}

#loglik.IC.beta.1covcont2(c(0,1,1,1))
#loglik.IC.beta.1covcont2(c(0,1,0.5,0.5))

#optim(par=c(0,1,1,1,1,1),fn=loglik.IC.beta.1covcont2,control = list(maxit=10000))





#En fait on peut ajouter d'autres covariate sas rien changer au code.
#On essaie ci-dessous en ajoutant une deuxième covariate discrete a la continue


loglik.IC.beta.general<-function(beta,theta,theta_tilde,Y,tau,X)
{
  n <- nrow(Y)
  X=cbind(rep(1,n),X)
  
  longY <- c(Y[,1],Y[,2])
  
  res = ELD.cdf2.cont(longY,tau,rep(X%*%beta,2),theta,theta_tilde)
  
  sum(log(res[(n+1):(2*n)]-res[1:n]))
}

#loglik.IC.beta.general(beta,theta,theta_tilde,Y,tau,X)
#loglik.IC.beta.general(beta=c(0,1,0),theta,theta_tilde,Y,tau,X)
#loglik.IC.beta.general(beta=c(0,0,1),theta,theta_tilde,Y,tau,X)
#loglik.IC.beta.general(beta=c(0,0,0),theta,theta_tilde,Y,tau,X)

loglik.IC.beta.general2<-function(param2){
  m=length(theta)-1;m_tilde=length(theta_tilde)-1
  
  lenbeta=ncol(as.matrix(X))+1
  beta=param2[1:lenbeta]
  
  if (m != 0 & m_tilde != 0) {
    theta=c(1,param2[(lenbeta+1):(lenbeta+1+m-1)])
    theta_tilde=c(1,param2[(lenbeta+1+m):(lenbeta+1+m+m_tilde-1)])
    return(-1*loglik.IC.beta.general(beta,theta,theta_tilde,Y,tau,X))
  } else if (m == 0 & m_tilde == 0) {
    theta=1
    theta_tilde=1
    return(-1*loglik.IC.beta.general(beta,theta,theta_tilde,Y,tau,X))
  } else if (m == 0 & m_tilde != 0) {
    theta=1
    theta_tilde=c(1,param2[(lenbeta+1):(lenbeta+m_tilde)])
    return(-1*loglik.IC.beta.general(beta,theta,theta_tilde,Y,tau,X))
  } else if (m != 0 & m_tilde == 0) {
    theta_tilde=1
    theta=c(1,param2[(lenbeta+1):(lenbeta+m)])
    return(-1*loglik.IC.beta.general(beta,theta,theta_tilde,Y,tau,X))
  }  
  
}


#loglik.IC.beta.general2(c(0,1,1,1,1,1,1))
#optim(par=c(0,1,1,1,1,1,1),fn=loglik.IC.beta.general2,control = list(maxit=10000))


################################################################################
################################################################################
################################################################################
##################### A partir d'ici on ajoute sigma ###########################
################################################################################
################################################################################
################################################################################

ELD.pdf.sigma<-function(y,tau,mu,sigma,theta,theta_tilde){
  
  #we separe the observations larger than mu from the ones smaller than mu
  yu<-y[y>=mu]
  yd<-y[y<mu]
  
  #we compute n+ (and n-) the number of observations larger (smaller) than mu
  nplus<-length(yu)
  nminus<-length(yd)
  
  n=nplus+nminus
  
  #correspond to the m in our equation
  m <- length(theta)-1
  m_tilde <- length(theta_tilde)-1
  
  #In this part, we compute the Laguerre polynomial :
  # The Laguerre polynomial of order k evaluated in x is given by :
  #\sum_k=0^m (k choose m) (-x)^k / k!
  #We compute this in 3 separated parts. In each part, the code is repeated for
  #observations smaller and larger than mu.
  
  #Part 1 computes the binomial coefficients needed for the computation of the 
  #Laguerre Polynomial. These coefficients are stored in a (lower) triangular
  #matrix with ones on the diagonal. In other words, Part 1 computes k choose m
  #and store this in a matrix called B.
  B <- matrix(0,ncol=m+1,nrow=m+1)
  for(k in 1:(m+1)) {
    B[k,1:k] <- choose(k-1,0:(k-1))
  }
  
  B_tilde <- matrix(0,ncol=m_tilde+1,nrow=m_tilde+1)
  for(k in 1:(m_tilde+1)) {
    B_tilde[k,1:k] <- choose(k-1,0:(k-1))
  }
  
  
  #In part 2, we compute the part " (-x)^k / k!" of the Laguerre polynomial.  
  #We store these results in a matrix with m columns and n lines (where n is the
  #number of observations). So the element of the i-th row and j-th column of Z
  #corresponds to (-x_i ^(j-1))/j-1!. The j-th column is given by the j-1th 
  #column multiplied by -x/j.
  #We also note that the Laguerre polynomial is evauated in tau*(y-mu) for y >mu
  #and in -(1-tau)(y-mu) for y < mu.
  Zu <- matrix(1,ncol=m+1,nrow=length(yu))
  if(m+1>=2) {
    for(i in 2:(m+1)) {
      Zu[,i] <- Zu[,i-1]*(-tau*(yu-mu)/sigma)/(i-1)
    }
  }
  
  Zd <- matrix(1,ncol=m_tilde+1,nrow=length(yd))
  if(m_tilde+1>=2) {
    for(i in 2:(m_tilde+1)) {
      Zd[,i] <- Zd[,i-1]*((1-tau)*(yd-mu)/sigma)/(i-1)
    }
  }
  
  #Part 3 computes the coefficients of the Laguerre Polynomial. These coeffi-
  #cients are stored in a (lower) triangular matrix. The i+1th-row of the matrix
  #contains the k+1 coefficients of i-th order Laguerre polynomial. In this
  #part, the results of part 1 and 2 are put together and the matrician product 
  #computes the sum in the definition of the Laguerre polynomial.
  Lu <- matrix(0,ncol=m+1,nrow=m+1)
  Lu <- Zu%*%t(B)
  
  Ld <- matrix(0,ncol=m_tilde+1,nrow=m_tilde+1)
  Ld <- Zd%*%t(B_tilde)
  
  
  #Now that the Laguerre polynomial is computed, we can construct the pdf of the
  #Enriched Laplace distribution.
  
  f<-rep(0,n)
  f[y>=mu]<-exp(-tau*(yu-mu)/sigma)/Norm(theta)^2*(Lu%*%theta)^2
  f[y<mu]<-exp((1-tau)*(yd-mu)/sigma)/Norm(theta_tilde)^2*(Ld%*%theta_tilde)^2
  f<-f*(1-tau)*tau/sigma
  return(f)
}

#x1=seq(from=-10,to=30,by=0.01)
#y1=ELD.pdf.sigma(y=x1,mu=0,tau=0.25,sigma=1,theta=c(1,1,1),theta_tilde = 1)
#y2=ELD.pdf.sigma(y=x1,mu=0,tau=0.25,sigma=1.25,theta=c(1,1,1),theta_tilde = 1)
#plot(x1,y1,type="l",lwd=1)
#lines(x1,y2,col=2)




ELD.cdf2.sigma<-function(y,tau,mu,sigma,theta,theta_tilde){
  m <- length(theta)-1
  m_tilde <- length(theta_tilde)-1
  
  ################################### D?but partie 1 #############################
  #Cette partie contient la partie o? on calcule tous les coefficients
  #Son temps de calcul n'est pas tres optimis? mais cette partie n'est ?x?cut?e
  #qu'une fois et les matrices sont souvent assez petites). Le code est r?p?t?
  #deux fois 
  Au <- matrix(0,ncol=m+1,nrow=m+1)
  for (k in 1:(m+1)) {
    Au[k,1:k] <- choose(k-1,0:(k-1)) *(-1)^(0:(k-1))/ factorial(0:(k-1))
  }
  
  Ad <- matrix(0,ncol=m_tilde+1,nrow=m_tilde+1)
  for (k in 1:(m_tilde+1)) {
    Ad[k,1:k] <- choose(k-1,0:(k-1)) *(-1)^(0:(k-1))/ factorial(0:(k-1))
  }
  
  
  #B=(j+j')! / l!
  Bu <- matrix(0,ncol=2*m+1,nrow=2*m+1)
  for(j.j_prime in 0:(2*m))
  {
    Bu[j.j_prime+1,(1:(j.j_prime+1))]<-factorial(j.j_prime)/factorial(0:j.j_prime)
  }
  #Bu
  
  
  Bd <- matrix(0,ncol=2*m_tilde+1,nrow=2*m_tilde+1)
  for(j.j_prime in 0:(2*m_tilde))
  {
    Bd[j.j_prime+1,(1:(j.j_prime+1))]<-factorial(j.j_prime)/factorial(0:j.j_prime)
  }
  #Bd
  
  
  #C contient toutes les combinaisons du produit de A j,k avec A jprime,kprime
  Cu <- array(0, dim=c(m+1,m+1,m+1,m+1))
  
  for (k in 0:m) {
    for (k_prime in 0:m) {
      for (j in 0:k) {
        for (j_prime in 0:k_prime) {
          Cu[j+1,j_prime+1,k+1,k_prime+1]<-Au[k+1,j+1] * Au[k_prime+1,j_prime+1]
        }
      }
    }
  }
  
  
  Cd <- array(0, dim=c(m_tilde+1,m_tilde+1,m_tilde+1,m_tilde+1))
  
  for (k in 0:m_tilde) {
    for (k_prime in 0:m_tilde) {
      for (j in 0:k) {
        for (j_prime in 0:k_prime) {
          Cd[j+1,j_prime+1,k+1,k_prime+1]<-Ad[k+1,j+1] * Ad[k_prime+1,j_prime+1]
        }
      }
    }
  }
  
  ###Cetta partie est un peu plus rapide que la pr?c?dente mais il y a un probleme avec la sym?trie. Bizarre :'(
  #  tic()
  
  #    for (j in 0:k) {
  #      for (j_prime in 0:k_prime) {
  #       my.array[j+1,j_prime+1,1:(m+1),1:(m+1)]<-A[(0:m)+1,j+1] * A[1:(m+1),j_prime+1]
  
  #      }
  #   }
  
  #toc()
  #my.array
  
  
  #Maintenant on veut connaitre le coefficient devant le terme en j* du produit
  #de vk et vkprime (avec j* in 0,1,...,k+k'). Pour ca on prend :
  #pour une valeur de k,k' la somme tel des ??ments tels que j+j' = j*
  
  #Dk kprime jpjprime = sum  Ajk*Ajprimekprime tel que j+jprime = jpjprime
  
  #Un as.mtrix a du etre rajout? car visiblement les fonctions row et col  ne 
  #marchent pas quand m et/ou mtilde = 1
  
  Du <- array(0, dim=c(m+1,m+1,2*m+1))
  for (k in 1:(m+1)) {
    for (k_prime in 1:(m+1)) {
      for (j_star in 1:(k+k_prime-1)) {
        Du[k,k_prime,j_star]<-
          sum(Cu[,,k,k_prime][(row(as.matrix(Cu[,,k,k_prime])) +col(as.matrix(Cu[,,k,k_prime])) -1)==j_star])
      }
    }
  }
  
  Dd <- array(0, dim=c(m_tilde+1,m_tilde+1,2*m_tilde+1))
  for (k in 1:(m_tilde+1)) {
    for (k_prime in 1:(m_tilde+1)) {
      for (j_star in 1:(k+k_prime-1)) {
        Dd[k,k_prime,j_star]<-
          sum(Cd[,,k,k_prime][(row(as.matrix(Cd[,,k,k_prime])) +col(as.matrix(Cd[,,k,k_prime])) -1)==j_star])
      }
    }
  }
  
  
  ################################### D?but partie 2 #############################
  
  index<-rank(y)
  y<-sort(y)
  
  yu <- y[y>=mu]
  su <- (yu-mu)*tau/sigma
  lenu<-length(su)
  resu <- rep(NA,lenu)#ce vecteur contiendra le r?sultat pour les valeurs 
  #demand?es sup?rieures ? ?.
  if(lenu == 0) {resu=numeric(0)} else
    if(lenu != 0) {
      #E contient -(s)^l * exp(-s) pour diff?rentes valeurs de l dans 0,...,2m
      Eu<-matrix(NA,2*m +1,lenu)
      for (l in 1:(2*m+1)) {
        Eu[l,] = (su)^(l-1) * exp(-su) #pas sur du moins ici, ? v?rifier.  
      }
      
      "Gj+jprime,l  contient - j+j!/l! s^l exp(-s) pour toutes les valeurs de l 
  #entre 0 et 2m et de j + j_prime entre 0 et 2m"
      Gu <- array(0,dim=c(2*m+1, 2*m+1,lenu))
      for (j.j_prime in 1:(2*m+1)) {
        for (l in 1:(j.j_prime)) { #-1 cark+kprime sont plus longs de 1
          Gu[j.j_prime,l,]<- Bu[j.j_prime,l] * Eu[l,]
        }
      }
      #Hj+jprime contient - sum l=0^j+jprime j+j!/l! s^l exp(-s)
      Hu<-array(0,dim=c(2*m+1,lenu))
      if(m!=0){
        for (i in 1:lenu) {
          Hu[,i] <- rowSums(Gu[,,i])
        }
      } else if (m==0) {  
        for (i in 1:lenu) {
          Hu[,i] <- Gu[,,i]
        }}
      #J contient 
      Ju<- array(0,dim=c(m+1, m+1, lenu))
      for (k in 1:(m+1)) {
        for (k_prime in 1:(m+1)) {
          Ju[k,k_prime,]<- (Du[k,k_prime,] %*% Hu)
        }
      }
      THETA <- (theta) %*% t(theta)
      
      #Le probleme avec J C'est que c'est difficile de faire des sommes 
      #sans utiliser de boucles for sur y. Ainsi, on transforme l'array en une
      #matrice avec m+1^2 colonnes et lenu lignes et on fait le sop?rations 
      #n?c?ssaires sur cette matrice K
      Ku <- matrix(Ju[,,1:lenu],ncol = (m+1)^2,nrow=lenu,byrow=T)*matrix(rep(THETA),byrow=T,ncol = (m+1)^2,nrow=lenu)
      
      if(m!=0){
        resu <- 1 + (1 -tau)* -rowSums(Ku)/Norm(theta)^2 
      } else if (m==0) {
        resu <- 1 + ((1-tau)* -Ku/Norm(theta)^2) 
      }
      
    }
  
  
  ##############################################################################
  #Pareil pour les observations < mu
  yd <- y[y<mu]
  sd <- -(yd-mu)*(1-tau)/sigma
  lend<-length(sd)
  if(lend == 0) {resd=numeric(0)} else
    if(lend != 0){
      
      resd<-rep(NA,lend)
      
      Ed<-matrix(NA,2*m_tilde +1,lend)
      for (l in 1:(2*m_tilde+1)) {
        Ed[l,] = (sd)^(l-1) * exp(-sd) #pas sur du moins ici, ? v?rifier.  
      }
      
      Gd <- array(0,dim=c(2*m_tilde+1, 2*m_tilde+1,lend))
      for (j.j_prime in 1:(2*m_tilde+1)) {
        for (l in 1:(j.j_prime)) { #-1 cark+kprime sont plus longs de 1
          Gd[j.j_prime,l,]<- Bd[j.j_prime,l] * Ed[l,]
        }
      }
      
      Hd<-array(0,dim=c(2*m_tilde+1,lend))
      if(m_tilde!=0){
        for (i in 1:lend) {
          Hd[,i] <- rowSums(Gd[,,i])
        }
      } else if (m_tilde==0) {  
        for (i in 1:lend) {
          Hd[,i] <- Gd[,,i]
        }}
      
      Jd<- array(0,dim=c(m_tilde+1, m_tilde+1, lend))
      for (k in 1:(m_tilde+1)) {
        for (k_prime in 1:(m_tilde+1)) {
          Jd[k,k_prime,]<- (Dd[k,k_prime,] %*% Hd)
        }
      }
      THETA_tilde <- (theta_tilde) %*% t(theta_tilde)
      
      Kd <- matrix(Jd[,,1:lend],ncol = (m_tilde+1)^2,nrow=lend,byrow=T)*matrix(rep(THETA_tilde),byrow=T,ncol = (m_tilde+1)^2,nrow=lend)
      
      if(m_tilde!=0){
        resd <- -tau* -rowSums(Kd)/Norm(theta_tilde)^2 
      } else if (m_tilde==0) {
        resd <- -tau* -Kd/Norm(theta_tilde)^2 
      }
      
    }
  res<-c(resd,resu)
  indexinfinite = is.nan(res) #Note : ces deux lignes ci sont présentes dans le 
  #cas où on introduit de la censure à droite. Si on a de la censure à droite,
  #le bord droit de l'intervalle est +infini. A l'étape du calcul de Eu cela 
  #génère des NaN car pour R exp(-infini) * infini^a (avec a >0) ne vaut pas 0 
  #mais vaut NaN. Note 2 cette modification ne convient pas si le bord gauche
  res[indexinfinite] = 1 #d'un intervalle est -infini
  res[index]
  
}






ELD.cdf2.cont.sigma<-function(y,tau,mu,sigma,theta,theta_tilde){
  m <- length(theta)-1
  m_tilde <- length(theta_tilde)-1
  n <- length(y)
  
  if(length(mu) == 1){print("plutot utiliser ELD.cdf2.sigma") 
  } else if (length(mu) != length(y)){print("erreur de dim dans mu ou dans y");return()}
  
  #pour gérer quand on a des infinis dans les données
  yfinite=y[is.finite(y)]
  yinfinite=y[!is.finite(y)]
  
  y=yfinite
  n=length(y)
  ################################### D?but partie 1 #############################
  #Cette partie contient la partie o? on calcule tous les coefficients
  #Son temps de calcul n'est pas tres optimis? mais cette partie n'est ?x?cut?e
  #qu'une fois et les matrices sont souvent assez petites). Le code est r?p?t?
  #deux fois 
  Au <- matrix(0,ncol=m+1,nrow=m+1)
  for (k in 1:(m+1)) {
    Au[k,1:k] <- choose(k-1,0:(k-1)) *(-1)^(0:(k-1))/ factorial(0:(k-1))
  }
  
  Ad <- matrix(0,ncol=m_tilde+1,nrow=m_tilde+1)
  for (k in 1:(m_tilde+1)) {
    Ad[k,1:k] <- choose(k-1,0:(k-1)) *(-1)^(0:(k-1))/ factorial(0:(k-1))
  }
  
  
  #B=(j+j')! / l!
  Bu <- matrix(0,ncol=2*m+1,nrow=2*m+1)
  for(j.j_prime in 0:(2*m))
  {
    Bu[j.j_prime+1,(1:(j.j_prime+1))]<-factorial(j.j_prime)/factorial(0:j.j_prime)
  }
  #Bu
  
  
  Bd <- matrix(0,ncol=2*m_tilde+1,nrow=2*m_tilde+1)
  for(j.j_prime in 0:(2*m_tilde))
  {
    Bd[j.j_prime+1,(1:(j.j_prime+1))]<-factorial(j.j_prime)/factorial(0:j.j_prime)
  }
  #Bd
  
  
  #C contient toutes les combinaisons du produit de A j,k avec A jprime,kprime
  Cu <- array(0, dim=c(m+1,m+1,m+1,m+1))
  
  for (k in 0:m) {
    for (k_prime in 0:m) {
      for (j in 0:k) {
        for (j_prime in 0:k_prime) {
          Cu[j+1,j_prime+1,k+1,k_prime+1]<-Au[k+1,j+1] * Au[k_prime+1,j_prime+1]
        }
      }
    }
  }
  
  
  Cd <- array(0, dim=c(m_tilde+1,m_tilde+1,m_tilde+1,m_tilde+1))
  
  for (k in 0:m_tilde) {
    for (k_prime in 0:m_tilde) {
      for (j in 0:k) {
        for (j_prime in 0:k_prime) {
          Cd[j+1,j_prime+1,k+1,k_prime+1]<-Ad[k+1,j+1] * Ad[k_prime+1,j_prime+1]
        }
      }
    }
  }
  
  ###Cetta partie est un peu plus rapide que la pr?c?dente mais il y a un probleme avec la sym?trie. Bizarre :'(
  #  tic()
  
  #    for (j in 0:k) {
  #      for (j_prime in 0:k_prime) {
  #       my.array[j+1,j_prime+1,1:(m+1),1:(m+1)]<-A[(0:m)+1,j+1] * A[1:(m+1),j_prime+1]
  
  #      }
  #   }
  
  #toc()
  #my.array
  
  
  #Maintenant on veut connaitre le coefficient devant le terme en j* du produit
  #de vk et vkprime (avec j* in 0,1,...,k+k'). Pour ca on prend :
  #pour une valeur de k,k' la somme tel des ??ments tels que j+j' = j*
  
  #Dk kprime jpjprime = sum  Ajk*Ajprimekprime tel que j+jprime = jpjprime
  
  #Un as.mtrix a du etre rajout? car visiblement les fonctions row et col  ne 
  #marchent pas quand m et/ou mtilde = 1
  
  Du <- array(0, dim=c(m+1,m+1,2*m+1))
  for (k in 1:(m+1)) {
    for (k_prime in 1:(m+1)) {
      for (j_star in 1:(k+k_prime-1)) {
        Du[k,k_prime,j_star]<-
          sum(Cu[,,k,k_prime][(row(as.matrix(Cu[,,k,k_prime])) +col(as.matrix(Cu[,,k,k_prime])) -1)==j_star])
      }
    }
  }
  
  Dd <- array(0, dim=c(m_tilde+1,m_tilde+1,2*m_tilde+1))
  for (k in 1:(m_tilde+1)) {
    for (k_prime in 1:(m_tilde+1)) {
      for (j_star in 1:(k+k_prime-1)) {
        Dd[k,k_prime,j_star]<-
          sum(Cd[,,k,k_prime][(row(as.matrix(Cd[,,k,k_prime])) +col(as.matrix(Cd[,,k,k_prime])) -1)==j_star])
      }
    }
  }
  
  
  ################################### D?but partie 2 #############################
  
  #C'est un peu embetant de devoir passer par des data frame :( Je fais ca pour
  #pouvoir rendre le vecteur réponse dans l'ordre correspondant aux vecteurs 
  #d'entrées)
  
  data.temp=data.frame(num=1:n,y=y,mu=mu[1:length(y)])
  
  data.up=subset(data.temp, y >= mu[1:length(y)])
  
  num.up=data.up$num
  
  yu <- data.up$y
  muu <- data.up$mu
  su <- (yu-muu)*tau/sigma
  lenu<-length(su)
  resu <- rep(NA,lenu)#ce vecteur contiendra le r?sultat pour les valeurs 
  #demand?es sup?rieures ? ?.
  if(lenu == 0) {resu=numeric(0)} else
    if(lenu != 0) {
      #E contient -(s)^l * exp(-s) pour diff?rentes valeurs de l dans 0,...,2m
      Eu<-matrix(NA,2*m +1,lenu)
      for (l in 1:(2*m+1)) {
        Eu[l,] = (su)^(l-1) * exp(-su) #pas sur du moins ici, ? v?rifier.  
      }
      
      "Gj+jprime,l  contient - j+j!/l! s^l exp(-s) pour toutes les valeurs de l 
  #entre 0 et 2m et de j + j_prime entre 0 et 2m"
      Gu <- array(0,dim=c(2*m+1, 2*m+1,lenu))
      for (j.j_prime in 1:(2*m+1)) {
        for (l in 1:(j.j_prime)) { #-1 cark+kprime sont plus longs de 1
          Gu[j.j_prime,l,]<- Bu[j.j_prime,l] * Eu[l,]
        }
      }
      #Hj+jprime contient - sum l=0^j+jprime j+j!/l! s^l exp(-s)
      Hu<-array(0,dim=c(2*m+1,lenu))
      if(m!=0){
        for (i in 1:lenu) {
          Hu[,i] <- rowSums(Gu[,,i])
        }
      } else if (m==0) {  
        for (i in 1:lenu) {
          Hu[,i] <- Gu[,,i]
        }}
      #J contient 
      Ju<- array(0,dim=c(m+1, m+1, lenu))
      for (k in 1:(m+1)) {
        for (k_prime in 1:(m+1)) {
          Ju[k,k_prime,]<- (Du[k,k_prime,] %*% Hu)
        }
      }
      THETA <- (theta) %*% t(theta)
      
      #Le probleme avec J C'est que c'est difficile de faire des sommes 
      #sans utiliser de boucles for sur y. Ainsi, on transforme l'array en une
      #matrice avec m+1^2 colonnes et lenu lignes et on fait les operations 
      #necessaires sur cette matrice K
      Ku <- matrix(Ju[,,1:lenu],ncol = (m+1)^2,nrow=lenu,byrow=T)*matrix(rep(THETA),byrow=T,ncol = (m+1)^2,nrow=lenu)
      
      if(m!=0){
        resu <- 1 + (1 -tau)* -rowSums(Ku)/Norm(theta)^2 
      } else if (m==0) {
        resu <- 1 + ((1-tau)* -Ku/Norm(theta)^2) 
      }
      
    }
  
  
  ##############################################################################
  #Pareil pour les observations < mu
  
  data.down=subset(data.temp, y < mu[1:length(y)])
  
  num.down=data.down$num
  
  yd <- data.down$y
  mud <- data.down$mu
  sd <- -(yd-mud)*(1-tau)/sigma
  lend<-length(sd)
  if(lend == 0) {resd=numeric(0)} else
    if(lend != 0){
      
      resd<-rep(NA,lend)
      
      Ed<-matrix(NA,2*m_tilde +1,lend)
      for (l in 1:(2*m_tilde+1)) {
        Ed[l,] = (sd)^(l-1) * exp(-sd) #pas sur du moins ici, ? v?rifier.  
      }
      
      Gd <- array(0,dim=c(2*m_tilde+1, 2*m_tilde+1,lend))
      for (j.j_prime in 1:(2*m_tilde+1)) {
        for (l in 1:(j.j_prime)) { #-1 cark+kprime sont plus longs de 1
          Gd[j.j_prime,l,]<- Bd[j.j_prime,l] * Ed[l,]
        }
      }
      
      Hd<-array(0,dim=c(2*m_tilde+1,lend))
      if(m_tilde!=0){
        for (i in 1:lend) {
          Hd[,i] <- rowSums(Gd[,,i])
        }
      } else if (m_tilde==0) {  
        for (i in 1:lend) {
          Hd[,i] <- Gd[,,i]
        }}
      
      Jd<- array(0,dim=c(m_tilde+1, m_tilde+1, lend))
      for (k in 1:(m_tilde+1)) {
        for (k_prime in 1:(m_tilde+1)) {
          Jd[k,k_prime,]<- (Dd[k,k_prime,] %*% Hd)
        }
      }
      THETA_tilde <- (theta_tilde) %*% t(theta_tilde)
      
      Kd <- matrix(Jd[,,1:lend],ncol = (m_tilde+1)^2,nrow=lend,byrow=T)*matrix(rep(THETA_tilde),byrow=T,ncol = (m_tilde+1)^2,nrow=lend)
      
      if(m_tilde!=0){
        resd <- -tau* -rowSums(Kd)/Norm(theta_tilde)^2 
      } else if (m_tilde==0) {
        resd <- -tau* -Kd/Norm(theta_tilde)^2 
      }
      
    }
  
  data.temp2<-data.frame(num=c(num.down,num.up),res=c(resd,resu))
  c(data.temp2[order(data.temp2$num),]$res,rep(1,length(yinfinite)))#on ajoute à la fin du vecteur des 1 pour chaque foi qu'il y a un infini. Ca va fonctionner que si on a des inf uniquement à la fin du vecteur
  
}

#x1=seq(from=-10,to=30,by=0.01)
#y1=ELD.cdf2.cont.sigma(y=x1,mu=0,tau=0.25,sigma=1,theta=c(1,1,1),theta_tilde = 1)
#y2=ELD.cdf2.cont.sigma(y=x1,mu=0,tau=0.25,sigma=2,theta=c(1,1,1),theta_tilde = 1)
#plot(x1,y1,type="l",lwd=1)
#lines(x1,y2,col=2)




# beta=c(0,1,1);
# sigma=2
# theta=theta_tilde=c(1,1,1);
# eps=rnorm(100);
# tau=0.5;
# X1=runif(100)
# X2=rbinom(100,size=1,p=0.5)
# Y=1*X1 + 1*X2 +eps
# L=Y-1*runif(100)
# R=Y+1*runif(100)
# Y=cbind(L,R)
# X=cbind(X1,X2)

loglik.IC.beta.general.sigma<-function(beta, sigma, theta,theta_tilde,Y,tau,X)
{
  n <- nrow(Y)
  X=cbind(rep(1,n),X)
  
  longY <- c(Y[,1],Y[,2])
  
  res = ELD.cdf2.cont.sigma(longY,tau,rep(X%*%beta,2),sigma,theta,theta_tilde)
  
  sum(log(res[(n+1):(2*n)]-res[1:n]))
}

# loglik.IC.beta.general.sigma(beta,sigma=2,theta,theta_tilde,Y,tau,X)
# loglik.IC.beta.general.sigma(beta=c(0,1,0),sigma=2,theta,theta_tilde,Y,tau,X)
# loglik.IC.beta.general.sigma(beta=c(0,0,1),sigma=2,theta,theta_tilde,Y,tau,X)
# loglik.IC.beta.general.sigma(beta=c(0,0,0),sigma=2,theta,theta_tilde,Y,tau,X)
# loglik.IC.beta.general.sigma(beta,sigma=1,theta,theta_tilde,Y,tau,X)
# loglik.IC.beta.general.sigma(beta=c(0,1,0),sigma=1,theta,theta_tilde,Y,tau,X)
# loglik.IC.beta.general.sigma(beta=c(0,0,1),sigma=1,theta,theta_tilde,Y,tau,X)
# loglik.IC.beta.general.sigma(beta=c(0,0,0),sigma=1,theta,theta_tilde,Y,tau,X)

#cette deuxième fonction retourne simplement la première mais ne prend en argument que les parametres qui seront optimisés
loglik.IC.beta.general2.sigma<-function(param2){
  m=length(theta)-1;m_tilde=length(theta_tilde)-1
  
  lenbeta=ncol(as.matrix(X))+1
  beta=param2[1:lenbeta]
  sigma=param2[lenbeta+1]
  
  theta=c(1,param2[(lenbeta+2):(lenbeta+2+m-1)])
  theta_tilde=c(1,param2[(lenbeta+2+m):(lenbeta+2+m+m_tilde-1)])
  return(-1*loglik.IC.beta.general.sigma(beta,sigma,theta,theta_tilde,Y,tau,X))
}

loglik.IC.beta.general2.Laplace<-function(param2){
  m=length(theta)-1;m_tilde=length(theta_tilde)-1
  
  lenbeta=ncol(as.matrix(X))+1
  beta=param2[1:lenbeta]
  sigma=param2[lenbeta+1]
  
  theta=c(1)
  theta_tilde=c(1)
  return(-1*loglik.IC.beta.general.sigma(beta,sigma,theta,theta_tilde,Y,tau,X))
}


#set.seed(1)
#beta=c(0,1,1);
#theta=theta_tilde=c(1,1,1);
#eps=rnorm(1000);
#tau=0.5;
#X1=runif(1000)
#X2=rbinom(1000,size=1,p=0.5)
#Y=1*X1 + 1*X2 +eps
#L=Y-1*runif(1000)
#R=Y+1*runif(1000)
#Y=cbind(L,R)
#X=cbind(X1,X2)
#loglik.IC.beta.general2.sigma(c(0,1,1,2,1,1,1,1))
#optim(par=c(0,1,1,2,1,1,1,1),fn=loglik.IC.beta.general2.sigma,control = list(maxit=10000))
#optim(par=c(0,1,1,0.5,1,1,1,1),fn=loglik.IC.beta.general2.sigma,control = list(maxit=10000))

#set.seed(2)
#beta=c(0,1,1);
#theta=theta_tilde=c(1,1,1);
#eps=rnorm(1000);
#tau=0.5;
#X1=runif(1000)
#X2=rbinom(1000,size=1,p=0.5)
#Y=1*X1 + 1*X2 +eps
#L=Y-1*runif(1000)
#R=Y+1*runif(1000)
#Y=cbind(L,R)
#X=cbind(X1,X2)
#loglik.IC.beta.general2.sigma(c(0,1,1,2,1,1,1,1))
#optim(par=c(0,1,1,2,1,1,1,1),fn=loglik.IC.beta.general2.sigma,control = list(maxit=10000))
#optim(par=c(0,1,1,0.5,1,1,1,1),fn=loglik.IC.beta.general2.sigma,control = list(maxit=10000))

#Attention visiblement il faudra mettre bcp de points de départ des que l'ordres des polynomes est grand
#(en tout cas si on utilise Nelder-Mead)

#set.seed(2)
#beta=c(0,1,1);
#theta=theta_tilde=c(1);
#eps=rnorm(1000);
#tau=0.5;
#X1=runif(1000)
#X2=rbinom(1000,size=1,p=0.5)
#Y=1*X1 + 1*X2 +eps
#L=Y-1*runif(1000)
#R=Y+1*runif(1000)
#Y=cbind(L,R)
#X=cbind(X1,X2)
#loglik.IC.beta.general2.sigma(c(0,1,1,2))
#optim(par=c(0,1,1,2),fn=loglik.IC.beta.general2.sigma,control = list(maxit=10000))
#optim(par=c(0,1,1,0.5),fn=loglik.IC.beta.general2.sigma,control = list(maxit=10000))


#SIGMA=seq(from=0.01,to=20,by=0.01)
#testsigma=rep(NA,length(SIGMA))
#for (i in 1:length(SIGMA)) {
#  testsigma[i] =  loglik.IC.beta.general2.sigma(c(0,1,1,SIGMA[i])) 
#}
#plot(SIGMA,testsigma)
#Pas convexe




#On va donner une valeur à sigma (donc le fixer à une autre valeur que 1 e faire
#l'optimisation) Ca a pour objectif de répondre à la question : "est ce que
#fixer une valeur un peu mauvaise pour sigma mène à de mauvais résultats
loglik.IC.beta.general2.sigmafixea2<-function(param2){
  m=length(theta)-1;m_tilde=length(theta_tilde)-1
  
  lenbeta=ncol(as.matrix(X))+1
  beta=param2[1:lenbeta]
  sigma=2
  
  if (m != 0 & m_tilde != 0) {
    theta=c(1,param2[(lenbeta+1):(lenbeta+1+m-1)])
    theta_tilde=c(1,param2[(lenbeta+1+m):(lenbeta+1+m+m_tilde-1)])
    return(-1*loglik.IC.beta.general.sigma(beta,sigma,theta,theta_tilde,Y,tau,X))
  } else if (m == 0 & m_tilde == 0) {
    theta=1
    theta_tilde=1
    return(-1*loglik.IC.beta.general.sigma(beta,sigma,theta,theta_tilde,Y,tau,X))
  } else if (m == 0 & m_tilde != 0) {
    theta=1
    theta_tilde=c(1,param2[(lenbeta+1):(lenbeta+m_tilde+0)])
    return(-1*loglik.IC.beta.general.sigma(beta,sigma,theta,theta_tilde,Y,tau,X))
  } else if (m != 0 & m_tilde == 0) {
    theta_tilde=1
    theta=c(1,param2[(lenbeta+1):(lenbeta+m+0)])
    return(-1*loglik.IC.beta.general.sigma(beta,sigma,theta,theta_tilde,Y,tau,X))
  }  
  
}


#Dans le cas où il n'y a pas dee covariate (voir ci dessous)
loglik.IC.beta.general.sigma.nocov<-function(beta, sigma, theta,theta_tilde,Y,tau)
{
  n <- nrow(Y)
  X=cbind(rep(1,n))
  
  longY <- c(Y[,1],Y[,2])
  
  res = ELD.cdf2.cont.sigma(longY,tau,rep(X%*%beta,2),sigma,theta,theta_tilde)
  
  sum(log(res[(n+1):(2*n)]-res[1:n]))
}

loglik.IC.beta.general2.sigma.nocov<-function(param2){
  m=length(theta)-1;m_tilde=length(theta_tilde)-1
  
  lenbeta=1
  beta=param2[1:lenbeta]
  sigma=param2[lenbeta+1]
  
  if (m != 0 & m_tilde != 0) {
    theta=c(1,param2[(lenbeta+1):(lenbeta+1+m-1)])
    theta_tilde=c(1,param2[(lenbeta+1+m):(lenbeta+1+m+m_tilde-1)])
    return(-1*loglik.IC.beta.general.sigma.nocov(beta,sigma,theta,theta_tilde,Y,tau))
  } else if (m == 0 & m_tilde == 0) {
    theta=1
    theta_tilde=1
    return(-1*loglik.IC.beta.general.sigma.nocov(beta,sigma,theta,theta_tilde,Y,tau))
  } else if (m == 0 & m_tilde != 0) {
    theta=1
    theta_tilde=c(1,param2[(lenbeta+1):(lenbeta+m_tilde+0)])
    return(-1*loglik.IC.beta.general.sigma.nocov(beta,sigma,theta,theta_tilde,Y,tau))
  } else if (m != 0 & m_tilde == 0) {
    theta_tilde=1
    theta=c(1,param2[(lenbeta+1):(lenbeta+m+0)])
    return(-1*loglik.IC.beta.general.sigma.nocov(beta,sigma,theta,theta_tilde,Y,tau))
  }  
  
}





simul.WeibIC.FU <- function(n=100,p=1,insptime=seq(from=-0,to=10,by=0.05),seed,shape=1.5,scale=2){
  set.seed(seed)
  Time=rweibull(n,shape=shape,scale=scale)
  L=R=rep(NA,n)
  for (i in 1:n) {
    presence=rbinom(n=length(insptime),size=1,prob=p)  
    deltaR=(Time[i]<insptime)&presence
    deltaL=(Time[i]>insptime)&presence
    if(sum(deltaR)==0){R[i]=10} else {
      R[i]=insptime[which.max(deltaR)]}
    if(sum(deltaL)==0){L[i]=0} else {
      L[i]=insptime[which.max2(deltaL)]
    }
  } 
  return(cbind(L,R))
}






simul.NormIC.FU <- function(n=100,p=1,insptime=seq(from=-5,to=5,by=0.5),seed){
  set.seed(seed)
  Time=rnorm(n)
  L=R=rep(NA,n)
  for (i in 1:n) {
    presence=rbinom(n=length(insptime),size=1,prob=p)  
    deltaR=(Time[i]<insptime)&presence
    deltaL=(Time[i]>insptime)&presence
    if(sum(deltaR)==0){R[i]=5} else {
      R[i]=insptime[which.max(deltaR)]}
    if(sum(deltaL)==0){L[i]=-5} else {
      L[i]=insptime[which.max2(deltaL)]
    }
  } 
  return(cbind(L,R))
}








simul.WeibIC.FU2 <- function(n=100,p=1,seed,shape=1.5,scale=2){
  set.seed(seed)
  Time=rweibull(n,shape=shape,scale=scale)
  L=R=rep(NA,n)
  for (i in 1:n) {
    insptime=seq(from=-6+runif(1,min=0,max=0.5),to=10,by=0.05)
    presence=rbinom(n=length(insptime),size=1,prob=p)  
    deltaR=(Time[i]<insptime)&presence
    deltaL=(Time[i]>insptime)&presence
    if(sum(deltaR)==0){R[i]=10} else {
      R[i]=insptime[which.max(deltaR)]}
    if(sum(deltaL)==0){L[i]=-6} else {
      L[i]=insptime[which.max2(deltaL)]
    }
  } 
  return(cbind(L,R))
}

# test=simul.WeibIC.FU2(n=20,p=1,seed=1)
# hist(test[,1])
# hist(test[,2])


simul.NormIC.FU2 <- function(n=100,p=1,seed){
  set.seed(seed)
  Time=rnorm(n)
  L=R=rep(NA,n)
  for (i in 1:n) {
    incr=runif(1000,min=0,max=1)
    insptimetemp=rep(NA,1000)
    insptimetemp[1]=-8+incr[1]
    dd=1
    while(insptimetemp[dd] < 8){
      insptimetemp[dd+1]=insptimetemp[dd]+incr[dd]
      dd = dd+1
    }
    insptime=insptimetemp[1:dd]
    #insptime=seq(from=-5+runif(n=1,min=0,max=0.5),to=5,by=0.5)
    presence=rbinom(n=length(insptime),size=1,prob=p)  
    deltaR=(Time[i]<insptime)&presence
    deltaL=(Time[i]>insptime)&presence
    if(sum(deltaR)==0){R[i]=5} else {
      R[i]=insptime[which.max(deltaR)]}
    if(sum(deltaL)==0){L[i]=-5} else {
      L[i]=insptime[which.max2(deltaL)]
    }
  } 
  return(cbind(L,R))
}





simul.WeibIC.FU3 <- function(n=100,p=1,seed){
  set.seed(seed)
  Time=rweibull(n,shape=1.5,scale=2)
  L=R=rep(NA,n)
  for (i in 1:n) {
    incr=runif(200,min=0,max=1)
    insptimetemp=rep(NA,100)
    insptimetemp[1]=-2+incr[1]
    dd=1
    while(insptimetemp[dd] < 15){
      insptimetemp[dd+1]=insptimetemp[dd]+incr[dd]
      dd = dd+1
    }
    insptime=insptimetemp[1:dd]
    #insptime=seq(from=-5+runif(n=1,min=0,max=0.5),to=5,by=0.5)
    presence=rbinom(n=length(insptime),size=1,prob=p)  
    deltaR=(Time[i]<insptime)&presence
    deltaL=(Time[i]>insptime)&presence
    if(sum(deltaR)==0){R[i]=20} else {
      R[i]=insptime[which.max(deltaR)]}
    if(sum(deltaL)==0){L[i]=-3} else {
      L[i]=insptime[which.max2(deltaL)]
    }
  } 
  return(cbind(L,R))
}





simul.WeibIC4 <- function(n=100,l=1,seed,shape=1.5,scale=2){
  set.seed(seed)
  Time=rweibull(n,shape=1.5,scale=2)
  L=R=rep(NA,n)
  for (i in 1:n) {
    incr=runif(500,min=0,max=2*l)
    insptimetemp=rep(NA,500)
    insptimetemp[1]=-3+incr[1]
    dd=1
    while(insptimetemp[dd] < 15){
      insptimetemp[dd+1]=insptimetemp[dd]+incr[dd]
      dd = dd+1
      if(dd == 499){insptimetemp[500] = 15}
    }
    insptime=insptimetemp[1:dd]
    R[i]=insptime[which.max(Time[i]<insptime)]
    L[i]=insptime[which.max2(Time[i]>insptime)]
  } 
  return(cbind(L,R))
}





##Transforme des observations d'une variable aleatoire (univariee) en des donnees
##censurees par intervalles. Les parametres sont adaptés a une Weibull

tranform.IC <- function(Time,l=1,seed){
  set.seed(seed)
  n=length(Time)
  L=R=rep(NA,n)
  for (i in 1:n) {
    set.seed(n*seed+i)
    incr=runif(5000,min=0,max=2*l)
    insptimetemp=rep(NA,5000)
    insptimetemp[1]=-15+incr[1]
    dd=1
    while(insptimetemp[dd] < 50){
      insptimetemp[dd+1]=insptimetemp[dd]+incr[dd]
      dd = dd+1
      if(dd == 4999){insptimetemp[5000] = 50}
    }
    insptime=insptimetemp[1:dd]
    R[i]=insptime[which.max(Time[i]<insptime)]
    L[i]=insptime[which.max2(Time[i]>insptime)]
  } 
  return(cbind(L,R))
}



##Transforme des observations non censurées en des donnees
##censurees par intervalles. 
##
##Les parametres sont adaptés a une Weibull. 
##Le temps entre 2 inspections est simulé à partir d'une distribution uniforme de parametres lmin et lmax
##La premiere inspection est en ximin et la derniere en ximax
##On assume qu'il n'y a aucune visite ratée

transform.IC.extrval <- function(Time,lmin=0.1,lmax=1.9,seed,ximin=-30,ximax=15){
  set.seed(seed)
  n=length(Time)
  L=R=rep(NA,n)
  for (i in 1:n) {
    set.seed(n*seed+i)
    incr=runif((ximax-ximin)/lmin,min=lmin,max=lmax)
    insptimetemp=rep(NA,(ximax-ximin)/lmin)
    insptimetemp[1]=ximin+incr[1]
    dd=1
    while(insptimetemp[dd] < ximax){
      insptimetemp[dd+1]=insptimetemp[dd]+incr[dd]
      dd = dd+1
      if(dd == ((ximax-ximin)/lmin)-1 ){insptimetemp[(ximax-ximin)/lmin] = ximax}
    }
    insptime=insptimetemp[1:dd]
    R[i]=insptime[which.max(Time[i]<insptime)]
    L[i]=insptime[which.max2(Time[i]>insptime)]
    #which.max(1) prend la premiere qui satisfait la condition
    #tandis que which.max2 prend la dernière qui satisfait la condition
  } 
  return(cbind(L,R))
}

#transform.IC.extrval(Yfulluncen,seed=1)


##Tentative d adapter les fonction loglik.IC pour sigma =gamma0+gamma1


ELD.cdf2.cont.gamma<-function(y,tau,mu,sigma,theta,theta_tilde){
  m <- length(theta)-1
  m_tilde <- length(theta_tilde)-1
  n <- length(y)
  
  if(length(mu) == 1){print("plutot utiliser ELD.cdf2.sigma") 
  } else if (length(mu) != length(y)){print("erreur de dim dans mu ou dans y");return()}
  
  #pour gérer quand on a des infinis dans les données
  yfinite=y[is.finite(y)]
  yinfinite=y[!is.finite(y)]
  
  y=yfinite
  n=length(y)
  ################################### D?but partie 1 #############################
  #Cette partie contient la partie o? on calcule tous les coefficients
  #Son temps de calcul n'est pas tres optimis? mais cette partie n'est ?x?cut?e
  #qu'une fois et les matrices sont souvent assez petites). Le code est r?p?t?
  #deux fois 
  Au <- matrix(0,ncol=m+1,nrow=m+1)
  for (k in 1:(m+1)) {
    Au[k,1:k] <- choose(k-1,0:(k-1)) *(-1)^(0:(k-1))/ factorial(0:(k-1))
  }
  
  Ad <- matrix(0,ncol=m_tilde+1,nrow=m_tilde+1)
  for (k in 1:(m_tilde+1)) {
    Ad[k,1:k] <- choose(k-1,0:(k-1)) *(-1)^(0:(k-1))/ factorial(0:(k-1))
  }
  
  
  Bu <- matrix(0,ncol=2*m+1,nrow=2*m+1)
  for(j.j_prime in 0:(2*m))
  {
    Bu[j.j_prime+1,(1:(j.j_prime+1))]<-factorial(j.j_prime)/factorial(0:j.j_prime)
  }
  
  
  Bd <- matrix(0,ncol=2*m_tilde+1,nrow=2*m_tilde+1)
  for(j.j_prime in 0:(2*m_tilde))
  {
    Bd[j.j_prime+1,(1:(j.j_prime+1))]<-factorial(j.j_prime)/factorial(0:j.j_prime)
  }
  
  #C contient toutes les combinaisons du produit de A j,k avec A jprime,kprime
  Cu <- array(0, dim=c(m+1,m+1,m+1,m+1))
  
  for (k in 0:m) {
    for (k_prime in 0:m) {
      for (j in 0:k) {
        for (j_prime in 0:k_prime) {
          Cu[j+1,j_prime+1,k+1,k_prime+1]<-Au[k+1,j+1] * Au[k_prime+1,j_prime+1]
        }
      }
    }
  }
  
  
  Cd <- array(0, dim=c(m_tilde+1,m_tilde+1,m_tilde+1,m_tilde+1))
  
  for (k in 0:m_tilde) {
    for (k_prime in 0:m_tilde) {
      for (j in 0:k) {
        for (j_prime in 0:k_prime) {
          Cd[j+1,j_prime+1,k+1,k_prime+1]<-Ad[k+1,j+1] * Ad[k_prime+1,j_prime+1]
        }
      }
    }
  }
  
  
  
  #Maintenant on veut connaitre le coefficient devant le terme en j* du produit
  #de vk et vkprime (avec j* in 0,1,...,k+k'). Pour ca on prend :
  #pour une valeur de k,k' la somme tel des ??ments tels que j+j' = j*
  
  #Dk kprime jpjprime = sum  Ajk*Ajprimekprime tel que j+jprime = jpjprime
  
  #Un as.mtrix a du etre rajout? car visiblement les fonctions row et col  ne 
  #marchent pas quand m et/ou mtilde = 1
  
  Du <- array(0, dim=c(m+1,m+1,2*m+1))
  for (k in 1:(m+1)) {
    for (k_prime in 1:(m+1)) {
      for (j_star in 1:(k+k_prime-1)) {
        Du[k,k_prime,j_star]<-
          sum(Cu[,,k,k_prime][(row(as.matrix(Cu[,,k,k_prime])) +col(as.matrix(Cu[,,k,k_prime])) -1)==j_star])
      }
    }
  }
  
  Dd <- array(0, dim=c(m_tilde+1,m_tilde+1,2*m_tilde+1))
  for (k in 1:(m_tilde+1)) {
    for (k_prime in 1:(m_tilde+1)) {
      for (j_star in 1:(k+k_prime-1)) {
        Dd[k,k_prime,j_star]<-
          sum(Cd[,,k,k_prime][(row(as.matrix(Cd[,,k,k_prime])) +col(as.matrix(Cd[,,k,k_prime])) -1)==j_star])
      }
    }
  }
  
  
  ################################### D?but partie 2 #############################
  
  #C'est un peu embetant de devoir passer par des data frame :( Je fais ca pour
  #pouvoir rendre le vecteur réponse dans l'ordre correspondant aux vecteurs 
  #d'entrées)
  
  data.temp=data.frame(num=1:n,y=y,mu=mu[1:length(y)],sigma=sigma[1:length(y)])
  
  data.up=subset(data.temp, y >= mu[1:length(y)])
  
  num.up=data.up$num
  
  yu <- data.up$y
  muu <- data.up$mu
  sigmau <- data.up$sigma
  su <- (yu-muu)*tau/sigmau
  lenu<-length(su)
  resu <- rep(NA,lenu)#ce vecteur contiendra le r?sultat pour les valeurs 
  #demand?es sup?rieures ? ?.
  if(lenu == 0) {resu=numeric(0)} else
    if(lenu != 0) {
      #E contient -(s)^l * exp(-s) pour diff?rentes valeurs de l dans 0,...,2m
      Eu<-matrix(NA,2*m +1,lenu)
      for (l in 1:(2*m+1)) {
        Eu[l,] = (su)^(l-1) * exp(-su) #pas sur du moins ici, ? v?rifier.  
      }
      
      "Gj+jprime,l  contient - j+j!/l! s^l exp(-s) pour toutes les valeurs de l 
  #entre 0 et 2m et de j + j_prime entre 0 et 2m"
      Gu <- array(0,dim=c(2*m+1, 2*m+1,lenu))
      for (j.j_prime in 1:(2*m+1)) {
        for (l in 1:(j.j_prime)) { #-1 cark+kprime sont plus longs de 1
          Gu[j.j_prime,l,]<- Bu[j.j_prime,l] * Eu[l,]
        }
      }
      #Hj+jprime contient - sum l=0^j+jprime j+j!/l! s^l exp(-s)
      Hu<-array(0,dim=c(2*m+1,lenu))
      if(m!=0){
        for (i in 1:lenu) {
          Hu[,i] <- rowSums(Gu[,,i])
        }
      } else if (m==0) {  
        for (i in 1:lenu) {
          Hu[,i] <- Gu[,,i]
        }}
      #J contient 
      Ju<- array(0,dim=c(m+1, m+1, lenu))
      for (k in 1:(m+1)) {
        for (k_prime in 1:(m+1)) {
          Ju[k,k_prime,]<- (Du[k,k_prime,] %*% Hu)
        }
      }
      THETA <- (theta) %*% t(theta)
      
      #Le probleme avec J C'est que c'est difficile de faire des sommes 
      #sans utiliser de boucles for sur y. Ainsi, on transforme l'array en une
      #matrice avec m+1^2 colonnes et lenu lignes et on fait les operations 
      #necessaires sur cette matrice K
      Ku <- matrix(Ju[,,1:lenu],ncol = (m+1)^2,nrow=lenu,byrow=T)*matrix(rep(THETA),byrow=T,ncol = (m+1)^2,nrow=lenu)
      
      if(m!=0){
        resu <- 1 + (1 -tau)* -rowSums(Ku)/Norm(theta)^2 
      } else if (m==0) {
        resu <- 1 + ((1-tau)* -Ku/Norm(theta)^2) 
      }
      
    }
  
  
  ##############################################################################
  #Pareil pour les observations < mu
  
  data.down=subset(data.temp, y < mu[1:length(y)])
  
  num.down=data.down$num
  
  yd <- data.down$y
  mud <- data.down$mu
  sigmad <- data.down$sigma
  sd <- -(yd-mud)*(1-tau)/sigmad
  lend<-length(sd)
  if(lend == 0) {resd=numeric(0)} else
    if(lend != 0){
      
      resd<-rep(NA,lend)
      
      Ed<-matrix(NA,2*m_tilde +1,lend)
      for (l in 1:(2*m_tilde+1)) {
        Ed[l,] = (sd)^(l-1) * exp(-sd) #pas sur du moins ici, ? v?rifier.  
      }
      
      Gd <- array(0,dim=c(2*m_tilde+1, 2*m_tilde+1,lend))
      for (j.j_prime in 1:(2*m_tilde+1)) {
        for (l in 1:(j.j_prime)) { #-1 cark+kprime sont plus longs de 1
          Gd[j.j_prime,l,]<- Bd[j.j_prime,l] * Ed[l,]
        }
      }
      
      Hd<-array(0,dim=c(2*m_tilde+1,lend))
      if(m_tilde!=0){
        for (i in 1:lend) {
          Hd[,i] <- rowSums(Gd[,,i])
        }
      } else if (m_tilde==0) {  
        for (i in 1:lend) {
          Hd[,i] <- Gd[,,i]
        }}
      
      Jd<- array(0,dim=c(m_tilde+1, m_tilde+1, lend))
      for (k in 1:(m_tilde+1)) {
        for (k_prime in 1:(m_tilde+1)) {
          Jd[k,k_prime,]<- (Dd[k,k_prime,] %*% Hd)
        }
      }
      THETA_tilde <- (theta_tilde) %*% t(theta_tilde)
      
      Kd <- matrix(Jd[,,1:lend],ncol = (m_tilde+1)^2,nrow=lend,byrow=T)*matrix(rep(THETA_tilde),byrow=T,ncol = (m_tilde+1)^2,nrow=lend)
      
      if(m_tilde!=0){
        resd <- -tau* -rowSums(Kd)/Norm(theta_tilde)^2 
      } else if (m_tilde==0) {
        resd <- -tau* -Kd/Norm(theta_tilde)^2 
      }
      
    }
  
  data.temp2<-data.frame(num=c(num.down,num.up),res=c(resd,resu))
  c(data.temp2[order(data.temp2$num),]$res,rep(1,length(yinfinite)))#on ajoute à la fin du vecteur des 1 pour chaque foi qu'il y a un infini. Ca va fonctionner que si on a des inf uniquement à la fin du vecteur
  
}



loglik.IC.beta.general.gammav1<-function(beta, gamma, theta,theta_tilde,Y,tau,X)
{
  if(gamma[1] < gamma[2]){return(-10000000)}
  
  n <- nrow(Y)
  X=cbind(rep(1,n),X)
  
  longY <- c(Y[,1],Y[,2])
  
  res = ELD.cdf2.cont.gamma(longY,tau,rep(X%*%beta,2),rep(X%*%gamma,2),theta,theta_tilde)#put exponential here in v2
  
  indexLverybig = res[1:n] > 1 - 1e-10
  indexRverysmall = res[(n+1):(2*n)] < 1e-10
  indexna = is.na(res[1:n]) | is.na(res[(n+1):(2*n)])
  indexaretirer= as.logical(indexLverybig+indexRverysmall+indexna)
  ntohide = sum(indexaretirer,na.rm = T)
  if(ntohide>= 1){
    res=res[-c(indexaretirer,n+indexaretirer)]
    n = n - ntohide  
  }
  
  sum(log(res[(n+1):(2*n)]-res[1:n]))
}


loglik.IC.beta.general2.gammav1<-function(param2){
  m=length(theta)-1;m_tilde=length(theta_tilde)-1
  
  lenbeta=ncol(as.matrix(X))+1
  beta=param2[1:lenbeta]
  
  gamma0=param2[lenbeta+1]
  gamma1=param2[lenbeta+2]
  gamma=c(gamma0,gamma1)
  
  if (m != 0 & m_tilde != 0) {
    theta=c(1,param2[(lenbeta+2):(lenbeta+2+m-1)])
    theta_tilde=c(1,param2[(lenbeta+2+m):(lenbeta+2+m+m_tilde-1)])
    return(-1*loglik.IC.beta.general.gammav1(beta,gamma,theta,theta_tilde,Y,tau,X))
  } else if (m == 0 & m_tilde == 0) {
    theta=1
    theta_tilde=1
    return(-1*loglik.IC.beta.general.gammav1(beta,gamma,theta,theta_tilde,Y,tau,X))
  } else if (m == 0 & m_tilde != 0) {
    theta=1
    theta_tilde=c(1,param2[(lenbeta+2):(lenbeta+m_tilde+1)])
    return(-1*loglik.IC.beta.general.gammav1(beta,gamma,theta,theta_tilde,Y,tau,X))
  } else if (m != 0 & m_tilde == 0) {
    theta_tilde=1
    theta=c(1,param2[(lenbeta+2):(lenbeta+m+1)])
    return(-1*loglik.IC.beta.general.gammav1(beta,gamma,theta,theta_tilde,Y,tau,X))
  }  
  
}




loglik.IC.beta.general.gammav2<-function(beta, gamma, theta,theta_tilde,Y,tau,X)
{
  n <- nrow(Y)
  X=cbind(rep(1,n),X)
  
  longY <- c(Y[,1],Y[,2])
  
  res = ELD.cdf2.cont.gamma(longY,tau,rep(X%*%beta,2),rep(exp(X%*%gamma),2),theta,theta_tilde)#put exponential here in v2
  
  indexLverybig = res[1:n] > 1 - 1e-12
  indexRverysmall = res[(n+1):(2*n)] < 1e-12
  indexna = is.na(res[1:n]) | is.na(res[(n+1):(2*n)])
  indexaretirer=  indexLverybig | indexRverysmall | indexna #as.logical(indexLverybig+indexRverysmall+indexna)
  which(indexaretirer)
  ntohide = sum(indexaretirer,na.rm = T)
  if(ntohide>= 1){
    res=res[-c(which(indexaretirer),n+which(indexaretirer))]
    #res[1:n][-indexaretirer]
    #res[(n+1):(2*n)][-indexaretirer]
    n = n - ntohide  
  }
  
  sum(log(res[(n+1):(2*n)]-res[1:n]))
}


loglik.IC.beta.general2.gammav2<-function(param2){
  m=length(theta)-1;m_tilde=length(theta_tilde)-1
  
  lenbeta=ncol(as.matrix(X))+1
  beta=param2[1:lenbeta]
  
  gamma=param2[(lenbeta+1):(2*lenbeta)]
  #3 lignes ci-dessous fonctionnent pour les simulations
  #gamma0=param2[lenbeta+1]
  #gamma1=param2[lenbeta+2]
  #gamma=c(gamma0,gamma1)
  
  if (m != 0 & m_tilde != 0) {
    theta=c(1,param2[(lenbeta+2):(lenbeta+2+m-1)])
    theta_tilde=c(1,param2[(lenbeta+2+m):(lenbeta+2+m+m_tilde-1)])
    return(-1*loglik.IC.beta.general.gammav2(beta,gamma,theta,theta_tilde,Y,tau,X))
  } else if (m == 0 & m_tilde == 0) {
    theta=1
    theta_tilde=1
    return(-1*loglik.IC.beta.general.gammav2(beta,gamma,theta,theta_tilde,Y,tau,X))
  } else if (m == 0 & m_tilde != 0) {
    theta=1
    theta_tilde=c(1,param2[(lenbeta+2):(lenbeta+m_tilde+1)])
    return(-1*loglik.IC.beta.general.gammav2(beta,gamma,theta,theta_tilde,Y,tau,X))
  } else if (m != 0 & m_tilde == 0) {
    theta_tilde=1
    theta=c(1,param2[(lenbeta+2):(lenbeta+m+1)])
    return(-1*loglik.IC.beta.general.gammav2(beta,gamma,theta,theta_tilde,Y,tau,X))
  }  
  
}






















#Dans le cas homoskedastique, modif de la fct optimisée pour qu'on puisse faire quelque chose quand pas de covariate

loglik.IC.beta.general.sigmanocov<-function(beta, sigma, theta,theta_tilde,Y,tau)
{
  n <- nrow(Y)
  X=cbind(rep(1,n))
  
  longY <- c(Y[,1],Y[,2])
  
  res = ELD.cdf2.cont.sigma(longY,tau,rep(X%*%beta,2),sigma,theta,theta_tilde)
  
  sum(log(res[(n+1):(2*n)]-res[1:n]))
}

loglik.IC.beta.general2.sigmanocov<-function(param2){
  m=length(theta)-1;m_tilde=length(theta_tilde)-1
  
  lenbeta=1
  beta=param2[1:lenbeta]
  sigma=param2[lenbeta+1]
  
  if (m != 0 & m_tilde != 0) {
    theta=c(1,param2[(lenbeta+2):(lenbeta+2+m-1)])
    theta_tilde=c(1,param2[(lenbeta+2+m):(lenbeta+2+m+m_tilde-1)])
    return(-1*loglik.IC.beta.general.sigmanocov(beta,sigma,theta,theta_tilde,Y,tau))
  } else if (m == 0 & m_tilde == 0) {
    theta=1
    theta_tilde=1
    return(-1*loglik.IC.beta.general.sigmanocov(beta,sigma,theta,theta_tilde,Y,tau))
  } else if (m == 0 & m_tilde != 0) {
    theta=1
    theta_tilde=c(1,param2[(lenbeta+2):(lenbeta+m_tilde+1)])
    return(-1*loglik.IC.beta.general.sigmanocov(beta,sigma,theta,theta_tilde,Y,tau))
  } else if (m != 0 & m_tilde == 0) {
    theta_tilde=1
    theta=c(1,param2[(lenbeta+2):(lenbeta+m+1)])
    return(-1*loglik.IC.beta.general.sigmanocov(beta,sigma,theta,theta_tilde,Y,tau))
  }  
  
}



CVtobestmGPT <- function(CV.crit) {
  # Vérification que la matrice est carrée
  if (dim(CV.crit)[1] != dim(CV.crit)[2]) {
    stop("La matrice doit être carrée")
  }
  
  # Calcul de l'indice de la valeur minimale
  ind_min <- which(CV.crit == min(CV.crit), arr.ind = TRUE)
  
  # Renvoi de l'indice sous forme de vecteur
  return(c(ind_min[1]-1, ind_min[2]-1))
}

CVtobestmGPT2 <- function(CV.crit,mmax) {
  # Vérification que la matrice est carrée
  if (dim(CV.crit)[1] != dim(CV.crit)[2]) {
    stop("La matrice doit être carrée")
  }
  
  # Exclusion des éléments de la première ligne et de la première colonne
  CV.crit[(1:(mmax+1)),1] <- Inf
  CV.crit[1,(1:(mmax+1))] <- Inf
  
  # Calcul de l'indice de la valeur minimale
  ind_min <- which(CV.crit == min(CV.crit), arr.ind = TRUE)
  
  # Renvoi de l'indice sous forme de vecteur
  return(c(ind_min[1]-1, ind_min[2]-1))
}


CVtobestmGPT3 <- function(CV.crit,mmax) {
  # Vérification que la matrice est carrée
  if (dim(CV.crit)[1] != dim(CV.crit)[2]) {
    stop("La matrice doit être carrée")
  }
  
  # Exclusion des éléments de la première ligne et de la première colonne
  CV.crit[(1:(mmax+1)),1] <- Inf
  CV.crit[(1:(mmax+1)),2] <- Inf
  CV.crit[1,(1:(mmax+1))] <- Inf
  CV.crit[2,(1:(mmax+1))] <- Inf
  
  # Calcul de l'indice de la valeur minimale
  ind_min <- which(CV.crit == min(CV.crit), arr.ind = TRUE)
  
  # Renvoi de l'indice sous forme de vecteur
  return(c(ind_min[1]-1, ind_min[2]-1))
}














simuler_points_contrainte <- function(n,m,m_tilde,beta) {
  # Initialiser une matrice pour stocker les points
  points <- matrix(NA, nrow=n, ncol=length(beta)+1+m+m_tilde)
  
  # Compter le nombre de points satisfaisant la contrainte
  nb_points_satisfaisants <- 0
  
  # Tant que nous n'avons pas suffisamment de points
  while (nb_points_satisfaisants < n) {
    # Générer aléatoirement des valeurs pour x1, x2 et x3
    x1 <- runif(length(beta), min=-0.5, max=0.5) + beta
    x2 <- runif(1, min=0.05,max=1)#Si sigma trop petit probleme de minimisation
    x3 <- runif(m, min=0, max=1)
    x4 <- runif(m_tilde, min=0, max=1)
    
    # Vérifier si les valeurs satisfont la contrainte
    if ( abs(double_constraint_fun_sigma(c(x1,x2,x3,x4)))[1] < 0.05) {
      # Si oui, stocker les valeurs dans la matrice de points
      nb_points_satisfaisants <- nb_points_satisfaisants + 1
      points[nb_points_satisfaisants,] <- c(x1, x2, x3,x4)
    }
  }
  
  return(points)
}


simuler_points_contrainte_beta <- function(n,m,m_tilde) {
  # Initialiser une matrice pour stocker les points
  points <- matrix(NA, nrow=n, ncol=length(beta)+1+m+m_tilde)
  
  # Compter le nombre de points satisfaisant la contrainte
  nb_points_satisfaisants <- 0
  
  # Tant que nous n'avons pas suffisamment de points
  while (nb_points_satisfaisants < n) {
    # Générer aléatoirement des valeurs pour x1, x2 et x3
    x1 <- runif(1, min=-0.5, max=0.5)
    x2 <- runif(1, min=0.5, max=1.5)
    x3 <- runif(1, min=0.05,max=1)#Si sigma trop petit probleme de minimisation
    x4 <- runif(m, min=0, max=1)
    x5 <- runif(m_tilde, min=0, max=1)
    
    # Vérifier si les valeurs satisfont la contrainte
    if ( abs(double_constraint_fun_sigma(c(x1,x2,x3,x4,x5)))[1] < 0.05) {
      # Si oui, stocker les valeurs dans la matrice de points
      nb_points_satisfaisants <- nb_points_satisfaisants + 1
      points[nb_points_satisfaisants,] <- c(x1, x2, x3,x4,x5)
    }
  }
  
  return(points)
}


simuler_points_contrainte_beta2 <- function(n,m,m_tilde) {
  # Initialiser une matrice pour stocker les points
  points <- matrix(NA, nrow=n, ncol=length(beta)+1+m+m_tilde)
  
  # Compter le nombre de points satisfaisant la contrainte
  nb_points_satisfaisants <- 0
  
  # Tant que nous n'avons pas suffisamment de points
  while (nb_points_satisfaisants < n) {
    # Générer aléatoirement des valeurs pour x1, x2 et x3
    x1 <- runif(1, min=-0.5, max=0.5)
    x2 <- runif(2, min=0.5, max=1.5)
    x3 <- runif(1, min=0.05,max=1)#Si sigma trop petit probleme de minimisation
    x4 <- runif(m, min=0, max=1)
    x5 <- runif(m_tilde, min=0, max=1)
    
    # Vérifier si les valeurs satisfont la contrainte
    if ( abs(double_constraint_fun_sigma(c(x1,x2,x3,x4,x5)))[1] < 0.05) {
      # Si oui, stocker les valeurs dans la matrice de points
      nb_points_satisfaisants <- nb_points_satisfaisants + 1
      points[nb_points_satisfaisants,] <- c(x1, x2, x3,x4,x5)
    }
  }
  
  return(points)
}


#Amélioration de simuler_points_contrainte_beta2 apres experience avce le 2eme papier
simuler_points_contrainte_beta2V2 <- function(n,m,m_tilde,sigmastart) {
  # Initialiser une matrice pour stocker les points
  points <- matrix(NA, nrow=n, ncol=length(beta)+1+m+m_tilde)
  
  # Compter le nombre de points satisfaisant la contrainte
  nb_points_satisfaisants <- 0
  
  #sera utilisé apres pour les theta qui sont généralement de plu en plus petits
  suite_geom=numeric(m)
  for (j in 1:m) {
    suite_geom[j] <- 1 * 0.75^(j-1)
  }
  suite_geom_2=numeric(m_tilde)
  for (j in 1:m_tilde) {
    suite_geom_2[j] <- 1 * 0.75^(j-1)
  }
  
  # Tant que nous n'avons pas suffisamment de points
  while (nb_points_satisfaisants < n) {
    # Générer aléatoirement des valeurs pour x1, x2 et x3
    x1 <- runif(1, min=-0.25, max=0.25)
    x2 <- runif(2, min=0.75, max=1.25)
    x3 <- sigmastart#Si sigma trop petit probleme de minimisation
    x4 <- runif(m, min=-0.75, max=1)* suite_geom
    x5 <- runif(m_tilde, min=-0.75, max=1) * suite_geom_2
    
    # Vérifier si les valeurs satisfont la contrainte
    if ( abs(double_constraint_fun_sigma.step1(c(x1,x2,x3,x4,x5),m,m_tilde))[1] < 0.05) {
      # Si oui, stocker les valeurs dans la matrice de points
      nb_points_satisfaisants <- nb_points_satisfaisants + 1
      points[nb_points_satisfaisants,] <- c(x1, x2, x3,x4,x5)
    }
  }
  
  return(points)
}




simuler_points_contrainte_Salary <- function(n,m,m_tilde) {
  # Initialiser une matrice pour stocker les points
  points <- matrix(NA, nrow=n, ncol=length(beta)+length(beta)+m+m)#m fixé a m tilde
  
  # Compter le nombre de points satisfaisant la contrainte
  nb_points_satisfaisants <- 0
  
  # Tant que nous n'avons pas suffisamment de points
  while (nb_points_satisfaisants < n) {
    # Générer aléatoirement des valeurs pour x1, x2 et x3
    x10 <- runif(1, min=-0.2, max=0.2)
    x11 <- runif(1, min=-0.05, max=0.05)
    x12 <- runif(1, min=-0.01, max=0.01)
    x13 <- runif(1, min=-0.05, max=0.05)
    x14 <- runif(1, min=-0.05, max=0.05)
    x15 <- runif(1, min=-0.03, max=0.03)
    x16 <- runif(1, min=-0.03, max=0.03)
    x1=c(x10,x11,x12,x13,x14,x15,x16)
    x2 <- runif(length(beta), min=0.05,max=1)#Si sigma trop petit probleme de minimisation
    x3 <- runif(m, min=0, max=1)
    x4 <- runif(m, min=0, max=1)
    
    # Vérifier si les valeurs satisfont la contrainte
    if ( abs(double_constraint_fun_sigma_nocov(c(x1,x2,x3,x4)))[1] < 0.05) {
      # Si oui, stocker les valeurs dans la matrice de points
      nb_points_satisfaisants <- nb_points_satisfaisants + 1
      points[nb_points_satisfaisants,] <- c(x1, x2, x3,x4)
    }
  }
  
  return(points)
}


# m=5
# m_tilde=4
# simuler_points_contrainte(5,m,m_tilde)
# 
# m=1
# m_tilde=1
# 
# simuler_points_contrainte(5,m,m_tilde)

double_constraint_fun_sigma_nocov <- function(x){
  if (m != 0 & m_tilde != 0) {
    part1 <- ((1+sum(x[3:(m+2)]))^2)/(1+sum(x[3:(m+2)]^2))
    part2 <- ((1+sum(x[(m+3):(m+2+m_tilde)]))^2)/(1+sum(x[(m+3):(m+2+m_tilde)]^2))
    common_term <- part1 - part2
    return(c(common_term,- common_term))
  } else if (m == 0 & m_tilde == 0) {
    return(c(0,0))
  } else if (m == 0 & m_tilde != 0) {
    part1 <- 1
    part2 <- ((1+sum(x[(m+3):(m+2+m_tilde)]))^2)/(1+sum(x[(m+3):(m+2+m_tilde)]^2))
    common_term <- part1 - part2
    return(c(common_term,- common_term))
  } else if (m != 0 & m_tilde == 0) {
    part1 <- ((1+sum(x[3:(m+2)]))^2)/(1+sum(x[3:(m+2)]^2))
    part2 <- 1
    common_term <- part1 - part2
    return(c(common_term,- common_term))
  }
}


double_constraint_fun_sigma <- function(x){
  if (m != 0 & m_tilde != 0) {
    part1 <- ((1+sum(x[(lenbeta+lengamma+1):(lenbeta+lengamma+m)]))^2)/(1+sum(x[(lenbeta+lengamma+1):(lenbeta+lengamma+m)]^2))
    part2 <- ((1+sum(x[(lenbeta+lengamma+m+1):(lenbeta+lengamma+m+m_tilde)]))^2)/(1+sum(x[(lenbeta+lengamma+m+1):(lenbeta+lengamma+m+m_tilde)]^2))
    common_term <- part1 - part2
    return(c(common_term,- common_term))
  } else if (m == 0 & m_tilde == 0) {
    return(c(0,0))
  } else if (m == 0 & m_tilde != 0) {
    part1 <- 1
    part2 <- ((1+sum(x[(lenbeta+lengamma+m+1):(lenbeta+lengamma+m+m_tilde)]))^2)/(1+sum(x[(lenbeta+lengamma+m+1):(lenbeta+lengamma+m+m_tilde)]^2))
    common_term <- part1 - part2
    return(c(common_term,- common_term))
  } else if (m != 0 & m_tilde == 0) {
    part1 <- ((1+sum(x[(lenbeta+lengamma+1):(lenbeta+lengamma+m)]))^2)/(1+sum(x[(lenbeta+lengamma+1):(lenbeta+lengamma+m)]^2))
    part2 <- 1
    common_term <- part1 - part2
    return(c(common_term,- common_term))
  }
}




ELD.cdf4 <- function(x){
  ELD.cdf2.cont.sigma(y=x,tau=tau,mu=beta0,sigma=gamma0,theta=theta0,theta_tilde=theta_tilde0)
}


# ELD.cdf4(beta0)
# 
# 
# 
# qELD4<-function(q){
#   ELD.inversecdf<-inverse(ELD.cdf4)
#   return(ELD.inversecdf(q))
# }
# 
# qELD4(0,tau=0.5,mu=beta0,sigma=sigma0,theta=theta0, theta_tilde = theta_tilde0)
# 
# qELD5<-function(q,tau,mu,theta,theta_tilde){
#   ELD.inversecdf<-myinverse(ELD.cdf2.cont.sigma)
#   return(ELD.inversecdf(q))
# }


ELD.cdf6 <- function(x){
  ELD.cdf2.sigma(y=x,tau=tau,mu=beta0,sigma=gamma0,theta=theta0,theta_tilde=theta_tilde0)
}


#ELD.cdf6(beta0)



qELD6<-function(q){
  ELD.inversecdf<-inverse(ELD.cdf6)
  return(ELD.inversecdf(q))
}

# qELD6(0.5)
# qELD6(0.75)








double_constraint_fun_sigma.step1 <- function(x,m,m_tilde){
  if (m != 0 & m_tilde != 0) {
    part1 <- ((1+sum(x[(lenbeta+lengamma+1):(lenbeta+lengamma+m)]))^2)/(1+sum(x[(lenbeta+lengamma+1):(lenbeta+lengamma+m)]^2))
    part2 <- ((1+sum(x[(lenbeta+lengamma+m+1):(lenbeta+lengamma+m+m_tilde)]))^2)/(1+sum(x[(lenbeta+lengamma+m+1):(lenbeta+lengamma+m+m_tilde)]^2))
    common_term <- part1 - part2
    return(c(common_term,- common_term))
  } else if (m == 0 & m_tilde == 0) {
    return(c(0,0))
  } else if (m == 0 & m_tilde != 0) {
    part1 <- 1
    part2 <- ((1+sum(x[(lenbeta+lengamma+m+1):(lenbeta+lengamma+m+m_tilde)]))^2)/(1+sum(x[(lenbeta+lengamma+m+1):(lenbeta+lengamma+m+m_tilde)]^2))
    common_term <- part1 - part2
    return(c(common_term,- common_term))
  } else if (m != 0 & m_tilde == 0) {
    part1 <- ((1+sum(x[(lenbeta+lengamma+1):(lenbeta+lengamma+m)]))^2)/(1+sum(x[(lenbeta+lengamma+1):(lenbeta+lengamma+m)]^2))
    part2 <- 1
    common_term <- part1 - part2
    return(c(common_term,- common_term))
  }
}
