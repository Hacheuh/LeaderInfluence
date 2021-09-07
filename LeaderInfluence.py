#coding=UTF8
import csv
import numpy as np
import math as m
import matplotlib.pyplot as plt
import os
import time
import pylab
from random import uniform
from random import randint
##########################################################################################
class ConditionsInitiales:
    """Cette classe permet de gerer les conditions initiales"""
    
    def __init__(self):
        """Par défaut on ne veut pas générer de nouvelles conditions initiales"""
        self.question='N'
        self.x=[]
        self.y=[]
        self.vx=[]
        self.vy=[]
        self.vini=0
        self.tetavini=[]
        self.fichier=0
        self.fichiercsv=0
            
    def ecrire(self,population):
        """Methode permettant de générer de nouvelles conditions initiales dans le fichier correspondant"""
        self.fichier = open ('pop.csv','wt')
        for i in range(population):
            self.x.append(uniform(-1,1))
            self.y.append(uniform(-1,1))
            self.vini=uniform(0,2)
            self.tetavini.append(uniform(-np.pi,np.pi))
            self.vx.append(self.vini*np.cos(self.tetavini[i]))
            self.vy.append(self.vini*np.sin(self.tetavini[i]))

            if i==N-1:
                self.fichier.write(str(self.x[i])+','+str(self.y[i])+','+str(self.vx[i])+','+str(self.vy[i])+','+str(self.tetavini[i]))
            else:
                self.fichier.write(str(self.x[i])+','+str(self.y[i])+','+str(self.vx[i])+','+str(self.vy[i])+','+str(self.tetavini[i])+'\n')
    
        self.fichier.close()
    
    def lire(self):
        """Methode permettant de lire le fichier de conditions initiales"""
        self.fichier = open ('pop.csv','rt')
        self.fichiercsv = csv.reader(self.fichier, delimiter=' ')
        for ligne in self.fichiercsv:
            self.x.append(float(ligne[0]))
            self.y.append(float(ligne[1]))
            self.vx.append(float(ligne[2]))
            self.vy.append(float(ligne[3]))
        self.fichier.close()
##########################################################################################
time1=time.clock()

#Topologie de l'arène
L=10.0      #taille de la demi largeur du domaine

#Paramètres modèle
sigma=5.0   #ecart-type de la force stochastique
Kp=1.0      #constante de positionnement
Kal=0.3     #constante d'alignement
k_f=30.0    #coefficient de frottement
w=10.0      #sensibilité à la préférence

#Capacités et préférences des individus
d0=0.6      #centre de la zone de positionnement #valable qu'en intéraction de paire
a=1.0       #largeur du potentiel de morse
r_perc=2.0  #rayon de perception pour les individus
normv_lim=5.0   #limite de vitesse
tetapref1=np.pi/4   #choix de la majorité
tetapref2=0.0


#Fonctionnalités
mpl=0           # 1 si création d'une série de figure matplotlib, 0 sinon
Progression='Y' #Affichage de la progression dans l'invite de commande Y/N
#compteur matplotlib
uni=0
diz=0
cen=0
pro_av=0        #initialisation du parametre de test pour l'affichage de la progression de manière plus propre


#Temporalité
h=0.01          #pas de temps
         
N_simu=input("Combien d'expériences voulez-vous simuler ? ")
N_simu=int(N_simu)

#Sur le même format possibilité de varier le coef leader ou la taille de groupe, ici choix du coef leader
N=30            #Taille de groupe
coef=[1.0,3.0,10.0,30.0,100.0]     #Liste des coefficients leader à parcourir
N_leader=input("Combien de leaders doit-il y avoir ? ")
N_leader=int(N_leader)

#Boucle en paramètre significatif
for coef_lead in coef:
    simu=0      #initialisation du nombre de simulations
    if N_leader==0:
        coef_lead=1
        numero_lead1,numero_lead2=113,114

    elif N_leader==1:
        numero_lead1,numero_lead2=3,114

    elif N_leader==2:
        numero_lead1,numero_lead2=3,4

    else:
        N_simu=0
        print("Le nombre de leader n'est pas valide, veuillez sélectionner parmi 0,1 ou 2.")

    #fichier de test
    DPIcsv=open("distpremindiv.csv","wt")       #Distance au premier individu
    tetacsv=open("direction.csv","wt")          #Directions prises à tous les temps par les individus
    DIIcsv=open("distinterindiv.csv","wt")      #Distance interindividuelles
    tetafocsv=open("force.csv","wt")            #Orientation des forces à chaque instant
    tpocsv=open("transparam.csv","wt")          #Translational parameter order
    rcpocsv=open("rotcoheparam.csv","wt")       #Rotational cohesiveness parameter order
    vitmoy=open("vecteurvitessemoyen"+str(coef_lead)+"_"+str(N_leader)+".csv","wt") #Orientation du vecteur vitesse moyen à chaque instant
    effi=open("efficacité"+str(coef_lead)+"_"+str(N_leader)+".csv","wt")            #Orientation du vecteur vitesse individuel à la fin de la simulation
    temps=open("temps"+str(coef_lead)+"_"+str(N_leader)+".csv","wt")                #temps mis pour réaliser la tâche
    final=open("Vfinalemoy"+str(coef_lead)+"_"+str(N_leader)+".csv","wt")           #image de l'orientation finale moyenne du groupe
    
    #boucle en simulations
    while simu<N_simu:
        print(simu)
        
        #génération de conditions initiales
        Cini=ConditionsInitiales()
        i=0.0       #initialisation du temps
        pro_av=0    #initialisation du parametre de test pour l'affichage de la progression de manière plus propre
        Quest='Y'   #Génération de nouvelles conditions initiales Y/N
        Cini.question=Quest
        if Cini.question=='Y':
            Cini.ecrire(N)
        else:
            Cini.lire()
        
        #Ecriture des conditions initiales dans les variables
        x=Cini.x
        y=Cini.y
        vx=Cini.vx
        vy=Cini.vy

        pop=len(x) #nombre d'individus dans la population
        DPI=[]
        DII=[]

        #Boucle de tirage des directions
        N_tirer=N
        tirage=[tir for tir in range(pop)]
        tirage.remove(3)        #le leader est placé directement en désaccord de la majorité
        
        dirpref1=[]
        dirpref1.append(3)
        dirpref2=[]
        tirer=2
        borne=pop-2
        while tirer<=N_tirer:
            elu=randint(0,borne)
            borne-=1
            choix=randint(1,100)
            if choix<=20:
                dirpref1.append(tirage[elu])
                tirage.remove(tirage[elu])
            elif choix>20:
                dirpref2.append(tirage[elu])
                tirage.remove(tirage[elu])
            tirer+=1
        print(dirpref1)
        print(dirpref2)


        #Main script

        out=open("mouv"+str(coef_lead)+"_"+str(N_leader)+".xyz","wt") #fichier de visualisation
        out.write(str(pop)+'\n')


        #Boucle temporelle
        stop=0      #initialisation de la condition d'arrêt du programme 
        obj=0       #variable de comptage pour l'arrêt
        
        #Boucle de simu individuelle
        while stop==0:
            
            out.write(str(i)+'\n')

            #Affichage de la progression du programme
            pro=int(i)
            if  pro!=pro_av and Progression=='Y':
                print ('Progression : ' +str(int(i)))
                pro_av=pro

            N_vu=[]         #nombre d'individus vus par chacun
            N_voit=[]       #liste des listes des individus vus par chaque individus
            for k in range(pop):
                N_vu.append(0)
                N_voit.append([])             
           
            

            #creation des fantomes
            fantx=[]
            fanty=[]
            fantvx=[]
            fantvy=[]
            for k in range(pop):
                fantx.append(x[k])
                fantx.append(x[k]+2*L)
                fantx.append(x[k])
                fantx.append(x[k]-2*L)
                fantx.append(x[k])
                fantx.append(x[k]+2*L)
                fantx.append(x[k]-2*L)
                fantx.append(x[k]-2*L)
                fantx.append(x[k]+2*L)

                fanty.append(y[k])
                fanty.append(y[k])
                fanty.append(y[k]+2*L)
                fanty.append(y[k])
                fanty.append(y[k]-2*L)
                fanty.append(y[k]+2*L)
                fanty.append(y[k]+2*L)
                fanty.append(y[k]-2*L)
                fanty.append(y[k]-2*L)

                for l in range(9):
                    fantvx.append(vx[k])
                    fantvy.append(vy[k])

            Ftot=[] #tableau des forces totales pour chaque individu 
            r_l =[] #vecteur des normes des positions à chaque instant


            for l in range(pop):
                Ftot.append(0)
                r_l.append(0)


            RCPO=0      #rotational cohesiveness parameter order
            TPOs=[0,0]  #translational parameter order
            TPO=0       #translational parameter order
            effilist=[] #initialisation de la liste des tetas de chaque individus
            
            #Boucle sur chaque individus
            for l in range(pop):
                
                filtreeffilist=1 #le filtre permet de remplir une fois le teta de l'individu l et pas 29 fois pour chaque (boucle j)
                
                barycentre=[0,0]    #position du barycentre du groupe
                for l_1 in range(pop):
                    barycentre[0]+=(x[l_1]/pop)
                    barycentre[1]+=(y[l_1]/pop)
                TPOs[0]+=(vx[l]/np.sqrt(vx[l]**2+vy[l]**2))
                TPOs[1]+=(vy[l]/np.sqrt(vx[l]**2+vy[l]**2))

                ny=1/(np.sqrt(1+((barycentre[1]-y[l])/(barycentre[0]-x[l]))**2))
                nx=np.sqrt(1-ny**2)
                RCPO+=((vx[l]/np.sqrt(vx[l]**2+vy[l]**2)*nx+vy[l]/np.sqrt(vx[l]**2+vy[l]**2)*ny)/pop)

                ftot=[0,0] #force totale pour un individu

                #force stochastique
                Fs=sigma*np.sqrt(12)*uniform(0,0.5)
                tetasto=uniform(-np.pi,np.pi)

                #force de frottement
                Ff=[0,0]
                if np.sqrt((vx[l])**2+(vy[l])**2)>normv_lim-1:
                    Ff[0]=-k_f*vx[l]
                    Ff[1]=-k_f*vy[l]

                #Boucle sur tous les autres individus possiblement en interaction
                DP=[] 
                DI=[]
                for j in range(9*pop):
                    if int(j/9)!=l:
                        K_leadp,K_leada=1,1 

                        d=np.sqrt((fantx[j]-x[l])**2+(fanty[j]-y[l])**2)

                        if (int(j/9)==numero_lead1 and l!=numero_lead2) or (int(j/9)==numero_lead2 and l!=numero_lead1) :
                            K_leadp=coef_lead/2
                            K_leada=coef_lead
                        #Condition suivante permet si deux leaders qu'ils n'intergissent pas entre eux
                        elif (int(j/9)==numero_lead1 and l==numero_lead2) or (int(j/9)==numero_lead2 and l==numero_lead1):
                            K_leadp,K_leada=0,0


                        if d<r_perc:    #condition pour uniquement prendre les individus dans le champs de vision
                            
                            N_vu[l]=N_vu[l]+1           #Ajoute l'indiv j à la liste des membres perçus
                            N_voit[l].append(int(j/9))
                            
                            DP.append(d)
                            if int(j/9)>=l+1:
                                DI.append(d)

                            r_l[l]=np.sqrt((x[l])**2+(y[l])**2) #position dans l'espace de l'individu l
                            normv_l=np.sqrt((vx[l])**2+(vy[l])**2) #norme de la vitesse de l'individu l
                            normv_j=np.sqrt((fantvx[j])**2+(fantvy[j])**2) #norme de la vitesse de l'individu j 


                            if normv_l!=0 and normv_j!=0:   #les individus doivent être en mouvement
                                cosavit=(vx[l]*fantvx[j]+vy[l]*fantvy[j])/(normv_l*normv_j) #cosinus de l'angle entre les directions des deux individus
                                sinavit=(fantvy[j]*vx[l]-fantvx[j]*vy[l])/(normv_l*normv_j) #sinus de l'angle entre les directions des deux individus
                                cosperc=(vx[l]*(fantx[j]-x[l])+vy[l]*(fanty[j]-y[l]))/(d*normv_l)   #cosinus de l'angle entre la direction dans laquelle l va et la position où j est
                                sinperc=(vy[l]*(fantx[j]-x[l])-vx[l]*(fanty[j]-y[l]))/(d*normv_l)   #sinus de l'angle entre la direction dans laquelle l va et la position où j est
                                cosvl=vx[l]/normv_l     #cos du vecteur vitesse de l
                                sinvl=vy[l]/normv_l     #sin du vecteur vitesse de l
                                cosvj=fantvx[j]/normv_j #cos du vecteur vitesse de j
                                sinvj=fantvy[j]/normv_j #sin du vecteur vitesse de j
                                cospref1=np.cos(tetapref1)  #cos du vecteur direction dans le cas où la préférence est l'angle 1
                                sinpref1=np.sin(tetapref1)  #sin du vecteur direction dans le cas où la préférence est l'angle 1
                                cospref2=np.cos(tetapref2)  #cos du vecteur direction dans le cas où la préférence est l'angle 2
                                sinpref2=np.sin(tetapref2)  #sin du vecteur direction dans le cas où la préférence est l'angle 2
                            else:
                                cosavit=0
                                sinavit=0
                                cosperc=0
                                sinperc=0
                                cosvl=0
                                sinvl=0
                                cosvj=0
                                sinvj=0
                                cospref1=np.cos(tetapref1)
                                sinpref1=np.sin(tetapref1)
                                cospref2=np.cos(tetapref2)
                                sinpref2=np.sin(tetapref2)

                            teta=m.acos(cosvl)
                            if sinvl<0:
                                    teta=-teta
                            if i>100*h:
                                tetacsv.write(str(teta)+'\n')
                                if filtreeffilist==1:
                                    effilist.append(teta)
                                    filtreeffilist=0
                                    
                            #Force de positionnement
                            Fp=K_leadp*2*Kp*a*m.exp(-a*(d-d0))*(1-m.exp(-a*(d-d0))) #potentiel de Morse

                            #Force d'alignement    
                            if sinavit!=0 :
                                Fa=K_leada*Kal
                            else:
                                Fa=0

                            #projeté des forces d'intéractions + ajout de la force cible      
                            ftot[0]=ftot[0]+(Fp*((fantx[j]-x[l])/d)+Fa*cosvj)
                            ftot[1]=ftot[1]+(Fp*((fanty[j]-y[l])/d)+Fa*sinvj)

                            #Calcul de la force de préférence
                            if l in dirpref1:
                                Fpref=[w*cospref1,w*sinpref1]
                            elif l in dirpref2:
                                Fpref=[w*cospref2,w*sinpref2]
                            else:
                                Fpref=[0,0]


                    #on rajoute la force stochastique et de frottement               
                ftot[0]=ftot[0]+Fs*np.cos(tetasto)+Ff[0]+Fpref[0]
                ftot[1]=ftot[1]+Fs*np.sin(tetasto)+Ff[1]+Fpref[1]


                #Matrice des forces totales
                Ftot[l]=ftot
                
            #calcul des test de force, de DII, de TPO et de RCPO
                if i>100*h:
                    tetafo=m.acos(ftot[0]/(np.sqrt(ftot[0]**2+ftot[1]**2)))
                    if (ftot[1]/(np.sqrt(ftot[0]**2+ftot[1]**2)))<0:
                        tetafo=-tetafo
                    tetafocsv.write(str(tetafo)+'\n')
                    try:
                        DPI.append(min(DP))
                    except ValueError:
                        1==1
                    for dist in DI:
                        DII.append(dist)
                        
            TPO=np.sqrt(np.sqrt(TPOs[0]**2+TPOs[1]**2))/pop
            tpocsv.write(str(TPO)+'\n')
            rcpocsv.write(str(RCPO)+'\n')
            
            #initialisation des nouvelles coordonnées
            newvx=[]
            newvy=[]                    
            newx=[]                    
            newy=[]    
            for p in range(pop):
                newvx.append(0.0)
                newvy.append(0.0)
                newx.append(0.0)
                newy.append(0.0)  

            #Calcul des next coord
            for q in range(pop):
                newvx[q]=vx[q]+Ftot[q][0]*h
                newvy[q]=vy[q]+Ftot[q][1]*h
                newx[q]=x[q]+vx[q]*h+Ftot[q][0]*h**2
                newy[q]=y[q]+vy[q]*h+Ftot[q][1]*h**2

                #Passage à travers les bords périodiques
                if newx[q]>L:
                    newx[q]-=L*2
                elif newx[q]<-L:
                    newx[q]+=L*2

                if newy[q]>L:
                    newy[q]-=L*2
                elif newy[q]<-L:
                    newy[q]+=L*2
                    
                #écriture des nouvelles coordonnées dans le fichier visu    
                if q==numero_lead1 or q==numero_lead2:
                    out.write('N'+' '+str(newx[q])+' '+str(newy[q])+' 0'+'\n')
                else:
                    out.write('C'+' '+str(newx[q])+' '+str(newy[q])+' 0'+'\n')    

                #Actualisation des coord
                x[q]=newx[q]
                y[q]=newy[q]
                vx[q]=newvx[q]
                vy[q]=newvy[q]

            #condition d'arrêt    
            Vmoy=[0,0]
            for j in range(N):
                Vmoy[0]+=(vx[j]/np.sqrt((vx[j])**2+(vy[j])**2))/pop
                Vmoy[1]+=(vy[j]/np.sqrt((vx[j])**2+(vy[j])**2))/pop 
            if np.sqrt((Vmoy[0])**2+(Vmoy[1])**2)>0.99939 or i>=150: #si l'orientation moyenne est de norme>0.9 i.e. le faisceau angulaire est au plus de 15° alors on stop la simu 
                obj+=1
                

                if obj>=100 or i>=150:#on attends d'etre sûr que les fluctuations restent dans le faisceau pendant au moins 100 pas de temps
                    stop=1
                    for k in effilist:
                        effi.write(str(k)+'\n')
            else:
                obj=0
            vitmoy.write(str(Vmoy[0])+','+str(Vmoy[1])+'\n')
            i+=h
                        

               
            out.write('\n')
        out.close()
    ###################################################################################
            #module matplotlib
            if mpl==1:
                if round(i,2)*100%4==0:
                    xmin, ymin, xmax, ymax = -10, -10, 10, 10
                    plt.xlim(xmin, xmax) # Limites des axes du repère
                    plt.ylim(ymin, ymax)
                    plt.quiver(x,y,vx,vy,width=0.005, pivot='mid')
                    #plt.show()
                    pylab.savefig('C:/Users/Hadrien/Desktop/Leader/Video/Creation/test'+str(cen)+str(diz)+str(uni)+'.png', bbox_inches='tight')
                    plt.close()
                    uni+=1
                    if uni==10:
                        uni=0
                        diz+=1
                        if diz==10:
                            diz=0
                            cen+=1
    #####################################################################################           
        
        #Ecriture dans les fichiers de test
        for dist in DPI:
            DPIcsv.write(str(dist)+'\n')
        for dist in DII:
            DIIcsv.write(str(dist)+'\n')
        effi.write('\n')
        temps.write(str(round(i,2))+' ')
        final.write(str(Vmoy[0])+','+str(Vmoy[1])+'\n')
        print("Le temps mis pour se regrouper est de "+str(round(i,2))+' s.')
        simu+=1
        
    final.close()
    vitmoy.close()    
    DPIcsv.close()
    tetacsv.close()
    DIIcsv.close()
    tetafocsv.close()
    tpocsv.close()
    rcpocsv.close()
    effi.close()
    temps.close()
    time2=time.clock()
    
    
    print("\nTemps d'execution: " +str(round(time2-time1,2))+' s')
    print(N,N_leader)