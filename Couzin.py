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
            self.vini=uniform(-2,2)
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
        self.fichiercsv = csv.reader(self.fichier, delimiter=',')
        for ligne in self.fichiercsv:
            self.x.append(float(ligne[0]))
            self.y.append(float(ligne[1]))
            self.vx.append(float(ligne[2]))
            self.vy.append(float(ligne[3]))
        self.fichier.close()
##########################################################################################
time1=time.clock()

#Topologie de l'arène
L=10 #taille de la demi largeur du domaine

#Paramètres modèle
sigma=1.5 #ecart-type de la force stochastique
Kp=0.1 #constante de positionnement
Kal=2.7 #constante d'alignement
k_f=30.0 #coefficient de frottement
coef_lead=10.0 #doit être mis à 1 pour les simu sans leaders ou pour la condition XP4
w=1 #sensibilité à la préférence

#Capacités et préférences des individus
d0=0.4 #centre de la zone de positionnement
a=0.1
r_perc=2.0 #rayon de perception pour les individus
normv_lim=5.0
tetapref=0

#Fonctionnalités
mpl=0 # 1 si création d'une série de figure matplotlib 0 sinon
Progression='Y'
uni=0#compteur matplotlib
diz=0
cen=0
pro_av=0 #initialisation du parametre de test pour l'affichage de la progression de manière plus propre

#Temporalité
h=0.01 #pas de temps
stop=0


simu=0
N_simu=input("Combien d'expériences voulez-vous simuler ? ")
N_simu=int(N_simu)
temps=open("temps.csv","wt")
N= input("Combien d'individus compte la population ?  ")
N=int(N)
N_leader=input("Combien de leaders doit-il y avoir ? ")
N_leader=int(N_leader)
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


while simu<N_simu:
    print(simu)
    
    #génération de conditions initiales
    Cini=ConditionsInitiales()
    i=0.0
    Quest='N'
    Cini.question=Quest
    if Cini.question=='Y':
        Cini.ecrire(N)
    else:
        Cini.lire()

    x=Cini.x
    y=Cini.y
    vx=Cini.vx
    vy=Cini.vy

    pop=len(x) #nombre d'individus dans la population
    
    N_tirer=10
    tirage=[i for i in range(pop)]
    dirpref=[]
    tirer=1
    borne=pop-1
    while tirer<=N_tirer:
        elu=randint(0,borne)
        borne-=1
        dirpref.append(tirage[elu])
        tirage.remove(tirage[elu])
        tirer+=1
    print(dirpref)
    
    
    
    #Main script

    out=open("mouv.xyz","wt")
    out.write(str(pop)+'\n')


    #Boucle temporelle
    while stop==0:
        out.write(str(i)+'\n')

        #Affichage de la progression du programme
        pro=int(i)
        if  pro!=pro_av and Progression=='Y':
            print ('Progression : ' +str(int(i)))
            pro_av=pro

        N_vu=[] #nombre d'individus vus par chacun
        N_voit=[]
        for k in range(pop):
            N_vu.append(0)
            N_voit.append([])
        
        #creation des fantomes
        fantx=[]
        fanty=[]
        fantvx=[]
        fantvy=[]
        #une série de 9 comprend tous les fantomes d'un individus (plus précisément sa position et les 8 fantomes)
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

            for l in range(9): #on met dans les listes de vitesse 9x la même comme elle ne change pas pour les fantomes
                fantvx.append(vx[k])
                fantvy.append(vy[k])

        Ftot=[] #tableau des forces totales pour chaque individu 
        r_l =[] #vecteur des normes des positions à chaque instant
        

        for l in range(pop):
            Ftot.append(0)
            r_l.append(0)
            

        #Boucle sur chaque individus
        for l in range(pop):
            
            ftot=[0,0] #force totale pour un individu

            #force stochastique
            Fs=sigma*np.sqrt(12)*uniform(0,0.5)
            tetasto=uniform(-np.pi,np.pi)

            #force de frottement
            Ff=[0,0]
            if np.sqrt((vx[l])**2+(vy[l])**2)>normv_lim-1:
                Ff[0]=-k_f*vx[l]
                Ff[1]=-k_f*vy[l]

            #Boucle sur tous les autres individus
            for j in range(9*pop):
                if j/9<l or j/9>=l+1: #on regarde pas les interactions avec soi même
                    K_lead=1 #valeur de base
                    
                    d=float(str(np.sqrt((fantx[j]-x[l])**2+(fanty[j]-y[l])**2)))
                    
                    if (j/9>=numero_lead1 and j/9<numero_lead1+1 and l!=numero_lead2) or (j/9>=numero_lead2 and j/9<numero_lead2+1 and l!=numero_lead1) : #cas où on a l un individu et j un leader
                        K_lead=coef_lead
                    #Condition suivante permet si deux leaders qu'ils n'intergissent pas entre eux
                    elif (j/9>=numero_lead1 and j/9<numero_lead1+1 and l==numero_lead2) or (j/9>=numero_lead2 and j/9<numero_lead2+1 and l==numero_lead1): #cas où l un leader et j l'autre leader
                        K_lead=0
                    

                    if d<r_perc or (int(j/9)==numero_lead1 and d<r_perc) or (int(j/9)==numero_lead2 and d<r_perc): #mettre 2* r_perc pour se mettre en condition XP4 #sinon représente condition pour être dans le rayon de perception
                        
                        N_vu[l]=N_vu[l]+1 #Ajoute l'indiv j à la liste des membres perçus
                        N_voit[l].append(int(j/9))
                        
                        r_l[l]=np.sqrt((x[l])**2+(y[l])**2) #position dans l'espace de l'individu l
                        normv_l=np.sqrt((vx[l])**2+(vy[l])**2) #norme de la vitesse de l'individu l
                        normv_j=np.sqrt((fantvx[j])**2+(fantvy[j])**2) #norme de la vitesse de l'individu j 


                        if normv_l!=0 and normv_j!=0:
                            cosavit=(vx[l]*fantvx[j]+vy[l]*fantvy[j])/(normv_l*normv_j) #cosinus de l'angle entre les directions des deux individus
                            sinavit=(fantvy[j]*vx[l]-fantvx[j]*vy[l])/(normv_l*normv_j) #sinus de l'angle entre les directions des deux individus
                            cosperc=(vx[l]*(fantx[j]-x[l])+vy[l]*(fanty[j]-y[l]))/(d*normv_l)
                            sinperc=(vy[l]*(fantx[j]-x[l])-vx[l]*(fanty[j]-y[l]))/(d*normv_l)
                            cosvl=vx[l]/normv_l
                            sinvl=vy[l]/normv_l
                            cosvj=fantvx[j]/normv_j
                            sinvj=fantvy[j]/normv_j
                            cospref=np.cos(tetapref)
                            sinpref=np.sin(tetapref)
                            cosalpref=(vx[l]*cospref+vy[l]*sinpref)/(normv_l)
                            sinalpref=(sinpref*vx[l]-cospref*vy[l])/(normv_l)
                        else:
                            cosavit=0
                            sinavit=0
                            cosperc=0
                            sinperc=0
                            cosvl=0
                            sinvl=0
                            cosvj=0
                            sinvj=0
                            cosalpref=0
                            sinalpref=0

                        #Force de positionnement

                        #Fp=K_lead*Kp*(d-d0)#*(1+cosperc) #champ de vision frontal et latéral en diminuant
                        #Fp=Kp*(d-d0)*(sinperc) #champ de vision latéral #calovi gautrais
                        #Fp=-4*Kp*(-2*(d0/2)**2/d**3+(d0/2)/d**2)#*sinperc #potentiel de lennard Jones
                        #Fp=-Kp*(12*d**11-(d0*2**13)/d**2)
                        Fp=2*Kp*a*m.exp(-a*(d-d0))*(1-m.exp(-a*(d-d0)))

                        #Force d'alignement    
                        if sinavit!=0 :
                            Fa=K_lead*Kal*abs(sinavit) #permission de l'anti parallèle #Calovi/Gautrais
                            #Fa=Kal*(1-cosavit)*(sinavit)/abs(sinavit) 
                        else:
                            Fa=0


                        #composantes de la force individuelle totale
                        #projeté des forces d'intéractions      
                        ftot[0]+=Fp*((fantx[j]-x[l])/d)+Fa*cosvj
                        ftot[1]+=Fp*((fanty[j]-y[l])/d)+Fa*sinvj
                        
                        #Calcul de la force de préférence
                        if l in dirpref:
                            Fpref=w*(1-cosalpref)

            #on rajoute la force stochastique et la force de frottement              
                    
            ftot[0]=ftot[0]+Fs*np.cos(tetasto)+Ff[0]+Fpref*cospref
            ftot[1]=ftot[1]+Fs*np.sin(tetasto)+Ff[1]+Fpref*sinpref

            #Matrice des forces totales
            Ftot[l]=ftot
         

        #initialisation des next coord
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

            if newx[q]>L:
                newx[q]-=L*2
            elif newx[q]<-L:
                newx[q]+=L*2

            if newy[q]>L:
                newy[q]-=L*2
            elif newy[q]<-L:
                newy[q]+=L*2

            if q==numero_lead1 or q==numero_lead2:
                out.write('N'+' '+str(newx[q])+' '+str(newy[q])+' 0'+'\n')    #écriture des nouvelles coordonnées pour le(s) leader(s)
            else:
                out.write('C'+' '+str(newx[q])+' '+str(newy[q])+' 0'+'\n')    #écriture des nouvelles coordonnées pour les autres 



            #Actualisation des coord
            x[q]=newx[q]
            y[q]=newy[q]
            vx[q]=newvx[q]
            vy[q]=newvy[q]
            



            
        out.write('\n')

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
        #condition d'arrêt de la boucle
        obj=0
        
        
        if i>20: 
            stop=1
                      

        i=i+h

    out.close()
    
    stop=0
    
    temps.write(str(round(i,2))+' ')
    print("Le temps mis pour se regrouper est de "+str(round(i,2))+' s.')
    simu+=1
    
temps.close()    
time2=time.clock()
print("\nTemps d'execution: " +str(round(time2-time1,2))+' s')

print(N,N_leader)