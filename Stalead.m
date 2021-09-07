% CE FICHIER CONTIENT TOUS LES TRAITEMENTS MATLAB QUI ONT ETE UTILISE POUR
% L'ANALYSE DES SIMULATIONS. NECESSITE DE COMMENTER LES PARTIES INUTILES ET
% ACTIVER EN DECOMMENTANT LES PARTIES D'INTERET. IL EST RECOMMANDE DE FAIRE
% ATTENTION LORS DE L'ACTIVATION DE PLUSIEURS MODULES EN MEME TEMPS

%%
%%%%%%%%%%% MODULE 1 : ANALYSE DU COMPORTEMENT DU MODELE SANS TACHES %%%%%
%%%Affichage sous forme d'histogramme de l'orientation de la force, de la
%%%distance inter_individuelle, de la distance au premier individu, de la
%%%direction du mouvement (orientation vitesse), et des param�tres d'ordre
%%%translationnel et rotationnel (pas utilis� jusqu'a pr�sent)

% cd 'C:\Users\*******' 

% Force=csvread('force.csv');
% DII=csvread('distinterindiv.csv');
% DPI=csvread('distpremindiv.csv');
% Direction=csvread('direction.csv');
% TPO=csvread('transparam.csv');
% RCPO=csvread('rotcoheparam.csv');

% %figure;
% %hist(DII(:,1),50)

% %figure;
% %hist(DPI(:,1),50)

% figure;
% xscal=floor(50*(max(Direction)-min(Direction))/(2*pi));
% hist(Direction(:,1),50)% utiliser rose() pour un histogramme polaire
% V=axis;
% axis([-pi pi V(3) V(4)]);

% figure;
% hist(Force(:,1),50)
% V=axis;
% axis([-pi pi V(3) V(4)]);

% %figure;
% %hist(TPO(:,1),50)

% %figure;
% %hist(RCPO(:,1),50)
%
%%% Quantification permettant de voir � quel point les forces sont
%%% uniformes en direction ou pas 
% % counts=hist(Force(:,1),50);
% % (max(counts)-min(counts))/mean(counts)

%%
%%%%%%%%%%%%%%%%%%%%% MODULE 2: TRAITEMENT DES SIMULATIONS TACHE ORIENTATION %%%%%%%%%%%%%%%%%%%%
% cd 'C:\Users\******'
% clear all

% %%% Choisir l'indicatif relatif aux tailles de groupes ou au coefficients
% %%% leaders, l'indicatif2 est un indicatif normalis� par la taille de
% %%% groupe
% indicatif={'10.0','11.5','13.0','15.4','18.0','20.5','23.7','27.4','30.0','36.5','42.1','48.6','56.0','64.9','75.0','86.6','100.0'};
% indicatif={'5','10','15','20','25','30','35','40','45'}
% indicatif2={num2str(round(10/N,2)),num2str(round(11.5/N,2)),num2str(round(13.0/N,2)),num2str(round(15.4/N,2)),num2str(round(18.0/N,2)),num2str(round(20.5/N,2)),num2str(round(23.7/N,2)),num2str(round(27.4/N,2)),num2str(round(30/N,2)),num2str(round(36.5/N,2)),num2str(round(42.1/N,2)),num2str(round(48.6/N,2)),num2str(round(56.0/N,2)),num2str(round(64.9/N,2)),num2str(round(75.0/N,2)),num2str(round(86.6/N,2)),num2str(round(100/N,2))};
% ind=[10,11.5,13.0,15.4,18.0,20.5,23.7,27.4,30,36.5,42.1,48.6,56.0,64.9,75.0,86.6,100];

% %%% Handles pour le trac� par r�currence de chaque violin relatif � un
% %%% param�tre
% h1=figure('name','Violin plot des distributions d �cart-type � param�tre fix�');
% h2=figure('name','Violin plot des distributions de temps � param�tre fix�');
% h3=figure('name','Violin plot des distributions d orientations de groupe � param�tre fix�');

% clear effi* EFFI* final* temps*
% N=30; % Taille de groupe
% i=1; % Initialisation du num�ro de fichier

% %%% Boucle de chargement des fichiers orientation
% for j = indicatif
%     j1=strcat('efficacit�',j,'_0.csv');
%     j1=j1{1};
%     eval(sprintf('effi%g=csvread(%c%s%c);',i,39,j1,39));
%     i=i+1;
% end
% 
% k=size(indicatif,2); % k est le nombre de fichier qui ont �t� charg�s

% %%% Boucle sur les orientations et les ecart-types associ�s. La partie
% %%% triplement comment�e permet le trac� d'un histogramme complet des
% %%% orientations individuelles, la partie doublement comment�e trace
% %%% directement les distributions d'ecart-types associ�es sous forme
% %%% d'histogrammes, la partie simplement comment�e g�n�re les fichier EFFI
% %%% contenant les ecart-types associ�s
% for i=1:k
% % %     figure('name',sprintf('Distribution des orientations finales individuelles (compil�es) n� %s',num2str(i)),'NumberTitle','off');
% % %     xscal=floor(50*(max(eval(sprintf('effi%g',i)))-min(eval(sprintf('effi%g',i))))/(2*pi));
% % %     hist(eval(sprintf('effi%g(:,1)',i)),50); % utiliser rose() pour un histogramme polaire
% % %     V=axis;
% % %     axis([-0.2 pi/4+0.2 V(3) V(4)]);
%     for j=1:199%(size(eval(sprintf('effi%g',i)),1)/30)
%         eval(sprintf('EFFI%g(%g)=std(effi%g((%g-1)*N+1:(%g*N)));',i,j,i,j,j));
%     end
% %     figure('name',sprintf('Distribution des ecart-type � param�tre fix� n� %s',num2str(k)),'NumberTitle','off');
% %     hist(eval(sprintf('EFFI%g',i)),10);
% %
% %     %%% Quantification des distributions
% %     Moment.effi(i,1)=mean(eval(sprintf('EFFI%g',i)));
% %     Moment.effi(i,2)=moment(eval(sprintf('EFFI%g',i)),2);
% %     Moment.effi(i,3)=moment(eval(sprintf('EFFI%g',i)),3);
% %     Moment.effi(i,4)=moment(eval(sprintf('EFFI%g',i)),4);
% end
% %
% %%% Cette partie est � pr�f�rer � a partie trac� de la pr�c�dente, n�cessite la
% %%% conservation de la g�n�ration des fichiers 'EFFI' de la partie
% %%% pr�c�dente(partie simplement comment�e)
% vio=[];
% for i =1:k
%     vio=[vio eval(sprintf('EFFI%g',i)).'];
% end
% 
% % figure(h1); % seulement si lancement de la totalit� du module
% violin(vio,'xlabel',indicatif2,'facecolor','y','facealpha',0.4);
% set(gca,'XTickMode','auto');
% set(gca,'XTickLabel',indicatif2);
% set(gca,'XTick',[1:1:34]);
% xlabel('Coefficient Leader normalis�');
% ylabel('Efficacit� (std)');
% 
% i=1; %r�initialisation du numero de fichier
% %%%Boucle de chargement des fichiers temps
% for j = indicatif
%     j1=strcat('temps',j,'_1.csv');
%     j1=j1{1};
%     eval(sprintf('temps%g=csvread(%c%s%c);',i,39,j1,39));
%     i=i+1;
% end
% 
% k=size(indicatif,2); % nombre de fichiers
% 
% %%% Trac� et quantification des fichiers temps sous forme d'histogrammes,
% %%% Pr�f�rez le trac� par violin plot plus bas 
% for i=1:k
%     %     figure('name',sprintf('Distribution des temps � param�tre fix� n� %s',num2str(k),'NumberTitle','off');
%     %     xscal=floor(30*(max(eval(sprintf('temps%g(1,:)',i)))-min(eval(sprintf('temps%g(1,:)',i))))/15);
%     %     hist(eval(sprintf('temps%g(1,:)',i)),xscal);
%     %     V=axis;
%     %     axis([0 15 0 V(4)]);
%     Moment.temps(i,1)=mean(eval(sprintf('temps%g',i)));
%     Moment.temps(i,2)=moment(eval(sprintf('temps%g',i)),2);
%     Moment.temps(i,3)=moment(eval(sprintf('temps%g',i)),3);
%     Moment.temps(i,4)=moment(eval(sprintf('temps%g',i)),4);
% end
%
% %figure(h2); % si utilisation du module complet
% %figure('name',sprintf('Violin plot des distributions de temps � param�tre fix�, groupe %g',N),'NumberTitle','off'); % si utilisation de cette partie de module uniquement
% vio=[];
% for i =1:k
%     vio=[vio eval(sprintf('temps%g',i)).'];
% end
% violin(vio,'xlabel',indicatif2,'facecolor','y','facealpha',0.4);
%
% set(gca,'XTickMode','auto');
% set(gca,'XTickLabel',indicatif2);
% set(gca,'XTick',[1:1:34]);
% xlabel('Coefficient Leader');
% ylabel('Temps (secondes)');
% V=axis;
% axis([0 V(2) 0 5]);
% 
% i=1; %r�initialisation du numero de fichier
% %%% Boucle de chargement des fichiers d'orientation FINALE
% for j = indicatif
%     j1=strcat('Vfinalemoy',j,'_0.csv');
%     j1=j1{1};
%     eval(sprintf('final%g=csvread(%c%s%c);',i,39,j1,39));
%     i=i+1;
% end
%
% k=size(indicatif,2); %nombre de fichier
%
% vio=[];
% for i=1:k
%     for j=1:size(eval(sprintf('final%g',i)),1)
%         %%% Calcul du vecteur d'orientation moyen
%         teta(j)=acos(eval(sprintf('final%g(j,1)',i))/sqrt(eval(sprintf('final%g(j,1)',i))^2+eval(sprintf('final%g(j,2)',i))^2));
%         if asin(eval(sprintf('final%g(j,2)',i))/sqrt(eval(sprintf('final%g(j,1)',i))^2+eval(sprintf('final%g(j,2)',i))^2))<0
%             teta(j)=-teta(j);
%         end
%     end
%     % %%% Trac� par histogramme, pr�f�rez les violin plus bas
%     % figure('name',sprintf('Distribution des orientations finales de groupe n� %s',num2str(k),'NumberTitle','off');
%     % xscal=floor(30*max(teta)-min(teta)/(pi/4));
%     % hist(teta,xscal);
%     % V=axis;
%     % axis([-pi/8 pi/2 V(3) V(4)]);
%     vio=[vio teta.'];
% end
% 
% %figure(h3); % si utilisation du module en entier
% %figure('name',sprintf('Violin plot des distributions d orientations de groupe � param�tre fix�, groupe %g', N),'NumberTitle','off'); %si utilisation de cette partie du module uniquement
% violin(vio,'xlabel',indicatif,'facecolor','y','facealpha',0.4);
% set(gca,'XTickMode','auto');
% set(gca,'XTickLabel',indicatif);
% set(gca,'XTick',[1:1:34]);
% xlabel('Leading coefficient normalized');
% ylabel('Orientation (radian)');
% V=axis;
% axis([V(1) V(2) -pi/8 pi/3]);

%%
%%%%%%%%%%%%%% MODULE 3 : Animation dynamique du vecteur vitesse moyen %%%%%%%%%%
%%% Ce module permet l'obtention d'une animation de l'orientation prise au
%%% cours de chaque simulations par le groupe dans son ensemble
% cd 'C:\Users\******'
% vit=csvread('vecteurvitessemoyen10_0.csv');
% figure;
% x=-pi:0.01:pi;
% plot(cos(x),sin(x))
% hold on;
% for i=1:size(vit)
%     plot(vit(i,1),vit(i,2),'b.');
%     axis([-1.1 1.1 -1.1 1.1]);
%     pause(.01)
% end
%%
%%%%%%%%%%%%%%%%%%% MODULE 4 : Violin plot concat�n� %%%%%%%%%%%%%%%%%%%%%%%
%%% Ce module a �t� utilis� pour l'obtention de la corr�lation entre d�rive
%%% du groupe et influence du leader normalis�e par la taille de groupe, il
%%% n�cessite les donn�es de simulations pour des tailles de groupes
%%% diff�rentes (ici 20, 30, 40, 50)
% ind=[10,11.5,13.0,15.4,18.0,20.5,23.7,27.4,30,36.5,42.1,48.6,56.0,64.9,75.0,86.6,100];
% ind2={'10','11.5','13.0','15.4','18.0','20.5','23.7','27.4','30','36.5','42.1','48.6','56.0','64.9','75.0','86.6','100'};
% ind20=[];
% ind30=[];
% ind40=[];
% ind50=[];
%
% %%% Boucle de normalisation des coefficient leader
% for i = ind
%     ind20=[ind20 i/20]; 
%     ind30=[ind30 i/30];
%     ind40=[ind40 i/40];
%     ind50=[ind50 i/50];
% end
%
% indtot=[ind20 ind30 ind40 ind50]; %r�union de tous les coefficients
% normalis�s
% IND=sort(indtot); %mise dans l'ordre croissant des coefficients
% %%% Retrait des coefficients en double vu � l'oeil, n�cessite
% l'amelioration du syst�me par automatisation, i.e. ATTENTION en fonction
% des tailles mises pas forc�ment les m�mes indices � enlever de la liste,
% ATTENTION 2 : le retrait d'un indice i d�cale tous les suivants de -1,
% les enlever dans l'ordre d�croissant pour �viter de reregarder le
% nouvel indice des suivants � enlever (ce qui a �t� fait ici)
% IND(23)=[]
% IND(47)=[]
% IND(59)=[]
% 
% %%% On charge ici les fichiers efficacit�, temps et orientation finale
% dans l'ordre de la liste totale dans chaque dossier associ�s aux tailles
% de groupes
% for i = ind
%     j=ind2{find(ind==i)};
%     cd C:\Users\*******\N20lf
%     i2=strcat('efficacit+�',j,'_1.csv');
%     eval(sprintf('effi%g=csvread(%c%s%c);',find(IND==i/20),39,i2,39))
%     i2=strcat('temps',j,'_1.csv');
%     eval(sprintf('temps%g=csvread(%c%s%c);',find(IND==i/20),39,i2,39))
%     i2=strcat('Vfinalemoy',j,'_1.csv');
%     eval(sprintf('final%g=csvread(%c%s%c);',find(IND==i/20),39,i2,39))
% 
%     cd C:\Users\********\N30lf
%     i2=strcat('efficacit+�',j,'_1.csv');
%     eval(sprintf('effi%g=csvread(%c%s%c);',find(IND==i/30),39,i2,39))
%     i2=strcat('temps',j,'_1.csv');
%     eval(sprintf('temps%g=csvread(%c%s%c);',find(IND==i/30),39,i2,39))
%     i2=strcat('Vfinalemoy',j,'_1.csv');
%     eval(sprintf('final%g=csvread(%c%s%c);',find(IND==i/30),39,i2,39))
% 
%     cd C:\Users\*********\N40lf
%     i2=strcat('efficacit+�',j,'_1.csv');
%     eval(sprintf('effi%g=csvread(%c%s%c);',find(IND==i/40),39,i2,39))
%     i2=strcat('temps',j,'_1.csv');
%     eval(sprintf('temps%g=csvread(%c%s%c);',find(IND==i/40),39,i2,39))
%     i2=strcat('Vfinalemoy',j,'_1.csv');
%     eval(sprintf('final%g=csvread(%c%s%c);',find(IND==i/40),39,i2,39))
%
%     cd C:\Users\**********\N50lf
%     i2=strcat('efficacit+�',j,'_1.csv');
%     eval(sprintf('effi%g=csvread(%c%s%c);',find(IND==i/50),39,i2,39))
%     i2=strcat('temps',j,'_1.csv');
%     eval(sprintf('temps%g=csvread(%c%s%c);',find(IND==i/50),39,i2,39))
%     i2=strcat('Vfinalemoy',j,'_1.csv');
%     eval(sprintf('final%g=csvread(%c%s%c);',find(IND==i/50),39,i2,39))
% end
% indicatif={'0.20' '0.23' '0.25' '0.26' '0.29' '0.31' '0.33' '0.36' '0.38' '0.39' '0.41' '0.43' '0.45' '0.47' '0.50' '0.51' '0.55' '0.58' '0.59' '0.60' '0.60' '0.65' '0.68' '0.68' '0.73' '0.75' '0.77' '0.79' '0.84' '0.90' '0.91' '0.97' '1' '1.03' '1.05' '1.12' '1.19' '1.22' '1.22' '1.30' '1.37' '1.40' '1.40' '1.50' '1.50' '1.62' '1.62' '1.73' '1.83' '1.87' '1.88' '2' '2.11' '2.16' '2.17' '2.43' '2.50' '2.50' '2.80' '2.89' '3.25' '3.33' '3.75' '4.33' '5'};
% 
% k=size(indicatif,2); % nombre de fichier
% vio=[];
%
% %%% Calcul de l'orientation moyenne finale du groupe
% for i=1:k
%     for j=1:size(eval(sprintf('final%g',i)),1)
%         teta(j)=acos(eval(sprintf('final%g(j,1)',i))/sqrt(eval(sprintf('final%g(j,1)',i))^2+eval(sprintf('final%g(j,2)',i))^2));
%         if asin(eval(sprintf('final%g(j,2)',i))/sqrt(eval(sprintf('final%g(j,1)',i))^2+eval(sprintf('final%g(j,2)',i))^2))<0
%             teta(j)=-teta(j);
%         end
%         %%% On prend la moyenne � part pour pouvoir tracer la courbes
%         qu'on cherchera � ajuster
%         resultat.moy(i,1)=mean(teta);
%         resultat.moy(i,2)=IND(i);
%     end
%     vio=[vio teta.'];
% end
% 
% %%% Trac� des violin 
% figure('name','Violin plot des distributions d orientations de groupe � param�tre fix�','NumberTitle','off');
% violin(vio,'xlabel',indicatif,'facecolor','b','facealpha',0.4);
% set(gca,'XTickMode','auto');
% set(gca,'XTickLabel',indicatif,'Fontsize', 5);
% set(gca,'XTick',[1:1:130]);
% xlabel('Coefficient Leader','Fontsize', 12);
% ylabel('Orientation (radian)','Fontsize', 12);
% V=axis;
% axis([V(1) V(2) -pi/8 pi/3]);
% 
% %%%Trac� de la courbe � ajuster
% figure;
% plot(resultat.moy(:,2),resultat.moy(:,1),'b+')
% set(gca,'XTickMode','auto');
% 
% xlabel('Coefficient Leader','Fontsize', 12);
% ylabel('Orientation (radian)','Fontsize', 12);
% 
%%% Pour l'ajustement utiliser l'appli curve fitting dans l'onglet apps

%%
%%%%%%%%%%%%%%%%%% MODULE 4 : VIOLIN INDIV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%notice: cette partie de script permet la creation de violin en
%subplot(donc pas sur le m�me plot) � partir des orientations individuelles
%prises sur l'ensemble des simulations, puis de tracer les courbes avec et sans leader. Ne pas faire en simultan�
%pr�f�rentiellement (soucis d'affichage dans les 'set gca' et 'handles')
%Param�tres : veiller � changer le repertoire // activer la premiere fois
%puis commenter pour toutes les autres fois les handles de figure,

cd 'C:\Users\Hadrien\Desktop\Leader\xpinf1lead\lead13_100'
%indicatif={'5','10','15','20','25','30','35','40','45'};
%indicatif={'5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'};
indicatif={'50','60','70','80','90','100'};
% h1=figure('name','Courbe moyenne+std')
% h2=figure('name','Courbe std')

i=1;
for j = indicatif
    j1=strcat('efficacit�',j,'_1.csv');
    j1=j1{1};
    eval(sprintf('effi%g=csvread(%c%s%c);',i,39,j1,39));
    i=i+1;
end

set(gca,'XTickMode','auto');
set(gca,'XTickLabel',indicatif);
set(gca,'XTick',[1:1:34]);
xlabel('Taille de groupe');
ylabel('Orientation');

k=size(indicatif,2);

for i =1:k
    subplot(2,3,i)
    violin(eval(sprintf('effi%g',i)),'xlabel',indicatif(i),'facecolor','c','facealpha',0.4);
    V=axis;
    axis([V(1) V(2) -pi/8 pi/3]);
    legend('hide')
end

%%% Trac� des courbes avec et sans leader moyenne+ecart-types ou juste ecart-types
% for j =1:k
%     EFFI(j,1)=mean(eval(sprintf('effi%g',j)))
%     EFFI(j,2)=std(eval(sprintf('effi%g',j)))
% end
% figure(h1);
% errorbar([5,10,15,20,25,30,35,40,45],EFFI(:,1),EFFI(:,2),'color','r');
% hold on;
% legend('Sans leader','Avec leader')
% figure(h2);
% plot([5,10,15,20,25,30,35,40,45],EFFI(:,2),'color','r') 
% hold on;
% legend('Sans leader','Avec leader')
