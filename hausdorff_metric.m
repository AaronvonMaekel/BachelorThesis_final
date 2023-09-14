 fntsize=11;
n=1000000;
%wir nehmen den standardmäßigen cutoff von 100
%% EQUATION 1: LOGISTIC EQUATION


%% classify theta in theta c and theta o
pmesh = get_model_attributes("log");
parfor A=1:length(pmesh)   
    x_i=create_time_series("log",pmesh(A),0,n,NaN,NaN);
    
    %wir benutzen methode 1
    lyaexp(A)=(1/n)*sum(log(abs(pmesh(A)-x_i.*(2*pmesh(A)))))
   
end
theta_c= lyaexp>0;
theta_o = ~theta_c;

%%  Algorithmus 1.
empirical_borders= nan(2,length(pmesh));
i=1;
 for p=pmesh
    xt = create_time_series("log",p,0,n,NaN,NaN);
    empirical_borders(:,i)=[max(xt);min(xt)];
    i=i+1;
 end


%hier berechnen wir die durchschnittlichen und worst case fehler der beiden
%grenzfunktionen, dies wird in theta_c und theta_o getrennt.
"chaotischer Attraktor"
"f1"
avg_error_upper_boundary_c =mean(abs(empirical_borders(1,theta_c)-0.25*pmesh(theta_c)))
max_error_upper_boundary_c =max(abs(empirical_borders(1,theta_c)-0.25*pmesh(theta_c)))
"f2"
avg_error_lesser_boundary_c =mean(abs(empirical_borders(2,theta_c)-0.25*(pmesh(theta_c).^2)+1/16*(pmesh(theta_c).^3)))
max_error_lesser_boundary_c =max(abs(empirical_borders(2,theta_c)-0.25*(pmesh(theta_c).^2)+1/16*(pmesh(theta_c).^3)))

"periodische Orbits"
"f1"
avg_error_upper_boundary_o =mean(abs(empirical_borders(1,theta_o)-0.25*pmesh(theta_o)))
max_error_upper_boundary_o =max(abs(empirical_borders(1,theta_o)-0.25*pmesh(theta_o)))

"f2"
avg_error_lesser_boundary_o =mean(abs(empirical_borders(2,theta_o)-0.25*(pmesh(theta_o).^2)+1/16*(pmesh(theta_o).^3)))
max_error_lesser_boundary_o =max(abs(empirical_borders(2,theta_o)-0.25*(pmesh(theta_o).^2)+1/16*(pmesh(theta_o).^3)))

%% fehler an f1 kriterium

avg_error_0 = mean(abs(pmesh(theta_c)-4* empirical_borders(1,theta_c)))
worst_error_0 = max(abs(pmesh(theta_c)-4* empirical_borders(1,theta_c)))

%% histogramme von theta
chaosgut = 3.68;
ordnungs_p = 3.84;
chaosschlecht = 3.767;

%hiermit wurden geeignete Fälle ausgesucht.
% figure(3)
% subplot(2,1,1)
% xt=create_time_series("log",chaosgut,0,n,NaN,NaN);
% histogram(xt,200,BinLimits=[0 1])
% subplot(2,1,2)
% xt=create_time_series("log",chaosgut,100,n,NaN,NaN);
% histogram(xt,200,BinLimits=[0 1])
% figure(4)
% subplot(2,1,1)
% xt=create_time_series("log",chaosschlecht,0,n,NaN,NaN);
% histogram(xt,400,BinLimits=[0 1])
% subplot(2,1,2)
% xt=create_time_series("log",chaosschlecht,100,n,NaN,NaN);
% histogram(xt,400,BinLimits=[0 1])
% figure(5)
% subplot(2,1,1)
% xt=create_time_series("log",ordnungs_p,0,n,NaN,NaN);
% histogram(xt,200,BinLimits=[0 1])
% subplot(2,1,2)
% xt=create_time_series("log",ordnungs_p,100,n,NaN,NaN);
% histogram(xt,200,BinLimits=[0 1])
%% 
mode=2; %1=gut konditioniertes chaos,  2= schlecht konditioniertes Chaos  3=periodischer Orbit;

param=[chaosgut chaosschlecht ordnungs_p];
upper_bound=+0.25*param(mode);
lower_bound=0.25*param(mode)^2-1/16*(param(mode)^3);
tries=100;

 
noisemesh=[round(logspace(3,1,3))]; %wir wenden dies nur auf gestörte Fälle an.
binlength_mesh=logspace(-1,-4,40);
Erfolgsmatrix=zeros(length(noisemesh),length(binlength_mesh));
 
    
for noise=1:length(noisemesh)
    tic
    runde=1;     
    while runde <= tries
        xt = create_time_series("log",param(mode),noisemesh(noise),n,NaN,NaN);
        bin_border =[min(xt) max(xt)];
        jj=1;
        for bn=binlength_mesh
        
            bn_fixed= (bin_border(2)-bin_border(1))/round((bin_border(2)-bin_border(1))/bn);
            [hvalue,hedges]= histcounts(xt,"BinWidth",bn_fixed,"BinLimits",bin_border);

            grenzen=islocalmax([ 0 hvalue 0]);
            grenzen=grenzen(2:end-1);
           
            grosser_rand_bin = find(grenzen,1,'last');
            Erfolgsmatrix(noise,jj)=Erfolgsmatrix(noise,jj)+and(upper_bound>= hedges(grosser_rand_bin),upper_bound < hedges(grosser_rand_bin+1));
            jj=jj+1;     
        end
        runde=runde+1;
    end 
    toc
    
end

%PROBLEM WENN bin limit min(xt)-max(xt): ohne noise wird die echte grenze
%ausgeschlossen da diese nie erreicht wird. 

%%
fntsize=11;
figure(2)

imagesc(1:length(noisemesh),binlength_mesh,flip(Erfolgsmatrix,2)')
clim([0,100])

set(gca,'YScale','log');
ax = gca;
ax.FontSize = fntsize;
ax.XTick = 1:length(noisemesh);
ax.XTickLabel = string(noisemesh);   

yticks(logspace(-5,-1,5))
yticklabels(cellstr(num2str(logspace(-1,-5,5)', '%.0e')))
xlabel('SNR') 
ylabel('Partitionsbreite') 
h = colorbar;
h.Label.String = "Erfolgsrate";
h.Label.FontSize =fntsize;
h.Label.Rotation = 270;
h.Label.VerticalAlignment = "bottom";
h.Ruler.TickLabelFormat =  '%g%%';
h.Ruler.FontSize = fntsize;

%% EQUATION 2: TENT MAP
pmesh = get_model_attributes("tent");

empirical_borders= nan(2,length(pmesh));


i=1;
 for p=pmesh
    xt = create_time_series("tent",p,0,n,100,NaN);
    empirical_borders(:,i)=[max(xt);min(xt)];
    i=i+1;
 end
    

 
figure(1)
mksize=12;
semilogy(pmesh,abs(empirical_borders(1,:)-0.5*pmesh),Marker=".",MarkerSize=mksize)
hold on;
semilogy(pmesh,abs(empirical_borders(2,:)+0.5*(pmesh.^2)-pmesh),Marker=".",MarkerSize=mksize)
hold off;
xlabel("\theta")
ylabel("Abweichung")
ax = gca;
ax.FontSize = fntsize;
legend("e_b","e_a",Location="southeast")
"Zeltabbilung"
avg_error_upper_boundary_c =mean(abs(empirical_borders(1,:)-0.5*pmesh))
avg_error_lesser_boundary_c =mean(abs(empirical_borders(2,:)+0.5*(pmesh.^2)-pmesh))

max_error_upper_boundary_c =max(abs(empirical_borders(1,:)-0.5*pmesh))
max_error_lesser_boundary_c =max(abs(empirical_borders(2,:)+0.5*(pmesh.^2)-pmesh))

%% histogramme an ausgewählten Parametern.
param=[1.1 1.55 1.8];
for mode=1:3
    xt = create_time_series("tent",param(mode),0,n,100,NaN);
    figure(mode+10)
    histogram(xt,100,BinLimits=[min(xt) max(xt)])
    
end
% wir sehen dass sich unser Algorithmus aufgrund der Dichte hier nicht geeignet
% ist.


%% EQUATION 3: BENCHMARK FUNCTION ( IRRELEVANT)

%% EQUATION 4: 2D MAP
n=1000000;

pmesh = get_model_attributes("2D");

smesh = unique(pmesh(1,:));
qmesh = pmesh(2,1:length(smesh));



%% schätze omega durch R und R'
noisemesh= [0 1000 100 10];
error1=nan(length(noisemesh),length(pmesh));
error2=nan(length(noisemesh),length(pmesh));
errortotal=nan(length(noisemesh),length(pmesh));

for Noise=1:length(noisemesh)
    tic;
    Noise
    dimmi = nan(1,length(pmesh));
    for A=1:length(pmesh) %parfor später       
        A
        x_i=create_time_series("2D",[pmesh(1,A) pmesh(2,A)],noisemesh(Noise),n,NaN,NaN);
        diag = x_i(2,:)-x_i(1,:);
        [~,indexR]= max(diag);
        [~,indexRstrich]= min(diag);
        R=x_i(:,indexR);
        Rstrich=x_i(:,indexRstrich);
        omega1= 0.5*R(1);
        omega2= 0.5*Rstrich(2);
        error1(Noise,A) = abs(omega1-pmesh(2,A));
        error2(Noise,A) =abs( omega2-pmesh(2,A));
        if noisemesh(Noise)==0
         omega= min(omega1,omega2);
        else
          omega = mean([omega1 omega2]);
        end
        errortotal(Noise,A) = abs(omega- pmesh(2,A));
        
    end
   
    toc
    
end
%% 
avg_error1= mean(error1,2)
max_error1 = max(error1,[],2)
avg_error2 = mean(error2,2)
max_error2 = max(error2,[],2)
avg_errortotal = mean(errortotal,2)
max_errortotal = max(errortotal,[],2)

figure(134)
semilogy(avg_errortotal)
hold on;
semilogy(max_errortotal)
hold off;
legend("\langle e_\theta \rangle","E_\theta",location="southeast")
xticks(1:length(noisemesh))
xticklabels(["0" "1000" "100" "10"])
xlabel("SNR")
ax = gca;
ax.FontSize = fntsize;

    
%% plotting 

for i=1%:length(noisemesh)
    line_ratio=18;
    fntsize=11;
    figure(9+i)
  
    img_data= reshape( errortotal(i,:),length(smesh),[]);
  
   
    t = tiledlayout(line_ratio,1);
    lim=[min(img_data,[],"all") max(img_data,[],"all")];
    nexttile
  
    t.TileSpacing = 'compact';
   
    %t.Padding = 'compact';

    imagesc(smesh,qmesh(1),img_data(1,:))
    ax = gca;
    xticks([])
    % Set where ticks will be
    ax.YTick = [qmesh(1)];
    % Set TickLabels;
    ax.YTickLabel = {'0'};
    ax.FontSize = fntsize;
    clim(lim)
    nexttile(2,[line_ratio-1 1]);
    imagesc(smesh,qmesh(2:end),img_data(2:end,:))
    ax = gca;
    ax.FontSize = fntsize;
    ax.YTick = [qmesh(2:2:end)];
    ytickformat('%.2f')
   
    xlabel('s') 
    ylabel('ω') 
    h = colorbar;
    h.Label.String = "e_\omega";
    h.Label.FontSize = fntsize;
    h.Label.Rotation = 360;
    h.Label.VerticalAlignment = "bottom";
    %h.Ruler.TickLabelFormat =  
    clim(lim)
    h.Label.Position = [2.8017 0.0047];
    %set(h.Label,{'Position'},{[3.7067 0.0284]})
    h.Label
end


%% approximierte wahrscheinlichkeitsverteilung im Falle des Edgecases
figure(14)
x_i=create_time_series("2D",[smesh(end) qmesh(end)],0,n,NaN,NaN);
       histogram2(x_i(1,:),x_i(2,:),100)
smesh(end)
qmesh(end)
ax = gca;
ax.FontSize = fntsize;
grid on;
xlabel("x_1(t)")
ylabel("x_2(t)")