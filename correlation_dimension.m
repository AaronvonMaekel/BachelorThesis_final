fntsize=11;

%%  EQUATION 1: LOGISTIC MAP
%we have a interval for corDIm for param=3.5699456, vergleichen präzision
%in abhängig von N.
%param=3.5699456;
param= 3.56994567187094490184200515;%hier werden nur 16 stellen genommen.
nmesh=round(linspace(1000,45000,20));
i=1;
corDim=nan(1,length(nmesh));
for n=nmesh
    n 
    tic;
    xt=create_time_series("log",param,0,n,NaN,NaN);
    corDim(i)= correlationDimension(xt,1,1);
    i=i+1;
    toc;
end
interval= [0.4926 0.5024];
%% berechne Abstand zum Intervall
figure(1)
error=nan(1,length(nmesh));
for i=1:length(nmesh)
    if corDim(i) < interval(1)
        error(i)= interval(1) - corDim(i);
    elseif corDim(i) > interval(2)
        error(i)= corDim(i)-interval(2);
    else
        error(i)=0;
    end
end
semilogy(nmesh,error)
xlabel("N")
ylabel("e_{D_2}")
ax = gca;
ax.FontSize = fntsize;
grid on;

%N \approx 8000 ist ausreichend.

%% linear regression analysis 
n=10000;
i=1;
numpoint=20;
dimmi =nan(length(noisemesh),1);
rrange =nan(numpoint,length(noisemesh));
corrint=nan(numpoint,length(noisemesh));

for Noise=noisemesh
    x_i=create_time_series("log",param,Noise,n,NaN,NaN);
    [dimmi(i),rrange(:,i),corrint(:,i)]=correlationDimension(x_i,1,1,"NumPoints",numpoint) ;     
    i=i+1;
end
%% plotting
figure(2)
for i=1:length(noisemesh)
    loglog(rrange(:,i),corrint(:,i),LineWidth=2)
    hold on;
    
end
ax = gca;
ax.FontSize = fntsize;
grid on;
hold off;
legend("SNR=0","SNR=1000","SNR=100","SNR=10","Location","southeast")
xlabel("\epsilon")
ylabel("C(\epsilon)")





%% EQUATION 2 und 3: Zeltabbildung und Benchmark funktion

%% variable initialisation
n=10000;
%chaotic parameter is fixed
pmesh1 = get_model_attributes("benchmark");%parametergitter
pmesh2 = get_model_attributes("tent");



%% berechnung von D_2 beider Funktionen auf dem SNR spektrum 

i=1;
t=nan(1,length(nmesh));
noisemesh= [0 1000 100 10];
corDim1=nan(length(noisemesh),length(pmesh1));
corDim2=nan(length(noisemesh),length(pmesh2));

for Noise=noisemesh
    j=1;
    
    for p1=pmesh1
    
        xt=create_time_series("benchmark",p1,Noise,n,NaN,NaN);
        
        corDim1(i,j)= correlationDimension(xt,1,1);
        
        
        j=j+1;
    end
    i=i+1
end

i=1;
for Noise=noisemesh
    j=1;
    for p2=pmesh2
    
        xt=create_time_series("tent",p2,Noise,n,NaN,NaN);
        
        corDim2(i,j)= correlationDimension(xt,1,1);
        
        
        j=j+1;
    end
    i=i+1
end
%%
mean(corDim1,2)%durchschnittliche dimension des attraktors
mean(corDim2,2)


figure(3)
plot(pmesh1,corDim1(1,:));
%yline(mean(corDim1(1,:)))
xlabel("\theta")
ylabel("D_2")
ax = gca;
ax.FontSize = fntsize;

figure(4)

xlabel("\theta")
ylabel("D_2")
ax = gca;
ax.FontSize = fntsize;

plot(pmesh2,corDim2(1,:),"g")
hold on;
plot(pmesh2,corDim2(2,:),"r")
hold on;
plot(pmesh2,corDim2(3,:),"b")
hold on;
plot(pmesh2,corDim2(4,:),"k")
xlabel("\theta")
ylabel("D_2")
%yline(1)
legend("SNR=0","SNR=1000","SNR=100","SNR=10",location="east")


mean(corDim2(1,6:end))




%% EQUATION 4: 2D MAP
  n=10000;

pmesh = get_model_attributes("2D");

smesh = unique(pmesh(1,:));
qmesh = pmesh(2,1:length(smesh));
edim=2;
%%

noisemesh=[0 1000 100 10];


cordim2D=nan(length(noisemesh),length(smesh),length(qmesh));
i=1;
for Noise=noisemesh
    tic;
    Noise
    dimmi = nan(1,length(pmesh));
    parfor A=1:length(pmesh)        
        x_i=create_time_series("2D",[pmesh(1,A) pmesh(2,A)],Noise,n,NaN,NaN);
        x_i=x_i(1,:);
            dimmi(A)= correlationDimension(x_i,1,edim);
        %else
        %[~,timelag,~]= phaseSpaceReconstruction(x_i,[],edim);
         %    x_i= x_i(1:n);
         %   
        %    dimmi(A)= correlationDimension(x_i,timelag,edim);    
       % end
    end
    
    cordim2D(i,:,:)= reshape(dimmi,length(smesh),[]);

       
    i=i+1;
    toc
    
end
%% plotting
mean(cordim2D,[2 3])
for i=1:4
    line_ratio=18;
    fntsize=11;
    figure(4+i) 
    img_data= reshape(cordim2D(i,:,:),length(smesh),[]);%is das korrekte dimension?
    lim=[min(img_data,[],"all") max(img_data,[],"all")];
    
    t = tiledlayout(line_ratio,1);
    
    nexttile
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    
    imagesc(smesh,qmesh(1),img_data(1,:))
    clim(lim)
    ax = gca;
    xticks([])
    % Set where ticks will be
    ax.YTick = [qmesh(1)];
    % Set TickLabels;
    ax.YTickLabel = {'0'};
    ax.FontSize = fntsize;
    nexttile(2,[line_ratio-1 1]);
    imagesc(smesh,qmesh(2:end),img_data(2:end,:))
    ax = gca;
    ax.YTick = [qmesh(2:2:end)];
    ax.FontSize = fntsize;
    ytickformat('%.2f')
    
    clim(lim)
    xlabel('s') 
    ylabel('ω') 

    h = colorbar;
    h.Label.String = "D_2";
    h.Label.FontSize = fntsize;
    h.Label.Position(1) = [2.98];
    h.Label.Rotation = 360;
    h.Label.VerticalAlignment = "bottom";
    
end
%% error comparison
figure(9)
cordim0= reshape(cordim2D(1,:,:),length(smesh),[]);
cordim1000= reshape(cordim2D(2,:,:),length(smesh),[]);
cordim100= reshape(cordim2D(3,:,:),length(smesh),[]);
cordim10= reshape(cordim2D(4,:,:),length(smesh),[]);
plot([ mean(cordim0(1,:)) mean(cordim1000(1,:)) mean(cordim100(1,:)) mean(cordim10(1,:))])
hold on;
plot([mean(cordim0(2:end,:),"all") mean(cordim1000(2:end,:),"all") mean(cordim100(2:end,:),"all") mean(cordim10(2:end,:),"all")]);
hold off;
legend("Entkopplung","normal",Location="southeast")
ylabel("\langle D_2 \rangle")
%yline(1)
xticklabels(["0" "1000" "100" "10"])
xlabel("SNR")
xticks(1:length(noisemesh))
ax = gca;
ax.FontSize = fntsize;
