 fntsize=11;

%% Histogramm der logistischen Gleichung

xt = create_time_series("log",3.89,0,1000000,NaN,NaN);
figure(1)
histogram(xt,100)
ax = gca;
ax.FontSize = fntsize;
xlabel("x")
ylabel("M(I_k)")


%% Auswirkung von N auf die Schnelligkeit des Rosenstein-Algorithmus 
pmesh= [1.9 1.5 1.1 1.0002];
nmesh=round(logspace(3,6,30)) ;
t=nan(length(pmesh),length(nmesh));

j=1;
for p=pmesh
    p
    i=1;
    xt=create_time_series("tent",p,0,nmesh(end),NaN,NaN);
    for n=nmesh
        n
        tic
        [~] = lyapunovExponent(xt(1:n),1,1,2);
        t(j,i) = toc
        if t(j,i) > 60
            break;
        end
    
        i=i+1;
    end
    j=j+1;
end
%%
 for i=1:length(pmesh)
    semilogx(nmesh,t(i,:))
    hold on;
end
hold off;
xlabel("N")
ylabel("Sekunden")
legend("p=1.9","p=1.5","p=1.1","p=1.0002",Location="northwest")


%% Vorher Nachher Bilder einer Phasenraumrekonstruktion  für das 2-dimensionale System.
%maybe choose a better parameter
pmesh=get_model_attributes("2D");
xt = create_time_series("2D",pmesh(:,end-1) ,0,10000,NaN,NaN);

[xtt,lag,edim] = phaseSpaceReconstruction(xt(1,:),1,2);
% figure(1)
% subplot(1,2,1)
% scatter(xt(1,:),xt(2,:),".")
% xlabel("x_1(t)")
% ylabel("x_2(t)")
% subplot(1,2,2)
% %scatter(xt(1,1:end-1),xt(1,2:end),".")
% %xlabel("x_1(t)")
% %ylabel("x_1(t+1)")
% %figure(2)
% scatter3(xtt(:,1)',xtt(:,2)',xtt(:,3)',".")
% xlabel("x_1(t)")
% ylabel("x_1(t+1)")
% zlabel("x_1(t+2)")
figure(2)
scatter(xt(1,:),xt(2,:),".")
xlabel("x_1(t)")
ylabel("x_2(t)")
ax = gca;
ax.FontSize = fntsize;
grid on;
figure(3)
scatter(xtt(:,1)',xtt(:,2)',".")
xlabel("x_1(t)")
ylabel("x_1(t+1)")
zlabel("x_1(t+2)")
grid on;
ax = gca;
ax.FontSize = fntsize;

%% Test der Breite des Attraktors des 2D-Systems
pmesh = get_model_attributes("2D");
breite = nan(1,length(pmesh));
i=1;
for p=pmesh
    xt = create_time_series("2D",p,0,100000,NaN,NaN);
    breite(i) = mean(abs(xt(1,:)-xt(2,:)));
    i=i+1;
end

smesh = unique(pmesh(1,:));
qmesh = pmesh(2,1:length(smesh));
breite= reshape(breite,length(smesh),[]);

imagesc(smesh,qmesh,breite)
colorbar
%% Veranschaulichung des Intermittenzeffekts für hohe omega?
pmesh = get_model_attributes("2D");
smesh = unique(pmesh(1,:));
qmesh = pmesh(2,1:length(smesh));
s_start = smesh(1)
s_end = smesh(end)
q_start = qmesh(1)
q_end= qmesh(end)
figure(234)
% subplot(1,2,1)
xt = create_time_series("2D",[smesh(1) qmesh(end)],0,10000,NaN,NaN);
scatter(xt(1,:),xt(2,:),".")
% subplot(1,2,2)
% xt = create_time_series("2D",[smesh(end)  qmesh(end)],0,10000,NaN,NaN);
% scatter(xt(1,:),xt(2,:),".")
ax = gca;
ax.FontSize = fntsize;
grid on;
xlabel("x_1(t)")
ylabel("x_2(t)")
%% Approximation des Wahrscheinlichkeitmaßes

figure(1)
 subplot(1,2,1)
xt = create_time_series("2D",[smesh(1) qmesh(10)],0,10000,NaN,NaN);
histogram2(xt(1,:),xt(2,:),100)
 subplot(1,2,2)
 xt = create_time_series("2D",[smesh(1)  qmesh(end)],0,10000,NaN,NaN);
histogram2(xt(1,:),xt(2,:),100)
ax = gca;
ax.FontSize = fntsize;
%%   kleinste Quadrate Schätzung mangelhaft durch "grobes" Parametergitter.
n=1000;
%pmesh= get_model_attributes("log");
pmesh=3:0.001:4;
x_orig = create_time_series("log",3.9015,0,n,0,NaN);
x_orignoise = create_time_series("log",3.9015,10,n,0,NaN);
mse=nan(1,length(pmesh));
mse1000=nan(1,length(pmesh));
i=1;
for p=pmesh
    xt=create_time_series("log",p,0,n,0,NaN);
    mse(i)= mean(abs(xt(2:end)-x_orig(2:end)).^2);
    mse1000(i) = mean(abs(xt-x_orignoise).^2);
    i=i+1;
end
figure(1)

plot(pmesh,mse,"b")
hold on;

%plot(pmesh,mse1000,"r")
hold off;
ylim([0 0.3])
xlabel("\theta")
ylabel("\phi_{MSE}")
ax = gca;
ax.FontSize = fntsize;


%% plotte alle drei funktionen

xmesh=0:0.005:1;
param_log = 3.9;
param_tent = 1.95;
param_bench = 2.2;
ymesh=nan(3,length(xmesh))% log tent benchmark
i=1
for x=xmesh
ymesh(1,i) = param_log*x*(1-x);
ymesh(2,i) = param_tent*min([x 1-x]);
ymesh(3,i) = (1-2^(1-param_bench))^-1 *(1 - x^param_bench- (1- x)^param_bench);
i=i+1;
end
figure(12)

plot(xmesh,ymesh(1,:),"b",LineWidth=1.5);
hold on;
plot(xmesh,ymesh(2,:),"g",LineWidth=1.5);
hold on;
plot(xmesh,ymesh(3,:),"r",LineWidth=1.5);
hold on
plot(xmesh,xmesh,"k");
hold off;
legend("logistische Gleichung","Zeltabbildung","Benchmark-System",location="south")
%% how does the benchmark look for all theta?
% figure(13)
% for ii=2:0.1:3
% i=1;
% for x=xmesh
% imesh(i) = (1-2^(1-ii))^-1 *(1 - x^ii- (1- x)^ii);
% i=i+1;
% end
% plot(xmesh,imesh)
% hold on;
% end
% plot(xmesh,xmesh)