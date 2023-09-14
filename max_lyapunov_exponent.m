 fntsize=11;

%% EQUATION 1: LOGISTIC EQUATION

%% checking stability while changing p
n=1000000;
pmesh = 3.4:0.0001:4;
lyaexp=nan(length(pmesh),1);

parfor A=1:length(pmesh)
    
    x_i=create_time_series("log",pmesh(A),0,n,NaN,NaN);
    %method1
    lyaexp(A)=(1/n)*sum(log(abs(pmesh(A)-x_i.*(2*pmesh(A)))))
   
end

%% markiere bereich zum reinzoomen.

figure(1)
area([3.7185 3.7367 ],[ 0.5 0.5 ],FaceColor="g",LineStyle="none")
hold on;
area([3.7185 3.7367 ],[ -0.2 -0.2 ],FaceColor="g",LineStyle="none")
plot(pmesh,lyaexp,"b")
hold on;
yline(0)
xlabel("\theta")
ylabel("λ_{max}")
ax = gca;
ax.FontSize = fntsize;
xlim([3.4 4])
hold off;
%% zooming in.

pmesh = 3.7185:3e-6:3.7367;
lyaexp_zoom=nan(length(pmesh),1);

parfor A=1:length(pmesh)
    
    x_i=create_time_series("log",pmesh(A),0,n,NaN,NaN);
    %method1
    lyaexp_zoom(A)=(1/n)*sum(log(abs(pmesh(A)-x_i.*(2*pmesh(A)))))
   
end


figure(2)
plot(pmesh,lyaexp_zoom,"b")
hold on;
yline(0)
xlabel("\theta")
ylabel("λ_{max}")
hold off;
ax = gca;
ax.FontSize = fntsize;
xlim([pmesh(1) pmesh(end)])


%% EQUATION 3: BENCHMARK FUNCTION

%%  variable initialisation
n=1000000;

p=2.56;% wir fixieren einen Parameter 
d0 =1e-15;%bestimmt durch präzision der matlab dezimaldarstellung plus Überprüfung, dass initialwerte verschieden sind
Tau = 48;
pmesh = get_model_attributes("benchmark");


%% checking parameters Tau and p

 Taumesh = 40:60;%round(linspace(10,50,40));

 tau_p_err=nan(length(Taumesh),length(pmesh));

 parfor pm=1:length(pmesh)
 
    tauerr= nan(length(Taumesh),1);
    i=1;
    for Tau_p=Taumesh
        x_i=create_time_series("benchmark",p,0,n,NaN,NaN);    
        
        %method2
        init_value=x_i(1);
        terr=0;
        %wir indexieren die Implementierung von diesem Algorithmus anderst,
        %indem wir j mit tau schritten implementieren, dies ändert aber
        %nichts am Algorithmus
        for j = 1:Tau_p:n
            if j+Tau_p > n
                break;
            end 
            %WICHTIG, KEIN CUTOFF, DA WIR ZWEI TRAJEKTORIEN VERGLEICHEN
            if (init_value + d0)>=1
                xi_disturbed=create_time_series("benchmark",pmesh(pm),0,Tau_p,0,(init_value - d0)) ;     
                
            else
                xi_disturbed=create_time_series("benchmark",pmesh(pm),0,Tau_p,0,(init_value + d0)) ;     
            end
            abstand_tau = abs(x_i(:,j+Tau_p-1)-xi_disturbed(:,end));    
            terr = terr + log(norm(abstand_tau./d0));               
            init_value = x_i(:,j);
           

        end
        tauerr(i) =(1/j)*terr;
        i=i+1;
    end
    tau_p_err(:,pm)=tauerr;

end
tau_p_err= abs(tau_p_err-log(2));
%% seems consistent
figure(3)
imagesc(pmesh,Taumesh,tau_p_err)
set(gca, 'Colorscale','log')
colorbar
xlabel("\theta")
ylabel("\tau")
ax = gca;
ax.FontSize = fntsize;
%deshalb wählen wir tau=48

%% checking parameter n



nmesh=round(logspace(2,6,100));
nerror1=nan(length(nmesh),1);
nerror2=nan(length(nmesh),1);
nerror3=nan(length(nmesh),1);

parfor N=1:length(nmesh)

    x_i=create_time_series("benchmark",p,0,nmesh(N),NaN,NaN);
    
    %method1 
    nerror1(N)=(1/nmesh(N))*sum(log(abs((1-2^(1-p))^-1 *(-p*x_i.^(p-1) + p*(-x_i + 1).^(p-1)))));
    
    %method2
    nerr=0;
    init_value=x_i(1);
    for j = 1:Tau:nmesh(N)
            if j+Tau > nmesh(N)
                break;
            end 
            if (init_value + d0)>=1
                xi_disturbed=create_time_series("benchmark",p,0,Tau,0,(init_value - d0)) ;     

            else
                xi_disturbed=create_time_series("benchmark",p,0,Tau,0,(init_value + d0)) ;     
            end
            abstand_tau = abs(x_i(:,j+Tau-1)-xi_disturbed(:,end));    
            nerr = nerr + log(norm(abstand_tau./d0));               
            init_value = x_i(:,j);
           

     end
    nerror2(N) =(1/j)*nerr;
    

    

    %method3
    nerror3(N)= lyapunovExponent(x_i,1,1,1);
   

end
%% plotte comparison der drei funktionen in abhängigkeit zu N
figure(4)
loglog(nmesh,abs(nerror1-log(2)),"g")
hold on;
loglog(nmesh,abs(nerror2-log(2)),"r")
hold on;
loglog(nmesh,abs(nerror3-log(2)),"b")
hold off;
legend("Berechnung über die Ableitung ","Berechnung über Störungsvektoren","Rosensteinalgorithmus",location="southwest")
xlabel("N")
ylabel("e_{\lambda_{max}}")
ax = gca;
ax.FontSize = fntsize;

%% adding noise
noisemesh=[0 1000 100 10];
nerror1=nan(length(noisemesh),1);
nerror2=nan(length(noisemesh),1);
nerror3=nan(length(noisemesh),1);
ldiv=nan(5,length(noisemesh));%TESTING


for Noise=1:length(noisemesh)
    
  
    x_i=create_time_series("benchmark",p,noisemesh(Noise),n,NaN,NaN);
    
    %method1
    
    nerror1(Noise)=(1/n)*sum(log(abs((1-2^(1-p))^-1 *(-p*x_i.^(p-1) + p*(-x_i + 1).^(p-1)))));
    
    %method2
    init_value=x_i(1);
    nerr=0;
    for j = 1:Tau:n
            if j+Tau > n
                break;
            end 
            if (init_value + d0)>=1
                xi_disturbed=create_time_series("benchmark",p,0,Tau,0,(init_value - d0)) ;     

            else
                xi_disturbed=create_time_series("benchmark",p,0,Tau,0,(init_value + d0)) ;     
            end
            abstand_tau = abs(x_i(:,j+Tau-1)-xi_disturbed(:,end));    
           
            nerr = nerr + log(norm(abstand_tau./d0));               
            init_value = x_i(:,j);
           

     end
    nerror2(Noise) =(1/j)*nerr;
    

    %method3
        [nerror3(Noise),~,ldiv(:,Noise)]= lyapunovExponent(x_i,1,1,1);
   

end
%%
figure(5)
semilogy(abs(nerror1-log(2)),"g",Marker="+")
hold on;
semilogy(abs(nerror2-log(2)),"r",Marker="+")
hold on;
semilogy(abs(nerror3-log(2)),"b",Marker="+")
hold off;
legend("Ableitungsverfahren","Abstandsverfahren","Rosenstein",location="northwest")
xlabel("SNR")
xticks(1:length(noisemesh))
xticklabels(["0" "1000" "100" "10"])
ylabel("e_{\lambda_{max}}")
ax = gca;
ax.FontSize = fntsize;


% %% APPENDIX :checking stability while changing p
% 
% 
% lyaexp1=nan(length(pmesh),1);
% lyaexp2=nan(length(pmesh),1);
% lyaexp3=nan(length(pmesh),1);
% 
% 
% 
% 
% 
% n=100000;
% parfor A=1:length(pmesh)
    %     x_i=create_time_series("benchmark",pmesh(A),0,n,0,NaN);
% 
%     %method1
%     lyaexp1(A)=(1/n)*sum(log(abs((1-2^(1-pmesh(A)))^-1 *(-pmesh(A)*x_i.^(pmesh(A)-1) + pmesh(A)*(-x_i + 1).^(pmesh(A)-1)))));
% 
%     %method2
%     init_value=x_i(1);
%     lya=0;
%     for j = 1:Tau:n
%         xi_disturbed=create_time_series("benchmark",pmesh(A),0,Tau,0,(init_value + abstand_zero)) ;     
%         abstand_tau = abs(x_i(:,end)-xi_disturbed(:,end));    
%         lya = lya + 1/(n)*log(norm(abstand_tau./abstand_zero));    
%         if j*Tau <n
%             init_value = x_i(:,j*Tau);
%         end
%     end
%     lyaexp2(A) =lya;
% 
%     %method3
%     lyaexp3(A)= lyapunovExponent(x_i,1,1,2);
% 
% 
% end

% figure(6)
% semilogy(pmesh,abs(lyaexp1-log(2)).^2,"g")
% hold on;
% semilogy(pmesh,abs(lyaexp2-log(2)).^2,"r")
% hold on;
% semilogy(pmesh,abs(lyaexp3-log(2)).^2,"b")
% 
% legend("Ableitungsverfahren","Abstandsverfahren","Rosenstein","Location","southwest")
% xlabel("\theta")
% ylabel("e_\theta")
% 
%bekommen stabile ergebnisse, nicht weiter relevant

%% EQUATION 2: TENT MAP
n=100000;
pmesh = get_model_attributes("tent");
%% calculate lyapunov exactly 

lyaexp = zeros(length(pmesh),1);
parfor A=1:length(pmesh)
    lyaexp(A)=log(abs(pmesh(A)));%exakt berechenbar, da unabhängig von x_i
end

figure(6)
plot(pmesh,lyaexp,"g")
hold on;
xlabel("\theta")
ylabel("λ_{max}")
ax = gca;
ax.FontSize = fntsize;
%% use Rosenstein at different NoiseLevels


noisemesh=[0 1000 100 10];

lyaexp3=nan(length(noisemesh),length(pmesh));

i=1;
for Noise=noisemesh
    tic;
    Noise
    lya3 = nan(1,length(pmesh));
    parfor A=1:length(pmesh)      %PARFOR HIER  
        x_i=create_time_series("tent",pmesh(A),Noise,n,NaN,NaN);        

            lya3(A)= lyapunovExponent(x_i,1,1,1);

    end

    lyaexp3(i,:)=lya3;
    i=i+1;
    toc
end
schaetzer_theta=exp(lyaexp3);

%%
figure(7)
semilogy(pmesh,abs(lyaexp3(1,:)-lyaexp'),"g")
hold on;
semilogy(pmesh,abs(lyaexp3(2,:)-lyaexp'),Color="r")
hold on;
semilogy(pmesh,abs(lyaexp3(3,:)-lyaexp'),"b")
hold on;
semilogy(pmesh,abs(lyaexp3(4,:)-lyaexp'),"k")
hold on;
xlabel("\theta")
ylabel("e_{\lambda_{max}}")
legend("SNR=0","SNR=1000","SNR=100","SNR=10",location="southeast")
ax = gca;
ax.FontSize = fntsize;
xlim([pmesh(1) pmesh(end)])

figure(8)

avg_e(1) = mean(abs(lyaexp3(1,:)-log(pmesh)));
avg_e(2) = mean(abs(lyaexp3(2,:)-log(pmesh)));
avg_e(3) = mean(abs(lyaexp3(3,:)-log(pmesh)));
avg_e(4) = mean(abs(lyaexp3(4,:)-log(pmesh)));

semilogy(avg_e,"b",Marker="+")
hold off;
xlabel("SNR")
xticks(1:length(noisemesh))
xticklabels(["0" "1000" "100" "10"])
ylabel("\langle e_{\lambda_{max}} \rangle")
ax = gca;
ax.FontSize = fntsize;



%% Berechnung des Schätzer \hat\theta
figure(9)
semilogy(pmesh,abs(schaetzer_theta(1,:)-pmesh),"g")
hold on;
semilogy(pmesh,abs(schaetzer_theta(2,:)-pmesh),Color="r")
hold on;
semilogy(pmesh,abs(schaetzer_theta(3,:)-pmesh),"b")
hold on;
semilogy(pmesh,abs(schaetzer_theta(4,:)-pmesh),"k")
hold on;
xlabel("\theta")
ylabel("e_\theta")
legend("SNR=0","SNR=1000","SNR=100","SNR=10",location="southeast")
"worst case"
max(abs(schaetzer_theta(1,:)-pmesh))
max(abs(schaetzer_theta(2,:)-pmesh))
max(abs(schaetzer_theta(3,:)-pmesh))
max(abs(schaetzer_theta(4,:)-pmesh))
"average case"
mean(abs(schaetzer_theta(1,:)-pmesh))
mean(abs(schaetzer_theta(2,:)-pmesh))
mean(abs(schaetzer_theta(3,:)-pmesh))
mean(abs(schaetzer_theta(4,:)-pmesh))


%% EQUATION 4: 2D MAP

n=1000000;

pmesh = get_model_attributes("2D");

smesh = unique(pmesh(1,:));
qmesh = pmesh(2,1:length(smesh));
edim=2;
         

noisemesh=[0 1000 100 10];

%23min for 20x20 and one noise level.
lya2D=nan(length(noisemesh),length(smesh),length(qmesh));
i=1;
for Noise=noisemesh
    tic;
    Noise
    lya3 = nan(1,length(pmesh));
    parfor A=1:length(pmesh)        
        A
        x_i=create_time_series("2D",[pmesh(1,A) pmesh(2,A)],Noise,n,NaN,NaN);
        x_i= x_i(1,:);
       lya3(A)= lyapunovExponent(x_i,1,1,edim);

    end
    
    lya2D(i,:,:)= reshape(lya3,length(smesh),[]);

       
    i=i+1;
    toc
    
end
%% plotting
for i=1:4
    line_ratio=18;
    
    figure(9+i)
    if i==1
    img_data= reshape(lya2D(i,:,:),length(smesh),[]);
    else
        img_data= (reshape(lya2D(1,:,:),length(smesh),[])-reshape(lya2D(i,:,:),length(smesh),[]))^2;
    end
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
    ax.FontSize = fntsize;

    ax.YTick = [qmesh(2:2:end)];
    ytickformat('%.2f')
    colorbar
    clim(lim)
    xlabel('s') 
    ylabel('ω') 
    
end


%% error comparison
figure(14)
lya2D0= reshape(lya2D(1,:,:),length(smesh),[]);
lya2D1000=abs(reshape(lya2D(2,:,:),length(smesh),[])-lya2D0);
lya2D100= abs(reshape(lya2D(3,:,:),length(smesh),[])-lya2D0);
lya2D10= abs(reshape(lya2D(4,:,:),length(smesh),[])-lya2D0);
semilogy(1:3,[ mean(lya2D1000(1,:)) mean(lya2D100(1,:)) mean(lya2D10(1,:))],Marker="+")
hold on;
semilogy(1:3,[ mean(lya2D1000(2:end,:),"all") mean(lya2D100(2:end,:),"all") mean(lya2D10(2:end,:),"all")],Marker="+");
hold off;
xlabel("SNR")
ylabel("Abweichung")
legend("Entkopplung","Normalfall",Location="southeast")
xticklabels([ "1000" "100" "10"])
xticks(1:(length(noisemesh)-1))


% %% linear regression analysis
% i=1;
% regrange=[1 5];
% dimmi =nan(length(noisemesh),1);
% estep =nan(regrange(2)-regrange(1)+1,length(noisemesh));
% ldiv=nan(regrange(2)-regrange(1)+1,length(noisemesh));
% for Noise=noisemesh
%     x_i=create_time_series("2D",[pmesh(1,56) pmesh(2,56)],Noise,n,NaN,NaN);
%         x_i= x_i(1,:);
%         %method3
%         if Noise==0
%              [lya(i),estep(:,i),ldiv(:,i)]= lyapunovExponent(x_i,1,1,edim(56),"ExpansionRange",[regrange(1) regrange(2)])
%         else
%             [lya(i),estep(:,i),ldiv(:,i)]= lyapunovExponent(x_i,1,[],edim(56),"ExpansionRange",[regrange(1) regrange(2)])
% 
%         end
% 
%         i=i+1;
% end
% %%
% figure(56)
% for i=1:length(noisemesh)
%     plot(estep(:,i),ldiv(:,i))
%     hold on;
% 
% end
% grid on;
% hold off;
% legend("Noisefree","SNR1000","SNR100","SNR10","Location","southeast")
