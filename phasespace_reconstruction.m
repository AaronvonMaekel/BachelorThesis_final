%% Konvergenz der Einbettungsdimension
nmesh=round(logspace(2,5,10));
edim=nan(1,length(nmesh));
i=1;
for n=nmesh
    xt=create_time_series("benchmark",1.1,0,n,0,NaN);
[~,~,edim(i)] = phaseSpaceReconstruction(xt,1);
i=i+1;
end

figure(1)
semilogx(nmesh,edim)
%konvergiert schon für kleine Datensätze.


%% EQUATION 4: 2-dimensionales Chaos

n=10000;

pmesh = get_model_attributes("2D");
smesh = unique(pmesh(1,:));
omegamesh = pmesh(2,1:length(smesh));

edim=nan(1,length(pmesh));
i=1;
for p=pmesh
    xt=create_time_series("2D",p,0,n,0,NaN);
    [~,~,edim(i)] = phaseSpaceReconstruction(xt(1,:),1);
    i=i+1;
end

unique(edim)
edim = reshape(edim,length(smesh),[])
%edim=2 wenn omega=0, ansonsten edim=3;

%% TIMELAG Analyse

%% Auswirkung von N auf die Konvergenz des AMI-Algorithmus
pmesh= [1.9 1.5 1.3 1.1];
nmesh=round(logspace(3,6.5,30)) ;
edim=nan(length(pmesh),length(nmesh));
j=1;
for p=pmesh
    
    i=1;
    xt=create_time_series("tent",p,10,nmesh(end),0,NaN);
    for n=nmesh
        
        tic
        [~,lag(j,i),~] = phaseSpaceReconstruction(xt(1:n),[],2);
       
        i=i+1;
    end
    j=j+1;
end

%%
figure(2)
for i=1:length(pmesh)
    semilogx(nmesh,lag(i,:))
    hold on;
end
hold off;
hold off;
xlabel("N")
ylabel("\tau")
legend("\theta=1.9","\theta=1.5","\theta=1.3","\theta=1.1",Location="northeast")
ylim([0 11])
ax = gca;
ax.FontSize = fntsize;
grid on;
%konvergenz unter einer Million

%%
% %%  testing hyptothese
% 
% pmesh = [1.2 1.3 1.6 1.7 1.9];
% n=1000000;
% noisemesh= [0 1000 100 10];
% lya_real=log(abs(pmesh));
% lya_edim1=nan(length(noisemesh),length(pmesh));
% lya_edim2=nan(length(noisemesh),length(pmesh));
% i=1;
% for noise=noisemesh
%     j=1;
%     noise
%     for p = pmesh
%            p
%            x_t=create_time_series("tent",p,noise,n,NaN,NaN);
%            if noise==0
% 
%            lya_edim1(i,j)=lyapunovExponent(x_t,1,1,1);
%            lya_edim2(i,j)=lyapunovExponent(x_t,1,1,2);
%            else
%             lya_edim1(i,j)=lyapunovExponent(x_t,1,[],1);
%             lya_edim2(i,j)=lyapunovExponent(x_t,1,[],2);
%            end
%     j=j+1;
%     end
% i=i+1;
% end
% 
% for i=1:4
%     figure(2+i)
%     semilogy(abs(lya_edim1(i,:)-lya_real))
%     hold on;
%     semilogy(abs(lya_edim2(i,:)-lya_real))
%     hold off;
%     legend("edim1","edim2")
% end
% %without noise, edim1 > edim2 (timelag 1 ofc)
% %with noise: with timelag 1 : edim2> edim1
% %with noise: with timelag calculated : edim2> edim1

