function xt = create_time_series(func,param_A,SNR,datalength,cutoff,init_v)
    
     if ~isnan(cutoff)
      cut = cutoff;
     else %default ist 100 
      cut = 100;
     end

     switch func
        case "log"
        
            init_val = 0.41;
            xt = zeros(1,cut+datalength);
            f= @(x) param_A*x*(1-x);
            lowerborder=0;%beschränkung des Zustandsraumes, falls durch einen Observationsfehler ein Datenpunkt 
            % den Zustandsraum verlässt, müssen wir den Fehler neu generieren
            upperborder=1;
       
        case "tent"
        
            init_val = 0.52;
            xt = zeros(1,cut+datalength);
            f= @(x) param_A* min(x,1-x);
            lowerborder=0;
            upperborder=1;
        
        case "benchmark"
        
            init_val = 0.41;
            xt = zeros(1,cut+datalength);
            f= @(x) (1-2^(1-param_A))^-1 *(1 - x^param_A - (1- x)^param_A);
            lowerborder=0;
            upperborder=1;
            
        case "2D"
            
            init_val = [0.46;0.51];
            xt = zeros(2,cut+datalength);
            s= param_A(1);
            w= param_A(2);
            ts1 = @(x1) (x1<= 1/s).* s * x1 + (x1> 1/s)*s/(s-1) * (1- x1);
            ts2 = @(x2) (x2<= 1/s).* s * x2 + (x2> 1/s)*s/(s-1) * (1- x2);
            f = @(x)  [(1-w)*ts1(x(1))+w*ts2(x(2))  w*ts1(x(1))+(1-w)*ts2(x(2))] ;  
            lowerborder=0;
            upperborder=1;
        end
    
     if ~isnan(init_v)
      init_val=init_v;
     end
        
     

    


    xt(:,1) = init_val;
        for i = 2:cut+datalength
            xt(:,i) = f( xt(:,i-1));
        end
    xt = xt(:,(cut+1):end);%schneide den Anfang ab.
    
    if SNR ~= 0 
        standardabw = sqrt(var(xt(1,:)));
        noise_standardabw = standardabw/SNR;
        xt_neu=nan(size(xt));%2dim
        for  k = 1:size(xt,1)%berechnen die observationsfehler für jede dimension einzeln, code nimmt an dass zustandsraum entlang jeder dimension gleich ist
            xt_neu(k,:)= xt(k,:) + normrnd(0,noise_standardabw,1,datalength); %generiere observationsfehler für alle observationen     
            while (nnz(xt_neu(k,:)<lowerborder) +nnz(xt_neu(k,:)>upperborder))>0 %index bisschen überfüllt, wir betrachten immer nur observationen welche sich nach 
                %der hinzugabe von Fehlern sich ausßerhalb des
                %Zustandsraumes befinden und generieren diese neu, bis alle
                %passen
                xt_neu(k,logical((xt_neu(k,:)<lowerborder) + (xt_neu(k,:)>upperborder)))= xt(k,logical((xt_neu(k,:)<lowerborder) + (xt_neu(k,:)>upperborder)))+...
                normrnd(0,noise_standardabw,1,nnz(xt_neu(k,:)<lowerborder) +nnz(xt_neu(k,:)>upperborder));          
            end
        end
        xt=xt_neu;  
        if cutoff==0%da wir den initialwert auch gestört haben,dieser bleibt immer Störungsfrei.
            xt(:,1)=init_val;
        end
    end

end

