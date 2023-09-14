function pmesh = get_model_attributes(func)%speichert die ausgewählten Parametergitter der Systeme
epsilon= 10e-8;%damit wir nicht die Ränder der Parameterbereiche benutzen
twoD= false;%dient zur klassifizierung des 2D Beispiels
    switch func
        case "log"
            p_start=3.57;
            prec=0.001;
            p_end=4;
        case "tent"
            p_start=1.1;% der Attraktor wird sonst zu klein um ihn zu untersuchen
            prec=0.02;
            p_end=2;
        case "benchmark"
            p_start=2;
            prec=0.01;
            p_end=3;
        case "2D"
            twoD= true;
            p_start1=1.6181;
            p_end1=2;
            prec1=(p_end1-p_start1)/20;                       
            p_start2=0-+epsilon;  %wir möchten den Spezialfall w=0 , (Entkopplung)     
            
            p_end2=1/(2*p_end1); %0.2500        
            prec2=(p_end2-p_start2)/20;
            
           

    
    end
    if twoD ==false
        pmesh=(p_start+epsilon):prec:(p_end-epsilon);
    else
        pmesh1 = (p_start1+epsilon):prec1:(p_end1-epsilon);
        pmesh2 =(p_start2+epsilon):prec2:(p_end2-epsilon);               
        pmesh=[repelem(pmesh1,length(pmesh2));repmat(pmesh2,1,length(pmesh1))];
    end
end