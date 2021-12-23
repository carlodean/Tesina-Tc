close all
clear all
clc

%% Dati problema
    sm = 2.36e-3;   % m
    sc=sm;
    s = 14.16e-3;   % m
    rho_m = 8e3;    %kg/m^3
    rho_c = 1e3;    %kg/m^3
    Cm = 480;       %J/KgK
    Cc = 1.5e3;     %J/KgK
    km = 12;        %W/mK
    kc = 0.3;       %W/mK
    alpha_m=km/(rho_m*Cm);
    alpha_c=kc/(rho_c*Cc);
    nodi_interni_layer=1;       %commentare se stai testando consistenza della griglia  
    layer=6;
    
    cond_var=false;             % se è false la condizione al contorno è fissa, viceversa con true
    instabilita=false;          % false non valuta l'instabilità, viceversa con true
    consistenza_griglia=false;  % false non valuta la consistenza della griglia, viceversa con true
    
%    for j=0:1                   %commentare questa linea se non si deve testare la consistenza della griglia + end alla fine
%    nodi_interni_layer=4^j;     %commentare questa linea se non si deve testare la consistenza della griglia + end alla fine
    
%% Condizioni
    Ts = 190;       % °C
    T0 = 15;        % °C   
    
%% Numero di nodi e griglia 
    N = (layer-1) + layer*nodi_interni_layer + 2;    %nodi totali presi

    x=linspace(0,s,N);
    dx= x(2)-x(1);
    
%% Matrice dei coefficienti e vettore dei termini noti
    A = zeros(N);
    b = zeros(N,1);
    
%% Condizione metodo esplicito 
    dt_ins = (dx^2)/(2*alpha_m);

%% Passo temporale
    dt=1;

%% Costanti
    phi = (2*dt)/(dx^2*(rho_m*Cm+rho_c*Cc));
    Foxm = alpha_m*dt/dx^2;
    Foxc = alpha_c*dt/dx^2;
    
%% Inizializzazione vettore termini noti
    if(~cond_var)
        b(1) = -190;
        b(end) = -190;
    end   
    
%% Nodi interni
tipo_equazione=0;
cnt=0;
for ir=2:N-1
    for ic = 1:N
        if ir == ic && tipo_equazione == 0
            A(ir,ic)=-2*Foxm;
            A(ir,ic-1) =Foxm;
            A(ir,ic+1) =Foxm;
            cnt=cnt+1;
            
        elseif ir == ic && tipo_equazione==1
            A(ir,ic)=-phi*(kc+km);
            A(ir,ic-1) =phi*(kc);
            A(ir,ic+1) =phi*(km);
            tipo_equazione=0;
            
        elseif ir == ic && tipo_equazione==2
            A(ir,ic)=-2*Foxc;
            A(ir,ic-1) =Foxc;
            A(ir,ic+1) =Foxc;
            cnt=cnt+1;   
            
        elseif ir == ic && tipo_equazione==3
            A(ir,ic)=-phi*(kc+km);
            A(ir,ic-1) =phi*(km);
            A(ir,ic+1) =phi*(kc);
            tipo_equazione=2;
        end
    end
            
    if tipo_equazione==0 && cnt==nodi_interni_layer
        tipo_equazione=3;
        cnt=0;
    elseif tipo_equazione==2 && cnt==nodi_interni_layer
        tipo_equazione=1;
        cnt=0; 
    end         
end

%% Problema in transitorio
%Condizione iniziale
    Tn = ones(N,1) * T0;       %Tn è un vettore con tante tighe quanti sono i nodi
                               % inizializzate al valore T0                         
%Tempo totale
    t_tot = 500; %s

%Parametro per metodo applicato (0 exp, 1/2 C-N, 1 imp)
    par = 1;    

%Numero step (arrotondato per eccesso)
    nt = ceil(t_tot/dt);      %ceil => arrotonda al numero intero superiore

%Vettore tempo
    t = 0:dt:nt*dt;    %vai da 0 a nt*dt con uno step di dt

%Matrice identità 
    I = speye(size(A));
    Imod=I;
    
%Condizoni del primo tipo agli estremi della griglia
    Imod(1,:) = 0; 
    Imod(end,:) = 0;   
    
%Variabili utili
    Ts_var=zeros(nt,1);
    mezzeria = (N+1)/2;
    pol = 0;      %variabile che diventa 1 se è stata raggiunta la temperatura di pol.
    raffr = 0;    %variabile che diventa 1 se è stata raggiunta la temperatura di raffr.
    
for it = 1:nt
    if(cond_var)
        m = 175/1200;   %coefficiente angolare retta condizione variabile
        if dt*it < 1200                         %andamento lineare crescente fino a 20 minuti
            temp_var = 15 + m*(dt*(it-1));      
            b(1) = -temp_var;
            b(end) = -temp_var;
        elseif dt*it >= 1200 && dt*it <= 3000   %andamento costante a 190 °C fino a 50 minuti
            b(1) = -190;
            b(end) = -190;
        elseif dt*it > 3000                     %andamento costante a 15 °C dopo 50 minuti
            b(1) = -15;
            b(end) = -15;
        end 
        Ts_var(it) = -b(1);
    end
    
    bm = (( Imod + (1-par)*A ) *  Tn(:,it)) - b;
    Am =  I - ( par*A );
    % Risolvo il problema moltiplicando l'inversa di Am con bm
    Tn(:,it+1) = Am\bm;
    
    %Tempo affinchè la mezzeria raggiunga 170 °C
    if Tn(mezzeria,it+1) > 170 && pol == 0      
        pol = 1;
        tempo_necessario_pol = dt*it;   
        disp('Il tempo di polimerizzazione per la mezzeria è [s]: ');
        disp(tempo_necessario_pol);
    end
    
    %Tempo affinchè la mezzeria si raffreddi
    if(cond_var)
        if Tn(mezzeria,it+1) < 37 && raffr == 0 && dt*it > 3000 
            raffr = 1;
            tempo_necessario_raffr = (dt*it) - 3000;
            disp('Il tempo di raffreddamento a 37 °C per la mezzeria è [s]: ');
            disp(tempo_necessario_raffr);
        end
    end
        
    figure(1)
    plot(x,Tn(:,it+1),'LineWidth',3);
    grid on
    axis([0 s 0 200]);
    legend('T(t)');
    xlabel('x (m)','FontSize', 12);
    ylabel('T (°C)','FontSize', 12);
    set(gca,'FontSize',14);
    title({'t (s)',(it * dt)},'FontSize',14);   %non c'è¨ it+1 perchè l'istante iniziale t=0 è¨ associato ad it=1.
end

%Andamento della temperatura nel tempo, ai bordi e alla mezzeria
pause
figure(2)
    plot(t,Tn(2,:),t,Tn(mezzeria,:),t,Tn(N-1,:),'LineWidth',3); 
    grid on
    axis([0 t_tot 0 200]);
    legend('T(t) in x = dx','T(t) in x = mezzeria','T(t) in x = s-dx');
    xlabel('t (s)','FontSize', 12);
    ylabel('T (°C)','FontSize', 12);
    set(gca,'FontSize',14);
    
%Andamento temperatura nello spazio in 1/4, 2/4, 3/4 e 4/4 di t_totale
if(~consistenza_griglia)
    pause
    figure(3)
    plot(x,Tn(:,ceil(1/4*nt)),x,Tn(:,ceil(2/4*nt)),x,Tn(:,ceil(3/4*nt)),x,Tn(:,ceil(4/4*nt)),'LineWidth',3);
    grid on
    axis([0 s 0 200]);    
    legend('1/4*t_tot','2/4*t_tot','3/4*t_tot','4/4*t_tot');
    xlabel('x (m)','FontSize', 12);
    ylabel('T (°C)','FontSize', 12);
    set(gca,'FontSize',14); 
end
    
%Andamento condizione al contorno variabile nel tempo
if(cond_var)
    pause
    figure(4)
    plot(t(1:nt),Ts_var,'LineWidth',3);        %t(nt+1) NON serve perchè la condizione al contorno varia in nt passi
    grid on
    axis([0 t_tot 0 200]);
    legend('Temperatura al contorno [°C]');
    xlabel('t (s)','FontSize', 12);
    ylabel('T (°C)','FontSize', 12);
    set(gca,'FontSize',14);
end

%% Instabilità
if(instabilita)
  % Prendiamo un punto nel dominio ed evidenziamo l'instabilità  del
  % Metodo esplicito  
    tic
    for it = 1:nt
    bm = (( Imod + (1-par)*A ) *  Tn(:,it) ) - b;
    Am =  I - ( par*A );
    Tn(:,it+1) = Am\bm;
    end
    disp ('Tempo di simulazione del metodo esplicito: ')
    toc
    
    pause
    figure(5)
    plot(t,Tn(6,:),'r','LineWidth',3);
    grid on
    axis([0 t_tot 0 200]);
    xlabel('t (s)','FontSize', 12);
    ylabel('T (°C)','FontSize', 12);
    set(gca,'FontSize',14);     

    % Metodo implicito
    par=1;
    tic
    for it = 1:nt
    bm = (( Imod + (1-par)*A ) *  Tn(:,it) ) - b;
    Am =  I - ( par*A );
    Tn(:,it+1) = Am\bm;
    end
    disp ('Tempo di simulazione del metodo implicito: ')
    toc
    
    hold on
    plot(t,Tn(6,:),'b','LineWidth',3);
    grid on
    axis([0 t_tot 0 200]);
    xlabel('t (s)','FontSize', 12);
    ylabel('T (°C)','FontSize', 12);
    set(gca,'FontSize',14);
    

    %Metodo C-N
    par=1/2;
    tic
    for it = 1:nt
    bm = (( Imod + (1-par)*A ) *  Tn(:,it) ) - b;
    Am =  I - ( par*A );
    Tn(:,it+1) = Am\bm;
    end
    disp ('Tempo di simulazione del metodo C-N: ')
    toc
        
    hold on
    plot(t,Tn(6,:),'g','LineWidth',3);
    grid on
    axis([0 t_tot 0 200]);
    xlabel('t (s)','FontSize', 12);
    ylabel('T (°C)','FontSize', 12);
    set(gca,'FontSize',14);
    legend('Metodo esplicito','Metodo implicito','Metodo C-N')
    hold off
end

%% Consistenza griglia
if(consistenza_griglia)
   if(j == 0)
       x1=x;
       Tn1=Tn;
   elseif (j==1)
       x2=x;
       Tn2=Tn;
       
       pause
       figure(6)
       plot(x1,Tn1(:,ceil(1/4*nt)),x1,Tn1(:,ceil(2/4*nt)),x1,Tn1(:,ceil(3/4*nt)),'LineWidth',3);
       grid on
       axis([0 s 0 200]);
       xlabel('x (m)','FontSize', 12);
       ylabel('T (°C)','FontSize', 12);
       set(gca,'FontSize',14); 
       
       hold on
       plot(x2,Tn2(:,ceil(1/4*nt)),x2,Tn2(:,ceil(2/4*nt)),x2,Tn2(:,ceil(3/4*nt)),'LineWidth',3);
       grid on
       axis([0 s 0 200]);
       legend('1/4*t_tot(N=13)','2/4*t_tot(N=13)','3/4*t_tot(N=13)','1/4*t_tot(N=31)','2/4*t_tot(N=31)','3/4*t_tot(N=31)');
       xlabel('x (m)','FontSize', 12);
       ylabel('T (°C)','FontSize', 12);
       set(gca,'FontSize',14); 
       hold off
   end
   end
%end