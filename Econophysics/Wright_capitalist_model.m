close all
clc
clear

M = 1e7;                    %Cantidad de Dinero
N = 1000;                    %Cantidad de personas
C = 100000   ;                   %Clases            
w = [10,90];              %Intervalo de salarios
w_avg = (w(2) + w(1))/2;    %Salario promedio
e = zeros(1,N);             %Indice de empleador
W = zeros(1,N);             %Conjunto de empleados
M_n = M/N;                  %Dinero promedio
M_t = 1000;                 %Constante de transacción
t = 10000;                    %Tiempo
M_c = 1e4;                  %Cantidad de dinero Mc (Aún por determinar su significado)
V = 0;                      %Market Value
MT = (0);
clases = zeros(1,C);
c = 0;
MT = (0);
P = (0);            %Potential employers  
p = (0);
A = (0);
B = (0);
X = (0);
S = (0);
sumNk = 0;
sumO = 0;
objetivo = (0);
a = 1/M_c; %a = 1/M_c;
caso = 1                        ;
FSM = zeros(N,1);

IM = zeros(1,N);
IMC = zeros(1,N);
IMW = zeros(1,N);

% Distribución incial del dinero
n = ones(1,N) .* M_n;               %Dinero de cada agente i


tic
    for l = 1:1
        n_0 = n;            %Guarda dinero al principio del año
        % Reglas del sistema SIMULATION RULE SR_1
        for i = 1:N*12
            
            
            a_1 = int16(ceil(rand*N));      %Actor selection rule
            [e,p,W] = H_1(W,a_1,e,n,w_avg,N,P,p);  %Hiring rule
            b_1 = int16(ceil(rand*N));
            [V,n] = E_1(a_1,b_1,n,V,N);         %Expediture rule
            [V,n] = M_1(V,n,W,a_1,e);             %Market sample rule
            [e,W] = F_1(n,N,a_1,e,W,w_avg) ;        %Firing rulle
            n = w_1(N,a_1,e,n,W,w);                %Wage payment rule
            MT(i) = V + sum(n) ;                    %Suma del dinero
               
            nk = histcounts(n,C);
            [counts,centers] = hist(n,C);
            S(i) = entropia(nk,N);                  %Entropia 
            objetivo(i) = o_objetivo(centers,C,nk,a,caso,sumO);
            [A(i),B(i),X(i)] = countclass(e,W,N);   %[Workers,Employers,Unemployed]


    
        end
       for k = 1:N
           IM(l,k) = n(k);   %  - n_0(k);
           if IM(l,k) < 0 
               IM(l,k) = 0;
           end
       end
       for k = 1:N
           if e(k) == 0 && W(k) ~= 0
               IMC(l,k) = IM(l,k);
           else
               IMW(l,k) = IM(l,k);
           end
       end
               
    end
    
       
        


  
    

figure(1)
%tiledlayout(2,2)
%nexttile
plot(MT)
title('Dinero en el sistema (CUWUUUV)')
xlabel('Tiempo')
ylabel('Dinero en el sistema')

%nexttile
figure(2)
plot(A)
hold on
plot(B)
hold on
plot(X)
title('Evolución de clases')
xlabel('Meses')
ylabel('Clases')
legend('Workers','Employers','Unemployed')

%nexttile
figure(3)
plot(S)
ylabel("Entropía")
xlabel("Tiempo")
title("Evolución temporal de la entropía")

%nexttile
figure(4)
plot(objetivo)
grid on 
title("Función objetivo del sistema")
xlabel("Tiempo")
ylabel("O(n_1,n_2,...,n_c)")

figure(5)
histogram(n, 1000)
xlabel("Dinero")
ylabel("Peronas")
title("Dinero distribuido en clases")


C = 20;
D = linspace(1,M,C);
HC = hist(IMC(1,:),C);
HC = HC ./ sum(HC);
HW = hist(IMW(1,:),C);
HW = HW ./ sum(HW);

figure(6)
H = hist(IM(1,:),C);
H = H ./ sum(H);
plot(D,H,'-.')
set(gca,'YScale','log')
xlabel("Dinero")
ylabel("Distribucióon del dinero normalizada")
title("Distribuciones agregadas")
%legend('Año 1','Año 3', 'Año 5')
grid on

figure(7)
plot(D,HC,'-.')
set(gca,'YScale','log')
hold on
plot(D,HW,'-.')
set(gca,'YScale','log')
xlabel("Dinero")
ylabel("Distribucióon del dinero normalizada")
title("Distribuciones desagregadas")
legend('Caapitalistas','No capitalistas')
grid on
%set(gca,'XScale','log')
%set(gca,'YScale','log')

figure(9)
histogram(IMC(1,:),C);
xlabel("Dinero")
ylabel("Peronas")
title("Dinero distribuido en clase empresaria")

figure(10);
histogram(IMW(1,:),C);
xlabel("Dinero")
ylabel("Peronas")
title("Dinero distribuido en clase trabajadora y desempleada")



toc
%%%%%% Funciones %%%%%%
function [e,p,W] = H_1(W,a_1,e,n,w_avg,N,P,p) %Hiring rule
        if e(a_1) == 0 && W(a_1) == 0
            for k = 1:N
                if e(k) == 0
                    P(k) = n(k);                    
                else
                    P(k) = 0;
                end
            end
            for j = 1:length(P)
                    p(j) = P(j) / sum(P);
                    if p(j) < 0
                        p(j) = 0;
                    %elseif p(j) < 0
                    %    p(j) = 0;  
                    end
                    
            end
            %pipol = 1:1:N;
            %[~, c] = max(p);
            %c = int16(ceil(max(p)*randi(length(P),1,1))); 
            %c = int16(ceil(max(p)*N)); 
            %c = find(n==c,1);
            %c = randi(length(P),1,1);
            %c = int16(ceil(randsample(P,1)));
            c = randsample(N,1,true,p);
    
            if n(c) > w_avg 
                e(a_1) = c;
                W(c) = W(c) + 1;
            else
                e(a_1) = e(a_1);
            end
        end
end

function [V,n] = E_1(a_1,b_1,n,V,N)   %Expediture rule
    while b_1 == a_1
        b_1 = int16(ceil(rand*N));
    end
    M_b = [0,n(b_1)];
    m = (M_b(2)-M_b(1))*rand(1) + M_b(1); %Transacción aleatoria en el intervalo M_b
    V = V + m;
    n(b_1) = n(b_1) - m;
end

function [V,n] = M_1(V,n,W,a_1,e)   %Market sample rule
    if e(a_1) ~= 0 && W(a_1) == 0
         M_v = [0,V];                            
         m_v = (M_v(2)-M_v(1))*rand(1) + M_v(1);
        n(e(a_1)) = n(e(a_1)) + m_v;
        V = V - m_v;
    elseif e(a_1) == 0 && W(a_1) ~= 0
        M_v = [0,V];                            
        m_v = (M_v(2)-M_v(1))*rand(1) + M_v(1);
        n(a_1) = n(a_1) + m_v;
        V = V - m_v;
    end
end

function [e,W] = F_1(n,N,a_1,e,W,w_avg) % Firing rule
    if e(a_1) == 0 && W(a_1) ~= 0
        u = max(W(a_1) - (n(a_1)/w_avg),0);
        u = ceil(u);
        for j = 1:u
            for k = 1:N
                if e(k) == a_1
                    e(k) = 0;
                end
            end
        end
        W(a_1) = W(a_1) - u;
         
    end
end

function n = w_1(N,a_1,e,n,W,w)        %Wage payment
    if e(a_1) == 0 && W(a_1) ~= 0
        for k = 1:W(a_1)
            wage = (w(2)-w(1))*rand(1) + w(1);
            if n(a_1) >= wage
                for j = 1:N
                    if e(j) == a_1
                        n(j) = n(j) + wage;
                        n(a_1) = n(a_1) - wage;
                    end
                end
            else 
                wage = (n(a_1))*rand(1);
                for j = 1:N
                    if e(j) == a_1
                        n(j) = n(j) + wage;
                        n(a_1) = n(a_1) - wage;
                    end
                end
            end
        end
    end
end

function [workers,employers,unemployed] = countclass(e,W,N)
workers = 0;
employers = 0;
unemployed = 0;
for k = 1:N
    if e(k) ~= 0 && W(k) == 0
        workers = workers +1;
    elseif e(k) == 0 && W(k) ~= 0
        employers = employers + 1;
    elseif e(k) == 0 
        unemployed = unemployed + 1;
    end
end
end

function S = entropia(nk,N) %Entropia
            sumNk= 0;
            for k = 1:length(nk)
                if nk(k) ~= 0
                    sumNk = (nk(k) * log(nk(k))) + sumNk;
                end
        
            end
            S = N * log(N) - sumNk; 
end

function o_objetivo = o_objetivo(centers,C,nk,a,caso,sumO)
    
    if caso == 0
        for k = 1:C
            if nk(k) ~= 0
                sumO = (nk(k) * (a*centers(k))) + sumO;
            end
            o_objetivo = sumO;
            
        
        end
    else

        for k = 1:C
            if nk(k) ~= 0 
                sumO = (nk(k) * (1-exp(-a*centers(k)))) + sumO;
            end
            o_objetivo = sumO;
           
        
        end
    end
end
