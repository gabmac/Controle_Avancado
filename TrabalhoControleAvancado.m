    close all
    clear 
    clc
    
    %variaveis que serao usadas ao longo do projeto
    syms J m l M g u k c1 c2 c3 c4 k
    % x_n = c_n - cequilibrio_n
    syms x1 x2 x3 x4
    
    %% Matrizes na forma csi_ponto*M1 = M2
    % matrizes de entrada para o sistema nao linear

    M1 = [(J+m*l^2) m*l*cos(c1);m*l*cos(c1) (M+m) ]; %M1*csi_ponto = M2
    M2 = [u - m*g*l*sin(c1); -k*c3+m*l*c2^2*sin(c1)];

    equations = simplify(M1\M2);%equacoes de csi ponto
    csi = [c1 c2 c3 c4]; %vetor com o espaco de estado pra derivacao
    
    %% Parte 1&2 - Representacao de Estado nao linear
    funcoes = [c2 equations(1) c4 equations(2)]; %espaço de Estado csi_ponto = M*csi
    
    %% Parte 3 - Pontos Equilibrio para entrada de controle nula
    
    % encontra os pontos de equilibrio, sendo sol a solucao para tal
    ptoequi = funcoes';
    ptoequi = ptoequi == 0; %iguala as funcoes a 0
    ptoequi = subs(ptoequi,u,0); %faz a entrada de controle nula
    sol = solve(ptoequi,csi); %encontra os pontos de equilibrio
    
    %% parte 4&5 - Encontrando as funcoes linearizadas em torno do ponto [0,0,0,0] 
    %             
    % c1 c2 c3 c4 sao as variaveis de espaço de estado, referentes aos csi's
    
    equilibrio_k1 = [0 0 0 0]; %linearizacao em volta dos pontos pedidos
    
    X = [x1 x2 x3 x4];% vetor equivalente a x_n = (c_n - c_n_equilibrio)
    
    
    %lineariza para os valores da matriz A
    multA = sym(zeros(4));
    multB = sym(zeros(1,4));
    for i = 1:1:4
        for j = 1:1:4
            derA = diff(funcoes(i),csi(j)); %deriva cada funcao para cada variavel
            substA = subs(derA,csi,equilibrio_k1); %substitui no ponto de equilibrio
            multA(i,j) = substA; %monta a matriz A literal
        end
    end

    A = multA;%Matriz A
    
    %lineariza para os valores da matriz B
    for i = 1:1:4
        derB = diff(funcoes(i),u); %deriva cada funcao para cada variavel
        substB = subs(derB,csi,equilibrio_k1); %substitui no ponto de equilibrio
        multB(1,i) = substB; % monta a matriz B literal
    end
    
    B = multB'; %MatrizB

    %% parte 6 - Analise da estabilidade de A 
    
    % dados fisicos do problema
    prop_fisicas_literal = [M m k J l g];
    prop_fisicas_numericas = [4 2 1 1 0.3 9.8];
    
    A = subs(A,prop_fisicas_literal, prop_fisicas_numericas);%substituicoes das propriedades
    B = subs(B,prop_fisicas_literal, prop_fisicas_numericas);%substitui as propriedads
    
    [Avetor,Avalor] = eig(A);
    Avalor = simplify(vpa(Avalor)); % autovalor de A
    
    
    %% parte 7 - Calculo numerico das respostas no tempo
   
    sim('trabControle17a.slx',[0,150]);
    % condicoes iniciais
    Xo = [pi/3 0 0.1 0]'; 
    
    % Plot das variaveis angulares x1, x2 (theta, theta ponto) 
    % e x3, x4 (posicaoo, velocidade).
    
    figure(1);
    title('u(t) = 0');
    hold on
    plot(tout,x1_t);
    plot(tout,x2_t);
    legend('Posição Angular','Velocidade Angular');
    saveas(gcf,'x1x2.png')
    hold off;
    figure(2);
    hold on
    title('u(t) = 0');
    plot(tout,x3_t);
    plot(tout,x4_t);
    legend('Posição Linear','Velocidade Linear');
    saveas(gcf,'x3x4.png')
    hold off;
    
    %% parte 8 - encontrando a nova matriz A
    
    %para u = - Kx => xdot = (A-BK)x // newA = (A-BK) xdot = (newA)x
    
    K = [8.75 31.50 -14.50 66.50]; % valores dados de ganho da matriz K
    newA = eval(A-B*K); 
    
    % Plot das variaveis angulares nx1, nx2 (theta, theta ponto)
    % e nx3, nx4 (posicao, velocidade).
    
    figure(3);
    hold on
    title('u(t) = -Kx')
    plot(tout,nx1_t);
    plot(tout,nx2_t);
    legend('Posição Angular','Velocidade Angular');
    saveas(gcf,'nx1nx2.png')
    hold off;
    figure(4);
    hold on
    title('u(t) = -Kx');
    plot(tout,nx3_t);
    plot(tout,nx4_t);
    legend('Posição Linear','Velocidade Linear');
    saveas(gcf,'nx3nx4.png')
    hold off;
    
    %% parte 9 - avaliando a estabilidade de newA e encontrando P com laypunov
    
    % encontrando os autovalores e autovetores de newA
    [newAvetor,newAvalor] = eig(newA);
    newAvalor = simplify(vpa(newAvalor)); %autovalor de newA
    
    % usando lyap para encontrar P
    Q = eye(4);
    P = lyap(newA,Q);
    [Pvetor,Pvalor] = eig(P);
    Pvalor = simplify(vpa(Pvalor));
    
    %% parte 10 - Matriz de controlabilidade e forma canonica controlavel
    
    % a matriz de controlabilidade Cfresco = [ B AB A^2B A^3B]
    Cfresco = [ B A*B A^2*B A^3*B];
    
    % Calculando a matriz na forma canonica controlavel
    C = [0 0 1 0]; %definindo y=xc
    B = zeros(4,1);
    [Num,Div] = ss2tf(newA,B,C,0); % encontrando a funcao transferencia e 
    [newA,B,C,D] = tf2ss(Num,Div); % voltando para para a forma de espaco de estado
    
    % arrumando as formas das matrizes newA, B e C 
    newAint = zeros(4);
    Bint = zeros(4,1);
    Cint = zeros(1,4);
    for i=1:1:4
        for j=1:1:4
            newAint(j,i) = newA(5-j,i);
        end
            Bint(i,1) = B(5-i,1);
            Cint(1,i) = C(1,5-i);
    end
    for i=1:1:4
        for j=1:1:4
            newA(i,j) = newAint(i,5-j); %newA esta na forma canonica controlavel adequada
        end
    end
    B = Bint; % B esta na forma canonica controlavel adequada
    C = Cint; % C esta na forma canonica controlavel adequada
    
    