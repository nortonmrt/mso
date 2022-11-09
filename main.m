%% Exemplo de programa para análise de treliça 2D
% adaptado em 26/04/2021 por Norton Trennepohl
clear;

%% Dados de entrada
E = [2e11;70e9]; %em Pa
C = E;

L = [1;1]; %em m 

A = [0.1;0.2]; %em m²

%% Construção da malha
%coordenadas dos nós
x1 = [0 0];
x2 = [L(1)*cos(pi/4) L(1)*sin(pi/4)];
x3 = [0,(L(1)+L(2))*sin(pi/4)];

xy = [x1; x2; x3];

%número de nós e número de elementos
number_nodes = size(xy,1);
number_elements = number_nodes - 1;

%matriz de conectividade
eNodes = [1 2;2 3];

%% Matriz de rigidez global
nDof = 2*number_nodes;
K = zeros(nDof);

for i=1%:number_elements
    eNode = eNodes(i,:);
    eDof = [eNode;eNode + number_nodes];
    %obs: lembrar que ele convenciona os GDLs 1 2 3 para a horizontal e 4 5
    %6 para a vertical, por isso o algoritmo acima funciona.
    nnDof = length(eNode);
    
    %Comprimento do elemento
    L = sqrt((xy(eNode(2),2)-xy(eNode(1),1))^2 + (xy(eNode(2),2)-xy(eNode(1),2))^2);
    cos = (xy(eNode(2),1)-xy(eNode(1),1))/L;
    sin = (xy(eNode(2),2)-xy(eNode(1),2))/L;

    %Quadratura de Gauss
    detJ = L/2;
    invJ = 1/detJ;
    
    [shape,nDeriv] = shapeFunct_Truss(0);
    gaussWt = 2;
    
    %Matriz B
    Xderiv = nDeriv*invJ;
    B = zeros(1,nnDof); B(1:nnDof) = Xderiv(:);
    
    %Matriz de rigidez local
    Ke = A(i,1)*(B'*C(i,1)*B)*detJ*gaussWt;
    
    %%%%% Transformation matrix
    Te = [cos, sin, 0,   0;
           0,   0, cos, sin];
    KeG = Te'*Ke*Te;
    
    K(eDof,eDof) = K(eDof,eDof) + KeG;
end

%% Condições de contorno
fixedNodes = [1;3];

fixedDof = [fixedNodes; fixedNodes + number_nodes];

%% Carregamentos
force = zeros(nDof,1);

force(2) = 5e9; force(5) = 3e9;

%% Solução
displacements = solution(nDof,fixedDof,K,force);

%% Gráfico
scale = 1; %pode ser alterada posteriormente

u_Dof = (1:number_nodes)';
v_Dof = u_Dof + number_nodes;
u = displacements(u_Dof,1); %desolcamentos em x
v = displacements(v_Dof,1); %deslocamentos em y

figure
original_config = xy;
deformed_config = xy + scale*[u,v];
plot(original_config(:,1),original_config(:,2));
hold on
plot(deformed_config(:,1),deformed_config(:,2));
xlabel('x (m)'); %legenda em x
ylabel('displacement: u (m)'); %legenda em y

%% Tensão e deformação
stress = zeros(number_elements,1);
strain = zeros(number_elements,1);

%strain: e = B*displacements
%stress: sigma = C*e

for i=1:number_elements
    eNode = eNodes(i,:);
    eDof = [eNode;eNode + number_nodes];
    %obs: lembrar que ele convenciona os GDLs 1 2 3 para a horizontal e 4 5
    %6 para a vertical, por isso o algoritmo acima funciona.
    nnDof = length(eNode);
    
    %Comprimento do elemento
    L = sqrt((xy(eNode(2),2)-xy(eNode(1),1))^2 + (xy(eNode(2),2)-xy(eNode(1),2))^2);
    cos = (xy(eNode(2),1)-xy(eNode(1),1))/L;
    sin = (xy(eNode(2),2)-xy(eNode(1),2))/L;

    %Quadratura de Gauss
    detJ = L/2;
    invJ = 1/detJ;
    
    [shape,nDeriv] = shapeFunct_Truss(0);
    gaussWt = 2;
    
    %Matriz B
    Xderiv = nDeriv*invJ;
    B = zeros(1,nnDof); B(1:nnDof) = Xderiv(:);
    
    Te = [cos, sin, 0,   0;
           0,   0, cos, sin];
    
    eDof_v = [eDof(:,1);eDof(:,2)];
    local_disp = Te * displacements(eDof_v,1);
    local_strain = B*local_disp;
    
    %deformação
    strain(i,:) = local_strain;
    
    %tensão
    Ci = C(i,1);
    local_stress = Ci * local_strain;
    
    stress(i,:) = local_stress;
    
    disp("Deformações:");
    disp(strain);
    disp("Tensões:");
    disp(stress);
end

%% Fim
