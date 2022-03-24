function parabolic1
%--------------------------------------------------------------------------
% Soluçao da equação de Petroleos, usando Matlab    
% Pontificia Universidade Catolica do Rio de Janeiro  pdepe
% Jaime Andres Castañeda Barbosa 
% Feito pelo metodo pdepe do matlab
%--------------------------------------------------------------------------

rw = 0.5; % m
Rinf = 100; % m
k = 5e-10; % condutividade hidraulica m^2
Pi = 10e6;   % Pressao inicial do poço (Pa)
Qw = 10e-3;   % vazao m^3/s
n = 0.25;     % porosidade do reservatorio
ni = 1;    % viscosidade (Pa*s)
Beta = 1e-8;  % compressibilidade do meio poroso 1/
Time = 1000  % tempo em (seg)

% EDP  FUNCTION PARABOLIC 


m = 1;
x = linspace(rw,Rinf,100);
t = linspace(0,Time,50);

options=odeset('RelTol',1e-4,'AbsTol',1e-4,'NormControl','off','InitialStep',1e-7)

u = pdepe(m,@prespde,@presic,@presbc,x,t,options,k,Pi,Qw,n,ni,Beta);
Pressao = u(:,:,1);


plot (x,Pressao(end,:))
xlabel('Profundidade [L]');
ylabel('Pressao [kPa]');





