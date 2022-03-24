%--------------------------------------------------------------------------
% Soluçao da equação de Petroleos, usando Matlab    
% Pontificia Universidade Catolica do Rio de Janeiro  pdepe
% Data: 11-03-2018
% Jaime Andres Castañeda Barbosa 
% Finite Difference Methods (FDM)
%--------------------------------------------------------------------------


N =  input('digite o numero de nós  N = :');
disp('');
dt = input('digite o numero de intervalos de tempo dt = :');
disp('');
Tmax = input('digite o tempo maximo da simulacao Tmax = :');
disp('');

% dados do problema

Rw = 0.5;
Rinf = 100;
K = 5e-10;
Pini = 1e+6;
Qw = 1e-3;
beta = 1e-8;
mu = 1;
phi = 0.25;

D = K/(mu*beta*phi);

L = Rinf - Rw;

IT = Tmax/dt;   % numero de passos de tempo
dr = L/(N-1);   % distancia entre os nos

r(1) = Rw;
for i=2:N
    r(i) = r(i-1)+dr;
end

% condiciones iniciais para t=0
t(1) = 0;
for i=1:N;
    p(i,1) = Pini;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%

for j=2:IT
    
    t(j) = t(j-1)+dt;
    
    % condiciones de contorno
    
    A(1,1) = -1/dr;
    A(1,2) = 1/dr;
    f(1,1) = mu*Qw/K;
    
    A(N, N-1)=-1;
    A(N,N)=1;
    f(N,1)=0;
    
    
    for i=2:N-1
        A(i,i+1)=-D/(2*r(i)*dr^2)*(r(i)+r(i+1));
        A(i,i)=1/dt+D/(2*r(i)*dr^2)*(2*r(i)+r(i+1)+r(i-1));
        A(i,i-1)=-D/(2*r(i)*dr^2)*(r(i)+r(i-1));
        f(i,1)=p(i,j-1)/dt;
    end
    
    Patual = inv(A)*f;
    
    for i=1:N
        p(i,j) = Patual(i);
    end;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(r,p(:,IT));
title(' Pressao final ao longo do reservatorio');
xlabel(' profundidade' );
ylabel(' pressao' );

%figure
%a=round((N-1)/2);
%plot(t,p(1,:),t,p(a+1,:),t,p(N,:))
%title(' variacao da pressao em determinados pontos' );
%%figure(2);
%%xlabel(' t' );
%%ylabel(' p' );
%%legend(' poco','L/2','L');

    








