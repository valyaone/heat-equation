clear all;
clc;

g0= @ (x,y) 0; % начальное условие
f= @ (x,y,t) x*t^2*sin(y);

a11=1;

a_x=0; % левая координата по оси x
b_x=1; % правая координата по оси x
N_x=100; % количество шагов по оси x
h_x=(b_x-a_x)/N_x; % шаг по координате x

a_y=0; % левая координата по оси y
b_y=pi/2; % правая координата по оси y
N_y=200; % количество шагов по оси y
h_y=(b_y-a_y)/N_y; % шаг по координате y

t0=0;% начальный момент времени
tcon=2; % до какого момента времени смотрим решение
M=100; % количество шагов по оси времени
tau=(tcon-t0)/M; % шаг по оси времени

%сетка по х
X=zeros(1,N_x+1);
for i=1:N_x+1
    X(i)=h_x*(i-1)+a_x; % сетка по оси X
end  

%сетка по у
Y=zeros(1,N_y+1);
for i=1:N_y+1
    Y(i)=h_y*(i-1)+a_y; % сетка по оси Y
end

%сетка по Т
T=zeros(1,M+1);
for i=1:M+1
    T(i)=tau*(i-1)+t0; % сетка по оси Y
end

U=zeros(N_x+1,N_y+1,M+1); % массив для численного решения
for i=1:N_x+1
    for k=1:N_y+1
        W(k,i,1)=g0(X(i),Y(k)); % начальные условия
    end
end

%Массив перехода со слоя j на j+1/2. 

c=zeros(1,N_y+1);
c(1)=1;
for i=2:N_y
    c(i)=1+2*a11*tau/(h_y)^2; % матрица значение ci
end

a=zeros(1,N_y+1);
for i=2:N_y
     a(i)=a11*tau/(h_y)^2; % матрица значение ai
end

b=zeros(1,N_y);
b(1)=0;
for i=2:N_y
     b(i)=a11*tau/(h_y)^2; % матрица значение bi
end

a(N_y+1)=4/3-c(N_y)/(3*a(N_y));
c(N_y+1)=1-b(N_y)/(3*a(N_y));

% второй переход (c j+1/2 на j+1) слой

a_2=zeros(1,N_x+1);
for i=2:N_x
    a_2(i)=a11*tau/(h_x)^2; % матрица значение ai
end

c_2=zeros(1,N_x+1);
c_2(1)=1;
for i=2:N_x
    c_2(i)=1+2*a11*tau/(h_x)^2; % матрица значение ci
end

b_2=zeros(1,N_x);
b_2(1)=0;
for i=2:N_x
    b_2(i)=a11*tau/(h_x)^2; % матрица значение bi
end

c_2(1)=3-a_2(2)/(b_2(2));
b_2(1)=4-c_2(2)/(b_2(2));
a_2(N_x+1)=4/3-c_2(N_x)/(3*a_2(N_x));
c_2(N_x+1)=1-b_2(N_x)/(3*a_2(N_x));
        
for j=1:M
    for n=2:N_x
        % первый переход (вспомогательный)
        F=zeros(1,N_y+1);
        F(1)=0;
        for i=2:N_y
            F(i)=(W(i,n+1,j)-2*W(i,n,j)+W(i,n-1,j))/(h_x)^2+(W(i+1,n,j)-2*W(i,n,j)+W(i-1,n,j))/(h_y)^2+f(X(n),Y(i),T(j)+tau/2); % матрица значение fi
        end
        F(N_y+1)=F(N_y)/(3*a(N_y));
        PROB=progonka(a,b,c,F);
        for i=1:N_y+1
            W_1(i,n)=PROB(i);% промежуточный Ошибка
        end
    end
    for m=2:N_y
        % второй переход (на j+1) слой
        F_2=zeros(1,N_x+1);
        for i=2:N_x
            F_2(i)=W(m,i,j)+tau*W_1(m,i)-a11*tau/(h_x)^2*(W(m,i+1,j)-2*W(m,i,j)+W(m,i-1,j)); % матрица значение fi
        end
        F_2(1)=F_2(2)/(b_2(2));
        F_2(N_x+1)=F_2(N_x)/(3*a_2(N_x));
        PROB=progonka(a_2,b_2,c_2,F_2);
        for i=1:N_x+1
            W(m,i,j+1)=PROB(i);
        end
    end
    for i=1:N_x+1
        W(1,i,j+1)=0;
        W(N_y+1,i,j+1)=4/3*W(N_y,i,j+1)-1/3*W(N_y-1,i,j+1);
    end
end
figure(2);
for j = 1:(M+1)
mesh(X,Y,W(:,:,j));
hold on; 
xlim([a_x b_x]);
ylim([a_y b_y]);
zlim([0 1]);
title('Численное решение');
xlabel ('x'); ylabel ('y' ); zlabel('U')
txt2 = num2str(j); txt1 = 'Итерация = ';
txt = [txt1 ' ' txt2] ;
text(-0.4,0.85,txt);
hold off; 
drawnow; 
pause (0.01);
end
                