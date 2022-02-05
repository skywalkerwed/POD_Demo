clear;clc;close all;
a=1;

l=0.2;
x=linspace(0, 0.2,2e2);
dx=x(2)-x(1);
t=0:0.0001:0.02;

u0=20;
u1=40;
u2=100;

N=500;
[X, T]=meshgrid(x, t);
u=u1*ones(length(t), length(x))+(u2-u1)/l*X.*ones(length(t), length(x));

for n=1:500
    u=u+2/pi*1/n*((u0-u1)*(1-(-1)^n)+(u2-u1)*(-1)^n)*sin(n*pi*X/l).*exp(-n^2*pi^2*a^2*T/l^2);
end
%% demo
% figure()
% i=1;
% plot(x, u(i, :));
% figure()
% i=2;
% plot(x, u(i, :));
% for i=1:length(t)
%     clf;
%     plot(x, u(i, :));
%     pause(0.05);
% end
%% numeric
usol0 = u0*ones(1, length(x));
usol0(1) = u1;
usol0(end) = u2; 
A = zeros(1, length(x));
for i = 2:length(x)-1
    A(i, i-1) = 1;
    A(i, i) = -2;
    A(i, i+1)=1;
end
A(end+1, :) = zeros(1, length(x));
[t, usol] = ode45('myTranseq', t, usol0,[], a, dx, A);

%% 1-d animation
err=abs(usol-u)./abs(u);
figure()
for i=1:length(t)
    clf;
    subplot(2,2,1)
    plot(x, u(i, :));
    subplot(2,2,2)
    plot(x, usol(i, :));
    subplot(2,2,[3 4])
    plot(x, err(i, :));
    pause(0.05)
end
%% POD
[U,S,V] = svd(u');
Eng = diag(S).^2/norm(diag(S),2)^2;
N = find(Eng < 1e-6);N=N(1)-1;%cut-of engery lower than 0.0001%
phi = U(:, 1:N);
a0 = (phi'*usol0')';
A_r = phi'*A*phi;
[t, asol] = ode45('myTranseq', t, a0, [], a, dx, A_r);
u_r = (phi*asol')';
%% 1-d animation
err = abs(u_r-u)./u;
figure()
pic_num=1;
for i=1:5:length(t)
    clf;
    subplot(2,2,1)
    plot(x, u(i, :), 'linewidth', [2]);
    title('exact')
    subplot(2,2,2)
    plot(x, u_r(i, :), 'linewidth', [2]);
    title('POD')
    subplot(2,2,[3 4])
    plot(x, err(i, :), 'linewidth', [2]);
    title('error')
    
    F=getframe(gcf);
    I=frame2im(F);
    [I,map]=rgb2ind(I,256);
    
    if (pic_num == 1)
    imwrite(I,map,'1-d.gif','gif','Loopcount',inf,'DelayTime',1e-4);
    else
    imwrite(I,map,'1-d.gif','gif','WriteMode','append','DelayTime',1e-4);
    end
    pic_num = pic_num + 1;
%     pause()
end
%% 2-d animation
[X, Y] = meshgrid(x, x);
figure()
set(gcf,'position',[500, 500, 720,240])
pic_num = 1;
for i=1:10:length(t)
    clf;
    u2d = repmat(u(i, :), length(x), 1);
    u_r2d = repmat(u_r(i, :), length(x), 1);
    subplot(1, 2, 1)
    pcolor(X, Y, u2d),shading interp,colormap('hot')
    title('exact')
    subplot(1, 2, 2)
    pcolor(X, Y, u2d),shading interp,colormap('hot')
    title('POD')
    %suptitle(strcat('T=',num2str(t(i))))
    
    F=getframe(gcf);
    I=frame2im(F);
    [I,map]=rgb2ind(I,256);
    
    if (pic_num == 1)
    imwrite(I,map,'2-d.gif','gif','Loopcount',inf,'DelayTime',1e-4);
    else
    imwrite(I,map,'2-d.gif','gif','WriteMode','append','DelayTime',1e-4);
    end
    pic_num = pic_num + 1;
    
%     pause()
end
    