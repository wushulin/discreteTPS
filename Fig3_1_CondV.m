clc;
clear;
global omega dt alpha nu theta
theta=2/2;
JJ=[8,16,32,64,128,256];
omega=0;
dt=1e-2;
alpha=1;
nuu=10.^(-9:0);
condV=zeros(length(JJ),length(nuu),2);

for j1=1:length(JJ)
    for j2=1:length(nuu)
        nu=nuu(j2);
        J=JJ(j1);
        C1=toeplitz([1;-1;zeros(J-2,1)],[1,zeros(1,J-1)]);
        C1(1,J)=-omega;
        C2=toeplitz([theta;1-theta;zeros(J-2,1)],[theta,zeros(1,J-1)]);
        C2(1,J)=(1-theta)*omega;
        C3=toeplitz([1-theta;zeros(J-1,1)],[1-theta,theta,zeros(1,J-2)]);
        C3(J,1)=theta*omega;
        Q=[inv(C2),zeros(J,J);zeros(J,J),inv(transpose(C2))]*[C1,-(dt/sqrt(nu))*C3;dt*alpha/sqrt(nu)*transpose(C3),transpose(C1)];
        [V,D]=eig(Q);
        condV(j1,j2,1)=cond(V,2);
        Q=[inv(C2),zeros(J,J);zeros(J,J),inv(transpose(C2))]*[C1,-(dt/nu)*C3;dt*alpha*transpose(C3),transpose(C1)];
        [V,D]=eig(Q);
        condV(j1,j2,2)=cond(V,2);
    end
end
j=1;
loglog(nuu,condV(6,:,j),'r-.o',nuu,condV(5,:,j),'b-.h',nuu,condV(4,:,j),'k-.s',...
    nuu,condV(3,:,j),'c-.+',nuu,condV(2,:,j),'y-.*',nuu,condV(1,:,j),'m-.d','MarkerSize',12,'LineWidth',1.5);
hold on;
j=2;
loglog(nuu,condV(6,:,j),'r-.o',nuu,condV(5,:,j),'b-.h',nuu,condV(4,:,j),'k-.s',...
    nuu,condV(3,:,j),'c-.+',nuu,condV(2,:,j),'y-.*',nuu,condV(1,:,j),'m-.d','MarkerSize',12,'LineWidth',1.5);
hold off;
set(gca,'fontsize',16)
shg
xlim([0.25e-9,3]);
set(gca,'xtick',[1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1]);
%set(gca,'ytick',[200,1000,3000]);
xlabel('$\nu$','Interpreter','latex','FontSize',23);
ylabel('condition number','FontSize',23);
if theta==0
    title('explicit Euler ($\theta=0$)','Interpreter','latex','FontSize',23);
elseif theta==0.5
    title('trapezoidal rule ($\theta=0.5$)','Interpreter','latex','FontSize',23);
else
    title('implicit  Euler ($\theta=1$)','Interpreter','latex','FontSize',23);
end
leg=legend('$J=256$','$J=128$','$J=64$','$J=32$','$J=16$','$J=8$');
set(leg,'fontsize',16,'interpreter','latex')
 %ylim([200,3500]);
% mesh(nuu,JJ,condV);
% set(gca,'fontsize',16)
%set(gca, 'XScale','log', 'YScale','log', 'ZScale','log');
% set(gca,'ytick',JJ);
% set(gca,'xtick',nuu);
%set(gca,'ztick',10.^(1:3));
%zlim([100,5e+3]);