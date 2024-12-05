  % simulation
clc, clear all, close all

n=2;
tstep = 0.01;
tau = 1;
pi=0.6;
tmin = -tau;
tmax = 2;
T=tmin:tstep:tmax;
L=length(T);
h=(-1)^0.5;j=(-1)^0.5;k=(-1)^0.5;
hj=k;jh=-k;jk=h;kj=-h;kh=j;hk=-j;

% initial condition
E1=[0.6;-1.5];
F1=[-1.0;2.2];
G1=[-0.9;2.5];
H1=[2;-1.2];

x = zerosq(n, L);
tempL=length(-tau:tstep:0);
for i=1:tempL
    x(:, i) = quaternion(E1,F1,G1,H1);
end
[M,N]=cd(x);
s=sym('s');
for oi=1:tmax  %1,2,3,4,5共5秒时间
    for ii=1:(round(1/tstep-1))  %细分，每秒的节点个数
        t=(oi-1)+ii*tstep;
        index=round((t+tau)/tstep)+1;  %加上延迟一共的时间(-0.5--5s)节点个数
        %Modified parameters        
        dM=-1.5*M(:,index-1)+A1*f1(M(:,index-1))+B1*f1(M(:,index-1-round(tau/tstep)))+C1*int(f1(M(:,index-1))*theta(index-1-s),s,index-1-round(pi/tstep),index-1);
        dN=-1.5*N(:,index-1)+A2*f2(N(:,index-1))+B2*f2(N(:,index-1-round(tau/tstep)))+C2*int(f2(N(:,index-1))*theta(index-1-s),s,index-1-round(pi/tstep),index-1);
        M(:,index)=M(:,index-1)+dM*tstep;
        N(:,index)=N(:,index-1)+dN*tstep;
    end
end

E=[-2.0;2.6];
F=[1.6;-2.1];
G=[2.2;2.2];
H=[-1.2;1.0];
y = zerosq(n, L);
tempL=length(-tau:tstep:0);
for i=1:tempL
    y(:, i) = quaternion(E,F,G,H);
end
[Y1,Y2]=cd(y);
% CD: returns two complex numbers which are the Cayley-Dickson components
% of the quaternion.
e = zerosq(n, L);
u = zerosq(n, L);
[U1,U2]=cd(u);
[E1,E2]=cd(e);

 for oi=1:tmax %index=tau:tstep:tmax
     for ii=1:(round(1/tstep-1))
         t=(oi-1)+ii*tstep;
         index=round((t+tau)/tstep)+1;
        
        E1(1,index-1)=Y1(1,index-1)-M(1,index-1);
        E1(2,index-1)=Y1(2,index-1)-M(2,index-1);
        E2(1,index-1)=Y2(1,index-1)-N(1,index-1);
        E2(2,index-1)=Y2(2,index-1)-N(2,index-1);
        %controller case1=33,35 and 32, 36; case2=22,25 and 20, 26. 
        if max(abs(E1(:,index-1)))>=1 && max(abs(E2(:,index-1)))>=1
           U1(:,index-1)=-22*E1(:,index-1)-20*sign(E1(:,index-1));
           U2(:,index-1)=-25*E2(:,index-1)-26*sign(E2(:,index-1));
        %   break;
        %end
        elseif max(abs(E1(:,index-1)))<1 && max(abs(E2(:,index-1)))<1
           % case1=31.3,and 35; case2=16.3,and 26. 
           U1(:,index-1)=-26.3*E1(:,index-1)-1.0*sign(E1(:,index-1))-1.5*(E1(:,index-1)).^0.5-1.4*(E1(:,index-1)).^1.5;
           U2(:,index-1)=-26*E2(:,index-1)-1.2*sign(E2(:,index-1))-1.6*(E2(:,index-1)).^0.5-1.5*(E2(:,index-1)).^1.5;
        %   break;
        end
        %drive system    
        dY1=-1.5*Y1(:,index-1)+A1*f1(Y1(:,index-1))+B1*f1(Y1(:,index-1-round(tau/tstep)))+C1*int(f1(Y1(:,index-1))*theta(index-1-s),s,index-1-round(pi/tstep),index-1)+U1(:,index-1);
        dY2=-1.5*Y2(:,index-1)+A2*f2(Y2(:,index-1))+B2*f2(Y2(:,index-1-round(tau/tstep)))+C2*int(f2(Y2(:,index-1))*theta(index-1-s),s,index-1-round(pi/tstep),index-1)+U2(:,index-1);
        Y1(:,index)=Y1(:,index-1)+dY1*tstep;
        Y2(:,index)=Y2(:,index-1)+dY2*tstep;
    end
end

figure(1)
grid on
h1=plot(T((tau/tstep+1):L), real(E1(1,(tau/tstep+1):L)),'r');
hold on
h2=plot(T((tau/tstep+1):L), imag(E1(1,(tau/tstep+1):L)),'b');
hold on
h3=plot(T((tau/tstep+1):L), real(E2(1,(tau/tstep+1):L)),'k');
hold on
h4=plot(T((tau/tstep+1):L), imag(E2(1,(tau/tstep+1):L)),'g');
hold on
h5=plot(T((tau/tstep+1):L), real(E1(2,(tau/tstep+1):L)),'Color',[0.5 0.5 0.5]);
hold on
h6=plot(T((tau/tstep+1):L), imag(E1(2,(tau/tstep+1):L)),'m');
hold on
h7=plot(T((tau/tstep+1):L), real(E2(2,(tau/tstep+1):L)),'c');
hold on
h8=plot(T((tau/tstep+1):L), imag(E2(2,(tau/tstep+1):L)),'Color',[0.67 0 1]);
hold off
xlabel('$t$', 'Interpreter','latex', 'FontSize', 12)
ylabel('state error')
lgd1=legend([h1,h2,h3,h4],'\epsilon^{0}_{1}(t)','\epsilon^{1}_{1}(t)','\epsilon^{2}_{1}(t)','\epsilon^{12}_{1}(t)','orientation','horizontal','location','north');
ah=axes('position',get(gca,'position'),'visible','off');
lgd2=legend(ah,[h5,h6,h7,h8],'\epsilon^{0}_{2}(t)','\epsilon^{1}_{2}(t)','\epsilon^{2}_{2}(t)','\epsilon^{12}_{2}(t)','orientation','horizontal','location','north');

% figure(2)
% plot(T((tau/tstep+1):L), real(M(1,(tau/tstep+1):L)),'g')
% hold on
% plot(T((tau/tstep+1):L), real(Y1(1,(tau/tstep+1):L)),'r')
% hold off
% xlabel('$t$', 'Interpreter','latex', 'FontSize', 10)
% legend({'x^{R}_{1}(t)','y^{R}_{1}(t)'},'FontSize',10);
% 
% figure(3)
% plot(T((tau/tstep+1):L), imag(M(1,(tau/tstep+1):L)),'g')
% hold on
% plot(T((tau/tstep+1):L), imag(Y1(1,(tau/tstep+1):L)),'r')
% hold off
% xlabel('$t$', 'Interpreter','latex', 'FontSize', 10)
% legend({'x^{I}_{1}(t)','y^{I}_{1}(t)'},'FontSize',10);
% 
% figure(4)
% plot(T((tau/tstep+1):L), real(N(1,(tau/tstep+1):L)),'g')
% hold on
% plot(T((tau/tstep+1):L), real(Y2(1,(tau/tstep+1):L)),'r')
% hold off
% xlabel('$t$', 'Interpreter','latex', 'FontSize', 10)
% legend({'x^{J}_{1}(t)','y^{J}_{1}(t)'},'FontSize',10);
% 
% figure(5)
% plot(T((tau/tstep+1):L), imag(N(1,(tau/tstep+1):L)),'g')
% hold on
% plot(T((tau/tstep+1):L), imag(Y2(1,(tau/tstep+1):L)),'r')
% hold off
% xlabel('$t$', 'Interpreter','latex', 'FontSize', 10)
% legend({'x^{K}_{1}(t)','y^{K}_{1}(t)'},'FontSize',10);
% 
% figure(6)
% plot(T((tau/tstep+1):L), real(M(2,(tau/tstep+1):L)),'m')
% hold on
% plot(T((tau/tstep+1):L), real(Y1(2,(tau/tstep+1):L)),'b')
% hold off
% xlabel('$t$', 'Interpreter','latex', 'FontSize', 10)
% legend({'x^{R}_{2}(t)','y^{R}_{2}(t)'},'FontSize',10);
% 
% figure(7)
% plot(T((tau/tstep+1):L), imag(M(2,(tau/tstep+1):L)),'m')
% hold on
% plot(T((tau/tstep+1):L), imag(Y1(2,(tau/tstep+1):L)),'b')
% hold off
% xlabel('$t$', 'Interpreter','latex', 'FontSize', 10)
% legend({'x^{I}_{2}(t)','y^{I}_{2}(t)'},'FontSize',10);
% 
% figure(8)
% plot(T((tau/tstep+1):L), real(N(2,(tau/tstep+1):L)),'m')
% hold on
% plot(T((tau/tstep+1):L), real(Y2(2,(tau/tstep+1):L)),'b')
% hold off
% xlabel('$t$', 'Interpreter','latex', 'FontSize', 10)
% legend({'x^{J}_{2}(t)','y^{J}_{2}(t)'},'FontSize',10);
% 
% figure(9)
% plot(T((tau/tstep+1):L), imag(N(2,(tau/tstep+1):L)),'m')
% hold on
% plot(T((tau/tstep+1):L), imag(Y2(2,(tau/tstep+1):L)),'b')
% hold off
% xlabel('$t$', 'Interpreter','latex', 'FontSize', 10)
% legend({'x^{K}_{2}(t)','y^{K}_{2}(t)'},'FontSize',10);

% figure(10)
% p1=plot(T((tau/tstep+1):L), real(M(1,(tau/tstep+1):L)),'r');
% hold on
% p2=plot(T((tau/tstep+1):L), imag(M(1,(tau/tstep+1):L)),'b');
% hold on
% p3=plot(T((tau/tstep+1):L), real(N(1,(tau/tstep+1):L)),'g');
% hold on
% p4=plot(T((tau/tstep+1):L), imag(N(1,(tau/tstep+1):L)),'m');
% hold on
% p5=plot(T((tau/tstep+1):L), real(Y1(1,(tau/tstep+1):L)),'c');
% hold on
% p6=plot(T((tau/tstep+1):L), imag(Y1(1,(tau/tstep+1):L)),'g');
% hold on
% p7=plot(T((tau/tstep+1):L), real(Y2(1,(tau/tstep+1):L)),'Color',[0.67 0 1]);
% hold on
% p8=plot(T((tau/tstep+1):L), imag(Y2(1,(tau/tstep+1):L)),'r');
% hold off
% xlabel('$t$', 'Interpreter','latex', 'FontSize', 10)
% lgd1=legend([p1,p2,p3,p4],'x^{R}_{1}(t)','x^{I}_{1}(t)','x^{J}_{1}(t)','x^{K}_{1}(t)','orientation','horizontal','location','north');
% ah=axes('position',get(gca,'position'),'visible','off');
% lgd2=legend(ah,[p5,p6,p7,p8],'y^{R}_{1}(t)','y^{I}_{1}(t)','y^{J}_{1}(t)','y^{K}_{1}(t)','orientation','horizontal','location','north');
% 
% figure(11)
% l1=plot(T((tau/tstep+1):L), real(M(2,(tau/tstep+1):L)),'r');
% hold on
% l2=plot(T((tau/tstep+1):L), imag(M(2,(tau/tstep+1):L)),'b')
% hold on
% l3=plot(T((tau/tstep+1):L), real(N(2,(tau/tstep+1):L)),'g');
% hold on
% l4=plot(T((tau/tstep+1):L), imag(N(2,(tau/tstep+1):L)),'m');
% hold on
% l5=plot(T((tau/tstep+1):L), real(Y1(2,(tau/tstep+1):L)),'c');
% hold on
% l6=plot(T((tau/tstep+1):L), imag(Y1(2,(tau/tstep+1):L)),'g');
% hold on
% l7=plot(T((tau/tstep+1):L), real(Y2(2,(tau/tstep+1):L)),'Color',[0.67 0 1]);
% hold on
% l8=plot(T((tau/tstep+1):L), imag(Y2(2,(tau/tstep+1):L)),'r');
% hold off
% xlabel('$t$', 'Interpreter','latex', 'FontSize', 12)
% lgd1=legend([l1,l2,l3,l4],'x^{R}_{2}(t)','x^{I}_{2}(t)','x^{J}_{2}(t)','x^{K}_{2}(t)','orientation','horizontal','location','north');
% ah=axes('position',get(gca,'position'),'visible','off');
% lgd2=legend(ah,[l5,l6,l7,l8],'y^{R}_{2}(t)','y^{I}_{2}(t)','y^{J}_{2}(t)','y^{K}_{2}(t)','orientation','horizontal','location','north');

