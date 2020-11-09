%  close all
%  clear all
%  clc
%  set(groot,'defaultFigureVisible','on')
%  cd 'C:\DUKE\courses\FEM\project'
%% set parameters
% resonance_case=0; %set 0 for normal case
% if resonance_case
%     n_resonance=1;%resonance eigenmode in x
%     m_resonance=1;%resonance eigenmode in y
%     omega=sqrt(n_resonance*pi^2+m_resonance*pi^2); %frequency for eigen value lambda_nm
%     oscillation='resonance case';
% else
%     omega=1;
%     oscillation='normal case';
% end
% m=20;%mesh size
% dt=0.01;%time step
% Ndt=500;%maximum iteration time
% dt_save=10;%time step save interval
% savecount=0;
% t_iteration=0;
%% generate mesh
num_e=4*m^2;
num_n=(m+1)^2+m^2;
%node indices
nodes_all=1:num_n;
nodes_boundary=sort([1:m+1, 2*m+2:2*m+1:2*m*m-m,3*m+2:2*m+1:2*m*m,2*m*m+m+1:2*m*m+2*m+1]);
nodes_inner=setdiff(nodes_all,nodes_boundary);
for i=1:m^2
    %get nodes label
    row=ceil((i-0.5)/m);
    column=mod(i-1,m)+1;
    TL=column+(2*m+1)*row;
    TR=column+(2*m+1)*row+1;
    BL=column+(2*m+1)*(row-1);
    BR=column+(2*m+1)*(row-1)+1;
    C=column+(2*m+1)*(row-1)+m+1;
    %get nodes coordinates
    node_coord(:,BL)=[column-1;row-1]/m;
    node_coord(:,BR)=[column;row-1]/m;
    node_coord(:,TL)=[column-1;row]/m;
    node_coord(:,TR)=[column;row]/m;
    node_coord(:,C)=[column-0.5;row-0.5]/m;
    %get element connectivity
    ele_conn(:,4*i-3)=[BR;C;BL];
    ele_conn(:,4*i-2)=[TR;C;BR];
    ele_conn(:,4*i-1)=[TL;C;TR];
    ele_conn(:,4*i)=[BL;C;TL]; 
end
%% compute and visualize shape functions
for i=1:num_e
    syms x y;
    node_ele=node_coord(:,ele_conn(:,i));
    Se=polyarea(node_ele(1,:),node_ele(2,:));
    for j=1:3
    a(i,j)=(node_ele(1,mod(j,3)+1)*node_ele(2,mod(j+1,3)+1)-node_ele(1,mod(j+1,3)+1)*node_ele(2,mod(j,3)+1))/2/Se;
    b(i,j)=(node_ele(2,mod(j,3)+1)-node_ele(2,mod(j+1,3)+1))/2/Se;
    c(i,j)=(node_ele(1,mod(j+1,3)+1)-node_ele(1,mod(j,3)+1))/2/Se;
    phi(i,j)=a(i,j)+b(i,j).*x+c(i,j).*y;
    end
end
%% compute local K M b0
K=zeros(num_n,num_n);
M=zeros(num_n,num_n);
b0=zeros(num_n,1);
Me=zeros(4,3,3);
Ke=zeros(4,3,3);
b0e=zeros(3,1);
for i=1:4
    Ne=phi(i,:);
    Be=diff(Ne,x)+diff(Ne,y);
    NeTNe=transpose(Ne)*Ne;
    BeTBe=transpose(Be)*Be;
    NeT=transpose(Ne);
    node_ele=node_coord(:,ele_conn(:,i));
    triangleX = node_coord(1,ele_conn(:,i));
    triangleY = node_coord(2,ele_conn(:,i));
    for j=1:3
        for k=1:3
            fMe = matlabFunction(NeTNe(j,k),'vars',[x y]);
            fBe = matlabFunction(BeTBe(j,k),'vars',[x y]);
            Me(i,j,k)=integral2( @(x,y) inpolygon(x, y, triangleX, triangleY) .* fMe(x,y), min(triangleX), max(triangleX), min(triangleY), max(triangleY),'Method', 'iterated');
            Ke(i,j,k)=integral2( @(x,y) inpolygon(x, y, triangleX, triangleY) .* fBe(x,y), min(triangleX), max(triangleX), min(triangleY), max(triangleY),'Method', 'iterated');
        end
        fNe = matlabFunction(Ne(j),'vars',[x y]);
        b0e(j,i)=integral2( @(x,y) inpolygon(x, y, triangleX, triangleY) .* fNe(x,y), min(triangleX), max(triangleX), min(triangleY), max(triangleY),'Method', 'iterated');
    end
end
Me=vpa(Me,5);
Ke=vpa(Ke,5);
b0e=vpa(b0e,5);
%% compute global K M b0
for i =1:num_e
    M(ele_conn(:,i),ele_conn(:,i))=reshape(Me(mod(i-1,4)+1,:,:),3,3)+M(ele_conn(:,i),ele_conn(:,i));
    K(ele_conn(:,i),ele_conn(:,i))=reshape(Ke(mod(i-1,4)+1,:,:),3,3)+K(ele_conn(:,i),ele_conn(:,i));
    b0(ele_conn(:,i))=b0e(:,mod(i-1,4)+1)+b0(ele_conn(:,i));
end
%% initial conditions Dirchlet 0 Neumann 0
uh_new=zeros(num_n,1);
uh_old=zeros(num_n,1);
uh_old_old=uh_old;
t_iteration=t_iteration+1;
% save initial condition
savecount=savecount+1;
uh_storage(:,savecount)=uh_old_old;
%% solve equation by time iteration
while(t_iteration<Ndt)
    t_iteration=t_iteration+1;
    t=t_iteration*dt;
    uh_new(nodes_inner)=inv(M(nodes_inner,nodes_inner)/dt^2+K(nodes_inner,nodes_inner))...
        *(2*M(nodes_inner,nodes_inner)/dt^2*uh_old(nodes_inner)-...
        M(nodes_inner,nodes_inner)/dt^2*uh_old_old(nodes_inner)+b0(nodes_inner)*cos(omega*t));
    
    if mod(t_iteration,dt_save)==0
        savecount=savecount+1;
        uh_storage(:,savecount)=uh_new;
    end
    uh_old_old=uh_old;
    uh_old=uh_new;
end
%% post processing
% t_storage=linspace(0,t,savecount);
% for plttime=[2 5]
%     i=find(t_storage==plttime);
%     figure
%     [xi,yi] = meshgrid(0:1/m/10:1, 0:1/m/10:1);
%     zi = griddata(node_coord(1,:),node_coord(2,:),uh_storage(:,i),xi,yi);
%     s=surf(xi,yi,zi,'CDataMapping','scaled');
%     set(gca,'xlim',[0 1],'ylim',[0 1]);
%     grid on
%     plttitle={'simulated vertical velocity surface plot', ['t=',num2str(plttime),'s  ',oscillation]};
%     title(plttitle)
%     xlabel('x')
%     ylabel('y')
%     zlabel('temperature')
%     %shading interp
%     s.EdgeColor = 'none';
%     hold on
%     spacing = 10;  % play around so it fits the size of your data set
%     for i = 1 : spacing : length(xi(:,1))
%         plot3(xi(:,i), yi(:,i), zi(:,i),'-g');
%         plot3(yi(:,i), xi(:,i), zi(:,i),'-g');
%     end
%     colorbar
%     saveas(gcf,['simulated vertical velocity surface plot at t=',num2str(plttime),'s ',oscillation,'.png'])
% 
%     figure
%     x = [0 1];
%     y = [1 0];
%     image(x,y,zi,'CDataMapping','scaled')
%     set(gca,'xlim',[0 1],'ylim',[0 1]);
%     ax = gca;
%     ax.YDir = 'normal';
%     grid on
%     plttitle={'simulated vertical velocity heatmap', ['t=',num2str(plttime),'s  ',oscillation]};
%     title(plttitle)
%     xlabel('x')
%     ylabel('y')
%     colorbar
%     saveas(gcf,['simulated vertical velocity heatmap at t=',num2str(plttime),'s ',oscillation,'.png'])
% end
% 
% figure
% plot(t_storage,max(abs(uh_storage)),'-o')
% plttitle=['simulated oscillation magnitude over time in ',oscillation];
% title(plttitle)
% xlabel('time/s')
% ylabel('oscillation magnitude/ms^{(-1)}')
% grid on
% saveas(gcf,[plttitle,'.png'])
