close
clear
clc
set(groot,'defaultFigureVisible','on')
cd 'C:\DUKE\courses\FEM\project'
%% set parameters
m=5;%mesh size
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
    if i==1
        [X,Y] = meshgrid(0:0.01/m:1/m,0:0.01/m:1/2/m);
        figure
        for j=1:3
            Z=double(subs(phi(i,j), {x,y}, {X, Y}));
            Z(X+Y>1/m)=nan;
            Z(X<Y)=nan;
            surf(X,Y,Z)
            hold on
        end
        set(gca,'xlim',[0 1/m],'ylim',[0 1/2/m]);
        grid on
        title('elementwise shape functions')
        xlabel('x')
        ylabel('y')
        zlabel('value')
        saveas(gcf,'elementwise shape functions.png')
    end
end
%% patch test with heat conduction problem
tem_analytical_fun=@(x,y)-pi*exp(1)+pi*x+exp(1)*y;
tem_analytical=zeros(num_n,1);
for i=1:num_n
    tem_analytical(i)=tem_analytical_fun(node_coord(1,i),node_coord(2,i));
end
%% compute local K 
K=zeros(num_n,num_n);
Ke=zeros(4,3,3);
for i=1:4
    Ne=phi(i,:);
    Be=diff(Ne,x)+diff(Ne,y);
    BeTBe=transpose(Be)*Be;
    node_ele=node_coord(:,ele_conn(:,i));
    triangleX = node_coord(1,ele_conn(:,i));
    triangleY = node_coord(2,ele_conn(:,i));
    for j=1:3
        for k=1:3
            fBe = matlabFunction(BeTBe(j,k),'vars',[x y]);
            Ke(i,j,k)=integral2( @(x,y) inpolygon(x, y, triangleX, triangleY) .* fBe(x,y), min(triangleX), max(triangleX), min(triangleY), max(triangleY),'Method', 'iterated');
        end
    end
end
Ke=vpa(Ke,5);
%% solve pde
for i =1:num_e
    K(ele_conn(:,i),ele_conn(:,i))=reshape(Ke(mod(i-1,4)+1,:,:),3,3)+K(ele_conn(:,i),ele_conn(:,i));
end
tem_tilde=zeros(num_n,1);
tem_tilde(nodes_boundary)=tem_analytical(nodes_boundary);
F_modified=-K*tem_tilde;
temh=zeros(num_n,1);
temh(nodes_boundary)=tem_analytical(nodes_boundary);
temh(nodes_inner)=inv(K(nodes_inner,nodes_inner))*F_modified(nodes_inner);
max_diff_fem_between_analytical=abs(max(round(temh-tem_analytical,8)))% The difference between simulated and analytical solution is zero. Pass the patch test!
%% visualize result
figure
[xi,yi] = meshgrid(0:1/m/10:1, 0:1/m/10:1);
zi = griddata(node_coord(1,:),node_coord(2,:),temh,xi,yi);
surf(xi,yi,zi);
set(gca,'xlim',[0 1],'ylim',[0 1]);
grid on
title('simulated temperature surface plot')
xlabel('x')
ylabel('y')
zlabel('temperature')
shading interp
colorbar
saveas(gcf,'simulated temperature surface plot.png')

figure
x = [0 1];
y = [1 0];
image(x,y,zi,'CDataMapping','scaled')
set(gca,'xlim',[0 1],'ylim',[0 1]);
ax = gca;
ax.YDir = 'normal';
grid off
title('simulated temperature heatmap')
xlabel('x')
ylabel('y')
colorbar
saveas(gcf,'simulated temperature heatmap.png')
