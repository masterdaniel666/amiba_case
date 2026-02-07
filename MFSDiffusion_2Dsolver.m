function  MFSDiffusion_2Dsolver
clear;clc;format long; warning('off');

Z = [];

for b = 9         %1:1:6    %形狀 6
for c = 4         %1~4 解析解
for d = 5         %0.5:0.5:10  zo 源點隨機參數 7
for e = 1/6         %0:1:9    計算時間 5
for f = e+1/6       %6,10     結束時間 5
for g = 25        %25         內部點密度 
for h = 84       %2:1:100     外框圖形完整度 60 (邊界點數) sp=h*kt+kt
for i = 1       %0.5:0.1:6  波速  (0.1:0.1:2)*60
for j = 1         % 1:1:3     正向、反向、反算(缺象限) 反算(缺初始)
k = [2,4];      % 反算問題保留象限 [2, 4]



para.nbp_buttom = g;     % 布點參數22 ip & bp_2D 內部點密度
para.nbp = h;     % 布點參數60 bp & sp 外框圖形完整度
para.mag = 1;     % 形狀大小
zo = para.mag*d;
% para.ori = [0, 0, 1, zo]; % 2D圖形中心點位
% para.zw = 2;     % 結束時間
para.ori = [0, 0, e, zo]; % 2D圖形中心點位
para.zw = f;     % 結束時間
para.kt = 6;      % 初始至結束時間計算之節點數(層數)
% 形狀選擇 1 Cassani shape 2 Amoeba shape 3 Cardioid shape 
para.sh = b;  % 4 Peanut shape 5 Star shape 6 齒輪 7 cross
para.q = j;          % 1 正向問題 2 反向問題 3 反算問題(缺象限)  4 反算問題(缺初始邊界)
Exact.e = c;       % 解析解選擇 
% Exact.tr = 50;   % time scale factor 
Exact.ap = i;       % 波速
para.keepquad = k; % 反算問題保留象限 [2, 4]
%-------------------------------------------------------------------
[sp,bp,ip,bp2] = CollocationPTs(para);
[A1] = MFSMatrix(bp,sp,para.zw,Exact.ap);
[B1] = Exactsol(bp,Exact.ap,Exact.e);
a=A1\B1;
    %=======Validate
    [A_ip] = MFSMatrix(ip,sp,para.zw,Exact.ap);
    [AS]    = Exactsol(ip,Exact.ap,Exact.e);
    NS = A_ip*a;
    save('dist_data.mat','ip','bp','AS','NS');
%-------------------------------------------------------------------
    ae = (abs(NS - AS))';
    re = (abs(NS - AS)./abs(AS))';
    rmse = sqrt((sum(NS - AS).^2./length(AS)));
    L = sqrt(sum(NS - AS).^2./sum(AS));
    max(ae);
    max(re);
    rmse;
    Plot(bp,bp2,ip,sp,ae,re);
    Z = [Z;max(ae),max(re),rmse,b,c,length(bp),length(sp)]


end
end
end
end
end
end
end
end
end


save Z


5

%=========================================
function [sp,bp,ip,bp2] = CollocationPTs(params)
    nbp_buttom = params.nbp_buttom;
    nbp = params.nbp;
    sh = params.sh;
    mag = params.mag;
    ori = params.ori;
    zw = params.zw;
    kt = params.kt;
    q = params.q;
    %======================================
    %建立邊界 上源點
    [ths,thos] = shape(nbp,sh,mag);
        x_sp(:,1)=ori(1)+thos.*cos(ths);
        y_sp(:,1)=ori(2)+thos.*sin(ths);
        z_sp(1:length(x_sp),1)=ori(3);
        zcs = linspace(ori(4),zw,kt);
        tv2 = length(x_sp);
        for i = 1:kt
            xs((i-1)*tv2+1:i*tv2,1) = x_sp(:);
            ys((i-1)*tv2+1:i*tv2,1) = y_sp(:);
            zs((i-1)*tv2+1:i*tv2,1) = zcs(i);
        end   
        zs = ori(4)*(1+randn(length(zs),1)); % the best
        x_sbp = [xs]; y_sbp = [ys]; z_sbp = [zs]; 
        sp = [x_sbp,y_sbp,z_sbp];



    %建立邊界點Constructing boundary points
    [th1,tho1] = shape(nbp,sh,mag);
        x_bp(:,1)=ori(1)+tho1.*cos(th1);
        y_bp(:,1)=ori(2)+tho1.*sin(th1);
        z_bp(1:length(x_bp),1)=ori(3);
        zc = linspace(ori(3),zw,kt);
        tv2 = length(x_bp);
        for i = 1:kt
            xs((i-1)*tv2+1:i*tv2,1) = x_bp(:);
            ys((i-1)*tv2+1:i*tv2,1) = y_bp(:);
            zs((i-1)*tv2+1:i*tv2,1) = zc(i);
        end   
        x_sbp = [xs]; y_sbp = [ys]; z_sbp = [zs]; 


        %頂部或底部邊界點 bottom or top boundary points
          sbp = [x_sbp,y_sbp];
         [xbm,ybm] = Innerpt(sbp,nbp_buttom);

          % N_ip = 400; % 或任何你想要固定的點數
          
          % [xbm,ybm] = Innerpt(sbp, N_ip);
         zb(1:length(xbm),1) = zc(1);
         zt(1:length(xbm),1) = zc(kt);
        switch q
            case 1 %----- bottom boundary points
         x_bp=[xbm;x_sbp]; y_bp=[ybm;y_sbp]; z_bp=[zb;z_sbp];
            case 2 %----- top boundary points
         x_bp=[xbm;x_sbp]; y_bp=[ybm;y_sbp]; z_bp=[zt;z_sbp];
            case 3 %----- bottom boundary points
         x_bp=[xbm;x_sbp]; y_bp=[ybm;y_sbp]; z_bp=[zb;z_sbp];
            case 4 %----- bottom boundary points
         x_bp=[x_sbp]; y_bp=[y_sbp]; z_bp=[z_sbp];
            otherwise
         error(' case error !!! ')
        end
     bp = [x_bp,y_bp,z_bp];
     bp2 = [x_sbp,y_sbp];
    %驗證點  inner points
        tv1 = length(xbm);
        zc = linspace(ori(3),zw,kt);
        for i = 1:kt
            x_ip(1+(i-1)*tv1:i*tv1,1) = xbm(:);
            y_ip(1+(i-1)*tv1:i*tv1,1) = ybm(:);
            z_ip(1+(i-1)*tv1:i*tv1,1) = zc(i);
        end    
     ip = [x_ip,y_ip,z_ip];

if params.q == 3   % 反算問題下才啟動篩選（正、反向不動）
     keep = params.keepquad;
%=============================================  % 定義每個點屬於哪個象限
    quadrant_bp = zeros(size(x_bp));
    quadrant_bp2 = zeros(size(x_sbp));

    % 篩選邊界點
    quadrant_bp((x_bp >= ori(1)) & (y_bp >= ori(2))) = 1;  % 第一象限
    quadrant_bp((x_bp <= ori(1)) & (y_bp >= ori(2))) = 2;  % 第二象限
    quadrant_bp((x_bp <= ori(1)) & (y_bp <= ori(2))) = 3;  % 第三象限
    quadrant_bp((x_bp >= ori(1)) & (y_bp <= ori(2))) = 4;  % 第四象限
    mask_bp = ismember(quadrant_bp, keep);
    x_bp = x_bp(mask_bp);
    y_bp = y_bp(mask_bp);
    z_bp = z_bp(mask_bp);
    bp = [x_bp, y_bp, z_bp];

    % 篩選邊界投影點（2D）
    quadrant_bp2((x_sbp >= ori(1)) & (y_sbp >= ori(2))) = 1;  % 第一象限
    quadrant_bp2((x_sbp <= ori(1)) & (y_sbp >= ori(2))) = 2;  % 第二象限
    quadrant_bp2((x_sbp <= ori(1)) & (y_sbp <= ori(2))) = 3;  % 第三象限
    quadrant_bp2((x_sbp >= ori(1)) & (y_sbp <= ori(2))) = 4;  % 第四象限
    mask_bp2 = ismember(quadrant_bp2, keep);
    x_sbp = x_sbp(mask_bp2);
    y_sbp = y_sbp(mask_bp2);
    bp2 = [x_sbp, y_sbp];
end
%=========================================
function [x_ip,y_ip] = Innerpt(bp,NN)
    x=linspace(1.0*min(bp(:,1)),1.0*max(bp(:,1)),NN);
    y=linspace(1.0*min(bp(:,2)),1.0*max(bp(:,2)),NN);
    [xx,yy]=meshgrid(x,y);
    [in1]=inpolygon(xx,yy,bp(:,1),bp(:,2));  %選擇小於邊界內之內點x,y
    x_ip = xx(in1);
    y_ip = yy(in1);

% function [x_ip, y_ip] = Innerpt(bp, N_ip)
%     % 建立一個比邊界略大的格點網格
%     x_grid = linspace(min(bp(:,1)), max(bp(:,1)), ceil(sqrt(N_ip)*1.5));
%     y_grid = linspace(min(bp(:,2)), max(bp(:,2)), ceil(sqrt(N_ip)*1.5));
%     [xx, yy] = meshgrid(x_grid, y_grid);
% 
%     % 判斷每個點是否在邊界內部
%     inside_mask = inpolygon(xx, yy, bp(:,1), bp(:,2));
%     x_all = xx(inside_mask);
%     y_all = yy(inside_mask);
% 
%     % 若候選點數 >= N_ip，從中挑 N_ip 個；否則全部用上
%     total_candidates = length(x_all);
%     if total_candidates >= N_ip
%         idx = randperm(total_candidates, N_ip);
%         x_ip = x_all(idx);
%         y_ip = y_all(idx);
%     else
%         warning('可用內部點數少於目標數量，將回傳 %d 個點。', total_candidates)
%         x_ip = x_all;
%         y_ip = y_all;
%     end


%=========================================
function [th,tho] = shape(nbp,sh,mag)    
    th(:,1)=linspace(0,2*pi,nbp+1);
    switch sh
       case 1
        tho=((cos(3*th)+sqrt(2-(sin(3*th)).^2)).^(1/3)); %Cassani shape
       case 2
         tho=(exp(sin(th)).*(sin(2*th)).^2 +exp(cos(th)).*(cos(2*th)).^2); %Amoeba shape
        %tho=mag*((cos(8*th)+sqrt(4-(cos(4*th)).^2)).^(1/3));
       case 3
        tho=sqrt( cos(2*th) + sqrt(1.1-(sin(2*th)).^2)); %Peanut shape  效果差
       case 4
        tho=(1+(cos(4.*th)).^2); %Star shape
       case 5
        tho=(1+0.1*tanh(10*sin(12*th)));%齒輪
       case 6
        tho=mag*(cos(4*th)+(3.6-(sin(4*th).^2)).^(1/2)).^(1/3); %cross
       case 7
        tho=mag*(cos(th).^2+sin(th).^2).^(1/2);  %圓形
       case 8
        tho = abs(cos(3*th)); %花瓣形
        case 9
          s = sin(th);
          c = cos(th);
          tho = 2 - 2*s + (s .* sqrt(abs(c))) ./ (s + 1.4 + eps);
        case 10
           n = 2/3;
    denom = (abs(cos(th)).^n + abs(sin(th)).^n).^(1/n);
    tho = 1 ./ max(denom, eps);
       otherwise
          error(' Boundary Shape error !!! ')
    end
    %
    %(2+1.3*cos(th));%Cardioid shape
    %
    %------------------------------------------
    % normalization
    tho = mag*tho./max(tho);



%=========================================
function [A] = MFSMatrix(op,sp,zw,ap)
    x_op = op(:,1); y_op = op(:,2); z_op = op(:,3);
    x_sp = sp(:,1); y_sp = sp(:,2); z_sp = sp(:,3);
    %====Constructing A matrix for the MFS
    for i = 1:length(x_op)
        for j = 1:length(x_sp)
            rs(i,j)   = ((x_op(i)-x_sp(j))^2 + (y_op(i)-y_sp(j))^2 ...
                            + (z_op(i)-z_sp(j))^2)^0.5;
            r(i,j)   =((x_op(i)-x_sp(j))^2 + (y_op(i)-y_sp(j))^2)^(0.5);
            t(i,j)   =z_op(i,1)-z_sp(j,1);
        end
    end
    %====Constructing A matrix 
    % Laplace Trefftz negative basis function
            for i = 1:length(x_op)
                for j = 1:length(x_sp)
                             A11(i,j) = 1/(2*pi)*log(rs(i,j));
                end
            end
        %====Constructing A matrix for for Time dependent Sol
            % for i = 1:length(x_op)
            %     for j = 1:length(x_sp)
            %             if t(i,j) == 0
            %                t(i,j) = t(i,j)+tr;
            %                apt(i,j) = ap*t(i,j)/tr;
            %             else
            %                apt(i,j) = ap*t(i,j);
            %             end
            %             if t(i,j) < 0
            %                 A12(i,j) = 0;
            %             else
            %                 A12(i,j) = (1 / (4 * pi * apt(i,j))) * exp((-r(i,j)^2) / (4 * apt(i,j)));
            %             end
            %     end
            % end
            for i = 1:length(x_op)
                for j = 1:length(x_sp)
                    H(i,j)=ap*t(i,j)-r(i,j);
                    if H(i,j) <= 0
                        A12(i,j)=0;
                    else
                        % apt(i,j) = ap*t(i,j);
                        % A12(i,j) = ap./(2*pi*sqrt(ap^2*t(i,j).^2-r(i,j).^2));
                        A12(i,j) = 2*ap./(sqrt(ap^2*(t(i,j).^2)-r(i,j).^2));
                    end
                end
            end
            % A12 = (1./(4*pi*apt)).*exp((-r.^2)./(4*apt));
        %====Constructing A matrix
             A = [A11  A12];
            %A=A12;



%===解析解============================================
function [BD] = Exactsol(op,ap,e)
 x = op(:,1);   y = op(:,2); z  = op(:,3);
 r = sqrt((x.^2+y.^2));
switch e
    case 1
      BD = 3+cos(pi*x./10).*cos(pi*y./10).*sin(2^0.5*pi*z./10);
      %The method of fundemental solutions for the multi-dimensional wave equations
    case 3
      BD = 2+sin(pi*x/4).*sin(pi*y/4).*cos(2^0.5*pi*z/4);
      %Domain Type Kernel-Based Meshless Methods for Solving Wave Equations
    % case 3
    %   BD = 3+sin(pi*(r.*x+1)/2).*sin(pi*(r.*y+1)/2).*cos(2^0.5*pi*z/2);
    case 2
      BD = 4+besselj(0,r).*(cos(5*z)+sin(5*z)); 
      %Solving Inverse Wave Problems Using Spacetime Radial Basis Functions in Neural Networks
    %case 4 
      %BD = 3+(x+y)+(x.^2+y.^2).*z+(2/3).*(z.^3);
      % Analytical Solution of Two Dimensional Wave equations by New Laplace Transform Variational Iterative Method
    case 4
      BD = 3+cos(pi*x/10).*cos(pi*y/10).*sin(2^0.5*pi*z/10);
    otherwise
      error(' Exact sol error !!! ')
end




%=========================================
function [] = Plot(bp,bp2,ip,sp,ae,re)
%==Plotting boundary inner source pts
     figure(2)
      plot3(bp(:,1),bp(:,2),bp(:,3),'ro','MarkerSize',9,'MarkerFaceColor',[1 0.6 0.6],'linewidth',1);
      hold on
      %plot3(bp(:,1),bp(:,2),bp(:,3),'-r','LineWidth',1);
      % plot3(ip(:,1),ip(:,2),ip(:,3),'ro', ...
      %  'MarkerFaceColor','b','MarkerSize',5);
      plot3(sp(:,1),sp(:,2),sp(:,3),'bs','MarkerFaceColor','c','MarkerSize',8,'linewidth',1);
      legend('boundary points','source points');
      xmax = max(sp(:,1)); xmin =  min(sp(:,1));
      ymax = max(sp(:,2)); ymin =  min(sp(:,2));
      zmax = max(sp(:,3)); zmin =  min(sp(:,3));
      xlim([floor(xmin), ceil(xmax)]);
      ylim([floor(ymin), ceil(ymax)]);
      zlim([floor(zmin), ceil(zmax)]);
      xticks(floor(xmin):(ceil(xmax)-floor(xmin))/4:ceil(xmax))
      yticks(floor(ymin):(ceil(ymax)-floor(ymin))/4:ceil(ymax))
      axis square
      xlabel('{\it{x}\rm (km)}'), ylabel('{\it{y}\rm (km)}'), zlabel('{\it t\rm (min)}')   % Italic
      set(gca,'fontsize',24,'linewidth',2)
      set(gca,'Fontname', 'Times New Roman')
      grid on
%==Plotting error
figure(5)
xlabel('X','FontSize',10,'FontName','Book Antiqua');
ylabel('Y','FontSize',10,'FontName','Book Antiqua');
plot3(bp(:,1),bp(:,2),bp(:,3),'r');
%hold on
scatter3(ip(:,1),ip(:,2),ip(:,3),[],ae,'filled')
view(-45,10)
title('Absolute error at inner pts')
colorbar
axis equal
xlim([-1,1]);
ylim([-1,1]);
xlabel('\it{x}\rm (km)'), ylabel('\it{y}\rm (km)'), zlabel('\it{t}\rm (min)')   % Italic
set(gca,'fontsize',24,'linewidth',2)
set(gca,'Fontname', 'Times New Roman')
xticks(-1:0.5:1);
yticks(-1:0.5:1);
grid on
% box on
f = gcf;  % 取得當前 figure
f.Position = [200, 200, 1000, 800];
%============================================================
figure(6)
xlabel('X','FontSize',10,'FontName','Book Antiqua');
ylabel('Y','FontSize',10,'FontName','Book Antiqua');
plot3(bp(:,1),bp(:,2),bp(:,3),'r');
%hold on
scatter3(ip(:,1),ip(:,2),ip(:,3),[],re,'filled')
view(-45,10)
title('Relative error at inner pts')
colorbar
axis equal
xlim([-1,1]);
ylim([-1,1]);
xlabel('\it{x}\rm (km)'), ylabel('\it{y}\rm (km)'), zlabel('\it{t}\rm (min)')   % Italic
set(gca,'fontsize',24,'linewidth',2)
set(gca,'Fontname', 'Times New Roman')
xticks(-1:0.5:1);
yticks(-1:0.5:1);
grid on
% box on
f = gcf;  % 取得當前 figure
f.Position = [200, 200, 1000, 800];


 figure(3)
      plot3(bp(:,1),bp(:,2),bp(:,3),'ro','MarkerSize',9,'MarkerFaceColor',[1 0.6 0.6],'linewidth',1);
      hold on
      legend('boundary points');
      xlim([floor(xmin), ceil(xmax)]);
      ylim([floor(ymin), ceil(ymax)]);
      zlim([0, 1]);
      axis square
      xlabel('{\it{x}\rm (km)}'), ylabel('{\it{y}\rm (km)}'), zlabel('{\it t\rm (min)}')   % Italic
      set(gca,'fontsize',24,'linewidth',2)
      set(gca,'Fontname', 'Times New Roman')
      grid on


 figure(4)
      plot(bp2(:,1),bp2(:,2),'ro','MarkerSize',9,'MarkerFaceColor',[1 0.6 0.6],'linewidth',1);
      hold on
      legend('boundary points');
      xlim([-1.5, 1.5]);
      ylim([-1.5, 1.5]);
      xticks(-1.5:0.5:1.5)
      yticks(-1.5:0.5:1.5)
      axis square
      xlabel('{\it{x}\rm (km)}'), ylabel('{\it{y}\rm (km)}'), zlabel('{\it t\rm (min)}')   % Italic
      set(gca,'fontsize',24,'linewidth',2)
      set(gca,'Fontname', 'Times New Roman')
      grid on

figure(50)  %單位  (s) 分鐘換秒鐘
xlabel('X','FontSize',10,'FontName','Book Antiqua');
ylabel('Y','FontSize',10,'FontName','Book Antiqua');
plot3(bp(:,1),bp(:,2),bp(:,3)*60,'r');  %   z 乘上 60
%hold on
scatter3(ip(:,1),ip(:,2),ip(:,3)*60,[],ae,'filled')  % z 乘上 60
view(-45,10)
 title('Absolute error at inner pts')
colorbar
axis equal
xlim([-1,1]);
ylim([-1,1]);
xlabel('\it{x}\rm (km)'), ylabel('\it{y}\rm (km)'), zlabel('\it{t}\rm (s)')   % 單位從 min 改為 s
set(gca,'fontsize',24,'linewidth',2)
set(gca,'Fontname', 'Times New Roman')
xticks(-1:0.5:1);
yticks(-1:0.5:1);
pbaspect([1 1 0.3]);
grid on
view(-45,40)

%============================================================
figure(60)
xlabel('X','FontSize',10,'FontName','Book Antiqua');
ylabel('Y','FontSize',10,'FontName','Book Antiqua');
plot3(bp(:,1),bp(:,2),bp(:,3)*60,'r');
%hold on
scatter3(ip(:,1),ip(:,2),ip(:,3)*60,[],re,'filled')
view(-45,10)
title('Relative error at inner pts')
colorbar
axis equal
xlim([-1,1]);
ylim([-1,1]);
xlabel('\it{x}\rm (km)'), ylabel('\it{y}\rm (km)'), zlabel('\it{t}\rm (s)')   % Italic
set(gca,'fontsize',24,'linewidth',2)
set(gca,'Fontname', 'Times New Roman')
xticks(-1:0.5:1);
yticks(-1:0.5:1);
grid on
view(-45,40)