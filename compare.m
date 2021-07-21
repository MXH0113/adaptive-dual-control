clear;
clc;
% close all;

load('experiment_data20170613.mat')
OriData = OriData(201:600,:);
InitTemp = OriData(:,2);         %  x1  unsterilized medium temperature
Speed = OriData(:,3);            %  x2  medium flow rate
SteamTemp = OriData(:,4);        %  x3  steam temperature
% Output = OriData(:,5);           
Multi = OriData(:,1);            
x1 = Speed;
x2 = InitTemp;
x3 = SteamTemp;

% x1(240) = x1(240)+16.8;
x1(329) = x1(329)+5.8;
x2(56) = x2(56)+36.8;
x2(200) = x2(200)-8.8;
% x2(350) = x1(350)-14.8;
% x3(90) = x3(90)+16.8;
x3(230) = x3(230)-6.8;

figure(1)
plot(x1,'b','linewidth',2,'color',[0.20,0.55,1.00]);hold on;
plot(x2,'r','linewidth',2,'color',[1.00,0.55,0.00]);hold on;
plot(x3,'g','linewidth',2,'color',[0.07,0.88,0.10]);hold on;

%% BDC no outlier detection
N = 400;

%real parameter
theta = [ 0.0049; 0.9840; 0.0091; 0.0035; 0.0042; 0.0016 ];

%  reference
y_r = 124*ones(400,1);

% noise
e = 0.01*randn(1,N);
% e(29) = -4;
e(88) = -4.4;
e(114) = 3.6;
e(146) = 8;
e(253) = 5.8;
e(256) = 2.4;

% parameter estimation
theta_hat(:,2) = [0.01; 0.01; 0.01; 0.01; 0.01; 0.01];

% initial value
u = 0*ones(1,400);
y(1) = 100;
y(2) = 100;
y_real(1) = 100;
y_real(2) = 100;
P(:,:,1) = 100*eye(6);
P(:,:,2) = 100*eye(6);

u(2) = 1;
phi(:,2) = [u(2); y(2); y(1); x1(2); x2(2); x3(2)];

lamda=0.01;
gamma = 0;

for k = 2:398
    y_real(k+1) = theta'*phi(:,k);
    y(k+1) = y_real(k+1)+e(k);
       
    K(:,k) = P(:,:,k)*phi(:,k)*inv(0.01+phi(:,k)'*P(:,:,k)*phi(:,k));
    theta_hat(:,k+1) = theta_hat(:,k)+K(:,k)*(y(k+1)-theta_hat(:,k)'*phi(:,k));
    P(:,:,k+1) = (eye(6)-K(:,k)*phi(:,k)')*P(:,:,k);

    %%%%%%  CE
%     u(k+1) = (y_r(k+2)-(theta_hat(2:6,k+1)'*[ y(k+1); y(k); x1(k+1); x2(k+1); x3(k+1)]))/theta_hat(1,k+1);
    
    %%%  BDC
    u(k+1) = -((P(1,2:6,k+1)+theta_hat(1,k+1)*theta_hat(2:6,k+1)')*[ y(k+1); y(k); x1(k+1); x2(k+1); x3(k+1)]-theta_hat(1,k+1)*y_r(k+2))/(P(1,1,k+1) +theta_hat(1,k+1)^2+gamma);
    flag = u(k+1)*P(1,1,k+1)+P(1,2:6,k+1)*[ y(k+1); y(k); x1(k+1); x2(k+1); x3(k+1)];
        if flag>=0
            u(k+1) = u(k+1)+lamda*trace(P(:,:,k+1));
        end
        if flag<0
            u(k+1) = u(k+1)-lamda*trace(P(:,:,k+1));
        end

    phi(:,k+1) = [u(k+1); y(k+1); y(k); x1(k+1); x2(k+1); x3(k+1)];
    
end

y_BDC_no = y;
% plot(y_CE,'b');hold on;
plot(y,'y','linewidth',2,'color',[1.00,0.85,0.00]);hold on;
outlier_k = [56; 88; 114; 146; 200; 253; 256; 230; 329];
outlierplot = zeros(1,400);
plot(outlier_k, outlierplot(outlier_k),'*r','linewidth',2);hold off;
axis([0 400 0 180])
set(gca,'FontSize',14,'Fontname', 'Times New Roman');
xlabel('k','fontsize',14); % x
ylabel('Data','fontsize',14); % y
legend('T_m','R_m','T_s','y(k)')
grid on;


%% BDC+online outlier detetion
clear u  y theta_hat phi K P;

N = 400;

%real parameter
theta = [ 0.0049; 0.9840; 0.0091; 0.0035; 0.0042; 0.0016 ];

%  reference
y_r = 124*ones(400,1);

% noise
%%%%% outlier %%%%%
e = 0.01*randn(1,N);
% e(29) = -4;
e(88) = -4.4;
e(114) = 3.6;
e(146) = 8;
e(253) = 5.8;
e(256) = 2.4;

% parameter estimation
theta_hat(:,2) = [0.01; 0.01; 0.01; 0.01; 0.01; 0.01];

% initial value
u = 0*ones(1,400);
y(1) = 100;
y(2) = 100;
y_real(1) = 100;
y_real(2) = 100;
P(:,:,1) = 100*eye(6);
P(:,:,2) = 100*eye(6);

u(2) = (y_r(3)-(theta_hat(2:6,2)'*[ y(2); y(1); x1(2); x2(2); x3(2)]))/theta_hat(1,2);
u(2) = 1;
phi(:,2) = [u(2); y(2); y(1); x1(2); x2(2); x3(2)];
% phi(:,2) = [60; y(2); y(1); x1(2); x2(2); x3(2)];
lamda=0.01;
gamma = 0;

%outlier detection set
Dimension = 4;              %

% eta
eta_dist = 3.8;
eta_direct =1.55;  %2

% 
dist_store(1) = 0.1;  %  
direct_store(:,1) = [0.1; 0.1; 0.1; 0.1];  % 
nonOutlier = 2;
flag_dist=0;
flag_direct=0;
flag_outlier = 0;

i_dis = 0;
i_dir_1 = 0;
i_dir_2 = 0;
i_dir_3 = 0;
i_dir_4 = 0;

a=0;

for k = 2:398
    y_real(k+1) = theta'*phi(:,k);
    y(k+1) = y_real(k+1)+e(k);
    
    if k>=21

%%%%%%%%%%%%%  online outlier detection   %%%%%%%%%
        dist_store(k) = 0.000000001+norm([y(k); x1(k); x2(k); x3(k)]-[y(k-1); x1(k-1); x2(k-1); x3(k-1)]);         %
        delta_d = std(dist_store(1:k));      % 
        
        % 
        END = 0;       
        for i = 1:k
            END = END +  2*i/k/(k+1)*dist_store(i);
        end          
        if k>=11        
                for i = k-10:k
                     END = END + 2*i/11/(2*k-10)*dist_store(i);
                end
        end
        if k<11   
                for i = 1:k
                    END = END +  2*i/k/(k+1)*dist_store(i);
                end     
        end      

         % 
        BoundUp_dist = END+eta_dist*delta_d;
        BoundDown_dist = END-eta_dist*delta_d;       
        
        dist_up(k) = BoundUp_dist;
        dist_do(k) = BoundDown_dist;
        
         % 
        dist_real = norm([y(k+1); x1(k+1); x2(k+1); x3(k+1)]-[y(k); x1(k); x2(k); x3(k)]);        
        dis_real_p(k) = dist_real;
        
        % 
        flag_dist=0;
        if dist_real<BoundDown_dist  ||  dist_real>BoundUp_dist
            flag_dist = 1;
            i_dis = i_dis+1;
            criti_dis(i_dis) = k;
        end      

    %%%%%   
             direct_store(:,k) = ([y(k); x1(k); x2(k); x3(k)]-[y(k-1); x1(k-1); x2(k-1); x3(k-1)])/(0.00000001+norm([y(k); x1(k); x2(k); x3(k)]-[y(k-1); x1(k-1); x2(k-1); x3(k-1)]));  %
             delta_v = std(direct_store(:,1:k),0,2);      %  

            % 
            ENDV=zeros(Dimension,1);
            if k>=11        
                    for i = k-10:k
                        for j=1:Dimension
                                ENDV(j) = ENDV(j) + 2*i/11/(2*k-10)*direct_store(j,i);
                        end
                    end
            end
            if k<11    
                    for i = 1:k
                        for j=1:Dimension
                                ENDV(j) = ENDV(j) + 2*i/k/(k+1)*direct_store(j,i);
                        end
                    end            
            end

            %
            BoundUp_direct = ENDV+eta_direct*delta_v;
            BoundDown_direct = ENDV-eta_direct*delta_v;  

            dir_up(:,k) = BoundUp_direct;
            dir_do(:,k) = BoundDown_direct;

            %
            direct_real = ([y(k+1); x1(k+1); x2(k+1); x3(k+1)]-[y(k); x1(k); x2(k); x3(3)])/norm([y(k+1); x1(k+1); x2(k+1); x3(k+1)]-[y(k); x1(k); x2(k); x3(k)]);
             dir_real_p(:,k) = direct_real;
        %
        flag_direct=0;
        for i=1:Dimension
            if direct_real(i)<BoundDown_direct(i)  || direct_real(i)>BoundUp_direct(i) 
                flag_direct = 1;
                if i==1
                        i_dir_1 = i_dir_1+1;
                        criti_dir_1(i_dir_1) = k;
                end
                if i==2
                        i_dir_2 = i_dir_2+1;
                        criti_dir_2(i_dir_2) = k;
                end
                if i==3
                        i_dir_3 = i_dir_3+1;
                        criti_dir_3(i_dir_3) = k;
                end
                if i==4
                        i_dir_4 = i_dir_4+1;
                        criti_dir_4(i_dir_4) = k;
                end
            end         
        end    

                % 
            if flag_dist==1 &&flag_direct==1
                a=a+1;
                outlier(a)=k+1;   
                 x1(k+1) = x1(k);
                 x2(k+1) = x2(k);
                 x3(k+1) = x3(k);
                  y(k+1) = y(k);

            else
                flag_outlier = 0;
                nonOutlier = k+1;
            end        

    end
     
    K(:,k) = P(:,:,k)*phi(:,k)*inv(0.01+phi(:,k)'*P(:,:,k)*phi(:,k));
    theta_hat(:,k+1) = theta_hat(:,k)+K(:,k)*(y(k+1)-theta_hat(:,k)'*phi(:,k));
    P(:,:,k+1) = (eye(6)-K(:,k)*phi(:,k)')*P(:,:,k);

%%%%  CE
%     u(k+1) = (y_r(k+2)-(theta_hat(2:6,k+1)'*[ y(k+1); y(k); x1(k+1); x2(k+1); x3(k+1)]))/theta_hat(1,k+1);

%%%  BDC
    u(k+1) = -((P(1,2:6,k+1)+theta_hat(1,k+1)*theta_hat(2:6,k+1)')*[ y(k+1); y(k); x1(k+1); x2(k+1); x3(k+1)]-theta_hat(1,k+1)*y_r(k+2))/(P(1,1,k+1) +theta_hat(1,k+1)^2+gamma);
    flag = u(k+1)*P(1,1,k+1)+P(1,2:6,k+1)*[ y(k+1); y(k); x1(k+1); x2(k+1); x3(k+1)];
        if flag>=0
            u(k+1) = u(k+1)+lamda*trace(P(:,:,k+1));
        end
        if flag<0
            u(k+1) = u(k+1)-lamda*trace(P(:,:,k+1));
        end

    phi(:,k+1) = [u(k+1); y(k+1); y(k); x1(k+1); x2(k+1); x3(k+1)];
    
end
outlier_x = 100*ones(1,a);
y_BDC_online = y;
figure(2)
ref = 124*ones(1,400);
plot(ref,'r','linewidth',2,'color',[1.00,0.55,0.00]);hold on;
plot(y_BDC_no,'-.b','linewidth',2,'color',[0.20,0.55,1.00]);hold on;
plot(y_BDC_online,'g','linewidth',2,'color',[0.07,0.88,0.10]);hold on;

%% BDC+MA
clear u  y theta_hat phi K P;

% x1(240) = x1(240)+16.8;
x1(329) = x1(329)+5.8;
x2(56) = x2(56)+36.8;
x2(200) = x2(200)-8.8;
% x2(350) = x1(350)-14.8;
% x3(90) = x3(90)+16.8;
x3(230) = x3(230)-6.8;

N = 400;

%real parameter
theta = [ 0.0049; 0.9840; 0.0091; 0.0035; 0.0042; 0.0016 ];

%  reference
y_r = 124*ones(400,1);

% noise
%%%%% outlier %%%%%
e = 0.01*randn(1,N);
% e(29) = -4;
e(88) = -4.4;
e(114) = 3.6;
e(146) = 8;
e(253) = 5.8;
e(256) = 2.4;

% parameter estimation
theta_hat(:,2) = [0.01; 0.01; 0.01; 0.01; 0.01; 0.01];
% theta_hat(:,2) = [ 0.0049; 0.9840; 0.0091; 0.0035; 0.0042; 0.0016 ];

% initial value
u = 0*ones(1,400);
y(1) = 100;
y(2) = 100;
y_real(1) = 100;
y_real(2) = 100;
P(:,:,1) = 100*eye(6);
P(:,:,2) = 100*eye(6);

u(2) = (y_r(3)-(theta_hat(2:6,2)'*[ y(2); y(1); x1(2); x2(2); x3(2)]))/theta_hat(1,2);
u(2) = 1;
phi(:,2) = [u(2); y(2); y(1); x1(2); x2(2); x3(2)];
% phi(:,2) = [60; y(2); y(1); x1(2); x2(2); x3(2)];
lamda=0.01;
gamma = 0;

%outlier detection set
Dimension = 4;              %

%
eta_dist = 3.8;
eta_direct =1.55;  %2

%
dist_store(1) = 0.1;  
direct_store(:,1) = [0.1; 0.1; 0.1; 0.1];  
nonOutlier = 2;
flag_dist=0;
flag_direct=0;
flag_outlier = 0;

i_dis = 0;
i_dir_1 = 0;
i_dir_2 = 0;
i_dir_3 = 0;
i_dir_4 = 0;

a=0;


for k = 2:398
    y_real(k+1) = theta'*phi(:,k);
    y(k+1) = y_real(k+1)+e(k);

%%%%%%%%%%%%%  Moving Average   %%%%%%%%%%%%%
    if k>=21
        for nn = 1:10
            y(k+1) = y(k+1-nn)+y(k+1);
            x1(k+1) = x1(k+1-nn)+x1(k+1);
            x2(k+1) = x2(k+1-nn)+x2(k+1);
            x3(k+1) = x3(k+1-nn)+x3(k+1);
        end
        y(k+1) = 1/11*y(k+1);
        x1(k+1) = 1/11*x1(k+1);
        x2(k+1) = 1/11*x2(k+1);
        x3(k+1) = 1/11*x3(k+1);
    end

     
    K(:,k) = P(:,:,k)*phi(:,k)*inv(0.01+phi(:,k)'*P(:,:,k)*phi(:,k));
    theta_hat(:,k+1) = theta_hat(:,k)+K(:,k)*(y(k+1)-theta_hat(:,k)'*phi(:,k));
    P(:,:,k+1) = (eye(6)-K(:,k)*phi(:,k)')*P(:,:,k);

%%%%  CE
%     u(k+1) = (y_r(k+2)-(theta_hat(2:6,k+1)'*[ y(k+1); y(k); x1(k+1); x2(k+1); x3(k+1)]))/theta_hat(1,k+1);

%%%  BDC
    u(k+1) = -((P(1,2:6,k+1)+theta_hat(1,k+1)*theta_hat(2:6,k+1)')*[ y(k+1); y(k); x1(k+1); x2(k+1); x3(k+1)]-theta_hat(1,k+1)*y_r(k+2))/(P(1,1,k+1) +theta_hat(1,k+1)^2+gamma);
    flag = u(k+1)*P(1,1,k+1)+P(1,2:6,k+1)*[ y(k+1); y(k); x1(k+1); x2(k+1); x3(k+1)];
        if flag>=0
            u(k+1) = u(k+1)+lamda*trace(P(:,:,k+1));
        end
        if flag<0
            u(k+1) = u(k+1)-lamda*trace(P(:,:,k+1));
        end

    phi(:,k+1) = [u(k+1); y(k+1); y(k); x1(k+1); x2(k+1); x3(k+1)];
    
end
y_BDC_MA = y;
plot(y_BDC_MA,'--y','linewidth',2,'color',[1.00,0.85,0.00]);hold on;

%% BDC+MH
clear u  y theta_hat phi K P;

% x1(240) = x1(240)+16.8;
x1(329) = x1(329)+5.8;
x2(56) = x2(56)+36.8;
x2(200) = x2(200)-8.8;
% x2(350) = x1(350)-14.8;
% x3(90) = x3(90)+16.8;
x3(230) = x3(230)-6.8;

N = 400;

%real parameter
theta = [ 0.0049; 0.9840; 0.0091; 0.0035; 0.0042; 0.0016 ];

%  reference
y_r = 124*ones(400,1);

% noise
%%%%% outlier %%%%%
e = 0.01*randn(1,N);
% e(29) = -4;
e(88) = -4.4;
e(114) = 3.6;
e(146) = 8;
e(253) = 5.8;
e(256) = 2.4;

% parameter estimation
theta_hat(:,2) = [0.01; 0.01; 0.01; 0.01; 0.01; 0.01];

% initial value
u = 0*ones(1,400);
y(1) = 100;
y(2) = 100;
y_real(1) = 100;
y_real(2) = 100;
P(:,:,1) = 100*eye(6);
P(:,:,2) = 100*eye(6);

u(2) = (y_r(3)-(theta_hat(2:6,2)'*[ y(2); y(1); x1(2); x2(2); x3(2)]))/theta_hat(1,2);
u(2) = 1;
phi(:,2) = [u(2); y(2); y(1); x1(2); x2(2); x3(2)];
% phi(:,2) = [60; y(2); y(1); x1(2); x2(2); x3(2)];
lamda=0.01;
gamma = 0;

%outlier detection set
Dimension = 4;              %

% 
eta_dist = 3.8;
eta_direct =1.55;  %2

%
dist_store(1) = 0.1;  % 
direct_store(:,1) = [0.1; 0.1; 0.1; 0.1];  % 
nonOutlier = 2;
flag_dist=0;
flag_direct=0;
flag_outlier = 0;

i_dis = 0;
i_dir_1 = 0;
i_dir_2 = 0;
i_dir_3 = 0;
i_dir_4 = 0;

a=0;


for k = 2:398
    y_real(k+1) = theta'*phi(:,k);
    y(k+1) = y_real(k+1)+e(k);
    delta_y(k+1) = (y(k+1)-theta_hat(:,k)'*phi(:,k))^2;

%%%%%%%%%%%%%  Moving Horizon   %%%%%%%%%%%%%
    J_MH(:,1) = zeros(5,1);
    if k>=21
        sum = delta_y(k+1)+delta_y(k)+delta_y(k-1)+delta_y(k-2)+delta_y(k-3);
        for nn = 1:5
            J_MH(nn,1) = sum - delta_y(k+nn-4);
        end
        [value, iter] = min(J_MH);
        if iter==5
            y(k+1) = y(k);
        end  
    end

     
    K(:,k) = P(:,:,k)*phi(:,k)*inv(0.01+phi(:,k)'*P(:,:,k)*phi(:,k));
    theta_hat(:,k+1) = theta_hat(:,k)+K(:,k)*(y(k+1)-theta_hat(:,k)'*phi(:,k));
    P(:,:,k+1) = (eye(6)-K(:,k)*phi(:,k)')*P(:,:,k);

%%%%  CE
%     u(k+1) = (y_r(k+2)-(theta_hat(2:6,k+1)'*[ y(k+1); y(k); x1(k+1); x2(k+1); x3(k+1)]))/theta_hat(1,k+1);

%%%  BDC
    u(k+1) = -((P(1,2:6,k+1)+theta_hat(1,k+1)*theta_hat(2:6,k+1)')*[ y(k+1); y(k); x1(k+1); x2(k+1); x3(k+1)]-theta_hat(1,k+1)*y_r(k+2))/(P(1,1,k+1) +theta_hat(1,k+1)^2+gamma);
    flag = u(k+1)*P(1,1,k+1)+P(1,2:6,k+1)*[ y(k+1); y(k); x1(k+1); x2(k+1); x3(k+1)];
        if flag>=0
            u(k+1) = u(k+1)+lamda*trace(P(:,:,k+1));
        end
        if flag<0
            u(k+1) = u(k+1)-lamda*trace(P(:,:,k+1));
        end

    phi(:,k+1) = [u(k+1); y(k+1); y(k); x1(k+1); x2(k+1); x3(k+1)];
    
end
y_BDC_MH = y;
plot(y_BDC_MH,':m','linewidth',2);hold on;
plot(outlier, outlier_x,'*r','linewidth',2);hold off;
axis([0 400 100 145])
set(gca,'FontSize',14,'Fontname', 'Times New Roman');
xlabel('k','fontsize',14); % 
ylabel('y(k)','fontsize',14); % 
legend('reference','without outlier detection','online outlier detection','moving average', 'robust moving horizon','outliers')
grid on;



