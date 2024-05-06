clc;clear;close all;                                                % 确保本文件可靠运行
addpath(genpath('.\CalledFunctions'));                          	% 获取自定义函数     
[~,name,~] = fileparts(mfilename('fullpath'));                      % 获取当前文件名
delta_mat = [[0,0,0];[10,0,0];[0,10,0];[0,0,10];[10,10,10];
    [20,0,0];[0,20,0];[0,0,20];[20,20,20];
    [30,0,0];[0,30,0];[0,0,30];[30,30,30]]./100;

%% 选用频点
tic;
G = 1.55;                                                           % z方向梯度T/m

r_Max = 15;                                                        	% 垂直于FFL方向最大正向移动范围：mm
a_Max = pi;                                                         % 最大正向旋转范围
z_Max = 15;                                                         % z方向最大正向移动范围：mm
Nr = 8;                                                          	% 径向像素点数
Na = 15;                                                          	% 旋转点数
Nz = 31;                                                            % z方向像素点数，2D用不到
r = linspace(1,r_Max,Nr).';                                         % 径向采样点
a = linspace(0,a_Max,Na+1).'; 
z = linspace(-z_Max,z_Max,Nz).';                                    % 采样层面
a2 = linspace(0,2*a_Max,2*Na+1).'; a2 = a2(1:end-1);                % 采样角度
[rr,aa] = ndgrid(r,a2);                                           	% 极坐标系网格剖分
[xq,yq] = pol2cart(aa,rr);                                          % 坐标系变换：极坐标系→正交坐标系

x_Max = r_Max;                                                      % x方向最大正向移动范围：mm
y_Max = r_Max;                                                      % y方向最大正向移动范围：mm
Nx = 31;                                                          	% x方向像素点数
Ny = 31;                                                         	% y方向像素点数
x = linspace(-x_Max,x_Max,Nx).';                                    % x方向采样位置
y = linspace(-y_Max,y_Max,Ny).';                                    % y方向采样位置
[xx,yy] = ndgrid(x,y);                                              % 正交坐标系网格剖分
[aq,rq] = cart2pol(xx,yy);                                          % 坐标系变换：正交坐标系→极坐标系，用于构建

%% calibration sample
dx = gpuArray(0);
dy = gpuArray((-12:12).');
dz = gpuArray(0);
[dxx,dyy,dzz] = ndgrid(dx,dy,dz);

%% sampling
Fs= 3e6;                                                            % 采样率
Fz1 = 3e3;                                                          % 激励频率
Fz2 = 120;                                                          % 扫描频率
totaltime = 0.1;                                                    % 总采样时长
t = linspace(0,totaltime,totaltime*Fs+1).';                       	% 采样时间点，0时刻采样点在磁化矢量微分时去掉
Nt = length(t);                                                     % 采样数

%% drive and rotate field                                        
B_Dz1_Amp = 7;                                                      % 激励场幅值
B_Dz2_Amp = 14;                                                     % 扫描场幅值
B_Dz = B_Dz1_Amp*cos(2*pi*Fz1*t)+B_Dz2_Amp*cos(2*pi*Fz2*t);         % z方向激励场强分布：mT
B_Dz = gpuArray(B_Dz);

%% Frequency components
[F_Vec_Idx,F_Vec] = MixFreqs([Fz1,Fz2],[8,12],Fs,totaltime,30e3); 	% 选用频点（索引）
Nf = numel(F_Vec_Idx);                                              % 选用频点数

tic;
for count_file = 1:size(delta_mat,1)
    delta_G = delta_mat(count_file,1);
    delta_B = delta_mat(count_file,2);
    delta_p = delta_mat(count_file,3);
    
    SM_DataPath = ['./data/',name,'/delta_',...
        num2str(delta_G.*100),'_',num2str(delta_B.*100),'_',num2str(delta_p.*100),'/'];
    if ~exist(SM_DataPath,'dir');mkdir(SM_DataPath);end                 % 确保创建数据存储路径
    Phantom = importdata('./Phantoms/Phantom_3D.mat');                  % 仿体图片路径+文件名
	PH_DataPath = [SM_DataPath,'PH/'];                                  % 图片存储路径
    if ~exist(PH_DataPath,'dir');mkdir(PH_DataPath);end                 % 确保创建图片存储路径

    dim_reco = [Nx,Ny,Nz];                                             	% 重建图片像素
    save([SM_DataPath, 'dim_reco.mat'],'dim_reco','-v7.3');          	% 存储dim_reco
    Phi_reco = SDF(dim_reco);
    save([SM_DataPath, 'Phi_reco.mat'],'Phi_reco','-v7.3');           	% 存储Phi_reco
    size_reco = [2*Nf*Na*Nz,Nx*Ny*Nz];
    save([SM_DataPath, 'size_reco.mat'],'size_reco','-v7.3');          	% 存储size_reco

    %% system matrix
    tic;
    for iz = 1:Nz                                                       % 激励层索引
        jz = 1;                                                         % 利用平面场强Br_P均匀性，仅采样最底/高层
        Sz_2D = zeros(Nf*(Na+1),2*Nr*Na);                               % 创建Sz_2D
        Sz_2D = gpuArray(Sz_2D);
        for jr = 1:Nr                                               	% 采样半径
            ja = 1;                                                   	% 采样角度
            for ia = 1:(Na+1) %                                       	% FFL角度
                drr = r(jr).*cos(a(ia)-a(ja))+dxx(:).';
                [a_j,r_j] = cart2pol(drr,dyy(:).');
                z_j = z(jz)+dzz(:).';

                B_ratio = 1-(z_j./z_Max).^2.*delta_B;
                G_ratio = 1+((z_j-z(iz))./(2*z_Max)).^2.*delta_G;
                p_ratio = 1-(z_j./z_Max).^2.*delta_p;

                Ba_P = repmat((G*G_ratio).*(z(iz)-z_j),[Nt,1]);         % 中心平面的平面场强，正方向为a(j1)+pi/2
                Bz_P = -(G_ratio*G).*drr+B_ratio.*B_Dz;                 % z方向场强，正方向为z
                B_P = sqrt(Ba_P.^2+Bz_P.^2);                            % 总场强
                B_P(B_P==0) = eps;                                      % 消除0场强，避免求导问题
                M_P = OSLangevin(B_P,1);                                % 总磁化矢量
                Mz_P = M_P.*Bz_P./B_P;                                  % z方向磁化矢量
                Uz_T_r = p_ratio.*(diff(Mz_P,1,1)*Fs);                  % z方向接收电压（假设灵敏度均匀且为1）
                Uz_F_r = fft(Uz_T_r,[],1);                              % 傅里叶变换

                idx_row = (ia-1)*Nf+1:ia*Nf;                           	% 填充行
                Sz_2D(idx_row,jr) = mean(Uz_F_r(F_Vec_Idx,:),2);        % 单次测量数据填充
            end
        end
        for ja = 1:2*Na                                                	% 生成其他采样角度
            idx_col = (ja-1)*Nr+1:ja*Nr;                               	% 批量生成所有采样半径数据
            for ia = 1:(Na+1)                                         	% FFL角度
                % 其他采样角度矩阵数据映射到采样角度1
                if abs(ia-ja)<Na
                    ka = abs(ia-ja)+1;
                else
                    ka = 2*Na+1-abs(ia-ja);
                end
                idx_row = (ia-1)*Nf+1:ia*Nf;                         	% 采样角度ja填充行
                idx_row1 = (ka-1)*Nf+1:ka*Nf;                         	% 采样角度1填充行
                Sz_2D(idx_row,idx_col) = Sz_2D(idx_row1,1:Nr);         	% 矩阵填充
            end
        end
        Sz_2D = gather(Sz_2D);
        save([SM_DataPath, 'Sz_2D_',num2str(iz),'.mat'],'Sz_2D','-v7.3');% 存储S_2D
    end
    fprintf('Progress1 = %% %.2f\t toc = %.2f\n',100*iz/Nz,toc);     % 进度

    %% Sz_reco generation
    tic;
    for iz = 1:Nz
        Sz_2D = importdata([SM_DataPath, 'Sz_2D_',num2str(iz),'.mat']);
        Sz_2D_cubic = zeros(Nf*Na,Nx*Ny);                               % 网格插值平面系统矩阵，180度角度数据删除
        for i = 1:Nf*Na                                                	% 网格插值系统矩阵行索引
            F = scatteredInterpolant(xq(:),yq(:),Sz_2D(i,:).');       	% 生成插值函数
            Sz_2D_cubic(i,:) = F(xx(:),yy(:));                         	% 生成插值：正交坐标系平面系统矩阵
        end
        Sz_reco = vertcat(real(Sz_2D_cubic),imag(Sz_2D_cubic));         % 实部、虚部分离    
        save([SM_DataPath, 'Sz_reco_',num2str(iz),'.mat'],'Sz_reco','-v7.3');% 存储Sz_reco
        fprintf('Progress2 = %% %.2f\t toc = %.2f\n',100*iz/Nz,toc);  	% 进度
    end
    clearvars Sz_2D Sz_2D_cubic

    %% signal generation   
    tic;
	idx_z = find(sum(Phantom,[1,2]));                                   % 非零z层索引集
    Uz_mat1 = zeros(2*Nf*Na*Nz,numel(idx_z));                           % 接收信号矩阵
    for iz = 1:Nz                                                       % 激励层索引
        for jjz = 1:numel(idx_z)
            jz = idx_z(jjz);                                          	% 非零z层索引
            B_ratio = 1-(z(jz)/z_Max)^2*delta_B;
            G_ratio = 1+((z(jz)-z(iz))/(2*z_Max))^2*delta_G;
            p_ratio = 1-(z(jz)/z_Max)^2*delta_p;
            Phantom_z = Phantom(:,:,jz);                               	% 非零z层仿体
            idx_xy = find(Phantom_z);                                 	% 非0值像素点索引集
            Uz_mat2 = zeros(Nf*Na,numel(idx_xy));                       % 接收信号矩阵
            for jj_xy = 1:numel(idx_xy)                                	% 非0值像素点索引集索引
                j_xy = idx_xy(jj_xy);                                	% 非0值像素点索引
                for ia = 1:Na                                          	% FFL旋转角度索引
                    Ba_P = (G*G_ratio)*(z(iz)-z(jz));                	% 中心平面的平面场强，正方向为a(j1)+pi/2       
                    Bz_P = -(G_ratio*G)*(rq(j_xy).*cos(a(ia)-aq(j_xy)))+B_ratio.*B_Dz;% z方向场强，正方向为z 
                    B_P = sqrt(Ba_P.^2+Bz_P.^2);                      	% 总场强
                    B_P(B_P==0) = eps;                                 	% 消除0场强，避免求导问题
                    M_P = OSLangevin(B_P,Phantom_z(j_xy));            	% 总磁化矢量
                    Mz_P = M_P.*Bz_P./B_P;                             	% z方向磁化矢量
                    Uz_T_r = p_ratio.*(diff(Mz_P)*Fs);                	% z方向接收电压（假设灵敏度均匀且为1）
                    Uz_F_r = fft(Uz_T_r);                              	% 傅里叶变换
                    idx_row = (ia-1)*Nf+1:ia*Nf;                        % 接收信号子矩阵行索引
                    Uz_mat2(idx_row,jj_xy) = Uz_F_r(F_Vec_Idx);         % 接收信号子矩阵赋值
                end
            end
            idx_row = 2*Nf*Na*(iz-1)+1:2*Nf*Na*iz;                     	% 接收信号子矩阵行索引
            Uz_mat1(idx_row,jjz) = vertcat(real(sum(Uz_mat2,2)),imag(sum(Uz_mat2,2)));% 接收信号子矩阵赋值
            fprintf('SubProgress = %% %f\t toc = %f\n',...
                100*(jjz+(iz-1)*numel(idx_z))/(numel(idx_z)*Nz),toc);   % 进度
        end
        fprintf('Progress3 = %% %f\t toc = %f\n',100*iz/Nz,toc);      % 进度
    end
    Uz = sum(Uz_mat1,2);                                                % 接收信号数组
    save([PH_DataPath,'Uz.mat'],'Uz','-v7.3');                          % 存储接收信号
    Uz_reco = awgn(Uz,40,'measured');                                   % 接收信号加噪
    dz_reco = Uz_reco - Uz;                                             % 重建用底噪
    save([PH_DataPath,'Uz_reco.mat'],'Uz_reco','-v7.3');                % 存储Uz_reco
    save([PH_DataPath,'dz_reco.mat'],'dz_reco','-v7.3');                % 存储dz_reco
    fprintf('\nBigProgress = %% %f\t toc = %f\n\n',100*count_file/size(delta_mat,1),toc);      % 进度
end