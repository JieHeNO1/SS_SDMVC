clc;clear;close all;                                                % 确保本文件可靠运行
addpath('.\CalledFunctions');                                       % 获取自定义函数     
[~,name,~] = fileparts(mfilename('fullpath'));                      % 获取当前文件名% 数据存储路径
delta_mat = [[0,0,0];[10,0,0];[0,10,0];[0,0,10];[10,10,10];
    [20,0,0];[0,20,0];[0,0,20];[20,20,20];
    [30,0,0];[0,30,0];[0,0,30];[30,30,30]]./100;

%% 选用频点
tic;
G = 1.55;                                                           % z方向梯度T/m

r_Max = 15;                                                        	% 垂直于FFL方向最大正向移动范围：mm
a_Max = pi;                                                         % 最大正向旋转范围
z_Max = 15;                                                         % z方向最大正向移动范围：mm
Nr = 15;                                                          	% 径向像素点数
Na = 30;                                                          	% 旋转点数
Nz = 15;                                                            % z方向像素点数，2D用不到
r = linspace(0,r_Max,Nr+1).'; r = r(2:end);                         % 径向采样点
a = linspace(0,a_Max,Na+1).'; 
z = linspace(-z_Max,z_Max,Nz).';                                    % 采样层面
a2 = linspace(0,2*a_Max,2*Na+1).'; a2 = a2(1:end-1);                % 采样角度
[rr,aa] = ndgrid(r,a2);                                           	% 极坐标系网格剖分
[xq,yq] = pol2cart(aa,rr);                                          % 坐标系变换：极坐标系→正交坐标系

x_Max = r_Max;                                                      % x方向最大正向移动范围：mm
y_Max = r_Max;                                                      % y方向最大正向移动范围：mm
Nx = 61;                                                          	% x方向像素点数
Ny = 61;                                                         	% y方向像素点数
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
    
    SM_Datapath = ['./data/',name,'/delta_',...
        num2str(delta_G.*100),'_',num2str(delta_B.*100),'_',num2str(delta_p.*100),'/'];
    if ~exist(SM_Datapath,'dir');mkdir(SM_Datapath);end      	        % 确保创建数据存储路径

    dim_reco = [Nx,Ny];                                              	% 重建图片像素
    save([SM_Datapath, 'dim_reco.mat'],'dim_reco','-v7.3');             % 存储dim_reco
    Phi_reco = SDF(dim_reco);
    save([SM_Datapath, 'Phi_reco.mat'],'Phi_reco','-v7.3');             % 存储Phi_reco
    size_reco = [2*Nf*Na,8*Nr*Na];
    save([SM_Datapath, 'size_reco.mat'],'size_reco','-v7.3');          	% 存储size_reco

    %% system matrix
    tic;
    Sz_2D = zeros(Nf*(Na+1),2*Nr*Na);                                   % 创建Sz_2D
    for jr = 1:Nr                                                       % 采样半径
        ja = 1;                                                         % 采样角度

        fprintf('Progress = %% %f\t toc = %f\n',100*jr/Nr,toc);         % 进度
        for ia = 1:(Na+1)                                               % FFL角度
            drr = r(jr).*cos(a(ia)-a(ja))+dxx(:).';
            [a_ja,r_jr] = cart2pol(drr,dyy(:).');
            z_jz = dzz(:).';

            G_ratio = 1+(r_jr./r_Max).^2.*delta_G;
            B_ratio = 1-(r_jr./r_Max).^2.*delta_B;
            p_ratio = 1-(r_jr./r_Max).^2.*delta_p;

            Ba_P = repmat((G*G_ratio).*(-z_jz),[Nt,1]);                 % 中心平面的平面场强，正方向为a(j1)+pi/2
            Bz_P = -(G_ratio*G).*drr+B_ratio.*B_Dz;                     % z方向场强，正方向为z
            B_P = sqrt(Ba_P.^2+Bz_P.^2);                                % 总场强
            B_P(B_P==0) = eps;                                          % 消除0场强，避免求导问题
            M_P = OSLangevin(B_P,1);                                    % 总磁化矢量
            Mz_P = M_P.*Bz_P./B_P;                                      % z方向磁化矢量
            Uz_T_r = p_ratio.*(diff(Mz_P,1,1)*Fs);                      % z方向接收电压（假设灵敏度均匀且为1）
            Uz_F_r = fft(Uz_T_r,[],1);                                  % 傅里叶变换

            idx_row = (ia-1)*Nf+1:ia*Nf;                                % 填充行
            Sz_2D(idx_row,jr) = mean(Uz_F_r(F_Vec_Idx,:),2);            % 单次测量数据填充
        end
    end
    for ja = 1:2*Na                                                     % 生成其他采样角度
        idx_col = (ja-1)*Nr+1:ja*Nr;                                    % 批量生成所有采样半径数据
        for ia = 1:(Na+1)                                               % FFL角度
            % 其他采样角度矩阵数据映射到采样角度1
            if abs(ia-ja)<Na
                ka = abs(ia-ja)+1;
            else
                ka = 2*Na+1-abs(ia-ja);
            end
            idx_row = (ia-1)*Nf+1:ia*Nf;                                % 采样角度ja填充行
            idx_row1 = (ka-1)*Nf+1:ka*Nf;                               % 采样角度1填充行
            Sz_2D(idx_row,idx_col) = Sz_2D(idx_row1,1:Nr);              % 矩阵填充
        end
    end
    Sz_2D = gather(Sz_2D);
    save([SM_Datapath, 'Sz_2D.mat'],'Sz_2D','-v7.3');                  	% 存储S_2D
    fprintf('Progress = S_2D \t toc = %f\n',toc);                       % 进度

    %% Sz_reco generation
    Sz_2D_cubic = zeros(Nf*Na,Nx*Ny);                                 	% 网格插值平面系统矩阵，180度角度数据删除
    for i = 1:Nf*Na                                                     % 网格插值系统矩阵行索引
        F = scatteredInterpolant(xq(:),yq(:),Sz_2D(i,:).');          	% 生成插值函数
        Sz_2D_cubic(i,:) = F(xx(:),yy(:));                              % 生成插值：正交坐标系平面系统矩阵
    end
    Sz_reco = vertcat(real(Sz_2D_cubic),imag(Sz_2D_cubic));             % 实部、虚部分离

    clearvars Sz_2D Sz_2D_cubic
    save([SM_Datapath, 'Sz_reco.mat'],'Sz_reco','-v7.3');            	% 存储Sz_reco
    fprintf('Progress = Sz_reco \t toc = %f\n',toc);                    % 进度

    %% System Matrix Visualization
    Sz_tmp = Sz_reco(1:end/2,:)+1j*Sz_reco(end/2+1:end,:);          % 还原系统矩阵
    hfig = figure;
    hfig.Position = [962 42 958 953];
    for i = 1:Nf
        count = 0*numel(F_Vec)+i;                                   % FFL角度1对应行
        hax = subplot(14,15,i);
        imagesc(reshape(abs(Sz_tmp(count,:)),dim_reco));            % FFL角度1的PSF        
        hax.XTick = [];hax.YTick = [];
        xlabel(num2str(F_Vec(i)));
    end
    saveas(hfig,[SM_Datapath, 'freqs.png']);                        % 存储PSF图片

    %% signal generation    
    Phantom = importdata('./Phantoms/Phantom_2D.mat');              % 仿体图片路径+文件名
    PH_Datapath = [SM_Datapath,'PH/'];                      	    % 图片存储路径
    if ~exist(PH_Datapath,'dir');mkdir(PH_Datapath);end             % 确保创建图片存储路径
    idx_P = find(Phantom);                                       	% 非0值像素点索引集
    Uz_mat = zeros(Nf*Na,numel(dz));                               	% 接收信号矩阵
    for iz = 1:numel(dz)
        for ia = 1:Na                                             	% FFL旋转角度索引
            [a_jr,r_jr] = cart2pol(xx(idx_P(:)).',yy(idx_P(:)).');
            drr = gpuArray(r_jr.*cos(a(ia)-a_jr));
            z_jz = gpuArray(repmat(dz(iz),[1,numel(drr)]));

            B_ratio = 1-(r_jr./r_Max).^2.*delta_B;
            G_ratio = 1+(r_jr./r_Max).^2.*delta_G;
            p_ratio = 1-(r_jr./r_Max).^2.*delta_p;

            Ba_P = repmat((G*G_ratio).*(-z_jz),[Nt,1]);           	% 中心平面的平面场强，正方向为a(j1)+pi/2
            Bz_P = -(G_ratio*G).*drr+B_ratio.*B_Dz;                	% z方向场强，正方向为z
            B_P = sqrt(Ba_P.^2+Bz_P.^2);                          	% 总场强
            B_P(B_P==0) = eps;                                     	% 消除0场强，避免求导问题
            M_P = OSLangevin(B_P,Phantom(idx_P(:)).');            	% 总磁化矢量
            Mz_P = M_P.*Bz_P./B_P;                                 	% z方向磁化矢量
            Uz_T_r = p_ratio.*(diff(Mz_P,1,1)*Fs);                 	% z方向接收电压（假设灵敏度均匀且为1）
            Uz_F_r = fft(Uz_T_r,[],1);                           	% 傅里叶变换

            idx_row = (ia-1)*Nf+1:ia*Nf;                          	% 填充行
            Uz_mat(idx_row,iz) = sum(Uz_F_r(F_Vec_Idx,:),2);     	% 单次测量数据填充
            fprintf('SubProgress = %% %f\t toc = %f\n',100*ia/Na,toc);	% 进度
        end  
        fprintf('Progress = %% %f\t toc = %f\n',100*iz/numel(dz),toc);	% 进度
    end
    Uz = gather(sum(Uz_mat,2));                                    	% 接收信号数组
    save([PH_Datapath,'Uz.mat'],'Uz','-v7.3');                    	% 存储接收信号
    Uz_noise = awgn(Uz,40,'measured');                              % 接收信号加噪

    Uz_reco = vertcat(real(Uz_noise),imag(Uz_noise));             	% 重建用接收信号
    dz_reco = Uz_reco - vertcat(real(Uz),imag(Uz));                	% 重建用底噪
    save([PH_Datapath,'Uz_reco.mat'],'Uz_reco','-v7.3');         	% 存储Uz_reco
    save([PH_Datapath,'dz_reco.mat'],'dz_reco','-v7.3');         	% 存储dz_reco
    fprintf('BigProgress = %% %f\t toc = %f\n',100*count_file/size(delta_mat,1),toc);	% 进度
end