clc;clear;close all;                                            	% 确保本文件可靠运行
addpath('.\CalledFunctions');                                       % 获取自定义函数
delta_mat = [[0,0,0];[10,0,0];[0,10,0];[0,0,10];[10,10,10];
    [20,0,0];[0,20,0];[0,0,20];[20,20,20];
    [30,0,0];[0,30,0];[0,0,30];[30,30,30]]./100;
Ori_data_name = '\DataGeneration_SS_SDMVC_RobustnessToGBp_2D';
    
for count_file = 1:size(delta_mat,1)
    delta_G = delta_mat(count_file,1);
    delta_B = delta_mat(count_file,2);
    delta_p = delta_mat(count_file,3);
    DataPath = ['.\data',Ori_data_name,'\delta_',...
            num2str(delta_G.*100),'_',num2str(delta_B.*100),...
            '_',num2str(delta_p.*100),'\'];
    PHPath = ['.\data',Ori_data_name,'\delta_',...
            num2str(delta_G.*100),'_',num2str(delta_B.*100),...
            '_',num2str(delta_p.*100),'\PH\'];
    FigPath = [PHPath,'plot\'];
    if ~exist(FigPath,'dir');mkdir(FigPath);end                 	% 确保创建图片存储路径
    Phantom = importdata('./Phantoms/Phantom_2D.mat');             	% 仿体图片路径+文件名
    
    %% generate noisy data
    dim_reco = importdata([DataPath, 'dim_reco.mat']);
    Phi_reco = importdata([DataPath, 'Phi_reco.mat']);

    Sz = importdata([DataPath,'Sz_reco.mat']);                	    % 加载Sz
    Uz = importdata([PHPath,'Uz.mat']);                             % 加载Uz
    Uz = vertcat(real(Uz),imag(Uz));                                % 重建用接收信号
    dz1 = importdata('.\MeasuredNoise\dz1.mat');                    % 加载dz1
    dz2 = importdata('.\MeasuredNoise\dz2.mat');                    % 加载dz2
    idx_true = false(200,1); idx_true(26:175) = true;               % 2-7 Harmonics with side-bands
    idx_row = repmat(idx_true,[numel(Uz)/200,1]);
    Sz2 = 25*gpuArray(Sz(idx_row,:));                               % * calibration size
    Uz2 = 0.1*gpuArray(Uz(idx_row));                                % * dilution ratio
    dz12 = 0.1*gpuArray(dz1(idx_row));                              % * dilution ratio
    dz22 = 0.1*gpuArray(dz1(idx_row));                              % * dilution ratio

    % The real and imaginary parts are denoised separately
    Uz2_mat = [Uz2(1:end/2),Uz2(1+end/2:end)];                      
    dz12_mat = [dz12(1:end/2),dz12(1+end/2:end)];
    kappa = (vecnorm(Uz2_mat)./vecnorm(dz12_mat))./(10^(40/20));
	dz12_reco = kappa.*dz12_mat;
    Uz_reco = dz12_reco+Uz2_mat;             	                   
    dz12_reco = dz12_reco(:); Uz_reco = Uz_reco(:);
    dz22_reco = kappa.*dz22_mat;
    Sz_reco = Sz2+repmat(dz22_reco(:),[1,size(Sz2,2)]);
    
    %%
    hcf = figure('Name','colorbar','Visible','off'); 
    hax = imagesc(Phantom);axis('square');axis('off');
    hcb = colormap('gray');
    clim([0 1]);
    savefig(hcf,[FigPath,'Phantom.fig']);
    exportgraphics(hcf,[FigPath,'Phantom.png'],'Resolution',600);
    hcb = colorbar;
    hcb.Position(3)=2*hcb.Position(3);
    hcb.Position(1)=hcb.Position(1)-hcb.Position(3);
    hcb.Ticks = [0 0.25 0.5 0.75 1];
    hcb.TickLength = 0.12;
    hcb.FontName = 'Times New Roman';
    hcb.FontSize = 30;
    hax.Visible = 'off';
    cla;
    exportgraphics(hcf,[FigPath,'colorbar.png'],'Resolution',600);
    
    %% reconstruction 
    lambda_vec = logspace(-5,-2,16);
    lambda_vec(13:16) = [];
    Indicator = zeros(3,numel(lambda_vec));
    c_reco_mat = zeros([dim_reco(1),dim_reco(2),numel(lambda_vec)]);
    for j = 1:numel(lambda_vec)   
        lambda = lambda_vec(j);
        opts = struct();                                            	% 算法参数
        opts.maxt = inf;
        opts.SNR = 3;
        opts.tol = 1e-3;
        [c_reco,opts] = TSD_FISTA(Sz_reco,Uz_reco,dz12_reco,Phi_reco,lambda,opts);% TSD_FISTA
        c_reshape = 250*gather(reshape(c_reco,dim_reco));
        c_reco_mat(:,:,j) = c_reshape;
        %% result visualization 0
        error = abs(c_reshape-Phantom);
        PSNR = psnr(c_reshape,Phantom);
        SSIM = ssim(c_reshape,Phantom);
        RMSE = sqrt(sum((c_reshape-Phantom).^2,'all')/numel(Phantom));
        Indicator(:,j) = [PSNR;SSIM;RMSE];
        
        hcf = figure('Name','Conc','Visible','off'); 
        hax = imagesc(c_reshape);axis('square');axis('off');
        hcb = colormap('gray');
        clim([0 1]);
        savefig(hcf,[FigPath,'Conc_',num2str(j),'.fig']);
        exportgraphics(hcf,[FigPath,'Conc_',num2str(j),'.png'],'Resolution',600);

        hcf = figure('Name','Error','Visible','off'); 
        hax = imagesc(error);axis('square');axis('off');
        hcb = colormap('jet');
        clim([0 1]);
        savefig(hcf,[FigPath,'Error_',num2str(j),'.fig']);
        exportgraphics(hcf,[FigPath,'Error_',num2str(j),'.png'],'Resolution',600);
        
        fprintf('Progress1 = %% %f\t Progress2 = %% %f\t toc = %f\n',...
            100*count_file/size(delta_mat,1),100*j/numel(lambda_vec),toc);% 进度
        
    end
    hcb = colorbar;
    hcb.Position(3)=2*hcb.Position(3);
    hcb.Position(1)=hcb.Position(1)-hcb.Position(3);
    hcb.Ticks = [0 0.25 0.5 0.75 1];
    hcb.TickLength = 0.12;
    hcb.FontName = 'Times New Roman';
    hcb.FontSize = 30;
    hax.Visible = 'off';
    cla;
    exportgraphics(hcf,[FigPath,'colorbar_error.png'],'Resolution',600);
    
    save([FigPath,'Indicator.mat'],'Indicator');
    save([FigPath,'c_reco_mat.mat'],'c_reco_mat');
end