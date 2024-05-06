function [oper,oper_cell] = SDF(dim)          % sum-difference frame
    dim = reshape(dim,1,length(dim));
    dim(dim==1) = [];
    if length(dim) == 1
        % direction (1) [rows]
        dim_ex = dim+2;
        n_ex = dim_ex;
        idx_ex = (1:n_ex).';
        idx_ori = idx_ex(2:end-1);
        idx_ori = idx_ori(:);
        i = idx_ex(1:end-1);
        j = i+1;
        diff_ex = sparse(i,i,1,n_ex,n_ex)+sparse(i,j,-1,n_ex,n_ex);
        sum_ex = abs(diff_ex);
        oper_ex = [diff_ex;sum_ex];
        oper_ori = oper_ex(:,idx_ori);
        oper_ori(all(oper_ori==0,2),:) = [];
        oper = normalize(oper_ori,1,'norm',2);

        oper_cell = cell(1);
        oper_cell{1} = oper;
    end
    
    if length(dim) == 2
        dim1_ex = dim(1)+2;
        dim2_ex = dim(2)+2;
        n_ex = dim1_ex*dim2_ex;
        idx_ex = reshape(1:n_ex,dim1_ex,dim2_ex);  
        idx_ori = idx_ex(2:end-1,2:end-1);
        idx_ori = idx_ori(:);
        % [rows,cols] direction: (1,0);(-1,1);(0,1);(1,1) 
        d_set = [[1, 0]; [-1, 1]; [0, 1]; [1,1]];
        d = d_set(:,1)+dim1_ex*d_set(:,2);
        omega = ones(4,1);
        omega = normalize(omega,'norm',2)./2;	% normalization of weighted parameters
        % diff_ex = [];
        
        oper_cell = cell(4,1);
        for k = 1:4
            switch k
                case 1  % [rows,cols] direction: (1,0)          
                    i = idx_ex(1:end-1,:);
                case 2  % [rows,cols] direction: (-1,1)
                    i = idx_ex(2:end,1:end-1); 
                case 3  % [rows,cols] direction: (0,1)
                    i = idx_ex(:,1:end-1); 
                case 4  % [rows,cols] direction: (1,1)
                    i = idx_ex(1:end-1,1:end-1);
            end
            i = i(:);
            j = i+d(k);
            diff_ex_temp = sparse(i,i,1,n_ex,n_ex)+sparse(i,j,-1,n_ex,n_ex);  
            sum_ex_temp = abs(diff_ex_temp);
            oper_ex_temp = [diff_ex_temp;sum_ex_temp];
            oper_ori_temp = oper_ex_temp(:,idx_ori);
            oper_ori_temp(all(oper_ori_temp==0,2),:) = [];
            oper_cell{k} = oper_ori_temp;
        end
        oper = cat(1,omega(1)*oper_cell{1},omega(2)*oper_cell{2},...
                    omega(3)*oper_cell{3},omega(3)*oper_cell{3});
%         oper = normalize(oper_ori,1,'norm',2);
    end
    
    if length(dim) == 3
        dim1_ex = dim(1)+2;
        dim2_ex = dim(2)+2;
        dim3_ex = dim(3)+2;
        n_ex = dim1_ex*dim2_ex*dim3_ex;
        idx_ex = reshape(1:n_ex,dim1_ex,dim2_ex,dim3_ex);  
        idx_ori = idx_ex(2:end-1,2:end-1,2:end-1);
        idx_ori = idx_ori(:);
        % [rows,cols,pages] direction: (1,0,0);(-1,1,0);(0,1,0);(1,1,0);
        % [rows,cols,pages] direction: (-1,-1,1);(0,-1,1);(1,-1,1);
        % [rows,cols,pages] direction: (-1,0,1);(0,0,1);(1,0,1);
        % [rows,cols,pages] direction: (-1,1,1);(0,1,1);(1,1,1);
        d_set = [[1, 0, 0]; [-1, 1, 0]; [0, 1, 0]; [1, 1, 0];
                    [-1, -1, 1]; [0, -1, 1]; [1, -1, 1];
                    [-1, 0, 1]; [0, 0, 1]; [1, 0, 1];
                    [-1, 1, 1]; [0, 1, 1]; [1, 1, 1]];
        d = d_set(:,1)+dim1_ex*d_set(:,2)+dim1_ex*dim2_ex*d_set(:,3); 
        omega = ones(13,1);
        omega = normalize(omega,'norm',2)./2;	% normalization of weighted parameters
        
        oper_cell = cell(13,1);
        for k = 1:13
            switch k
                case 1  % [rows,cols] direction: (1,0,0)          
                    i = idx_ex(1:end-1,:,:);
                case 2  % [rows,cols] direction: (-1,1,0)
                    i = idx_ex(2:end,1:end-1,:); 
                case 3  % [rows,cols] direction: (0,1,0)
                    i = idx_ex(:,1:end-1,:); 
                case 4  % [rows,cols] direction: (1,1,0)
                    i = idx_ex(1:end-1,1:end-1,:);
                case 5  % [rows,cols] direction: (-1,-1,1)
                    i = idx_ex(2:end,2:end,1:end-1); 
                case 6  % [rows,cols] direction: (0,-1,1)
                    i = idx_ex(:,2:end,1:end-1); 
                case 7  % [rows,cols] direction: (1,-1,1)
                    i = idx_ex(1:end-1,2:end,1:end-1);
                case 8  % [rows,cols] direction: (-1,0,1)
                    i = idx_ex(2:end,:,1:end-1); 
                case 9  % [rows,cols] direction: (0,0,1)
                    i = idx_ex(:,:,1:end-1); 
                case 10  % [rows,cols] direction: (1,0,1)
                    i = idx_ex(1:end-1,:,1:end-1);
                case 11  % [rows,cols] direction: (-1,1,1)
                    i = idx_ex(2:end,1:end-1,1:end-1); 
                case 12  % [rows,cols] direction: (0,1,1)
                    i = idx_ex(:,1:end-1,1:end-1); 
                case 13  % [rows,cols] direction: (1,1,1)
                    i = idx_ex(1:end-1,1:end-1,1:end-1);
            end
            i = i(:);
            j = i+d(k);
            fprintf('k = %d\n',k);
            diff_ex_temp = sparse(i,i,1,n_ex,n_ex)+sparse(i,j,-1,n_ex,n_ex);
            sum_ex_temp = abs(diff_ex_temp);
            oper_ex_temp = [diff_ex_temp;sum_ex_temp];
            oper_ori_temp = oper_ex_temp(:,idx_ori);
            oper_ori_temp(all(oper_ori_temp==0,2),:) = [];
            oper_cell{k} = oper_ori_temp;
        end
        oper = cat(1,omega(1)*oper_cell{1},omega(2)*oper_cell{2},...
                    omega(3)*oper_cell{3},omega(4)*oper_cell{4},...
                    omega(5)*oper_cell{5},omega(6)*oper_cell{6},...
                    omega(7)*oper_cell{7},omega(8)*oper_cell{8},...
                    omega(9)*oper_cell{9},omega(10)*oper_cell{10},...
                    omega(11)*oper_cell{11},omega(12)*oper_cell{12},...
                    omega(13)*oper_cell{13});
%         oper = normalize(oper_ori,1,'norm',2);
    end
end