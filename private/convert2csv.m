% Convert MAT file data to CSV (preprocessed data)

% Get original MAT file from here: 
% https://github.com/lacerbi/visvest-causinf/blob/master/VestBMS_ModelWork/VestBMS_data.mat

temp = load('VestBMS_data.mat');
data = temp.data;
Ns = numel(data);

data_mat = [];
for iSubj = 1:Ns
    data_one = [];
    for iNoise = 1:3
       for iTask = 2:3
           X = data{iSubj}.X.bimodal{iNoise}{iTask};
           Nx = size(X,1);
           data_one = [data_one; [iSubj*ones(Nx,1), X, iNoise*ones(Nx,1)]];
       end
    end
    
    % Sort by trial number
    [~,idx] = sort(data_one(:, 2));
    data_one = data_one(idx,:);
    
    data_mat = [data_mat; data_one];
end

csvwrite('bisensory_data.csv',data_mat);