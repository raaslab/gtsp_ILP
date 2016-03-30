    tic;
    targets_x = guard_target_struct.targets_x; % 
    targets_y = guard_target_struct.targets_y; % 
    
    
    guards_x = guard_target_struct.guards_x;
    guards_y = guard_target_struct.guards_y;
    target_mat = zeros(length(targets_x), length(guards_x));

    for i = 1:length(targets_x)

       target_mat(i,:)  = (visibility_adjacency_matrix((length(guards_x)+i) , (1:length(guards_x)))==1)*i;


    end

    V_Cluster = cell(length(guards_x), 1);

    for i = 1:length(guards_x)
        V_Cluster{i} = target_mat(find(target_mat(:,i)), i)';

    end

    %%
    % Eucledean Distance
    guard_mat = visibility_adjacency_matrix(1:length(guards_x),1:length(guards_x));

    guard_mat = guard_mat - diag(diag(guard_mat));

    guard_mat_weight = zeros(size(guard_mat));

    [r_vis c_vis]= find(guard_mat ~=0); % row column of visiblility 

    for i = 1:length(r_vis)

        guard_mat_weight(r_vis(i), c_vis(i)) = pdist2( [guards_x(r_vis(i)), guards_y(r_vis(i))], [guards_x(c_vis(i)), guards_y(c_vis(i))]);    

    end
    
    %scaled_guard_mat_weight = round(10*guard_mat_weight);


    V_adj  = double(guard_mat_weight);
    
       
    [x_reshape, G_final,fval,exitflag,output] = call_gtsp_func(V_Cluster, V_adj);
    toc;