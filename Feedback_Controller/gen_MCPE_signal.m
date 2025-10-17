function u_all = gen_MCPE_signal(m, T_i, L, N, para)
%% 3. Calculate Δ_i
Delta = zeros(N, 1); % Initialize Delta as an N+1 element vector
Delta(1) = L - 1;
T_i_s(1) = 0; 
for i = 2:N
    Delta(i) = L - 1 - mod(sum(T_i(1:i-1) - L + 1), L);
    T_i_s(i) = sum(T_i(1:i-1) - L + 1);
end

for i = 1:N
    % Current trajectory length
    T_current = T_i(i);
    Delta_i = Delta(i);
    
    % Initialize the control sequence u_all{i} as m x T_current
    u_all{i} = zeros(m, T_current);
    
    if i == 1
        % For the first trajectory
        if L-1 >= 1
            u_all{i}(:, 1:L-1) = randn(m, L-1); % Randomly generate the first L-1 control inputs
        end
        % Set u_all{i}(:, L) as a specific vector to ensure rank=1
        u_all{i}(:, L) = para * rand(m,1); % Example: first input is 1, others are 0
    else
        % For i > 1 trajectories
        if Delta_i >= 1
            u_all{i}(:, 1:Delta_i) = randn(m, Delta_i); % Randomly generate the first Delta_i control inputs
        end
        
        % Set u_all{i}(:, Delta_i + 1)
        if (T_i(i-1) >= Delta(i-1) + L) && (T_i_s(i) < (m+1)*L - 1)
            % Condition1: u_all{i}(:, Delta_i + 1) should not be in the span of [u_all{1}(:, Delta(1) + 1), ..., u_all{i-1}(:, Delta(i-1) + 1)]
            existing_z = [];
            for j = 1:i-1
                existing_z = [existing_z, u_all{j}(:, Delta(j) + 1)];
            end
            valid = false;
            attempts = 0;
            max_attempts = 1000;  % Prevent infinite loop
            while ~valid && attempts < max_attempts
                candidate = randn(m,1); % Generate a random candidate vector
                % Check if candidate is in the column space of existing_z
                alpha = existing_z \ candidate;
                residual = candidate - existing_z * alpha;
                if norm(residual) > 1e-6 % Not in the column space
                    valid = true;
                end
                attempts = attempts + 1;
            end
            if ~valid
                warning('Could not find a valid u_all{%d}(:, %d). Using a random vector.', i, Delta_i + 1);
                candidate = para * randn(m,1);
            end
            u_all{i}(:, Delta_i + 1) = candidate;
        else
            % Otherwise, set u_all{i}(:, Delta_i + 1) = u_all{i-1}(:, Delta(i-1) + 1)
            u_all{i}(:, Delta_i + 1) = u_all{i-1}(:, Delta(i-1) + 1);
        end
        
        % If Delta_i < L-1, set u_all{i}(:, Delta_i + 2) to u_all{i}(:, L) = 0
        if Delta_i + 1 < L
            u_all{i}(:, Delta_i + 2:L) = 0;
        end
    end
end

for i = 1:N
    % Ensure that control inputs not in positions {L-1, 2L-1, ..., mL-1} are zero
    for k = L+1:T_i(i)
        if ismember(T_i_s(i) + k - 1, (L-1):L:((m+1)*L - 1))
            remaining_steps = T_i(i) - k;
           
            for j = 1:N
                zBar_i = [];
                step = 1;
                while true
                    idx = step * L;
                    if idx < T_i(j) + 1
                         zBar_i = [zBar_i, u_all{j}(:, idx)]; % Concatenate horizontally
                        step = step + 1;
                   else
                       break;
                    end
                end
                 zBar{j} = zBar_i;
            end
            existing_z = [];
            for j = 1:N
                existing_z = [existing_z, u_all{j}(:, Delta(j) + 1), zBar{j}];
            end
            % Store \bar{z}_i (i.e., zBar{i})
    
            if remaining_steps >= L
                % Set z_i(k) not in [z1(Delta1), ..., zp(Deltap), \bar{z}_1, ..., \bar{z}_p]
                valid = false;
                attempts = 0;
                max_attempts = 1000;  % Prevent infinite loop
                candidate = randn(m,1); % Generate a random candidate vector
                % Check if candidate is in the column space of existing_z
                alpha = existing_z \ candidate;
                residual = candidate - existing_z * alpha;
                if norm(residual) > 1e-6 % Not in the column space
                    valid = true;
                end
                attempts = attempts + 1;
                if ~valid
                    warning('Could not find a valid z_i(%d). Using a random vector.', k + 1);
                    candidate = randn(m,1);
                end
                if k <= T_current
                    u_all{i}(:, k) = candidate;
                end
            else
                u_all{i}(:, k) = u_all{i+1}(:, Delta(i+1) + 1);
            end
            % Only set to 0 outside these positions to avoid overwriting non-zero values
        end
    end
end

end
