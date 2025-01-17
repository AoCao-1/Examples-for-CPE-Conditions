function u_all = gen_rank1_sum0_signal(m, T_i, L, p, alpha, sv)
for i = 1:p
    % Length of the trajectory for the i-th controller
    T_current = T_i(i);
    % Initialize the input matrix for this controller: size m x T_current
    u_all{i} = zeros(m, T_current);

    % ========== 4.1. Construct the first mL control inputs ==========
    % Define key time steps k ∈ {L-1, 2L-1, ..., mL-1}
    excitation_steps = (L-1):L:(m*L-1);  
    existing_excitation_sums = [];  % Store the "weighted sum" vectors of previous key time steps

    for k = 0 : (m*L -1)
        % If k >= T_current, it means the trajectory is too short, break the loop
        if k >= T_current
            break;
        end

        if ismember(k, excitation_steps)
            % Key time step: Ensure that sum_{i=1:p} alpha_i * u_i(k) != 0 
            % and is not linearly dependent on previous excitation vectors
            valid = false;
            attempts = 0; 
            max_attempts = 1000;

            while ~valid && attempts < max_attempts
                % Randomly generate for this controller
                u_candidate = sv*randn(m,1);

                % Compute sum(alpha_i * u_i(k)) for all controllers
                % Note: This is for the "i-th trajectory", should we consider all p controllers?
                % If it’s a "multi-controller" case, we need to generate them together.
                % Assuming 1-to-1 generation for simplicity here, 
                % other controllers’ inputs are set to 0 by default.
                % ----------------------------------
                % For multi-controller simultaneous generation:
                % Use a loop `for ctrl=1:p` and merge checks.
                % ----------------------------------
                % For now, simplifying: generate only for the i-th controller.
                u_all{i}(:, k+1) = u_candidate;

                sum_u = alpha(i) * u_candidate; 

                if norm(sum_u) > 1e-6
                    % Check the linear independence with existing_excitation_sums
                    if isempty(existing_excitation_sums)
                        valid = true;
                    else
                        new_matrix = [existing_excitation_sums, sum_u];
                        if rank(new_matrix) == size(new_matrix, 2)
                            valid = true;
                        end
                    end
                end
                attempts = attempts + 1;
            end

            if ~valid
                error('Trajectory %d: Could not generate valid input for k=%d after %d attempts.', ...
                      i, k, attempts);
            end

            % Record the sum vector for this key time step
            existing_excitation_sums = [existing_excitation_sums, alpha(i)*u_all{i}(:,k+1)];

        else
            for ctrl = 1:(p-1)
                u_all{ctrl}(:,k+1) = sv*randn(m,1);
            end
            % Compute the sum of alpha_i * u_all{i}(:,k+1) for the first p-1 controllers
            sum_prev = zeros(m,1);
            for ctrl = 1:(p-1)
                sum_prev = sum_prev + alpha(ctrl) * u_all{ctrl}(:,k+1);
            end
            % Set the input for the p-th controller to -sum_prev / alpha(p)
            u_all{p}(:,k+1) = -sum_prev / alpha(p);
        end
    end

    % If T_len > mL, fill in the remaining part
    for k = m*L:(T_i(i)-1)
        if k >= T_current
            break;
        end
        % Custom compensation method
        u_all{p}(:,k+1) = zeros(m,1);
    end
end
