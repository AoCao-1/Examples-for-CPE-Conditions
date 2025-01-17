function u_all = gen_CCPE_signal(m, T_i, L, p, alpha, sv)
for i = 1:p
    % Length of the i-th controller's trajectory
    T_current = T_i(i);
    % Initialize the control input matrix for the i-th controller: size m x T_current
    u_all{i} = zeros(m, T_current);

    % ========== 4.1. Construct the first mL control inputs ==========
    % Define key time steps k ∈ {L-1, 2L-1, ..., mL-1}
    excitation_steps = (L-1):L:(m*L-1);  
    existing_excitation_sums = [];  % To store the weighted sum vectors of previous key time steps

    for k = 0 : (m*L - 1)
        % If k >= T_current, the trajectory is not long enough, so exit
        if k >= T_current
            break;
        end

        if ismember(k, excitation_steps)
            % Key time steps: Ensure that sum_{i=1:p} alpha_i * u_i(k) != 0 
            % and it is not linearly dependent on previous excitation vectors
            valid = false;
            attempts = 0; 
            max_attempts = 1000;

            while ~valid && attempts < max_attempts
                % Generate a random candidate input for this controller
                u_candidate = sv * randn(m,1);

                % Calculate sum(alpha_i * u_i(k)) for all controllers
                % Note: we are simplifying to generate only for controller i,
                % others are assumed to be zero.
                % For simultaneous generation of p trajectories, iterate over all controllers
                % and merge the checks.
                % This code assumes a 1-to-1 correspondence (only for controller i).
                % -----------------------------------
                % If multiple controllers are being generated simultaneously,
                % this needs to be adapted to: for ctrl = 1:p, u_all{ctrl}(:, k+1) = ...
                % and merge the sum_u checks.
                % -----------------------------------
                % For now, simplify by only generating for controller i:
                u_all{i}(:, k + 1) = u_candidate;

                sum_u = alpha(i) * u_candidate;

                if norm(sum_u) > 1e-6
                    % Check for linear independence from existing excitation sums
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
                error('Trajectory %d: Unable to generate a valid input for k=%d after %d attempts', i, k, attempts);
            end

            % Record the weighted sum vector
            existing_excitation_sums = [existing_excitation_sums, alpha(i) * u_all{i}(:, k + 1)];

        else
            % For non-key time steps (not in excitation_steps)
            for ctrl = 1:(p-1)
                u_all{ctrl}(:, k + 1) = sv * randn(m,1);
            end
            % Calculate the weighted sum for the first p-1 controllers
            sum_prev = zeros(m,1);
            for ctrl = 1:(p-1)
                sum_prev = sum_prev + alpha(ctrl) * u_all{ctrl}(:, k + 1);
            end
            % Set the input for the p-th controller to balance the sum
            u_all{p}(:, k + 1) = -sum_prev / alpha(p);
        end
    end

    % ========== 4.2. If T_current > mL, generate subsequent inputs (k > mL-1) ==========
    for k = m * L : (T_current - 1)
        j = k - m * L + 1; 
        % Construct Z_i^j = [u_i(j), u_i(j+L), ..., u_i(j+(m-1)L)]
        % Here m=2 => [u_i(j), u_i(j+L)]
        % Ensure that j + (m-1) * L <= T_current
        if (j + (m-1) * L) > T_current
            warning('Trajectory %d: Unable to construct input for k=%d due to insufficient length.', i, k);
            break;
        end

        % Z_i^j
        Z_i_j = zeros(m, m);  % m x m matrix
        for col = 1:m
            idx_col = j + (col - 1) * L;
            Z_i_j(:, col) = u_all{i}(:, idx_col);
        end

        % Z_i^0 = [u_i(0), u_i(L)] 
        Z_i_0 = [u_all{i}(:, 1), u_all{i}(:, L + 1)];

        if rank(Z_i_j) < m || rank(Z_i_0) < m
            warning('Trajectory %d: Z_i_j or Z_i^0 is not invertible for k=%d', i, k);
            break;
        end

        % Calculate the control input for u_all{i}(:, k+1) based on previous inputs
        u_all{i}(:, k + 1) = Z_i_j * (Z_i_0 \ u_all{i}(:, m * L));
    end
end
end
