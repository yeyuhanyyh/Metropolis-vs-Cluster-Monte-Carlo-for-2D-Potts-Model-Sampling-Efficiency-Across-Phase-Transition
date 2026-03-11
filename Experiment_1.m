%% scan_equal_work.m
% Scan and plot u(T), c(T), comparing Metropolis vs Wolff
% Standardize "work" = spin flips. Each run takes the k-th item from all parameter lists.
% Dependencies: calculateH, calculateDeltaH, and local helper functions at the end of this file.
clear; clc; rng(1);

%% ----------------------- Experiment Parameters (Paired k-th items) -----------------------
% Default parameter values are chosen to be good for quick experiments and sanity checks,
% although results will not be that accurate.
% For more accurate result, refer to the parameter setting in the experiment section of the paper.
N_list               = [32];    % List of lattice sizes (N)
q_list               = [3];      % List of q values (must align with N_list or be scalar)
W_perT_metro_list    = [300];  % Metropolis: target sweeps per temperature (1 sweep = N^2 flips)
W_perT_wolff_list    = [200];  % Wolff: target sweeps per temperature (1 sweep = N^2 flips)
repeats              = 1;           % Number of repeats for each parameter set
% Temperature grid and burn-in config (Annealed warm-start: from high T to low T)
T_values         = linspace(1.4, 0.6, 25);
scan_burn_sweeps = 50;              % Number of sweeps to discard at each temperature point
% Align parameter list lengths: allow some to be scalar, others to have length L
L = max([numel(N_list), numel(q_list), numel(W_perT_metro_list), numel(W_perT_wolff_list)]);
assert(all(ismember([numel(N_list), numel(q_list), numel(W_perT_metro_list), numel(W_perT_wolff_list)], [1, L])), ...
  'All parameter lists must have length 1 or the same length L.');
%% ----------------------- Main Loop (Paired k-th items) -----------------------
run_id = 0;
for k = 1:L
  N            = pick_k(N_list, k);
  q            = pick_k(q_list, k);
  W_perT_metro = pick_k(W_perT_metro_list, k);
  W_perT_wolff = pick_k(W_perT_wolff_list, k);
  M  = N*N;
  Tc = 1/log(1+sqrt(q));   % kB=J=1
  % Statistics containers
  uM_all = zeros(repeats, numel(T_values));
  cM_all = zeros(repeats, numel(T_values));
  uW_all = zeros(repeats, numel(T_values));
  cW_all = zeros(repeats, numel(T_values));
  for rep = 1:repeats
    run_id = run_id + 1;
    fprintf('Run %d | Scan: N=%d, q=%d, Wm/N^2=%d, Ww/N^2=%d, rep=%d\n', ...
            run_id, N, q, W_perT_metro, W_perT_wolff, rep);
    % Annealed warm-start: from high T to low T; functions must return sigma_* to preserve it for the next T point
    sigma_m = randi([1 q], N, N);
    sigma_w = randi([1 q], N, N);
    u_m = zeros(1, numel(T_values)); c_m = u_m;
    u_w = zeros(1, numel(T_values)); c_w = u_w;
    for t_index = 1:numel(T_values)
      T    = T_values(t_index);
      beta = 1/T; 
      % --- Metropolis: Collect samples for equivalent work ---
      W_target_m = W_perT_metro * M;
      [Hs_m, ~, sigma_m] = collect_sweeps_equal_work_metro(sigma_m, T, W_target_m, q);
      H_use = drop_burn(Hs_m, scan_burn_sweeps);
      u_m(t_index) = mean(H_use)/M;
      c_m(t_index) = ( (1/T)^2 * var(H_use) )/M;
      % --- Wolff: Collect samples for equivalent work ---
      W_target_w = W_perT_wolff * M;
      [Hs_w, ~, sigma_w] = collect_sweeps_equal_work_wolff(sigma_w, T, W_target_w, q);
      H_use = drop_burn(Hs_w, scan_burn_sweeps);
      u_w(t_index) = mean(H_use)/M;
      c_w(t_index) = ( (1/T)^2 * var(H_use) )/M;
    end
    uM_all(rep,:) = u_m;  cM_all(rep,:) = c_m;
    uW_all(rep,:) = u_w;  cW_all(rep,:) = c_w;
  end
  % Calculate mean and standard error
  mu_uM = mean(uM_all,1); se_uM = std(uM_all,0,1)/sqrt(repeats);
  mu_cM = mean(cM_all,1); se_cM = std(cM_all,0,1)/sqrt(repeats);
  mu_uW = mean(uW_all,1); se_uW = std(uW_all,0,1)/sqrt(repeats);
  mu_cW = mean(cW_all,1); se_cW = std(cW_all,0,1)/sqrt(repeats);
  % --- Plot: u(T) ---
  figure('Name', sprintf('u(T): N=%d, q=%d, Wm/N^2=%d, Ww/N^2=%d',N,q,W_perT_metro,W_perT_wolff)); hold on;
  errorbar(T_values, mu_uM, se_uM, '-o', 'DisplayName','Metropolis');
  errorbar(T_values, mu_uW, se_uW, '-s', 'DisplayName','Wolff');
  xline(Tc,'r--','T_c'); xlabel('T'); ylabel('u(T)');
  title(sprintf('u(T), N=%d, q=%d, Wm/N^2=%d, Ww/N^2=%d', N, q, W_perT_metro, W_perT_wolff));
  legend('Location','best'); grid on;
  % --- Plot: c(T) ---
  figure('Name', sprintf('c(T): N=%d, q=%d, Wm/N^2=%d, Ww/N^2=%d',N,q,W_perT_metro,W_perT_wolff)); hold on;
  errorbar(T_values, mu_cM, se_cM, '-o', 'DisplayName','Metropolis');
  errorbar(T_values, mu_cW, se_cW, '-s', 'DisplayName','Wolff');
  xline(Tc,'r--','T_c'); xlabel('T'); ylabel('c(T)');
  title(sprintf('c(T), N=%d, q=%d, Wm/N^2=%d, Ww/N^2=%d', N, q, W_perT_metro, W_perT_wolff));
  legend('Location','best'); grid on;
end
%% ======================== Local Helper Functions =========================
function v = pick_k(lst, k)
  if numel(lst) == 1
    v = lst(1);
  else
    v = lst(k);
  end
end
function H_use = drop_burn(Hs, burn)
  if numel(Hs)>burn, H_use = Hs(burn+1:end);
  else, H_use = Hs;
  end
end
function [H_sweeps, flips, sigma] = collect_sweeps_equal_work_metro(sigma, T, W_target, q)
  % Single-spin proposal + Metropolis acceptance; 1 sweep per N^2 flips
  [N,~] = size(sigma); beta = 1/T; M = N*N;
  H_cur = calculateH(sigma, 1, 0);
  flips = 0; in_sweep = 0; H_sweeps = [];
  while flips < W_target
    i = randi([1 N]); j = randi([1 N]);
    newv = mod(sigma(i,j)+randi([1 q-1])-1, q)+1;
    dH = calculateDeltaH(newv, sigma, i, j, 1, 0);
    if dH < 0 || rand < exp(-beta*dH)
      sigma(i,j) = newv;
      H_cur = H_cur + dH;
    end
    flips = flips + 1;
    in_sweep = in_sweep + 1;
    if in_sweep >= M
      H_sweeps(end+1,1) = H_cur; 
      in_sweep = 0;
    end
  end
end
function [H_sweeps, flips, sigma] = collect_sweeps_equal_work_wolff(sigma, T, W_target, q)
  % Cluster growth + immediate recolor; accumulate cluster size as flips; 1 sweep per >= N^2 accumulated flips
  [N,~] = size(sigma); beta = 1/T; M=N*N; p = 1 - exp(-beta*1);
  H_cur = calculateH(sigma, 1, 0);
  flips = 0; carry = 0; H_sweeps = [];
  while flips < W_target
    x = randi(N); y = randi(N);
    old = sigma(x,y);
    new = mod(old + randi([1 q-1]) - 1, q) + 1;
    cluster = [x, y]; sigma(x,y) = new; tail = 1;
    while ~isempty(cluster)
      [cx, cy] = deal(cluster(end,1), cluster(end,2));
      cluster(end,:) = [];
      neighbors = [mod(cx-2,N)+1, cy;
                   mod(cx,  N)+1, cy;
                   cx,            mod(cy-2,N)+1;
                   cx,            mod(cy,  N)+1];
      for ii=1:4
        xn = neighbors(ii,1); yn = neighbors(ii,2);
        if sigma(xn,yn)==old && rand<p
          sigma(xn,yn)=new; cluster=[cluster; xn, yn];
          tail = tail + 1;
        end
      end
    end
    H_cur = calculateH(sigma, 1, 0);
    flips = flips + tail; carry = carry + tail;
    while carry >= M
      H_sweeps(end+1,1) = H_cur;
      carry = carry - M;
    end
  end
end