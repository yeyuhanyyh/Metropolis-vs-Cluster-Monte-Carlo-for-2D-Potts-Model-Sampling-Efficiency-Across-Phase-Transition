%% dyn_compare_at_Tc.m
% Compare Metropolis vs Wolff's IACT, ESS, Cost/ESS at T=Tc
% Paired k-th item: N_list and q_list can be scalar or lists of equal length.
% Record time series by sweep (N^2 flips), calculate tau_int, ESS, Cost/ESS.
clear; clc; rng(2);
%% ----------------------- Experiment Parameters (Paired k-th items) -----------------------
N_list  = [32, 64];   % Can be scalar or list of equal length
q_list  = [3, 3];     % Can be scalar or list of equal length
repeats = 2;
% Sampling settings (shared by both algorithms; add new parameters if separate settings are needed)
burn_sweeps = 100;
keep_sweeps = 1000;
maxlag_dyn  = floor(keep_sweeps/5);
% Align list lengths: allow some to be scalar, others to have length L
L = max([numel(N_list), numel(q_list)]);
assert(ismember(numel(N_list), [1 L]) && ismember(numel(q_list), [1 L]), ...
  'N_list and q_list must have the same length or be scalar');
%% ----------------------- Main Loop (Paired k-th items) -----------------------
run_id = 0;
for k = 1:L
  N = pick_k(N_list, k);
  q = pick_k(q_list, k);
  M = N*N;
  Tc = 1/log(1+sqrt(q));  % J=kB=1
  fprintf('\n=== Dynamics @Tc: N=%d, q=%d (Tc=%.6f) ===\n', N, q, Tc);
  % stats columns: [tau, ESS, flips, Cost/ESS]
  statsM = zeros(repeats, 4);
  statsW = zeros(repeats, 4);
  for rep = 1:repeats
    run_id = run_id + 1;
    % --- Metropolis ---
    sigma0 = randi([1 q], N, N);
    [H_m, flips_m] = metropolis_sweeps_series(sigma0, Tc, burn_sweeps, keep_sweeps, 1, 0, q);
    [tau_m, ~, K_m] = iact_ips(H_m, maxlag_dyn);
    ESS_m = max(numel(H_m)/max(tau_m, eps), 0);   % Protect against division by zero
    cost_m = flips_m / max(ESS_m, eps);
    statsM(rep,:) = [tau_m, ESS_m, flips_m, cost_m];
    % --- Wolff ---
    sigma0 = randi([1 q], N, N);
    [H_w, flips_w] = wolff_sweeps_series(sigma0, Tc, burn_sweeps, keep_sweeps, 1, 0, q);
    [tau_w, ~, K_w] = iact_ips(H_w, maxlag_dyn);
    ESS_w = max(numel(H_w)/max(tau_w, eps), 0);
    cost_w = flips_w / max(ESS_w, eps);
    statsW(rep,:) = [tau_w, ESS_w, flips_w, cost_w];
    fprintf('Run %d | rep %d | Metro: IACT=%.1f (K~%d), ESS=%.0f, Cost/ESS=%.2g | Wolff: IACT=%.1f (K~%d), ESS=%.0f, Cost/ESS=%.2g\n',...
      run_id, rep, tau_m, K_m, ESS_m, cost_m, tau_w, K_w, ESS_w, cost_w);
  end
  muM = mean(statsM,1); seM = std(statsM,0,1)/sqrt(repeats);
  muW = mean(statsW,1); seW = std(statsW,0,1)/sqrt(repeats);
  fprintf('Avg  Metro: IACT=%.1f±%.1f, ESS=%.0f±%.0f, flips=%.2g±%.2g, Cost/ESS=%.2g±%.2g\n',...
    muM(1), seM(1), muM(2), seM(2), muM(3), seM(3), muM(4), seM(4));
  fprintf('Avg  Wolff: IACT=%.1f±%.1f, ESS=%.0f±%.0f, flips=%.2g±%.2g, Cost/ESS=%.2g±%.2g\n',...
    muW(1), seW(1), muW(2), seW(2), muW(3), seW(3), muW(4), seW(4));
end
%% ======================== Local Helper Functions =====================
function v = pick_k(lst, k)
  if numel(lst) == 1, v = lst(1); else, v = lst(k); end
end
function [H_sweeps, flips_total] = metropolis_sweeps_series(sigma, T, burn_sweeps, keep_sweeps, J, h, q)
  [N,~] = size(sigma); beta = 1/T; M=N*N;
  H_cur = calculateH(sigma, J, h);
  total_sweeps = burn_sweeps + keep_sweeps;
  flips_total = 0; rec = 0; H_sweeps = zeros(keep_sweeps,1); sweeps_done = 0;
  while sweeps_done < total_sweeps
    for t = 1:M
      i=randi([1 N]); j=randi([1 N]);
      newv = mod(sigma(i,j)+randi([1 q-1])-1, q)+1;
      dH = calculateDeltaH(newv, sigma, i, j, J, h);
      if dH < 0 || rand < exp(-beta*dH)
        sigma(i,j) = newv; H_cur = H_cur + dH;
      end
      flips_total = flips_total + 1;
    end
    sweeps_done = sweeps_done + 1;
    if sweeps_done > burn_sweeps
      rec = rec + 1; H_sweeps(rec) = H_cur;
      if rec == keep_sweeps, return; end
    end
  end
end
function [H_sweeps, flips_total] = wolff_sweeps_series(sigma, T, burn_sweeps, keep_sweeps, J, h, q)
  [N,~] = size(sigma); beta = 1/T; M=N*N;
  H_cur = calculateH(sigma, J, h);
  p = 1 - exp(-beta*J);
  flips_total = 0; H_sweeps = zeros(keep_sweeps,1);
  sweeps_done = 0; carry = 0; rec = 0;
  while sweeps_done < burn_sweeps + keep_sweeps
    x=randi(N); y=randi(N);
    old = sigma(x,y);
    new = mod(old + randi([1 q-1]) - 1, q) + 1;
    cluster=[x,y]; sigma(x,y)=new; tail=1;
    while ~isempty(cluster)
      [cx,cy] = deal(cluster(end,1), cluster(end,2)); cluster(end,:)=[];
      neighbors = [mod(cx-2,N)+1, cy;
                   mod(cx,  N)+1, cy;
                   cx,            mod(cy-2,N)+1;
                   cx,            mod(cy,  N)+1];
      for ii=1:4
        xn=neighbors(ii,1); yn=neighbors(ii,2);
        if sigma(xn,yn)==old && rand < p
          sigma(xn,yn)=new; cluster=[cluster; xn, yn]; %#ok<AGROW>
          tail = tail + 1;
        end
      end
    end
    H_cur = calculateH(sigma, J, h);
    flips_total = flips_total + tail; carry = carry + tail;
    while carry >= M
      sweeps_done = sweeps_done + 1;
      if sweeps_done > burn_sweeps
        rec = rec + 1; H_sweeps(rec) = H_cur;
        if rec == keep_sweeps, return; end
      end
      carry = carry - M;
    end
  end
end
function [tau_int, rho, K] = iact_ips(x, maxlag)
  x = x(:); n = numel(x);
  x = x - mean(x);
  c0 = sum(x.^2) / n;
  rho = zeros(maxlag,1);
  for k = 1:maxlag
    rho(k) = sum( x(1:n-k).*x(1+k:n) ) / ((n-k)*c0);
    if rho(k) <= 0
      rho = rho(1:k-1); break;
    end
  end
  if isempty(rho), rho = 0; end
  K = numel(rho);
  tau_int = 1 + 2*sum(rho);
end