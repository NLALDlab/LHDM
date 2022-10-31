function [resnorm,true_nnz,x_nnz,true_xl1norm,x_l1norm, opt_dist, time_el, iter_hist, maxcounter] = test_marco_FM(A,x,b,itest,NNEG,runsolver,solver_options,do_plots)


iter_hist = zeros(1);

maxcounter = 0;

n_of_repetitions = 5;
velt = zeros(1,n_of_repetitions);

resnorm = -1.0;
time_el = -1.0;
x_nnz = -1;
opt_dist = -1.0;
x_l1norm = -1.0;

true_nnz = -1;
true_xl1norm = -1.0;

true_nnz = nnz(x);
Ic_true = find(x);
true_xl1norm = norm(x,1);
[m,n] = size(A);


% METODO: L1 MAGIC
if runsolver.l1_magic % l1-magic: https://statweb.stanford.edu/~candes/software/l1magic/



    if issparse(A)
        return
    end
    if i==3 | i==26 | i==51 | i==55 | i==61 | i==65 | i==74 | i==77 | i==78 | i==104 | i==112 | i==136 | i==137 | i==147  | i==234 | i==277 | i==278 | i==332 | i==379    
        return 
    end 

    
    if issparse(A)
        return
    end
    
    disp("l1-magic:")
    % initial guess = min energy
    x0 = A' * b;
    % solve the LP
    tic
        x_comp = l1eq_pd(x0, A, [], b, 1e-11); % era 1e-3
    elapsed = toc



% METODO: L1 HOMOTOPY
elseif runsolver.l1_homotopy


    if issparse(A)
        return
    end
    if i==130 %| i==26
        return 
    end 
    disp("l1-homotopy:")
    err_fun = @(z) (norm(x-z)/norm(x))^2;
    SNR = 100; %40;       % additive Gaussian noise
    sigma = sqrt(norm(A'*b)^2/10^(SNR/10)/m);
    tau = 1.e-11;  %max(1e-4*max(abs(A'*b)),sigma*sqrt(log(n)));
    % rank-1 update mode 
    delx_mode = 'qr'; % mil or qr
    in = [];   
    in.tau = tau;
    in.delx_mode = delx_mode;
    in.debias = 0;
    in.verbose = 0;
    in.plots = 0;
    in.record = 1;
    in.err_fun = err_fun;
    in.early_terminate = false;
    
    tic
    out = l1homotopy(A,b,in);
    elapsed = toc
    
    x_comp = out.x_out;
    iter_bpdn = out.iter;
    time_bpdn = out.time;
    gamma_bpdn = out.gamma;
    %err_bpdn = out.error_table;
    
    
% METODO: L1 HOMOTOPY
elseif runsolver.l1_homotopy_early



    if issparse(A)
        return
    end

    disp("l1-homotopy early_terminate:")
    err_fun = @(z) (norm(x-z)/norm(x))^2;
    SNR = 40;       % additive Gaussian noise
    sigma = sqrt(norm(A'*b)^2/10^(SNR/10)/m);
    tau = max(1e-4*max(abs(A'*b)),sigma*sqrt(log(n)));
    % rank-1 update mode 
    delx_mode = 'mil'; % mil or qr
    in = [];   
    in.tau = tau;
    in.delx_mode = delx_mode;
    in.debias = 0;
    in.verbose = 0;
    in.plots = 0;
    in.record = 1;
    in.err_fun = err_fun;
    in.early_terminate = true;
    
    tic
    out = l1homotopy(A,b,in);
    elapsed = toc
    
    x_comp = out.x_out;
    iter_bpdn = out.iter;
    time_bpdn = out.time;
    gamma_bpdn = out.gamma;
    %err_bpdn = out.error_table;

% METODO: SolveBP
elseif runsolver.SolveBP

    if issparse(A)
        return
    end


    disp("SolveBP:")
    maxIters = 500;
    lambda = 0; % era 0;
    OptTol = 1e-11; % era 1e-3;
    tic
    x_comp = SolveBP(A, b, n, maxIters, lambda, OptTol);
    elapsed = toc


% METODO: SolveOMP
elseif runsolver.SolveOMP

    if issparse(A)
        return
    end


    disp("SolveOMP:")
    maxIters = 2000;
    lambdaStop = 0; % era 0;
    solFreq = 0;
    verbose = 0;
    OptTol = 1e-14; %1e-11; % era 1e-3;
    for ir = 1:n_of_repetitions
        tic
        x_comp = SolveOMP(A, b, n, maxIters, lambdaStop, solFreq, verbose, OptTol);
        elapsed = toc;
        velt(ir) = elapsed;
    end
    elapsed = min(velt);
    disp(['elapsed_time = ',num2str(elapsed)]);


% METODO: SPGL1
elseif runsolver.SPGL1

    if issparse(A)
        return
    end

    disp("SPGL1:")
    opts = spgSetParms('verbosity',0,'optTol',1e-9);
    tic
    x_comp = spg_bp(A, b, opts);
    elapsed = toc
    
% METODO: SPGL1 tol
elseif runsolver.SPGL1_tol
    
    % tol = 1e-08
    %if i==154 | i==155 | i==238 | i==239 | i==251 | i==252 | i==254 | i==264 | i==265 | i==266
    %    return
    %end
    
    % tol = 1e-11
    %if i<=8 || (i>=10 && i<=40)          % e altre
    %    return
    %end
    
    
    if issparse(A)
        return
    end
    
    disp("SPGL1 tol:")
  
    opts = spgSetParms('verbosity',0, 'optTol', 1e-8, 'decTol', 1e-8);
    tic
    x_comp = spg_bp(A, b, opts);
    elapsed = toc


% METODO: ISAL1
elseif runsolver.ISAL1

    %if itest==3 || itest==4 || itest==6
    %    return
    %end
    
    if issparse(A)
        return
    end


    disp("ISAL1:")
    tic
    [x_comp,fval,err,exfl,it] = ISAL1(A,b,[]); %,varargin)
    elapsed = toc


% METODO: PBP
elseif runsolver.PBP

    %if itest==3 || itest==4 || itest==6
    %    return
    %end
    
    if issparse(A)
        return
    end


    disp("PBP:")
    tic
    [x_comp] = blocknnls(A,b); %,varargin)
    elapsed = toc



% METODO: LHDM
elseif runsolver.LHDM
    if (issparse(A))
        return
    end
    disp("LHDM:")
    LHDM_options = solver_options;
    LHDM_options.x_true = x;
    verbose = 0;
    for ir = 1:n_of_repetitions
        if NNEG
            tic
            [tmpx,resnorm_LHDM_v4,exitflag, outeriter_LHDM_v4, itervec_LHDM_v4] = LHDM(A,b,LHDM_options,verbose);
            elapsed = toc;
            x_comp = tmpx;
        else
            tic
            [tmpx,resnorm_LHDM_v4,exitflag, outeriter_LHDM_v4, itervec_LHDM_v4] = LHDM([A -A],b,LHDM_options,verbose);
            elapsed = toc;
            x_comp = tmpx(1:n) - tmpx(n+1:end);
        end
        velt(ir) = elapsed;
    end
    elapsed = min(velt);
    disp(['elapsed_time = ',num2str(elapsed)]);
    iter_hist = itervec_LHDM_v4;

end


resnorm = norm(b - A*x_comp); disp(['resnorm = ',num2str(resnorm)]);
time_el = elapsed;
x_nnz = nnz(x_comp); disp(['x_comp_nnz = ',num2str(x_nnz),' , nnz(x) = ',num2str(nnz(x))]);
opt_dist = norm(x_comp - x); %/norm(x); 
disp(['opt_dist = ',num2str(opt_dist)]); tmpI = find(abs(x)+abs(x_comp)); fullx=full(x); %[tmpI'; x_comp(tmpI)'; fullx(tmpI)']
x_l1norm = norm(x_comp, 1); disp(['x_l1norm = ',num2str(x_l1norm),' , ||x||_1 = ',num2str(norm(x,1))]);
