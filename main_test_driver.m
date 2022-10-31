addpath('L1_Testset_ascii')
addpath('L1_Testset_ascii/Dense')
addpath('L1_Testset_ascii/Sparse')
addpath('L1_homotopy_v2.0')
addpath('L1_homotopy_v2.0/utils')
addpath('L1_homotopy_v2.0/solvers')
addpath('l1magic')
addpath('l1magic/Optimization')
addpath('SparseLab2.1-Core');
addpath('SparseLab2.1-Core/Solvers');
addpath('SparseLab2.1-Core/Utilities');
addpath('spgl1-2.1');
addpath('spgl1-2.1/spgl_in');
addpath('ISAL1_v1.0');
addpath('ISAL1_v1.0/ISAL1bpdn');
addpath('Principal_block_pivoting');

matrix_type=1;
case_study = 1;

% number of matrices in the dataset.
nmat = 548;
nsparse = 0;

l1_magicnnz_vec = -1 + zeros(nmat, 1);
l1_magictime_vec = -1 + zeros(nmat, 1);
l1_magicopt_dist_vec = -1 + zeros(nmat, 1);
l1_magicresnorm_vec = -1 + zeros(nmat, 1);
l1_magicx_l1norm_vec = -1 + zeros(nmat, 1);

l1_homotopynnz_vec = -1 + zeros(nmat, 1);
l1_homotopytime_vec = -1 + zeros(nmat, 1);
l1_homotopyopt_dist_vec = -1 + zeros(nmat, 1);
l1_homotopyresnorm_vec = -1 + zeros(nmat, 1);
l1_homotopyx_l1norm_vec = -1 + zeros(nmat, 1);

l1_homotopy_earlynnz_vec = -1 + zeros(nmat, 1);
l1_homotopy_earlytime_vec = -1 + zeros(nmat, 1);
l1_homotopy_earlyopt_dist_vec = -1 + zeros(nmat, 1);
l1_homotopy_earlyresnorm_vec = -1 + zeros(nmat, 1);
l1_homotopy_earlyx_l1norm_vec = -1 + zeros(nmat, 1);

SolveBPnnz_vec = -1 + zeros(nmat, 1);
SolveBPtime_vec = -1 + zeros(nmat, 1);
SolveBPopt_dist_vec = -1 + zeros(nmat, 1);
SolveBPresnorm_vec = -1 + zeros(nmat, 1);
SolveBPx_l1norm_vec = -1 + zeros(nmat, 1);

SolveOMPnnz_vec = -1 + zeros(nmat, 1);
SolveOMPtime_vec = -1 + zeros(nmat, 1);
SolveOMPopt_dist_vec = -1 + zeros(nmat, 1);
SolveOMPresnorm_vec = -1 + zeros(nmat, 1);
SolveOMPx_l1norm_vec = -1 + zeros(nmat, 1);

SPGL1nnz_vec = -1 + zeros(nmat, 1);
SPGL1time_vec = -1 + zeros(nmat, 1);
SPGL1opt_dist_vec = -1 + zeros(nmat, 1);
SPGL1resnorm_vec = -1 + zeros(nmat, 1);
SPGL1x_l1norm_vec = -1 + zeros(nmat, 1);

SPGL1_tolnnz_vec = -1 + zeros(nmat, 1);
SPGL1_toltime_vec = -1 + zeros(nmat, 1);
SPGL1_tolopt_dist_vec = -1 + zeros(nmat, 1);
SPGL1_tolresnorm_vec = -1 + zeros(nmat, 1);
SPGL1_tolx_l1norm_vec = -1 + zeros(nmat, 1);

YALL1nnz_vec = -1 + zeros(nmat, 1);
YALL1time_vec = -1 + zeros(nmat, 1);
YALL1opt_dist_vec = -1 + zeros(nmat, 1);
YALL1resnorm_vec = -1 + zeros(nmat, 1);
YALL1x_l1norm_vec = -1 + zeros(nmat, 1);

ISAL1nnz_vec = -1 + zeros(nmat, 1);
ISAL1time_vec = -1 + zeros(nmat, 1);
ISAL1opt_dist_vec = -1 + zeros(nmat, 1);
ISAL1resnorm_vec = -1 + zeros(nmat, 1);
ISAL1x_l1norm_vec = -1 + zeros(nmat, 1);

PBPnnz_vec = -1 + zeros(nmat, 1);
PBPtime_vec = -1 + zeros(nmat, 1);
PBPopt_dist_vec = -1 + zeros(nmat, 1);
PBPresnorm_vec = -1 + zeros(nmat, 1);
PBPx_l1norm_vec = -1 + zeros(nmat, 1);

LHDM_v1nnz_vec = -1 + zeros(nmat, 1);
LHDM_v1time_vec = -1 + zeros(nmat, 1);
LHDM_v1opt_dist_vec = -1 + zeros(nmat, 1);
LHDM_v1resnorm_vec = -1 + zeros(nmat, 1);
LHDM_v1x_l1norm_vec = -1 + zeros(nmat, 1);

LHDM_v3_NNREGnnz_vec = -1 + zeros(nmat, 1);
LHDM_v3_NNREGtime_vec = -1 + zeros(nmat, 1);
LHDM_v3_NNREGopt_dist_vec = -1 + zeros(nmat, 1);
LHDM_v3_NNREGresnorm_vec = -1 + zeros(nmat, 1);
LHDM_v3_NNREGx_l1norm_vec = -1 + zeros(nmat, 1);

LHDM_v1_hard_thresholdednnz_vec = -1 + zeros(nmat, 1);
LHDM_v1_hard_thresholdedtime_vec = -1 + zeros(nmat, 1);
LHDM_v1_hard_thresholdedopt_dist_vec = -1 + zeros(nmat, 1);
LHDM_v1_hard_thresholdedresnorm_vec = -1 + zeros(nmat, 1);
LHDM_v1_hard_thresholdedx_l1norm_vec = -1 + zeros(nmat, 1);

OLSnnz_vec = -1 + zeros(nmat, 1);
OLStime_vec = -1 + zeros(nmat, 1);
OLSopt_dist_vec = -1 + zeros(nmat, 1);
OLSresnorm_vec = -1 + zeros(nmat, 1);
OLSx_l1norm_vec = -1 + zeros(nmat, 1);

LH_v3nnz_vec = -1 + zeros(nmat, 1);
LH_v3time_vec = -1 + zeros(nmat, 1);
LH_v3opt_dist_vec = -1 + zeros(nmat, 1);
LH_v3resnorm_vec = -1 + zeros(nmat, 1);
LH_v3x_l1norm_vec = -1 + zeros(nmat, 1);

LHDM_v3nnz_vec = -1 + zeros(nmat, 1);
LHDM_v3time_vec = -1 + zeros(nmat, 1);
LHDM_v3opt_dist_vec = -1 + zeros(nmat, 1);
LHDM_v3resnorm_vec = -1 + zeros(nmat, 1);
LHDM_v3x_l1norm_vec = -1 + zeros(nmat, 1);

LHDM_v4nnz_vec = -1 + zeros(nmat, 1);
LHDM_v4time_vec = -1 + zeros(nmat, 1);
LHDM_v4opt_dist_vec = -1 + zeros(nmat, 1);
LHDM_v4resnorm_vec = -1 + zeros(nmat, 1);
LHDM_v4x_l1norm_vec = -1 + zeros(nmat, 1);

LH_Cnnz_vec = -1 + zeros(nmat, 1);
LH_Ctime_vec = -1 + zeros(nmat, 1);
LH_Copt_dist_vec = -1 + zeros(nmat, 1);
LH_Cresnorm_vec = -1 + zeros(nmat, 1);
LH_Cx_l1norm_vec = -1 + zeros(nmat, 1);

LHDM_Cnnz_vec = -1 + zeros(nmat, 1);
LHDM_Ctime_vec = -1 + zeros(nmat, 1);
LHDM_Copt_dist_vec = -1 + zeros(nmat, 1);
LHDM_Cresnorm_vec = -1 + zeros(nmat, 1);
LHDM_Cx_l1norm_vec = -1 + zeros(nmat, 1);

iter_vec = zeros(nmat,1);

% choose here the solvers to execute (1=execute, 0=not-execute):
runsolver = struct('l1_magic',0, ...
                   'l1_homotopy',0, ...
                   'l1_homotopy_early', 0, ...
                   'SolveBP',0, ...
                   'SolveOMP',0, ...
                   'SPGL1',0, ...
                   'SPGL1_tol', 0,... 
                   'YALL1',0, ...
                   'ISAL1',0, ...
                   'PBP',0, ...
                   'LHDM',1);

nrepeats = nmat;
vrepeats = [1:nrepeats];

% choose to run a modified dataset, respectively, non-negative, denser (with more non-zeros) or ill-conditioned, by setting to "1" one of the following variables:
NNEG = 0;
modify_nnz = 0;
ILLCOND = 0;


% LHDM setup:
var_cos = 4*0.1;
c_thres_nrm = 2;
var_w = 6*0.1;

% NB: "k" is the block size; with "k=1" runs the LH algorithm
k = 1 %16;

positrick = 0;
ensure_descend_direction = 1;
k_adattativo = 0; % not used

t = 0.0;
solver_options = struct('init',false, ...
                        'itmax',2500, ...
                        'k',k, ...
                        'positrick',positrick, ...
                        't',t,...
                        'thres_cos', var_cos, ... 
                        'thres_nrm',c_thres_nrm*0.1, ... 
                        'thres_w',var_w, ... 
                        'NNEG',NNEG, ...
                        'tol',2.e-15, ...
                        'k_adattativo',k_adattativo, ...
                        'ensure_descend_direction',ensure_descend_direction ...
                        );  
                        %thres_cos = 0.8, thres_nrm = 0.8, thres_w = 0.8

ticStart = tic;

for i=vrepeats
    itest = i; 
    fprintf("itest = %d \n",itest)
    
    if matrix_type == 1
        i = itest;   
        [A,x,b,specs] = read_ascii_instance(itest);

        if NNEG && issparse(A)==0 
            x = abs(x);
            continua = 1;
            diag_term = 1;
            [m,n] = size(A);
            I = find(abs(x) > 0);
            tmpn = length(I)
            while continua
                A(:,I) = A(:,I) + diag_term*eye(m,tmpn);
                % test the ERC condition:
                pinvAST = pinv(A(:,I));
                vn1AS = zeros(1,n-tmpn);
                Iz = find(x == 0);
                for jc = Iz'
                    vn1AS(jc) = norm(pinvAST*A(:,jc),1);
                end
                disp(['ERC condition: ',num2str(max(vn1AS)),' < 1 ?']);
                diag_term = 1.1;
                if max(vn1AS) < 1
                    continua = 0;
                end
            end
            b = A * x;
        end
        if modify_nnz && issparse(A)==0 
            [m,n] = size(A);
            I = find(abs(x) > 0);
            tmpn = length(I)
            % test the ERC condition:
            pinvAST = pinv(A(:,I)');
            vn1AS = zeros(1,n-tmpn);
            Iz = find(x == 0);
            for jc = Iz'
                vn1AS(jc) = norm(A(:,jc)'*pinvAST,1);
            end
            disp(['ERC condition: ',num2str(max(vn1AS)),' < 1 ?']);

            new_nnz = floor(m*1.0);
            continua = 1;
            diag_term = 1;
            I = floor(rand(1,new_nnz)*n)+1;
            tmpn = length(I)
            maxx = max(abs(x));
            x = zeros(n,1);
            x(I) = maxx*(rand(tmpn,1) - 0.5);
            while continua
                A(:,I) = A(:,I) + diag_term*eye(m,tmpn);
                % test the ERC condition:
                pinvAST = pinv(A(:,I));
                vn1AS = zeros(1,n-tmpn);
                Iz = find(x == 0);
                for jc = Iz'
                    vn1AS(jc) = norm(pinvAST*A(:,jc),1);
                end
                disp(['ERC condition: ',num2str(max(vn1AS)),' < 1 ?']);
                diag_term = 1.1;
                if max(vn1AS) < 1
                    continua = 0;
                end
            end
            b = A * x;
            disp(['cond(A) = ',num2str(cond(A))]);
        end
        if ILLCOND && issparse(A)==0 
            [m,n] = size(A);
            I = find(abs(x) > 0);
            tmpn = length(I);
            condAI_prima = cond(A(:,I))
            % test the ERC condition:
            pinvAST = pinv(A(:,I)');
            vn1AS = zeros(1,n-tmpn);
            Iz = find(x == 0);
            for jc = Iz'
                vn1AS(jc) = norm(A(:,jc)'*pinvAST,1);
            end
            disp(['ERC condition: ',num2str(max(vn1AS)),' < 1 ?']);
            [U,dummy] = qr(rand(m,tmpn));
            [V,dummy] = qr(rand(tmpn,tmpn));
            tmpc = 20/tmpn; % max 45
            vs = 0.5.^(tmpc*[1:tmpn]);
            A(:,I) = U(:,1:tmpn) * diag(vs) * V';
            for jc = I'
                %nc_prima=norm(A(:,jc))
                A(:,jc) = A(:,jc) / norm(A(:,jc));
                %nc_dopo=norm(A(:,jc))
            end
            condAI_dopo = cond(A(:,I))
            A(:,Iz) = A(:,Iz) ./ (condAI_dopo/4);
            % test the ERC condition:
            pinvAST = pinv(A(:,I)');
            vn1AS = zeros(1,n-tmpn);
            for jc = Iz'
                vn1AS(jc) = norm(A(:,jc)'*pinvAST,1);
            end
            %vn1AS
            disp(['min vn1AS = ',num2str(min(vn1AS))]);
            disp(['ERC condition: ',num2str(max(vn1AS)),' < 1 ?']);
            b = A * x;
            x_prova = zeros(size(x));
            x_prova(I) = A(:,I) \ b;
            disp(['relerr x_prova = ',num2str(norm(x_prova-x)/norm(x))]);
            disp(['residual x_prova = ',num2str(norm(b - A*x_prova))]);
        end
    end
    
    disp(['size(A) = ',mat2str(size(A))])
    
    
    
    if runsolver.ISAL1 == 1
        tmp_runsolver = struct('l1_magic',0,'l1_homotopy',0,'l1_homotopy_early',0,'SolveBP',0,'SolveOMP',0,'SPGL1',0,'SPGL1_tol',0,'YALL1',0,'ISAL1',0,'PBP',0,'LHDM_v1',0,'LHDM_v3_NNREG',0,'LHDM_v1_hard_thresholded',0,'OLS',0,'LH_v3',0,'LHDM_v3',0,'LHDM_v4',0,'LH_C',0,'LHDM_C',0);
        tmp_runsolver.ISAL1=1;

        [resnorm,true_nnz,x_nnz,true_xl1norm,x_l1norm, opt_dist, time_el, iter_hist] = ...
            runner(A,x,b,itest,NNEG,tmp_runsolver,solver_options);

        ISAL1resnorm_vec(i) = resnorm;
        ISAL1nnz_vec(i) = x_nnz;
        ISAL1time_vec(i) = time_el;
        ISAL1opt_dist_vec(i) = opt_dist;
        ISAL1x_l1norm_vec(i) = x_l1norm;
    end



    if runsolver.PBP == 1
        tmp_runsolver = struct('l1_magic',0,'l1_homotopy',0,'l1_homotopy_early',0,'SolveBP',0,'SolveOMP',0,'SPGL1',0,'SPGL1_tol',0,'YALL1',0,'ISAL1',0,'PBP',0,'LHDM_v1',0,'LHDM_v3_NNREG',0,'LHDM_v1_hard_thresholded',0,'OLS',0,'LH_v3',0,'LHDM_v3',0,'LHDM_v4',0,'LH_C',0,'LHDM_C',0);
        tmp_runsolver.PBP=1;

        [resnorm,true_nnz,x_nnz,true_xl1norm,x_l1norm, opt_dist, time_el, iter_hist] = ...
            runner(A,x,b,itest,NNEG,tmp_runsolver,solver_options);

        PBPresnorm_vec(i) = resnorm;
        PBPnnz_vec(i) = x_nnz;
        PBPtime_vec(i) = time_el;
        PBPopt_dist_vec(i) = opt_dist;
        PBPx_l1norm_vec(i) = x_l1norm;
    end



    if runsolver.SPGL1 == 1
        tmp_runsolver = struct('l1_magic',0,'l1_homotopy',0,'l1_homotopy_early',0,'SolveBP',0,'SolveOMP',0,'SPGL1',0,'SPGL1_tol',0,'YALL1',0,'ISAL1',0,'PBP',0,'LHDM_v1',0,'LHDM_v3_NNREG',0,'LHDM_v1_hard_thresholded',0,'OLS',0,'LH_v3',0,'LHDM_v3',0,'LHDM_v4',0,'LH_C',0,'LHDM_C',0);
        tmp_runsolver.SPGL1=1;

        [resnorm,true_nnz,x_nnz,true_xl1norm,x_l1norm, opt_dist, time_el, iter_hist] = ...
            runner(A,x,b,itest,NNEG,tmp_runsolver,solver_options);

        SPGL1resnorm_vec(i) = resnorm;
        SPGL1nnz_vec(i) = x_nnz;
        SPGL1time_vec(i) = time_el;
        SPGL1opt_dist_vec(i) = opt_dist;
        SPGL1x_l1norm_vec(i) = x_l1norm;
    end



    if runsolver.SolveBP == 1
        tmp_runsolver = struct('l1_magic',0,'l1_homotopy',0,'l1_homotopy_early',0,'SolveBP',0,'SolveOMP',0,'SPGL1',0,'SPGL1_tol',0,'YALL1',0,'ISAL1',0,'PBP',0,'LHDM_v1',0,'LHDM_v3_NNREG',0,'LHDM_v1_hard_thresholded',0,'OLS',0,'LH_v3',0,'LHDM_v3',0,'LHDM_v4',0,'LH_C',0,'LHDM_C',0);
        tmp_runsolver.SolveBP=1;

        [resnorm,true_nnz,x_nnz,true_xl1norm,x_l1norm, opt_dist, time_el, iter_hist] = ...
            runner(A,x,b,itest,NNEG,tmp_runsolver,solver_options);

        SolveBPresnorm_vec(i) = resnorm;
        SolveBPnnz_vec(i) = x_nnz;
        SolveBPtime_vec(i) = time_el;
        SolveBPopt_dist_vec(i) = opt_dist;
        SolveBPx_l1norm_vec(i) = x_l1norm;
    end



    if runsolver.SolveOMP == 1
        tmp_runsolver = struct('l1_magic',0,'l1_homotopy',0,'l1_homotopy_early',0,'SolveBP',0,'SolveOMP',0,'SPGL1',0,'SPGL1_tol',0,'YALL1',0,'ISAL1',0,'PBP',0,'LHDM_v1',0,'LHDM_v3_NNREG',0,'LHDM_v1_hard_thresholded',0,'OLS',0,'LH_v3',0,'LHDM_v3',0,'LHDM_v4',0,'LH_C',0,'LHDM_C',0);
        tmp_runsolver.SolveOMP=1;

        [resnorm,true_nnz,x_nnz,true_xl1norm,x_l1norm, opt_dist, time_el, iter_hist] = ...
            runner(A,x,b,itest,NNEG,tmp_runsolver,solver_options);

        SolveOMPresnorm_vec(i) = resnorm;
        SolveOMPnnz_vec(i) = x_nnz;
        SolveOMPtime_vec(i) = time_el;
        SolveOMPopt_dist_vec(i) = opt_dist;
        SolveOMPx_l1norm_vec(i) = x_l1norm;
    end



    if runsolver.l1_magic == 1
        tmp_runsolver = struct('l1_magic',0,'l1_homotopy',0,'l1_homotopy_early',0,'SolveBP',0,'SolveOMP',0,'SPGL1',0,'SPGL1_tol',0,'YALL1',0,'ISAL1',0,'PBP',0,'LHDM_v1',0,'LHDM_v3_NNREG',0,'LHDM_v1_hard_thresholded',0,'OLS',0,'LH_v3',0,'LHDM_v3',0,'LHDM_v4',0,'LH_C',0,'LHDM_C',0);
        tmp_runsolver.l1_magic=1;

        [resnorm,true_nnz,x_nnz,true_xl1norm,x_l1norm, opt_dist, time_el, iter_hist] = ...
            runner(A,x,b,itest,NNEG,tmp_runsolver,solver_options);

        l1_magicresnorm_vec(i) = resnorm;
        l1_magicnnz_vec(i) = x_nnz;
        l1_magictime_vec(i) = time_el;
        l1_magicopt_dist_vec(i) = opt_dist;
        l1_magicx_l1norm_vec(i) = x_l1norm;
    end



    if runsolver.l1_homotopy_early == 1
        tmp_runsolver = struct('l1_magic',0,'l1_homotopy',0,'l1_homotopy_early',0,'SolveBP',0,'SolveOMP',0,'SPGL1',0,'SPGL1_tol',0,'YALL1',0,'ISAL1',0,'PBP',0,'LHDM_v1',0,'LHDM_v3_NNREG',0,'LHDM_v1_hard_thresholded',0,'OLS',0,'LH_v3',0,'LHDM_v3',0,'LHDM_v4',0,'LH_C',0,'LHDM_C',0);
        tmp_runsolver.l1_homotopy_early=1;

        [resnorm,true_nnz,x_nnz,true_xl1norm,x_l1norm, opt_dist, time_el, iter_hist] = ...
            runner(A,x,b,itest,NNEG,tmp_runsolver,solver_options);

        l1_homotopy_earlyresnorm_vec(i) = resnorm;
        l1_homotopy_earlynnz_vec(i) = x_nnz;
        l1_homotopy_earlytime_vec(i) = time_el;
        l1_homotopy_earlyopt_dist_vec(i) = opt_dist;
        l1_homotopy_earlyx_l1norm_vec(i) = x_l1norm;
    end


    if runsolver.l1_homotopy == 1
        tmp_runsolver = struct('l1_magic',0,'l1_homotopy',0,'l1_homotopy_early',0,'SolveBP',0,'SolveOMP',0,'SPGL1',0,'SPGL1_tol',0,'YALL1',0,'ISAL1',0,'PBP',0,'LHDM_v1',0,'LHDM_v3_NNREG',0,'LHDM_v1_hard_thresholded',0,'OLS',0,'LH_v3',0,'LHDM_v3',0,'LHDM_v4',0,'LH_C',0,'LHDM_C',0);
        tmp_runsolver.l1_homotopy=1;

        [resnorm,true_nnz,x_nnz,true_xl1norm,x_l1norm, opt_dist, time_el, iter_hist] = ...
            runner(A,x,b,itest,NNEG,tmp_runsolver,solver_options);

        l1_homotopyresnorm_vec(i) = resnorm;
        l1_homotopynnz_vec(i) = x_nnz;
        l1_homotopytime_vec(i) = time_el;
        l1_homotopyopt_dist_vec(i) = opt_dist;
        l1_homotopyx_l1norm_vec(i) = x_l1norm;
    end


    if runsolver.SPGL1_tol == 1
        tmp_runsolver = struct('l1_magic',0,'l1_homotopy',0,'l1_homotopy_early',0,'SolveBP',0,'SolveOMP',0,'SPGL1',0,'SPGL1_tol',0,'YALL1',0,'ISAL1',0,'PBP',0,'LHDM_v1',0,'LHDM_v3_NNREG',0,'LHDM_v1_hard_thresholded',0,'OLS',0,'LH_v3',0,'LHDM_v3',0,'LHDM_v4',0,'LH_C',0,'LHDM_C',0);
        tmp_runsolver.SPGL1_tol=1;

        [resnorm,true_nnz,x_nnz,true_xl1norm,x_l1norm, opt_dist, time_el, iter_hist] = ...
            runner(A,x,b,itest,NNEG,tmp_runsolver,solver_options);

        SPGL1_tolresnorm_vec(i) = resnorm;
        SPGL1_tolnnz_vec(i) = x_nnz;
        SPGL1_toltime_vec(i) = time_el;
        SPGL1_tolopt_dist_vec(i) = opt_dist;
        SPGL1_tolx_l1norm_vec(i) = x_l1norm;
    end


    if runsolver.LHDM == 1
        tmp_runsolver = struct('l1_magic',0,'l1_homotopy',0,'l1_homotopy_early',0,'SolveBP',0,'SolveOMP',0,'SPGL1',0,'SPGL1_tol',0,'YALL1',0,'ISAL1',0,'PBP',0,'LHDM_v1',0,'LHDM_v3_NNREG',0,'LHDM_v1_hard_thresholded',0,'OLS',0,'LH_v3',0,'LHDM_v3',0,'LHDM_v4',0,'LH_C',0,'LHDM_C',0);
        tmp_runsolver.LHDM=1;
        iter_cell_LHDM_v4 = {};

        [resnorm,true_nnz,x_nnz,true_xl1norm,x_l1norm, opt_dist, time_el, iter_hist] = ...
            runner(A,x,b,itest,NNEG,tmp_runsolver,solver_options);

        LHDM_v4resnorm_vec(i) = resnorm;
        LHDM_v4nnz_vec(i) = x_nnz;
        LHDM_v4time_vec(i) = time_el;
        LHDM_v4opt_dist_vec(i) = opt_dist;
        LHDM_v4x_l1norm_vec(i) = x_l1norm;

        iter_cell_LHDM_v4{end+1} = iter_hist;
    end

end

if runsolver.l1_magic == 1, 
    nnz_vec=l1_magicnnz_vec; time_vec=l1_magictime_vec; opt_dist_vec=l1_magicopt_dist_vec; resnorm_vec=l1_magicresnorm_vec; x_l1norm_vec=l1_magicx_l1norm_vec;
    save('data_l1_magic.mat', 'nnz_vec', 'time_vec', 'opt_dist_vec', 'resnorm_vec', 'x_l1norm_vec'); 
end
if runsolver.l1_homotopy == 1, 
    nnz_vec=l1_homotopynnz_vec; time_vec=l1_homotopytime_vec; opt_dist_vec=l1_homotopyopt_dist_vec; resnorm_vec=l1_homotopyresnorm_vec; x_l1norm_vec=l1_homotopyx_l1norm_vec;
    save('data_l1_homotopy.mat', 'nnz_vec', 'time_vec', 'opt_dist_vec', 'resnorm_vec', 'x_l1norm_vec'); 
end
if runsolver.l1_homotopy_early == 1, 
    nnz_vec=l1_homotopy_earlynnz_vec; time_vec=l1_homotopy_earlytime_vec; opt_dist_vec=l1_homotopy_earlyopt_dist_vec; resnorm_vec=l1_homotopy_earlyresnorm_vec; x_l1norm_vec=l1_homotopy_earlyx_l1norm_vec;
    save('data_l1_homotopy_early.mat', 'nnz_vec', 'time_vec', 'opt_dist_vec', 'resnorm_vec', 'x_l1norm_vec'); 
end
if runsolver.SolveBP == 1, 
    nnz_vec=SolveBPnnz_vec; time_vec=SolveBPtime_vec; opt_dist_vec=SolveBPopt_dist_vec; resnorm_vec=SolveBPresnorm_vec; x_l1norm_vec=SolveBPx_l1norm_vec;
    save('data_SolveBP.mat', 'nnz_vec', 'time_vec', 'opt_dist_vec', 'resnorm_vec', 'x_l1norm_vec'); 
end
if runsolver.SolveOMP == 1, 
    nnz_vec=SolveOMPnnz_vec; time_vec=SolveOMPtime_vec; opt_dist_vec=SolveOMPopt_dist_vec; resnorm_vec=SolveOMPresnorm_vec; x_l1norm_vec=SolveOMPx_l1norm_vec;
    save('data_SolveOMP.mat', 'nnz_vec', 'time_vec', 'opt_dist_vec', 'resnorm_vec', 'x_l1norm_vec'); 
end
if runsolver.SPGL1 == 1, 
    nnz_vec=SPGL1nnz_vec; time_vec=SPGL1time_vec; opt_dist_vec=SPGL1opt_dist_vec; resnorm_vec=SPGL1resnorm_vec; x_l1norm_vec=SPGL1x_l1norm_vec;
    save('data_SPGL1.mat', 'nnz_vec', 'time_vec', 'opt_dist_vec', 'resnorm_vec', 'x_l1norm_vec'); 
end
if runsolver.SPGL1_tol==1, 
    nnz_vec=SPGL1_tolnnz_vec; time_vec=SPGL1_toltime_vec; opt_dist_vec=SPGL1_tolopt_dist_vec; resnorm_vec=SPGL1_tolresnorm_vec; x_l1norm_vec=SPGL1_tolx_l1norm_vec;
    save('data_spgl1_tol.mat', 'nnz_vec', 'time_vec', 'opt_dist_vec', 'resnorm_vec', 'x_l1norm_vec'); 
end
if runsolver.ISAL1 == 1, 
    nnz_vec=ISAL1nnz_vec; time_vec=ISAL1time_vec; opt_dist_vec=ISAL1opt_dist_vec; resnorm_vec=ISAL1resnorm_vec; x_l1norm_vec=ISAL1x_l1norm_vec;
    save('data_ISAL1.mat', 'nnz_vec', 'time_vec', 'opt_dist_vec', 'resnorm_vec', 'x_l1norm_vec'); 
end
if runsolver.PBP == 1, 
    nnz_vec=PBPnnz_vec; time_vec=PBPtime_vec; opt_dist_vec=PBPopt_dist_vec; resnorm_vec=PBPresnorm_vec; x_l1norm_vec=PBPx_l1norm_vec;
    save('data_PBP.mat', 'nnz_vec', 'time_vec', 'opt_dist_vec', 'resnorm_vec', 'x_l1norm_vec'); 
end
if runsolver.LHDM==1, 
    nnz_vec=LHDM_v4nnz_vec; time_vec=LHDM_v4time_vec; opt_dist_vec=LHDM_v4opt_dist_vec; resnorm_vec=LHDM_v4resnorm_vec; x_l1norm_vec=LHDM_v4x_l1norm_vec;
    save('data_LHDM_v4.mat', 'nnz_vec', 'time_vec', 'opt_dist_vec', 'resnorm_vec', 'x_l1norm_vec', 'iter_cell_LHDM_v4'); 
end

fprintf("time = %f \n", time_el);
fprintf("resnorm = %f \n", resnorm);
fprintf("nnz = %d \n", x_nnz);
fprintf("opt_dist = %f \n", opt_dist);
fprintf("x_l1norm = %f \n", x_l1norm);
fprintf("true_xl1norm = %f \n", true_xl1norm);

