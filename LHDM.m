function [x,resnorm,exitflag, outeriter, iter_vec] = LHDM_v4(C,d,options,verbose) 

% LHDM solves underdetermined linear least squares with nonnegativity constraints
% using Lawson-Hanson with Deviation Maximization algorithm, described in 
%
% M.Dessole, M.Dall'Orto, F.Marcuzzi, 
% "The Lawson-Hanson Algorithm with Deviation Maximization: Finite Convergence and Sparse Recovery", 2022
%
% INPUT
% C: underdetermined matrix
% d: column array of set points coordinates
% options: structure containing the values of optimization parameters:
%    init       true if ULS initialization of Passive set is desired, false otherwise
%    tol        tolerance on the projected residual, stop criterion
%    itmax      maximum number of iterations 
%    k          n. of columns considered by the DM
%    positrick  solve an arbitrary signed problem with the 'positivity trick'
%    thres_cos  threshold on cosine between new columns 
%    thres_nrm  threshold on relative length of columns
%    thres_w    threshold on dual variables
%    NNEG       look for a nonnegative solution
%    tol        tolerance on dual variables positivity
%    ensure_descend_direction  ensure that DM finds a descend direction
%
% OUTPUT
% x: sparse vector that minimizes NORM(d-C*x)
% resnorm: squared 2-norm of the residual: norm(d-C*X)^2
% exitflag: exit condition, possible value are:
%    1  LH converged with a solution X
%    0  Iteration count was exceeded, increasing the tolerance
%       (OPTIONS.Tol) may lead to a solution

% Dates:
% Written by M. Dessole, F. Marcuzzi - April 2022

if nargin < 2
    error('MATLAB:LHDM:NotEnoughInputs',...
        getString(message('MATLAB:LHDM:NotEnoughInputs')));
end

positrick = options.positrick; 
ensure_descend_direction = options.ensure_descend_direction; 

[m,n] = size(C);
nmezzi = n/2;

C_colnorms = vecnorm(C,2,1);

ls_opts.UT = true;

% Initialize vector of n zeros and Infs (to be used later)
nZeros = zeros(n,1);
if positrick
    wz = zeros(nmezzi,1);
else
    wz = zeros(n,1);
end

% Initialize set of non-active columns to null
if positrick
    P = false(nmezzi,1);
    Cr = C(:,1:nmezzi);
else
    P = false(n,1);
end
cardP = 0;
IP = [];
jmax_valid_VH = 0;
% Initialize set of active columns to all and the initial point to zeros
if positrick
    Z = true(nmezzi,1);
else
    Z = true(n,1);
end
x = nZeros;

tmptime = 0;

% Check if options was created with optimoptions
if ~isempty(options) 
    if ~isa(options,'struct')
    error('MATLAB:LDHM:ArgNotStruct',...
        getString(message('MATLAB:LHDM:commonMessages:ArgNotStruct')));
    end
    if isfield(options,'itmax')
        itmax = options.itmax;
    else
        itmax = 2*m;
    end
    
    if isfield(options,'thres_cos')
        thres_cos = options.thres_cos;
        if (thres_cos <= eps)
            LHDMflag = 0;
        end
    else
        thres_cos = 0.2;
    end
    if isfield(options,'thres_w')
        thres_w = options.thres_w;
    else
        thres_w = 0.8;
    end
    if isfield(options,'thres_nrm')
        thres_nrm = options.thres_nrm;
    else
        thres_nrm = 0.05;
    end
    if isfield(options,'k')
        k = options.k;
        if (k == 1)
            LHDMflag = 0;
        else
            LHDMflag = 1;
        end
    else
        k = ceil(m/20); 
    end
    if isfield(options,'ncMax')
        ncMax = options.ncMax;
    else
        ncMax = k;
    end
    if isfield(options,'tol')
        tol = options.tol;
        %tol = options.tol * norm(d);
    else
        %deprecated (for its cost) tol = 10*eps*norm(C,1)*m*n;
        tol = 1.e-5;
    end
    if isfield(options,'init') 
        if options.init
            xtmp = C\d;
            Idx = find(xtmp>0);    
            % Initialize set of non-active columns 
            P(Idx) = true; 
            % Initialize set of active columns 
            Z(Idx) = false;
            % Initialize starting point
            x(P) = C(:,P)\d;
            tmp = find(x<0);
            if(size(tmp,1)>0)
               x(tmp) = 0; 
               P(tmp) = false;
               Z(tmp) = true;
            end
            cardP = sum(P);
            resid = d - C*x; %disp(['norm(resid) = ',num2str(norm(resid))]);
        else
            resid = d; % x==0! - C*x;
        end   
    end
else
    thres_cos   = 0.2222;
    thres_w = 0.8;
    thres_nrm = 0.05;
    k = ceil(m/20);
    ncMax = k;
    itmax = 3*n;
    tol = 1.e-5; %deprecated (for its cost)  tol = 10*eps*norm(C,1)*m*n;
    LHDMflag = 1;
    resid = d; % x==0! - C*x;
end

old_resnorm = Inf;

if positrick
    wa = Cr'*resid; %w = [wa; -wa];
    w = abs(wa); 
else
    w = C'*resid;
end

% Set up iteration criterion
outeriter = 0;
totiter = 0;
iter = 0;
iter_vec = zeros(itmax,1);
flag_lin_dep = 0;
Q = eye(m);
R = zeros(m);
iVHG = 0;
VHG = zeros(m,2*m);
whichVHG = zeros(1,2*m); % 0=H, 1=G
fromVHG = zeros(1,2*m); % indice da cui parte "v_k"
n2VHG = zeros(1,2*m); % norm-2 di "v_k"
iVH = 0;
VH = zeros(m);
W = zeros(m,m);
invalid_Gcoeff = 99;
VG = zeros(3,2*m);
iVG = 0;
QTd = d; % Q'd
MQTd = zeros(m,m);
exclude = [];
after_Givens = 0;

old_t = 0;

k0 = k;

%disp('-----------------------------------------')
% Outer loop to put variables into set to hold positive coefficients
ok_descend = 1;
go_on = 1;
while ( go_on && (totiter < itmax) )  
    %increase total iteration counter
    totiter = totiter+1;
    % Create wz, a Lagrange multiplier vector of variables in the zero set.
    % wz must have the same size as w to preserve the correct indices, so
    % set multipliers to -Inf for variables outside of the zero set.
    wz(P) = -Inf;
    wz(Z) = w(Z);
    wz(exclude) = -Inf;

    t = [];
    % Deviation Maximization    
    [wzI, tmpIDM] = sort(wz, 'descend');
    % initialize set of selected indices with the index selected by standard LH:
    thres_wloc = thres_w*wzI(1);
    % candidates indices
    ncDM_0 = min([length(find(wzI>thres_wloc)),ncMax,m-cardP]);
    ncDM = ncDM_0;
    % restrict number of candidates if necessary
    IDM = tmpIDM(1:ncDM)';
    if positrick
        jj = 0;
        for jc=tmpIDM(1:ncDM)'
            jj = jj + 1;
            if wa(jc) < 0
                IDM(jj) = jc + nmezzi;
            else
                IDM(jj) = jc;
            end
        end
    end
    % number of indices added to P
   
    %increase outer iteration counter
    outeriter = outeriter + 1;
   
    % Reset intermediate solution z
    z = nZeros; 

    % Compute intermediate solution using only variables in positive set
    j = cardP + 1;
    R2add = C(:,IDM);
    for jj = 1:iVHG
       if whichVHG(jj)==0
           iH = fromVHG(jj);
           v_k = VHG(iH:m,jj);
           n2v_k = n2VHG(jj);
           tmpr = 2/n2v_k/n2v_k;
           R2add(iH:m,1:ncDM) = R2add(iH:m,1:ncDM) - tmpr*v_k*(v_k'*R2add(iH:m,1:ncDM)); %tmp(iH:m,j:j2);
       else   
           cG = VHG(1,jj); sG = VHG(2,jj); tmpj = VHG(3,jj);
           tmpr = cG*R2add(tmpj,1:ncDM) -sG*R2add(tmpj+1,1:ncDM);
           R2add(tmpj+1,1:ncDM) = sG*R2add(tmpj,1:ncDM) +cG*R2add(tmpj+1,1:ncDM);
           R2add(tmpj,1:ncDM) = tmpr;
       end
    end

    vnorm_tmpc = vecnorm(R2add(j:m,1:ncDM),2,1); %disp(['vnorm_tmpc = ',mat2str(vnorm_tmpc')]);
    maxnorm_tmpc = 0;
    thres_nrmloc = 0;
    viIDMok = [];
    iIDM = 0;
    iH = 1;
    while iIDM < ncDM
       iIDM = iIDM + 1;
       ir = cardP + iH;
       n2iH = vnorm_tmpc(iH);
       if n2iH >= thres_nrmloc
           tmp_s1 = max(abs(R2add(j:m,iH)'*R2add(j:m,1:iH-1))./vnorm_tmpc(iH)./vnorm_tmpc(viIDMok));
           if iH==1 || (tmp_s1 < thres_cos)
               signfc = sign(R2add(ir,iH)); if signfc==0, signfc = 1; end; 
               v_k = [signfc.*norm(R2add(ir:m,iH)) - R2add(ir,iH); - R2add(ir+1:m,iH)]; %disp(['v_k = ',mat2str(v_k')]);
               n2v_k = norm(v_k);  %disp(['n2v_k = ',num2str(n2v_k)]);
               if n2v_k == 0, n2v_k = 1; end
               tmpr = 2/n2v_k/n2v_k;
               R2add(ir:m,iH:ncDM) = R2add(ir:m,iH:ncDM) - tmpr*v_k*(v_k'*R2add(ir:m,iH:ncDM));
               if abs(R2add(ir,iH)) >= 1.e-15
                   if n2iH > maxnorm_tmpc, maxnorm_tmpc = n2iH; thres_nrmloc = thres_nrm * maxnorm_tmpc; end
                   QTd(ir:m) = QTd(ir:m) - tmpr*v_k*(v_k'*QTd(ir:m));
                   iVHG = iVHG + 1;
                   VHG(ir:m,iVHG) = v_k;
                   whichVHG(iVHG) = 0;
                   fromVHG(iVHG) = ir;
                   n2VHG(iVHG) = n2v_k;
                   viIDMok = [viIDMok iIDM];
                   iH = iH + 1;
               else
                   exclude = [exclude IDM(iIDM)];
               end
           else
               R2add = [R2add(:,1:iH-1) R2add(:,iH+1:end)];
               ncDM = ncDM - 1;
           end
       else
           R2add = [R2add(:,1:iH-1) R2add(:,iH+1:end)];
           ncDM = ncDM - 1;
       end
    end
    R(:,j:cardP+ncDM) = R2add;
    t = IDM(viIDMok);
    addedP = length(t);
    if positrick
       tr = t; tmpI = find(tr > nmezzi); tr(tmpI) = tr(tmpI) - nmezzi; 
    end

   % Move variable t from zero set to positive set
   IP = [IP t];
   old_cardP = cardP;
   cardP = cardP + addedP;
   iter_vec(iter+outeriter) = cardP;
   
   z(IP) = linsolve(R(1:cardP,1:cardP),QTd(1:cardP),ls_opts);
   if ensure_descend_direction
       if positrick
           while sum(z(tr).*w(tr)) < 0 
               cardP = cardP - 1;
               IP = IP(1:end-1);
               z(IP) = linsolve(R(1:cardP,1:cardP),QTd(1:cardP),ls_opts);
               t = t(1:end-1);
               tr = t; tmpI = find(tr > nmezzi); tr(tmpI) = tr(tmpI) - nmezzi; 
           end
       else
           while (sum(z(t).*w(t)) < 0)
               cardP = cardP - 1;
               IP = IP(1:end-1);
               z(IP) = linsolve(R(1:cardP,1:cardP),QTd(1:cardP),ls_opts);
               t = t(1:end-1);
           end
       end   
       if length(t) == 0, ok_descend = 0; end
   end

   if positrick
       P(tr) = true; 
       Z(tr) = false;
   else
       P(t) = true;
       Z(t) = false;
   end
   
   removedP = 0;
   t_removed = [];
   jmax_valid_VH = cardP;
   if positrick
       for ip=1:length(IP)
           jj = IP(ip);
           if z(jj) < 0
               if jj <= nmezzi
                   z(jj+nmezzi) = -z(jj); z(jj) = 0;
                   IP(ip) = IP(ip) + nmezzi;
               else
                   z(jj-nmezzi) = -z(jj); z(jj) = 0;
                   IP(ip) = IP(ip) - nmezzi;
               end
               R(:,ip) = -R(:,ip);
           end
       end
   else
       % inner loop to remove elements from the positive set which no longer belong
       while (any(z(P) <= 0.) && (totiter < itmax))
           totiter = totiter +1;
           iter = iter+1;
           % Find indices where intermediate solution z is approximately negative
           IPneg = (z <= 0.) & P;
           % Choose new x subject to keeping new x nonnegative
           b = x(IPneg)./(x(IPneg) - z(IPneg));
           alpha = min(b);
           x = x + alpha*(z - x);
           % number of indices removed from P
           tr = find((x <= 0.) & P)'; 
           if (length(tr) > 0)
               for t_j = tr
                   if t_j == t(1)
                       continue
                   end
                   t_removed = [t_removed, t_j];
                   it = find(IP == t_j); 
                   removedP = removedP + 1; %length(tr);
                   % left circular shift of R columns from "it":
                   if 1
                       tmpc = R(:,it); 
                       R(:,it:cardP) = [R(:,it+1:cardP) tmpc];
                   else
                       tmpc = R(1:cardP,it); 
                       R(1:cardP,it:cardP) = [R(1:cardP,it+1:cardP) tmpc];
                   end
                   % ri-triangularize R:
                   for j = it:cardP-1
                       [cG,sG] = Givens_matrix(R(j,j),R(j+1,j));
                       iVHG = iVHG + 1;
                       VHG(1,iVHG)=cG; VHG(2,iVHG)=sG; VHG(3,iVHG)=j;
                       whichVHG(iVHG) = 1;
                       tmpr = cG*R(j,j:cardP) -sG*R(j+1,j:cardP);
                       R(j+1,j:cardP) = sG*R(j,j:cardP) +cG*R(j+1,j:cardP);
                       R(j,j:cardP) = tmpr;
                       tmpd = cG*QTd(j) - sG*QTd(j+1);
                       QTd(j+1) = sG*QTd(j) + cG*QTd(j+1);
                       QTd(j) = tmpd;
                   end
                   IP = [IP(1:it-1) IP(it+1:cardP)];
                   cardP = cardP - 1;
               end
           end
           % Reset Z and P given intermediate values of x
           Z = ((x <= 0.) & P) | Z;
           P = ~Z; 
           iter_vec(iter+outeriter) = cardP;
           % Reset z 
           z = nZeros;
           % Compute intermediate solution using only variables in positive set
           z(IP) = linsolve(R(1:cardP,1:cardP),QTd(1:cardP),ls_opts);
       end
       % update cardinality of P
   end

   x=z;

   resid = d - C(:,IP)*x(IP); 
    if positrick
        wa = Cr'*resid; %w = [wa; -wa];
        w = abs(wa); 
    else
        w = C'*resid;
    end

   if removedP > 0
       w(t_removed) = -inf;
   end

   % break if the algorithm is stagnating
   if 0
       if ( (iter+outeriter>10) && ( abs(iter_vec(iter+outeriter) - iter_vec(iter+outeriter-10) )/10 < 1e-8 ))
           break;
       end
   end
   
   if positrick
       go_on = any(Z) && any(w(Z) > tol) && ok_descend;
   else
       go_on = any(Z) && (any(w(Z) > tol) || any(x(P) <= 0.)) && ok_descend;
   end
end

if (outeriter < itmax)
    exitflag = 1;
else
    exitflag = 0;
end

iter_vec = iter_vec(1:iter+outeriter);

resnorm = resid'*resid;

if totiter == itmax, 
    disp("error totiter too small!"); 
    if itmax < 500, keyboard, end
end

end



function [c,s] = Givens_matrix(xp,xq)

if xq==0
  c=1; s=0;
else
  if abs(xq)>abs(xp)
    r=xp/xq; s=-1/sqrt(1+r^2); c=-s*r;
  else
    r=xq/xp; c=1/sqrt(1+r^2); s=-c*r;
  end;
end;

end
