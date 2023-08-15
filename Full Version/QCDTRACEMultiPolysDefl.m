

%  Version 2 tries to automate how many noises for each part.

%  QCDTRACEMultiPolys uses a sequence of poly's. 

%  This program is to compute the Trace of the inverse of a QCD matrix.  It
%  first finds a poly, then used many PP-GMRES's to compute Tr(A^{-1}-p(A)), then
%  does Tr(p(A)) stochastically.  2021
%
%  For QCD, need to adjust to deal with complex harm Ritz val's not in
%  pairs !

%  Poly Precondioned GMRES - 2017 version
%   
%
clear all;
format long e

global A
rand('seed',0)

% tic

% load qcd4; A = qcd;  n = size(A,1);  A = A + 3e-3*speye(n,n);   
load h01_4444; n = size(B,1); spI = speye(n,n); A = spI - 0.1570*B;
% load h01_8888; n = size(B,1); spI = speye(n,n); A = spI - 0.1570*B;
% load h01_12121216; n = size(B,1); spI = speye(n,n); A = spI - 0.1570*B;
% load h01_16161624; n = size(B,1); spI = speye(n,n); A = spI - 0.1570*B;
% load h02_12121216; n = size(B,1); spI = speye(n,n); A = spI - 0.1570*B;
clear B
% load matmarket888000; A = qcd;
%                       n = 49152;
%                       spidentity = spdiags(ones(n,1),0,n,n);
%                       A = A + 6.4*spidentity;
% load h01; n = size(B,1); spI = speye(n,n); A = spI - 0.1570*B;
% n = size(A,1);
%
% n = 49152  % (for qcd matrix matmarket888000)
%                                           A = A - .3*eye(n,n);


tic
% j

trtol = 0.005*(4^4)
% trtol = 0.002*(8^4)  %5.0e-1
% trtol = 0.001*(12^3*16) 
npolys = 3;
ndefl = 16
deflate = 1;  % deflate = 1 means: yes, deflate
deg = [200 50 7]'  % Degree of poly  (degree of alpha*p(alpha))
% deg = [800 100 7]'  
% deg = [ 700 300 7]' 
% nnva = [2 100 100 100]'  % Number Z4 noise vects for each part (there are npolys+1 parts)
% nnva = [5000 40000 80000 120000 2]'
% nnva = [ 40000 120000 120000 2]'
nnva = [ 20 1000 3000 2]'
freqcherr = 5
freqcherr2 = 5
thqa = zeros(max(deg)+1,npolys);

m = 50 % Dimension of Krylov subspace
mperm = m;


% k = 0    % Number of approximate eigenvectors saved at the restart
%numev = 1
rtol = 1.e-8;
cyclim = 100
istab = 1;

% kperm = k;
%-------------- RHS vectors.--------------------------%
cpu0 = cputime
% randn('seed',0)
rng('default')
% rng(3)
% b = randn(n,1);
b = randn(n,1); b = b/norm(b);
%                 b = ones(n,1) + sqrt(-1)*ones(n,1); b = b/norm(b);
                                     
% % related rhs's:
% epsi = 1.e-4
% temp = randn(n,nrhs);
% b(:,1) = temp(:,1);
% for i=2:nrhs
%     b(:,i) = b(:,1) + epsi*temp(:,i);
% end
% %
% load burakrhs;
% b(:,1) = burakrhs;
% b(:,2) = b(:,1);
% b = ones(n,nrhs);
%-------------------------------------------------------%
%  Begin with generating a polynomial for the poly preconditioning.
%  First find coeff's of p.
%
x = zeros(n,1);
cycle = 1;
mvp = 0;
iter = 0;
vops = 0;
dps = 0;

%%%%%%%%%%%%%%%%%%%%% 
% Loop to generate all the polys

ipoly = 1;
pp = deg(ipoly,1) + 1;

j = 1;
rninit = norm(b(:,1));
rn = rninit;
vops = vops + 1;
dps = dps + 1;
% v(:,1) = b(:,1)/rninit;
% %
% for j=1:p 
%     v(:,j+1) = op(n,v(:,j));
%     mvp = mvp + 1;
% end
% lsmat = v(:,2:p+1)'*v(:,2:p+1);
%     vops = vops + p*p;
% cls = v(:,2:p+1)'*v(:,1);
%     vops = vops + p;
% dls = lsmat\cls;
%                                         dls
%                                         co = dls;
                                        
%  Now find roots of q.
v(:,1) = b(:,1)/rn;
c(1,1) = rn;
%                             % For damped poly:
%                             temp = A*b;
%                             tn = norm(temp);
%                             v(:,1) = temp/tn;
%                             c(1,1) = tn;
c(2:pp+1,1) = zeros(pp,1);
j = 1;
while ( j <= pp ) 
    f = op(n,v(:,j));
    iter = iter + 1;
    mvp = mvp + 1;
    vnf = norm(f);
    for i = 1:j
      h(i,j) = v(:,i)'*f;
      f = f - h(i,j) * v(:,i);
    end  %for
    vn = norm(f);
            vops = vops + 2*j+1;
            dps = dps + j+1;
%--------------------------------------------------%
%--------reorthogonalization section-------%
    if( vn < 0*vnf )
         % disp( 'did a reorthog')
      for i = 1:j
        dot = v(:,i)'*f;
        f = f - dot * v(:,i);
        h(i,j) = h(i,j) + dot;
      end  %for
      vn = norm(f);
             vops = vops + 2*j+1;
    end  %if
%--------------------------------------------------%
    h(j+1,j) = vn;
    v(:,j+1) = f/h(j+1,j);
        vops = vops + 1;
    j = j + 1;
end  %while

% hh = h(1:p,1:p);
%  ----Set up and solve linear equations.-----%
d = h \ c;

srv = c-h*d;   
srvnorm = norm(srv)
% x(:,1) = x(:,1) + v(:,1:pp)*d;
%        vops = vops + pp;
% %

for ipoly = 1:npolys
ipoly
pp = deg(ipoly,1) + 1

  hh = h(1:pp,1:pp);             
  em = zeros(pp,1);
  em(pp) = 1;
  ff = hh' \ em;
%
%FOM section
%Comment next line for FOM instead of GMRES:
          hh(:,pp) = hh(:,pp) + h(pp+1,pp)^2 * ff;
  [g,dd] = eig(hh); 
  ddd = diag(dd);
  dabs = abs(ddd);
  [thabs,ind] = sort(dabs);
  g = g(:,ind);
  th = ddd(ind);
     
  % Section for finding right and left eigenvectors for Deflation
  if(deflate == 1)
  if(ipoly == 1)
  nh = n/2;
  ii = sqrt(-1);
  for id = 1:2:ndefl
%       thdefl = th(1:id);
      thdefl = th(1:ndefl,1);
      zr(:,id) = v(:,1:pp)*g(:,id);
      zr(:,id) = zr(:,id)/norm(zr(:,id));
      rnar(id,1) = norm(A*zr(:,id) - th(id)*zr(:,id));
      tt = gamma5(zr(:,id),n,1);
      zl(:,id+1) = tt;
      rnal(id+1,1) = norm(A'*zl(:,id+1) - conj(th(id+1))*zl(:,id+1));
      
      zr(:,id+1) = v(:,1:pp)*g(:,id+1);
      zr(:,id+1) = zr(:,id+1)/norm(zr(:,id+1));
      rnar(id+1,1) = norm(A*zr(:,id+1) - th(id+1)*zr(:,id+1));
      tt = gamma5(zr(:,id+1),n,1);
      zl(:,id) = tt;
      rnal(id,1) = norm(A'*zl(:,id) - conj(th(id))*zl(:,id));
  end

% % This section is for first eval real
%       id = 1;
%       thdefl = th(1:ndefl,1);
%       zr(:,id) = v(:,1:pp)*g(:,id);
%       zr(:,id) = zr(:,id)/norm(zr(:,id));
%       rnar(id,1) = norm(A*zr(:,id) - th(id)*zr(:,id));
%       tt = Travisgamma5(zr(:,id),n,1);
%       zl(:,id) = tt;
%       rnal(id,1) = norm(A'*zl(:,id) - conj(th(id))*zl(:,id));
% for id = 2:2:ndefl
% %       thdefl = th(1:ndefl,1);
%       zr(:,id) = v(:,1:pp)*g(:,id);
%       zr(:,id) = zr(:,id)/norm(zr(:,id));
%       rnar(id,1) = norm(A*zr(:,id) - th(id)*zr(:,id));
%       tt = Travisgamma5(zr(:,id),n,1);
%       zl(:,id+1) = tt;
%       rnal(id+1,1) = norm(A'*zl(:,id+1) - conj(th(id+1))*zl(:,id+1));
%       
%       zr(:,id+1) = v(:,1:pp)*g(:,id+1);
%       zr(:,id+1) = zr(:,id+1)/norm(zr(:,id+1));
%       rnar(id+1,1) = norm(A*zr(:,id+1) - th(id+1)*zr(:,id+1));
%       tt = Travisgamma5(zr(:,id+1),n,1);
%       zl(:,id) = tt;
%       rnal(id,1) = norm(A'*zl(:,id) - conj(th(id))*zl(:,id));
%   end
                                thdefl 
                                rnar
                                rnal
  end
  end
  
  
%  Now the modified Leja ordering from Reichel.
%  Can make more efficient (according to Bai, Hu, Riechel paper, but they
%    don't give details there)    
                                            figure(7)
                                            plot(th,'r*')
    thr(1) = th(1);
    index(1) = 1;
    j = 2;
    for i=2:pp
        if( abs(th(i)) > abs(thr(1)) )
          thr(1) = th(i);
        end
    end
    iatendofloop = i;
    while( j <= pp )
        for i=1:pp
            pr(i) = 0;
            for ii=1:j-1
                pr(i) = pr(i) + log(abs(th(i)-thr(ii)));
            end
        end
        thr(j) = th(1);
        prj = pr(1);
        itemp = 1;
        for i=2:pp 
            if( pr(i) > prj )
                thr(j) = th(i);
                prj = pr(i);
            end
        end
        j = j + 1;
    end
   
    thq = transpose(thr); 
                    if( pp == 1)
                        thq = 1;
                    end
                                figure(8)
                                plot(thq,'*')
                                
% SECTION to find value of rest of poly at every root.  This tell us about
%       the conditioning.

        for j=1:pp
            prod = 1;
            for i=1:pp
                if(i~=j)
                    prod = prod*(1-thq(j)/thq(i));
                end
            end
            ropatr(j,1) = abs(prod);
        end
                                              ropatra(:,ipoly) = ropatr;     
% Set up an array with roots to add and merge it with the previous array.
        if(istab == 1)
        j = 1;
        iex = 0;
        pporig = pp;
        while (j <= pp)
            rootj = thq(j);
            extra = ceil((log10(ropatr(j,1)) - 4)/14);
            for ie = 1:extra
                iex = iex + 1;
                exa(iex,1) = rootj;
                place = j + ceil(ie*(pp-j)/extra);
                pla(iex,1) = place;
            end
            j = j + 1;
        end
        
        if(iex > 0 )
        pp = pp + iex;
        [pla,ind] = sort(pla);
        exa = exa(ind);
        
        j = 1;
        ir = 1;
        ie = 1;
        while( j <= pp ) 
            if( ir <= pla(ie) )
                rootj = thq(ir);
                tht(j) = rootj;
                j = j + 1;
                ir = ir + 1;
            else
                rootj = exa(ie);
                tht(j) = rootj;
                j = j + 1;
                ie = ie + 1;
            end
        end
        thq = tht';
        end
        end
            
cqzero = 1;
cq(1,1) =  - 1/thq(1);
for j=2:pp
%     num = -cqzero;
%     den = thq(j);
%     for i=1:j-1
%         num = num - cq(i,1)*den;             
%         den = den*(thq(j)-thq(i));
%     end
%     cq(j,1) = num/den;
    cq(j,1) = - cq(j-1,1)/thq(j); 
end
  
% Now for coeff's of p (Newton basis)
cp = - cq;

% % xt = oppppNewt(n,b,pp,cp,thq);
% xt = oppppNewtNew(n,b,pp,thq);
% r3 = b(:,1) - op(n,xt);
% r4 = oppprB(n,b,pp,thq);
% rd43 = norm(r4-r3)

thqa(1:pp,ipoly) = thq;
deg(ipoly,1) = pp - 1;

clear c
clear v
clear thq
clear thr
clear tht

end
clear h
toc
tic

% end of generating the poly
%                                     
                                        
%                                                                
%
%-------------------------------------------------------%
%

% Which poly to use for the POLY Prec? Here may not be highest degree poly.
ipoly = 1;
pp = deg(ipoly,1) + 1;
ipoly2 = 2;
pp2 = deg(ipoly2,1) + 1;
thq = thqa(:,ipoly);
thq2 = thqa(:,ipoly2);
nnv = nnva(1,1);
errsqtarget = trtol^2/(npolys+1);  errtarget = sqrt(errsqtarget)

in = 0;
errsq = 1.e+50;
while( (in < nnv) & (errsq > errsqtarget) )
for i=1:freqcherr
    in = in + 1;
    b = rand(n,1);
    b = z4noise2(n,b);
m = mperm;
% in
cycle = 1;   
x = zeros(n,1);
j = 1;          
rorig = b;
r = rorig;
rn = norm(rorig);  rninit = rn;
    mvp = mvp + 1;
    vops = vops + 1;
                                    
vn = rn;
v(:,1) = rorig/vn;
    vops = vops + 1;
c(1,1) = vn;
c(2:m+1,1) = zeros(m,1);
                                                    
while ( (rn/rninit > rtol) & (cycle <= cyclim) ) 
while ( (j <= m) & (rn/rninit > rtol) )
%     wv = op(n,v(:,j));
    f = oppprQCD(n,v(:,j),pp2,thq2);
    iter = iter + 1;
        mvp = mvp + pp2;
        vops = vops + pp2;
    vnf = norm(f);
    for i = 1:j
      h(i,j) = v(:,i)'*f;
      f = f - h(i,j) * v(:,i);
    end  %for
    vn = norm(f);
            vops = vops + 2*j+1;
            dps = dps + j+1;
%--------------------------------------------------%
%--------reorthogonalization section-------%
    if( vn < 0*vnf )
         % disp( 'did a reorthog')
      for i = 1:j
        dot = v(:,i)'*f;
        f = f - dot * v(:,i);
        h(i,j) = h(i,j) + dot;
      end  %for
      vn = norm(f);
             vops = vops + 2*j+1;
    end  %if
%--------------------------------------------------%
    h(j+1,j) = vn;
    v(:,j+1) = f/h(j+1,j);
        vops = vops + 1;
%subsection for res.norm at every iteration
	d = h(1:j+1,1:j) \ c(1:j+1);  
    srv = c(1:j+1)-h(1:j+1,1:j)*d;
    gdr(iter,1) = norm(srv);
    rn = gdr(iter,1);
%now subsection for res.norm if preconditioning used
 %     xw = x(:,1) + v(:,1:j)*d;
 %     wv = b(:,1) - op(n,xw);   %--may remove this and next line----%
 %     gdr(mvp) = norm(wv);
%end subsection
    j = j + 1;
end  %while
m = j-1;
hh = h(1:m,1:m);

  %----Set up and solve linear equations.-----%
%  Should already be done:  d = h \ c;
%
%FOM section.
%Next line is for FOM instead of GMRES:
%  d = h(1:m,1:m) \ c(1:m);
%End FOM section
%
% srv = c-h(1:m+1,1:m)*d;
cycle;
%                                                 xbef = x(:,1)
% rn = norm(srv)
rna(cycle,1) = rn;
x(:,1) = x(:,1) + v(:,1:m)*d;
    vops = vops + m;

r = v(:,1:m+1)*srv;
       vops = vops + m+1;
%                                                  srv
rnwithpp = norm(r);

  j = 1;
  cycle = cycle + 1;
  vn = norm(r);
    vops = vops + 1;
  v(:,1) = r/vn;
  c(1,1) = vn;
end  %while
% cpu1 = cputime - cpu0
rntrwoppdivbyrninit = rn/rninit;

xreal(:,1) = oppppNewtNewQCD(n,x(:,1),pp2,thq2); %  NEED TO check counts for this! (mvp,vops)
rreal = rorig(:,1) - op(n,xreal);
rrealn = norm(rreal);
x = xreal;
if( pp > 1)
    mvp = mvp + pp2;  % pp-1 for oppppNewtNew 
    vops = vops + 2*pp2; % 2pp-1 for oppppNewtNew 
else
    mvp = mvp + 1;
end

pb = oppppNewtNewQCD(n,b,pp,thq);
bAib(in,1) = b'*(x-pb);
%                                         norm(b-A*x)
%                                         norm(b-A*pb)
% mvp

end

mvp
tra = sum(bAib(1:in,1))/in;
var = 0;
for i=1:in
    diff = bAib(i,1) - tra;
    var = var + diff'*diff;
end
var = var/(in-1);  % check this
error = sqrt(var/in);
errsq = error^2 ;

in
currenterr = error

end
errsum = error;
numnoia(1,1) = in;
%
% %%%%%%%%%
% %TRACE and Variance
% tra = sum(bAib)/nnv;
% var = 0;
% for i=1:nnv
%     diff = bAib(i,1) - tra;
%     var = var + diff'*diff;
% end
% var = var/(nnv-1);  % check this
% error1 = sqrt(var/nnv)
% 
% errsum = error1;
% % stdev = sqrt(var)
% numnoia(1,1) = nnv;

%%%%%%%%%%%%%%%
%Phase 2, deal with traces of p_i(A) - p_i+1(A)

for ipoly=1:npolys-1
    
nnv2 = nnva(ipoly+1);
pp1 = deg(ipoly,1) + 1;  deg1 = pp1-1
pp2 = deg(ipoly+1,1) + 1; deg2 = pp2-1
thq1 = thqa(1:pp1,ipoly);
thq2 = thqa(1:pp2,ipoly+1);

% denomfactor = 2^(npolys-ipoly) - 1
% errsqtarget = (trtol^2 - errsum^2)/denomfactor
errsqtarget = (trtol^2 - errsum^2)/(npolys-ipoly); errtarget = sqrt(errsqtarget)

j = 0;
errsq = 1.e+50;
while( (j < nnv2) & (errsq > errsqtarget) )
for i=1:freqcherr2
    j = j + 1;
    b = rand(n,1);
    b = z4noise2(n,b);
    pb1 = oppppNewtNewQCD(n,b,pp1,thq1); 
    pb2 = oppppNewtNewQCD(n,b,pp2,thq2);
    bpAb(j,1) = b'*(pb1-pb2);
    if(deflate == 1)
    adjdefl = 0;
    for id=1:ndefl
        temp = bpAb(j,1);
        lamb = thdefl(id);
        plam1 = fpolypNewtNew(lamb,pp1,thq1);
        plam2 = fpolypNewtNew(lamb,pp2,thq2);
        adjdefl = adjdefl + (plam1 - plam2)*(b'*zr(:,id))*(zl(:,id)'*b)/(zl(:,id)'*zr(:,id));
    end
%                                                     bpAbj = bpAb(j,1)
%                                                     adjdefl
    temp = bpAb(j,1) - adjdefl;
%                                                     ab = abs(bpAbj)
%                                                     at = abs(temp)
    bpAb(j,1) = temp;
    end
end
tra2 = sum(bpAb(1:j,1))/j;
var2 = 0;
for i=1:j
    diff = bpAb(i,1) - tra2;
    var2 = var2 + diff'*diff;
end
var2 = var2/(j-1);  % check this
error2 = sqrt(var2/j);
errsq = error2^2 ;

j
currenterr = error2

end

numnoise = j
numnoia(ipoly+1,1) = j;
error2
% errsq
errsum = sqrt(errsum^2 + error2^2)

end

numnoia

% [S,U,V]=GolubKahQCD(n,pp1,pp2,thq1,thq2);

% %%%%%%%%%%%%%%%%%%%%
% %PHASE 3, deal with trace of p_npolys(A)
% 
% pp = deg(npolys,1) + 1;
% nnv3 = nnva(npolys+1,1)
% for i=1:nnv3
%     b = rand(n,1);
%     b = z4noise2(n,b);
%     pb = oppppNewtNew(n,b,pp,thq);
%     bpAb(i,1) = b'*pb;
% end
% %
% %%%%%%%%%
% %TRACE and Variance
% tra3 = sum(bpAb)/nnv3;
% var3 = 0;
% for i=1:nnv3
%     diff = bpAb(i,1) - tra3;
%     var3 = var3 + diff'*diff;
% end
% var3 = var3/(nnv3-1);  % check this
% error3 = sqrt(var3/nnv3)
% % stdev = sqrt(var2)


% its = iter - pp
% % itertimesp = iter*pp
% mvp
% vops
% dps
% % rnfinal = norm(b(:,1) - op(n,xreal(:,1)))  
% pp
%
toc
%
figure(1)
clf reset
semilogy(gdr(:,1),'linewidth',3);
% xlim([0 250])
ylim([1.e-10 1])
xlabel('Matrix-vector products','Fontsize',13)
ylabel('Residual Norm','Fontsize',13)

% figure(2)
% clf reset
% semilogy(rreala)
% %
% figure(3)
% clf reset
% semilogy(rreala,'linewidth',3)
% hold on
% semilogy(rna,'-.','linewidth',3)


%
%
%
% %-------------SECTION FOR PLOTTING POLYNOMIAL---------
% 

% %  TWO-DIM
% xx = [0:.1:5000];
% for i=1:50001
%     yy(i) = fpoly(xx(i),p,co);
% end
% figure(10)
% clf reset
% plot(xx,yy)

%  TWO-DIM
% xx = [0:.1:322];
% for i=1:3221
% xx = [-2000:.1:2000];
% for i=1:40001
%     yy(i) = 1-fpolyr(xx(i),pp,thq);
% end
% figure(10)
% clf reset
% plot(xx,yy,'g','linewidth',3)
% hold on
% xx = [-2000:-1,1:2000];
% clear yy
% for i=1:n
%     yy(i) = 1-fpolyr(xx(i),pp,thq);
% end
% plot(xx,yy,'r+','linewidth',1)

% %  TWO-DIM for simple 1-D Helmholtz, gamma = 100
% %  NOTE: MOVED this to gppf4 program.
% xx = [-1000:.01:1000];
% for i=1:200001
%     yy(i) = 1-fpolyr(xx(i),pp,thq);
% end
% figure(10)
% clf reset
% plot(xx,yy,'g','linewidth',3)


%  TWO-DIM for simple 1-D Helmholtz, gamma = 100
%  NOTE: MOVED this to gppf4 program.

% xx = [-101:.01:20000];
% for i=1:2010101
%     yy(i) = 1-fpolyr(xx(i),pp,thq);
% end
% figure(10)
% clf reset
% plot(xx,yy,'g','linewidth',3)
% % hold on
% % xx = [-100,.1:.1:1,2:20000]';
% % % xx = [.01,.1,1:998];
% % clear yy
% % for i=1:n 
% %     yy(i) = 1-fpolyr(xx(i),pp,thq);
% % end
% % plot(xx,yy,'r+','linewidth',1)

% ev = eig(full(A));
% for i=1:n
%     evpp(i,1) = fpoly(ev(i,1),p,co);
% end
% hold on
% plot(ev,evpp,'go')
% xlim([-1 3.002])

%  THREE-DIM for matrix with complex eigenvalues in circle
iii = sqrt(-1);
xx = [0:.01:2];
yy = [-1:.01:1];
for i=1:201
    for j=1:201
        xy = xx(i) + iii*yy(j);
        zz(j,i) = fpolyr(xy,pp,thq);
        zza(j,i) = abs(fpolyr(xy,pp,thq));
%         zz(j,i) = abs(fpoly(xy,p,co));
%         zz(j,i) = 1-abs(1-fpoly(xy,p,co));
    end
end
figure(11)
clf reset
zzr = real(zz);
mesh(xx,yy,zzr)
% zlim([0 2.5])         % 
% colormap(gray); 
zlim([0 1.25])
xlabel('Real Axis','Fontsize',12)
ylabel('Imaginary Axis','Fontsize',12)
zlabel('Real Val of Polynomial','Fontsize',12)

figure(12)
clf reset
zzi = imag(zz);
mesh(xx,yy,zzi)
% zlim([0 2.5])         % 
% colormap(gray); 
zlim([-1 .25])
xlabel('Real Axis','Fontsize',12)
ylabel('Imaginary Axis','Fontsize',12)
zlabel('Imag Val of Polynomial','Fontsize',12)

figure(13)
clf reset
mesh(xx,yy,zza)
% zlim([0 2.5])         % 
% colormap(gray); 
zlim([0 1.25])
xlabel('Real Axis','Fontsize',12)
ylabel('Imaginary Axis','Fontsize',12)
zlabel('Absolute Val of Polynomial','Fontsize',12)

% 
% ev = eig(full(A));
% for i=1:n
%     evpp(i,1) = fpoly(ev(i,1),p,co);
% end
% figure(12)
% clf reset
% plot(ev,'go')
% hold on
% plot(evpp,'r*')
% xlabel('Real Axis','Fontsize',13)
% ylabel('Imaginary Axis','Fontsize',13)
% text(1.58,0.4,'Spectrum','Fontsize',12)
% text(1.68,0.3,'of A','Fontsize',12)
% text(1.15,0.,'Spectrum','Fontsize',12)
% text(1.2,-0.1,'of p(A)A','Fontsize',12)
% figure(13)
% clf reset
% plot(ev,'go','linewidth',2)
% hold on
% plot(evpp,'r*','linewidth',3)
% xlim([0 .005])
% xlabel('Real Axis','Fontsize',13)
% ylabel('Imaginary Axis','Fontsize',13)
% text(0.003,-0.05,'Eigenvalues of A:','Fontsize',12)
% text(0.00335,-0.06,'asterisks','Fontsize',12)
% text(0.0003,-0.072,'Eigenvalues of p(A)A:','Fontsize',12)
% text(0.0008,-0.082,'circles','Fontsize',12)
% % 

% %
% % p = 4
% p = 6
% % 
% imag = sqrt(-1); 
% clear co
% % % co(1,1) = 8.23*10^(-2)-(2.24*10^(-4))*imag;
% % % co(2,1) = -2.93*10^(-3)+(1.95*10^(-6))*imag;
% % % co(3,1) = 5.02*10^(-5)+(1.42*10^(-7))*imag;
% % % co(4,1) = -3.37*10^(-7)-(2.53*10^(-9))*imag;
% temp =   [0.121684976415471       1.265914029070125E-003
%  -7.478541098766387E-003  -1.429600307203536E-004
%   2.742738430053471E-004  7.157288747416652E-006
%  -6.061929180859466E-006  -1.885287631710066E-007
%   7.463097191631308E-008  2.546266347524936E-009
%  -3.935788023600409E-010  -1.396667139038465E-011];
% for i=1:6
%     co(i,1) = temp(i,1) + temp(i,2)*imag;
% end
% 
% %  TWO-DIM
% xx = [0:.1:80];
% for i=1:801
%     yy(i) = fpoly(xx(i),p,co);
% end
% figure(10)
% clf reset
% % plot(xx,yy)
% plot(xx,abs(yy))
% %
% %  THREE-DIM
% % xx = [0:.1:80];
% % yy = [-40:.1:40];
% xx = [0:.1:60];
% yy = [-30:.1:30];
%         ylim([-30 30])
% for i=1:601
%     for j=1:601
%         xy = xx(i) + imag*yy(j);
% %         zz(j,i) = 1 - abs(1-fpoly(xy,p,co));
%         zz(j,i) = abs(fpoly(xy,p,co));
%     end
% end
% figure(11)
% clf reset
% mesh(xx,yy,zz)
% % zlim([0 2])
% % evsep = [4.22764881336477       -10.7774842112737        
% %     7.77219798867201       -16.7111738989533        
% %     3.28638001620863        8.53199499999199        
% %     13.9104152088469       -20.6188732667854        
% %     20.1136163069666       -23.2568602485243        
% %     7.44434009358470        15.1682076530793        
% %     27.9383975351494       -23.5081929297608        
% %     12.2662117879156        20.3912730697859        
% %     35.0206287365598       -22.1711830862555        
% %     18.6701794854557        22.5330356170241        
% %     41.9301428584312       -18.5885106795300        
% %     47.4518703334573       -14.0464612618020        
% %     51.7047311863208       -8.69912899189167        
% %     53.3423154093707       -2.37616080011424        
% %     25.9117598295351        23.2825631718495        
% %     53.0158742129377        4.02625383953611        
% %     50.5344521668820        10.0667861530487        
% %     46.0304054142224        15.0103783491194        
% %     33.6573628143367        22.4081873475982        
% %     40.3832049573956        19.0615142497454];
% % for i=1:20
% %     ev(i) = evsep(i,1) + imag*evsep(i,2);
% % end
% %     for j=1:20
% % %         zval(j) = 1 - abs(1-fpoly(ev(i),p,co));
% %         zval(j) = abs(fpoly(ev(i),p,co));
% %     end
% % hold on
% % plot3(evsep(:,1),evsep(:,2),zval,'k*')
% % figure(21)
% % clf reset
% % plot(ev,'*')
