%
%  For QCD, need to adjust to deal with complex harm Ritz val's not in
%  pairs !

%  Poly Precondioned GMRES - 2017 version
%   
%
clear all;
format long e

global A

% tic

%load saylor4;
%load pores2;
%load sher5;
%load sher2;
%load sher1;
%load fid18;
%load sher4;
% load burak; A = burak;
% load qcd4; A = qcd;    
% load matmarket888000; A = qcd;
%                       n = 49152;
%                       spidentity = spdiags(ones(n,1),0,n,n);
%                       A = A + 6.4*spidentity;

% load s1rmq4m1
% load bwm2000
% load olm1000
% load af23560
% % %%%%%%%%%%%%%%%%%%%%%%%%%  Don't know if right:  !!!!!!!!!
% A = af23560;  A(23560,23560) = -71.155520800000005;
% load e20r0100; A = testMatrix.A;
% load quarkmatrix_12_16; A = M; clear M;
% load h01; n = size(B,1); spI = speye(n,n); A = spI - 0.1570*B;
% load h01_4444; n = size(B,1); spI = speye(n,n); %A = spI - 0.1570*B;
load h01_8888; n = size(B,1); spI = speye(n,n); %A = spI - 0.1570*B;
                        A = spI - .175*B;
% load XyceTest3M6
n = size(A,1);
% global L U
%global D

%[L,U] = luinc(A,'0');
%D = diag(A);
%T = diag( diag(A) ) + diag( diag(A,1), 1) + diag( diag(A,-1), -1);
%[L,U] = luinc(T,'0');
%
%n = 3564
%n = 1224
%n = 3312
%n = 3079
%n=5773
%n = 1080
%n = 1104
%n = 2400 (for Burak)
% n = 10000
% n = 1536
% n = 49152  % (for qcd matrix matmarket888000)
%                                           A = A - .3*eye(n,n);
%n = 15552
% %%randn('seed',0)
% %A = diag(randn(n,1));
% xxx = [1:n]';
% xxx = [.01,.1,1:n-2]';
% xxx = [.1:.1:.9,1:n-9]';
% xxx = [.01:.01:.09,.1:.1:.9,1:n-18]';
% % xxx = [1:n-1,2*n]';
% xxx = [.1:.1:.9,1:n-10,5200]';
% xxx = [.1:.1:.9,1:n-12,7500,8000,8300]';
%     xxx = [-1000:-11,-1.e-8,11:1000]'; n = size(xxx,1);
%     xxx = [-1000:-1,1:1000]'; n = size(xxx,1);
%     xxx = [-100,.1:.1:1,2:20000]'; n = size(xxx,1);
%     xxx = [-100,1:20000]'; n = size(xxx,1);
% xxx = [-n/2:-1,1:n/2]'; n = size(xxx,1);
% % % % xxx = [10.1:.1:100,101:n-800]';
% p = 2.0 ; %1.75;
% xxx = [1:n]';
% xxx = (xxx.^p)/n^(p-1);
% A = spdiags(xxx,0,n,n); % + .2*spdiags(ones(n,1),1,n,n);
% A = spdiags(xxx,0,n,n) + spdiags(ones(n,1),1,n,n);
%                                                         A = A/10;
%                                                         b = b/10;
% e = ones(n,1);
% A = spdiags([e -2*e e], -1:1, n, n)

% rng('default')
% A = randn(n,n);
% A = A + 40.*eye(n,n);

% %

tic
j
pp = 15 % Degree of poly  (degree of alpha*p(alpha))
m = 50   % Dimension of Krylov subspace
k = 0    % Number of approximate eigenvectors saved at the restart
%numev = 1
rtol = 1.e-10;
cyclim = 10000
istab = 1;

kperm = k;
%-------------- RHS vectors.--------------------------%
cpu0 = cputime
% randn('seed',0)
rng('default')
% rng(3)
% b = randn(n,1);
b = randn(n,1); b = b/norm(b);
                b = ones(n,1) + sqrt(-1)*ones(n,1); b = b/norm(b);
                                     
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
j = 1;
rninit = norm(b(:,1))
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
    wv = op(n,v(:,j));
    f = minv(n,wv);
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
% x(:,1) = x(:,1) + v(:,1:pp)*d;
%        vops = vops + pp;
% %
  hh = h(1:pp,1:pp);                             %   hhhhh = hh;
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
  th = ddd(ind);
                               
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
    iatendofloop = i
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

%%%%%%%%%%%%% Look over next dozen lines - not sure what
% w = randn(n,1); w = w/norm(w);
% %                                         w = b;
% qawq = w - opppr(n,w,pp,thq);
% temp = oppppNewt(n,w,pp,cp,thq);
% qawp = w - op(n,temp);
% qawd = norm(qawq - qawp)
%                                         dxtt = norm(xt-temp)
%                                         opxt = op(n,xt);
%                                         optemp = op(n,temp);
%                                         dop = norm(opxt-optemp)
%                                         bmoxt = b - opxt;
%                                         bmote = b - optemp;
%                                         dbm = norm(bmoxt-bmote)
clear h
clear c
clear v

% end of generating the poly
%                                     
                                        
%                                                                
%
%-------------------------------------------------------%
%

%  Now may use results from the run for gen the poly instead of starting
%  over.   old:  
x = zeros(n,1);
j = 1;
%                             r = b;
%                             [x,flag,relres,iterbi,resvec] = bicgstab(A,r,1.e-10,2000000);
%                             figure(66)
%                             clf reset
%                             semilogy(resvec)
%                             iterbi
%                             relres
%                             flag
                            
%                             maxit = 50000
%                             S = 4
%                             r = b;
%                             [X,FLAG,RELRES,iterIDR,resvec,replacements] = IDRS(A,r,S,rtol,maxit);
%                             figure(66)
%                             clf reset
%                             semilogy(resvec)
%                             iterIDR
%                             RELRES
%                             FLAG
%                             realrIDR = b-op(n,X);
%                             realrIDRn = norm(realrIDR)
                            
                            
% r = opppp(n,b(:,1),pp,co); 
%     mvp = mvp + pp - 1;
%     vops = vops + pp;
rorig = b(:,1) - op(n,x);
rn = norm(rorig);
    mvp = mvp + 1;
    vops = vops + 1;
                                                   
% rninitofppsystem = norm(r);
%     vops = vops + 1;

% note minv stuff is not done here, so don't use an minv
vn = rn;
v(:,1) = rorig/vn;
    vops = vops + 1;
c(1,1) = vn;
c(2:m+1,1) = zeros(m,1);
                                                    
while ( (rn/rninit > rtol) & (cycle <= cyclim) ) 
while ( (j <= m) & (rn/rninit > rtol) )
%     wv = op(n,v(:,j));
    wv = oppprQCD(n,v(:,j),pp,thq);
                                     
    f = minv(n,wv);
    iter = iter + 1;
        mvp = mvp + pp;
        vops = vops + pp;
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
	d = h(1:j+1,1:j) \ c(1:j+1) ;  
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
cycle
%                                                 xbef = x(:,1)
% rn = norm(srv)
rna(cycle,1) = rn;
x(:,1) = x(:,1) + v(:,1:m)*d;
    vops = vops + m;
% xreal(:,1) = oppppNewt(n,x(:,1),pp,cp,thq);  
% xreal(:,1) = oppppNewtNew(n,x(:,1),pp,thq); %  NEED TO check counts for this! (mvp,vops)
% rreal = rorig(:,1) - op(n,xreal);
% rrealn = norm(rreal)
% rreala(cycle,1) = rrealn;
% rtruewithpp = b(:,1) - oppprQCD(n,x,pp,thq);
% rtruewithppnorm = norm(rtruewithpp)
%        vops = vops + m;
%                                                 xaft = x(:,1)
% wv = rorig - opppr(n,x(:,1),pp,thq); %-may remove this and next line%
% rnale(cycle) = norm(wv);
%                                                 oppp(n,x(:,1),p,th)
%  r = minv(n,wv); etc.
r = v(:,1:m+1)*srv;
       vops = vops + m+1;
%                                                  srv
rnwithpp = norm(r)

%
%   hh = h(1:m,1:m);
%   em = zeros(m,1);
%   em(m) = 1;
%   ff = hh' \ em;
%
% %FOM section
% %Comment next line for FOM instead of GMRES:
%   hh(:,m) = hh(:,m) + h(m+1,m)^2 * ff;
%   [g,dd] = eig(hh,'nobalance'); 
%   ddd = diag(dd);
%   dabs = abs(ddd);
%   [thabs,ind] = sort(dabs);
%   th = ddd(ind);
%      thk = th(1:k)
%   gg = g(:,ind(1:k));
%   for i=1:k
%     rho(i) = gg(:,i)'*h(1:m,:)*gg(:,i);
%     tv = h(1:m,:)*gg(:,i)-rho(i)*gg(:,i);
%     tvn = norm(tv);
%     rna(cycle,i) = sqrt( tvn*tvn+ abs(h(m+1,m))^2*abs(gg(m,i))^2 );
%     tha(cycle,i) = th(i);
%     rhoa(cycle,i) = rho(i);
%   end  %for
%   %
%   greal = gg;
%   % Comment this section for complex matrices:
%   
% %   i = 1;
% %   while( i <= k )
% %     if( imag(th(i)) ~= 0 )
% %       if( i ~= k )
% %         greal(:,i+1) = imag(gg(:,i));
% %         greal(:,i) = real(gg(:,i));
% % disp('split complex vector')
% %         i = i + 1;
% %       else
% %         k = k - 1
% %      end  %if
% %     end  %if
% %     i = i + 1;
% %   end  %while
%   %
%   greal(m+1,1:k) = zeros(1,k);
%   greal(:,k+1) = srv;
%   [gon,rr] = qr(greal(:,1:k+1),0);
%   hnew = gon'*h*gon(1:m,1:k);
%   h(kperm+1,:) = zeros(1,m);
%   h(1:k+1,1:k) = hnew;
%   c(1:k+1,1) = gon(:,1:k+1)'*srv(:,1);
%   c(k+2:m+1,1) = zeros(m-k,1);    
%   work = v*gon;
%         vops = vops + m*(k+1) + m+1;
% %                                                      gon
%   v(:,1:k+1) = work;
%   %
%   %section for just reorthog. one vector, v_{k+1}
%   for i = 1:k
%     dot = v(:,i)'*v(:,k+1) ;
%     v(:,k+1) = v(:,k+1) - dot * v(:,i);
%   end  %for
%   v(:,k+1) = v(:,k+1)/norm(v(:,k+1));
%          vops = vops + 2*k+1;
%   % Don't recompute the k+1 line of h.  Is that OK?
%   %
%   % end section
%   %
% %   rv = b(:,1) - op(n,x(:,1));
% %         mvp = mvp + 1;
% %   rntruewithoutpp = norm(rv)
% %         vops = vops + 1;
% %   rn = rntruewithoutpp;
% %   rnaletrwo(cycle,1) = rn;
%   j = k + 1;
%   kold = k;
%   k = kperm;
  j = 1;
  cycle = cycle + 1;
  vn = norm(r);
    vops = vops + 1;
  v(:,1) = r/vn;
  c(1,1) = vn;
end  %while
cpu1 = cputime - cpu0
rntrwoppdivbyrninit = rn/rninit


xreal(:,1) = oppppNewtNew(n,x(:,1),pp,thq); %  NEED TO check counts for this! (mvp,vops)
rreal = rorig(:,1) - op(n,xreal);
rrealn = norm(rreal)
if( pp > 1)
    mvp = mvp + pp;  % pp-1 for oppppNewtNew 
    vops = vops + 2*pp; % 2pp-1 for oppppNewtNew 
else
    mvp = mvp + 1;
end
%
% k = kold;
% vk = v(:,1:k+1);
% hk = h(1:k+1,1:k);

its = iter - pp
% itertimesp = iter*pp
mvp
vops
dps
% rnfinal = norm(b(:,1) - op(n,xreal(:,1)))  
pp
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
xx = [-.5:.01:2.5];
yy = [-1:.01:1];
for i=1:301
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
% % 
% % evsep = [51.6691236689541       -7.04196886099235        
% %     52.8299963489828        3.82126366203520        
% %     47.4864887830559        14.4564777773204        
% %     44.0845179819575       -16.9179058260804        
% %     36.1687404771749        21.9817739600345        
% %     32.9434813080091       -22.7258809628284        
% %     21.7631830609789        23.3319955878390        
% %     18.9423594394479       -22.7364165584338        
% %     12.1405348046398        19.7073418521863        
% %     9.31417740790867       -17.0704437206335        
% %     4.23932551397388        10.9079978972671        
% %     3.55667471432984       -8.25670719315009        
% %     2.70704136537297       -6.84889776659215        
% %     2.07363033517860       -7.17553990558425        
% %     1.90367788778560       -4.22459583925108        
% %     1.75583359164472       -5.19032285832277        
% %     2.14421105171973        7.06778815081143        
% %     1.90438243043914        4.22443987637251        
% %     2.59566654817047        5.63398110638954        
% %     1.75630177076622        5.18864345710068];
% % 
% % % sent by Quan(Andy) Liu on 2-12-13  First is after first cycle of
% % % g-dr(20,8), then at last cycle.  Next are similar results with pp.
% % evsep = [1.91835364718298       -10.0208979172880        
% %     5.48761562018846       -15.9800167171717        
% %     10.5672732073321       -20.5374788180325        
% %     2.04218258063942        10.6921123396977        
% %     17.3545566464748       -23.2775539877165        
% %     6.10777268593536        16.4958993102093        
% %     24.5981078706375       -23.7387789831597        
% %     11.5054130369249        20.7730541466863        
% %     32.0976472579974       -22.5592556115600        
% %     18.3368325887324        23.3095803460793        
% %     39.0178003384613       -19.4501516877685        
% %     44.8626013936720       -15.1858393365609        
% %     49.1118626977745       -9.69033408887553        
% %     25.5680403487075        23.8834075262292        
% %     51.3290110952449       -3.46044560428488        
% %     51.5113979753517        2.77575703218370        
% %     49.4121188281727        9.11489805612238        
% %     45.4960183468947        14.6224011070231        
% %     33.1960961835481        22.2311196393907        
% %     39.9669149565358        19.1977894854749];    
% % 
% % evsep = [47.8379792820664        5.23377973402470        
% %     48.3907658235861       -4.17622948020382        
% %     42.7967110699765       -13.4590245410070        
% %     41.4914861444724        13.9282553886422        
% %     31.6009349863660       -19.9351304571706        
% %     30.8566832849145        19.4408510853652        
% %     18.6856746793456        20.9489514689298        
% %     19.1083285823252       -19.7963913197891        
% %     10.8917699746469       -17.3645117305850        
% %     8.52352745654333        16.1975528411675        
% %     3.34858435991237       -9.04401410216039        
% %     2.97040602805294        7.55803690805435      
% %     -0.106900684269875       -4.20273341579507      
% %     0.111346328895684       -3.87128492884973      
% %     -0.200136463295853       -2.70957191072258       
% %     0.137198642639356       -3.33425737036532       
% %     0.386337243308860        4.79286552472616      
% %     -0.200176059775461        2.70957590620462      
% %     -0.106654211621985        4.20310692852228       
% %     0.132416251347769        3.33924936946822];
% % 
% % % evsep = [0.250257060635886       -1.24025292396856       
% % %     0.700215697847320       -1.97479401593052        
% % %     1.33573703196873       -2.53366012105129       
% % %     0.240032630572939        1.32837664251180        
% % %     2.18076480857523       -2.86505584851465       
% % %     0.737010837513613        2.05308479679660        
% % %     3.07955691334501       -2.91330371287933        
% % %     1.40105266946183        2.59013091532163        
% % %     4.00806608623279       -2.75777762646162        
% % %     2.24503710374330        2.91310440747475        
% % %     4.86234664290095       -2.36369316292698        
% % %     5.58185363013305       -1.82768732935969        
% % %     6.10198834114566       -1.14098022657645        
% % %     3.14102034027396        2.99319093603766        
% % %     6.36922790559315      -0.365716919740048        
% % %     6.38414301883037       0.407816526160548        
% % %     6.11599743981855        1.19129689646652       
% % %     5.62358749070382        1.86940792453279       
% % %     4.08896309999927        2.79772247423883        
% % %     4.93231136378679        2.42994189842298] ;       
% % % 
% % % evsep = [6.03984360550008      -0.598486937928023        
% % %     6.05617720790522       0.642135407531598        
% % %     5.27822748568857       -1.75016780413585       
% % %     5.28511427210362        1.81858853147034        
% % %     3.83170006839252       -2.59475448839326       
% % %     3.87875427021469        2.52592658518503       
% % %     2.40014008372036        2.64477244935197       
% % %     2.20135985375332       -2.54871600147248       
% % %     1.15427283946536       -2.05367612834619       
% % %     1.10858044298878        2.12700650844587       
% % %     0.281711015411004       -1.01294676027552      
% % %     0.363583539628893        1.03873939844053      
% % %     -8.069324289204094E-003 -0.521281642678995     
% % %     1.873237234390856E-002 -0.479720473244520      
% % %     -2.147135001896407E-002 -0.336242350732565     
% % %     2.106694777255480E-002 -0.413265850606190     
% % %     -2.816864496005418E-002  0.335748457770950    
% % %     1.302885053906668E-002  0.414351320901198     
% % %     -1.843943170915891E-002  0.521011702463749    
% % %     1.841543757969072E-002  0.495988009237250] ;       
% %     
 % sent by Quan(Andy) Liu on 2-14-13  First are h.Ritz vals at end of run
 % for non-PP, then for with PP.
%  evsep = [
% 50.0819597368169       -6.58605364874819        
% 50.4631893077417        5.56959580498800       
% 43.7354350632714        16.7056714290946       
% 42.2924339159846       -17.9084622076601       
% 32.7190928110312       -23.6836832844508        
% 31.2745018385909        24.0925627247128        
% 17.1601141923820        24.3560146630134        
% 17.1868305878498       -24.4150948236098       
% 5.42518525620335        18.7375462810975       
% 6.47163368806274       -19.1236241860240     
% -0.597064192383743        7.86635413634008  
% 0.266245546152205       -8.68137354243959     
% -1.29015453953825        4.19673051378464      
% -1.03316004930448        3.33438091915404      
% -1.36986450259367        2.70957561055731      
% 0.114168478808071      -4.100562266980041E-004 
% -1.36982494237620       -2.70957178713119      
% -1.03312056434207       -3.33437508002253      
% -1.27659916830863       -4.20273837582945      
% -1.05556622094986       -3.86951726218983];
% 
%  evsep2 = [
% 6.07029180735340      -0.771996457465142      
% 6.13725157448738       0.704279299486516     
% 5.30098967598444        2.05805676185901     
% 5.07055277304148       -2.16415766341890     
% 3.99092820396511       -2.85243271811937      
% 3.76634141327604        2.93454437518109      
% 2.08767402244614        3.00501397212764     
% 2.08382615384855       -2.92188434620229      
% 0.676003579387512        2.26722789016635    
% 0.785061661877054       -2.29664450653947     
% -6.945056064155296E-002  0.976686616836749   
% 4.548657606256257E-002  -1.06266193642860    
% -0.160739921500043       0.503540330408797    
% -0.129940866246647       0.404436724174846    
% -0.170122019468733       0.327980513574849    
% 1.389310767461823E-002  9.462987252407037E-005 
% -0.163257030868272      -0.331448259625016      
% -0.121494219340682      -0.407051194846051     
% -0.150022652988462      -0.513026196466926      
% -0.123548010113018      -0.472198420938924 ] 
% 
% for i=1:20
%     ev(i) = evsep(i,1) + imag*evsep(i,2);
% end
% for i=1:20
%     ev2(i) = evsep2(i,1) + imag*evsep2(i,2);
% end
% figure(21)
% clf reset
% plot(ev,'*')
% 
% for j=1:20
%     polyateval(j,1) = fpoly(ev(j),p,co);
% end
% figure(40)
% clf reset
% plot(polyateval,'o')
% hold on
% plot(ev2,'rx')
% 
