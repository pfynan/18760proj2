infile = 'benchmarks/toy/toy1';
outfile = 'outfile.txt';
fd = fopen(infile);

chipx = fscanf(fd,'%f',1);
chipy = fscanf(fd,'%f',1);
fundunit = fscanf(fd,'%f',1);
ngates = fscanf(fd,'%f',1);
nnets = fscanf(fd,'%f',1);

b = zeros(ngates,nnets);
nconn = zeros(ngates,1);
for i = 1:ngates
    gnum = fscanf(fd,'%f',1);
    assert(i==gnum);
    nconn(i) = fscanf(fd,'%f',1);
    c_net = fscanf(fd,'%f', nconn(i)); 
    b(gnum,c_net) = 1; % e.g. Column gives net 1(y_coord) connects to gate 2(x_coord) and gate3(x_coord)
end

npins = fscanf(fd,'%f',1);
pinc = zeros(npins,nnets);
pinp = zeros(npins,2);

for i = 1:npins
    pnum = fscanf(fd,'%f',1);
    assert(i==pnum);
    net = fscanf(fd,'%f',1);
    pinc(pnum,net) = 1;  % e.g. Column gives net 2(y_coord) connects to gate 4(x_coord); FIXED
    pos = fscanf(fd,'%f',2);
    pinp(pnum,:) = pos;
end

fclose(fd);
% Update via Gradient Descent

cellx = randi(chipx, ngates,1); %Random X coord of gates
celly = randi(chipy, ngates,1); % Random Y coord

state.b = b;
state.ngates = ngates;
state.chipx = chipx;
state.chipy = chipy;
state.nnets = nnets;
state.pinc = pinc;
state.pinp = pinp;
state.nconn = nconn;
state.fundunit = fundunit;


% Martin King, ICTP, 29 March 2005
%1. Conjugate Gradient Method with Flecther-Reeves (or Polak-Ribiere)
%   to find a vector x that gives a MINIMUM of a function (a scalar).
%2. Ideas taken from J.R. Shewchuk and Numerical Recipes.
%3. You must modify your own function to minimise in a Matlab 
%   function file called func.m (scalar output) and the first derivative of that
%   function in a matlab function file called dfunc.m, which has a vector output in 
%   (del/del(x1) del/del(x2) ... etc)'. 
%4. If you want to MAXIMISE a function, multiply -1 to the output of func.m
%   and dfunc.m (be careful here, dfunc.m may use func.m; to be safe, give dfunc.m
%   the original output of func.m and then multiply -1 to dfunc.m at the end). 
%   If you are brave, reverse the search direction r and d for maximisation.
%5. If the method is not converging or is giving you a solution that doesn't make sense, 
%   change the initial guess.
%6. As an example, a simple function is given in func.m and its gradient vector in dfunc.m. 
%   Change the initial guess to x = [1 ; 1] for example, the solution it gives is incorrect. 
%   The reason is obvious if you plot the function (it is the saddle points). 
%7. I have used these scripts to optimise a fairly complicated function. They 
%   seem to work well. If you notice any bug or have any comment, please email me
%   king at ictp at it

%put initial guess here (x is an n-dimensional vector).
%x=[-2 ; -3];
x = ones(ngates*2, 1);
x = [rand(ngates,1) * chipx;
     rand(ngates,1) * chipy];

[F] = func(x,state);
[F_prime] = dfunc(x,state);

%From here in this m-file, I usually follow the notations in the 'pseudocode' of B4 
%in Shewchuk's note.

r=-1.*F_prime;
d=r;
%g is for use in Polak-Ribiere
g=r;
delta_new=r'*r; 
delta_0=delta_new;
fp = F;
%ftol is a convergence tolerance
ftol=1.e-7;
ITERMAX = 10000;
%Don't worry too much about this
EPS=1.e-10;

for iter = 1 : ITERMAX 

%Doing the line search here. First bracket a minimum, then use Golden section to find it.
%Not using Newton-Raphson as in Shewchuk. So you don't need the second derivative.
  [ax,bx,cx,fa,fb,fc] = func_mnbrak(0,1,x,d,state);
  [xt,golden] = func_golden(ax,bx,cx,x,d,state);
%To recover vector x, which is along d at xt away from initial x.
   x = x + xt.*d;
%The function value at x is golden as returned by func_golden.
   F = golden;
   [F_prime] = dfunc(x,state);

   r = -1.*F_prime;
   delta_old = delta_new;
%This is Fletcher-Reeves
   delta_new = r'*r;
%This is Polak-Ribiere
%   delta_new = (F_prime+g)'*F_prime;
   beta = delta_new/delta_old;
   d = r + beta * d ; 
   g = r;
   if r'*d <= 0
      d=r;
   end
%this convergence criterion is taken from NR.
   if 2.*abs(F-fp) < ftol.*(abs(F)+abs(fp)+EPS)
      iter
      x
      F
      break;
   end
   fp = F;

   if mod(iter,50) == 0
       fprintf(1,'Iteration %d\n',iter);
   end
end

x_soln = x(1:ngates);
y_soln = x((ngates+1):(2*ngates));


fd = fopen(outfile,'w');
fprintf(fd,'%d %f %f\n',[ 1:size(x_soln,1) ; x_soln'; y_soln']);
fclose(fd);


