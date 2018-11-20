syms Fx Fxw Fy Fyw Tapp Tresb Tresc B C m a real
syms th dth ddth ph dph ddph x dx ddx y dy ddy real
syms k1 real

assume(B~=0); assume(C~=0); assume(m~=0); assume(k1==0);

constraint = dx*sin(th) - dy*cos(th) + a*dth;

eq1 = (Fx + Fxw - m*ddx)*cos(th) + (Fy + Fyw - m*ddy)*sin(th) == 0;
eq2 = Fxw*a*sin(th) - Fyw*a*cos(th) + Tresc + Tresb - (C+B)*ddth - B*ddph == 0;
eq3 = Tapp + Tresb - B*ddth - B*ddph == 0;
eq4 = dx*dth*cos(th)+ddx*sin(th)+dy*dth*sin(th) - ddy*cos(th) + a*ddth == 0;
eq4 = dx*dth*cos(th)+ddx*sin(th)+dy*dth*sin(th) - ddy*cos(th) + a*ddth == -k1*constraint;

eqns = [eq1;eq2;eq3;eq4];
vars = [ddx;ddy;ddth;ddph];

solutions = solve(eqns,vars);

fields = fieldnames(solutions)

for i=1:numel(fields)
  solutions.(fields{i}) = simplify(solutions.(fields{i}),100);
  solutions.(fields{i})
end

%% LQR

f1 = [dx;dy;dth;dph;solutions.ddx;solutions.ddy;solutions.ddth;solutions.ddph]

dfdx = jacobian(f1,[x,y,th,ph,dx,dy,dth,dph])

dfdu = jacobian(f1,[Tapp])


%[A,b] = equationsToMatrix(eqns,vars)