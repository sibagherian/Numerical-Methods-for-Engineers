% Compute the Feigenbaum delta
% Store approximate values in the row vector delta for assessment, where length(delta)= num_doublings and 
% delta(2:num_doublings) are computed from the algorithm described in Lectures 21-23.
num_doublings=11; delta=zeros(1,num_doublings); delta(1)=5;
% Write your code here
m = zeros(1,num_doublings);
m(1)= 2; 
m(2)= 1 + sqrt(5); 
x0  = 1/2; 
xp0 = 0; 
for n = 2: num_doublings
  period = 2^n ;
  m(n+1) = m(n) + ( m(n) - m(n-1))/delta(n-1);
  for j = 1:10
    for k = 1:period
      x  = m(n+1) * x0 * (1 - x0);
      xp = x0 * (1 - x0) + m(n+1) * xp0 * (1 - 2*x0) ;
      x0 = x ;
      xp0 = xp ;
    end
    u(j) = m(n+1) - (x -1/2 )/xp;
    m(n+1) = u(j) ;
  end
%  m(n+1) = u(j);
  delta(n) = (m(n) - m(n-1))/(m(n+1) - m(n));
end
% Output your results
fprintf('n        delta(n)\n');
for n=1:num_doublings
    fprintf('%2g %18.15f\n',n,delta(n));
end
