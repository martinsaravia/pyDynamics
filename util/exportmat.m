

A1 = [9.9, 9900];
A2 = [8.8,  7.7 ; ...
      8800, 7700];
formatSpec = 'X is %4.2f meters or %8.3f mm\n';
fprintf(formatSpec,A1,A2)

node = 1;
dofs = [100, 101, 102]
steps = len(U(:,1))
X0 = zeros(steps, 3)
V0 = zeros(steps, 3)
A0 = zeros(steps, 3)
VO = zeros(steps, 3)

for i=1:steps;
    X0(i,:) = U(i,dofs);
end


