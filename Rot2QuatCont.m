% Optimal Continuous Unit Quaternions from Rotation Matrices
% 
% Reference:
% Wu, J. (2018) Optimal Continuous Unit Quaternions from Rotation Matrices.
%        Journal of Guidance, Control and Dynamics
%
% (c) Jin Wu
% e-mail: jin_wu_uestc@hotmail.com


function q = Rot2QuatCont(C, last_q)
    z = [C(2, 3) - C(3, 2); 
         C(3, 1) - C(1, 3); 
         C(1, 2) - C(2, 1)];
    K = [trace(C), z';
         z, C + C' - trace(C) * eye(3)];
    [V, D] = eig(K);
    q = V(:, 4);
    
    if(last_q ~= 0)
        q = sign(q' * last_q) * q;
    end
end