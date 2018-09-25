clear all
close all
clc

dt = 1 / 500;
len = 10000;
time = dt * (1 : len);

beta = 1.9535;

quaternion_true = zeros(len, 4);
quaternion = zeros(len, 4);
quaternion_lasso = zeros(len, 4);
gains = zeros(len, 1);
eigen = zeros(len, 1);

% base1 = randn(4, 1);
% base2 = randn(4, 1);

base1 = [
  -0.833410380239347;
  -1.583349855149863;
   3.003802968285059;
  -1.120003595374574
  ];
base2 = [
   1.367881096733592;
  -0.147989676456610;
   2.006142475708269;
  -0.017958829584065
];
q = [1; 0; 0; 0];
qq = q;



for i = 1 : len
    quat = sin(base1 * i * 0.001 + base2) + 0.02 * randn(4, 1);
    quat = quat ./ norm(quat);
    quaternion_true(i, :) = quat';
    C_true = quat2dcm(quat');
    
    C_true = C_true / 3 + 0.01 * randn(3, 3);
    
    
    z = [C_true(2, 3) - C_true(3, 2); 
         C_true(3, 1) - C_true(1, 3); 
         C_true(1, 2) - C_true(2, 1)];
    K = [trace(C_true), z';
         z, C_true + C_true' - trace(C_true) * eye(3)];
    
    [V, D] = eig(K);
    q = V(:, 4);
    if(q(1) < 0)
        q = - q;
    end
    quaternion(i, :) = q';
    quaternion(i, :) = quaternion(i, :) ./ norm(quaternion(i, :));
     
    [V, D] = eig(K);
    eigen(i) = D(4, 4);
    qm = V(:, 4);
    F = qm * qm';
    qq = F * qq;
    qq = qq ./ norm(qq);
    quaternion_lasso(i, :) = qq';
    
end

figure(1);
subplot(2, 2, 1);
plot(time, quaternion(:, 1), 'LineWidth', 2); hold on
plot(time, quaternion_true(:, 1), 'LineWidth', 1); hold off
xlim([0 max(time)]);
xlabel('Time (s)');
ylabel('q_0');
title('q_0 Bar-Itzhack');
legend('Computed', 'Reference');

subplot(2, 2, 2);
plot(time, quaternion(:, 2), 'LineWidth', 2); hold on
plot(time, quaternion_true(:, 2), 'LineWidth', 1); hold off
xlim([0 max(time)]);
xlabel('Time (s)');
ylabel('q_1');
title('q_1 Bar-Itzhack');
legend('Computed', 'Reference');

subplot(2, 2, 3);
plot(time, quaternion(:, 3), 'LineWidth', 2); hold on
plot(time, quaternion_true(:, 3), 'LineWidth', 1); hold off
legend('Computed', 'Reference');
xlim([0 max(time)]);
xlabel('Time (s)');
ylabel('q_2');
title('q_2 Bar-Itzhack');

subplot(2, 2, 4);
plot(time, quaternion(:, 4), 'LineWidth', 2); hold on
plot(time, quaternion_true(:, 4), 'LineWidth', 1); hold off
legend('Computed', 'Reference');
xlim([0 max(time)]);
xlabel('Time (s)');
ylabel('q_3');
title('q_3 Bar-Itzhack');


figure(4);
subplot(2, 2, 1);
plot(time, quaternion_lasso(:, 1), 'LineWidth', 2); hold on
plot(time, quaternion_true(:, 1), '--', 'LineWidth', 1); hold off
legend('Proposed', 'Reference');
xlim([0 max(time)]);
xlabel('Time (s)');
ylabel('q_0');
title('q_0 Proposed');

subplot(2, 2, 2);
plot(time, quaternion_lasso(:, 2), 'LineWidth', 2); hold on
plot(time, quaternion_true(:, 2), '--', 'LineWidth', 1); hold off
legend('Proposed', 'Reference');
xlim([0 max(time)]);
xlabel('Time (s)');
ylabel('q_1');
title('q_1 Proposed');

subplot(2, 2, 3);
plot(time, quaternion_lasso(:, 3), 'LineWidth', 2); hold on
plot(time, quaternion_true(:, 3), '--', 'LineWidth', 1); hold off
legend('Proposed', 'Reference');
xlim([0 max(time)]);
xlabel('Time (s)');
ylabel('q_2');
title('q_2 Proposed');

subplot(2, 2, 4);
plot(time, quaternion_lasso(:, 4), 'LineWidth', 2); hold on
plot(time, quaternion_true(:, 4), '--', 'LineWidth', 1); hold off
legend('Proposed', 'Reference');
xlim([0 max(time)]);
xlabel('Time (s)');
ylabel('q_3');
title('q_3 Proposed');


figure(6);
subplot(2, 2, 1);
plot(time, quaternion(:, 1) - quaternion_true(:, 1), 'LineWidth', 1);
xlim([0 max(time)]);
xlabel('Time (s)');
ylabel('Error q_0');
title('Error q_0 Bar-Itzhack');

subplot(2, 2, 2);
plot(time, quaternion(:, 2) - quaternion_true(:, 2), 'LineWidth', 1);
xlim([0 max(time)]);
xlabel('Time (s)');
ylabel('Error q_1');
title('Error q_1 Bar-Itzhack');

subplot(2, 2, 3);
plot(time, quaternion(:, 3) - quaternion_true(:, 3), 'LineWidth', 1);
xlim([0 max(time)]);
xlabel('Time (s)');
ylabel('Error q_2');
title('Error q_2 Bar-Itzhack');

subplot(2, 2, 4);
plot(time, quaternion(:, 4) - quaternion_true(:, 4), 'LineWidth', 1);
xlim([0 max(time)]);
xlabel('Time (s)');
ylabel('Error q_3');
title('Error q_3 Bar-Itzhack');



figure(7);
subplot(2, 2, 1);
plot(time, quaternion_lasso(:, 1) - quaternion_true(:, 1), 'LineWidth', 1);
xlim([0 max(time)]);
xlabel('Time (s)');
ylabel('Error q_0');
title('Error q_0 Proposed');

subplot(2, 2, 2);
plot(time, quaternion_lasso(:, 2) - quaternion_true(:, 2), 'LineWidth', 1);
xlim([0 max(time)]);
xlabel('Time (s)');
ylabel('Error q_1');
title('Error q_1 Proposed');

subplot(2, 2, 3);
plot(time, quaternion_lasso(:, 3) - quaternion_true(:, 3), 'LineWidth', 1);
xlim([0 max(time)]);
xlabel('Time (s)');
ylabel('Error q_2');
title('Error q_2 Proposed');clear all
close all
clc

dt = 1 / 500;
len = 10000;
time = dt * (1 : len);

beta = 1.9535;

quaternion_true = zeros(len, 4);
quaternion = zeros(len, 4);
quaternion_lasso = zeros(len, 4);
gains = zeros(len, 1);
eigen = zeros(len, 1);

% base1 = randn(4, 1);
% base2 = randn(4, 1);

base1 = [
  -0.833410380239347;
  -1.583349855149863;
   3.003802968285059;
  -1.120003595374574
  ];
base2 = [
   1.367881096733592;
  -0.147989676456610;
   2.006142475708269;
  -0.017958829584065
];
q = [1; 0; 0; 0];
qq = q;



for i = 1 : len
    quat = sin(base1 * i * 0.001 + base2) + 0.02 * randn(4, 1);
    quat = quat ./ norm(quat);
    quaternion_true(i, :) = quat';
    C_true = quat2dcm(quat');
    
    C_true = C_true / 3 + 0.01 * randn(3, 3);
    
    
    z = [C_true(2, 3) - C_true(3, 2); 
         C_true(3, 1) - C_true(1, 3); 
         C_true(1, 2) - C_true(2, 1)];
    K = [trace(C_true), z';
         z, C_true + C_true' - trace(C_true) * eye(3)];
    
    [V, D] = eig(K);
    q = V(:, 4);
    if(q(1) < 0)
        q = - q;
    end
    quaternion(i, :) = q';
    quaternion(i, :) = quaternion(i, :) ./ norm(quaternion(i, :));
     
    [V, D] = eig(K);
    eigen(i) = D(4, 4);
    qm = V(:, 4);
    F = qm * qm';
    qq = F * qq;
    qq = qq ./ norm(qq);
    quaternion_lasso(i, :) = qq';
    
end

figure(1);
subplot(2, 2, 1);
plot(time, quaternion(:, 1), 'LineWidth', 2); hold on
plot(time, quaternion_true(:, 1), 'LineWidth', 1); hold off
xlim([0 max(time)]);
xlabel('Time (s)');
ylabel('q_0');
title('q_0 Bar-Itzhack');
legend('Computed', 'Reference');

subplot(2, 2, 2);
plot(time, quaternion(:, 2), 'LineWidth', 2); hold on
plot(time, quaternion_true(:, 2), 'LineWidth', 1); hold off
xlim([0 max(time)]);
xlabel('Time (s)');
ylabel('q_1');
title('q_1 Bar-Itzhack');
legend('Computed', 'Reference');

subplot(2, 2, 3);
plot(time, quaternion(:, 3), 'LineWidth', 2); hold on
plot(time, quaternion_true(:, 3), 'LineWidth', 1); hold off
legend('Computed', 'Reference');
xlim([0 max(time)]);
xlabel('Time (s)');
ylabel('q_2');
title('q_2 Bar-Itzhack');

subplot(2, 2, 4);
plot(time, quaternion(:, 4), 'LineWidth', 2); hold on
plot(time, quaternion_true(:, 4), 'LineWidth', 1); hold off
legend('Computed', 'Reference');
xlim([0 max(time)]);
xlabel('Time (s)');
ylabel('q_3');
title('q_3 Bar-Itzhack');


figure(4);
subplot(2, 2, 1);
plot(time, quaternion_lasso(:, 1), 'LineWidth', 2); hold on
plot(time, quaternion_true(:, 1), '--', 'LineWidth', 1); hold off
legend('Proposed', 'Reference');
xlim([0 max(time)]);
xlabel('Time (s)');
ylabel('q_0');
title('q_0 Proposed');

subplot(2, 2, 2);
plot(time, quaternion_lasso(:, 2), 'LineWidth', 2); hold on
plot(time, quaternion_true(:, 2), '--', 'LineWidth', 1); hold off
legend('Proposed', 'Reference');
xlim([0 max(time)]);
xlabel('Time (s)');
ylabel('q_1');
title('q_1 Proposed');

subplot(2, 2, 3);
plot(time, quaternion_lasso(:, 3), 'LineWidth', 2); hold on
plot(time, quaternion_true(:, 3), '--', 'LineWidth', 1); hold off
legend('Proposed', 'Reference');
xlim([0 max(time)]);
xlabel('Time (s)');
ylabel('q_2');
title('q_2 Proposed');

subplot(2, 2, 4);
plot(time, quaternion_lasso(:, 4), 'LineWidth', 2); hold on
plot(time, quaternion_true(:, 4), '--', 'LineWidth', 1); hold off
legend('Proposed', 'Reference');
xlim([0 max(time)]);
xlabel('Time (s)');
ylabel('q_3');
title('q_3 Proposed');


figure(6);
subplot(2, 2, 1);
plot(time, quaternion(:, 1) - quaternion_true(:, 1), 'LineWidth', 1);
xlim([0 max(time)]);
xlabel('Time (s)');
ylabel('Error q_0');
title('Error q_0 Bar-Itzhack');

subplot(2, 2, 2);
plot(time, quaternion(:, 2) - quaternion_true(:, 2), 'LineWidth', 1);
xlim([0 max(time)]);
xlabel('Time (s)');
ylabel('Error q_1');
title('Error q_1 Bar-Itzhack');

subplot(2, 2, 3);
plot(time, quaternion(:, 3) - quaternion_true(:, 3), 'LineWidth', 1);
xlim([0 max(time)]);
xlabel('Time (s)');
ylabel('Error q_2');
title('Error q_2 Bar-Itzhack');

subplot(2, 2, 4);
plot(time, quaternion(:, 4) - quaternion_true(:, 4), 'LineWidth', 1);
xlim([0 max(time)]);
xlabel('Time (s)');
ylabel('Error q_3');
title('Error q_3 Bar-Itzhack');



figure(7);
subplot(2, 2, 1);
plot(time, quaternion_lasso(:, 1) - quaternion_true(:, 1), 'LineWidth', 1);
xlim([0 max(time)]);
xlabel('Time (s)');
ylabel('Error q_0');
title('Error q_0 Proposed');

subplot(2, 2, 2);
plot(time, quaternion_lasso(:, 2) - quaternion_true(:, 2), 'LineWidth', 1);
xlim([0 max(time)]);
xlabel('Time (s)');
ylabel('Error q_1');
title('Error q_1 Proposed');

subplot(2, 2, 3);
plot(time, quaternion_lasso(:, 3) - quaternion_true(:, 3), 'LineWidth', 1);
xlim([0 max(time)]);
xlabel('Time (s)');
ylabel('Error q_2');
title('Error q_2 Proposed');

subplot(2, 2, 4);
plot(time, quaternion_lasso(:, 4) - quaternion_true(:, 4), 'LineWidth', 1);
xlim([0 max(time)]);
xlabel('Time (s)');
ylabel('Error q_3');
title('Error q_3 Proposed');

figure(8);
plot(time, eigen, 'LineWidth', 1); 
xlabel('Time (s)');
ylabel('Eigenvalue ${\lambda}_{\max}$', 'Interpreter', 'latex');



subplot(2, 2, 4);
plot(time, quaternion_lasso(:, 4) - quaternion_true(:, 4), 'LineWidth', 1);
xlim([0 max(time)]);
xlabel('Time (s)');
ylabel('Error q_3');
title('Error q_3 Proposed');

figure(8);
plot(time, eigen, 'LineWidth', 1); 
xlabel('Time (s)');
ylabel('Eigenvalue ${\lambda}_{\max}$', 'Interpreter', 'latex');

