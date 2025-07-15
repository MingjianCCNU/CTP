% 示例：计算四维H_{2,1,1,0}(x)
n = [2, 1, 1, 0];
x = [1, 2, 3, 4];
M = [1, 0.5, 0.2, 0;
     0.5, 2, 0.3, 0;
     0.2, 0.3, 3, 0;
     0, 0, 0, 4];

H = multiHermite_ascending(n, x, M);
fprintf('H_{%d,%d,%d,%d}(%.1f,%.1f,%.1f,%.1f) = %.4f\n', n, x, H);

% 验证正确性（手动计算）
expected = x(1)^2*x(2)*x(3) + 2*M(1,1)*x(2)*x(3) + 2*M(1,2)*x(1)*x(3) + 2*M(1,3)*x(1)*x(2);
assert(abs(H - expected) < 1e-10, '计算结果错误');

function result = multiHermite_ascending(n, x, M)
% MULTIHERMITE_ASCENDING 按总阶数递增递归计算多维厄米特多项式
%   n: 多重指标向量 [n1, n2, ..., nd]
%   x: 变量向量 [x1, x2, ..., xd]
%   M: 协方差矩阵 (d×d)
%   返回值: H_n(x) 的值

    % 检查输入维度
    d = length(n);
    assert(length(x) == d, 'x的维度必须与n一致');
    assert(size(M,1) == d && size(M,2) == d, 'M必须是d×d的矩阵');
    
    % 计算目标总阶数
    target_order = sum(n);
    
    % 初始化缓存：H{order}(i) 存储总阶数为order的第i个多项式值
    H = cell(target_order + 1, 1);
    
    % 生成所有阶数从0到target_order的多项式
    for order = 0:target_order
        % 生成当前阶数的所有可能多重指标
        indices = generate_indices(d, order);
        
        % 初始化当前阶数的结果数组
        H{order+1} = zeros(size(indices, 1), 1); % 修正：MATLAB索引从1开始
        
        % 计算每个多重指标对应的多项式值
        for i = 1:size(indices, 1)
            current_n = indices(i, :);
            
            % 零阶多项式处理
            if order == 0
                H{1}(i) = 1; % 修正：H{1}对应阶数0
                continue;
            end
            
            % 选择第一个非零的k
            k = find(current_n > 0, 1, 'first');
            
            % 构造 n - e_k
            n_minus_ek = current_n;
            n_minus_ek(k) = n_minus_ek(k) - 1;
            
            % 查找 H_{n - e_k}(x)
            prev_order = order - 1;
            prev_indices = H{prev_order+1}; % 修正：索引+1
            prev_index = find_index(indices_prev, n_minus_ek);
            H_prev = prev_indices(prev_index);
            
            % 计算主导项
            leading_term = x(k) * H_prev;
            
            % 计算求和项
            sum_term = 0;
            for j = 1:d
                if current_n(j) > 0
                    % 构造 n - e_j
                    n_minus_ej = current_n;
                    n_minus_ej(j) = n_minus_ej(j) - 1;
                    
                    % 查找 H_{n - e_j}(x)
                    prev_index_ej = find_index(indices_prev, n_minus_ej);
                    H_prev_ej = prev_indices(prev_index_ej);
                    
                    % 累加项
                    if j == k
                        sum_term = sum_term + 2 * M(k,j) * (current_n(j)-1) * H_prev_ej;
                    else
                        sum_term = sum_term + 2 * M(k,j) * current_n(j) * H_prev_ej;
                    end
                end
            end
            
            % 存储当前多项式值
            H{order+1}(i) = leading_term + sum_term; % 修正：索引+1
        end
        
        % 保存当前指标集，供下一阶使用
        indices_prev = indices;
    end
    
    % 查找目标指标的结果
    target_index = find_index(indices, n);
    result = H{target_order+1}(target_index); % 修正：索引+1
end

function indices = generate_indices(d, order)
% 生成所有d维、总阶数为order的非负整数多重指标
% 例如：d=2, order=1 -> [1 0; 0 1]

    if order == 0
        indices = zeros(1, d);
        return;
    end
    
    if d == 1
        indices = order;
        return;
    end
    
    indices = [];
    for i = 0:order
        sub_indices = generate_indices(d-1, order-i);
        current_indices = [i * ones(size(sub_indices, 1), 1), sub_indices];
        indices = [indices; current_indices];
    end
end

function idx = find_index(indices, target)
% 在指标矩阵中查找目标指标的位置
    [rows, ~] = size(indices);
    for i = 1:rows
        if all(indices(i, :) == target)
            idx = i;
            return;
        end
    end
    error('未找到目标指标');
end