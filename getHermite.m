function [out, norm] = getHermite(Va, a, n, T)
if nargin < 4
    T = 1;
end
ac = a;
% nbar不等于零：对应TMST

global nbar
nbar =0;

% 适用于PPC： 计算衰减后的Va
r = acosh(Va)/2;
l = tanh(r);

% l:更新过后的lambda
l = l*sqrt(T);
r = atanh(l);
Va = cosh(2*r);

% 后面一个1/2:如果是+则说明是PS，-则说明是PA
A = -(2*nbar+1)*Va/2 - 1/2;
B = A;

% C 若取反：g = -1的抵消作用
C = (2*nbar+1)*sqrt(Va^2-1)/2;

% A = -Va/2 + 1/2;
% B = A;
% C = sqrt(Va^2-1)/2;

k_out = zeros(1, numel(a));

for n5 = 0:max(n)
    for n6 = 0:max(n)
        for n7 = 0:max(n)
            for n8 = 0:max(n)
                if (n(1)-n5-n7)>=0 && (n(2)-n5-n8)>=0 && (n(3)-n6-n7)>=0 && (n(4)-n6-n8)>=0
                    k_out = k_out + ...
                        a.^(n(2)+n(3)-n5-n6-n7-n8).*...
                        ac.^(n(1)+n(4)-n5-n6-n7-n8)...
                        * (A + C) ^ (sum(n)-2*(n5+n6+n7+n8)) ...
                        * A^n5 * B^n6 * C^ (n7 + n8)...
                        * prod([(n5+1):n(1), (n6+1):n(2), (n7+1):n(3), (n8+1):n(4)]) ...
                        /factorial(n(1)-n5-n7)/factorial(n(2)-n5-n8)/factorial(n(3)-n6-n7)/factorial(n(4)-n6-n8);
                end
            end
        end
    end
end

% for n5 = 0:max(n)
%     for n6 = 0:max(n)
%         for n7 = 0:max(n)
%             for n8 = 0:max(n)
%                 if (n(1)-n5-n7)>=0 && (n(2)-n5-n8)>=0 && (n(3)-n6-n7)>=0 && (n(4)-n6-n8)>=0
%                     k_out = k_out + a.^(n(2)+n(3)-n5-n6-n7-n8).*ac.^(n(1)+n(4)-n5-n6-n7-n8)...
%                         * (A + C) ^ (sum(n)-2*(n5+n6+n7+n8)) ...
%                         * A^n5 * B^n6 * C^ (n7 + n8)...
%                         * exp(gammaln(n(1)+1)-gammaln(n5+1))*exp(gammaln(n(2)+1)-gammaln(n6+1))...
%                         * exp(gammaln(n(3)+1)-gammaln(n7+1))*exp(gammaln(n(4)+1)-gammaln(n8+1))...
%                         /gamma(n(1)-n5-n7+1)/gamma(n(2)-n5-n8+1)/gamma(n(3)-n6-n7+1)/gamma(n(4)-n6-n8+1);
%                 end
%             end
%         end
%     end
% end

out = k_out;
norm = k_out(1);

