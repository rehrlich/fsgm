function like = likelyhood( params )

global S K histogram;

% S = 17;
% K = 6;
% histogram = [ 104 576 329 93 102 68 45 20 40 32 16 25 19 41 50 81 179 1454 ];


condP = zeros(S+1,K); 
for i = 0:S
    for k=1:K
        
        condP(i+1,k) = nchoosek(S,i) * params(K+k)^i * ( 1-params(K+k) )^(S-i);
        
    end
end


like = 0;

for i = 0:S
    temp = 0;
    for k = 1:K
        temp = temp + params(k)*condP(i+1,k);
    end
    like = like + histogram(i+1) * log(temp);
end

% we need to maximize using the minimize function, so take negative.
if ( imag(like) ~= 0 )
    params
    histogram(1)
end
    
like = -real(like);