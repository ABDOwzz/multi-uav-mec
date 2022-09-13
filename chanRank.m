% This function ranks the user in ascending order of their channel gains.
function [gRank, rank, Psi] = chanRank(g, U, N)
    %  g is a U*N matrix denoting the channel gains for all links throughout the interval.
    gRank = zeros(U, N);
    rank = zeros(U, N);
    Psi = zeros(U, N);
    for n = 1:N
        [gRank(:, n), rank(:, n)] = sort(g(:, n));
        % rank gives a list of user index. Its value denotes the user index with that ranking, 
        % i.e. rank(k, n) means the index of user with the k-th smallest channel gain at time n.
        [B2, Psi(:, n)] = sort(rank(:, n));
        % Psi helps track the rankings of users, 
        % i.e. Psi(i, n) denotes the ranking of user i at time n.
    end
end