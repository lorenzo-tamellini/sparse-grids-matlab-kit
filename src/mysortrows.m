function [A,sorter]=mysortrows_new(A,Tol)

[n_row,N] = size(A);

n_blocks = 1;
beg_of_blocks = 1;
end_of_blocks = n_row; % this is supposed to become a vector of n_blocks components
sorter = 1:n_row;

for n = 1:N
    %next_beg_of_blocks = [];
    %next_end_of_blocks = [];
    next_beg_of_blocks = zeros(n_row,1);
    next_end_of_blocks = zeros(n_row,1);
    k = 1;
    r = 1;
    while k <= n_blocks
        b_o_b = beg_of_blocks(k);
        e_o_b = end_of_blocks(k);
        % I could use sortrows here, but that's slower  <---------- need sorter as well
        i = b_o_b:e_o_b;
        [~,ii]=sort(A(i,n));
        A(i,:)=A(i(ii),:);
        sorter(i) = sorter(i(ii));
        j=[1;diff(A(i,n))>Tol;1];
        %j=[1;diff(A(sorter(i),n))>Tol;1];
        %next_beg_of_blocks=[next_beg_of_blocks (b_o_b + find(diff(j)==-1) - 1)'];
        %next_end_of_blocks=[next_end_of_blocks (b_o_b + find(diff(j)==1) -1)'];
        v1 = (b_o_b + find(diff(j)==-1) - 1)';
        v2 = (b_o_b + find(diff(j)==1) -1)';
        %v1 = (beg_of_blocks(k) + find(diff(j)==-1) - 1)';
        %v2 = (beg_of_blocks(k) + find(diff(j)==1) -1)';
        Lv = length(v1); %% equal to length(v2) too
        next_beg_of_blocks(r:r+Lv-1) = v1;
        next_end_of_blocks(r:r+Lv-1) = v2;
        r = r+Lv;
        k=k+1;
    end

    beg_of_blocks=next_beg_of_blocks(1:r-1);
    end_of_blocks=next_end_of_blocks(1:r-1);
    %beg_of_blocks=next_beg_of_blocks;
    %end_of_blocks=next_end_of_blocks;
    n_blocks = length(beg_of_blocks);
    
end
