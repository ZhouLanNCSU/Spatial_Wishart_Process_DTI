function AI = inverse_lower(A)
    len = length(A);
    I  = eye(len);
    M  = [A I];
    for row = 1:len
        M(row,:) = M(row,:)/M(row,row);
        for idx = 1:row-1
            M(row,:) = M(row,:) - M(idx,:)*M(row,idx);
        end
    end
    AI = M(:,len+1:end);
end