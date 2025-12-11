function const = constants_ae4205()
% constants_ae4205
% Reads constants from a CSV file and returns a struct.

    % Load the CSV file as a table
    T = readtable('constants_ae4205.csv', 'Delimiter', ',', 'ReadVariableNames', true);

    % Convert table to struct: Name column becomes field names, Value column becomes field values
    const = struct();

    for i = 1:height(T)
        fieldName = string(T.Name{i});
        const.(fieldName) = T.Value(i);
    end

end
