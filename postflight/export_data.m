clear
clc

load("FTISxprt-20200310_flight2.mat");
names = strings([1, 49]);
data = zeros(size(flightdata.time.data, 2), 49);

fn = fieldnames(flightdata);
for k=1:numel(fn)
   cell = flightdata.(fn{k});
   names(k) = append(cell.description, ' [', cell.units, ']');
   data(:,k) = cell.data;
end

output = 'flightdata.csv';
writematrix(names, output);
writematrix(data, output, 'WriteMode', 'append');


