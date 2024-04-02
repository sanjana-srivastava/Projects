function [] = DesignVariable_Final()

global u

Data_1 = {'x(1)', 'x(2)', 'x(3)', 'x(4)', 'x(5)', 'x(6)', 'x(7)','x(8)','x(9)','x(10)','x(11)','x(12)',...
          'x(13)', 'x(14)', 'x(15)', 'x(16)'};
Data_2 = u';
fid = fopen('DesignVariables.txt','w+t');
if fid < 0
    fprintf('error opening file\n')
    return;
end

for i = 1:length(Data_1)
    fprintf(fid,'%s = '  ,Data_1{i});
    fprintf(fid,'%1.5f\n',Data_2(i,:));
end
fclose(fid);

