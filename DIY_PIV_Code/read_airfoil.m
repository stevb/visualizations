function [x,y]=read_airfoil(filename)

fid=fopen([filename],'r');

for i=1:2
    line=fgetl(fid);
end

[A,count]=fscanf(fid,'%g',[2,inf]);

x = A(1,:);
y = A(2,:);
