function []=Constellation_plot(input)
DataLen = length(input);
xyz(:,1)=real(input);
xyz(:,2)=imag(input);
z3 = round(rand(1,DataLen)*4+1);
xyz(:,3)=reshape(z3,DataLen,1);
xyz(:,3)=xyz(:,3)-xyz(:,3);
sizeXYZ = size(xyz);
searchR = 0.2;%����Բ�뾶
for i=1:sizeXYZ(1)
    index_i = find(xyz(:,1)>(xyz(i,1)-searchR) & xyz(:,1)<xyz(i,1)+searchR & xyz(:,2) > xyz(i,2) -searchR & xyz(:,2)<xyz(i,2)+searchR );
    sizeIndexI = size(index_i);
    xyz(i,3)=sizeIndexI(1); % �ҵ�����Χ�뾶��0.2�ĵ�ĸ�����Խ������ɫ��RGBֵԽ��
end
[sortXYZ,sortI]=sort(xyz(:,3));
sz = 5;% ԲȦ�Ĵ�С
scatter(xyz(:,1),xyz(:,2),sz,xyz(:,3),'filled');%% filled��ʾԲȦΪʵ��
% scatter(xyz(sortI,1),xyz(sortI,2),sz,xyz(sortI,3),'filled');%% filled��ʾԲȦΪʵ��
hold on ;
    






end