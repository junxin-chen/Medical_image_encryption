% ����ͼ�����������
% img Ϊ����ͼ�����ȡn�����ص�
function img_cov(img,n)

% M,N ���ͼ���С
[M,N]=size(img);

% gray_val_base,gray_val_hor,gray_val_ver,gray_val_dia�ֱ�洢���㣬ˮƽ���ڵ㣬��ֱ���ڵ㣬б���ڵ�����ػҶ�ֵ
gray_val_base = zeros(1,n);
gray_val_hor = zeros(1,n);
gray_val_ver = zeros(1,n);
gray_val_dia = zeros(1,n);

% ��������ԣ����ȡn�Ե㣬ˮƽ���ڵ��Ĺ�ϵΪ(x,y)-(x+1,y)����ֱ���ڵ��Ĺ�ϵΪ(x,y)-(x,y+1)��б���ڵ��Ĺ�ϵΪ(x
% ,y)-(x+1,y+1)
% ѭ�� n ��
for i = 1:n
    % ����һ�������(x,y)����
    x_coor = randint(1,1,[1,N-1]); % ȡN-1��M-1,�������ڵ㲻�ᳬ��
    y_coor = randint(1,1,[1,M-1]);
    % ȡ�Ҷ�ֵ��������Ӧ������
    gray_val_base(i) = img(x_coor,y_coor);
    gray_val_hor(i) = img(x_coor+1,y_coor);
    gray_val_ver(i) = img(x_coor,y_coor+1);
    gray_val_dia(i) = img(x_coor+1,y_coor+1);
end
% ������ضȹ�ʽ����
sum_base = 0;
sum_hor = 0;
sum_ver = 0;
sum_dia = 0;
for i=1:n
    sum_base = sum_base + gray_val_base(i);
    sum_hor = sum_hor + gray_val_hor(i);
    sum_ver = sum_ver + gray_val_ver(i);
    sum_dia = sum_dia + gray_val_dia(i);
end
mean_base = sum_base/n;
mean_hor = sum_hor/n;
mean_ver = sum_ver/n;
mean_dia = sum_dia/n;

sum_base_var = 0;
sum_hor_var = 0;
sum_ver_var = 0;
sum_dia_var = 0;
for i=1:n
    sum_base_var = sum_base_var + (gray_val_base(i) - mean_base)^2;
    sum_hor_var = sum_hor_var + (gray_val_hor(i) - mean_hor)^2;
    sum_ver_var = sum_ver_var + (gray_val_ver(i) - mean_ver)^2;
    sum_dia_var = sum_dia_var + (gray_val_dia(i) - mean_dia)^2;
end
variance_base = sum_base_var/n;
variance_hor = sum_hor_var/n;
variance_ver = sum_ver_var/n;
variance_dia = sum_dia_var/n;

sum_cov_base_hor = 0;
sum_cov_base_ver = 0;
sum_cov_base_dia = 0;
for i=1:n
    sum_cov_base_hor = sum_cov_base_hor + (gray_val_base(i) - mean_base)*(gray_val_hor(i) - mean_hor);
    sum_cov_base_ver = sum_cov_base_ver + (gray_val_base(i) - mean_base)*(gray_val_ver(i) - mean_ver);
    sum_cov_base_dia = sum_cov_base_dia + (gray_val_base(i) - mean_base)*(gray_val_dia(i) - mean_dia);
end
cov_base_hor = sum_cov_base_hor/n;
cov_base_ver = sum_cov_base_ver/n;
cov_base_dia = sum_cov_base_dia/n;

r_base_hor = cov_base_hor/(sqrt(variance_base)*sqrt(variance_hor));
r_base_ver = cov_base_ver/(sqrt(variance_base)*sqrt(variance_ver));
r_base_dia = cov_base_dia/(sqrt(variance_base)*sqrt(variance_dia));

% ������
digits(6)
vpa(r_base_hor)
vpa(r_base_ver)
vpa(r_base_dia)

%��ͼ
figure;
plot(gray_val_base,gray_val_hor,'.','Markersize',3);
axis([0,255,0,255]);
figure;
plot(gray_val_base,gray_val_ver,'.','Markersize',3);
axis([0,255,0,255]);
figure;
plot(gray_val_base,gray_val_dia,'.','Markersize',3);
axis([0,255,0,255]);






 