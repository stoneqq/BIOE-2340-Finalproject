%% calculate the final results
% each row is an image, the images included are 1~8, 13 and 14, the first column is the cell soma count, the second column is the total
% neurite length, the first five images are tunning set and the last five
% images are testing set.

% declare the final raw results
Manual_results={
15	8753
24	11297.409
24	5255.865
28	11139.504
29	11138.214
39	11145.665
28	3909.534
28	10240.888
13	10948.691
16	7316.06
};

NeuriteTracer_results={
16	9889.574
27	14034.603
27	6930.561
33	12468.42
37	12124.986
44	14565.284
35	5742.697
28	11413.247
14	10958.938
17	11177.011
};

GAIN_results={
16	3661.75477
25	5380.43724
23	3973.453262
28	6219.333773
33	6105.610912
41	10376.27441
33	4722.060471
30	7026.147049
12	3044.243007
17	2449.596254
};

GAIN_Hessian_results={
16	3751.565951
25	5652.674829
23	4058.578151
28	5951.711407
33	6106.865748
41	10600.6003
33	4797.374178
29	6613.488107
13	3032.314075
17	2422.424681
};

%% calculate the normalization factor and normalize the results to manual results
norm_factor_NeuriteTracer=zeros(1,2);
norm_factor_GAIN=zeros(1,2);
norm_factor_GAIN_Hessian=zeros(1,2);

norm_factor_NeuriteTracer(1)=mean(cell2mat(NeuriteTracer_results(1:5,1))./cell2mat(Manual_results(1:5,1)));
norm_factor_NeuriteTracer(2)=mean(cell2mat(NeuriteTracer_results(1:5,2))./cell2mat(Manual_results(1:5,2)));
norm_factor_GAIN(1)=mean(cell2mat(GAIN_results(1:5,1))./cell2mat(Manual_results(1:5,1)));
norm_factor_GAIN(2)=mean(cell2mat(GAIN_results(1:5,2))./cell2mat(Manual_results(1:5,2)));
norm_factor_GAIN_Hessian(1)=mean(cell2mat(GAIN_Hessian_results(1:5,1))./cell2mat(Manual_results(1:5,1)));
norm_factor_GAIN_Hessian(2)=mean(cell2mat(GAIN_Hessian_results(1:5,2))./cell2mat(Manual_results(1:5,2)));

NeuriteTracer_norm=[cell2mat(NeuriteTracer_results(:,1))/norm_factor_NeuriteTracer(1),cell2mat(NeuriteTracer_results(:,2))/norm_factor_NeuriteTracer(2)];
GAIN_norm=[cell2mat(GAIN_results(:,1))/norm_factor_GAIN(1),cell2mat(GAIN_results(:,2))/norm_factor_GAIN(2)];
GAIN_Hessian_norm=[cell2mat(GAIN_Hessian_results(:,1))/norm_factor_GAIN_Hessian(1),cell2mat(GAIN_Hessian_results(:,2))/norm_factor_GAIN_Hessian(2)];

NeuriteTracer_ssd=(cell2mat(NeuriteTracer_results)-cell2mat(Manual_results)).^2;
GAIN_ssd=(cell2mat(GAIN_results)-cell2mat(Manual_results)).^2;
GAIN_Hessian_ssd=(cell2mat(GAIN_Hessian_results)-cell2mat(Manual_results)).^2;

NeuriteTracer_ssd_sum=[sum(NeuriteTracer_ssd(1:5,1)),sum(NeuriteTracer_ssd(1:5,2));sum(NeuriteTracer_ssd(6:10,1)),sum(NeuriteTracer_ssd(6:10,2))];
GAIN_ssd_sum=[sum(GAIN_ssd(1:5,1)),sum(GAIN_ssd(1:5,2));sum(GAIN_ssd(6:10,1)),sum(GAIN_ssd(6:10,2))];
GAIN_Hessian_ssd_sum=[sum(GAIN_Hessian_ssd(1:5,1)),sum(GAIN_Hessian_ssd(1:5,2));sum(GAIN_Hessian_ssd(6:10,1)),sum(GAIN_Hessian_ssd(6:10,2))];

NeuriteTracer_norm_ssd=(NeuriteTracer_norm-cell2mat(Manual_results)).^2;
GAIN_norm_ssd=(GAIN_norm-cell2mat(Manual_results)).^2;
GAIN_Hessian_norm_ssd=(GAIN_Hessian_norm-cell2mat(Manual_results)).^2;

NeuriteTracer_norm_ssd_sum=[sum(NeuriteTracer_norm_ssd(1:5,1)),sum(NeuriteTracer_norm_ssd(1:5,2));sum(NeuriteTracer_norm_ssd(6:10,1)),sum(NeuriteTracer_norm_ssd(6:10,2))];
GAIN_norm_ssd_sum=[sum(GAIN_norm_ssd(1:5,1)),sum(GAIN_norm_ssd(1:5,2));sum(GAIN_norm_ssd(6:10,1)),sum(GAIN_norm_ssd(6:10,2))];
GAIN_Hessian_norm_ssd_sum=[sum(GAIN_Hessian_norm_ssd(1:5,1)),sum(GAIN_Hessian_norm_ssd(1:5,2));sum(GAIN_Hessian_norm_ssd(6:10,1)),sum(GAIN_Hessian_norm_ssd(6:10,2))];

NeuriteTracer_norm_ssd_sum_percent=NeuriteTracer_norm_ssd_sum(2,:)./NeuriteTracer_norm_ssd_sum(1,:);
GAIN_norm_ssd_sum_percent=GAIN_norm_ssd_sum(2,:)./GAIN_norm_ssd_sum(1,:);
GAIN_Hessian_norm_ssd_sum_percent=GAIN_Hessian_norm_ssd_sum(2,:)./GAIN_Hessian_norm_ssd_sum(1,:);
GAIN_Hessian_compare=GAIN_Hessian_norm_ssd_sum./GAIN_norm_ssd_sum;
% the final figures are plotted in GRAPHPAD PRISM
