%% clear
clc; clear all; close all;

%% load image and calculate vesselness
sigma_range=[0.01,0.1,0.2,0.5,0.8,1,1.5,2,2.5,3,3.5,4,5,10];
gamma_range=3;

im=cell(1,length(sigma_range)*length(gamma_range));
im_c=cell(1,length(sigma_range)*length(gamma_range));

histeq_im=cell(1,length(sigma_range)*length(gamma_range));
histeq_im_c=cell(1,length(sigma_range)*length(gamma_range));
% im{1} = imread ('jellyfish.png');
im{1} = imread('stackedImage_2.tif');
im_c{1} = imread('stackedImage_2.tif');
histeq_im{1} = imread('stackedImage_2.tif');
histeq_im_c{1} = imread('stackedImage_2.tif');
nuclear=imread('stackedImage_2.tif',2);
for sigma=1:length(sigma_range)
    for gamma = 1:length(gamma_range)
    [im_c{(sigma-1)*length(gamma_range)+gamma+1},vx,vy] = neuriteness2d(imcomplement(im_c{1}),sigma_range(sigma),gamma_range(gamma));
    [im{(sigma-1)*length(gamma_range)+gamma+1},vx,vy] = neuriteness2d((im{1}),sigma_range(sigma),gamma_range(gamma));
    [histeq_im_c{(sigma-1)*length(gamma_range)+gamma+1},vx,vy] = neuriteness2d(imcomplement(histeq_im_c{1}),sigma_range(sigma),gamma_range(gamma));
    [histeq_im{(sigma-1)*length(gamma_range)+gamma+1},vx,vy] = neuriteness2d((histeq_im{1}),sigma_range(sigma),gamma_range(gamma));
    end
end
% gamma doesn't really do a lot, the numerical difference results are on
% the level of 1e-17, which is negligable, choose gamma=3.
%% plot
figure(1)
for sigma=1:length(sigma_range)
    for gamma = 1:length(gamma_range)
%         figure(1)
%         subplot(length(sigma_range),length(gamma_range),(sigma-1)*length(gamma_range)+gamma)
%         imagesc(im_c{(sigma-1)*length(gamma_range)+gamma+1}); colormap gray; 
%         set(gca,'ytick',[]); set(gca,'xtick',[]); axis image; axis tight;
%         title(['s=', num2str(sigma_range(sigma)),', g=',num2str(gamma_range(gamma))])
%         fig_handle=figure;
%         imagesc(im{(sigma-1)*length(gamma_range)+gamma+1}); colormap gray; 
%         set(gca,'ytick',[]); set(gca,'xtick',[]); axis image; axis tight;
%         title(['s=', num2str(sigma_range(sigma)),', g=',num2str(gamma_range(gamma))])
        savename=['s', num2str(sigma),'g',num2str(3)];
        savename_series={['im', num2str(sigma*2-1)],['im', num2str(sigma*2)]};
%         print(fig_handle,savename,'-dtiffn')
%         imwrite(im{(sigma-1)*length(gamma_range)+gamma+1},[savename,'.tif'],'tif')
%         imwrite(im_c{(sigma-1)*length(gamma_range)+gamma+1},[savename,'_complement.tif'],'tif')
%         imwrite(histeq_im{(sigma-1)*length(gamma_range)+gamma+1},['histeq_',savename,'.tif'],'tif')
%         imwrite(histeq_im_c{(sigma-1)*length(gamma_range)+gamma+1},['histeq_',savename,'_complement.tif'],'tif')
%         imwrite(nuclear,[savename,'.tif'],'tif')
%         imwrite(nuclear,[savename,'_complement.tif'],'tif')
%         imwrite(nuclear,['histeq_',savename,'.tif'],'tif')
%         imwrite(nuclear,['histeq_',savename,'_complement.tif'],'tif')   
        
        imwrite(im{(sigma-1)*length(gamma_range)+gamma+1},[savename_series{1},'.tif'],'tif')
        imwrite(im_c{(sigma-1)*length(gamma_range)+gamma+1},['complement',savename_series{1},'.tif'],'tif')
        imwrite(histeq_im{(sigma-1)*length(gamma_range)+gamma+1},['histeq_',savename_series{1},'.tif'],'tif')
        imwrite(histeq_im_c{(sigma-1)*length(gamma_range)+gamma+1},['complement','histeq_',savename_series{1},'.tif'],'tif')
        imwrite(nuclear,[savename_series{2},'.tif'],'tif')
        imwrite(nuclear,['complement',savename_series{2},'.tif'],'tif')
        imwrite(nuclear,['histeq_',savename_series{2},'.tif'],'tif')
        imwrite(nuclear,['complement','histeq_',savename_series{2},'.tif'],'tif')  
    end
end

%% Sigma range from 0.01~0.2 works fine, adapthisteq does not make a difference, apply to all other images
tic
sigma=0.1;
gamma=3;

im=cell(1,40);
nuclear=cell(1,40);
neurite=cell(1,40);
[files_list,path]=uigetfile('*.*','Image File', 'MultiSelect','on');
for i=1:40
%     im{i} = imread(['stackedImage_',num2str(i),'.tif']);
%     nuclear{i}=imread(['stackedImage_',num2str(i),'.tif'],2);
    im{i} = imread([path,'stackedImage_',num2str(i),'.tif']);
    nuclear{i}=imread([path,'stackedImage_',num2str(i),'.tif'],2);
    [neurite{i},vx,vy] = neuriteness2d((im{i}),sigma,gamma);
    savename_series={['im_sig01gamma3_', num2str(i*2-1)],['im_sig01gamma3_', num2str(i*2)]};
    imwrite(neurite{i},[savename_series{1},'.tif'],'tif')
    imwrite(nuclear{i},[savename_series{2},'.tif'],'tif')
end

toc

