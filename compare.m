%% test data - proj2compare.mat

clear; close all

load proj2compare.mat
%% plot initial figure and sonogram

figure
subplot(1,2,1)
imagesc(A)
title('Test Image 1')
axis off
subplot(1,2,2)
imagesc(S2)
title('Sonogram S2')
axis off
colormap('Gray')
%% Compute and plot DFI of proj2compare.mat data

DFI_S1 = dirFourierInv(S1,theta1,t1);
DFI_S2 = dirFourierInv(S2,theta2,t2);



figure
subplot(2,2,1)
imagesc(DFI_S1)
title('DFI of S1')
axis off
subplot(2,2,2)
imagesc(DFI_S2)
title('DFI of S2')
axis off
subplot(2,2,3)
imagesc(1-DFI_S1)
title('Negative Image of DFI of S1')
axis off
subplot(2,2,4)
imagesc(1-DFI_S2)
title('Negative Image of DFI of S2')
axis off
colormap('Gray')

L1_DFI_S1 = norm(A - DFI_S1,1);
Linf_DFI_S1 = max(abs(A(:) - DFI_S1(:)));
L1_DFI_S2 = norm(A - DFI_S2,1);
Linf_DFI_S2 = max(abs(A(:) - DFI_S2(:)));
%% Compute and Plot FBP of proj2compare.mat data for different B's

for i = 1:3
    frl1(:,:,i) = filteredBackProj(S1,theta1,t1,'Ram-Lak',20*i);
    fh1(:,:,i) = filteredBackProj(S1,theta1,t1,'Hanning',30*i);
    fsl1(:,:,i) = filteredBackProj(S1,theta1,t1,'Shepp-Logan',60*i*pi);
end

figure
montage(frl1)
title('FBP on S1, Ram-Lak filter, B = 20:20:180')
figure
montage(fh1)
title('FBP on S1, Hanning filter, B = 30:30:270')
figure
montage(fsl1)
title('FBP on S1, Shepp Logan filter, B = 60\pi :60\pi :540\pi')

for i = 1:9
    frl2(:,:,i) = filteredBackProj(S2,theta2,t2,'Ram-Lak',20*i);
    fh2(:,:,i) = filteredBackProj(S2,theta2,t2,'Hanning',30*i);
    fsl2(:,:,i) = filteredBackProj(S2,theta2,t2,'Shepp-Logan',60*i*pi);
end

figure
montage(frl1)
title('FBP on S2, Ram-Lak filter, B = 20:20:180')
figure
montage(fh1)
title('FBP on S2, Hanning filter, B = 30:30:270')
figure
montage(fsl1)
title('FBP on S2, Shepp Logan filter, B = 60\pi :60\pi :540\pi')
%% Compare best FBP stuff wit DFI stuff

FBP_RL_S1 = frl1(:,:,end);
FBP_RL_S2 = frl2(:,:,end);
FBP_H_S1 = fh1(:,:,end);
FBP_H_S2 = fh2(:,:,end);
FBP_SL_S1 = fsl1(:,:,end);
FBP_SL_S2 = fsl2(:,:,end);
%% compare L-1 and L-inf norms


L1_FBP_RL_S1 = norm(A - FBP_RL_S1,1);
Linf_FBP_RL_S1 = max(abs(A(:) - FBP_RL_S1(:)));

L1_FBP_RL_S2 = norm(A - FBP_RL_S2,1);
Linf_FBP_RL_S2 = max(abs(A(:) - FBP_RL_S2(:)));

L1_FBP_H_S1 = norm(A - FBP_H_S1,1);
Linf_FBP_H_S1 = max(abs(A(:) - FBP_H_S1(:)));

L1_FBP_H_S2 = norm(A - FBP_H_S2,1);
Linf_FBP_H_S2 = max(abs(A(:) - FBP_H_S2(:)));

L1_FBP_SL_S1 = norm(A - FBP_SL_S1,1);
Linf_FBP_SL_S1 = max(abs(A(:) - FBP_SL_S1(:)));

L1_FBP_SL_S2 = norm(A - FBP_SL_S2,1);
Linf_FBP_SL_S2 = max(abs(A(:) - FBP_SL_S2(:)));



figure
subplot(2,2,1)
surf1 = surf(abs(A - DFI_S1));
surf1.EdgeColor = 'none';
title(['DFI Error: L^1 = ', num2str(L1_DFI_S1), ', L^\infty = ', num2str(Linf_DFI_S1)])
subplot(2,2,2)
surf2 = surf(abs(A - FBP_RL_S1));
surf2.EdgeColor = 'none';
title(['FBP RL Error: L^1 = ', num2str(L1_FBP_RL_S1), ', L^\infty = ', num2str(Linf_FBP_RL_S1)])
subplot(2,2,3)
surf3 = surf(abs(A - FBP_H_S1));
surf3.EdgeColor = 'none';
title(['FBP H Error: L^1 = ', num2str(L1_FBP_H_S1), ', L^\infty = ', num2str(Linf_FBP_H_S1)])
subplot(2,2,4)
surf4 = surf(abs(A - FBP_SL_S1));
surf4.EdgeColor = 'none';
title(['FBP SL Error: L^1 = ', num2str(L1_FBP_SL_S1), ', L^\infty = ', num2str(Linf_FBP_SL_S1)])

figure
subplot(2,2,1)
surf5 = surf(abs(A - DFI_S2));
surf5.EdgeColor = 'none';
title(['DFI Error: L^1 = ', num2str(L1_DFI_S2), ', L^\infty = ', num2str(Linf_DFI_S2)])
subplot(2,2,2)
surf6 = surf(abs(A - FBP_RL_S2));
surf6.EdgeColor = 'none';
title(['FBP RL Error: L^1 = ', num2str(L1_FBP_RL_S2), ', L^\infty = ', num2str(Linf_FBP_RL_S2)])
subplot(2,2,3)
surf7 = surf(abs(A - FBP_H_S2));
surf7.EdgeColor = 'none';
title(['FBP H Error: L^1 = ', num2str(L1_FBP_H_S2), ', L^\infty = ', num2str(Linf_FBP_H_S2)])
subplot(2,2,4)
surf8 = surf(abs(A - FBP_SL_S2));
surf8.EdgeColor = 'none';
title(['FBP SL Error: L^1 = ', num2str(L1_FBP_SL_S2), ', L^\infty = ', num2str(Linf_FBP_SL_S2)])




