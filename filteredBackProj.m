function [f] = filteredBackProj(S,theta,t,filter,B)
[n,~]   = size(t);
[m,~]   = size(theta);
N       = (n-1)/2;
L       = max(t);
r       = (-N:N)*pi/L; 

%% Compute phi_hat (the DFT of the apodizing function) for designated filter
if strcmp(filter,'Ram-Lak')
    phi_hat = RamLak(ifftshift(r),B)';
elseif strcmp(filter,'Hanning')
    phi_hat = Hanning(ifftshift(r),B)';
elseif strcmp(filter,'Shepp-Logan')
    phi_hat = SheppLogan(ifftshift(r),B)';
elseif strcmp(filter,'None')
    phi_hat = ones(n,1);
else
    error('Invalid filter selction. Choose from: "Ram-Lak", "Hanning", "Shepp-Logan", or "None"')
end

%% Convolving with filter
S_shifted   = ifftshift(S);         %   Shift the domain to [0,2L) for use 
                                    %   with fft function (by periodic-ness)

RF_tilde    = fft(S_shifted);       %  Compute the DFT of the radon transform

conv        = RF_tilde.*phi_hat;    %   Convolution of RF and phi is equivalent 
                                    %   to this element-wise multiplication of
                                    %   the DFTs.

S_filtered  = fftshift(ifft(conv,'symmetric'));% Computing the filtered sinogram
f = backProj(S_filtered,theta,t);   %   Computing the backprojection of the 
                                    %   filtered data to get the attenution
    
f(f<0) = 0;                         %   Since atteunation constants are nonnegative,
                                    %   any values less than zero are 
                                    %   numerical/aliasing errors and can be     
                                    %   corrected

%% Displaying data

figure(1);
subplot(3,1,1);plot(t,S(:,1)); title('Initial RF value for \Theta_0')
subplot(3,1,2);plot(r, fftshift(phi_hat)); title([filter ' filter'])
subplot(3,1,3);plot(t,S_filtered(:,1)); title('Filtered RF for \Theta_0')

end
%% Filter Functions

function [phi_hat] = RamLak(r, B)
% Applies a Ram-Lak Filter to a vector of frequencies, r, with bandlimit B
phi_hat = zeros(size(r));
idx = find(abs(r)>B);
phi_hat = abs(r);
phi_hat(idx) = 0;

end

function [phi_hat] = Hanning(r,B)
% Applies a Hanning Filter to a vector of frequencies, r, with bandlimit B
phi_hat         = zeros(size(r));
idx             = find(abs(r)>B);
phi_hat         = abs(r).*(cos((pi.*r)/(2*B)).^2);
phi_hat(idx)    = 0;

end

function [phi_hat] = SheppLogan(r,B)
% Applies B-low pass Shepp-Logan filter over the frequency domain r
phi_hat         = zeros(size(r));
d = pi/B;
phi_hat         = abs(r).*abs(sinc(r*d/2)).^3;
end

%% Backprojection Function
function [A] = backProj(S, theta, t)

% This function performs Back Projection on a sinogram and
% reconstructs the image of a two-dimensional region with non-constant 
% attenuation coefficients. 
%
% ----------INPUTS-------------
% S:        A M (rows) x N (columns) array corresponding to a sinogram with
%           N theta measurement
% theta:    An N X 1 array of the measurement angles (domain of sinogram),
%           [0,pi-e] **NOTE: the measurment angles are prependicular to the
%           incoming rays, so the computed A matrix will need to be rotated
%           90 degrees for the final image** 
% t :       An M X 1 array of the distance values (range of the sinogram)
%
% ----------OUTPUTS -----------
% A:        M X M matrix of attenutation coefficients associated with the
%           two-dimensional region
 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [1]   This function will apply the values of the radon transform data to 
%       all entries of the matrix temp matrix corresponding to the line 
%       used to compute the the transform values. 
%
% [2]   To do this for all measurments, the temp matrix will be rotated
%       according to the values of the measurement angles.
%
% [3]   The temp matrix values will be accumulated to the output A
%
% [4]   To attain the attenuation coefficients at each discrete entry in
%       the region, the values of A are averaged over the N measurements
%
% [5]   Since our measurement angles in theta are perpendicular to the
%       incoming waves, the output A is rotated 90 degrees to match the 
%       measured region
%
%temp    = S(:,1);
S       = [S flipud(S)]; % simulate measurements from 0 to 2*pi
theta   = [theta; theta+pi];
[M,N]   = size(S);
A_old   = zeros(M);
A       = zeros(M);
%for j = 1:2 
figure(2)
for i = 1:N
    temp = repmat(S(:,i),[1,M]);                     %[1]
    temp = imrotate(temp,theta(i)*(180/pi)+90,'bilinear','crop'); %[2] [5]
%       Uncomment/comment below to display graphics    
%     subplot(1,2,1);imagesc(flipud(A_old));colormap('gray'); 
%     axis('equal'); axis off;
%     subplot(1,2,2);spy(temp);
%     %subplot(1,2,2);surf(t,t,flipud(A_old),'edgecolor','none');
%     drawnow;
    
    A       = A_old + temp;                                 %[3]
    if i > 1 % zero out unaffected entries
        A(A==A_old) = 0;
        A(A==temp)  = 0;
    end 
    A_old = A;
end
A   = (1/(2*N))*A;   %[4]
A   = flipud(A);
% Displaying the results


% a = num2str(t(1)); b = num2str(t(end));
% figure(3);
% imagesc(f;colormap('gray');axis('equal')
% axis off;
% title(['Attenuation Coefficients over [' a ' , ' b ']\times[' a ' , ' b ']']);
% figure(4);
% subplot(1,2,1);surf(t,t,f,'edgecolor','none');
% subplot(1,2,2);surf(t
end 