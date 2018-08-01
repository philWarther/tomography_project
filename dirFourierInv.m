function [f] = dirFourierInv(R, theta, t)
%{
Inverse Radon Transform by Direct Fourier Inversion

R: Sinogram data on meshgrid(t,theta)
t: vector of shift values from discrete radon transform
theta: vector of angle values from discrete radon transform

F: function approximations over meshgrid(x,y)

Steps:
1) 1-D Fourier Coef. R
2) Convert from Polar to Cartesian Coordinates
2.5) Enforce symetry condition?
3) Inverse 2-D Fourier Transform
%}

%initial stuff

if mod(size(R,1),2) == 1
   R = R(1:end-1,:);%using left Reiman sum, so the largest t val row of R wont be used
end
L = (max(t) - min(t))/2;
[N,M] = size(R);
N = (N)/2;

%Step 1

% Find fft of sinogram in polar coords:
%1) shift DC for t to (1,1)
%2) fft-1D
%3) shift DC for r to center
%4) scale by timestep

F_pol = L*fftshift(fft(ifftshift(R,1)),1)/N;
r = -N:N;
r = r*pi/L;

%Step 2
[theta_m, r_m] = meshgrid(theta,r);
[cartx, carty] = pol2cart(theta_m, r_m);
F_pol = [F_pol ; conj(F_pol(1,:))]; %extend F_pol to include its r_n = Npi/L row

[x,y] = meshgrid(r,r);

F_cart = griddata( carty, cartx , F_pol, x, y, 'cubic');

dummy = zeros(2*N+1);%replace nans from griddata with zeros
b = isnan(F_cart);
F_cart(b) = dummy(b);

%Step 2.5
F_cart = (F_cart + conj(rot90(F_cart,2)))/2;%impose forier symmetry
F_cart = ifftshift(F_cart); %shift DC of r back to (1,1)


%Step 3
%ifft
f = N*N/(L*L)*ifft2(F_cart, 'symmetric');
f = fftshift(f); %undo fftshift
f = rot90(f,-1);%fix weird rotation issues
f = fliplr(f);

end

