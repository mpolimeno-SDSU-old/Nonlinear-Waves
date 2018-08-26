function dwdt = adv_diff(t,wF,dummy,KT,v,Dx,Dy,Dx2,Dy2)

sfF = wF./(Dx2+Dy2);

psix = real(ifft2(reshape((1i*Dx.*sfF).',KT,KT).'));
psiy = real(ifft2(reshape((1i*Dy.*sfF).',KT,KT).'));
wy = real(ifft2(reshape((1i*Dy.*wF).',KT,KT).'));
wx = real(ifft2(reshape((1i*Dx.*wF).',KT,KT).'));

% dwdt = v*(Dx2+Dy2).*w - (((Dx.*sf).*(Dy.*w))-((Dy.*sf).*(Dx.*w)));

dwdt = -v*(Dx2+Dy2).*wF+reshape(fft2(wx.*psiy-wy.*psix)',KT^2,1);

% dwdtF = reshape(fft2(dwdt)',length(Dx),1);

end