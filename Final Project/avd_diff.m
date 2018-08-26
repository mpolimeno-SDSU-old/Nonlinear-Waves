function dwdtF = avd_diff(t,w,RK,sf,v,Dx,Dy,Dx2,Dy2)

dx2w = real(ifft2(reshape((Dx2.*w).',RK,RK).'));
dy2w = real(ifft2(reshape((Dy2.*w).',RK,RK).'));
dxsf = real(ifft2(reshape((Dx.*sf).',RK,RK).'));
dyw = real(ifft2(reshape((Dy.*w).',RK,RK).'));
dysf = real(ifft2(reshape((Dy.*sf).',RK,RK).'));
dxw = real(ifft2(reshape((Dx.*w).',RK,RK).'));

% dwdt = v*(Dx2+Dy2).*w - (((Dx.*sf).*(Dy.*w))-((Dy.*sf).*(Dx.*w)));

dwdt = v*(dx2w+dy2w) - (dxsf.*dyw) + (dysf.*dxw);

dwdtF = reshape(fft2(dwdt)',length(Dx),1);

end