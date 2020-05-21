function sig = vanGelderen(bigDel,litDel,Gmax,R,D)
% approximation to DWI signal inside a cylinder as published by P Van
% Gelderen et al JMRB 1994

%nuclear gyromagnetic ratio in rad/G*ms, 42.576 MHz/T or 267.513 x106 rad s-1 T-1
gamma = 0.0267513; 

% the first 200 or so bessel derivative zeros, as output by 
% BessDerivZerosBisect2(1,1:200,1e-12); (found on file exchange)
alphaR = [1.841 5.331 8.536 11.706 14.864 18.016 21.164 24.311 27.457 ...
    30.602 33.746 36.890 40.033 43.177 46.320 49.462 52.605 55.748 ...
    58.890 62.032 65.175 68.317 71.459 74.601 77.743 80.885 84.027 ...
    87.169 90.311 93.453 96.595 99.737 102.879 106.020 109.162 112.304 ...
    115.446 118.588 121.730 124.871 128.013 131.155 134.297 137.438 ...
    140.580 143.722 146.863 150.005 153.147 156.289 159.430 162.572 ...
    165.714 168.855 171.997 175.139 178.280 181.422 184.564 187.705 ...
    190.847 193.989 197.131 200.272 203.414 206.555 209.697 212.839 ...
    215.980 219.122 222.264 225.405 228.547 231.689 234.830 237.972 ...
    241.114 244.255 247.397 250.539 253.680 256.822 259.963 263.105 ...
    266.247 269.388 272.530 275.672 278.813 281.955 285.096 288.238 ...
    291.380 294.521 297.663 300.805 303.946 307.088 310.229 313.371 ...
    316.513 319.654 322.796 325.938 329.079 332.221 335.362 338.504 ...
    341.646 344.787 347.929 351.070 354.212 357.354 360.495 363.637 ...
    366.779 369.920 373.062 376.203 379.345 382.487 385.628 388.770 ...
    391.911 395.053 398.195 401.336 404.478 407.620 410.761 413.903 ...
    417.044 420.186 423.328 426.469 429.611 432.752 435.894 439.036 ...
    442.177 445.319 448.460 451.602 454.744 457.885 461.027 464.168 ...
    467.310 470.452 473.593 476.735 479.876 483.018 486.160 489.301 ...
    492.443 495.584 498.726 501.868 505.009 508.151 511.292 514.434 ...
    517.576 520.717 523.859 527.001 530.142 533.284 536.425 539.567 ...
    542.709 545.850 548.992 552.133 555.275 558.417 561.558 564.700 ...
    567.841 570.983 574.125 577.266 580.408 583.549 586.691 589.833 ...
    592.974 596.116 599.257 602.399 605.541 608.682 611.824 614.965 ...
    618.107 621.249 624.390 627.532];

alpha = alphaR/R;

S = (2*D*alpha.^2*litDel - 2 + 2*exp(-D*alpha.^2*litDel) + ...
    2*exp(-D*alpha.^2*bigDel) - exp(-D*alpha.^2*(bigDel-litDel)) - ...
    exp(-D*alpha.^2*(bigDel+litDel)))./(alpha.^6.*(R.^2*alpha.^2-1));
sumS = sum(S);

sig = exp(-2*gamma^2*Gmax.^2/D.^2*sumS);
