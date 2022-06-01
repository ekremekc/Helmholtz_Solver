a = [ 1+2i,2+4i,-1+5i,2+3i,1+i]
norm_a = norm(a)
normalized = a/norm_a
normalized = normalized/abs(max(normalized))

figure(1)
subplot(2,1,1);
x = linspace(0,1,5);
y1 = real(normalized);
plot(x,y1)

subplot(2,1,2); 
y2 = imag(normalized);
plot(x,y2)

%% Seperate Normalization for real and imag parts
b = [ 1+2i,2+4i,-1+5i,2+3i,1+i];
norm_b_real = norm(real(b))
norm_b_imag = norm(imag(b))

normalized_real = real(b)/norm_b_real
normalized_real = normalized_real/abs(max(normalized_real))

normalized_imag = imag(b)/norm_b_imag
normalized_imag = normalized_imag/abs(max(normalized_imag))

figure(2)
subplot(2,1,1);
y3 = normalized_real;
plot(x,y3)

subplot(2,1,2); 
y4 = normalized_imag;
plot(x,y4)