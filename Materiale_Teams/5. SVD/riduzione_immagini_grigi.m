% Compressione immagine tramite SVD  

% 1. Lettura dell’immagine e conversione in scala di grigi
I_orig = imread('landscape.jpg');
I_gray = rgb2gray(I_orig);
I = double(I_gray);   % converti in double per i calcoli

% 2. Visualizza immagine originale
figure;
imshow(uint8(I_gray));
title('Immagine originale in scala di grigi');

% 3. Calcola la SVD
[U, S, V] = svd(I, 'econ');   % 'econ' per dimensione ridotta

% 4. Mostra i valori singolari (opzionale)
singvals = diag(S);
figure;
plot(singvals, 'o-');
title('Valori singolari');
xlabel('i');
ylabel('\sigma_i');

% 5. Scegli il numero di componenti k da mantenere
k = 50;   % per esempio, puoi variare: 10, 30, 100 …
fprintf('Usando k = %d componenti\n', k);

% 6. Ricostruisci l’immagine usando solo le prime k componenti
Uk = U(:,1:k);
Sk = S(1:k,1:k);
Vk = V(:,1:k);
I_approx = Uk * Sk * Vk';

% 7. Visualizza l’immagine ricostruita
figure;
imshow(uint8(I_approx));
title(sprintf('Ricostruzione con k = %d componenti', k));

% 8. Calcola l’errore di ricostruzione (norma Frobenius relativa)
err = norm(I - I_approx, 'fro') / norm(I, 'fro');
fprintf('Errore relativo (Frobenius) = %.4f\n', err);

% 9. Visualizza affiancate originale e ricostruita
figure;
subplot(1,2,1), imshow(uint8(I));
title('Originale');
subplot(1,2,2), imshow(uint8(I_approx));
title(sprintf('Ricostruita (k = %d)', k));
