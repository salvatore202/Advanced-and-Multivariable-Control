% ================================================
%  COMPRESSIONE IMMAGINE A COLORI CON SVD
% ================================================

clc; clear; close all;

% Lettura dell'immagine
I_orig = imread('Landscape.jpg'); % matrice a 3 dimensioni mxn (pixel) x3 (rgb)
I_orig = im2double(I_orig);  % Normalizza valori [0,1]

figure;
imshow(I_orig);
title('Immagine originale');

% Dimensioni immagine
[m, n, ~] = size(I_orig);
fprintf('Dimensioni immagine: %d x %d pixel\n', m, n);

% Numero di componenti da mantenere 
k = round(min([m n])/2);   % prova anche 10, 30, 100
k = 30;   % prova anche 10, 30, 100
fprintf('Numero di componenti mantenute: k = %d\n', k);

% Inizializzazione immagine compressa
I_compressa = zeros(m, n, 3);

% Applica SVD separatamente a ciascun canale RGB
for c = 1:3
    A = I_orig(:,:,c);
    [U, S, V] = svd(A, 'econ');
    Uk = U(:,1:k);
    Sk = S(1:k,1:k);
    Vk = V(:,1:k);
    I_compressa(:,:,c) = Uk * Sk * Vk';
end

% Visualizza risultato
figure;
imshow(I_compressa);
title(sprintf('Immagine ricostruita con k = %d', k));

% Confronto visivo
figure;
subplot(1,2,1);
imshow(I_orig);
title('Originale');
subplot(1,2,2);
imshow(I_compressa);
title(sprintf('Compressa (k = %d)', k));

%  errore medio
errore = norm(I_orig(:) - I_compressa(:)) / norm(I_orig(:));
fprintf('Errore relativo di ricostruzione = %.4f\n', errore);

% 9️⃣ Calcola tasso di compressione
% Numero di parametri salvati rispetto all'immagine originale
num_originale = m * n * 3;
num_compressi = k * (m + n + 1) * 3;   % U, S, V per 3 canali
compressione = 100 * (1 - num_compressi / num_originale);
fprintf('Compressione ≈ %.2f%% di spazio risparmiato\n', compressione);

% 🔟 (Opzionale) Salva il risultato
imwrite(I_compressa, sprintf('landscape_compressa_k%d.jpg', k));
