function html_color = wl2htmlcolor(wavelength, gamma)
    % Costanti per la conversione delle lunghezze d'onda in RGB
    intensity_max = 255;
    factor = intensity_max - 1;

    % Conversione della lunghezza d'onda in RGB
    if wavelength >= 380 && wavelength < 440
        R = -(wavelength - 440) / (440 - 380);
        G = 0.0;
        B = 1.0;
    elseif wavelength >= 440 && wavelength < 490
        R = 0.0;
        G = (wavelength - 440) / (490 - 440);
        B = 1.0;
    elseif wavelength >= 490 && wavelength < 510
        R = 0.0;
        G = 1.0;
        B = -(wavelength - 510) / (510 - 490);
    elseif wavelength >= 510 && wavelength < 580
        R = (wavelength - 510) / (580 - 510);
        G = 1.0;
        B = 0.0;
    elseif wavelength >= 580 && wavelength < 645
        R = 1.0;
        G = -(wavelength - 645) / (645 - 580);
        B = 0.0;
    elseif wavelength >= 645 && wavelength <= 780
        R = 1.0;
        G = 0.0;
        B = 0.0;
    else
        R = 0.0;
        G = 0.0;
        B = 0.0;
    end

    % Correzione gamma
    R = (R <= 0.0031308) * 12.92 * R .^ gamma + (R > 0.0031308) * (1.055 * R .^ (1.0 / gamma) - 0.055);
    G = (G <= 0.0031308) * 12.92 * G .^ gamma + (G > 0.0031308) * (1.055 * G .^ (1.0 / gamma) - 0.055);
    B = (B <= 0.0031308) * 12.92 * B .^ gamma + (B > 0.0031308) * (1.055 * B .^ (1.0 / gamma) - 0.055);

    % Normalizzazione dei valori RGB
    R = max(0, min(1, R));
    G = max(0, min(1, G));
    B = max(0, min(1, B));

    % Conversione RGB in valori esadecimali per HTML
    html_color = sprintf('#%02X%02X%02X', round(R * factor), round(G * factor), round(B * factor));
end

