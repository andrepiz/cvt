function rgb = wl2rgb(wl_nm, gamma)

    % Conversione della lunghezza d'onda in RGB
    if wl_nm >= 380 && wl_nm < 440
        R = -(wl_nm - 440) / (440 - 380);
        G = 0.0;
        B = 1.0;
    elseif wl_nm >= 440 && wl_nm < 490
        R = 0.0;
        G = (wl_nm - 440) / (490 - 440);
        B = 1.0;
    elseif wl_nm >= 490 && wl_nm < 510
        R = 0.0;
        G = 1.0;
        B = -(wl_nm - 510) / (510 - 490);
    elseif wl_nm >= 510 && wl_nm < 580
        R = (wl_nm - 510) / (580 - 510);
        G = 1.0;
        B = 0.0;
    elseif wl_nm >= 580 && wl_nm < 645
        R = 1.0;
        G = -(wl_nm - 645) / (645 - 580);
        B = 0.0;
    elseif wl_nm >= 645 && wl_nm <= 780
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

    rgb =  [R, G, B];
end

