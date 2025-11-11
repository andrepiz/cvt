function html_color = wl2htmlcolor(wavelength, gamma)
    % Costanti per la conversione delle lunghezze d'onda in RGB
    intensity_max = 255;
    factor = intensity_max - 1;

    rgb = wl2rgb(wavelength, gamma);

    % Conversione RGB in valori esadecimali per HTML
    html_color = sprintf('#%02X%02X%02X', round(rgb(1) * factor), round(rgb(2) * factor), round(rgb(3) * factor));
end

