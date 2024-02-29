function [Enc_Mat,Modes] = CalcEncMat(Enc_Scheme)

% Calculate Encoding Matrix
        if strcmp(Enc_Scheme,'FE')
            Modes = input('How many Fourier Encoding Modes? (8-16) ');
        elseif strcmp(Enc_Scheme,'B1TIAMO')
            Modes = 10; % 8 individual relative maps + 2 prepared shims (CP+ and CP2+)
        else
                    Modes = 8; % Modes
        end
        Channels = 8; % Channels
        Enc_Mat = zeros(Modes,Channels); % Pre-allocate Encoding Matrix
        for mode = 1:Modes
            for channel = 1:Channels
                
                if strcmp(Enc_Scheme,'FE')
                    Enc_Mat(mode,channel) = exp((2*pi*1i*(mode-1)*channel)/Modes);
                    
                elseif strcmp(Enc_Scheme,'B1TIAMO')
                    if mode <= 8 
                    if channel == mode
                        % Individual channels
                        Enc_Mat(mode,channel) =  1;
                    end    
                    elseif mode == 9
                        % CP+ 
                        Enc_Mat(mode,channel) = exp((2*pi*1i*(2-1)*channel)/8);
                    elseif mode == 10
                        % CP2+
                        Enc_Mat(mode,channel) = exp((2*pi*1i*(3-1)*channel)/8);
                    end
                end
            end
        end
end

