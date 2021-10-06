% funcion hom que regresa las coordenadas homogeneas (Hs) de Y con escala s
function [Hs] = hom(Y,S)
    if nargin < 2
        S = 1;
    end 
    
    Hs = [Y ; ones(1,size(Y,2))*S];
end
