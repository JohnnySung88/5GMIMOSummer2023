function [LLR] = LLRcul(y,No,a1,a2,b1,b2)
    LLR = (1/(2*No)) .* ( min( (y-[ a1; a2]).^2 ) - min( (y-[b1;b2]).^2 ));
end