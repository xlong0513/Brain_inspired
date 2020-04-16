function y = rectify(x)

% rectify is a linear threshold unit with threshold 0.

y = (x>0);   
y = y.*x;
