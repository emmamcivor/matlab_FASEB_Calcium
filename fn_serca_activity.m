function serca_activity = fn_serca_activity(cj,cs,x_serca_i,y_serca_i,Kmf,Kmr_Bers1998,n_H)
% This function calculates the activity of a SERCA pump at a specified
% position on the ER membrane.
% returns a number 
% returns a number whose absolute value is between 0 and 1 at baseline
% cytoplasmic calcium concentration
% returns a number whose absolute value is between 0 and 1 at high
% cytoplasmic calcium concentration
% returns a number that is of a sensible size at baseline cytoplasmic calcium
% concentration
% returns a number that is greater at high cytoplasmic calcium concentrations than
% baseline cytoplasmic calcium concentrations. Curve is a saturating function at high
% cytoplasmic calcium, as expected.
% inputs: 
% cj: a matrix of size (|x| x |y| x |z|)
% cs: a matrix of size (|x| x |y| x |z|)
% x_serca_i: a number representing the x co-ordinate on the ER membrane
% y_serca_i: a number representing the x co-ordinate on the ER membrane
% Kmf: a number representing the cytoplasmic calcium concentration resulting in half
% maximal SERCA pump activation (forward pumping)
% Kmr_Bers1998: a number representing the sub-PM ER calcium concentration resulting in half
% maximal SERCA pump activation (reverse pumping)
% n_H: a number representing the Hill coefficient used in the Hill equation


J_SERCA_num=(cj(x_serca_i,y_serca_i,1)/Kmf)^n_H - (cs(x_serca_i,y_serca_i,end)/Kmr_Bers1998)^n_H;
J_SERCA_den=1+(cj(x_serca_i,y_serca_i,1)/Kmf)^n_H +(cs(x_serca_i,y_serca_i,end)/Kmr_Bers1998)^n_H;
serca_activity=J_SERCA_num/J_SERCA_den;

end