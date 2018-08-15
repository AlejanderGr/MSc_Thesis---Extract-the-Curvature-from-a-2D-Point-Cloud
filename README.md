# MSc_Thesis---Extract-the-Curvature-from-a-2D-Point-Cloud
The first step is to extract the feature that we are going to use. This was chosen to be the curvature due to it’s very important characteristic of the independence from the Coordinate System. We present 4 techniques . The best one is this on the file calcCurvParam2nd and this is the one that we used on our final experiments.
cacalcCurv2nd: Estimate the function y(x) using 2nd order Polynomials in every step
calcCurv6th: Estimate the function y(x) using 5th order Polynomials. We devide the cloud into proper regions and do a new estimate for every region 
calcCurvParam5th: We extract the parametrical functions x(s) y(s) where s it the length from the starting point. Afterwards we estimate this two functions using one 5th order polynomial for each of them.
calcCurvParam2nd:  We extract the parametrical functions x(s) y(s) where s it the length from the starting point. Afterwards we estimate this two functions using 2nd order polynomials and doing new estimation in every step.
