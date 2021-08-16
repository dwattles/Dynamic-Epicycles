# Dynamic-Epicycles

Program demonstrates how to draw an image in the complex plane by using the discrete Fourier transform(DFT).

First, converts an RGB image to grayscale so that the Sobel operator can be used to extract a set of points on the edges of the image.
Creates a subset of the edge points such that all points in the subset are at least DISTANCE away from eachother.
Moves points to the complex plane and centers image.
Sorts this subset using a nearest neighbour approach; this set of points will serve as the discrete samples for the DFT.
Applys a naive implementation of the DFT to the subset of edge points in order to determine approximate Fourier coefficients(c_n).
Normalizes the coefficients by the length of the subset.
Interprets the real and imaginary part of each coefficient as the magnitude and phase of a vector that rotates with frequency given by it's order in the array.
Note that coefficients with indexes greater than N/2 correspond to negative frequencies, or vectors that rotate clockwise.
Sums back the series to draw the image with the stacked vectors in an animation.
