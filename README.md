# Multi-view triangulation and Non-linear optimization

### Description

In this code, reconstructing a synthetic cube (made up of 56 points) using multi-view triangulation is implemented. Multi-view triangulation is a straight forward extension of 2-view triangulation which you have already coded in the https://github.com/Pravin854/Two-View-Spare-Reconstruction. Similar to the 2-view triangulation, projection matrices of all the 8 views and setup a least square system of the form Ax = b and then solve it using SVD is used. For example 3D point X3 must satisfy the following constraints P1*X3 = x13, P2*X3 = x23, ..., P8*X3 = x83, where x13 denotes the 2D projection of X3 in image 1, x23 denotes the 2D projection of X3 in image 2,..., x83 denotes the 2D projection of X3 in image 8.

The images of the synthetic cube is provided to you in the form of a 8 × 2 × 56 (numOfViews × 2 × 56) tensor named cube_imgs.mat. And the corresponding projection matrices are provided as 8 × 1 (numOf V iews × 1) MATLAB cell array named projMatrices.mat; Both stored in MatFilesQues1 folder.

Levenberg-Marquardt (LM) Algorithm for non-linear least square is used for non-linear least square problems.

Results are present in Report.pdf