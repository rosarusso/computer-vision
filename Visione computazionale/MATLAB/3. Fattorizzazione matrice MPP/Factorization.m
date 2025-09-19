% MPP factorization
% Rosa Russo VR445639

clear all
close all
clc

% Camera calibration:
% - obtain internal parameters (intrinsic to the camera itself),
%   which allow a mapping between camera coordinates and pixel coordinates
%   in the image frame;
% - obtaint extrinsic parameters, that define the location and orientation
%   of the camera wrt the world frame.

camera = load(['Calib_Results.mat']);

% Intrinsic parameters matrix
K = camera.KK;

% Rotation matrix
R = camera.Rc_5;

% Translation vector
t = camera.Tc_5;

P = K * [R t];

[K_fact, R_fact, t_fact] = krt(P)
