% run_mdo_optimization.m
% Main driver for AE4205 MDO assignment – Fokker 100

clear; clc; close all;

% Global constants and reference data
const = constants_ae4205();

% Design vector initial guess and bounds (Table 1.2)
[x0, lb, ub] = design_vector_init();

