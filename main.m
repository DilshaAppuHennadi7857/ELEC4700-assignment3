%% ELEC 4700 Assignment 3: Monte-Carlo/Finite Difference Method
%
% Dilsha Appu-Hennadi, 101107857
% Mar. 20, 2022
% Latest rev. Mar. 7, 2022

clear all
clearvars
clearvars -GLOBAL
close all
format shorte

global regL regW  nx ny L W

regL = 100e-9; % length of area
regW = 100e-9; % width of area

nx = 50; % number of space in x direction
ny = 50; % number of spaces in y direction

L = linspace(0, regL, nx); % nx spaces from 0 to region length
W = linspace(0, regW, ny); % ny spaces from 0 to region width

q = 1.60217662e-19; %C % elementary charge
m_e = 9.10938356e-31; %kg % rest mass of an electron
m_n = 0.26*m_e; %kg % mass of an electron
kb = 1.38064852e-23; %m^2 kg s^-2 K^-1 % Boltzmann's Constant

T = 300; %K % temperature

numElec = 500;%50;
numDispElec = 5;
t_mn = 0.2e-12; %s % mean time between collisions

numTimeStep = 800; %50;
dt = 1e-16;


% V_0 = 0.1; %V
% V_0 = 0.8; %V
V_0 = 1; %V

Cols = hsv(numDispElec);

%%
%
% This assignment will take the work done in the previous two assignments
% to simulate the effect of an electric field on the movement of electrons.
% For this first part, a potential of 0.1V is applied across the x
% direction of a semiconductor. We will determine the electric field
% induced by this potential. Once we know the electric field, we will know
% what the force on the electrons will be. Then, from Newton's law, we can
% determine the acceleration on the electrons and simulate the movement of
% the electrons over time.
%

assign3part1

%%
%
% The density and temperature plots also confirm that the electrons tend
% towards the right side of the region - i.e. the side with the lower
% potential.
%

%%
%
% In this next part, we will run a similar simulation, however we will now
% add a bottle neck to the region. This bottle neck will be created by
% adding in two boxed regions where the conductivity is far lower than the
% surrounding material. This will give us some idea of how these electrons
% will behave in a Si crystal with added impurities which have a lower
% conductance.
%

assign3part2

%%
%
% From the density plot, we see electrons are being guided toward the right
% side of the region by the applied voltage. We observe this as a high
% density of electrons colliding and reflecting off of the left sides of
% the two boxes. We can also vaguely see electrons spread out as they exit
% the bottleneck on the right side of the region. As a reminder, a voltage
% $V_0$ is applied at x = 0 and 0V at x = 100nm. The applied voltage
% ,$V_0$, describes an area with excess electrons, thus it makes sense that
% the electrons in the simulation would be repelled from x = 0 towards x =
% 100nm.
%

currVsBotNeck = [0.0016 0.8*regW; 0.0015 0.6*regW; 0.0013 0.4*regW; 0.0012 0.2*regW; 0.0011 0.1*regW];

%%
%
% We can observer the effects of the bottleneck on the system by taking the
% current as the bottleneck narrows. This yields the following graph:
%

figure(10)
plot(currVsBotNeck(:,2),currVsBotNeck(:,1))
title('Current vs Bottle Neck Size')
xlabel('Size of Bottle Neck')
ylabel('Current')

%%
%
% We see, as the bottle neck size increases - i.e. the pathway between the
% two boxes widens - current increases. This intuitively makes sense as the
% two boxes represent areas of higher resistivitivy. These areas would
% impede current flow, thus decreasing the current through the device.
%
