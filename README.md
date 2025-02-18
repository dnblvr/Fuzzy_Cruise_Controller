# Fuzzy Advanced Cruise Controller

The objective of this project is to develop a model of a car's cruise controller that automatically maintains a certain distance with the front car. This braking distance, which is a safe distance that allows the driver to notice in time and perform a fast-enough brake, is proportional to the speed of the car. When given a chart showing the approximate distances at certain discrete velocities, this is a near-perfect application of fuzzy logic as values in-between can be interpolated. 

Additionally, the speed of the car is given as an input and is regulated at this speed with a multiple-in multiple-out (MIMO) fuzzy-based cruise controller. This controller takes in magnitudes of current position and velocities. Direction of changes to output velocities and accelerations are opposite to the direction of their errors.
