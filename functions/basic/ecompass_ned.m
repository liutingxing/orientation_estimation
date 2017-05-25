function [Cbn] = ecompass_ned(g, mag)

vz = g';
vx = mag';
vy = cross(vz, vx);
vx = cross(vy, vz);

vx = vx/norm(vx);
vy = vy/norm(vy);
vz = vz/norm(vz);

Cbn = [vx, vy, vz];
