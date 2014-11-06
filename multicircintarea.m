## Copyright (C) 2014 Arun Iyer
## 
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## multicircintarea

## Author: Arun Iyer <aruniyer@localhost.localdomain>
## Created: 2014-11-05

function [ A ] = multicircintarea (X, Y, R)

% compute intersection for every pair and 
% add those that form part of the intersection area
intpts = [];
for i = 1:length(X)
  for j = i+1:length(X)
    c1x = X(i);
    c1y = Y(i);
    c1r = R(i);
    c2x = X(j);
    c2y = Y(j);
    c2r = R(j);

    sol = circcirc(c1x, c1y, c1r, c2x, c2y, c2r);

    if (~isnan(sol(1,1)))
      for k = 1:2
        flag = isPartOfAll(sol(k, 1), sol(k, 2), X, Y, R);
        if (flag)
          intpts = [intpts; sol(k, 1), sol(k, 2) i j];
        end
      end
    end
  end 
end

% lets arrange these intersection pts in clockwise/counter-clockwise fashion
refx = intpts(1,1); refy = intpts(1,2);
ref = repmat([refx, refy, 0, 0], size(intpts, 1), 1);
shiftint = intpts - ref;
[val, ind] = sort(atan2(shiftint(:, 2), shiftint(:, 1)));
intpts = intpts(ind, :);

% area of the polygon
A1 = polyarea(intpts(:, 1), intpts(:, 2));

% area of the circular segments
len = size(intpts, 1);
A2 = 0;
for i = 1:len
  i1 = intpts(i, 3); j1 = intpts(i, 4);
  i2 = intpts(mod(i, len) + 1, 3); j2 = intpts(mod(i, len) + 1, 4);

  if (i1 == i2)
    cc = i1;
  else
    cc = j1;
  end
  
  pt1 = [intpts(i, 1), intpts(i, 2)];
  pt2 = [intpts(mod(i, len) + 1, 1), intpts(mod(i, len) + 1, 2)];
  
  cl = norm(pt1 - pt2, 2);
  r = R(cc);
  centang = 2*asin(cl / (2*r));
  area = (r^2 / 2)*(centang - sin(centang));
  A2 = A2 + area;
end

A = A1 + A2;

end

function [ sol ] = circcirc (x1, y1, r1, x2, y2, r2)

% shift matrix
SM = [1 0 -x1; 0 1 -y1; 0 0 1];
SMinv = [1 0 x1; 0 1 y1; 0 0 1];

% rotation matrix
theta = -atan2(y2 - y1, x2 - x1);
RM = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
RMinv = RM';

% shift and rotate
tt = RM*SM*[x1; y1; 1];
x1 = tt(1); y1 = tt(2);
tt = RM*SM*[x2; y2; 1];
x2 = tt(1); y2 = tt(2);

% solve for intersection
d = x2;
r = r2;
R = r1;

num = (d^2 - r^2 + R^2);
xout1 = num / (2*d);
yout1 = sqrt((4*d^2*R^2 - num^2)/(4*d^2));
xout2 = xout1;
yout2 = -yout1;

% rotate and shift
tt = SMinv*RMinv*[xout1; yout1; 1];
xout1 = tt(1); yout1 = tt(2);
tt = SMinv*RMinv*[xout2; yout2; 1];
xout2 = tt(1); yout2 = tt(2);

% return the solution
sol = [xout1, yout1; xout2, yout2];

end

function [ flag ] = isPartOfAll(sx, sy, X, Y, R)

flag = 1;
for i = 1:length(X)
  if (~inCircle(sx, sy, X(i), Y(i), R(i)))
    flag = 0;
  end
end

end

function [ flag ] = inCircle(sx, sy, x, y, r)

d = (sx - x)^2 + (sy - y)^2;
flag = 1;
if (d - r^2 > eps)
  flag = 0;
end

end
