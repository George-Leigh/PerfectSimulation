# -*- Text -*-
##############################################################################
# R script to create figures for paper on perfect simulation from unbiased
#   simulation

# Author: George Leigh, 201908, updated, 202001, 202005, 202105, 202107,
#   202206, 202301, 202305, 202307
# Run under R 4.2.2, Linux version
# Copyright 2023 George Leigh
# Associated journal article: Leigh, G. M., Yang, W-H, Wickens, M. and
#   Northrop, A. R. (2023).  "Perfect simulation from unbiased simulation".
# Citation of the journal article in any publication that uses the algorithms
#   published in it would be greatly appreciated by the authors.

# The programs in this project are free software: you can redistribute them
#   and/or modify them under the terms of the GNU General Public License as
#   published by the Free Software Foundation, either version 2 of the
#   License, or (at your option) any later version.
# The programs are distributed in the hope that they will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
#   or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#   for more details.
# You should have received a copy of the GNU General Public License along with
#   this program.  If not, see <https://www.gnu.org/licenses/>.

######################################## Figure showing elimination of holes
#   when points are grouped into intervals: this figure is not currently used
#   in the paper.
# We'll use a N(0, 1) distribution.
setEPS(horizontal = FALSE, onefile = FALSE, paper = "special",
 onefile = FALSE, pointsize = 14)
postscript(width = 6.25, height = 6.25, file = "figgrouping.ps")

plot(c(0, 5), c(0.5, 5.5), asp = 1, type = "n")

lines(c(0, 2, 2, 0, 0), c(3, 3, 5, 5, 3))
lines(c(0, 2, 2, 0, 0), c(1, 1, 3, 3, 1))

lines(c(3, 5, 5, 3, 3), c(3, 3, 5, 5, 3))
lines(c(3, 5, 5, 3, 3), c(1, 1, 3, 3, 1))

headsize = 0.15 # Size of arrow heads
headangle = 20 # Angle of arrow heads
arrows(2.1, 3, 2.9, 3, length = headsize, angle = headangle, lwd = 2)

x1 = c(0.1, 0.5, 0.75, 1.4, 1.6)
y1 = c(4.1, 4.75, 3.5, 3.15, 4.25)
x2 = c(0.25, 0.5, 0.75, 0.9, 1.25, 1.5)
y2 = c(2.2, 1.25, 2.6, 1.8, 2.5, 1.4)
xoffs = -0.05
yoffs = -0.02
points(x1, y1, pch = 16, cex = 0.8)
points(x2, y2, pch = 16, cex = 0.8)
text(x1 + xoffs, y1 + yoffs, labels = as.character(c(1, -1, 1, 1, 1)), pos = 4)
text(x2 + xoffs, y2 + yoffs, labels = as.character(c(1, -1, 1, -1, 1, 1)),
 pos= 4)
x3 = c(4, 4)
y3 = c(4, 2)
points(x3, y3, pch = 16, cex = 1)
text(x3 + xoffs, y3 + yoffs, labels = as.character(c(3, 2)), pos= 4)

dev.off()

######################################## Figure showing indicator functions
#   that, as they become finer, establish that unbiased simulation produces
#   perfect simulation.
# This figure is currently not used.
# We'll use a N(0, 1) distribution.
setEPS(horizontal = FALSE, onefile = FALSE, paper = "special",
 onefile = FALSE, pointsize = 11)
postscript(width = 6.25, height = 6.25, file = "figindicator.ps")

x = seq(-3, 3, 0.001)
y = dnorm(x)
h = 0.1 # Interval width
umin = -0.5
umax = 1.5
u1 = seq(umin, umax - h, h)
u2 = u1 + h
v1 = (pnorm(u2) - pnorm(u1)) / h
v2 = v1[c(2:length(v1), length(v1))]

l = x >= umin & x <= umax
ymax = 1.04 * max(y[l])
plot(x[l], y[l], type = "l", yaxs = "i", ylim = c(0, ymax),
 xaxs = "i", xlim = c(umin, umax), xlab = "x", ylab = "Probability density",
 col = "black", lty = 3, lwd = 2)
# Draw histogram manually.
#segments(u1, 0, u1, v, lty = 2)
#segments(u1, v, u2, v, lty = 2)
#segments(u2, 0, u2, v, lty = 2)
#segments(u1, 0, u1, v, col = "grey80")
#segments(u1, v, u2, v, col = "grey80")
#segments(u2, 0, u2, v, col = "grey80")
segments(u1, 0, u1, v1, col = "black")
segments(u1, v1, u2, v1, col = "black")
segments(u2, v1, u2, v2, col = "black")

legend(0.62, ymax, lty = c(3, 1), lwd = c(2, 1), legend = c("Actual p.d.f.",
 "P.d.f. from indicator sets"))

dev.off()

######################################## Matrix showing algorithm design for
#   use of multiple perfect simulation points in each independent sample set
setEPS(horizontal = FALSE, onefile = FALSE, paper = "special", onefile =
 FALSE, pointsize = 8)
postscript(width = 5, height = 5, file = "figsampleset.ps")

#X11(width = 5, height = 2.5)
plot(c(-5, 95), c(5, 105), asp = 1, type = "n")
#plot(c(0, 100), c(66.4, 103), asp = 1, type = "n",
# xaxt = "n", yaxt = "n", ann = FALSE, axes = FALSE)
Slice = 0.8
delta = 0.001
theta1 = (pi / 2) * seq(0, Slice, delta)
theta2 = (pi / 2) * seq(2 - Slice, 2, delta)
theta3 = (pi / 2) * seq(2, 2 + Slice, delta)
theta4 = (pi / 2) * seq(-Slice, 0, delta)
r = 3
x1 = r * cos(theta1)
x2 = r * cos(theta2)
x3 = r * cos(theta3)
x4 = r * cos(theta4)
y1 = r * sin(theta1)
y2 = r * sin(theta2)
y3 = r * sin(theta3)
y4 = r * sin(theta4)

# Left bracket
x = r + 1.5 + c(x2, x3)
ytop = 95.5
ybottom = 69.3
y = c(ytop + y2, ybottom + y3)
lines(x, y)

for (i in 1:6)
 for (j in 1:6) {
  xj = 5 * j
  yi = 100 - 5 * i
  if (i == 4 | j == 4) {
v   text(xj, yi, expression(cdots))
  } else if (i == j) {
   text(xj, yi, "I")
  } else if ((i - j) %% 6 == 1) {
   text(xj, yi, "F")
  } else {
   arrows(xj - 2, yi, xj + 2, yi, length = 0.05)
  }
 }

# Right bracket
x = 30.1 + c(x4, x1)
y = c(ybottom + y4, ytop + y1)
lines(x, y)

# Arrows indicating origin and destination points
xoffs = -2
arrows(xoffs, ytop + r, xoffs, ybottom - r, length = 0.06)
text(xoffs - 1.5, 0.5 * (ytop + ybottom), "Chain number", adj = c(0.5, 1),
 srt = -90)
yoffs = ytop + r + 4
xoffs1 = 1.5
xoffs2 = 30.1 + r
arrows(xoffs1, yoffs, xoffs2, yoffs, length = 0.06)
text(0.5 * (xoffs1 + xoffs2), yoffs + 2, "Random number block", adj = c(0.5, 0))

dev.off()

######################################## Figure showing maximal coupling using
#   solid disc uniform distributions in 2-D
setEPS(horizontal = FALSE, onefile = FALSE, paper = "special",
 onefile = FALSE, pointsize = 14)
postscript(width = 8.5, height = 8.5, file = "figmaximal.ps")

plot(c(0.75, 4.25), c(0.75, 4.25), asp = 1, xaxs = "i", yaxs = "i", type = "n")

# Plot centre points (origins of Metropolis-Hastings jumps)
xc = c(2, 3)
yc = c(3, 2)
points(xc, yc, pch = 1)

# Label the centre points for processes X and Y.
xoffs = rep(0.015, length(xc))
#yoffs = rep(-0.04, length(xc))
yoffs = c(-0.04, -0.08)
text(xc + xoffs, yc + yoffs, labels = c("X", "Y"), pos = 2, cex = 1.3)

# Draw circles for the discs within which the jumps must lie.
r = 1
theta = 2 * pi * (0:1000) / 1000
x = list()
y = list()
for (i in 1:length(xc)) {
 x[[i]] = xc[i] + r * cos(theta)
 y[[i]] = yc[i] + r * sin(theta)
 lines(x[[i]], y[[i]])
}

# Shade the area of overlap.
thetaint = pi * c(1.5, 2, 0.5, 1)
l1 = theta >= thetaint[1] & theta <= thetaint[2]
l2 = theta >= thetaint[3] & theta <= thetaint[4]
xbound = c(x[[1]][l1], x[[2]][l2])
ybound = c(y[[1]][l1], y[[2]][l2])
polygon(xbound, ybound, col = "grey85")

# Draw the line bisecting the vector XY.
xbisect = c(1.5, 3.5)
ybisect = c(1.5, 3.5)
lines(xbisect, ybisect, lty = 2, lwd = 2)

# Draw a coalesced proposal point, with arrows to it.
xcoal = 2.19
ycoal = 2.40
points(xcoal, ycoal, pch = 16, cex = 1.3)

#xoffscoal1 = c(0.008, -0.035)
#yoffscoal1 = c(-0.04, 0.023)
#xoffscoal2 = c(-0.02, 0.035)
#yoffscoal2 = c(0.04, -0.026)
xoffscoal1 = c(0.008, -0.035)
yoffscoal1 = c(-0.04, 0.02)
xoffscoal2 = c(-0.02, 0.035)
yoffscoal2 = c(0.04, -0.02)
arrows(xc + xoffscoal1, yc + yoffscoal1, xcoal + xoffscoal2,
 ycoal + yoffscoal2, lwd = 2, length = 0.18, angle = 20)

# Set the coordinates of the non-coalesced proposal point for X, which we'll
#   call Xstar.
xXstar = 1.92
yXstar = 3.74

# Coordinates of point C on figure
xC = 2.83
yC = 2.83

# Working for coordinates of other points A, B, D and E: z = x-offset to get
#   from point C to point D, y-offset = negative of x-offset:
# (0.83 + z)^2 + (-0.17 - z)^2 = 1
# 2 * z^2 + 2 * z + 0.83^2 + 0.17^2 = 1
# 2 * z^2 + 2 * z + c = 0, c = 0.83^2 + 0.17^2 - 1
# z = {-1 +/- sqrt(1 - 2 * c)} / 2

c = 0.83^2 + 0.17^2 - 1
sDelta = sqrt(1 - 2 * c)
zplus = (-1 + sDelta) / 2
zminus = (-1 - sDelta) / 2

xA = xC + zminus
yA = yC - zminus
xB = xC - zplus
yB = yC + zplus
xD = xC + zplus
yD = yC - zplus
xE = xC - zminus
yE = yC + zminus

segments(c(xA, xD), c(yA, yD), c(xB, xE), c(yB, yE))

xAtoE = c(xA, xB, xC, xD, xE)
yAtoE = c(yA, yB, yC, yD, yE)
points(xAtoE, yAtoE, pch = 16, cex = 1.3)

# Label points A to E.
#xoffs = c(0.165, 0.165, 0.015, 0.015, 0.015)
#yoffs = c(0.085, 0.085, -0.01, -0.04, -0.04)
xoffs = c(0.07, 0.165, 0.015, 0.015, 0.015)
yoffs = c(-0.13, 0.085, -0.01, -0.04, -0.04)
text(xAtoE + xoffs, yAtoE + yoffs, labels = LETTERS[1:5], pos = 2, cex = 1.3)

# Now, by our system, the vector from A to Xstar is the same as from D to
#   Ystar, where Ystar is the noncoalsesced proposal from a start at point Y
#   on the figure.
xYstar = xD + (xXstar - xA)
yYstar = yD + (yXstar - yA)

# Draw the noncoalesced points, with arrows.
xnoncoal = c(xXstar, xYstar)
ynoncoal = c(yXstar, yYstar)

xoffs = rep(0.265, 3)
yoffs = rep(0.025, 3) + c(0.005, 0, 0)
text(c(xcoal, xXstar, xYstar) + xoffs, c(ycoal, yXstar, yYstar) + yoffs,
 labels = c("X*", "X*", "Y*"), pos = 2, cex = 1.3)

points(xnoncoal, ynoncoal, pch = 16, cex = 1.3)
#xoffsnoncoal1 = c(0.007, 0.049)
#yoffsnoncoal1 = c(0.049, 0.005)
#xoffsnoncoal2 = c(-0.007, -0.049)
#yoffsnoncoal2 = c(-0.049, -0.007)
xoffsnoncoal1 = c(-0.002, 0.01)
yoffsnoncoal1 = c(0.049, 0.04)
xoffsnoncoal2 = c(0.002, -0.01)
yoffsnoncoal2 = c(-0.049, -0.02)
arrows(xc + xoffsnoncoal1, yc + yoffsnoncoal1, xnoncoal + xoffsnoncoal2,
 ynoncoal + yoffsnoncoal2, lwd = 2, length = 0.18, angle = 20, lty = 3)

dev.off()

######################################## Figure showing transition matrix for
#   the two-state example
# I edited the Postscript output file to make the Greek letter theta italic
#   (or, more correctly, "oblique").  I don't know how to do it any other way.
setEPS(horizontal = FALSE, onefile = FALSE, paper = "special",
 onefile = FALSE, pointsize = 8)
postscript(width = 5, height = 5, file = "figtransitionsimple.ps")

#X11(width = 5, height = 2.5)
plot(c(-5, 95), c(5, 105), asp = 1, type = "n")
#plot(c(0, 100), c(66.4, 103), asp = 1, type = "n",
# xaxt = "n", yaxt = "n", ann = FALSE, axes = FALSE)
Slice = 0.8
delta = 0.001
theta1 = (pi / 2) * seq(0, Slice, delta)
theta2 = (pi / 2) * seq(2 - Slice, 2, delta)
theta3 = (pi / 2) * seq(2, 2 + Slice, delta)
theta4 = (pi / 2) * seq(-Slice, 0, delta)
r = 3
x1 = r * cos(theta1)
x2 = r * cos(theta2)
x3 = r * cos(theta3)
x4 = r * cos(theta4)
y1 = r * sin(theta1)
y2 = r * sin(theta2)
y3 = r * sin(theta3)
y4 = r * sin(theta4)

# Draw matrix
x = r + 1.5 + c(x2, x3)
ytop = 95.5
ybottom = 89.3
y = c(ytop + y2, ybottom + y3)
lines(x, y)

text(6, 91.5, expression(1 - theta * italic(p)), pos = 3)
#text(13.5, 91.5, expression(theta * italic(p)), pos = 3)
text(17, 91.5, expression(theta * italic(p)), pos = 3)
text(6, 86.5, expression(italic(p)), pos = 3)
#text(13.5, 86.5, expression(1 - italic(p)), pos = 3)
text(17, 86.5, expression(1 - italic(p)), pos = 3)

#for (i in seq(95, 90, -5))
# for (j in seq(5, 10, 5))
#  text(j, i, "1")

#x = 15.4 + c(x4, x1)
x = 18.9 + c(x4, x1)
y = c(ybottom + y4, ytop + y1)
lines(x, y)

dev.off()

######################################## Line plot of proportion of holes in
#   the two-state example.
load("Example_2state.RData")
#x11(width = 10, height = 6.25)
setEPS(horizontal = FALSE, onefile = FALSE, paper = "special",
 onefile = FALSE, pointsize = 12.95)
postscript(width = 10, height = 6.25, file = "fignumberofholes.ps")
par(cex = 1.4)
plot(as.numeric(rownames(Summ)), Summ[, "Total"], type = "l",
 xlab = "k", ylab = "Mean number of holes per simulation",
 xaxs = "i", yaxs = "i", ylim = c(0, 2.55), lwd = 2)
dev.off()
