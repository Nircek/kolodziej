#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

# inspirated by http://people.cas.uab.edu/~mosya/cl/CPPcircle.html
# copyright (c) 2011-2014 Nikolai Chernov
# edited by Nircek 2019

# MIT License

# Copyright (c) 2019 Nircek

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

from math import sqrt

class Data:
    def __init__(self):
        self.n = 0
        self.X = []
        self.Y = []
        for i in range(len(self.n)):
            self.X += [0.]
            self.Y += [0.]

    def __init__(self, n):
        self.n = n
        self.X = []
        self.Y = []
        for i in range(len(self.n)):
            self.X += [0.]
            self.Y += [0.]

    def __init__(self, n, xs, ys):
        if len(xs) != n or len(ys) != n:
            raise Error('not as much xs or ys as declared')
        self.n = n
        self.X = xs
        self.Y = ys

    @property
    def mX(self):
        return sum(self.X) / len(self.X)

    @property
    def mY(self):
        return sum(self.Y) / len(self.Y)

    def center(self):
        mX = self.mX
        mY = self.mY
        for e in self.X:
            e -= mX
        for e in self.Y:
            e -= mY

    def scale(self):
        sXX = 0.
        sYY = 0.
        for i in range(self.n):
            sXX += X[i]*X[i]
            sYY += Y[i]*Y[i]
        scaling = sqrt((sXX+sYY)/self.n/2)
        for i in range(self.n):
            X[i] /= scaling
            Y[i] /= scaling

    def __repr__(self):
        s = ''
        s += 'DATA(len:' + str(self.n) + ')['
        for i in range(len(self.n)):
            s += '(' + str(X[i]) + ',' + str(Y[i]) + '), '
        if self.n == 0:
            s = s[:-2]
        s += ']'
        return s

class Circle:
    def __init__(self):
        self.a=0.
        self.b=0.
        self.r=1.
        self.s=0.
        self.i=0
        self.j=0

    def __init__(self, aa, bb, rr) {
        self.a=aa
        self.b=bb
        self.r=rr
        self.s=0.
        self.i=0
        self.j=0

    def __repr__(self):
        s = 'CIRCLE('
        s += 'x:' + str(self.a)
        s += ', y:' + str(self.b)
        s += ', r:' + str(self.r)
        s += ', s:' + str(self.s)
        s += ', g:' + str(self.g)
        s += ', i:' + str(self.i)
        s += ', j:' + str(self.j)
        s += ')'

def SQR(t):
    return t*t

def Sigma(data, circle):
    sum = 0.
    for i in range(self.n):
        dx = data.X[i] - circle.a
        dy = data.Y[i] - circle.b
        sum += SQR(sqrt(dx*dx+dy*dy) - circle.r)
    return sqrt(sum/data.n)

# -----------------

def CircleFitByLevenbergMarquardtFull(data, circleIni, LambdaIni, circle):
    '''
        Algorithm:  Levenberg-Marquardt running over the full parameter space (a,b,r)

        See a detailed description in Section 4.5 of the book by Nikolai Chernov:
        "Circular and linear regression: Fitting circles and lines by least squares"
        Chapman & Hall/CRC, Monographs on Statistics and Applied Probability, volume 117, 2010.

        Nikolai Chernov,  February 2014
        see the unedited code: http://people.cas.uab.edu/~mosya/cl/CircleFitByLevenbergMarquardtFull.cpp
        port from C++ to Python by Nircek in January 2019
    '''
    IterMAX = 99
    factorUp=10.
    factorDown=0.04
    ParLimit=1.e+6
    epsilon=3.e-8
    New = circleIni
    New.s = Sigma(data,New)
    lambd = LambdaIni
    iterr = inner = code = 0
NextIteration:
    Old = New
    iterr += 1
    if iterr > IterMAX:
        code = 1
        goto enough
    Mu=Mv=Muu=Mvv=Muv=Mr=0.
    for i in range(data.n):
        dx = data.X[i] - Old.a
        dy = data.Y[i] - Old.b
        ri = sqrt(dx*dx + dy*dy)
        u = dx/ri
        v = dy/ri
        Mu += u
        Mv += v
        Muu += u*u
        Mvv += v*v
        Muv += u*v
        Mr += ri
    Mu /= data.n
    Mv /= data.n
    Muu /= data.n
    Mvv /= data.n
    Muv /= data.n
    Mr /= data.n
    F1 = Old.a + Old.r*Mu - data.meanX
    F2 = Old.b + Old.r*Mv - data.meanY
    F3 = Old.r - Mr
    Old.g = New.g = sqrt(F1*F1 + F2*F2 + F3*F3)
try_again:
    UUl = Muu + lambd
    VVl = Mvv + lambd
    Nl = 1.0 + lambd
    G11 = sqrt(UUl)
    G12 = Muv/G11
    G13 = Mu/G11
    G22 = sqrt(VVl - G12*G12)
    G23 = (Mv - G12*G13)/G22
    G33 = sqrt(Nl - G13*G13 - G23*G23)
    D1 = F1/G11
    D2 = (F2 - G12*D1)/G22
    D3 = (F3 - G13*D1 - G23*D2)/G33
    dR = D3/G33
    dY = (D2 - G23*dR)/G22
    dX = (D1 - G12*dY - G13*dR)/G11
    if (abs(dR)+abs(dX)+abs(dY))/(1.0+Old.r) < epsilon:
        goto enough
    New.a = Old.a - dX
    New.b = Old.b - dY
    if abs(New.a)>ParLimit or abs(New.b)>ParLimit:
        code = 3
        goto enough
    New.r = Old.r - dR
    if New.r <= 0.:
        lambd *= factorUp
        inner += 1
        if inner > IterMAX:
            code = 2
            goto enough
        goto try_again
    New.s = Sigma(data,New)
    if New.s < Old.s:
        lambd *= factorDown
        goto NextIteration
    else:
        inner += 1
        if inner > IterMAX:
            code = 2
            goto enough
        lambd *= factorUp
        goto try_again
enough:
    Old.i = iterr
    Old.j = inner
    circle = Old
    return code



int main(int argc, const char **argv) {
    int code;
    const size_t quantity = 4;
    reals Xs[quantity] = {0,1,-1,0};
    reals Ys[quantity] = {1,0,0,-1};
    reals LambdaIni = 0.001;

    Data data1(quantity, Xs, Ys);
    Circle circle, circleIni;
    cout.precision(7);

    code = CircleFitByLevenbergMarquardtFull(data1, circleIni, LambdaIni, circle);
    if ((code == 1)||(code==2))
        cerr << "\n Geometric circle by Levenberg-Marquardt (full) did not converge. Iterations maxed out.\n";
    else if (code == 3)
        cerr << "\n Geometric circle by Levenberg-Marquardt (full) did not converge. Fitting circle too big.\n";
    else if (code == 0) {
        cout << "X Y Radius Sigma Iterations" << endl;
        cout << circle.a << ' ' << circle.b << ' ' << circle.r << ' ' << circle.s << ' ' << circle.i << endl;
        return 0;
    }
    return 1;
}
