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
from sys import argv
from tkinter import filedialog, messagebox
from tkinter import *
import traceback
from PIL import ImageTk, Image, ImageDraw, ImageFont
from webbrowser import open_new

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
    def __init__(self, aa=0., bb=0., rr=1.):
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
    for i in range(data.n):
        dx = data.X[i] - circle.a
        dy = data.Y[i] - circle.b
        sum += SQR(sqrt(dx*dx+dy*dy) - circle.r)
    return sqrt(sum/data.n)

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
    iterr = inner = 0
    code = -1
    while 1:
        Old = New
        iterr += 1
        if iterr > IterMAX:
            code = 1
            break
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
        F1 = Old.a + Old.r*Mu - data.mX
        F2 = Old.b + Old.r*Mv - data.mY
        F3 = Old.r - Mr
        Old.g = New.g = sqrt(F1*F1 + F2*F2 + F3*F3)
        while 1:
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
            if (abs(dR)+abs(dX)+abs(dY))/(1.+Old.r) < epsilon:
                code = 0
                break
            New.a = Old.a - dX
            New.b = Old.b - dY
            if abs(New.a)>ParLimit or abs(New.b)>ParLimit:
                code = 3
                break
            New.r = Old.r - dR
            if New.r <= 0.:
                lambd *= factorUp
                inner += 1
                if inner > IterMAX:
                    code = 2
                break
            New.s = Sigma(data,New)
            if New.s < Old.s:
                lambd *= factorDown
                break
            else:
                inner += 1
                if inner > IterMAX:
                    code = 2
                else:
                    lambd *= factorUp
                break
        if code != -1:
            break
    Old.i = iterr
    Old.j = inner
    return code, Old

def show(img):
    root = Tk()
    root.title('Kołodziej')
    n = lambda:0
    menubar = Menu(root)
    filemenu = Menu(menubar, tearoff=0)
    saveas = lambda: (lambda x: img.save(x) if x else '')((lambda: filedialog.asksaveasfilename(filetypes=(('png files','*.png'),('all files','*.*'))))())
    filemenu.add_command(label="Save as...", command=saveas)
    filemenu.add_command(label="Exit", command=root.destroy)
    menubar.add_cascade(label="File", menu=filemenu)
    helpmenu = Menu(menubar, tearoff=0)
    helpmenu.add_command(label="Get source code", command=lambda:open_new('https://github.com/Nircek/kolodziej'))
    helpmenu.add_command(label="About...", command=lambda:messagebox.showinfo("Kołodziej", "Kołodziej by Nircek\nCopyright \N{COPYRIGHT SIGN} Nircek 2019"))
    menubar.add_cascade(label="Help", menu=helpmenu)
    root.config(menu=menubar)
    tkimage = ImageTk.PhotoImage(img)
    Label(root, image=tkimage).pack()
    root.bind('<Control-s>', lambda x:saveas())
    root.mainloop()

if __name__ == '__main__':
    try:
        log = ''
        files = argv[1:]
        if len(files) == 0:
            tk = Tk()
            tk.withdraw()
            files = filedialog.askopenfilenames(filetypes=(('txt files','*.txt'),('all files','*.*')))
            tk.destroy()
        if len(files) != 0:
            s = lambda i: str(i).replace('.', ',')
            b = lambda s, l: str(s)[:l] + (l-len(str(s)))*' '
            t = [
                    len(str(len(files)))+1,
                    max([len(x) for x in files])+1,
                    len(str(-0.0001/3))+1,
                    len(str(-0.0001/3))+1,
                    len(str(-0.0001/3))+1,
                    len(str(-0.0001/3))+1,
                    len('Iterations')
                ]
            log += b('I', t[0]) + b('Name', t[1]) + b('X', t[2]) + b('Y', t[3]) + b('Radius', t[4]) + b('Sigma', t[5]) + 'Iterations' + '\n'
        else:
            t=[100]
        for i in range(len(files)):
            f = open(files[i], 'r')
            Xs = []
            Ys = []
            for ff in f:
                ff = ff.replace(',', '.')
                ff = ff.split()
                if len(ff) < 2:
                    log += 'WARN: ignoring \'' + ''.join(ff) + '\'\n'
                    continue
                if len(ff) > 2:
                    log += 'WARN: ignoring \'' + ''.join(ff[2:]) +  '\'\n'
                x = float(ff[0])
                y = float(ff[1])
                Xs += [x]
                Ys += [y]
            LambdaIni = 0.001
            data1 = Data(len(Xs), Xs, Ys)
            circle, circleIni = Circle(), Circle()
            code, circle = CircleFitByLevenbergMarquardtFull(data1, circleIni, LambdaIni, circle)
            log += b(i+1, t[0]) + b(files[i], t[1])
            if code == 1 or code == 2:
                log += 'ERR: Iterations maxed out.\n'
            elif code == 3:
                log += 'ERR: Fitting circle too big.\n'
            elif code == 0:
                log += b(s(circle.a), t[2]) + b(s(circle.b), t[3]) + b(s(circle.r), t[4]) + b(s(circle.s), t[5]) + str(s(circle.i)) + '\n'
                W = 595-16 # in px
                Wc = 10000 # in chart units
                x = lambda x: x*W/Wc+W/2
                y = lambda x: x*W/Wc
                img = Image.new('RGBA', (W, W), 'white')
                imgd = ImageDraw.Draw(img)
                circ = lambda x, y, r, a={'outline': 'red'}: imgd.ellipse((x-r, y-r, x+r, y+r), **a)
                circ(x(circle.a), x(-circle.b), 3, {'fill': 'red'})
                circ(x(circle.a), x(-circle.b), y(circle.r))
                imgd.line((0, W/2, W, W/2), fill='black')
                imgd.polygon((W, W/2, W-9, W/2-9, W-9, W/2+9), outline='black')
                imgd.line((W/2, 0, W/2, W), fill='black')
                imgd.polygon((W/2, 0, W/2-9, 9, W/2+9, 9), outline='black')
                for i in range(data1.n):
                    circ(x(data1.X[i]), x(-data1.Y[i]), 2, {'fill':'black'})
                pap = Image.new('RGBA', (595, 842), 'white')
                pap.paste(img, ((pap.size[0]-W)//2, (pap.size[1]-W)//2))
                papd = ImageDraw.Draw(pap)
                fontsize = 20
                fnt = ImageFont.truetype('fonts/Roboto-Regular.ttf', fontsize)
                fnt2 = ImageFont.truetype('fonts/Roboto-Regular.ttf', fontsize//2)
                nr = '4'
                gl = '275,1'
                papd.text((16, 48), 'Przekrój nr: ' + nr, 'black', font=fnt)
                papd.text((pap.size[0]//2, 48), 'Głębokość: '+ gl + ' m', 'black', font=fnt)
                under = (pap.size[1]-W)//2+W
                papd.text((pap.size[0]//2, under),
                          'Współrzędne środka:\nXśr = ' + str(round(circle.a)) + ' mm Yśr = ' + str(round(circle.b)) + ' mm\nŚrednica przekroju:\nD = ' + str(round(circle.r*2)) + ' mm',
                          'black', font=fnt)
                scale = y(1000)
                papd.text((pap.size[0]//2-5,(pap.size[1]-W)//2-fontsize-2), 'N', 'black', font=fnt)
                papd.text((24, under-fontsize//2), '1000 m', 'black', font=fnt2)
                papd.line([16, under, 16, under+16, 16+scale, under+16, 16+scale, under], fill='black')
                show(pap)
            else:
                log += 'Unexpected code:' + str(code) + '\n'
        tk = Tk()
        tk.title('LOG')
        w = Text(tk, width=sum(t), foreground='black')
        w.insert(INSERT, log)
        w.pack()
        w.configure(state="disabled")
        w.configure(background='white')
        w.configure(inactiveselectbackground=w.cget("selectbackground"))
        mainloop()
    except Exception as e:
        messagebox.showerror("Fatal error", traceback.format_exc())
