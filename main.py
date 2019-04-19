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
import docx
from docx.shared import Cm
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_TAB_ALIGNMENT
import io
from docx.shared import Pt, RGBColor
import os
from math import sin, cos, pi, atan2, sqrt

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

# src: https://stackoverflow.com/a/39885430
def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    try:
        # PyInstaller creates a temp folder and stores path in _MEIPASS
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.environ.get("_MEIPASS2",os.path.abspath("."))
    return os.path.join(base_path, relative_path)

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
        doc = docx.Document()
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
                W = 2220 # in px
                Wc = 10000 # in chart units
                M = 36
                x = lambda x: x*W/Wc+W/2
                y = lambda x: x*W/Wc
                img = Image.new('RGBA', (W+M, W+M), (255, 255, 255, 0))
                imgd = ImageDraw.Draw(img)
                circ = lambda x, y, r, a={'outline': 'red'}, i=None: imgd.ellipse((x-r, y-r, x+r, y+r), **a) if i is None else imgd.arc((x-r, y-r, x+r, y+r), i[0], i[1], **a)
                cx, cy, cr = x(circle.a), x(-circle.b)+M, y(circle.r)
                circ(cx, cy, 8, {'fill': 'red'})
                circ(cx, cy, cr)
                imgd.line((0, W/2+M, W, W/2+M), fill='black')
                imgd.line((W-15, W/2-5+M, W, W/2+M, W-15, W/2+5+M), fill='black')
                imgd.line((W/2, M, W/2, W+M), fill='black')
                imgd.line((W/2-5, 15+M, W/2, M, W/2+5, 15+M), fill='black')
                fnt = ImageFont.truetype(resource_path('./fonts/Roboto-Regular.ttf'), 24)
                fntb = ImageFont.truetype(resource_path('./fonts/Roboto-Regular.ttf'), 48)
                imgd.text((W/2, 0), 'Y', font=fntb, fill='black')
                imgd.text((W, W/2), 'X', font=fntb, fill='black')
                scale = y(1000)
                for i in range(1, 10):
                  l = 12
                  w = 0
                  if i == 5:
                    l = 20
                    w = 2
                  imgd.line([16+i*scale/10, W-l+M, 16+i*scale/10, W-2+M], fill='black', width=w)
                imgd.line([16, W-18+M, 16, W-2+M, 16+scale, W-2+M, 16+scale, W-18+M], fill='black')
                imgd.text((16-6, W-18-42+M), '0', font=fnt, fill='black')
                imgd.text((16+scale-50, W-18-42+M), '1000 mm', font=fnt, fill='black')
                imgpx = img.load()
                arr = [(x(data1.X[i]), x(-data1.Y[i])+M, str(i+1), False) for i in range(data1.n)]
                for e in arr:
                    circ(e[0], e[1], 5, {'outline': 'black'})
                    circ(e[0], e[1], 4, {'outline': 'black'})
                    circ(e[0], e[1], 3, {'outline': 'black'})
                arr += [(cx, cy, 'S', True)]
                for e in arr:
                    r = 33 # radius from point
                    ex, ey = e[0], e[1]
                    angle = atan2(ex-cx, ey-cr)
                    if cr > sqrt((ex-cx)**2 + (ey-cy)**2):
                        angle += pi
                    for r in (22,33,44,55,66,77,88,99):
                        s = 54 if e[3] else 30 # size of font + 6
                        for a in range(24):
                            a = a/12*pi + angle
                            m = len(e[2])/2
                            ix = int(ex+r*sin(a)-m*s/2)
                            iy = int(ey+r*cos(a)-s/2)
                            good = True
                            for dx in range(int(m*s)):
                                for dy in range(s):
                                    good = good and imgpx[ix+dx, iy+dy] == (255, 255, 255, 0)
                                    if not good:
                                        break
                                if not good:
                                    break
                            if not good:
                                continue
                            imgd.text((ix, iy), e[2], font=(fntb if e[3] else fnt), fill=('red' if e[3] else 'black'))
                            break
                        if good:
                           break
                    if not good:
                        print('ERR: no place for "',e[2],'"',sep='')
                p = doc.add_paragraph('')
                p.paragraph_format.tab_stops.add_tab_stop(Cm(14), WD_TAB_ALIGNMENT.RIGHT)
                r = p.add_run('Przekrój: \t')
                r.bold = True
                r.font.size = Pt(14)
                r = p.add_run()
                r.alignment = WD_ALIGN_PARAGRAPH.RIGHT
                r.add_picture('img/arrow.png', height=Cm(3))
                p.add_run('\nGłębokość:  m')
                with io.BytesIO() as out:
                    img.save(out, format='PNG')
                    p = doc.add_paragraph()
                    p.add_run().add_picture(out, width=Cm(17))
                    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
                p = doc.add_paragraph()
                p.paragraph_format.tab_stops.add_tab_stop(Cm(1.5), WD_TAB_ALIGNMENT.CENTER)
                p = doc.add_paragraph('\t')
                r = p.add_run('\N{Bullet}')
                r.font.color.rgb = RGBColor(0xff, 0x00, 0x00)
                r.bold = True
                p.add_run(' - środek okręgu\n\n\tWspółrzędne środka okręgu:\n\tX')
                p.paragraph_format.tab_stops.add_tab_stop(Cm(9), WD_TAB_ALIGNMENT.LEFT)
                p.add_run('S').font.subscript = True
                #p.add_run(' = ' + str(round(circle.a)) + ' mm Y')
                p.add_run('S').font.subscript = True
                p.add_run(' = ' + str(round(circle.b)) + ' mm\n\n\tŚrednica okręgu:\n\tD = ' + str(round(circle.r*2)) + ' mm')
                doc.add_page_break()
                for section in doc.sections:
                    section.top_margin = Cm(0.25)
                    section.bottom_margin = Cm(0.25)
                    section.right_margin = Cm(0.25)
            else:
                log += 'Unexpected code:' + str(code) + '\n'
        tk = Tk()
        tk.title('Kołodziej')
        menubar = Menu(tk)
        filemenu = Menu(menubar, tearoff=0)
        saveas = lambda: (lambda x: doc.save(x) if x else '')((lambda: filedialog.asksaveasfilename(defaultextension='.docx', filetypes=(('docx files','*.docx'),('all files','*.*'))))())
        filemenu.add_command(label="Save as...", command=saveas)
        filemenu.add_command(label="Exit", command=tk.destroy)
        menubar.add_cascade(label="File", menu=filemenu)
        helpmenu = Menu(menubar, tearoff=0)
        helpmenu.add_command(label="Get source code", command=lambda:open_new('https://github.com/Nircek/kolodziej'))
        helpmenu.add_command(label="About...", command=lambda:messagebox.showinfo("Kołodziej", "Kołodziej by Nircek\nCopyright \N{COPYRIGHT SIGN} Nircek 2019"))
        menubar.add_cascade(label="Help", menu=helpmenu)
        tk.config(menu=menubar)
        tk.bind('<Control-s>', lambda x:saveas())
        w = Text(tk, width=sum(t), foreground='black')
        w.insert(INSERT, log)
        w.pack()
        w.configure(state="disabled")
        w.configure(background='white')
        w.configure(inactiveselectbackground=w.cget("selectbackground"))
        mainloop()
    except Exception as e:
        messagebox.showerror("Fatal error", traceback.format_exc())
