
##############################################################################################
#User Interface
import Tkinter as tk
from Tkinter import *

class User:
    def __init__(self, master):
        self._name = "" 
        frame = Frame(master)
        frame.grid()
        master.title("AMERICAN OPTION PRICING INPUT")
        self.locator_label1 = Label(frame, text="Choose number of steps:",width=100, height=2)
        self.locator_label1.grid(row=0, sticky=E)
        self.entry1 = Entry(frame)
        self.entry1.grid(row=0, sticky=E)
        self.button1 = Button(frame, text="Confirm", command=self.context, pady=2)
        self.button1.grid(row=1)
        self.button2 = Button(frame, text="Exit", command=master.destroy, pady=2)
        self.button2.grid(row=2, sticky=E)
    def context(self):
        self._context1 = self.entry1.get()


root = Tk()
Input = User(root)
root.mainloop()

N=int(Input._context1)

##############################################################################################
#Main Code
import pandas as pd
from pandas import DataFrame, read_csv
from scipy.stats import norm
from scipy import optimize
from scipy.optimize import newton
import numpy as np
import math
import matplotlib.pyplot as plt

location1 = r'/Users/.../settlement.csv'
location2 = r'/Users/.../ICE-SBV2016.csv'
data1 = pd.read_csv(location1)
data1_1=data1[data1['Date']=='30-Jun-2016']
data2 = pd.read_csv(location2)
price=data2['Settle'].values[:]
price_call=data1['PRICE'].values[:]
strike_call=data1['STRIKE'].values[:]

T=float(90)/365
r=0.0123025
#1-year LIBOR rate

num=len(data1_1)
sigma_call=np.zeros((7,num))

##############################################################################################
#Implied Volatility
for i in range(0,7):
    for j in range(0,num):
        def eur_c(vol,r,f,k,t):
            d1 = (math.log(f/k)+(0.5*vol**2)*t)/(math.sqrt(t)*vol)
            d2 = d1-vol*math.sqrt(t)
            call = math.exp(-r*t)*(f*norm.cdf(d1)-k*norm.cdf(d2))
            return call
        eur_call = (lambda vol:eur_c(vol,r,price[i],strike_call[i*num+j],T+float(i)/365)-price_call[i*num+j])
        sigma_call[i][j] = newton(eur_call,0.9,tol=1.48e-08)


##############################################################################################
#American Call Pricing: binomial tree 
#Using control variate 
deltaT=float(T)/N
american_call=np.zeros((7,num))
european_call=np.zeros((7,num))
bs_call=np.zeros((7,num))
american_call_cv=np.zeros((7,num))
for k in range(0,7):
    for j in range(0,num):
        S0=price[k]
        K=strike_call[k*num+j]
        u=np.exp(sigma_call[k][j]*np.sqrt(deltaT))
        d=1/u
        p = (1 - d)/ (u - d)
        oneMinusP = 1.0 - p
        fs_a =  np.asarray([0.0 for i in xrange(N + 1)])
        fs_e =  np.asarray([0.0 for i in xrange(N + 1)])
        fs2 = np.asarray([(S0 * u**i * d**(N - i)) for i in xrange(N + 1)])
        fs3 =np.asarray( [float(K) for i in xrange(N + 1)])
        fs_a[:] = np.maximum(fs2-fs3, 0.0)
        fs_e[:] = np.maximum(fs2-fs3, 0.0)
        for i in xrange(N-1, -1, -1):
            fs_a[:-1]=np.exp(-r * deltaT) * (p * fs_a[1:] + oneMinusP * fs_a[:-1])
            fs2[:]=fs2[:]/d
            fs_a[:]=np.maximum(fs_a[:],fs2[:]-fs3[:])
            fs_e[:-1]=np.exp(-r * deltaT) * (p * fs_e[1:] + oneMinusP * fs_e[:-1])
        american_call[k][j]=fs_a[0]
        european_call[k][j]=fs_e[0]
        american_call_cv[k][j]=fs_a[0]+price_call[k*num+j]-fs_e[0]

##############################################################################################
#Output
from Tkinter import *

strike=strike_call[0:77].astype(str)
day=np.asarray(['JUN 30 2016','JUN 29 2016','JUN 28 2016','JUN 27 2016','JUN 24 2016','JUN 23 2016','JUN 22 2016'])
day=day.astype(str)
for i in range(0,7):
    lst = american_call_cv[i].astype(str)
    
    root = Tk()
    root.title("AMERICAN OPTION PRICING OUTPUT")
    def close_window (): 
        root.destroy()
        
    label_1=Label(root, text="American call price (binomial tree): ")
    label_1.grid(row=0, sticky=W)
    label_1.pack()
    scrollbar = Scrollbar(root)
    scrollbar.pack(side=RIGHT, fill=Y)
    t = Text(root,wrap=WORD, yscrollcommand=scrollbar.set)
    t.insert(END, 'Date: ' + day[i] + '\n' + '\n')
    t.pack()
    for j in range(0,num):
        t.insert(END, 'Strike Price ' + strike[j] +': ' + lst[j] + '\n' + '\n')
        t.pack()
        t.configure(font=("Times New Roman", 14, "bold"))
        t.configure(spacing2=10)
        
    frame = Frame(root)
    frame.pack()
    b=Button(frame, text="Next", command=close_window)
    b.pack()
    
    root.mainloop()

##############################################################################################
#Plotting
for i in range(0,7):
    x1=strike_call[0:77]
    y1=sigma_call[i]
    plt.plot(x1, y1,label='Call',color='g')
    plt.title('Volatility Skew')
    plt.xlabel('Strike')
    plt.ylabel('Implied Volatility')
    plt.legend()
    plt.show()
    plt.close('all')













