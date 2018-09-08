import numpy as np
import csv
import math
from scipy import stats
from scipy.stats import norm
import pandas as pd

def bs_call(S, K, r, sigma, t):
    d1 = (math.log(S/K) + (r  + 0.5 * sigma*sigma)*t)/(sigma*t**0.5)
    d2 = (math.log(S/K) + (r  - 0.5 * sigma*sigma)*t)/(sigma*t**0.5)
    return math.exp(-t)*S*norm.cdf(d1) - math.exp(-r * t) * K *norm.cdf(d2)

def bs_put(S, K, r, sigma, t):
    d1 = (math.log(S/K) + (r  + 0.5 * sigma*sigma)*t)/(sigma*t**0.5)
    d2 = (math.log(S/K) + (r  - 0.5 * sigma*sigma)*t)/(sigma*t**0.5)
    return -math.exp(-t)*S*norm.cdf(-d1) + math.exp(-r * t) * K *norm.cdf(-d2)

def dividends(S, N):
    sum = 0
    for n in range(N):
        sum = sum + 0.008306 * (1-0.008306)**(n-1)
    return sum

df1=pd.read_csv('call.csv', sep=',')
df=pd.read_csv('pro1_call.csv', sep=',')
df2 = pd.read_csv('pro1_put.csv', sep = ',')
print(df)
x = df['days']/365
y1 = df['45940']
y2 = df['45950']
y3 = df['45960']
# regress and predict implied volatilities on 20170825, then regard them as ivs on 20200825
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y1)
iv_1 = slope * 2 + intercept
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y2)
iv_2 = slope * 2 + intercept
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y3)

iv_3 = slope * 2 + intercept
x = df1['days']/365
slope, intercept, r_value, p_value, std_err = stats.linregress(x,df1['61260'])
iv_4 = slope * 2 + intercept

sigma = 0.45473
K1 = 45.94
K2 = 45.95
K3 = 45.96
N = 20
S = 61.26
t = 5
r = 1.631192132/100
d= dividends(S, N)
print(d)
digital_call_1 = bs_call(S-d, K1, r, iv_1, t)
digital_call_2 = bs_call(S-d, K1, r, iv_1, t)
digital_call_3 = bs_call(S-d, K1, r, iv_1, t)


sub = (digital_call_2-digital_call_3)
super = (digital_call_1-digital_call_2)

call = bs_call(S-d, 61.26, r, iv_4, t)

digital_call = (sub + super)/2

portfolio_super = 7.5/45.95*(61.26-d)-7.5/45.95* digital_call_1 + 2.5 * super * 100 + 16.379/61.26  * call
portfolio_sub = 7.5/45.95*(61.26-d)-7.5/45.95* digital_call_2 + 2.5 * sub * 100 + 16.379/61.26  * call
print(portfolio_super,portfolio_sub)

slope, intercept, r_value, p_value, std_err = stats.linregress(df2['days'],df2['45950'])
iv_5 = slope * 2 + intercept
put_super = bs_put(S-d, K1, r, iv_5, t)
put_sub = bs_put(S-d, K2, r, iv_5, t)

portfolio2_super = -7.5/45.95 * put_super +7.5/(1+1.631192132/100)**5 + 2.5 * super * 100 + 16.379/61.26 * call
portfolio2_sub = -7.5/45.95 * put_sub +7.5/(1+1.631192132/100)**5 + 2.5 * sub * 100 + 16.379/61.26 * call
print(portfolio2_super,portfolio2_sub)
