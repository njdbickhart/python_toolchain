# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 14:40:39 2022

@author: derek.bickhart-adm
"""
import string

def solution(n, b):
    prev = list()
    cycle = task(n, b, prev)
        
    return(cycle)
    
def task(n, b, prev):
    k = len(str(n))
    x = intSort(n, True)
    y = intSort(n, False)
    z = int(x, b) - int(y, b)
    z = intToBase(z, b)
    z = list(str(z))
    if len(z) != k:
        for i in range(k - len(z)):
            z = ['0'] + z
    z = ''.join(z)
    prev.append(z)
    
    cycle = 0
    if len(prev) > 1:
        for i in range(len(prev) - 2, -1, -1): 
            cycle += 1
            if int(prev[i], b) == int(prev[-1], b):
                return cycle
    return task(z, b, prev)
    
def intSort(n, reverse):
    l = list(str(n))
    l.sort(reverse=reverse)
    return ''.join(l)

def intToBase(n, b):
    base_n_digits = string.digits + string.ascii_lowercase + string.ascii_uppercase
    result = ""
    while n > 0:
        q, r = divmod(n,b)
        result += base_n_digits[r]
        n = q
    if result == "":
        result = "0"
    return "".join(reversed(result))
    

print(solution(210022, 3))
print(solution(1211, 10))