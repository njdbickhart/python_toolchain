# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 15:46:16 2022

@author: derek.bickhart-adm
"""

def solution(N):
    
    permutations = 0
    prevlen = 0
    for i in range(1, int(N / 10)):
       p, l = stepPermutator(N, i + 1)
       if l != prevlen:
           permutations += p
           prevlen = l    
    
            
    return permutations

def stepPermutator(N, maxS):
    stairs = widestStair(N, maxS)
    print(stairs)
    print(sum(stairs))
    prevWidest = len(stairs)
    permutations = 1
    val = stairs[-1] - stairs[-2] - 1
    
    while val > 0 and len(stairs) > 2:
        stairs = widestStair(val, len(stairs) - 2)
        print(f'i{len(stairs)} {stairs}')
        if len(stairs) == 1:
            permutations += stairs[0]
            break
        val = stairs[-1] - stairs[-2] - 1
        permutations += int(len(stairs) / (prevWidest - 2))
    print(permutations)
    return permutations, len(stairs)

def widestStair(N, maxS):
    total = 0
    stairs = list()
    for i in range(0, N):
        if len(stairs) == 0:
            stairs.append(1)
            total += 1
            continue
        if total + i + 1 <= N and i + 1 > stairs[-1] and len(stairs) < maxS:            
            stairs.append(i + 1)
            total += i + 1
        elif total < N:
            stairs[-1] += 1
            total += 1
    return stairs
    
print(solution(200))