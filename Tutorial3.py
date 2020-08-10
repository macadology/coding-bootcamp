#!/usr/bin/env python
# coding: utf-8

# # Functions in python
# 
# Functions help to organize a script into modular chunks
# 
# Functions are especially useful for repetitive, complex tasks
# 
# 
# ## Basic syntax
# 
# def function_name(parameters):
# 	statement(s)
#     #commonly print or return something
# 
# 

# What is the difference between the print and return commands in functions?
# 


def absolute_value(num):
    """This function returns the absolute
    value of the entered number"""

    if num >= 0:
        print(num)
    else:
        print(-num)

        
absolute_value(2)

absolute_value(-4)

test=absolute_value(2)

type(test)


def absolute_value(num):
    """This function returns the absolute
    value of the entered number"""

    if num >= 0:
        return num
    else:
        return -num

        
absolute_value(2)

absolute_value(-4)

test=absolute_value(2)

type(test)


# Function to count the number of bases in a string of nucleotides

def base_count(dna, base):
    i = 0 
    for c in dna:
        if c == base:
            i += 1
    return i


result=base_count("AATTAGCCTTA", "A")

print(result)

#More sophisticated printing?

print('%s appears %d times in %s' % (base, result, dna)) #Why does this break? Due to variable scoping.

# Function to count the number of bases in a string of nucleotides


def base_count(dna, base):
    i = 0 
    for c in dna:
        if c == base:
            i += 1
    result=i
    #More sophisticated printing
    print('%s appears %d times in %s' % (base, result, dna)) 
    
    print('{} appears {} times in {}'.format(base, result, dna))
    
    print(f'{base} appears {result} times in {dna}') #using fstrings


base_count("AATTAGCCTTA", "A")

# Function to count the number of times each unique base appears in a string of nucleotides


def base_count_v2(dna):
    all_freq = {} 
    for i in dna: 
        if i in all_freq: 
            all_freq[i] += 1
        else: 
            all_freq[i] = 1
    
    print(all_freq)
    

base_count_v2("AATTAGCCTTA")


# ## How to test if my function is working as intended? 
# 
# Positive controls, defined input should give expected output (e.g. assert in python)
# Printing variables and messages throughout function body can help in debugging
# 


#Model random mutations in a string of nucleotides

import random

random.seed(a=1000)

def mutate(dna):
    dna_list = list(dna)
    mutation_site = random.randint(0, len(dna_list) - 1)
    dna_list[mutation_site] = random.choice(list("ATCG"))
    return "".join(dna_list)


dna = "ACGGAGATTTCGGTATGCAT"
print("Starting DNA:", dna)
print(base_count_v2(dna))

nmutations = 10000

for i in range(nmutations):
    dna = mutate(dna)

print("DNA after %d mutations:" % nmutations, dna)
print(base_count_v2(dna))


# # List and dictionary comprehensions
# 
# Compact way of creating lists and dictionaries


##making a list by typing everything out is impractical

my_list=[1,2,3,4,5,6,6,8,9,10,11,12,13,14,15]

#a list comprehension is more compact

my_list2=[i for i in range(1,16,1)]

print(my_list)
print(my_list2)

#Alternating even and odd labels

my_list3 = ["Even" if i%2==0 else "Odd" for i in range(10)]
print(my_list3)



##Dictionary comprehensions

#Minimal syntax for dictionary comprehensions
#dictionary = {key: value for vars in iterable}



#Recall Dictionaries are data types in Python which allows us to store data in key/value pairs.

my_dict = {"Adenine":'A', "Thymine":'T', "Guanine":'G', "Cytosine":'C', "Uracil":"U"}

print(my_dict)


#Dictionary comprehensions allow us to generate new dictionaries in a compact manner

my_dict2 = {k: ('purine' if v in ("A","G")  else 'pyrimidine') for k,v in my_dict.items()}

print(my_dict2)

#An example of generating a dictionary matching each letter in the English alphabet to a number

import string

numbers=[i for i in range(1,27)]

#Python's zip() function creates an iterator that will aggregate elements from two or more iterables. 
object= zip(numbers, string.ascii_uppercase) 

my_dict3 = {k:v for k,v in object}

print(my_dict3)



# ## Exercises 
# 
# 
# Ex 1)
# Use list comprehension to obtain a list of the squares of the first 10 positive integers. 
# 
# Ex 2) 
# Write a Python function that that prints out the first n rows of Pascal's triangle.
# 
# Expected output for 6 rows:
# [1]                                                                                                           
# [1, 1]                                                                                                        
# [1, 2, 1]                                                                                                     
# [1, 3, 3, 1]                                                                                                  
# [1, 4, 6, 4, 1]                                                                                               
# [1, 5, 10, 10, 5, 1]  
# 
# 
# Ex 3)
# Given a nucleotide sequence: 
# ATGCAAGCGGATATGCATGGAAAACTTCACGCTGCCTTAGAAGATGGTTTCTTCCTCTTTTTTTGAACAGCAACAACAACCTAACATTTATTATGACACAACCACCGATCAAGAAGAC
# 
# Write a function count_substring(dna, substring), to count how many times a certain substring appears in the string of nucleotides.
# For example, the function returns 6 when called with the substring 'CAA'and the nucleotide string above. 
# 
# 
# 
