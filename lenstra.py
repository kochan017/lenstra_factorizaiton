'''
Purpose: Python Implementation of Lenstra's elliptic curve factoring algorithm.
@author: Ko Yat Chan
'''

'''
Primary function for Lenstra's factoring algorithm.
- input N, a, c, d, B
- N is a positive integer number to be factored,
- a,c,d are variables used in the function d^2 = c^3 + ac + b (mod N)
- B is a practical bound B.
'''
def lenstra(N, a, c, d, B):
    # I-V referring to the summary of the factoring algorithm in the paper.
    # I. Choose a positive integer N to be factorized
    if N < 1:
        print("n has to be a positive integer")
        return -1
    # II. Pick arbitrary numbers for a, c, d mod N
    a, c, d = a % N, c % N, d % N  # ensure the inputs are mod N.
    # III. Assign variable P = (c,d) and than calculate b = d^2 - c^3 - a * c (mod N).
    P = (c, d)
    b = get_b(N, a, c, d)
    # IV. at this point, we get all the variables for the elliptic curve E: d^2 = c^3 + ac + b (mod N)

    # V. init. variable i = 2, create a loop for i < some bound B
    p_list = [P]
    i = 2  # init. i = 2
    o = 0
    while 1:  # must enter the while loop
        print("\n>>>>>> Now computing", i, "P mod", N, "<<<<<<<")
        if i % 2 == 0:  # i is an even number for iP, we use the formula (3x^2+a)/(2y) to get lambda.
            print('---------{}P Computation Starts Here---------'.format(i))
            P = p_list[int(i / 2) - 1]
            x, y = P
            try:
                temp = 2 * y
                recip = pow(temp, -1, N)

                lamb = lambda_fun1(a, x, recip, N)
                print("Inverse of {} mod {} is {}".format(temp, N, recip))
                print("Lambda is {}".format(lamb))

                new_x = (pow(lamb, 2, N) - (2 * x % N)) % N
                new_y = (lamb * (x - new_x) % N - y) % N

                p_list = p_list + [(new_x, new_y)]
                print("{}P is ({},{})".format(i, new_x, new_y))
                print('---------{}P Computation Ends Here-----------'.format(i))

            except ValueError:  # to catch the value error when mod inverse cannot be found
                print("Inverse of {} mod {} cannot be found".format(temp % N, N))
                print("Now we compute d1 and d2")
                d1 = gcd(temp, N)
                if d1 == N:
                    print("d1 == N: You need to repick arbitrary numbers for a,c,d modulo N.")
                    break
                if d1 < N:
                    d2 = int(N / d1)
                    print("The two divisors are d1={} and d2={}, for N={}".format(d1, d2, N))
                    print("Done")
                break

        else:  # i is odd for iP, we use the formula (y2-y1)/(x2-x1) to compute lambda.
            print("in else")
            print('---------{}P Computation Starts Here---------'.format(i))

            P, P2 = p_list[len(p_list) - 2 - o], p_list[len(p_list) - 1 - o]
            x, y = P
            x2, y2 = P2

            print(
                'We use addtion of {}P and {}P for {}P: \nx{} = {}, y{} = {}; x{} = {}, y{} = {}'
                    .format(i - 2 - o,
                            i - 1 - o,
                            i,
                            i - 2 - o,
                            x,
                            i - 2 - o,
                            y,
                            i - 1 - o,
                            x2,
                            i - 1 - o,
                            y2))

            try:
                temp = x2 - x
                recip = pow(temp, -1, N)

                lamb = lambda_fun2(recip, y, y2, N)
                print("Inverse of {} mod {} is {}".format(x2 - x, N, recip))
                print("Lambda is {}".format(lamb))
                new_x = (pow(lamb, 2, N) - x2 - x) % N
                new_y = ((lamb * ((x - new_x) % N) % N - y) % N)
                print("{}P is ({},{})".format(i, new_x, new_y))
                p_list = p_list + [(new_x, new_y)]
                o += 1
                print('---------{}P Computation Ends Here-----------'.format(i))

            except ValueError:  # to catch the value error when mod inverse cannot be found
                print("Inverse of {} mod {} cannot be found".format(temp % N, N))
                print("Now we compute d1 and d2")
                d1 = gcd(temp, N)
                if d1 == N:
                    print("d1 == N: You need to repick arbitrary numbers for a,c,d modulo N.")
                    break
                if d1 < N:
                    d2 = int(N / d1)
                    print("The two divisors are d1={} and d2={}, for N={}".format(d1, d2, N))
                    # print("Done.")
                break
        i += 1
        if not (i < B + 1):
            print("factors cannot be found within the bound", B)
            break


# helper function to calculate b as d^2 - c^3 - a * c (mod N)
def get_b(N, a, c, d):
    return (pow(d, 2, N) - pow(c, 3, N) - a * c) % N


# helper function to evaluate the lambda by (3x^2 + a) / 2y (mod N)
def lambda_fun1(a, x, recip, N):
    return (((3 * pow(x, 2, N) + a) % N) * recip) % N


# helper function to use the formula (y2-y1)/(x2-x1) to compute lambda.
def lambda_fun2(recip, y, y2, N):
    a = (y2 - y) % N
    return (a * recip) % N


# a custom gcd helper function that uses Euclidean Algorithm
def gcd(a, b):
    while 1:
        if not b:
            break
        a, b = b, a % b
    return a


import time

st = time.time()  # start time
lenstra(1081, 1, 0, 1, 6) # simple test case 1
print("It takes {} milliseconds to complete the first test case.".format(time.time() * 1000 - st * 1000))
# print("It takes {} milliseconds to complete the first test case".format((time.perf_counter_ns()-st)/1000000))

st = time.time()  # update start time
lenstra(187, 3, 38, 112, 5) # simple test case 2
print("It takes {} milliseconds to complete the second test case.".format(time.time() * 1000 - st * 1000))
# print("It takes {} milliseconds to complete the second test case".format((time.perf_counter_ns()-st)/1000000))

st = time.time()  # update start time
lenstra(28102844557, 18, 7, 4, 50000) # advanced test case
print("It takes {} milliseconds to complete the advanced test case.".format(time.time() * 1000 - st * 1000))
# print("It takes {} milliseconds to complete the advanced test case".format((time.perf_counter_ns()-st)/1000000))

