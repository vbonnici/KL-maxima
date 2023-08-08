import math

def calculate_kl(d1, d2):
    """
    Calculate the entropic divergence betwen two multiplicity distributions defined on the same quantum, namely their equal in their sum.
    """
    s = sum(d1)
    k = 0.0
    for i in range(len(d1)):
        k += (d1[i]/s) * math.log( (d1[i]/d2[i]), 2 )
    return k

def calculate_nkl(v1, v2):
    """
    Calculate the normalized (in [0...1]) KL divergenge between two multiplicity distributions.
    It is given by NKL(v1|v2) = KL(v1|v2) / KL(v1|M), where M us the distribution for which the maximum divergence can be obtained.
    NKL is not symmetric.
    Distributions are defined as vectors of the same length such that for a generic event i, 
    the multiplcitiy of i in the two distribution is given by v1[i] and v2[i].
    To get a normalized KL, the two distributions must be defined on the same quantum, that is sum(v1)=sum(v2).
    In case the quantum is different, it is obtained as lcm(sum(v1),sum(v2)).
    """
    s1 = sum(v1)
    s2 = sum(v2)
    q = 1
    if s1 != s2:
        lcmv = (s1 * s2) // math.gcd(s1,s2)
        v1 = [ v1[i] * lcmv / s1 for i in range(len(v1)) ]
        v2 = [ v2[i] * lcmv / s2 for i in range(len(v2)) ]
        q = 1 / lcmv
    m = [ q for i in range(len(v1)) ]
    m[ v1.index(min(v1)) ] += sum(v1) - sum(m)
    if calculate_kl(v1,m) != 0:
        return calculate_kl(v1,v2) / calculate_kl(v1,m)
    else:
        return 0


if __name__ == "__main__":
    import random

    print("Testing normalized Kullback-Leibler")
    for i in range(10):
        print('-'*40)
        vecl = random.randint(2,10)
        maxv = random.randint(1,10)
        
        print("vectors length",vecl)
        print("max value",maxv)

        v1 = [ random.randint(1,maxv) for i in range(vecl) ]
        v2 = [ random.randint(1,maxv) for i in range(vecl) ]

        print('vectors')
        print(v1, sum(v1))
        print(v2, sum(v2))

        print("KL(v1|v2) = ", calculate_kl(v1,v2) )
        print("NKL(v1|v2) = ", calculate_nkl( v1,v2 ) )


    print('-'*40)
    print("Manual testing")

    print('-'*40)
    v1 = [5,4,3,2,1]
    v2 = [5,4,3,2,1]
    print('vectors')
    print(v1)
    print(v2)
    print("KL(v1|v2) = ", calculate_kl(v1,v2) )
    print("NKL(v1|v2) = ", calculate_nkl( v1,v2 ) )
    
    print('-'*40)
    v1 = [5,4,3,2,1]
    v2 = [1,2,3,4,5]
    print('vectors')
    print(v1)
    print(v2)
    print("KL(v1|v2) = ", calculate_kl(v1,v2) )
    print("NKL(v1|v2) = ", calculate_nkl( v1,v2 ) )


    print('-'*40)
    v1 = [5,4,3,2,1]
    v2 = [1,1,1,1,11]
    print('vectors')
    print(v1)
    print(v2)
    print("KL(v1|v2) = ", calculate_kl(v1,v2) )
    print("NKL(v1|v2) = ", calculate_nkl( v1,v2 ) )

    

