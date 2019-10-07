from math import *

def slope(b):
    return 1/(b-1)*log(((math.pi*b)**(1./b)*b/(2*math.pi*math.e)),2)


def svp_plausible(b):
	return b *log(sqrt(4./3))/log(2) +log(b,2)  # .2075 * b 

def svp_quantum(b):
	return b *log(sqrt(13./9))/log(2) +log(b,2) # .265 * b  

def svp_classical(b):
	return b * log(sqrt(3./2))/log(2)  +log(b,2) # .292 * b 

def rootHF(b):
    return 2 ** (0.5 * slope(b))
    
# def num_of_q(n, r1):
#     return n - r1
    
def dim_of_extracted(n, r1, r2):
    return 2 * n + 1 - r1 - r2


def term_of_alpha(j, b, first_term, r1):
    return first_term + slope(b) * (j - r1)

    
def hybrid_cost(n, q, coef=svp_classical):
    best_cost = 9999
    best_r1=1
    for b in range(200, 2 * n + 1):
        cb = coef(b)
        r2 = int(round(2 * cb / log(3, 2)))
        for r1 in range(1, 2 * n + 1 - r2):
            d = dim_of_extracted(n, r1, r2)
            delta = rootHF(b)
            last_alpha = (n-r1)/d-d*log(delta,q)
            if last_alpha > log(2., q):
                cost_wrt_r1 = log(2 ** cb + 3 ** (.5 * r2), 2)
                if cost_wrt_r1 < best_cost:
                    best_cost = cost_wrt_r1
                    best_r1 = r1
        if best_cost != 9999: break
    return b, dim_of_extracted(n, best_r1, r2), best_r1, best_cost

def dilithium_cost(n, q, a, fun=svp_quantum):
    best_cost = 9999
    best_b = 9999
    best_w = 9999
    best_r2 = 9999
    for w in range(n, 2 * n + 1):
        # forgetting q-vector, then r1=0
        for b in range(200, w):
            r2 = int(floor(-0.5 + sqrt(.25 + 2 * n * log(q,2) / slope(b)) + 1))
            l1 = slope(b) * (r2 - 1)
            # l1 = n * log(q, 2) / (r2 - 1) + 0.5 * (r2 - 2) * slope(b)
            if r2 - 1 > w: r2 = w + 1
            std = 2 ** l1 / sqrt(r2 - 1)
            log_p = log(erf(a / (sqrt(2) * std)), 2)
            cost = max(fun(b), fun(b) - (r2 - 1) * log_p - b*log(sqrt(4./3), 2))
            if cost < best_cost:
                best_cost = cost
                best_b = b
                best_w = w
                best_r2 = r2
    return best_b, best_w, best_r2, best_cost
           
for (n, q) in [(1024, 463873)]:
    b, N, r1, bitsec = hybrid_cost(n, q)
    print "blocksize\tN\tr1\t bitsec"
    print b, "&", N, "&", r1, "&", bitsec
    # a = (q - 1) / 8
    # b, w, r2, cost = dilithium_cost(n, q, a)
    # print "blocksize\tclassical\tquantum\tplausible"
    # print b, "&", svp_classical(b), "&", svp_quantum(b), "&", svp_plausible(b)



               
             
            
            

            
        
        

