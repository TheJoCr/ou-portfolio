from os import fork

def gen_summands( d, n=100 ):
    p = Primes().first()
    primes = []
    for _ in range(n):
        primes.append(p)
        p = Primes().next(p)

    return map(lambda p: floor(10**d * p**(1/3)), primes)

def target( d, n=100 ):
    return 2 * 113085286 / 10**8 * n * 10**d

def find_subset_sum( summands, M ):
    target_vec = copy(summands)
    target_vec = map(lambda a: a*-1, target_vec)
    target_vec.append(M)

    a = MatrixSpace(ZZ,len(summands) + 1)(1)
    a[-1,:] = 0
    a[:,-1] = vector(target_vec)

    r = a.BKZ()

    for row in r.rows():
        works = True
        l = 0
        for val in row:
            if val != 0:
                if l == 0:
                    l = val
                elif l != val:
                    works = False
        if works:
            return list(row / l)

    return None

def search_d(d):
    print("D = " + str(d))
    sys.stdout.flush()
    for i in range(1000000):
        summands = gen_summands(d, 100)
        targ = target(d, 100) + i
        sol = find_subset_sum( summands, targ )
        if sol != None:
            print( "Sol:", sol )
            print( "Ms:", vector(summands) )
            print( "Sum:", vector(summands) * vector(sol[:-1]) )
            print( "Target:", targ )
            break

if __name__ == '__main__':
        max_proc_num = 20
        init_d = 10
        print("Initalizing...")
        import sys
        sys.stdout.flush()
        while max_proc_num > 0 and fork() != 0:
            max_proc_num -= 1

        start = init_d + max_proc_num
        print("Starting proc number " + str(max_proc_num) + " at dim " + str(start))
        sys.stdout.flush()
        for d in range(start, 100, 20):
            search_d(d)
            sys.stdout.flush()
