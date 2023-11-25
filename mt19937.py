"""
A Python version of MT19937
Original C code from Takuji Nishimura and Makoto Matsumoto (mt19937ar.c)
SAP - JIL by JHL in November 22, 2023

Before using, initialize the state by using init_genrand(seed) or init_by_array(init_key, key_length).
"""
# 귀찮으니 이것도 마저 이식
from numpy import uint as unslong
from numpy import int_ as long
from numpy import double

# - Period parameters
N = 624
M = 397
MATRIX_A = unslong(0x9908b0df)  # constant vector a
UPPER_MASK = unslong(0x80000000)  # most significant w-r bits
LOWER_MASK = unslong(0x7fffffff)  # least significant r bits
mt = [unslong(0)] * N  # the array for the state vector - static!
mti: int = N + 1  # mti==N+1 means mt[N] is not initialized - static!


# - initializes mt[N] with a seed
def init_genrand(s: unslong):
    global N, M, MATRIX_A, UPPER_MASK, LOWER_MASK, mt, mti
    mt[0] = s & unslong(0xffffffff)  # s & 0xffffffffUL
    for mti in range(1, N):
        # - See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier.
        mt[mti] = (unslong(1812433253) * (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + mti)
        mt[mti] &= unslong(0xffffffff)  # for >32 bit machines


# - initialize by an array with array-length
def init_by_array(init_key: unslong, key_length: int):
    # init_key is the array for initializing keys
    # key_length is its length
    global N, M, MATRIX_A, UPPER_MASK, LOWER_MASK, mt, mti
    init_genrand(unslong(19650218))
    i = 1
    j = 0
    k = N if N > key_length else key_length  # (N>key_length ? N : key_length)

    for k in range(k, 0, -1):
        mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * unslong(1664525))) + init_key[j] + j  # non linear
        mt[i] &= unslong(0xffffffff)
        i += 1
        j += 1
        if i >= N:
            mt[0] = mt[N - 1]
            i = 1
        if j >= key_length:
            j = 0

    for k in range(N - 1, k, -1):
        mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * unslong(1566083941))) - i  # 위와 같음
        mt[i] &= unslong(0xffffffff)
        i += 1
        if i >= N:
            mt[0] = mt[N - 1]
            i = 1
        mt[0] = unslong(0x80000000)  # MSB is 1; assuring non-zero initial array


# - generates a random number on [0,0xffffffff]-interval
def genrand_int32():
    global N, M, MATRIX_A, UPPER_MASK, LOWER_MASK, mt, mti
    y: unslong
    mag01 = [unslong(0x0), MATRIX_A]  # mag01[x] = x * MATRIX_A  for x=0,1

    if mti >= N:  # generate N words at one time
        kk: int

        if mti == N + 1:  # if init_genrand() has not been called
            init_genrand(unslong(5489))  # a default initial seed is used

        for kk in range(0, N - M):
            y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK)
            mt[kk] = mt[kk + M] ^ (y >> 1) ^ mag01[y & unslong(0x1)]

        for kk in range(0, N - 1):
            y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK)
            mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ mag01[y & unslong(0x1)]

        y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK)
        mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ mag01[y & unslong(0x1)]

        mti = 0

    y = mt[mti]
    mti += 1

    # - Tempering
    y ^= (y >> 11)
    y ^= (y << 7) & unslong(0x9d2c5680)
    y ^= (y << 15) & unslong(0xefc60000)
    y ^= (y >> 18)

    return y


# - generates a random number on [0,0x7fffffff]-interval
def genrand_int31():
    return long(genrand_int32() >> 1)


# - generates a random number on [0,1]-real-interval
def genrand_real1():
    return double(genrand_int32() * (1.0 / 4294967295.0))  # divided by 2^32-1


# - generates a random number on [0,1)-real-interval
def genrand_real2():
    return double(genrand_int32() * (1.0 / 4294967296.0))  # divided by 2^32


# - generates a random number on (0,1)-real-interval
def genrand_real3():
    return double((double(genrand_int32())) + 0.5) * (1.0 / 4294967296.0)  # divided by 2^32


# - generates a random number on [0,1) with 53-bit resolution
def genrand_res53():
    a: unslong = genrand_int32() >> 5
    b: unslong = genrand_int32() >> 6
    return double((a * 67108864.0 + b) * (1.0 / 9007199254740992.0))

'''
# - These real versions are due to Isaku Wada, 2002/01/09 added
int main(void)
{
    int i;
    unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4;
    init_by_array(init, length);
    printf("1000 outputs of genrand_int32()\n");
    for (i=0; i<1000; i++) {
      printf("%10lu ", genrand_int32());
      if (i%5==4) printf("\n");
    }
    printf("\n1000 outputs of genrand_real2()\n");
    for (i=0; i<1000; i++) {
      printf("%10.8f ", genrand_real2());
      if (i%5==4) printf("\n");
    }
    return 0;
}
'''