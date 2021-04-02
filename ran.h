#ifndef _RAN_H_
#define _RAN_H_

#define Ullong unsigned long long
#define Uint unsigned int
#define Doub double

//Gan note 10/25/2019: this is the highest uality RNG from numerical recipe.
//recommended by Prof. Erik Luijten

struct Ran {
/*Implementation of the highest quality recommended generator. 
The constructor is called with an integer seed and creates an 
instance of the generator. The member functions int64, doub,
and int32 return the next values in the random sequence, 
as a variable type indicated by their
names. The period of the generator is 3.138^1057.*/
Ullong u,v,w;
Ran(Ullong j) : v(4101842887655102017LL), w(1) {
//Constructor. Call with any integer seed (except value of v above).
u = j ^ v; int64();
v = u; int64();
w = v; int64();
}
inline Ullong int64() {
u = u * 2862933555777941757LL + 7046029254386353087LL;
v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
w = 4294957665U*(w & 0xffffffff) + (w >> 32);
Ullong x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
return (x + v) ^ w;
}

inline Doub doub() { return 5.42101086242752217E-20 * int64(); }
//Return random double-precision floating value in the range 0. to 1.
inline Uint int32() { return (Uint)int64(); }
//Return 32-bit random integer.
};

#endif
