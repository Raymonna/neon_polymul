#include <stdint.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>

// Function to perform modular exponentiation
int16_t mod_exp(int16_t base, int16_t exp, int16_t mod) {
    int16_t result = 1;
    while (exp > 0) {
        if (exp % 2 == 1) {
            result = (result * base) % mod;
        }
        base = (base * base) % mod;
        exp = exp / 2;
    }
    return result;
}

// Function to check if a given candidate is a primitive root of unity
bool is_primitive_root_of_unity(int16_t candidate, int16_t n, int16_t q) {
    // Check if candidate^n == 1 (mod q)
    if (mod_exp(candidate, n, q) != 1) {
        return false;
    }

    // Check if candidate^(n/d) != 1 (mod q) for all proper divisors d of n
    for (int i = 2; i * i <= n; i++) {
        if (n % i == 0) {
            if (mod_exp(candidate, n / i, q) == 1) {
                return false;
            }
            if (i != 1 && i != n && mod_exp(candidate, i, q) == 1) {
                return false;
            }
        }
    }

    return true;
}

// Function to find a primitive root of unity
int16_t find_primitive_root_of_unity(int16_t n, int16_t q) {
    // Check the properties of q
    if ((q - 1) % n != 0) {
        return -1; // q-1 must be divisible by n
    }

    for (int16_t candidate = 2; candidate < q; candidate++) {
        if (is_primitive_root_of_unity(candidate, n, q)) {
            return candidate;
        }
    }
    return -1; // Return -1 if no primitive root of unity is found
}

int main() {
    int16_t n = 2048; // Example power of 2
    int16_t q = 4591; // Example modulus

    int16_t root_of_unity = find_primitive_root_of_unity(n, q);
    if (root_of_unity != -1) {
        printf("Primitive root of unity: %d\n", root_of_unity);
    } else {
        printf("No primitive root of unity found.\n");
    }

    return 0;
}

